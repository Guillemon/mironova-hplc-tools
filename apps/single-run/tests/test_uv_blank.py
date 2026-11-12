from pathlib import Path
import tempfile
import unittest
import json
import numpy as np
from hplc_autoquant.single_run import RunSettings
from hplc_autoquant.uv_blank import subtract_uv_trace,matching_blank,blank_corrected_analysis,uv_display_traces


class UVBlankTests(unittest.TestCase):
    def test_display_blank_is_recorded_response_even_with_zero_scale(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);x=np.linspace(0,4,5)
            np.savetxt(root/'sample.csv',np.column_stack([x,x+10,x*.1]),delimiter=',',header='t,y,b',comments='')
            np.savetxt(root/'corrected.csv',np.column_stack([x,x+10,x*.2]),delimiter=',',header='t,y,b',comments='')
            np.savetxt(root/'blank.csv',np.column_stack([x,x+20]),delimiter=',',header='t,y',comments='')
            channel={'trace_csv':'corrected.csv','blank_correction':{'applied':True,'original_trace_csv':'sample.csv','blank_channel':'UV','scale':0.,'shift_min':.2}}
            traces=uv_display_traces(root,channel,root,{'channels':[{'channel':'UV','trace_csv':'blank.csv'}]})
            self.assertEqual(list(traces),['Original sample','Blank','Baseline','Blank-corrected sample'])
            np.testing.assert_allclose(traces['Blank'][:,1],x+20)
            np.testing.assert_allclose(traces['Blank'][:,0],x+.2)
            np.testing.assert_allclose(traces['Baseline'][:,1],x*.2)
            np.testing.assert_allclose(traces['Original sample'][:,1],x+10)

    def test_pointwise_subtraction_retains_negative_values(self):
        x=np.linspace(0,4,401);background=5+20*np.exp(-.5*((x-1)/.08)**2);analyte=50*np.exp(-.5*((x-3)/.06)**2)
        xx,result,scaled,mask=subtract_uv_trace(x,background+analyte,x,background)
        np.testing.assert_allclose(result,analyte,atol=1e-12)
        _,negative,_,_=subtract_uv_trace(x,np.zeros_like(x),x,background)
        self.assertTrue(np.all(negative<0))

    def test_overlap_only_and_explicit_time_shift(self):
        x=np.linspace(0,10,101);bx=np.linspace(1,9,81)
        xx,result,background,_=subtract_uv_trace(x,x,bx,bx,shift_min=1.)
        self.assertGreaterEqual(xx[0],2);self.assertLessEqual(xx[-1],10)
        np.testing.assert_allclose(result,np.ones(len(xx)))
        with self.assertRaises(ValueError):subtract_uv_trace(x,x,bx,bx,scale=-1.)

    def test_wavelength_and_optics_must_match(self):
        sample={'metadata':{'signal':'DAD1C, Sig=210,4 Ref=360,100','bandwidth':4.}}
        wrong={'metadata':{'signal':'DAD1A, Sig=250,4 Ref=360,100','bandwidth':4.}}
        right={'metadata':{'signal':'DAD1B, Sig=210,4 Ref=360,100','bandwidth':4.}}
        self.assertIsNone(matching_blank(sample,[wrong]))
        self.assertIs(matching_blank(sample,[wrong,right]),right)
        with self.assertRaises(ValueError):matching_blank(sample,[right,right])

    def test_background_peak_is_removed_uv_reintegrated_ms_untouched(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);base=root/'base';blank=root/'blank'
            (base/'data').mkdir(parents=True);(base/'source_snapshots').mkdir()
            (blank/'data').mkdir(parents=True);(blank/'sources').mkdir()
            x=np.linspace(0,4,4001);solvent=10+100*np.exp(-.5*((x-1)/.04)**2);analyte=70*np.exp(-.5*((x-3)/.04)**2)
            np.savetxt(base/'data/uv.csv',np.column_stack([x,solvent+analyte,np.ones(len(x))*10]),delimiter=',',header='time,signal,baseline',comments='')
            np.savetxt(blank/'data/uv.csv',np.column_stack([x,solvent]),delimiter=',',header='time,signal',comments='')
            np.savetxt(base/'data/tic.csv',np.column_stack([x,solvent+analyte,np.zeros(len(x))]),delimiter=',',header='time,tic,baseline',comments='')
            np.savez_compressed(base/'data/ms.npz',time_min=x,mz=np.array([100.,200.]),intensity=np.column_stack([solvent+analyte,solvent]))
            original_ms=(base/'data/ms.npz').read_bytes()
            metadata={'signal':'DAD1C, Sig=210,4 Ref=360,100','unit':'mAU'}
            record={'source_path':str(root/'A.D'),'sources':[],'review_notes':[],
                    'uv_channels':[{'channel':'DAD1C','metadata':metadata,'trace_csv':'data/uv.csv','peaks':[],'area_unit':'mAU*s'}],
                    'ms_channels':[{'channel':'MSD1','trace_csv':'data/tic.csv','matrix_export':'data/ms.npz','peaks':[]}]}
            (base/'run_summary.json').write_text(json.dumps(record));(blank/'sources/DAD1C.ch').write_bytes(b'blank snapshot')
            info={'source_path':str(root/'B.D'),'injection':'B.D','signature':'test','channels':[{'channel':'DAD1C','file':'DAD1C.ch','metadata':metadata,'trace_csv':'data/uv.csv','sha256':'blank hash','source_path':str(root/'B.D/DAD1C.ch')} ]}
            out,result,_=blank_corrected_analysis(base,blank,info,root/'results',RunSettings())
            peaks=result['uv_channels'][0]['peaks']
            self.assertEqual(len(peaks),1);self.assertAlmostEqual(peaks[0]['rt_min'],3.)
            self.assertFalse(result['blank_correction']['ms_modified'])
            self.assertEqual((out/'data/ms.npz').read_bytes(),original_ms)
            self.assertEqual(result['ms_channels'][0]['peaks'],[])
            self.assertTrue((out/result['uv_channels'][0]['blank_correction']['original_trace_csv']).exists())
            self.assertTrue(result['mass_associations'])
            # Changing the window beyond the blank endpoint must still allow
            # valid corrected peaks, without inventing blank data past 4 min.
            out,extended,_=blank_corrected_analysis(base,blank,info,root/'results',RunSettings(integration_start_min=.5,integration_end_min=4.2))
            uv=extended['uv_channels'][0]
            self.assertEqual(uv['blank_correction']['effective_integration_min'],[.5,4.])
            self.assertEqual(len(uv['peaks']),1)
            self.assertAlmostEqual(uv['peaks'][0]['rt_min'],3.)
            self.assertEqual(np.loadtxt(out/uv['trace_csv'],delimiter=',',skiprows=1)[-1,0],4.)
            self.assertEqual((out/'data/ms.npz').read_bytes(),original_ms)
            _,outside,_=blank_corrected_analysis(base,blank,info,root/'results',RunSettings(integration_start_min=5.,integration_end_min=6.))
            self.assertEqual(outside['uv_channels'][0]['peaks'],[])
            self.assertIsNone(outside['uv_channels'][0]['blank_correction']['effective_integration_min'])


if __name__=='__main__':unittest.main()
