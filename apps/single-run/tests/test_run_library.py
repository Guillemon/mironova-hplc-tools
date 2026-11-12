from dataclasses import replace
from pathlib import Path
import json
import tempfile
import unittest
from unittest.mock import patch
import numpy as np

from hplc_autoquant.single_run import RunSettings,detect_integrate
from hplc_autoquant.run_library import scan_library,filter_library,cached_analysis,display_summary,saved_library,reintegrate_saved,selected_library_index


class RunLibraryTests(unittest.TestCase):
    def test_any_library_cell_selects_its_injection(self):
        self.assertEqual(selected_library_index({'selection':{'cells':[[2,'Sample']]}},4),2)
        self.assertEqual(selected_library_index({'selection':{'rows':[1]}},4),1)
        self.assertIsNone(selected_library_index({'selection':{'cells':[[8,'Day']]}},4))
        self.assertIsNone(selected_library_index({'selection':{}},4))
    def test_status_change_without_content_change_is_accepted(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=root/'A.D';source.mkdir();(source/'DAD1A.ch').write_bytes(b'stable data')
            def process(source,out,settings,generate_pdf):
                out.mkdir(parents=True);(out/'trace.csv').write_text('decoded')
                summary={'uv_channels':[],'ms_channels':[],'settings':{}}
                (out/'run_summary.json').write_text(json.dumps(summary));(out/'manifest.json').write_text(json.dumps({'files':{'trace.csv':'hash'}}))
                return out,summary
            before=[['DAD1A.ch',11,100,100]];after=[['DAD1A.ch',11,100,200]]
            with patch('hplc_autoquant.run_library.file_stamp',side_effect=[before,after,after]),patch('hplc_autoquant.run_library.process_injection',side_effect=process):
                folder,_,_=cached_analysis(source,root/'results')
            with patch('hplc_autoquant.run_library.file_stamp',return_value=after),patch('hplc_autoquant.run_library.source_fingerprint',side_effect=AssertionError('Timestamp was not recorded')):
                reopened,_,hit=cached_analysis(source,root/'results')
                self.assertTrue(hit);self.assertEqual(reopened,folder)

    def test_content_change_during_processing_is_still_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=root/'A.D';source.mkdir();(source/'DAD1A.ch').write_bytes(b'original')
            def process(source,out,settings,generate_pdf):
                out.mkdir(parents=True);(Path(source)/'DAD1A.ch').write_bytes(b'changed!')
                return out,{}
            with patch('hplc_autoquant.run_library.file_stamp',side_effect=[['before'],['after']]),patch('hplc_autoquant.run_library.process_injection',side_effect=process):
                with self.assertRaisesRegex(ValueError,'contents changed'):cached_analysis(source,root/'results')

    def test_nested_library_search_utf16_and_day(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);run=root/'sequence 2026-08-27'/'Sample 1.D';run.mkdir(parents=True)
            (run/'SAMPLE.XML').write_bytes('<Sample><Name>CEA lot ABC</Name><Description>MeOH blank</Description><ACQMethodPath>TAF-B01</ACQMethodPath></Sample>'.encode('utf-16'))
            (run/'acq.txt').write_bytes('Acq. Method: TAF-B01'.encode('utf-16'))
            (run/'ignored.D').mkdir()  # Do not descend into injections.
            other=root/'sequence 2026-11-06'/'Sample 2.d';other.mkdir(parents=True)
            rows,warnings=scan_library(root)
            self.assertEqual(len(rows),2)
            self.assertEqual(rows[0]['day'],'2026-11-06')
            self.assertEqual(len(filter_library(rows,'CEA MeOH TAF-B01')),1)
            self.assertEqual(len(filter_library(rows,day='2026-08-27')),1)
            self.assertEqual(len(filter_library(rows,folder='sequence 2026-08-27')),1)
            self.assertEqual(filter_library(rows,'unmatched'),[])
            self.assertEqual(rows[1]['date_source'],'Folder name')

    def test_listing_opens_only_sample_xml_and_filenames_mode_opens_nothing(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);run=root/'A.D';run.mkdir()
            for name in ['DAD1A.ch','MSD1.MS','DAD1.UV','acq.txt']:(run/name).write_bytes(b'Never read this content during listing')
            (run/'SAMPLE.XML').write_text('<Sample><Name>ABC</Name><Description>Comment</Description></Sample>')
            opened=[];original=Path.open
            def tracked(path,*args,**kwargs):
                opened.append(path.name)
                if path.name!='SAMPLE.XML':raise AssertionError('Listing opened instrument data')
                return original(path,*args,**kwargs)
            with patch.object(Path,'open',new=tracked):
                rows,_=scan_library(root)
                self.assertEqual(opened,['SAMPLE.XML'])
                self.assertEqual(rows[0]['sample'],'ABC')
                self.assertEqual(len(filter_library(rows,'Comment')),1)
                opened.clear()
                rows,_=scan_library(root,include_sample_metadata=False)
                self.assertEqual(opened,[])
                self.assertIn('MSD1.MS',rows[0]['files'])

    def test_cache_reopens_and_content_change_invalidates(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=root/'A.D';source.mkdir();(source/'DAD1A.ch').write_bytes(b'raw1')
            def fake_process(source,out,settings,generate_pdf):
                self.assertFalse(generate_pdf);out.mkdir(parents=True)
                (out/'trace.csv').write_text('trace')
                result={'uv_channels':[{'trace_csv':'trace.csv'}],'ms_channels':[],'settings':{}}
                (out/'run_summary.json').write_text(json.dumps(result));(out/'manifest.json').write_text(json.dumps({'files':{'trace.csv':'hash'}}))
                return out,result
            with patch('hplc_autoquant.run_library.process_injection',side_effect=fake_process) as process:
                first,_,hit=cached_analysis(source,root/'results');self.assertFalse(hit)
                second,_,hit=cached_analysis(source,root/'results',RunSettings(filter_report=True,report_min_area_pct=2));self.assertTrue(hit)
                with patch('hplc_autoquant.run_library.source_fingerprint',side_effect=AssertionError('Raw content reread')):
                    cached_analysis(source,root/'results')
                self.assertEqual(first,second);self.assertEqual(process.call_count,1)
                (source/'DAD1A.ch').write_bytes(b'raw2') # Same size; content hash still detects change.
                third,_,hit=cached_analysis(source,root/'results');self.assertFalse(hit);self.assertNotEqual(first,third)
                (third/'trace.csv').unlink()
                _,_,hit=cached_analysis(source,root/'results');self.assertFalse(hit)
                self.assertEqual(process.call_count,3)

    def test_saved_index_survives_sessions_and_refresh_skips_unchanged_metadata(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)/'raw';run=root/'Run.D';run.mkdir(parents=True);(run/'sample.xml').write_text('<Sample><Name>First</Name></Sample>')
            store=Path(tmp)/'saved'
            _,_,hit=saved_library(root,store);self.assertFalse(hit)
            with patch('hplc_autoquant.run_library.scan_library',side_effect=AssertionError('Unexpected recursive scan')):
                rows,_,hit=saved_library(root,store);self.assertTrue(hit);self.assertEqual(rows[0]['sample'],'First')
            with patch('hplc_autoquant.run_library.index_injection',side_effect=AssertionError('Unchanged metadata reread')):
                saved_library(root,store,refresh=True)
            (run/'sample.xml').write_text('<Sample><Name>Changed</Name></Sample>')
            rows,_,_=saved_library(root,store,refresh=True);self.assertEqual(rows[0]['sample'],'Changed')

    def test_reintegration_uses_saved_arrays_and_changes_window(self):
        from dataclasses import asdict
        with tempfile.TemporaryDirectory() as tmp:
            base=Path(tmp)/'base';(base/'data').mkdir(parents=True)
            x=np.linspace(0,4,4001);y=10+100*np.exp(-.5*((x-1)/.04)**2)+80*np.exp(-.5*((x-3)/.04)**2)
            np.savetxt(base/'data/uv.csv',np.column_stack([x,y,np.ones(len(x))*10]),delimiter=',',header='time,signal,baseline',comments='')
            np.savetxt(base/'data/tic.csv',np.column_stack([x,y*1.5,np.ones(len(x))*15]),delimiter=',',header='time,tic,baseline',comments='')
            np.savez_compressed(base/'data/ms.npz',time_min=x,mz=np.array([100.,200.]),intensity=np.column_stack([y,y*.5]))
            summary={'uv_channels':[{'channel':'UV','trace_csv':'data/uv.csv','peaks':[]}],
                     'ms_channels':[{'channel':'MSD1','trace_csv':'data/tic.csv','matrix_export':'data/ms.npz','peaks':[{'old_integrated_peak':True}]}],
                     'settings':asdict(RunSettings()),'review_notes':[]}
            (base/'run_summary.json').write_text(json.dumps(summary))
            with patch('hplc_autoquant.run_library.process_injection',side_effect=AssertionError('Unexpected raw decode')),patch('hplc_autoquant.single_run.detect_integrate',wraps=detect_integrate) as integrate:
                folder,result=reintegrate_saved(base,Path(tmp)/'new',RunSettings(integration_start_min=2.9,integration_end_min=3.1))
                self.assertEqual(integrate.call_count,1)  # UV only, despite a TIC containing the same peaks.
            peaks=result['uv_channels'][0]['peaks'];self.assertEqual(len(peaks),1);self.assertEqual(peaks[0]['area_pct_all_detected'],100.)
            self.assertAlmostEqual(peaks[0]['rt_min'],3.);self.assertTrue((folder/'manifest.json').exists())
            self.assertEqual(result['ms_channels'][0]['peaks'],[])
            self.assertFalse(result['ms_channels'][0]['integration_performed'])
            tic=np.loadtxt(folder/'data/tic.csv',delimiter=',',skiprows=1)
            np.testing.assert_allclose(tic[:,1],y*1.5)
            self.assertEqual(result['mass_associations'][0]['base_peak_mz'],100.)

    def test_range_changes_integral_and_denominator(self):
        x=np.linspace(0,4,4001);y=10+100*np.exp(-.5*((x-1)/.04)**2)+80*np.exp(-.5*((x-3)/.04)**2)
        _,full=detect_integrate(x,y,5)
        _,restricted=detect_integrate(x,y,5,RunSettings(integration_start_min=2.9,integration_end_min=3.1))
        self.assertEqual(len(full),2);self.assertEqual(len(restricted),1)
        self.assertAlmostEqual(restricted[0]['rt_min'],3)
        self.assertEqual(restricted[0]['area_pct_all_detected'],100.)
        self.assertGreaterEqual(restricted[0]['start_min'],2.9);self.assertLessEqual(restricted[0]['end_min'],3.1)
        _,clipped=detect_integrate(x,y,5,RunSettings(integration_start_min=3.,integration_end_min=3.05))
        self.assertLess(clipped[0]['area_response_seconds'],restricted[0]['area_response_seconds'])
        self.assertIn('INTEGRATION_RANGE_EDGE_REVIEW',clipped[0]['flags'])
        for start,end in [(2,1),(4,5),(float('nan'),3)]:
            with self.assertRaises(ValueError):detect_integrate(x,y,5,RunSettings(integration_start_min=start,integration_end_min=end))

    def test_live_cutoff_retains_export_rows_and_percentages(self):
        peaks=[{'peak_id':'P1','area_response_seconds':99.,'area_pct_all_detected':99.,'estimated_snr':20.,'report_included':True},
               {'peak_id':'P2','area_response_seconds':1.,'area_pct_all_detected':1.,'estimated_snr':None,'report_included':True}]
        raw={'settings':{},'review_notes':[],'uv_channels':[{'channel':'UV','peaks':peaks}],'ms_channels':[],
             'mass_associations':[{'uv_channel':'UV','uv_peak_id':'P2'}]}
        shown=display_summary(raw,2)
        self.assertEqual(len(shown['uv_channels'][0]['peaks']),2)
        self.assertEqual(shown['uv_channels'][0]['peaks'][0]['area_pct_all_detected'],99.)
        self.assertEqual(shown['uv_channels'][0]['peaks'][0]['area_pct_retained'],100.)
        self.assertIsNone(shown['uv_channels'][0]['peaks'][1]['area_pct_retained'])
        self.assertFalse(shown['mass_associations'][0]['uv_report_included'])
        self.assertTrue(raw['uv_channels'][0]['peaks'][1]['report_included'])
        restored=display_summary(raw,0)
        self.assertTrue(all(p['report_included'] for p in restored['uv_channels'][0]['peaks']))
        raw['uv_channels'][0]['peaks'][0]['area_response_seconds']=990.
        raw['uv_channels'][0]['peaks'][1]['area_response_seconds']=10.
        absolute=display_summary(raw,min_uv_area=20.)
        self.assertFalse(absolute['uv_channels'][0]['peaks'][1]['report_included'])
        self.assertEqual(absolute['uv_channels'][0]['peaks'][0]['area_pct_all_detected'],99.)

    def test_remaining_percentages_sum_to_100_and_all_excluded_is_defined(self):
        peaks=[{'peak_id':str(i),'area_response_seconds':area,'area_pct_all_detected':area,'estimated_snr':20.} for i,area in enumerate([60.,36.,4.])]
        raw={'settings':{},'review_notes':[],'uv_channels':[{'channel':'UV','peaks':peaks}],'ms_channels':[],'mass_associations':[]}
        result=display_summary(raw,5)['uv_channels'][0]
        self.assertEqual([p['area_pct_retained'] for p in result['peaks']],[62.5,37.5,None])
        self.assertEqual(result['area_reporting']['retained_area'],96.)
        self.assertEqual(result['area_reporting']['excluded_pct_original'],4.)
        self.assertEqual([p['area_response_seconds'] for p in result['peaks']],[60.,36.,4.])
        restored=display_summary(raw,0)['uv_channels'][0]
        self.assertEqual([p['area_pct_retained'] for p in restored['peaks']],[60.,36.,4.])
        empty=display_summary(raw,100)['uv_channels'][0]
        self.assertEqual(empty['area_reporting']['retained_area'],0.)
        self.assertTrue(all(p['area_pct_retained'] is None for p in empty['peaks']))
        self.assertTrue(all('area_pct_retained' not in p for p in peaks))


if __name__=='__main__':unittest.main()
