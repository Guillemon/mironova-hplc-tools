from pathlib import Path
import tempfile
import unittest
import numpy as np
from hplc_autoquant.single_run import RunSettings
from hplc_autoquant.manual_integration import integrate_region,apply_channel_actions,record_edit,undo_edits,manual_region,deletion_region,peaks_in_box
from hplc_autoquant.run_library import read_saved


class ManualIntegrationTests(unittest.TestCase):
    def fixture(self):
        x=np.linspace(0,4,4001);z=10+2*x;y=z+100*np.exp(-.5*((x-2)/.05)**2)
        trace=np.column_stack([x,y,z]);auto=integrate_region(trace,1.8,2.2,'P001');auto['flags']=[]
        return trace,[auto]

    def test_exact_area_units_and_interpolated_endpoints(self):
        x=np.linspace(0,4,401);trace=np.column_stack([x,10+2*x+5,10+2*x])
        p=integrate_region(trace,1.003,2.237,'M001')
        self.assertAlmostEqual(p['area_response_seconds'],5*(2.237-1.003)*60)
        self.assertEqual((p['start_min'],p['end_min']),(1.003,2.237))
        self.assertIn('MANUAL_INTEGRATION',p['flags'])

    def test_new_region_replaces_overlap_without_double_counting(self):
        trace,peaks=self.fixture()
        actions=[{'kind':'add','channel':'UV','start_min':1.9,'end_min':2.1,'peak_id':'M001'}]
        result,audit=apply_channel_actions(trace,peaks,actions,RunSettings())
        self.assertEqual(len(result),1);self.assertEqual(result[0]['peak_id'],'M001')
        self.assertLess(result[0]['area_response_seconds'],peaks[0]['area_response_seconds'])
        self.assertEqual(result[0]['area_pct_all_detected'],100.)
        self.assertEqual(audit[0]['replaced_integrals'][0]['peak_id'],'P001')
        self.assertEqual(peaks[0]['peak_id'],'P001')

    def test_delete_removes_integral_and_preserves_audit(self):
        trace,peaks=self.fixture()
        result,audit=apply_channel_actions(trace,peaks,[{'kind':'delete','channel':'UV','peak_id':'P001'}],RunSettings())
        self.assertEqual(result,[])
        self.assertEqual(audit[0]['removed_integrals'][0]['area_response_seconds'],peaks[0]['area_response_seconds'])

    def test_box_delete_uses_apex_in_both_axes_and_undo_restores_group(self):
        x=np.linspace(0,4,4001);z=np.zeros_like(x)
        y=sum(h*np.exp(-.5*((x-t)/.04)**2) for t,h in [(1,30),(2,50),(3,90)])
        trace=np.column_stack([x,y,z])
        peaks=[integrate_region(trace,t-.15,t+.15,f'P{i:03d}') for i,t in enumerate([1,2,3],1)]
        bounds=deletion_region({'selection':{'delete_region':{'Time (min)':[2.5,.5],'Response':[60,10]}}})
        self.assertEqual(peaks_in_box(trace,peaks,bounds),['P001','P002'])
        self.assertEqual(peaks_in_box(trace,peaks,(.5,3.5,60,100)),['P003'])
        self.assertEqual(peaks_in_box(trace,peaks,(.5,3.5,100,120)),[])
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'edits.json'
            record_edit(path,{'kind':'delete_many','channel':'UV','peak_ids':['P001','P002'],'selection_bounds':list(bounds)})
            result,audit=apply_channel_actions(trace,peaks,read_saved(path)['actions'],RunSettings())
            self.assertEqual([p['peak_id'] for p in result],['P003'])
            self.assertEqual(result[0]['area_pct_all_detected'],100.)
            self.assertEqual(len(audit[0]['removed_integrals']),2)
            undo_edits(path,'UV')
            restored,_=apply_channel_actions(trace,peaks,read_saved(path)['actions'],RunSettings())
            self.assertEqual([p['peak_id'] for p in restored],['P001','P002','P003'])
        self.assertIsNone(deletion_region({'selection':{'delete_region':{}}}))

    def test_undo_and_reset_are_scoped_to_channel(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'edits.json'
            record_edit(path,{'kind':'add','channel':'UV1','start_min':1.,'end_min':2.})
            record_edit(path,{'kind':'add','channel':'UV2','start_min':1.,'end_min':2.})
            record_edit(path,{'kind':'delete','channel':'UV1','peak_id':'M001'})
            undo_edits(path,'UV1')
            self.assertEqual([a['kind'] for a in read_saved(path)['actions']],['add','add'])
            undo_edits(path,'UV1',True)
            self.assertEqual([a['channel'] for a in read_saved(path)['actions']],['UV2'])
            self.assertEqual(len(read_saved(path)['history']),5)

    def test_invalid_negative_and_outside_active_window_rejected(self):
        trace,peaks=self.fixture()
        with self.assertRaises(ValueError):integrate_region(trace,-1,2,'M001')
        with self.assertRaises(ValueError):integrate_region(np.column_stack([trace[:,0],trace[:,2]-1,trace[:,2]]),1,2,'M001')
        with self.assertRaises(ValueError):apply_channel_actions(trace,peaks,[{'kind':'add','channel':'UV','start_min':1.,'end_min':2.,'peak_id':'M001'}],RunSettings(integration_start_min=1.5,integration_end_min=3.))
        self.assertEqual(manual_region({'selection':{'manual_region':{'Time (min)':[2,1]}}}),(1.,2.))
        self.assertIsNone(manual_region({'selection':{'manual_region':{}}}))


if __name__=='__main__':unittest.main()
