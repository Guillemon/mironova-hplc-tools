from pathlib import Path
import tempfile
import unittest
import numpy as np
from hplc_autoquant.comparison_views import comparison_traces,uv_spectrum


class ComparisonViewTests(unittest.TestCase):
    def test_overlays_match_wavelength_and_leave_signals_unchanged(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);x=np.linspace(0,4,41);trace=np.column_stack([x,x*3])
            np.savetxt(root/'trace.csv',trace,delimiter=',',header='time,response',comments='')
            uv={'channel':'DAD1B','metadata':{'signal':'DAD1B, Sig=210,4 Ref=360,100'},'trace_csv':'trace.csv'}
            channel={'channel':'DAD1C','metadata':{'signal':'DAD1C, Sig=210,4 Ref=360,100'}}
            summary={'injection':'other.D','source_path':'/other.D','uv_channels':[uv],'ms_channels':[]}
            overlay,missing=comparison_traces(channel,True,[(root,summary)])
            self.assertFalse(missing);np.testing.assert_allclose(next(iter(overlay.values())),trace)
            self.assertEqual(summary['uv_channels'][0],uv)
            self.assertEqual(comparison_traces(channel,False,[(root,summary)]),({},['other.D']))
            ambiguous={**summary,'uv_channels':[uv,dict(uv,channel='DAD1Z')]}
            self.assertEqual(comparison_traces(channel,True,[(root,ambiguous)]),({},['other.D (ambiguous UV match)']))
            wrong={'channel':'DAD1A','metadata':{'signal':'DAD1A, Sig=280,4 Ref=360,100'}}
            self.assertEqual(comparison_traces(wrong,True,[(root,summary)]),({},['other.D']))

    def test_uv_spectrum_uses_nearest_actual_scan_and_rejects_remote_rt(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);times=np.array([1.,1.01,1.02]);wave=np.array([210.,220.,280.]);intensity=np.arange(9).reshape(3,3)
            np.savez(root/'cube.npz',time_min=times,wavelength_nm=wave,intensity=intensity)
            rt,spectrum=uv_spectrum(root,{'export':'cube.npz'},1.009)
            self.assertEqual(rt,1.01);np.testing.assert_array_equal(spectrum[:,0],wave);np.testing.assert_array_equal(spectrum[:,1],intensity[1])
            self.assertEqual(uv_spectrum(root,{'export':'cube.npz'},2.),(None,None))
