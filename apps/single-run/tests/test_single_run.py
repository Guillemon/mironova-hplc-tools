import importlib.util
from pathlib import Path
import tempfile
import unittest
import numpy as np

AVAILABLE=all(importlib.util.find_spec(n) for n in ['scipy','rainbow'])


@unittest.skipUnless(AVAILABLE,'Install desktop extras')
class SingleRunTests(unittest.TestCase):
    def test_baseline_and_gaussian_area(self):
        from hplc_autoquant.single_run import detect_integrate,RunSettings
        x=np.linspace(0,4,4001);y=10+3*x+100*np.exp(-.5*((x-2)/.04)**2)
        z,peaks=detect_integrate(x,y,5,RunSettings())
        self.assertEqual(len(peaks),1)
        self.assertAlmostEqual(peaks[0]['rt_min'],2,places=3)
        expected=100*.04*np.sqrt(2*np.pi)*60
        self.assertLess(abs(peaks[0]['area_response_seconds']/expected-1),.05)
        np.testing.assert_allclose(y,10+3*x+100*np.exp(-.5*((x-2)/.04)**2))

    def test_two_peaks_have_distinct_boundaries(self):
        from hplc_autoquant.single_run import detect_integrate
        x=np.linspace(0,4,4001);y=10+100*np.exp(-.5*((x-1.5)/.04)**2)+80*np.exp(-.5*((x-2.5)/.04)**2)
        _,peaks=detect_integrate(x,y,5)
        self.assertEqual(len(peaks),2)
        self.assertLessEqual(peaks[0]['end_min'],peaks[1]['start_min'])

    def test_constant_signal_no_false_zero_peak(self):
        from hplc_autoquant.single_run import detect_integrate
        _,peaks=detect_integrate(np.linspace(0,4,100),np.ones(100)*5,1)
        self.assertEqual(peaks,[])

    def test_scan_association_preserves_mass_and_timing(self):
        from hplc_autoquant.single_run import spectrum_at
        t=np.array([1.,1.01,1.02]);mz=np.array([100.05,200.10,301.15]);m=np.array([[1,2,0],[1,20,3],[0,4,1]])
        info,s=spectrum_at(t,mz,m,1.011)
        self.assertEqual(info['base_peak_mz'],200.10)
        self.assertAlmostEqual(info['delta_min'],-.001)
        np.testing.assert_array_equal(s,m[1])
        info,s=spectrum_at(t,mz,m,20)
        self.assertEqual(info['status'],'NO_NEARBY_SCAN');self.assertIsNone(s)

    def test_invalid_trace_rejected(self):
        from hplc_autoquant.single_run import detect_integrate
        with self.assertRaises(ValueError):detect_integrate([0,1,1,2,3],[1,2,3,4,5],1)

    def test_upload_accepts_only_one_injection(self):
        from hplc_autoquant.single_run import save_injection_upload
        with tempfile.TemporaryDirectory() as tmp:
            p=save_injection_upload([('A.D/MSD1.MS',b'raw'),('A.D/DAD1B.ch',b'uv')],Path(tmp))
            self.assertEqual(p.name,'A.D');self.assertEqual((p/'MSD1.MS').read_bytes(),b'raw')
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(ValueError):save_injection_upload([('A.D/MSD1.MS',b'a'),('B.D/MSD1.MS',b'b')],Path(tmp))


@unittest.skipUnless(importlib.util.find_spec('streamlit'),'Install desktop extras')
class SingleRunInterfaceTests(unittest.TestCase):
    def test_simple_initial_screen(self):
        from streamlit.testing.v1 import AppTest
        app=Path(__file__).resolve().parents[1]/'single_run_app.py'
        at=AppTest.from_file(str(app),default_timeout=30).run()
        self.assertFalse(at.exception)
        self.assertEqual(at.title[0].value,'Single Run - UV / Mass Spec')
        self.assertTrue(any(b.label=='Open / refresh library' for b in at.button))
        self.assertFalse(any('weight' in t.label.lower() for t in at.text_input))


if __name__=='__main__':unittest.main()
