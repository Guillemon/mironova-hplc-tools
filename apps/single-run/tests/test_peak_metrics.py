import importlib.util
from dataclasses import replace
import unittest
import numpy as np


@unittest.skipUnless(importlib.util.find_spec('scipy') and importlib.util.find_spec('rainbow'),'Install single-run dependencies')
class PeakMetricsTests(unittest.TestCase):
    def fixture(self):
        x=np.linspace(0,10,10001)
        noise=np.random.default_rng(832).normal(0,.5,len(x))
        peaks=[dict(peak_id='P1',rt_min=3.,start_min=2.9,end_min=3.1,area_response_seconds=90.,height_above_baseline=20.,flags=[]),
               dict(peak_id='P2',rt_min=7.,start_min=6.9,end_min=7.1,area_response_seconds=10.,height_above_baseline=1.,flags=[])]
        return x,noise,peaks

    def test_percentage_and_independent_white_noise_estimate(self):
        from hplc_autoquant.single_run import peak_metrics
        x,y,peaks=self.fixture();result=peak_metrics(x,y,peaks)
        self.assertAlmostEqual(result[0]['area_pct_all_detected'],90.)
        self.assertAlmostEqual(result[1]['area_pct_all_detected'],10.)
        self.assertAlmostEqual(result[0]['noise_sigma_estimate'],.5,delta=.07)
        self.assertGreater(result[0]['estimated_snr'],30)
        self.assertLess(result[1]['estimated_snr'],3)
        self.assertEqual(peaks[0]['flags'],[])  # Input data are not mutated.

    def test_filter_does_not_renormalize_or_remove_export_rows(self):
        from hplc_autoquant.single_run import peak_metrics,RunSettings
        x,y,peaks=self.fixture()
        result=peak_metrics(x,y,peaks,RunSettings(filter_report=True,report_min_area_pct=15,report_min_estimated_snr=0))
        self.assertEqual(len(result),2)
        self.assertTrue(result[0]['report_included']);self.assertFalse(result[1]['report_included'])
        self.assertEqual(result[0]['area_pct_all_detected'],90.)
        self.assertIn('BELOW_REPORT_AREA_PERCENT',result[1]['report_filter_reasons'])

    def test_filter_off_preserves_low_snr_peak(self):
        from hplc_autoquant.single_run import peak_metrics,RunSettings
        x,y,peaks=self.fixture();result=peak_metrics(x,y,peaks,RunSettings(filter_report=False))
        self.assertTrue(all(p['report_included'] for p in result))
        self.assertIn('BELOW_REPORT_SNR',result[1]['report_filter_reasons'])

    def test_unknown_noise_is_flagged_and_kept(self):
        from hplc_autoquant.single_run import peak_metrics,RunSettings
        x,y,peaks=self.fixture();result=peak_metrics(x,np.zeros_like(y),peaks,RunSettings(filter_report=True,report_min_area_pct=15))
        self.assertIsNone(result[1]['estimated_snr'])
        self.assertTrue(result[1]['report_included'])
        self.assertIn('NOISE_NOT_ESTIMABLE',result[1]['flags'])

    def test_nonpositive_area_not_used_as_percentage(self):
        from hplc_autoquant.single_run import peak_metrics
        x,y,peaks=self.fixture();peaks[1]['area_response_seconds']=-10.
        result=peak_metrics(x,y,peaks)
        self.assertEqual(result[0]['area_pct_all_detected'],100)
        self.assertIsNone(result[1]['area_pct_all_detected'])
        self.assertIn('NONPOSITIVE_AREA',result[1]['flags'])

    def test_manual_noise_region_and_broad_feature_flag(self):
        from hplc_autoquant.single_run import peak_metrics,RunSettings
        x,y,peaks=self.fixture();peaks[0]['start_min']=1.;peaks[0]['end_min']=5.
        settings=RunSettings(noise_window_start_min=8.,noise_window_end_min=9.)
        result=peak_metrics(x,y,peaks,settings)
        self.assertEqual(result[0]['noise_source'],'user_noise_region_peak_free')
        self.assertIn('BROAD_FEATURE_REVIEW',result[0]['flags'])
        with self.assertRaises(ValueError):peak_metrics(x,y,peaks,replace(settings,noise_window_end_min=20.))

    def test_invalid_noise_window_rejected_even_without_peaks(self):
        from hplc_autoquant.single_run import detect_integrate,RunSettings
        with self.assertRaises(ValueError):
            detect_integrate(np.linspace(0,1,200),np.ones(200),1,RunSettings(noise_window_start_min=2,noise_window_end_min=3))


if __name__=='__main__':unittest.main()
