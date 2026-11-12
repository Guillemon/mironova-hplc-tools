import unittest
from hplc_autoquant.peak_interaction import annotation_peak_key,apply_peak_names,clean_peak_name,chart_peak_id


class PeakInteractionTests(unittest.TestCase):
    def test_names_follow_rt_after_peak_ids_are_renumbered(self):
        original={'peak_id':'P002','rt_min':3.,'area_response_seconds':12.3}
        names={annotation_peak_key('DAD1B',original):'My peak'}
        raw={'uv_channels':[{'channel':'DAD1B','peaks':[dict(original,peak_id='P001')]}],
             'mass_associations':[{'uv_channel':'DAD1B','uv_peak_id':'P001'}]}
        labelled=apply_peak_names(raw,names)
        peak=labelled['uv_channels'][0]['peaks'][0]
        self.assertEqual(peak['peak_name'],'My peak')
        self.assertEqual(peak['peak_id'],'P001');self.assertEqual(peak['area_response_seconds'],12.3)
        self.assertEqual(labelled['mass_associations'][0]['uv_peak_name'],'My peak')
        self.assertNotIn('peak_name',raw['uv_channels'][0]['peaks'][0])
        self.assertNotEqual(labelled['peak_name_revision'],apply_peak_names(raw,{})['peak_name_revision'])

    def test_chart_peak_and_empty_selection(self):
        self.assertEqual(chart_peak_id({'selection':{'peak_pick':[{'Peak ID':'P003'}]}}),'P003')
        self.assertEqual(chart_peak_id({'selection':{'peak_pick':{'vlPoint':{'or':[{'Peak ID':'P002'}]}}}}),'P002')
        self.assertIsNone(chart_peak_id({'selection':{'peak_pick':[]}}))
        self.assertIsNone(chart_peak_id({}))

    def test_names_are_plain_labels_with_bounded_length(self):
        self.assertEqual(clean_peak_name('  CEA\n candidate '),'CEA candidate')
        self.assertEqual(clean_peak_name(''),'')
        with self.assertRaises(ValueError):clean_peak_name('a'*81)


if __name__=='__main__':unittest.main()
