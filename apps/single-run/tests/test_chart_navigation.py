import unittest
from hplc_autoquant.chart_navigation import box_zoom_bounds,wider_view

class ChartNavigationTests(unittest.TestCase):
    def test_box_updates_both_axes_and_accepts_reverse_drag(self):
        self.assertEqual(box_zoom_bounds({'selection':{'zoom_region':{'Time (min)':[12,9],'Response':[1000,0]}}}),(9,12,0,1000))
        self.assertEqual(box_zoom_bounds({'selection':{'zoom_region':{'Time (min)':[4,5],'TIC response':[10,20]}}}),(4,5,10,20))
    def test_empty_invalid_clicks_do_not_zoom(self):
        for values in [{},{'Time (min)':[1,1],'Response':[0,1]},{'Time (min)':[1,2],'Response':[0,float('inf')]}]:
            self.assertIsNone(box_zoom_bounds({'selection':{'zoom_region':values}}))
    def test_zoom_out_widens_and_caps_time_to_full_run(self):
        self.assertEqual(wider_view((9,11,0,100),(0,20)),(8,12,-50,150))
        self.assertEqual(wider_view((0,20,0,100),(0,20)),(0,20,-50,150))

if __name__=='__main__':unittest.main()
