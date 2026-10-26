import copy
import unittest
import numpy as np
import workflow as w


class ReviewTests(unittest.TestCase):
    def setUp(self):
        self.cfg={'main_integral_ppm':[5,5.35],'end_integral_ppm':[4.1,4.3],
                  'baseline_regions':[[3.9,4]],'output_dir':'first'}
        self.review={'source_sha256':'source_A','config_sha256':w.config_hash(self.cfg),
                     'decision':'quantify_provisionally'}
        for key in ['solvent_and_reference_assessment','sample_identity_basis','phase_assessment','baseline_assessment',
                    'main_assignment_and_bounds','end_assignment_and_bounds','overlap_assessment',
                    'boundary_sensitivity_assessment','formula_applicability']:
            self.review[key]='Synthetic test review entry; not experimental evidence.'

    def test_new_file_rejected(self):
        with self.assertRaises(ValueError):w.assert_review(self.cfg,self.review,'source_B')

    def test_changed_bounds_or_baseline_rejected(self):
        for field,value in [('end_integral_ppm',[4.2,4.4]),('baseline_regions',[[3.8,3.9]])]:
            cfg=copy.deepcopy(self.cfg);cfg[field]=value
            with self.assertRaises(ValueError):w.assert_review(cfg,self.review,'source_A')

    def test_draft_rejected_output_change_allowed(self):
        cfg=copy.deepcopy(self.cfg);cfg['output_dir']='second'
        w.assert_review(cfg,self.review,'source_A')
        self.review['decision']='needs_review'
        with self.assertRaises(ValueError):w.assert_review(cfg,self.review,'source_A')

    def test_candidates_follow_shifted_data(self):
        x=np.linspace(3.8,4.8,10001)
        def peaks(shift):
            y=np.exp(-((x-4.2-shift)/.01)**2)+.5*np.exp(-((x-4.24-shift)/.008)**2)
            return w.candidate_peaks(x,y,[3.9,4.7],.001)
        a,b=peaks(0),peaks(.16)
        self.assertEqual(len(a),2);self.assertEqual(len(b),2)
        for pa,pb in zip(a,b):self.assertAlmostEqual(pb['ppm']-pa['ppm'],.16,places=4)


if __name__=='__main__':unittest.main(verbosity=2)
