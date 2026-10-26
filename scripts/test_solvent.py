import copy
import unittest
import numpy as np
import solvent as v


class SolventTests(unittest.TestCase):
    def setUp(self):
        self.x=np.linspace(0,10,1001);self.s=np.ones(1001,dtype=complex)
        self.ref={'mode':'reviewed_peak','kind':'solvent','solvent':'DMSO-d6',
            'observed_ppm':2.46,'target_ppm':2.5,'source_sha256':'A',
            'basis':'Synthetic reference identity review for test only.'}

    def test_aliases_and_unknowns(self):
        for alias in ['DMSO-D6','DMSO-d₆','d6 DMSO']:
            self.assertEqual(v.normalize(alias),'DMSO-d6')
        self.assertEqual(v.normalize('CDCl₃'),'CDCl3')
        self.assertIsNone(v.normalize('DMSO'))
        self.assertEqual(v.describe({}, {})['status'],'missing')
        self.assertEqual(v.describe({'solvent':'CDCl3/DMSO-d6'}, {})['status'],'mixed_or_ambiguous')

    def test_metadata_does_not_trigger_axis_shift(self):
        shift,meta=v.reference_shift(self.x,self.s,{'solvent':'DMSO-D6'}, {},'A')
        self.assertEqual(shift,0);self.assertEqual(meta['status'],'original_axis_preserved')

    def test_reviewed_reference_and_mismatch(self):
        cfg={'reference':self.ref}
        shift,_=v.reference_shift(self.x,self.s,{'solvent':'DMSO-D6'},cfg,'A')
        self.assertAlmostEqual(shift,.04)
        for params,h in [({'solvent':'CDCl3'},'A'),({'solvent':'DMSO-D6'},'B'),({},'A')]:
            with self.assertRaises(ValueError):v.reference_shift(self.x,self.s,params,cfg,h)

    def test_old_automatic_reference_rejected(self):
        with self.assertRaises(ValueError):
            v.reference_shift(self.x,self.s,{'solvent':'DMSO-D6'},
                {'reference':{'solvent':'DMSO-D6','search_ppm':[2.3,2.6],'target_ppm':2.5}},'A')

    def test_d2o_no_fixed_hod_and_explicit_standard(self):
        self.assertIsNone(v.describe({'solvent':'D2O'}, {})['suggested_residual_ppm'])
        ref=copy.deepcopy(self.ref);ref['solvent']='D2O'
        with self.assertRaises(ValueError):v.reference_shift(self.x,self.s,{'solvent':'D2O'},{'reference':ref},'A')
        ref['kind']='standard'
        shift,_=v.reference_shift(self.x,self.s,{'solvent':'D2O'},{'reference':ref},'A')
        self.assertAlmostEqual(shift,.04)

    def test_overlap_flags_do_not_modify_regions(self):
        cfg={'reference':self.ref,'baseline_regions':[[2.49,2.53],[3.9,4]],
             'main_integral_ppm':[5,5.35],'end_integral_ppm':[4.1,4.3],
             'observed_contaminant_regions':[{'label':'observed water candidate','bounds_ppm':[4.2,4.25]}]}
        before=copy.deepcopy(cfg);flags=v.region_conflicts(cfg)
        self.assertEqual(cfg,before)
        self.assertEqual({a['region_type'] for a in flags},{'baseline','end_integral'})


if __name__=='__main__':unittest.main(verbosity=2)
