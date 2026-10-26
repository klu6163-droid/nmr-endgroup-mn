"""Scientific regression checks: formula, integration, FFT axis, phase recovery."""
import copy
import json
import unittest
import numpy as np
import engine as q


class ProcessingTests(unittest.TestCase):
    def test_image_arithmetic(self):
        self.assertAlmostEqual(q.molecular_weight(17.08,1),7922.56)
        self.assertAlmostEqual(q.molecular_weight(34.16,2),7922.56)
        with self.assertRaises(ValueError): q.molecular_weight(1,0)

    def test_signed_integral_and_exact_edges(self):
        x=np.linspace(0,10,1001)
        self.assertAlmostEqual(q.integrate(x,2*x+3,[1.123,5.456]),
                               (5.456**2+3*5.456)-(1.123**2+3*1.123),places=10)
        self.assertAlmostEqual(q.integrate(x,-np.ones(len(x)),[2,4]),-2)

    def test_baseline_preserves_known_area_ratio(self):
        x=np.linspace(0,7,14001)
        signal=100*np.exp(-((x-5.16)/.025)**2)+6*np.exp(-((x-4.2)/.025)**2)
        b=.3*x+5
        fit,_=q.baseline(x,signal+b,[[0,1],[3,3.5],[4.5,4.7],[6,7]])
        ratio=q.integrate(x,signal+b-fit,[5,5.35])/q.integrate(x,signal+b-fit,[4.1,4.3])
        self.assertAlmostEqual(ratio,100/6,places=5)

    def test_fft_axis_and_wrong_domain_rejection(self):
        n=2048; sw=1000.;obs=100.;offset=6.;frequency=-125.
        fid=np.exp(2j*np.pi*frequency*np.arange(n)/sw)*np.exp(-np.arange(n)/300)
        d={"header":{"data_dimension_number":1,"units":[("none",1,"second")],
            "data_offset_start":[0],"data_axis_start":[0.],"data_axis_stop":[(n-1)/sw]},
            "parameters":{"x_sweep":sw,"x_freq":obs*1e6,"x_offset":offset}}
        c={"fid_drop_points":0,"line_broadening_hz":0,"zero_fill_factor":1,"reference":None}
        x,s,_=q.transform(d,fid,c)
        self.assertAlmostEqual(x[np.argmax(abs(s))],4.75,places=8)
        d["header"]["units"][0]=("none",1,"ppm")
        with self.assertRaises(ValueError): q.transform(d,fid,c)

    def test_recover_known_linear_phase(self):
        x=np.linspace(0,7,32768)
        s=sum(a*.001/(.001+1j*(x-c)) for a,c in [(1,2.5),(.3,5.445),(.4,5.759)])
        truth=np.array([-45.,-8.])
        distorted=s*np.exp(-1j*np.deg2rad(truth[0]+truth[1]*(x-3.5)))
        c={"phase_pivot_ppm":3.5,"manual_phase":None,"phase_anchors":[
            {"center_ppm":v,"half_width_ppm":.018} for v in [2.5,5.445,5.759]]}
        _,fit=q.phase_spectrum(x,distorted,c)
        self.assertAlmostEqual(fit['p0_degrees_at_pivot'],truth[0],delta=.15)
        self.assertAlmostEqual(fit['p1_degrees_per_ppm'],truth[1],delta=.15)


if __name__=='__main__': unittest.main(verbosity=2)
