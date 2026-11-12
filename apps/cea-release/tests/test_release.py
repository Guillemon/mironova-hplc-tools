import unittest,tempfile,json
from pathlib import Path
import numpy as np
from cea_release.raw import assignments,integrate_target,IntegrationSettings,import_sequence
from cea_release.calculation import analyze,ValidationError
from cea_release.review import calculate,default_tables,fingerprint


class ReleaseTests(unittest.TestCase):
    def fixture(self):
        records=[]
        for name,kind,area,c,pid in [('s1a','standard',1000.,1.,''),('s1b','standard',1000.,1.,''),('s2','standard',9999.,2.,''),('sample','sample',500.,None,'P1')]:
            records.append({'injection_id':name,'proposed_kind':kind,'proposed_concentration_mg_ml':c,'preparation_id':pid,'sample_id':'batch','integration':{'area_mau_s':area},'reader_status':'EXTRACTED'})
        imported={'run_id':'test','records':records};rows,preps=default_tables(imported)
        preps[0].update(formulation_mass_mg=100.,final_volume_ml=25.)
        payload={'assignments':rows,'preparations':preps,'calibration_model':'single_standard','reference_nominal_mg_ml':1.,'purity_fraction':1.,'purity_confirmed':False,'analyst':'Test','verified':True}
        return imported,payload

    def test_name_comments_and_blank_classification(self):
        standard=assignments('CEA Standard #1','12.36mg in 25mL\nCEA reference','STD.D')
        self.assertEqual(standard['proposed_kind'],'standard');self.assertAlmostEqual(standard['proposed_concentration_mg_ml'],.4944)
        sample=assignments('CEA_F_1026-01 Sample 1','','misleading Sample 11.D')
        self.assertEqual(sample['preparation_id'],'Sample-1')
        self.assertEqual(assignments('CEA Sample #1','163 mg in 25 mL phenoxyethanol','')['proposed_kind'],'sample')
        self.assertEqual(assignments('CEA Standard #1','12 mg in 25 mL phenoxyethanol','')['proposed_kind'],'standard')
        self.assertEqual(assignments('Phenoxyethanol Blank','','')['proposed_kind'],'matrix_blank')
        self.assertEqual(assignments('MeOH_H2O','','')['proposed_kind'],'solvent_blank')
        self.assertEqual(assignments('Unknown','','')['proposed_kind'],'unassigned')
        self.assertEqual(assignments('CEA 0.4368 mg per mL','','')['proposed_concentration_mg_ml'],.4368)

    def test_exact_integration_units_and_ambiguous_peak(self):
        x=np.linspace(0,4,4001);y=10+100*np.exp(-.5*((x-2)/.04)**2)
        settings=IntegrationSettings(rt_min=1,rt_max=3,min_prominence_mau=10,half_window_min=.3)
        result=integrate_target(x,y,settings);expected=100*.04*np.sqrt(2*np.pi)*60
        self.assertAlmostEqual(result['area_mau_s'],expected,delta=.01)
        self.assertEqual(integrate_target(x,y+100*np.exp(-.5*((x-2.7)/.04)**2),settings)['flags'],['MULTIPLE_TARGET_CANDIDATES'])
        self.assertIsNone(integrate_target(x,np.zeros_like(x),settings)['area_mau_s'])
        with self.assertRaises(ValidationError):integrate_target(x[::-1],y,settings)

    def test_single_reference_checks_do_not_change_assay_and_purity_once(self):
        imported,payload=self.fixture();result=calculate(imported,payload)
        self.assertEqual(result['result']['preparations'][0]['assay_pct_ww'],12.5)
        self.assertIsNone(result['result']['r_squared'])
        payload.update(purity_confirmed=True,purity_fraction=.98)
        result=calculate(imported,payload)
        self.assertAlmostEqual(result['result']['preparations'][0]['assay_pct_ww'],12.25)
        self.assertNotIn('PURITY_UNCONFIRMED',[r['code'] for r in result['result']['flags']])

    def test_standard_mass_volume_takes_priority(self):
        imported,payload=self.fixture()
        for r in payload['assignments']:
            if r['kind']=='standard' and r['nominal_concentration_mg_ml']==1.:
                r.update(standard_mass_mg=20.,standard_volume_ml=25.)
        payload['reference_nominal_mg_ml']=.8
        self.assertAlmostEqual(calculate(imported,payload)['result']['preparations'][0]['assay_pct_ww'],10.)

    def test_exclusion_override_and_missing_preparation_guards(self):
        imported,payload=self.fixture();payload['assignments'][2]['include']=False
        with self.assertRaisesRegex(ValidationError,'exclusion reason'):calculate(imported,payload)
        payload['assignments'][2]['reason']='Check injection excluded for validation'
        calculate(imported,payload)
        payload['assignments'][3]['area']=600.
        with self.assertRaisesRegex(ValidationError,'override'):calculate(imported,payload)
        payload['assignments'][3]['reason']='Reviewed integration override';calculate(imported,payload)
        payload['preparations'][0]['formulation_mass_mg']=None
        with self.assertRaises(ValidationError):calculate(imported,payload)

    def test_undetected_blank_is_not_reported_as_zero(self):
        imported,payload=self.fixture()
        imported['records'].append({'injection_id':'blank','reader_status':'EXTRACTED','integration':{'area_mau_s':None}})
        payload['assignments'].append({'injection_id':'blank','include':True,'kind':'solvent_blank','area':None,'reason':''})
        result=calculate(imported,payload)
        self.assertTrue(any(f['code']=='BLANK_NO_CANDIDATE' for f in result['result']['flags']))
        self.assertFalse(any(r['injection_id']=='blank' for r in result['result']['injections']))

    def test_import_cache_tracks_raw_contents_and_repairs_changed_snapshot(self):
        from unittest.mock import patch
        from types import SimpleNamespace
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=root/'sequence';d=source/'Sample 1.D';d.mkdir(parents=True)
            (d/'DAD1B.ch').write_bytes(b'raw trace')
            (d/'SAMPLE.XML').write_text('<Sample><Name>CEA Sample 1</Name></Sample>')
            x=np.linspace(9,13,4001);y=20+500*np.exp(-.5*((x-11)/.04)**2)
            trace=SimpleNamespace(detector='UV',metadata={'unit':'mAU','wavelength':208.},xlabels=x,data=y[:,None])
            with patch('cea_release.raw.parse_ch',return_value=trace):
                folder,first=import_sequence(source,root/'results')
                (folder/'source_snapshots/Sample 1.D/DAD1B.ch').write_bytes(b'corrupted snapshot')
                repaired,again=import_sequence(source,root/'results')
                self.assertEqual(repaired,folder)
                self.assertEqual((repaired/'source_snapshots/Sample 1.D/DAD1B.ch').read_bytes(),b'raw trace')
                (d/'SAMPLE.XML').write_text('<Sample><Name>CEA Sample 2</Name></Sample>')
                _,changed=import_sequence(source,root/'results')
                self.assertNotEqual(first['run_id'],changed['run_id'])
                self.assertEqual(changed['records'][0]['preparation_id'],'Sample-2')

    def test_verification_and_stale_fingerprint(self):
        imported,payload=self.fixture();prior=fingerprint(payload);payload['analyst']='Changed'
        self.assertNotEqual(fingerprint(payload),prior)
        payload['verified']=False
        with self.assertRaises(ValidationError):calculate(imported,payload)

    def test_ols_equal_weight_level_means_and_two_level_guard(self):
        rows=[{'injection_id':f's{i}', 'kind':'standard','area':1000*c+100,'nominal_concentration_mg_ml':c,'purity_fraction':1.} for i,c in enumerate([.2,.4,.6])]
        rows.append({'injection_id':'sample','kind':'sample','area':600.,'preparation_id':'P1'})
        preps=[{'preparation_id':'P1','sample_id':'batch','formulation_mass_mg':100.,'final_volume_ml':25.,'additional_dilution_factor':1.}]
        result=analyze(rows,preps);self.assertAlmostEqual(result['slope'],1000.);self.assertAlmostEqual(result['intercept'],100.);self.assertAlmostEqual(result['preparations'][0]['assay_pct_ww'],12.5)
        with self.assertRaises(ValidationError):analyze(rows[1:],preps)
        self.assertIn('TWO_LEVEL_CALIBRATION',[r['code'] for r in analyze(rows[1:],preps,allow_two_levels=True)['flags']])


if __name__=='__main__':unittest.main()
