"""Content-bound review and assay calculations using the preserved engine."""
import json,math,hashlib,sys
from datetime import datetime,timezone
from importlib.metadata import version
from .calculation import analyze,ValidationError
from . import __version__


def clean(value):
    if hasattr(value,'item'):value=value.item()
    if isinstance(value,float) and not math.isfinite(value):return None
    if isinstance(value,dict):return {str(k):clean(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [clean(v) for v in value]
    return value


def fingerprint(value):return hashlib.sha256(json.dumps(clean(value),sort_keys=True,allow_nan=False).encode()).hexdigest()


def default_tables(imported):
    assignments=[];preparations={}
    for r in imported['records']:
        prep=r.get('preparation_record',{});pid=r.get('preparation_id') or ''
        assignments.append({'include':True,'injection_id':r['injection_id'],'kind':r['proposed_kind'],'standard_label':r.get('standard_label',''),
            'area':r['integration'].get('area_mau_s'),'nominal_concentration_mg_ml':r.get('proposed_concentration_mg_ml'),
            'standard_mass_mg':prep.get('formulation_mass_mg') if r['proposed_kind']=='standard' else None,
            'standard_volume_ml':prep.get('final_volume_ml') if r['proposed_kind']=='standard' else None,
            'preparation_id':pid,'reason':''})
        if r['proposed_kind']=='sample' and pid:
            candidate={'preparation_id':pid,'sample_id':r.get('sample_id',''),'formulation_mass_mg':prep.get('formulation_mass_mg'),
                       'final_volume_ml':prep.get('final_volume_ml'),'additional_dilution_factor':1. if all(str(r.get('sample_metadata',{}).get(k,'1') or '1') in ['1','1.0','1.00'] for k in ['Dilution','Multiplier']) else None}
            if pid in preparations:
                prior=preparations[pid]
                for k in ['sample_id','formulation_mass_mg','final_volume_ml','additional_dilution_factor']:
                    if prior[k]!=candidate[k]:prior[k]=None
            else:preparations[pid]=candidate
    return assignments,list(preparations.values())


def calculate(imported,payload):
    payload=clean(payload)
    if not payload.get('verified') or not str(payload.get('analyst','')).strip():raise ValidationError('A named analyst must verify the assignments, target identity, preparation records and purity settings')
    original={r['injection_id']:r for r in imported['records']};assignments=payload['assignments']
    ids=[r['injection_id'] for r in assignments]
    if len(ids)!=len(set(ids)) or set(ids)!=set(original):raise ValidationError('Every imported injection must appear exactly once in the assignment table')
    purity=payload['purity_fraction'] if payload.get('purity_confirmed') else 1.
    if purity is None or not 0<purity<=1:raise ValidationError('Purity fraction must be greater than 0 and at most 1')
    rows=[];exclusions=[];controls=[];review_flags=[]
    for a in assignments:
        iid=a['injection_id'];r=original[iid];area=a.get('area');reason=str(a.get('reason') or '').strip()
        if not a['include']:
            if not reason:raise ValidationError(iid+': give an exclusion reason')
            exclusions.append(a);continue
        if a['kind']=='unassigned':raise ValidationError(iid+': confirm a role or exclude with a reason')
        raw_area=r['integration'].get('area_mau_s')
        changed=(area is not None and (raw_area is None or not math.isclose(area,raw_area,rel_tol=1e-12,abs_tol=1e-9)))
        if changed and not reason:raise ValidationError(iid+': explain the reviewed area override')
        if changed:review_flags.append({'code':'AREA_OVERRIDE','subject':iid,'message':reason})
        if a['kind'] in ['solvent_blank','matrix_blank','system_suitability']:
            controls.append(a)
            if a['kind']=='system_suitability':review_flags.append({'code':'SST_NOT_EVALUATED','subject':iid,'message':'Chromatographic suitability is not automatically qualified'})
            elif area is None:review_flags.append({'code':'BLANK_NO_CANDIDATE','subject':iid,'message':'No unique candidate was detected; this does not demonstrate zero blank response'})
            else:rows.append({'injection_id':iid,'kind':'blank','area':area})
            continue
        if r['reader_status']!='EXTRACTED' and not changed:raise ValidationError(iid+': raw trace could not be read; resolve or exclude it')
        if area is None:raise ValidationError(iid+': no unique target area; review integration or enter a justified area override')
        row={'injection_id':iid,'kind':a['kind'],'area':area,'preparation_id':a.get('preparation_id')}
        if a['kind'] in ['standard','qc']:
            c=a.get('nominal_concentration_mg_ml');mass=a.get('standard_mass_mg');volume=a.get('standard_volume_ml')
            if mass is not None or volume is not None:
                if mass is None or volume is None or mass<=0 or volume<=0:raise ValidationError(iid+': standard mass and volume must both be positive')
                c=mass/volume
            row.update(nominal_concentration_mg_ml=c,purity_fraction=purity)
        rows.append(row)
    used={r['preparation_id'] for r in rows if r['kind']=='sample'}
    preps=[p for p in payload['preparations'] if p['preparation_id'] in used]
    reference=payload.get('reference_nominal_mg_ml')
    result=analyze(rows,preps,payload.get('criteria'),allow_two_levels=payload.get('allow_two_levels',False),calibration_model=payload['calibration_model'],reference_concentration_mg_ml=reference*purity if reference is not None else None)
    result['flags'].extend(review_flags)
    if not payload.get('purity_confirmed'):result['flags'].append({'code':'PURITY_UNCONFIRMED','subject':'reference','message':'Uncorrected estimate; reference purity was not supplied'})
    result['flags'].append({'code':'PROVISIONAL_INTEGRATION','subject':'method','message':'Target identity, baseline and method suitability require analyst review; no final batch release approval is issued'})
    prep_lookup={p['preparation_id']:p for p in preps}
    for r in result['injections']:
        if r['kind']=='sample':
            p=prep_lookup[r['preparation_id']];r['assay_pct_ww']=100*r['calculated_concentration_mg_ml']*p['final_volume_ml']*p['additional_dilution_factor']/p['formulation_mass_mg']
    return {'created_at_utc':datetime.now(timezone.utc).isoformat(),'runtime':{'python':sys.version,'libraries':{p:version(p) for p in ['numpy','scipy','rainbow-api','reportlab','streamlit']}},'analysis_id':fingerprint([imported['run_id'],payload]),'app_version':__version__,'input_fingerprint':fingerprint(payload),'import_run_id':imported['run_id'],
        'status':'PROVISIONAL - ANALYST REVIEW; NOT A BATCH RELEASE APPROVAL','analyst':payload['analyst'],'payload':payload,'result':result,'exclusions':exclusions,'controls':controls}
