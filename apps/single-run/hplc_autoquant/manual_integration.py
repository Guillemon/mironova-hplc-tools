"""Reversible UV integration edits on saved traces; raw signals and MS stay intact."""
from pathlib import Path
import copy
import json
import math
import shutil
import tempfile
import uuid
import numpy as np
from scipy.integrate import trapezoid
from .single_run import peak_metrics,spectrum_at,write_table
from .run_library import read_saved,write_saved
from .workflow import digest,json_write,utc_now
from .chart_navigation import box_zoom_bounds


def deletion_region(event):
    return box_zoom_bounds({'selection':{'zoom_region':event.get('selection',{}).get('delete_region',{})}})


def peaks_in_box(trace,peaks,bounds):
    if bounds is None:return []
    xmin,xmax,ymin,ymax=bounds
    trace=np.asarray(trace,float)
    return [p['peak_id'] for p in peaks if xmin<=p['rt_min']<=xmax and
            ymin<=float(np.interp(p['rt_min'],trace[:,0],trace[:,1]))<=ymax]


def manual_region(event):
    selected=event.get('selection',{}).get('manual_region',{})
    if not isinstance(selected,dict):return None
    values=selected.get('Time (min)',selected.get('x'))
    if not isinstance(values,(list,tuple)) or len(values)!=2:return None
    try:lo,hi=sorted(float(v) for v in values)
    except (TypeError,ValueError):return None
    return (lo,hi) if math.isfinite(lo) and math.isfinite(hi) and hi-lo>1e-10 else None


def integrate_region(trace,lo,hi,peak_id):
    trace=np.asarray(trace,float)
    if trace.ndim!=2 or trace.shape[1]<3 or len(trace)<5 or not np.all(np.isfinite(trace[:,:3])):raise ValueError('Manual integration requires a finite saved UV trace and baseline')
    x,y,z=trace[:,:3].T
    if np.any(np.diff(x)<=0) or not math.isfinite(lo) or not math.isfinite(hi) or not x[0]<=lo<hi<=x[-1]:raise ValueError('Select start and end times inside the recorded UV trace')
    xx=np.r_[lo,x[(x>lo)&(x<hi)],hi]
    corrected=y-z;yy=np.interp(xx,x,corrected)
    area=float(trapezoid(yy,xx)*60)
    if area<=0:raise ValueError('The selected region has no positive net area above the displayed baseline')
    apex=int(np.argmax(yy));flags=['MANUAL_INTEGRATION','PROVISIONAL_INTEGRATION','IDENTITY_NOT_ASSIGNED']
    if apex in (0,len(xx)-1):flags.append('MANUAL_APEX_AT_BOUNDARY_REVIEW')
    if hi-lo<2*float(np.median(np.diff(x))):flags.append('SHORT_MANUAL_REGION_REVIEW')
    return {'peak_id':peak_id,'rt_min':float(xx[apex]),'start_min':lo,'end_min':hi,'area_response_seconds':area,
            'height_above_baseline':float(yy[apex]),'prominence':max(0.,float(yy[apex]-min(yy[0],yy[-1]))),
            'flags':flags,'integration_method':'Manual bounds above the displayed baseline; interpolated endpoints and signed trapezoids'}


def apply_channel_actions(trace,peaks,actions,settings):
    result=copy.deepcopy(peaks);audit=[]
    for action in actions:
        if action['kind']=='add':
            lo,hi=action['start_min'],action['end_min']
            lower=settings.integration_start_min;upper=settings.integration_end_min
            if (lower is not None and lo<lower) or (upper is not None and hi>upper):raise ValueError('Manual bounds must lie inside the active integration time range')
            added=integrate_region(trace,lo,hi,action['peak_id'])
            overlapping=[p for p in result if p['start_min']<hi and p['end_min']>lo]
            result=[p for p in result if p not in overlapping]+[added]
            audit.append({'action':action,'replaced_integrals':overlapping,'added_integral':added})
        elif action['kind'] in ('delete','delete_many'):
            ids={action['peak_id']} if action['kind']=='delete' else set(action['peak_ids'])
            removed=[p for p in result if p['peak_id'] in ids]
            if len(removed)!=len(ids):raise ValueError('Some selected integrals no longer exist. Undo or reset the edits.')
            if not removed:raise ValueError('That integral no longer exists in this analysis. Undo or reset the edits.')
            result=[p for p in result if p['peak_id'] not in ids]
            audit.append({'action':action,'removed_integrals':removed})
        else:raise ValueError('Unknown integration edit')
    result.sort(key=lambda p:(p['rt_min'],p['peak_id']))
    result=peak_metrics(trace[:,0],trace[:,1]-trace[:,2],result,settings)
    for p in result:p['flags']=list(dict.fromkeys(p['flags']))
    return result,audit


def record_edit(path,action):
    record=read_saved(path,{'actions':[],'history':[]})
    record.setdefault('actions',[]);record.setdefault('history',[])
    action=dict(action,id=uuid.uuid4().hex,at_utc=utc_now())
    if action['kind']=='add':
        used=[int(a['peak_id'][1:]) for a in record['actions'] if a['kind']=='add' and a.get('peak_id','').startswith('M')]
        action['peak_id']=f'M{max(used,default=0)+1:03d}'
    record['actions'].append(action);record['history'].append(copy.deepcopy(action));write_saved(path,record)
    return action


def undo_edits(path,channel,reset=False):
    record=read_saved(path,{'actions':[],'history':[]});actions=record.get('actions',[])
    indexes=[i for i,a in enumerate(actions) if a['channel']==channel]
    if not indexes:return
    remove=set(indexes if reset else [indexes[-1]])
    record['history'].append({'kind':'reset' if reset else 'undo','channel':channel,'at_utc':utc_now(),'removed_action_ids':[actions[i]['id'] for i in sorted(remove)]})
    record['actions']=[a for i,a in enumerate(actions) if i not in remove];write_saved(path,record)


def edited_analysis(basefolder,summary,record,store,settings):
    actions=record.get('actions',[])
    if not actions and not record.get('history'):return Path(basefolder),summary
    basefolder=Path(basefolder)
    key=digest(json.dumps(['manual-UV-v1',str(basefolder),record],sort_keys=True).encode())
    out=Path(store)/'manual_analysis_cache'/key
    saved=read_saved(out/'run_summary.json');manifest=read_saved(out/'manifest.json')
    if saved and manifest.get('files') and all((out/f).is_file() for f in manifest['files']):return out,saved
    out.parent.mkdir(parents=True,exist_ok=True);stage=Path(tempfile.mkdtemp(prefix='.manual-UV-',dir=out.parent))
    try:
        shutil.copytree(basefolder,stage,dirs_exist_ok=True)
        (stage/'spectra').mkdir(exist_ok=True)
        (stage/'single_run_report.pdf').unlink(missing_ok=True)
        result=copy.deepcopy(summary);all_audit=[]
        for c in result['uv_channels']:
            selected=[a for a in actions if a['channel']==c['channel']]
            if selected:
                trace=np.loadtxt(stage/c['trace_csv'],delimiter=',',skiprows=1)
                c['peaks'],audit=apply_channel_actions(trace,c['peaks'],selected,settings)
                all_audit.extend(dict(a,channel=c['channel']) for a in audit)
            write_table(stage/'data'/f"{c['channel']}_peaks.csv",c['peaks'])
        # Rebuild UV-to-MS links from unchanged, decoded mass data. No MS integration.
        result['mass_associations']=[]
        for c in result['ms_channels']:
            with np.load(stage/c['matrix_export']) as m:time,mz,intensity=m['time_min'],m['mz'],m['intensity']
            for uv in result['uv_channels']:
                for p in uv['peaks']:
                    info,spectrum=spectrum_at(time,mz,intensity,p['rt_min']+settings.ms_minus_uv_delay_min)
                    link={'uv_channel':uv['channel'],'uv_peak_id':p['peak_id'],'uv_rt_min':p['rt_min'],'uv_area_pct_all_detected':p['area_pct_all_detected'],
                          'uv_estimated_snr':p['estimated_snr'],'uv_report_included':p['report_included'],'ms_channel':c['channel'],
                          'ms_minus_uv_delay_min':settings.ms_minus_uv_delay_min,**info}
                    if spectrum is not None:
                        rel=f"spectra/{uv['channel']}_{p['peak_id']}_{c['channel']}_spectrum.csv";link['spectrum_csv']=rel
                        np.savetxt(stage/rel,np.column_stack([mz,spectrum]),delimiter=',',header='mz,intensity',comments='')
                    result['mass_associations'].append(link)
        result['created_at_utc']=utc_now()
        result['manual_integration']={'baseline':'Displayed estimated baseline; no endpoint-to-endpoint baseline construction',
            'base_analysis':str(basefolder),'active_actions':actions,'history':record.get('history',[]),'audit':all_audit,
            'overlap_policy':'A new manual region replaces every integral overlapping its time bounds; no double counting'}
        result['review_notes'].append('Manual UV edits are applied before area cutoffs and percentage normalization. Dragged regions replace overlapping integrals. Deleted/replaced integrals and undo history are retained in integration_edits.json. MS data is unchanged.')
        json_write(stage/'original_auto_integrals.json',{c['channel']:c['peaks'] for c in summary['uv_channels']})
        json_write(stage/'integration_edits.json',result['manual_integration'])
        write_table(stage/'mass_associations.csv',result['mass_associations']);json_write(stage/'run_summary.json',result)
        json_write(stage/'manifest.json',{'algorithm':'SHA-256','files':{p.relative_to(stage).as_posix():digest(p.read_bytes()) for p in stage.rglob('*') if p.is_file() and p.name!='manifest.json'}})
        if out.exists():shutil.rmtree(out)
        stage.rename(out);return out,result
    except Exception:
        shutil.rmtree(stage,ignore_errors=True);raise
