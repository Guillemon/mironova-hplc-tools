"""Explicit UV-only blank subtraction. MS arrays are copied without modification."""
from dataclasses import asdict
from pathlib import Path
import hashlib
import json
import math
import shutil
import tempfile
import numpy as np
from rainbow.agilent.chemstation import parse_ch
from .run_library import read_saved,write_saved,reintegrate_saved
from .single_run_report import signal_wavelength
from .single_run import baseline_asls
from .workflow import digest,json_write


def matching_blank(channel,candidates):
    wavelength=signal_wavelength(channel)
    if wavelength is None:return None
    matches=[]
    for candidate in candidates:
        other=signal_wavelength(candidate)
        if other is None or abs(wavelength-other)>.01:continue
        a,b=channel.get('metadata',{}),candidate.get('metadata',{})
        if any(k in a and k in b and a[k]!=b[k] for k in ['bandwidth','reference_wavelength','reference_bandwidth']):continue
        matches.append(candidate)
    if len(matches)>1:raise ValueError(f'More than one matching blank channel at {wavelength:g} nm; choose a blank with an unambiguous matching signal')
    return matches[0] if matches else None


def subtract_uv_trace(sample_x,sample_y,blank_x,blank_y,scale=1.,shift_min=0.):
    arrays=[np.asarray(v,dtype=float) for v in [sample_x,sample_y,blank_x,blank_y]]
    x,y,bx,by=arrays
    if x.ndim!=1 or bx.ndim!=1 or x.shape!=y.shape or bx.shape!=by.shape or min(len(x),len(bx))<5 or any(not np.all(np.isfinite(v)) for v in arrays) or np.any(np.diff(x)<=0) or np.any(np.diff(bx)<=0):raise ValueError('UV blank correction needs finite, increasing time axes and matching signal arrays')
    if not math.isfinite(scale) or scale<0 or not math.isfinite(shift_min):raise ValueError('Blank scale must be nonnegative and time shift finite')
    # Positive shift moves the blank later on the sample time axis.
    shifted=bx+shift_min
    mask=(x>=shifted[0])&(x<=shifted[-1])
    if np.count_nonzero(mask)<5:raise ValueError('The sample and blank have no usable shared time range')
    xx=x[mask];background=np.interp(xx,shifted,by)*scale
    return xx,y[mask]-background,background,mask


def load_blank_uv(source,store):
    source=Path(source).expanduser().resolve()
    if source.suffix.lower()!='.d' or not source.is_dir():raise ValueError('Choose an existing .D injection as the UV blank')
    paths=sorted(p for p in source.iterdir() if p.is_file() and p.suffix.lower()=='.ch')
    if not paths:raise ValueError('The selected blank has no .CH UV channel files')
    records=[{'file':p.name,'source_path':str(p),'sha256':digest(p.read_bytes())} for p in paths]
    key=digest(json.dumps([str(source),records],sort_keys=True).encode());target=Path(store)/'uv_blank_cache'/key
    saved=read_saved(target/'blank.json')
    if saved and all((target/c['trace_csv']).is_file() for c in saved['channels']):return target,saved
    target.parent.mkdir(parents=True,exist_ok=True);stage=Path(tempfile.mkdtemp(prefix='.uv-blank-',dir=target.parent))
    try:
        (stage/'sources').mkdir();(stage/'data').mkdir();channels=[];errors=[]
        for record in records:
            raw=Path(record['source_path']).read_bytes()
            if digest(raw)!=record['sha256']:raise ValueError('The blank changed while it was being read; try again after acquisition finishes')
            snapshot=stage/'sources'/record['file'];snapshot.write_bytes(raw)
            try:
                trace=parse_ch(str(snapshot))
                if trace is None or trace.detector!='UV' or trace.metadata.get('unit')!='mAU':continue
                rel='data/'+Path(record['file']).stem+'.csv'
                np.savetxt(stage/rel,np.column_stack([trace.xlabels,trace.data[:,0]]),delimiter=',',header='time_min,absorbance_mau',comments='')
                channels.append({'channel':Path(record['file']).stem,'file':record['file'],'metadata':trace.metadata,'trace_csv':rel,'sha256':record['sha256'],'source_path':record['source_path']})
            except Exception as exc:errors.append({'file':record['file'],'error':str(exc)})
        if not channels:raise ValueError('No supported UV/mAU trace could be read from this blank')
        for record in records:
            if digest(Path(record['source_path']).read_bytes())!=record['sha256']:raise ValueError('Blank data contents changed during reading')
        saved={'source_path':str(source),'injection':source.name,'signature':key,'channels':channels,'errors':errors}
        json_write(stage/'blank.json',saved)
        if target.exists():shutil.rmtree(target)
        stage.rename(target);return target,saved
    except Exception:
        shutil.rmtree(stage,ignore_errors=True);raise


def blank_corrected_analysis(basefolder,blankfolder,blank,store,settings,scale=1.,shift_min=0.):
    if not math.isfinite(scale) or scale<0 or not math.isfinite(shift_min):raise ValueError('Invalid UV blank scale or time shift')
    basefolder,blankfolder=Path(basefolder),Path(blankfolder)
    key=digest(json.dumps(['UV-blank-v2-independent-baselines',str(basefolder),blank['signature'],scale,shift_min,asdict(settings)],sort_keys=True).encode())
    target=Path(store)/'uv_corrected_cache'/key
    saved=read_saved(target/'run_summary.json');manifest=read_saved(target/'manifest.json')
    if saved and manifest.get('files') and all((target/f).is_file() for f in manifest['files']):return target,saved,True
    target.parent.mkdir(parents=True,exist_ok=True)
    temporary=Path(tempfile.mkdtemp(prefix='.uv-correct-input-',dir=target.parent))
    try:
        prepared=temporary/'input';shutil.copytree(basefolder,prepared)
        summary=read_saved(prepared/'run_summary.json');matched=[];unmatched=[]
        if Path(summary['source_path']).resolve()==Path(blank['source_path']).resolve():raise ValueError('Select a different injection as the blank')
        (prepared/'source_snapshots'/'uv_blank').mkdir(parents=True,exist_ok=True)
        for channel in summary['uv_channels']:
            background=matching_blank(channel,blank['channels'])
            if background is None:
                channel['blank_correction']={'applied':False,'reason':'No matching wavelength / bandwidth / reference settings'}
                unmatched.append(channel['channel']);continue
            xy=np.loadtxt(prepared/channel['trace_csv'],delimiter=',',skiprows=1)
            bxy=np.loadtxt(blankfolder/background['trace_csv'],delimiter=',',skiprows=1)
            x,corrected,scaled,mask=subtract_uv_trace(xy[:,0],xy[:,1],bxy[:,0],bxy[:,1],scale,shift_min)
            blank_lambda=settings.baseline_lambda*(settings.baseline_reference_step_min/float(np.median(np.diff(bxy[:,0]))))**4
            blank_baseline=baseline_asls(bxy[:,1],blank_lambda,settings.baseline_asymmetry,settings.baseline_iterations)
            # Do not fit a new baseline through a potentially negative subtraction
            # trough. Subtract the independently estimated baselines instead.
            corrected_baseline=xy[mask,2]-scale*np.interp(x,bxy[:,0]+shift_min,blank_baseline)
            minimum=float(np.min(corrected));amplitude=max(float(np.ptp(xy[:,1])),settings.uv_min_prominence_mau)
            review=minimum < -max(settings.uv_min_prominence_mau,.05*amplitude)
            # Integration uses the requested window intersected with this
            # corrected trace. UV and MS acquisition endpoints can differ;
            # never reject the whole run or extrapolate the blank for that.
            lower=max(float(x[0]),settings.integration_start_min if settings.integration_start_min is not None else float(x[0]))
            upper=min(float(x[-1]),settings.integration_end_min if settings.integration_end_min is not None else float(x[-1]))
            original='data/'+channel['channel']+'_original_sample_trace.csv'
            shutil.copyfile(prepared/channel['trace_csv'],prepared/original)
            blank_csv='data/'+channel['channel']+'_scaled_uv_blank.csv'
            np.savetxt(prepared/blank_csv,np.column_stack([x,scaled]),delimiter=',',header='time_min,scaled_blank_mau',comments='')
            np.savetxt(prepared/channel['trace_csv'],np.column_stack([x,corrected,corrected_baseline,corrected-corrected_baseline]),delimiter=',',header='time_min,blank_corrected_uv,baseline_difference,corrected',comments='')
            channel['points']=len(x)
            channel['blank_correction']={'applied':True,'blank_injection':blank['injection'],'blank_channel':background['channel'],
                'blank_file':background['file'],'blank_sha256':background['sha256'],'scale':scale,'shift_min':shift_min,
                'coverage_min':[float(x[0]),float(x[-1])],'original_trace_csv':original,'scaled_blank_csv':blank_csv,
                'effective_integration_min':[lower,upper] if lower<upper else None,
                'baseline_method':'Difference of independently estimated sample and blank AsLS baselines',
                'minimum_residual_mau':minimum,'review_negative_residual':review,
                'interpretation':'Sample minus time-aligned scaled blank. No extrapolation; negative residuals retained.'}
            shutil.copyfile(blankfolder/'sources'/background['file'],prepared/'source_snapshots'/'uv_blank'/background['file'])
            summary['sources'].append({'file':'uv_blank/'+background['file'],'source_path':background['source_path'],'sha256':background['sha256']})
            matched.append(channel['channel'])
        if not matched:raise ValueError('No UV wavelengths with matching bandwidth/reference settings are shared by this sample and blank')
        summary['blank_correction']={'method':'Time-aligned UV absorbance subtraction','blank_injection':blank['injection'],'blank_path':blank['source_path'],
            'scale':scale,'shift_min':shift_min,'matched_channels':matched,'uncorrected_channels':unmatched,'ms_modified':False}
        summary['review_notes']=[n.replace('Spectra are not background subtracted.','Mass spectra are not background subtracted.') for n in summary['review_notes']]
        summary['review_notes'].insert(0,f"UV blank correction: sample minus {scale:g} x {blank['injection']}, shifted {shift_min:+g} min. Only shared time coverage is corrected; no extrapolation or negative clipping. MS data and m/z spectra are unchanged. Unmatched UV channels remain uncorrected: {', '.join(unmatched) or 'none'}.")
        summary['review_notes'].insert(1,'Corrected UV baseline is sample AsLS baseline minus aligned, scaled blank AsLS baseline. Negative subtraction troughs are not fitted as a new baseline. Substantial negative residuals are flagged for blank/scale/timing review; the warning is a screening diagnostic, not a validated acceptance criterion.')
        summary['review_notes'].insert(2,'The requested integration window is intersected with each corrected UV trace’s shared sample/blank coverage. No blank signal is extrapolated. Channels with no time overlap have no integrated peaks.')
        json_write(prepared/'run_summary.json',summary)
        if target.exists():shutil.rmtree(target)
        output,result=reintegrate_saved(prepared,target,settings)
        return output,result,False
    finally:shutil.rmtree(temporary,ignore_errors=True)


def uv_display_traces(folder,channel,blankfolder=None,blank=None):
    """Independent plot choices; the displayed blank is its recorded response."""
    folder=Path(folder)
    xy=np.loadtxt(folder/channel['trace_csv'],delimiter=',',skiprows=1)
    correction=channel.get('blank_correction',{})
    if not correction.get('applied'):
        return {'Original sample':xy[:,:2],'Baseline':xy[:,[0,2]]}
    original=np.loadtxt(folder/correction['original_trace_csv'],delimiter=',',skiprows=1)
    traces={'Original sample':original[:,:2]}
    if blankfolder is not None and blank:
        matched=next((c for c in blank['channels'] if c['channel']==correction['blank_channel']),None)
        if matched:
            background=np.loadtxt(Path(blankfolder)/matched['trace_csv'],delimiter=',',skiprows=1)
            # Align time as used for correction, but never scale the displayed response.
            background=background[:,:2].copy();background[:,0]+=correction.get('shift_min',0.)
            traces['Blank']=background
    traces['Baseline']=xy[:,[0,2]]
    traces['Blank-corrected sample']=xy[:,:2]
    return traces
