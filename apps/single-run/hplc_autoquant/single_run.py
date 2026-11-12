"""Single-injection UV/MS report, with no assay or automatic compound identity."""
from dataclasses import asdict, dataclass
from datetime import datetime
from importlib.metadata import version
import csv
import gc
import io
import json
import math
from pathlib import Path
import re
import shutil
import tempfile
import uuid
import xml.etree.ElementTree as ET
import zipfile

import numpy as np
from scipy import sparse
from scipy.integrate import trapezoid
from scipy.sparse.linalg import spsolve
from scipy.signal import find_peaks, peak_widths
from rainbow.agilent.chemstation import parse_ch, parse_ms, parse_uv

from .engine import ValidationError
from .workflow import digest, json_write, utc_now


@dataclass(frozen=True)
class RunSettings:
    uv_min_prominence_mau: float = 5.0
    ms_relative_prominence: float = .02
    baseline_lambda: float = 1e6
    baseline_reference_step_min: float = 1/150
    baseline_asymmetry: float = .01
    baseline_iterations: int = 20
    boundary_prominence_fraction: float = .99
    ms_bin_width_da: float = .05
    ms_minus_uv_delay_min: float = 0.0
    spectra_plotted_per_channel: int = 6
    filter_report: bool = False
    report_min_area_pct: float = .1
    report_min_estimated_snr: float = 3.0
    report_min_uv_area: float = 0.0
    report_min_ms_area: float = 0.0
    noise_side_window_min: float = .5
    noise_window_start_min: float = None
    noise_window_end_min: float = None
    width_review_above_min: float = 1.0
    integration_start_min: float = None
    integration_end_min: float = None


def baseline_asls(y, lam=1e6, asymmetry=.01, iterations=20):
    y=np.asarray(y,dtype=float)
    if y.ndim!=1 or len(y)<5 or not np.all(np.isfinite(y)):
        raise ValidationError('Baseline requires a finite 1D trace with at least five points')
    if not math.isfinite(lam) or lam<=0 or not 0<asymmetry<1 or not isinstance(iterations,int) or not 1<=iterations<=100:
        raise ValidationError('Invalid baseline settings')
    d=sparse.diags([1.,-2.,1.],[0,1,2],shape=(len(y)-2,len(y)),format='csc')
    penalty=lam*(d.T@d)
    w=np.ones(len(y))
    for _ in range(iterations):
        z=spsolve(sparse.diags(w,format='csc')+penalty,w*y)
        w=np.where(y>z,asymmetry,1-asymmetry)
    return z


def detect_integrate(x,y,min_prominence,settings=RunSettings(),baseline_override=None):
    x,y=np.asarray(x,float),np.asarray(y,float)
    if x.ndim!=1 or y.shape!=x.shape or len(x)<5 or not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)) or np.any(np.diff(x)<=0):
        raise ValidationError('Trace requires matching finite arrays and strictly increasing time')
    if not math.isfinite(min_prominence) or min_prominence<=0 or not 0<settings.boundary_prominence_fraction<1:
        raise ValidationError('Invalid peak threshold or boundary setting')
    # A second-difference penalty scales with dt**4. Keep its physical time
    # scale comparable between UV and MS acquisition rates.
    step=float(np.median(np.diff(x)))
    if settings.baseline_reference_step_min<=0 or not math.isfinite(settings.baseline_reference_step_min):
        raise ValidationError('Baseline reference sampling step must be positive')
    effective_lambda=settings.baseline_lambda*(settings.baseline_reference_step_min/step)**4
    if baseline_override is None:z=baseline_asls(y,effective_lambda,settings.baseline_asymmetry,settings.baseline_iterations)
    else:
        z=np.asarray(baseline_override,dtype=float)
        if z.shape!=y.shape or not np.all(np.isfinite(z)):raise ValidationError('Supplied baseline must match the finite UV trace')
    corrected=y-z
    lower=x[0] if settings.integration_start_min is None else settings.integration_start_min
    upper=x[-1] if settings.integration_end_min is None else settings.integration_end_min
    if not math.isfinite(lower) or not math.isfinite(upper) or lower>=upper or upper<=x[0] or lower>=x[-1]:
        raise ValidationError('Integration range must overlap the acquired trace and start before it ends')
    lower,upper=max(float(x[0]),lower),min(float(x[-1]),upper)
    peaks,props=find_peaks(corrected,prominence=min_prominence)
    noise_exclusions=[]
    if len(peaks) and (settings.integration_start_min is not None or settings.integration_end_min is not None):
        all_widths=peak_widths(corrected,peaks,rel_height=settings.boundary_prominence_fraction)
        noise_exclusions=list(zip(np.interp(all_widths[2],np.arange(len(x)),x),np.interp(all_widths[3],np.arange(len(x)),x)))
    selected=(x[peaks]>=lower)&(x[peaks]<=upper)
    peaks=peaks[selected];props={key:value[selected] for key,value in props.items()}
    if not len(peaks):return z,peak_metrics(x,corrected,[],settings,noise_exclusions)
    widths=peak_widths(corrected,peaks,rel_height=settings.boundary_prominence_fraction)
    left=np.interp(widths[2],np.arange(len(x)),x);right=np.interp(widths[3],np.arange(len(x)),x)
    flags=[['PROVISIONAL_INTEGRATION','IDENTITY_NOT_ASSIGNED'] for _ in peaks]
    for j in range(len(peaks)):
        if left[j]<lower or right[j]>upper:flags[j].append('INTEGRATION_RANGE_EDGE_REVIEW')
    left=np.maximum(left,lower);right=np.minimum(right,upper)
    for j in range(len(peaks)-1):
        if right[j]>left[j+1]:
            valley=int(peaks[j]+np.argmin(corrected[peaks[j]:peaks[j+1]+1]))
            right[j]=left[j+1]=x[valley]
            flags[j].append('SHARED_VALLEY_REVIEW');flags[j+1].append('SHARED_VALLEY_REVIEW')
    output=[]
    for j,k in enumerate(peaks):
        lo,hi=float(left[j]),float(right[j])
        xx=np.r_[lo,x[(x>lo)&(x<hi)],hi]
        area=float(trapezoid(np.interp(xx,x,corrected),xx)*60)
        if lo<=x[0]+np.median(np.diff(x)) or hi>=x[-1]-np.median(np.diff(x)):flags[j].append('RUN_EDGE_REVIEW')
        output.append({'peak_id':f'P{j+1:03d}','rt_min':float(x[k]),'start_min':lo,'end_min':hi,
                       'area_response_seconds':area,'height_above_baseline':float(corrected[k]),
                       'prominence':float(props['prominences'][j]),'flags':flags[j]})
    return z,peak_metrics(x,corrected,output,settings,noise_exclusions)


def peak_metrics(x,corrected,peaks,settings=RunSettings(),noise_exclusions=()):
    """Percentages are per channel, before display filtering; S/N is a screening estimate."""
    import copy
    x,corrected=np.asarray(x,float),np.asarray(corrected,float)
    if x.ndim!=1 or corrected.shape!=x.shape or len(x)<5 or np.any(np.diff(x)<=0) or not np.all(np.isfinite(x)) or not np.all(np.isfinite(corrected)):
        raise ValidationError('Invalid trace for peak metrics')
    for value in [settings.report_min_area_pct,settings.report_min_estimated_snr,settings.noise_side_window_min,settings.width_review_above_min]:
        if not math.isfinite(value) or value<0:raise ValidationError('Report thresholds must be finite and nonnegative')
    if settings.report_min_area_pct>100 or settings.noise_side_window_min<=0 or settings.width_review_above_min<=0:
        raise ValidationError('Invalid reporting/noise-window setting')
    if not isinstance(settings.filter_report,bool):raise ValidationError('Report filtering must be explicitly enabled or disabled')
    manual=settings.noise_window_start_min is not None or settings.noise_window_end_min is not None
    if manual:
        lo,hi=settings.noise_window_start_min,settings.noise_window_end_min
        if lo is None or hi is None or not math.isfinite(lo) or not math.isfinite(hi) or not x[0]<=lo<hi<=x[-1]:
            raise ValidationError('Noise window must lie inside the acquired trace')
    result=copy.deepcopy(peaks)
    positive_area=sum(p['area_response_seconds'] for p in result if p['area_response_seconds']>0)
    free=np.ones(len(x),dtype=bool)
    guard=float(np.median(np.diff(x)))*2
    for p in result:free &= ~((x>=p['start_min']-guard)&(x<=p['end_min']+guard))
    for lo,hi in noise_exclusions:free &= ~((x>=lo-guard)&(x<=hi+guard))
    delta=np.diff(corrected)
    def noise(mask):
        pairs=mask[:-1]&mask[1:]
        differences=delta[pairs]
        if len(differences)<20:return None,len(differences)
        sigma=float(np.median(np.abs(differences-np.median(differences)))/.6744897501960817/np.sqrt(2))
        return (sigma if sigma>np.finfo(float).eps*max(1,float(np.max(np.abs(corrected)))) else None),len(differences)
    global_sigma,global_pairs=noise(free)
    for p in result:
        p['flags']=list(p.get('flags',[]))
        area=p['area_response_seconds']
        p['area_pct_all_detected']=100*area/positive_area if area>0 and positive_area>0 else None
        p['area_denominator_response_seconds']=positive_area
        if area<=0:p['flags'].append('NONPOSITIVE_AREA')
        if p['end_min']-p['start_min']>settings.width_review_above_min:p['flags'].append('BROAD_FEATURE_REVIEW')
        if manual:
            window=(x>=settings.noise_window_start_min)&(x<=settings.noise_window_end_min)
            sigma,pairs=noise(window&free)
            source='user_noise_region_peak_free'
            region=[settings.noise_window_start_min,settings.noise_window_end_min]
        else:
            low=max(float(x[0]),p['start_min']-settings.noise_side_window_min)
            high=min(float(x[-1]),p['end_min']+settings.noise_side_window_min)
            sigma,pairs=noise(free&(x>=low)&(x<=high));source='local_peak_free';region=[low,high]
            if sigma is None and global_sigma is not None:
                sigma,pairs=global_sigma,global_pairs;source='global_peak_free_fallback';region=[float(x[0]),float(x[-1])]
                p['flags'].append('GLOBAL_NOISE_FALLBACK')
        p['noise_sigma_estimate']=sigma;p['noise_source']=source;p['noise_difference_pairs']=pairs;p['noise_region_min']=region
        p['estimated_snr']=p['height_above_baseline']/sigma if sigma is not None and p['height_above_baseline']>0 else None
        snr=p['estimated_snr']
        if snr is None:p['flags'].append('NOISE_NOT_ESTIMABLE')
        p['snr_screening_band']='not estimated' if snr is None else 'below 3 - low signal' if snr<3 else '3 to <10 - review' if snr<10 else '10 or above - screening only'
        reasons=[]
        if p['area_pct_all_detected'] is not None and p['area_pct_all_detected']<settings.report_min_area_pct:reasons.append('BELOW_REPORT_AREA_PERCENT')
        if snr is not None and snr<settings.report_min_estimated_snr:reasons.append('BELOW_REPORT_SNR')
        p['report_filter_reasons']=reasons
        uncertain=snr is None or p['area_pct_all_detected'] is None
        p['report_included']=not settings.filter_report or not reasons or uncertain
        if settings.filter_report and reasons and uncertain:p['flags'].append('FILTER_METRICS_INCOMPLETE_REVIEW')
    return result


def spectrum_at(time,mz,matrix,requested_rt):
    time,mz,matrix=np.asarray(time),np.asarray(mz),np.asarray(matrix)
    if time.ndim!=1 or mz.ndim!=1 or matrix.shape!=(len(time),len(mz)) or not len(time) or not len(mz):
        raise ValidationError('Invalid MS array dimensions')
    if not math.isfinite(requested_rt) or not np.all(np.isfinite(time)) or not np.all(np.isfinite(mz)) or np.any(np.diff(time)<=0) or np.any(np.diff(mz)<=0):
        raise ValidationError('MS time/mass axes must be finite and increasing')
    index=int(np.argmin(np.abs(time-requested_rt)))
    delta=float(time[index]-requested_rt)
    tolerance=max(float(np.median(np.diff(time)))*1.5,.025) if len(time)>1 else .025
    if abs(delta)>tolerance:
        return {'status':'NO_NEARBY_SCAN','requested_rt_min':requested_rt,'nearest_scan_rt_min':float(time[index]),'delta_min':delta},None
    intensities=matrix[index].astype(float,copy=True)
    top=np.argsort(intensities)[::-1][:10]
    top=[int(i) for i in top if intensities[i]>0]
    return {'status':'APEX_SCAN_ONLY','requested_rt_min':requested_rt,'scan_rt_min':float(time[index]),
            'delta_min':delta,'base_peak_mz':float(mz[top[0]]) if top else None,
            'top_mz_bins':[{'mz':float(mz[i]),'intensity':float(intensities[i])} for i in top],
            'interpretation':'No background subtraction. Nearby scan association, not confirmed compound identity or neutral molecular mass.'},intensities


def write_table(path,rows):
    fields=list(dict.fromkeys(k for row in rows for k in row))
    with Path(path).open('w',newline='',encoding='utf-8') as handle:
        if not fields:handle.write('No peaks met the recorded detection threshold\n');return
        writer=csv.DictWriter(handle,fieldnames=fields);writer.writeheader()
        for row in rows:
            safe={}
            for k,v in row.items():
                if isinstance(v,(list,dict)):v=json.dumps(v)
                if isinstance(v,str) and v.lstrip().startswith(('=','+','-','@')):v="'"+v
                safe[k]=v
            writer.writerow(safe)


def save_injection_upload(files,destination):
    """Select exactly one injection; names and case collisions are checked for Windows."""
    from pathlib import PurePosixPath
    destination=Path(destination);selected=[];roots=set();seen=set()
    allowed={'.ch','.ms','.uv'}
    for name,raw in files:
        p=PurePosixPath(str(name).replace('\\','/'))
        if p.is_absolute() or any(v=='..' or ':' in v for v in p.parts):raise ValidationError('Unsafe upload path')
        if p.suffix.lower() not in allowed and p.name.lower() not in {'sample.xml','acq.txt'}:continue
        if len(p.parts)<2 or not p.parts[-2].lower().endswith('.d'):raise ValidationError('Select one .D folder')
        roots.add('/'.join(p.parts[:-1]))
        if p.name.casefold() in seen:raise ValidationError('Duplicate file in selected injection')
        seen.add(p.name.casefold());selected.append((p.name,raw))
    if len(roots)!=1 or not any(Path(n).suffix.lower() in allowed for n,_ in selected):raise ValidationError('Select exactly one .D injection with CH, UV or MS data')
    target=destination/Path(next(iter(roots))).name;target.mkdir(parents=True,exist_ok=False)
    for name,raw in selected:(target/name).write_bytes(raw)
    return target


def process_injection(source,out,settings=RunSettings(),generate_pdf=True):
    source,out=Path(source).resolve(),Path(out).resolve()
    if source.suffix.lower()!='.d' or not source.is_dir():raise ValidationError('Choose one Agilent .D injection folder')
    if out.exists() or source==out or source in out.parents:raise ValidationError('Choose a new results directory outside the raw injection')
    if not 0<settings.ms_relative_prominence<1 or not .05<=settings.ms_bin_width_da<=1 or not math.isfinite(settings.ms_minus_uv_delay_min):raise ValidationError('Invalid MS settings')
    if not isinstance(settings.spectra_plotted_per_channel,int) or not 1<=settings.spectra_plotted_per_channel<=30:raise ValidationError('Choose 1-30 spectrum figures per channel')
    out.parent.mkdir(parents=True,exist_ok=True)
    stage=Path(tempfile.mkdtemp(prefix='.single-run-',dir=out.parent))
    for d in ['source_snapshots','data','spectra','figures']:(stage/d).mkdir()
    summary={'status':'PROVISIONAL - SINGLE INJECTION REVIEW','created_at_utc':utc_now(),
             'injection':source.name,'source_path':str(source),'settings':asdict(settings),
             'reader_versions':{p:version(p) for p in ['rainbow-api','numpy','scipy']},
             'sources':[],'uv_channels':[],'ms_channels':[],'spectral_cubes':[],'errors':[],
             'sample_metadata':{},'mass_associations':[],
             'scope':'All available supported CH/UV/MS files are attempted. Peaks are threshold-dependent; m/z is measured ion mass-to-charge, not neutral mass or confirmed identity.'}
    try:
        files=sorted(p for p in source.iterdir() if p.is_file())
        for p in files:
            if p.suffix.lower() not in ('.ch','.ms','.uv') and p.name.lower() not in ('sample.xml','acq.txt'):continue
            raw=p.read_bytes();(stage/'source_snapshots'/p.name).write_bytes(raw)
            summary['sources'].append({'file':p.name,'source_path':str(p),'sha256':digest(raw),'bytes':len(raw)})
        xml=next((p for p in (stage/'source_snapshots').iterdir() if p.name.lower()=='sample.xml'),None)
        if xml:
            try:summary['sample_metadata']={n.tag:n.text for n in ET.fromstring(xml.read_bytes()) if not list(n)}
            except ET.ParseError as exc:summary['errors'].append({'file':xml.name,'error':str(exc)})
        # UV channels are processed first so each UV peak can be paired to an MS apex scan.
        for p in sorted((stage/'source_snapshots').iterdir()):
            if p.suffix.lower()!='.ch':continue
            try:
                trace=parse_ch(str(p))
                if trace is None or trace.detector!='UV' or trace.metadata.get('unit')!='mAU':raise ValidationError('Not a supported UV/mAU channel')
                x,y=trace.xlabels,trace.data[:,0]
                z,peaks=detect_integrate(x,y,settings.uv_min_prominence_mau,settings)
                name=p.stem
                np.savetxt(stage/'data'/f'{name}_trace.csv',np.column_stack([x,y,z,y-z]),delimiter=',',header='time_min,intensity_mau,baseline_mau,corrected_mau',comments='')
                write_table(stage/'data'/f'{name}_peaks.csv',peaks)
                summary['uv_channels'].append({'file':p.name,'channel':name,'metadata':trace.metadata,
                    'area_unit':'mAU*s','points':len(x),'trace_csv':f'data/{name}_trace.csv','peaks':peaks})
            except Exception as exc:summary['errors'].append({'file':p.name,'error':str(exc)})
        for p in sorted((stage/'source_snapshots').iterdir()):
            if p.suffix.lower()!='.uv':continue
            try:
                cube=parse_uv(str(p))
                if cube is None:raise ValidationError('Unsupported UV spectrum cube')
                np.savez_compressed(stage/'data'/f'{p.stem}_uv_cube.npz',time_min=cube.xlabels,wavelength_nm=cube.ylabels,intensity=cube.data)
                summary['spectral_cubes'].append({'file':p.name,'time_points':len(cube.xlabels),'wavelength_points':len(cube.ylabels),'metadata':cube.metadata,'export':f'data/{p.stem}_uv_cube.npz'})
                del cube;gc.collect()
            except Exception as exc:summary['errors'].append({'file':p.name,'error':str(exc)})
        for p in sorted((stage/'source_snapshots').iterdir()):
            if p.suffix.lower()!='.ms':continue
            try:
                ms=parse_ms(str(p),display_precision=2,bin_width=settings.ms_bin_width_da)
                if ms is None:raise ValidationError('Unsupported MS file')
                x,mz,matrix=ms.xlabels,ms.ylabels,ms.data
                if len(x)<5 or matrix.shape!=(len(x),len(mz)) or not np.all(np.isfinite(matrix)) or np.any(matrix<0) or not np.all(np.isfinite(mz)) or np.any(np.diff(mz)<=0):raise ValidationError('Invalid MS array')
                tic=np.sum(matrix,axis=1,dtype=float)
                z=np.zeros_like(tic);peaks=[]  # MS TIC is displayed only, never integrated.
                name=p.stem
                np.savetxt(stage/'data'/f'{name}_tic.csv',np.column_stack([x,tic,z,tic-z]),delimiter=',',header='time_min,tic,baseline,corrected',comments='')
                np.savez_compressed(stage/'data'/f'{name}_ms_matrix.npz',time_min=x,mz=mz,intensity=matrix)
                channel={'file':p.name,'channel':name,'metadata':ms.metadata,'scans':len(x),'mz_bins':len(mz),
                         'mass_axis_note':f'Reader bin width {settings.ms_bin_width_da:g} Da. Display precision is not mass accuracy.',
                         'intensity_unit':'instrument response units','area_unit':'instrument response units*s',
                         'trace_csv':f'data/{name}_tic.csv','matrix_export':f'data/{name}_ms_matrix.npz',
                         'integration_performed':False,'peaks':peaks}
                for uv in summary['uv_channels']:
                    for peak in uv['peaks']:
                        rt=peak['rt_min']+settings.ms_minus_uv_delay_min
                        info,spectrum=spectrum_at(x,mz,matrix,rt)
                        assoc={'uv_channel':uv['channel'],'uv_peak_id':peak['peak_id'],'uv_rt_min':peak['rt_min'],
                               'uv_area_pct_all_detected':peak['area_pct_all_detected'],'uv_estimated_snr':peak['estimated_snr'],
                               'uv_report_included':peak['report_included'],
                               'ms_channel':name,'ms_minus_uv_delay_min':settings.ms_minus_uv_delay_min,**info}
                        if spectrum is not None:
                            rel=f'spectra/{uv["channel"]}_{peak["peak_id"]}_{name}_spectrum.csv'
                            np.savetxt(stage/rel,np.column_stack([mz,spectrum]),delimiter=',',header='mz,intensity',comments='')
                            assoc['spectrum_csv']=rel
                        summary['mass_associations'].append(assoc)
                summary['ms_channels'].append(channel)
                del ms,matrix;gc.collect()
            except Exception as exc:summary['errors'].append({'file':p.name,'error':str(exc)})
        if not summary['uv_channels'] and not summary['ms_channels']:raise ValidationError('No supported UV or MS trace could be read')
        summary['integration_scope']='UV only; raw MS total-ion chromatograms are not integrated.'
        summary['review_notes']=['MS total-ion chromatograms are displayed as raw traces without baseline correction, peak picking or area integration.',
            'UV/MS associations use nearest acquisition scan after the recorded time offset; coelution and detector timing require review.',
            'Mass spectra are not background subtracted. Dominant m/z may be matrix, solvent, adduct or fragment; no compound or neutral molecular mass is assigned.',
            'Only peaks exceeding recorded thresholds are tabulated. Run-edge peaks and unresolved shoulders may be missed or require manual review.',
            'UV area percentages, if computed externally, are response fractions and do not establish mass purity.',
            'All supported numerical arrays are retained in CSV/NPZ exports. The PDF illustrates only the most prominent peaks at the recorded figure limit.']
        summary['review_notes'].extend([
            'Area % = 100 x positive peak area / sum of all positive detected peak areas in the same channel, before report filtering. Channels/wavelengths are never pooled. It is not mass purity or formulation assay.',
            'Estimated S/N = baseline-corrected peak height / robust noise sigma. Sigma is MAD of adjacent peak-free residual differences / 0.67448975 / sqrt(2); this is a weakly-correlated-noise screening estimate, not an established detection or quantitation limit.',
            'Report filtering changes display only. Every detected peak, percentage, noise estimate and suppression reason remains in CSV/JSON. Unknown noise or nonpositive area is flagged and kept for review rather than silently suppressed.',
            'Broad or run-edge features require baseline/integration review even when estimated S/N is high. A small genuine impurity is not excluded on identity grounds.'])
        if settings.integration_start_min is not None or settings.integration_end_min is not None:
            summary['review_notes'].append(f'Integration window: {settings.integration_start_min} to {settings.integration_end_min} min, intersected with each trace. Only peaks with apex inside this window are integrated; crossing boundaries are clipped and flagged. Area percentages use positive detected areas in this window. Full raw traces remain exported. Baseline is fitted to the full trace.')
        write_table(stage/'mass_associations.csv',summary['mass_associations'])
        json_write(stage/'run_summary.json',summary)
        if generate_pdf:
            from .single_run_report import create_pdf
            create_pdf(summary,stage)
        for s in summary['sources']:
            if digest(Path(s['source_path']).read_bytes())!=s['sha256']:raise ValidationError('Source changed during processing; rerun on a stable copy')
        manifest={p.relative_to(stage).as_posix():digest(p.read_bytes()) for p in sorted(stage.rglob('*')) if p.is_file()}
        json_write(stage/'manifest.json',{'algorithm':'SHA-256','files':manifest})
        stage.rename(out)
        return out,summary
    except Exception:
        shutil.rmtree(stage,ignore_errors=True);raise


def report_zip(folder):
    data=io.BytesIO();folder=Path(folder)
    with zipfile.ZipFile(data,'w',zipfile.ZIP_DEFLATED) as z:
        for p in sorted(folder.rglob('*')):
            if p.is_file():z.write(p,p.relative_to(folder))
    return data.getvalue()
