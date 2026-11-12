"""Read-only raw-data catalog, persistent analysis cache, and explicit report exports."""
from dataclasses import asdict, replace
from datetime import datetime
from pathlib import Path
import copy
import hashlib
import json
import os
import re
import shutil
import uuid
import xml.etree.ElementTree as ET

from .single_run import RunSettings, process_injection, write_table
from .workflow import digest, json_write

CACHE_VERSION='uv-only-integration-1'


def decode_text(raw):
    if raw.startswith((b'\xff\xfe',b'\xfe\xff')):return raw.decode('utf-16',errors='replace')
    if b'\x00' in raw[:100]:return raw.decode('utf-16-le',errors='replace')
    return raw.decode('utf-8-sig',errors='replace')


def index_injection(path,root,include_sample_metadata=True):
    path,root=Path(path),Path(root)
    with os.scandir(path) as entries:
        files=sorted(Path(e.path) for e in entries if e.is_file(follow_symlinks=False))
    metadata={};texts=[];warnings=[]
    for p in files:
        if not include_sample_metadata or p.name.lower()!='sample.xml':continue
        try:
            with p.open('rb') as f:raw=f.read(64*1024)
            text=decode_text(raw);texts.append(text)
            if p.name.lower()=='sample.xml':
                metadata.update({n.tag:n.text or '' for n in ET.fromstring(raw) if not list(n)})
        except (OSError,ET.ParseError) as exc:warnings.append(f'{p.name}: {exc}')
    header={}  # Binary instrument files are never opened during library listing.
    recorded=header.get('date','');day=None;date_source=''
    for fmt in ['%d-%b-%y, %H:%M:%S','%d-%b-%Y, %H:%M:%S','%Y-%m-%d %H:%M:%S','%m/%d/%Y %H:%M:%S']:
        try:day=datetime.strptime(recorded,fmt).date().isoformat();date_source='Acquisition header';break
        except ValueError:pass
    if day is None:
        match=re.search(r'(20\d{2}-\d{2}-\d{2})',str(path))
        if match:
            try:day=datetime.strptime(match.group(),'%Y-%m-%d').date().isoformat();date_source='Folder name'
            except ValueError:pass
    if day is None:
        day=datetime.fromtimestamp(path.stat().st_mtime).date().isoformat();date_source='Folder modified date'
    relative=path.relative_to(root).as_posix() if path!=root else path.name
    return {'path':str(path.resolve()),'name':path.name,'relative':relative,'folder':path.parent.relative_to(root).as_posix() if path!=root else '.',
            'sample':metadata.get('Name') or header.get('sample') or path.stem,'day':day,'date_source':date_source,'recorded_date':recorded,
            'include_sample_metadata':include_sample_metadata,'metadata':metadata,'header':header,'files':[p.name for p in files],'warnings':warnings,
            'search_text':' '.join([relative,json.dumps(metadata),json.dumps(header),*texts,*[p.name for p in files]]).casefold()}


def scan_library(root,max_runs=5000,previous_rows=None,include_sample_metadata=True):
    root=Path(root).expanduser().resolve()
    if not root.is_dir():raise ValueError('Choose an existing parent data folder on this computer')
    found=[];warnings=[];visited=0
    previous={r["path"]:r for r in (previous_rows or [])}
    def onerror(exc):warnings.append(str(exc))
    for folder,dirs,_ in os.walk(root,followlinks=False,onerror=onerror):
        visited+=1
        if visited>50000:warnings.append('Stopped at 50,000 folders. Choose a smaller parent folder.');break
        path=Path(folder)
        dirs[:]=sorted(d for d in dirs if not d.startswith('.') and not (path/d).is_symlink())
        if path.suffix.lower()=='.d':
            dirs[:]=[]
            try:
                if previous_rows is None:found.append(index_injection(path,root,include_sample_metadata))
                else:
                    stamp=file_stamp(path);old=previous.get(str(path.resolve()))
                    row=old if old and old.get('file_stamp')==stamp and old.get('include_sample_metadata',True)==include_sample_metadata else index_injection(path,root,include_sample_metadata)
                    row['file_stamp']=stamp;found.append(row)
            except OSError as exc:warnings.append(f'{path}: {exc}')
            if len(found)>=max_runs:warnings.append(f'Stopped at {max_runs} injections. Choose a smaller parent folder.');break
    return sorted(found,key=lambda r:(r['day'],r['relative']),reverse=True),warnings


def filter_library(rows,query='',folder=None,day=None):
    words=query.casefold().split()
    return [r for r in rows if (folder is None or r['folder']==folder or r['folder'].startswith(folder.rstrip('/')+'/'))
            and (day is None or r['day']==day) and all(w in r['search_text'] for w in words)]


def selected_library_index(event,count):
    selection=event.get('selection',{})
    cells=selection.get('cells',[]);rows=selection.get('rows',[])
    index=cells[0][0] if cells else rows[0] if rows else None
    return index if isinstance(index,int) and 0<=index<count else None


def source_fingerprint(source):
    source=Path(source).resolve();h=hashlib.sha256();h.update(str(source).encode())
    for p in sorted(source.iterdir()):
        if not p.is_file() or p.is_symlink():continue
        if p.suffix.lower() not in ('.ch','.ms','.uv') and p.name.lower() not in ('sample.xml','acq.txt'):continue
        h.update(p.name.encode());h.update(b'\0')
        with p.open('rb') as f:
            while chunk:=f.read(1024*1024):h.update(chunk)
    return h.hexdigest()


def cached_analysis(source,store,settings=RunSettings()):
    # View cutoffs never change the analysis cache or integrated-area denominator.
    settings=replace(settings,filter_report=False,report_min_area_pct=0.,report_min_estimated_snr=0.,report_min_uv_area=0.,report_min_ms_area=0.)
    source=str(Path(source).resolve());stamp=file_stamp(source)
    pointer_path=Path(store)/'analysis_catalog.json'
    catalog=read_saved(pointer_path)
    slot=digest(json.dumps([source,asdict(settings)],sort_keys=True).encode())
    pointer=catalog.get(slot,{})
    if pointer.get('stamp')==stamp and pointer.get('cache_version')==CACHE_VERSION:
        target=Path(store)/pointer['folder']
        summary=read_saved(target/'run_summary.json')
        manifest=read_saved(target/'manifest.json')
        if summary and manifest.get('files') and all((target/r).is_file() for r in manifest['files']):return target,summary,True
    def usable(entry):
        folder=Path(store)/entry['folder'];manifest=read_saved(folder/'manifest.json')
        return (folder/'run_summary.json').is_file() and bool(manifest.get('files')) and all((folder/p).is_file() for p in manifest['files'])
    known=next((r for r in catalog.values() if r.get('source')==source and r.get('stamp')==stamp and usable(r)),None)
    signature=known['signature'] if known else source_fingerprint(source)
    key=digest(json.dumps([CACHE_VERSION,signature,asdict(settings)],sort_keys=True).encode())
    target=Path(store)/'analysis_cache'/key
    summary=read_saved(target/'run_summary.json');manifest=read_saved(target/'manifest.json')
    hit=bool(summary and manifest.get('files') and all((target/r).is_file() for r in manifest['files']))
    if not hit:
        if target.exists():target.rename(target.with_name(target.name+'-incomplete-'+uuid.uuid4().hex[:8]))
        target.parent.mkdir(parents=True,exist_ok=True)
        base_settings=replace(settings,integration_start_min=None,integration_end_min=None)
        saved_base=read_saved(Path(store)/known['folder']/'run_summary.json') if known else {}
        if known and saved_base.get('settings',{}).get('ms_bin_width_da')==settings.ms_bin_width_da:
            _,summary=reintegrate_saved(Path(store)/known['folder'],target,settings)
        elif settings!=base_settings:
            basefolder,_,_=cached_analysis(source,store,base_settings)
            _,summary=reintegrate_saved(basefolder,target,settings)
        else:_,summary=process_injection(source,target,settings,generate_pdf=False)
    current_stamp=file_stamp(source)
    if current_stamp!=stamp:
        # OneDrive hydration may update timestamps/status while the bytes stay
        # identical. Only a changed content hash invalidates this analysis.
        if source_fingerprint(source)!=signature:
            raise ValueError('Raw data contents changed while opening this run. Wait for acquisition to finish and open it again.')
        stamp=file_stamp(source)
    # Merge in any base-analysis entry created by the recursive call.
    catalog=read_saved(pointer_path)
    catalog[slot]={'source':source,'stamp':stamp,'folder':target.relative_to(store).as_posix(),'signature':signature,'cache_version':CACHE_VERSION}
    write_saved(pointer_path,catalog)
    return target,summary,hit


def display_summary(summary,min_area_pct=0.,min_snr=0.,min_uv_area=0.,min_ms_area=0.):
    import math
    if not all(math.isfinite(v) for v in [min_area_pct,min_snr,min_uv_area,min_ms_area]) or not 0<=min_area_pct<=100 or min(min_snr,min_uv_area,min_ms_area)<0:raise ValueError('Invalid display cutoff')
    result=copy.deepcopy(summary)
    result['review_notes']=[n for n in result['review_notes'] if not n.startswith(('Report filtering changes display only.','Area % = ','UV area percentages, if computed externally'))]
    result['review_notes'].append('Area cutoffs are applied once to the original percentages within each UV channel and integration window. Displayed Area % is then recalculated using only retained positive peak areas, which sum to 100% when any remain. Measured areas are unchanged. Original percentages, excluded peaks and reasons remain in CSV/JSON. Area % is detector response, not mass purity.')
    result['settings'].update(filter_report=bool(min_area_pct or min_snr or min_uv_area or min_ms_area),report_min_area_pct=min_area_pct,report_min_estimated_snr=min_snr,report_min_uv_area=min_uv_area,report_min_ms_area=min_ms_area)
    result['settings']['area_percent_basis']='Retained positive UV peak areas after cutoff'
    lookup={}
    for channel in result['uv_channels']+result['ms_channels']:
        absolute_cutoff=min_uv_area if channel in result['uv_channels'] else min_ms_area
        for p in channel['peaks']:
            pct,snr=p.get('area_pct_all_detected'),p.get('estimated_snr');reasons=[]
            if pct is not None and pct<min_area_pct:reasons.append('BELOW_REPORT_AREA_PERCENT')
            if snr is not None and snr<min_snr:reasons.append('BELOW_REPORT_SNR')
            if absolute_cutoff and p.get('area_response_seconds',0)>0 and p['area_response_seconds']<absolute_cutoff:reasons.append('BELOW_REPORT_ABSOLUTE_AREA')
            p['report_filter_reasons']=reasons
            # An unknown S/N does not block an independently known area cutoff.
            p['report_included']=not reasons
            if pct is None:p['report_included']=True
            lookup[(channel['channel'],p['peak_id'])]=p
        original_total=sum(p.get('area_response_seconds',0) for p in channel['peaks'] if p.get('area_response_seconds',0)>0)
        retained_total=sum(p.get('area_response_seconds',0) for p in channel['peaks'] if p['report_included'] and p.get('area_response_seconds',0)>0)
        channel['area_reporting']={'original_area':original_total,'retained_area':retained_total,
                                   'excluded_area':original_total-retained_total,
                                   'excluded_pct_original':100*(original_total-retained_total)/original_total if original_total>0 else None}
        for p in channel['peaks']:
            area=p.get('area_response_seconds',0)
            p['area_pct_retained']=100*area/retained_total if p['report_included'] and area>0 and retained_total>0 else None
            p['retained_area_denominator_response_seconds']=retained_total
    for a in result['mass_associations']:
        a['uv_report_included']=lookup[(a['uv_channel'],a['uv_peak_id'])]['report_included']
        a['uv_area_pct_retained']=lookup[(a['uv_channel'],a['uv_peak_id'])]['area_pct_retained']
    return result


def export_report(folder,summary,store):
    from .single_run_report import create_pdf
    target=Path(store)/'reports'/('report-'+datetime.now().strftime('%Y%m%d-%H%M%S')+'-'+uuid.uuid4().hex[:8])
    target.parent.mkdir(parents=True,exist_ok=True)
    shutil.copytree(folder,target)
    for channel in summary['uv_channels']:
        write_table(target/'data'/f"{channel['channel']}_peaks.csv",channel['peaks'])
    write_table(target/'mass_associations.csv',summary['mass_associations'])
    json_write(target/'run_summary.json',summary)
    create_pdf(summary,target)
    manifest={p.relative_to(target).as_posix():digest(p.read_bytes()) for p in sorted(target.rglob('*')) if p.is_file() and p.name!='manifest.json'}
    json_write(target/'manifest.json',{'algorithm':'SHA-256','files':manifest})
    return target


def read_saved(path,default=None):
    try:return json.loads(Path(path).read_text(encoding='utf-8'))
    except (OSError,ValueError):return {} if default is None else default


def write_saved(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_name('.'+path.name+'-'+uuid.uuid4().hex+'.tmp')
    try:
        temporary.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False),encoding='utf-8')
        os.replace(temporary,path)
    finally:
        if temporary.exists():temporary.unlink()


def file_stamp(source):
    """Only directory entries and file stats; no metadata or signal bytes are read."""
    source=Path(source)
    return [[p.name,p.stat().st_size,p.stat().st_mtime_ns,p.stat().st_ctime_ns]
            for p in sorted(source.iterdir()) if p.is_file() and not p.is_symlink()]


def saved_library(root,store,refresh=False):
    root=str(Path(root).expanduser().resolve())
    target=Path(store)/'library_indexes'/(digest(root.encode())+'.json')
    previous=read_saved(target)
    if previous.get('root')==root and not refresh:return previous['rows'],previous['warnings'],True
    rows,warnings=scan_library(root,previous_rows=previous.get('rows',[]))
    write_saved(target,{'root':root,'rows':rows,'warnings':warnings})
    return rows,warnings,False


def reintegrate_saved(basefolder,out,settings):
    """Reintegrate decoded arrays without reading raw CH/MS files or their metadata."""
    import tempfile
    import numpy as np
    from .single_run import detect_integrate,spectrum_at
    from .workflow import utc_now
    basefolder,out=Path(basefolder),Path(out)
    out.parent.mkdir(parents=True,exist_ok=True)
    stage=Path(tempfile.mkdtemp(prefix='.reintegrate-',dir=out.parent))
    try:
        shutil.copytree(basefolder,stage,dirs_exist_ok=True)
        for name in ['figures','spectra']:
            shutil.rmtree(stage/name,ignore_errors=True);(stage/name).mkdir()
        (stage/'single_run_report.pdf').unlink(missing_ok=True)
        summary=read_saved(stage/'run_summary.json')
        summary['created_at_utc']=utc_now();summary['settings']=asdict(settings)
        summary['processing']='Reintegrated from saved decoded arrays; raw files and metadata were not decoded again.'
        summary['mass_associations']=[]
        # Previous per-channel read errors are retained. A range with no overlap is
        # represented by zero peaks without losing the decoded channel.
        for channel in summary['uv_channels']+summary['ms_channels']:
            xy=np.loadtxt(stage/channel['trace_csv'],delimiter=',',skiprows=1)
            x,y=xy[:,0],xy[:,1]
            uv=channel in summary['uv_channels']
            overlap=(settings.integration_start_min is None or settings.integration_start_min<x[-1]) and (settings.integration_end_min is None or settings.integration_end_min>x[0])
            if uv and overlap:
                supplied=xy[:,2] if channel.get('blank_correction',{}).get('applied') else None
                z,peaks=detect_integrate(x,y,settings.uv_min_prominence_mau,settings,baseline_override=supplied)
            elif uv:z,peaks=xy[:,2],[]
            else:z,peaks=np.zeros_like(y),[]
            channel['peaks']=peaks
            if not uv:
                channel['integration_performed']=False
                channel.pop('peak_prominence_threshold',None)
                (stage/'data'/f"{channel['channel']}_tic_peaks.csv").unlink(missing_ok=True)
            np.savetxt(stage/channel['trace_csv'],np.column_stack([x,y,z,y-z]),delimiter=',',header='time_min,signal,baseline,corrected',comments='')
            if uv:write_table(stage/'data'/f"{channel['channel']}_peaks.csv",peaks)
        for channel in summary['ms_channels']:
            with np.load(stage/channel['matrix_export']) as matrix:
                time,mz,intensity=matrix['time_min'],matrix['mz'],matrix['intensity']
            name=channel['channel']
            for uv in summary['uv_channels']:
                for peak in uv['peaks']:
                    info,spectrum=spectrum_at(time,mz,intensity,peak['rt_min']+settings.ms_minus_uv_delay_min)
                    assoc={'uv_channel':uv['channel'],'uv_peak_id':peak['peak_id'],'uv_rt_min':peak['rt_min'],
                           'uv_area_pct_all_detected':peak['area_pct_all_detected'],'uv_estimated_snr':peak['estimated_snr'],
                           'uv_report_included':peak['report_included'],'ms_channel':name,'ms_minus_uv_delay_min':settings.ms_minus_uv_delay_min,**info}
                    if spectrum is not None:
                        rel=f'spectra/{uv["channel"]}_{peak["peak_id"]}_{name}_spectrum.csv';assoc['spectrum_csv']=rel
                        np.savetxt(stage/rel,np.column_stack([mz,spectrum]),delimiter=',',header='mz,intensity',comments='')
                    summary['mass_associations'].append(assoc)
        summary['integration_scope']='UV only; raw MS total-ion chromatograms are not integrated.'
        summary['review_notes']=[n for n in summary['review_notes'] if not n.startswith(('Integration window:','MS total-ion chromatograms'))]
        summary['review_notes'].insert(0,'MS total-ion chromatograms are raw traces without baseline correction, peak picking or area integration.')
        if settings.integration_start_min is not None:
            summary['review_notes'].append(f'Integration window: {settings.integration_start_min} to {settings.integration_end_min} min. Only peaks with apex inside the window are integrated; crossing bounds are clipped. Percentages use positive detected areas in this window. Full traces remain available; the baseline is fitted to the full trace.')
        write_table(stage/'mass_associations.csv',summary['mass_associations']);json_write(stage/'run_summary.json',summary)
        json_write(stage/'manifest.json',{'algorithm':'SHA-256','files':{p.relative_to(stage).as_posix():digest(p.read_bytes()) for p in stage.rglob('*') if p.is_file() and p.name!='manifest.json'}})
        stage.rename(out);return out,summary
    except Exception:
        shutil.rmtree(stage,ignore_errors=True);raise
