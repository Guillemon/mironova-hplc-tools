"""Read-only Agilent sequence import and the recovered CEA integration method."""
from pathlib import Path
from dataclasses import dataclass,asdict
import csv,hashlib,json,re,os,shutil,tempfile
import xml.etree.ElementTree as ET
import numpy as np
from scipy.signal import find_peaks,peak_widths
from scipy.integrate import trapezoid
from rainbow.agilent.chemstation import parse_ch
from .calculation import ValidationError


def digest(data):return hashlib.sha256(data).hexdigest()
def save_json(path,data):Path(path).write_text(json.dumps(data,indent=2,allow_nan=False),encoding='utf-8')

def decode(data):
    if data[:2] in (b'\xff\xfe',b'\xfe\xff'):return data.decode('utf-16',errors='replace')
    if b'\x00' in data[:200]:return data.decode('utf-16-le',errors='replace')
    return data.decode('utf-8-sig',errors='replace')


@dataclass(frozen=True)
class IntegrationSettings:
    rt_min:float=10.3
    rt_max:float=11.4
    min_prominence_mau:float=100.
    half_window_min:float=.30
    baseline_guard_halfwidth_min:float=.025
    min_width_min:float=.04
    max_width_min:float=.30

    def validate(self):
        values=asdict(self)
        if not all(np.isfinite(v) for v in values.values()):raise ValidationError('Integration parameters must be finite')
        if not 0<=self.rt_min<self.rt_max or min(self.min_prominence_mau,self.half_window_min,self.baseline_guard_halfwidth_min,self.min_width_min)<=0 or self.max_width_min<self.min_width_min:raise ValidationError('Invalid integration RT, prominence, width or baseline settings')


def integrate_target(x,y,settings):
    settings.validate();x,y=np.asarray(x,float),np.asarray(y,float)
    if x.ndim!=1 or y.shape!=x.shape or len(x)<5 or not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)) or np.any(np.diff(x)<=0):raise ValidationError('Invalid finite/increasing UV trace')
    indices,properties=find_peaks(y,prominence=settings.min_prominence_mau)
    widths=peak_widths(y,indices,rel_height=.5)[0]*np.median(np.diff(x)) if len(indices) else []
    candidates=[{'index':int(i),'rt_min':float(x[i]),'prominence_mau':float(p),'width_min':float(w)} for i,p,w in zip(indices,properties['prominences'],widths) if settings.rt_min<=x[i]<=settings.rt_max and settings.min_width_min<=w<=settings.max_width_min]
    result={'status':'REVIEW_REQUIRED','flags':[],'candidates':candidates,'area_mau_s':None}
    if len(candidates)!=1:
        result['flags']=['NO_TARGET_CANDIDATE' if not candidates else 'MULTIPLE_TARGET_CANDIDATES'];return result
    rt=candidates[0]['rt_min'];lo,hi=rt-settings.half_window_min,rt+settings.half_window_min
    guard=settings.baseline_guard_halfwidth_min
    if lo-guard<x[0] or hi+guard>x[-1]:result['flags']=['TARGET_WINDOW_OUTSIDE_TRACE'];return result
    left=y[np.abs(x-lo)<=guard];right=y[np.abs(x-hi)<=guard]
    if not len(left) or not len(right):result['flags']=['NO_BASELINE_GUARD_POINTS'];return result
    a,b=float(np.median(left)),float(np.median(right))
    xx=np.r_[lo,x[(x>lo)&(x<hi)],hi];yy=np.interp(xx,x,y);baseline=np.interp(xx,[lo,hi],[a,b])
    area=float(trapezoid(yy-baseline,xx)*60)
    result.update(status='PROVISIONAL',rt_min=rt,start_min=lo,end_min=hi,baseline_start_mau=a,baseline_end_mau=b,area_mau_s=area,height_mau=float(np.interp(rt,x,y)-np.interp(rt,[lo,hi],[a,b])),flags=['EXPERIMENTAL_INTEGRATION','CEA_IDENTITY_REQUIRES_REVIEW','COELUTION_NOT_RULED_OUT'])
    if area<=0:result['status']='REVIEW_REQUIRED';result['flags'].append('NONPOSITIVE_TARGET_AREA')
    return result


def assignments(name,description,filename):
    text=(name or filename)+' '+description
    lower=text.casefold();prep=None;sample_id='';standard_label='';concentration=None
    mass=re.search(r'(\d+(?:\.\d+)?)\s*mg\s*(?:in|/|per)\s*(\d+(?:\.\d+)?)\s*m[lL]',description,re.I)
    preparation={'formulation_mass_mg':float(mass[1]),'final_volume_ml':float(mass[2])} if mass else {}
    direct=re.search(r'(\d+(?:\.\d+)?)\s*mg\s*(?:/|per)\s*m[lL]',text,re.I)
    sample=re.search(r'\bsample\s*#?\s*(\d+)\b',name or filename,re.I)
    standard=re.search(r'\b(?:std|standard)\s*#?\s*(\d+)\b',name or filename,re.I)
    named_control=any(v in (name or filename).casefold() for v in ['blank','meoh','methanol','phenoxy','solvent'])
    if sample:
        kind='sample';prep='Sample-'+sample[1];sample_id=(name or filename)[:sample.start()].strip(' _-')
    elif standard or (direct and not named_control):
        kind='standard';concentration=float(direct[1]) if direct else (preparation.get('formulation_mass_mg',0)/preparation['final_volume_ml'] if preparation.get('final_volume_ml',0)>0 else None)
        standard_label='Standard-'+standard[1] if standard else ''
    elif 'phenoxy' in lower:kind='matrix_blank'
    elif any(s in lower for s in ['blank','meoh','methanol','solvent']):kind='solvent_blank'
    else:kind='unassigned'
    return {'proposed_kind':kind,'proposed_concentration_mg_ml':concentration,'preparation_id':prep,'sample_id':sample_id,'standard_label':standard_label,'preparation_record':preparation,'assignment_basis':'Sample name and comments; analyst verification required'}


def wavelength(metadata):
    if metadata.get('wavelength') is not None:return float(metadata['wavelength'])
    m=re.search(r'Sig\s*=\s*([\d.]+)',metadata.get('signal',''),re.I)
    return float(m[1]) if m else None


def sequence_folders(source):
    source=Path(source).expanduser().resolve()
    if not source.is_dir() or source.suffix.casefold()=='.d':raise ValidationError('Choose the sequence folder containing the .D injections')
    found=[]
    for parent,dirs,_ in os.walk(source,followlinks=False):
        dirs[:]=sorted(d for d in dirs if not d.startswith('.') and not (Path(parent)/d).is_symlink())
        if Path(parent).suffix.casefold()=='.d':found.append(Path(parent));dirs[:]=[]
        if len(found)>2000:raise ValidationError('Choose a smaller sequence folder (maximum 2,000 injections)')
    if not found:raise ValidationError('No .D injection folders were found')
    return source,sorted(found)


def import_sequence(source,store,settings=IntegrationSettings(),requested_wavelength=None):
    settings.validate();source,folders=sequence_folders(source)
    inputs=[];h=hashlib.sha256(str(source).encode())
    for d in folders:h.update(d.relative_to(source).as_posix().encode())
    for d in folders:
        for p in sorted(d.iterdir()):
            if p.is_file() and not p.is_symlink() and (p.suffix.casefold()=='.ch' or p.name.casefold() in ['sample.xml','report.txt']):
                data=p.read_bytes();rel=p.relative_to(source).as_posix();sha=digest(data);h.update(rel.encode());h.update(sha.encode());inputs.append((p,rel,sha,data))
    signature=h.hexdigest();key=digest(json.dumps(['CEA-release-raw-v1',signature,asdict(settings),requested_wavelength],sort_keys=True).encode())
    target=Path(store)/'imports'/key
    if (target/'manifest.json').is_file():
        manifest=json.loads((target/'manifest.json').read_text())
        if all((target/f).is_file() and digest((target/f).read_bytes())==sha for f,sha in manifest['files'].items()):return target,json.loads((target/'import.json').read_text(encoding='utf-8'))
    target.parent.mkdir(parents=True,exist_ok=True);stage=Path(tempfile.mkdtemp(prefix='.import-',dir=target.parent))
    try:
        (stage/'traces').mkdir();snap=stage/'source_snapshots';snap.mkdir()
        sources=[]
        for p,rel,sha,data in inputs:
            dest=snap/rel;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(data)
            sources.append({'source_path':str(p),'relative_file':rel,'sha256':sha,'bytes':len(data)})
        selected=requested_wavelength;records=[]
        for index,d in enumerate(folders):
            rel=d.relative_to(source).as_posix();local=snap/rel
            meta={};xml=next((p for p in local.iterdir() if p.name.casefold()=='sample.xml'),None) if local.exists() else None
            if xml:
                try:meta={n.tag:n.text or '' for n in ET.fromstring(xml.read_bytes()) if not list(n)}
                except ET.ParseError:meta={'metadata_error':'SAMPLE.XML could not be parsed'}
            rec={'injection_id':rel,'sample_metadata':meta,**assignments(meta.get('Name',''),meta.get('Description',''),d.stem),'sources':[s for s in sources if str(Path(s['relative_file']).parent)==rel]}
            try:
                choices=[]
                for p in sorted(local.iterdir()):
                    if p.suffix.casefold()!='.ch':continue
                    trace=parse_ch(str(p))
                    if trace is not None and trace.detector=='UV' and trace.metadata.get('unit')=='mAU' and wavelength(trace.metadata) is not None:choices.append((p,trace))
                available={wavelength(t.metadata) for _,t in choices}
                if selected is None and available:selected=208. if 208. in available else (210. if 210. in available else min(available))
                matching=[(p,t) for p,t in choices if selected is not None and abs(wavelength(t.metadata)-selected)<.01]
                if len(matching)!=1:raise ValidationError(f'Expected one UV channel at {selected} nm; found {len(matching)}')
                p,trace=matching[0];x,y=trace.xlabels,trace.data[:,0]
                rec.update(reader_status='EXTRACTED',channel_metadata=trace.metadata,channel_file=p.name,integration=integrate_target(x,y,settings),trace_npz=f'traces/{index:04d}.npz')
                np.savez_compressed(stage/rec['trace_npz'],time_min=x,intensity_mau=y)
            except Exception as exc:rec.update(reader_status='ERROR',error=str(exc),integration={'area_mau_s':None,'status':'REVIEW_REQUIRED','flags':['READ_ERROR'],'candidates':[]})
            records.append(rec)
        result={'run_id':key,'source':str(source),'source_signature':signature,'wavelength_nm':selected,'settings':asdict(settings),'records':records,'sources':sources,'integration_method':'Recovered local median-guard linear baseline; signed trapezoids, interpolated endpoints, response in mAU*s'}
        for p,_,sha,_ in inputs:
            if digest(p.read_bytes())!=sha:raise ValidationError('Raw data changed during import; wait for acquisition to finish and retry')
        if not any(r['reader_status']=='EXTRACTED' for r in records):raise ValidationError('No supported matching UV/mAU channel could be read. Check the wavelength and files.')
        save_json(stage/'import.json',result)
        save_json(stage/'manifest.json',{'algorithm':'SHA-256','files':{p.relative_to(stage).as_posix():digest(p.read_bytes()) for p in stage.rglob('*') if p.is_file()}})
        if target.exists():shutil.rmtree(target)
        stage.rename(target);return target,result
    except Exception:shutil.rmtree(stage,ignore_errors=True);raise
