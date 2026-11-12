"""Display-only chromatogram comparisons and raw DAD apex spectra."""
from pathlib import Path
import numpy as np
from .uv_blank import matching_blank


def comparison_traces(channel,uv,comparisons):
    traces={};missing=[]
    for folder,summary in comparisons:
        candidates=summary['uv_channels'] if uv else summary['ms_channels']
        try:
            match=matching_blank(channel,candidates) if uv else next((c for c in candidates if c['channel']==channel['channel']),None)
        except ValueError:
            missing.append(summary['injection']+' (ambiguous UV match)');continue
        if match is None:
            missing.append(summary['injection']);continue
        label='Overlay: '+summary['injection']
        if label in traces:label+=' · '+str(Path(summary['source_path']).parent)
        xy=np.loadtxt(Path(folder)/match['trace_csv'],delimiter=',',skiprows=1)
        traces[label]=xy[:,:2]
    return traces,missing


def uv_spectrum(folder,cube,rt):
    with np.load(Path(folder)/cube['export']) as data:
        times=data['time_min'];wavelength=data['wavelength_nm'];intensity=data['intensity']
        if times.ndim!=1 or wavelength.ndim!=1 or intensity.shape!=(len(times),len(wavelength)) or not len(times) or not len(wavelength):raise ValueError('Invalid saved UV spectrum dimensions')
        if not np.all(np.isfinite(times)) or np.any(np.diff(times)<=0) or not np.all(np.isfinite(wavelength)) or np.any(np.diff(wavelength)<=0):raise ValueError('Invalid saved UV spectrum axes')
        index=int(np.argmin(np.abs(times-rt)))
        tolerance=max(1.5*float(np.median(np.diff(times))),.025) if len(times)>1 else .025
        if abs(float(times[index])-rt)>tolerance:return None,None
        result=np.column_stack([wavelength,intensity[index]])
        if not np.all(np.isfinite(result)):raise ValueError('UV spectrum contains nonfinite values')
        return float(times[index]),result
