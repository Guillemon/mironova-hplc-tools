"""User peak labels are independent of integration and transient highlighting."""
import copy
import hashlib
import json


def annotation_run_key(summary):
    signals=sorted((s['file'],s['sha256']) for s in summary.get('sources',[]) if s['file'].lower().endswith(('.ch','.ms','.uv')))
    return hashlib.sha256(json.dumps([summary.get('source_path'),signals],sort_keys=True).encode()).hexdigest()


def annotation_peak_key(channel,peak):
    return f"{channel}:{peak['rt_min']:.10f}"


def apply_peak_names(summary,names):
    result=copy.deepcopy(summary);lookup={}
    for channel in result['uv_channels']:
        for peak in channel['peaks']:
            name=names.get(annotation_peak_key(channel['channel'],peak),'')
            peak['peak_name']=name.strip() if isinstance(name,str) and name.strip() else peak['peak_id']
            lookup[(channel['channel'],peak['peak_id'])]=peak['peak_name']
    for assoc in result['mass_associations']:
        assoc['uv_peak_name']=lookup.get((assoc['uv_channel'],assoc['uv_peak_id']),assoc['uv_peak_id'])
    result['peak_name_revision']=hashlib.sha256(json.dumps(names,sort_keys=True).encode()).hexdigest()
    return result


def clean_peak_name(value):
    value=' '.join(str(value).split())
    if len(value)>80:raise ValueError('Use a peak name of 80 characters or fewer')
    return value


def chart_peak_id(event):
    value=event.get('selection',{}).get('peak_pick',[])
    if isinstance(value,list):return next((v.get('peak_key',v.get('Peak ID')) for v in value if isinstance(v,dict) and (v.get('peak_key') or v.get('Peak ID'))),None)
    if isinstance(value,dict):
        if value.get('peak_key') or value.get('Peak ID'):return value.get('peak_key',value.get('Peak ID'))
        return next((v.get('peak_key',v.get('Peak ID')) for v in value.get('vlPoint',{}).get('or',[]) if v.get('peak_key') or v.get('Peak ID')),None)
    return None


# Trusted static code only. It reports an outside click; no user text becomes JS.
OUTSIDE_CLICK_JS="""
export default function(component) {
  const {data, setTriggerValue} = component;
  if (!data.active) return;
  const handler = (event) => {
    const path = event.composedPath();
    const keep = path.some(el => el && el.matches && (
      el.matches('[role="dialog"], [data-testid="stDataFrame"], [data-testid="stVegaLiteChart"], [data-testid="stVegaLiteChartV4"]') ||
      [...el.classList].some(name => name.startsWith('st-key-peak-name-editor'))
    ));
    if (!keep) setTriggerValue('clear', Date.now());
  };
  document.addEventListener('click', handler);
  return () => document.removeEventListener('click', handler);
}
"""
