"""Plot-view navigation, independent of peak picking and integration."""
import math


def box_zoom_bounds(event):
    selected=event.get('selection',{}).get('zoom_region',{})
    if not isinstance(selected,dict):return None
    def interval(names):
        for name in names:
            values=selected.get(name)
            if isinstance(values,(list,tuple)) and len(values)==2:
                try:lo,hi=sorted(float(v) for v in values)
                except (TypeError,ValueError):continue
                if math.isfinite(lo) and math.isfinite(hi) and hi>lo:return lo,hi
        return None
    x=interval(['Time (min)','x']);y=interval(['Response','TIC response','Raw','y'])
    return (*x,*y) if x and y else None


def wider_view(bounds,full_x):
    xmin,xmax,ymin,ymax=bounds;cx=(xmin+xmax)/2;cy=(ymin+ymax)/2
    return max(full_x[0],cx-(xmax-xmin)),min(full_x[1],cx+(xmax-xmin)),cy-(ymax-ymin),cy+(ymax-ymin)
