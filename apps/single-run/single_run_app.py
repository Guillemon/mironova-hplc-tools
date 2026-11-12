"""Local Agilent run library with cached analysis, live tables and optional PDF export."""
from pathlib import Path
from dataclasses import replace
import base64
import hashlib
import json
import math
import os
import uuid

import altair as alt
import numpy as np
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components
from hplc_autoquant.single_run import RunSettings, save_injection_upload, report_zip
from hplc_autoquant.single_run_report import signal_wavelength, uv_label
from hplc_autoquant.run_library import scan_library,filter_library,cached_analysis,display_summary,export_report,selected_library_index,read_saved,write_saved

from hplc_autoquant.comparison_views import comparison_traces,uv_spectrum
from hplc_autoquant.manual_integration import manual_region,deletion_region,peaks_in_box,apply_channel_actions,record_edit,undo_edits,edited_analysis
from hplc_autoquant.chart_navigation import box_zoom_bounds,wider_view
from hplc_autoquant.uv_blank import load_blank_uv,blank_corrected_analysis,uv_display_traces
from hplc_autoquant.single_point_assay import calculate_single_point
from hplc_autoquant.peak_interaction import annotation_run_key,annotation_peak_key,apply_peak_names,clean_peak_name,chart_peak_id,OUTSIDE_CLICK_JS

st.set_page_config(page_title='Mironova Run Library',page_icon='🧪',layout='wide')
st.markdown('<style>.block-container {padding-top: 2rem;}</style>',unsafe_allow_html=True)
st.title('Single Run - UV / Mass Spec')
st.caption('Browse your local data, open an injection, and review chromatograms and tables. Export a PDF only when you need one.')
store=Path(os.environ.get('AUTOQUANT_DATA',str(Path.cwd()/'single_run_results'))).expanduser().resolve()
store.mkdir(parents=True,exist_ok=True)
preferences=store/'library_settings.json'
try:remembered=json.loads(preferences.read_text())
except (OSError,ValueError):remembered={}

@st.cache_data(show_spinner=False,max_entries=128)
def array_csv(path,modified):
    return np.loadtxt(path,delimiter=',',skiprows=1)

def read_array(path):
    path=Path(path)
    return array_csv(str(path),path.stat().st_mtime_ns)

with st.sidebar:
    st.header('Data library')
    mode=st.radio('Open data',['Browse parent folder','Open one .D folder','Upload one .D folder'])
    chosen=None
    if mode=='Browse parent folder':
        with st.form('library_location'):
            root=st.text_input('Parent folder containing your .D runs',value=remembered.get('root',''),placeholder=r'C:\Users\Public\Documents\ChemStation\1\Data')
            include_metadata=st.checkbox('Read sample names and comments',value=True,help='Reads only SAMPLE.XML. Uncheck for filenames-only listing. Instrument signals and acquisition headers are never opened during listing.')
            scan=st.form_submit_button('Open / refresh library',type='primary')
        if scan:
            try:
                with st.spinner('Listing folders and sample records…'):rows,warnings=scan_library(root.strip().strip('"'),include_sample_metadata=include_metadata)
                st.session_state['library_rows']=rows;st.session_state['library_warnings']=warnings
                st.session_state['active_run']=None
                preferences.write_text(json.dumps({'root':str(Path(root.strip().strip('"')).expanduser().resolve())}))
            except Exception as exc:st.error(str(exc))
        rows=st.session_state.get('library_rows',[])
        query=st.text_input('Search filenames, sample names and comments',placeholder='Sample, lot, method, date…')
        folders=['All folders']+sorted(set(r['folder'] for r in rows))
        subfolder=st.selectbox('Sequence / subfolder',folders)
        days=['All days']+sorted({r['day'] for r in rows},reverse=True)
        day=st.selectbox('Day',days)
        visible=filter_library(rows,query,None if subfolder=='All folders' else subfolder,None if day=='All days' else day)
        st.caption(f'{len(visible)} of {len(rows)} injections · newest day first')
        if st.session_state.get('library_warnings'):
            with st.expander('Folders requiring attention'):st.write(st.session_state['library_warnings'])
        if visible:
            paths=[r['path'] for r in visible];lookup={r['path']:r for r in visible}
            if st.session_state.get('active_run') not in paths:st.session_state['active_run']=None
            def move(delta):
                index=paths.index(st.session_state['active_run']) if st.session_state['active_run'] in paths else -1;st.session_state['active_run']=paths[max(0,min(len(paths)-1,index+delta))]
            left,right=st.columns(2);position=paths.index(st.session_state['active_run']) if st.session_state['active_run'] in paths else -1
            left.button('← Previous',on_click=move,args=(-1,),disabled=position<=0)
            right.button('Next →',on_click=move,args=(1,),disabled=position==len(paths)-1)
            chosen=st.selectbox('Injection',paths,format_func=lambda p:lookup[p]['day']+' | '+lookup[p]['relative'],key='active_run',index=None,placeholder='Select an injection to view its data')
            if chosen:st.caption(lookup[chosen]['sample'])
        elif rows:st.info('No injections match these filters.')
    elif mode=='Open one .D folder':
        with st.form('one_path'):
            path=st.text_input('Path to one .D injection folder',placeholder=r'C:\ChemStation\Data\Sample 11.D')
            if st.form_submit_button('Open injection',type='primary'):st.session_state['direct_path']=path.strip().strip('"')
        chosen=st.session_state.get('direct_path')
    else:
        files=st.file_uploader('Select one .D injection folder',type=['ch','uv','ms','xml','txt'],accept_multiple_files='directory')
        if st.button('Open uploaded injection',disabled=not files):
            destination=store/('loaded-'+uuid.uuid4().hex[:12]);destination.mkdir()
            try:st.session_state['upload_path']=str(save_injection_upload([(f.name,f.getvalue()) for f in files],destination))
            except Exception as exc:st.error(str(exc))
        chosen=st.session_state.get('upload_path')
    st.caption('Original instrument files are read only. Results stay on this computer.')

if mode=='Browse parent folder' and st.session_state.get('library_rows'):
    with st.expander('Browse injections — click any cell to open',expanded=True):
        if visible:
            tablekey='catalog-cell-'+hashlib.sha256(json.dumps([[r['path'] for r in visible],chosen]).encode()).hexdigest()[:12]
            def choose_row():
                event=st.session_state.get(tablekey,{})
                index=selected_library_index(event,len(visible))
                if index is not None:st.session_state['active_run']=visible[index]['path']
            st.dataframe(pd.DataFrame([{'Day':r['day'],'Date source':r['date_source'],'Sample':r['sample'],'Injection':r['name'],'Folder':r['folder']} for r in visible]),hide_index=True,on_select=choose_row,selection_mode='single-cell',key=tablekey,width='stretch')


if not chosen:
    st.info('Open the parent folder in the sidebar, then select a .D injection or click its row. Folder listing does not load chromatograms or mass spectra until you choose a run.')
    st.stop()

try:
    with st.spinner('Opening this injection…'):
        basefolder,base,base_hit=cached_analysis(chosen,store)
except Exception as exc:
    st.error(str(exc));st.stop()

runkey=hashlib.sha256(str(chosen).encode()).hexdigest()[:12]
limits=[]
for c in base['uv_channels']+base['ms_channels']:
    xy=read_array(basefolder/c['trace_csv']);limits.append((float(xy[0,0]),float(xy[-1,0])))
start=min(v[0] for v in limits);end=max(v[1] for v in limits)
settings=RunSettings(**st.session_state.get('settings-'+runkey,{}))
st.subheader(base['sample_metadata'].get('Name') or base['injection'])
st.caption(base['injection']+' · '+str(chosen))
analysis_cols=st.columns(3)
integration_panel=analysis_cols[0].expander('Integration settings',expanded=False)
with integration_panel:
    with st.form('integration-'+runkey):
        left,right=st.columns(2)
        lo=left.number_input('Integration start (min)',value=float(settings.integration_start_min if settings.integration_start_min is not None else start),min_value=0.,format='%.3f',key='integration_start-'+runkey)
        hi=right.number_input('Integration end (min)',value=float(settings.integration_end_min if settings.integration_end_min is not None else end),min_value=0.,format='%.3f',key='integration_end-'+runkey)
        st.caption('Enter the time window, then click Reintegrate with these settings. Peaks outside it leave the tables; crossing boundaries are clipped. The raw trace stays visible, with the active window highlighted.')
        with st.expander('Advanced analysis settings'):
            c1,c2=st.columns(2)
            uv=c1.number_input('Minimum UV peak prominence (mAU)',value=settings.uv_min_prominence_mau,min_value=.01)
            delay=c2.number_input('MS time minus UV time (min)',value=settings.ms_minus_uv_delay_min,step=.01)
            figures=st.number_input('Spectrum figures per channel in PDF',value=settings.spectra_plotted_per_channel,min_value=1,max_value=30)
            manual=st.checkbox('Specify a baseline noise region',value=settings.noise_window_start_min is not None)
            c1,c2=st.columns(2)
            noise_lo=c1.number_input('Noise region start (min)',value=float(settings.noise_window_start_min if settings.noise_window_start_min is not None else 8.),min_value=0.)
            noise_hi=c2.number_input('Noise region end (min)',value=float(settings.noise_window_end_min if settings.noise_window_end_min is not None else 9.),min_value=0.)
        apply=st.form_submit_button('Reintegrate with these settings',type='primary')
        def reset_window():
            st.session_state['integration_start-'+runkey]=start
            st.session_state['integration_end-'+runkey]=end
        reset=st.form_submit_button('Reintegrate full run',on_click=reset_window)

if apply or reset:
    full=reset or (abs(lo-start)<1e-8 and abs(hi-end)<1e-8)
    if not full and not start<=lo<hi<=end:st.error('Choose a start and end inside the acquired run.');st.stop()
    settings=replace(settings,integration_start_min=None if full else lo,integration_end_min=None if full else hi,
                     uv_min_prominence_mau=uv,ms_minus_uv_delay_min=delay,spectra_plotted_per_channel=int(figures),
                     noise_window_start_min=noise_lo if manual else None,noise_window_end_min=noise_hi if manual else None)
    from dataclasses import asdict
    st.session_state['settings-'+runkey]=asdict(settings)


try:
    if settings==RunSettings():folder,raw,hit=basefolder,base,base_hit
    else:
        with st.spinner('Applying integration settings…'):folder,raw,hit=cached_analysis(chosen,store,settings)
except Exception as exc:st.error(str(exc));st.stop()
def library_picker_filters(rows,prefix):
    if not rows:
        st.caption('Load a parent folder in Data library to search its injections here, or paste a .D path below.')
        return []
    query=st.text_input('Search names and comments',key=prefix+'-query',placeholder='Sample, blank, lot, filename…')
    folders=['All folders']+sorted({r['folder'] for r in rows})
    if st.session_state.get(prefix+'-folder') not in folders:st.session_state[prefix+'-folder']='All folders'
    sequence=st.selectbox('Folder / sequence',folders,key=prefix+'-folder')
    folder_rows=filter_library(rows,query,None if sequence=='All folders' else sequence)
    days=['All days']+sorted({r['day'] for r in folder_rows},reverse=True)
    if st.session_state.get(prefix+'-day') not in days:st.session_state[prefix+'-day']='All days'
    day=st.selectbox('Day',days,key=prefix+'-day')
    visible=filter_library(folder_rows,'',None,None if day=='All days' else day)
    st.caption(f'{len(visible)} of {len(rows)} injections match. Filters update the list immediately.')
    return visible


with analysis_cols[1].expander('UV blank correction',expanded=False):
    st.caption('Select a blank measured with the same UV method. Subtraction is sample minus scaled blank at matched wavelengths and times. Mass data is not changed.')
    blank_rows=library_picker_filters(st.session_state.get('library_rows',[]),'blank-picker-'+runkey)
    with st.form('blank-'+runkey):
        blank_choice=None
        if blank_rows:
            indexed={r['path']:r for r in blank_rows}
            choices=sorted(indexed,key=lambda p:(not any(w in (indexed[p]['name']+' '+indexed[p]['sample']).casefold() for w in ['blank','meoh','phenoxy']),indexed[p]['relative']))
            blank_choice=st.selectbox('Blank injection from library',choices,index=None,placeholder='Choose the actual UV blank',format_func=lambda p:indexed[p]['relative']+' | '+indexed[p]['sample'])
        blank_path=st.text_input('Or paste a blank .D folder path',placeholder=r'C:\ChemStation\Data\Solvent blank.D')
        left,right=st.columns(2)
        blank_scale=left.number_input('Blank scale factor',value=1.,min_value=0.,step=.1,help='1.0 subtracts the recorded blank response directly. Adjust only with an analytical basis.')
        blank_shift=right.number_input('Blank time shift (min)',value=0.,step=.01,help='Positive moves the blank later. There is no automatic peak alignment.')
        apply_blank=st.form_submit_button('Apply UV blank correction')
        undo_blank=st.form_submit_button('Show uncorrected UV')
    if undo_blank:st.session_state.pop('uv-blank-'+runkey,None)
    if apply_blank:
        selected_blank=blank_path.strip().strip('"') or blank_choice
        if not selected_blank:st.error('Choose the actual blank injection first.')
        else:st.session_state['uv-blank-'+runkey]={'path':selected_blank,'scale':blank_scale,'shift':blank_shift}
    blank_config=st.session_state.get('uv-blank-'+runkey)
    if blank_config:
        try:
            with st.spinner('Reading the selected UV blank and correcting UV traces…'):
                blank_folder,blank_record=load_blank_uv(blank_config['path'],store)
                folder,raw,hit=blank_corrected_analysis(folder,blank_folder,blank_record,store,settings,blank_config['scale'],blank_config['shift'])
            st.success('UV blank correction active: '+blank_record['injection'])
            if raw['blank_correction']['uncorrected_channels']:st.warning('No matching blank for: '+', '.join(raw['blank_correction']['uncorrected_channels'])+'. Those UV channels remain uncorrected.')
            if blank_record['errors']:st.warning('Some blank files could not be read.');st.json(blank_record['errors'])
        except Exception as exc:
            st.error(str(exc));st.stop()
comparisons=[]
with analysis_cols[2].expander('Overlay chromatograms',expanded=False):
    st.caption('Compare original chromatograms at matching UV wavelengths and MS channels. Overlays keep their recorded times and response units; they do not change integration or blank correction. Overlays are for on-screen comparison; reports use the current injection.')
    overlay_rows=library_picker_filters([r for r in st.session_state.get('library_rows',[]) if r['path']!=str(chosen)],'overlay-picker-'+runkey)
    st.caption('Add selections from different folders or days without losing existing overlays. Clear overlays removes all comparisons.')
    with st.form('overlay-form-'+runkey):
        overlay_lookup={r['path']:r for r in overlay_rows if r['path']!=str(chosen)}
        overlay_choices=st.multiselect('Compare injections from library',list(overlay_lookup),format_func=lambda p:overlay_lookup[p]['relative']+' | '+overlay_lookup[p]['sample'])
        overlay_paths=st.text_area('Or paste other .D folder paths — one per line',placeholder=r'C:\ChemStation\Data\Other sample.D')
        apply_overlays=st.form_submit_button('Apply overlays')
        clear_overlays=st.form_submit_button('Clear overlays')
    if clear_overlays:st.session_state['overlays-'+runkey]=[]
    if apply_overlays:
        requested=list(dict.fromkeys(overlay_choices+[p.strip().strip('"') for p in overlay_paths.splitlines() if p.strip()]))
        st.session_state['overlays-'+runkey]=list(dict.fromkeys(st.session_state.get('overlays-'+runkey,[])+[p for p in requested if Path(p).expanduser().resolve()!=Path(chosen).expanduser().resolve()]))
    for comparison_path in st.session_state.get('overlays-'+runkey,[]):
        try:
            with st.spinner('Opening comparison injection…'):
                cf,cs,_=cached_analysis(comparison_path,store)
            comparisons.append((cf,cs));st.caption(cs['injection']+' · '+str(comparison_path))
        except Exception as exc:st.error('Could not open overlay '+comparison_path+': '+str(exc))

automatic_folder=Path(folder)
edits_path=store/'manual_integration'/(hashlib.sha256(str(automatic_folder).encode()).hexdigest()+'.json')
integration_edits=read_saved(edits_path,{'actions':[],'history':[]})
try:
    folder,raw=edited_analysis(automatic_folder,raw,integration_edits,store,settings)
except Exception as exc:
    st.error(str(exc))
    if st.button('Reset invalid integration edits'):
        write_saved(edits_path,{'actions':[],'history':integration_edits.get('history',[])+[{'kind':'reset_invalid'}]});st.rerun()
    st.stop()
if st.session_state.get('integration_edit_error'):st.error(st.session_state['integration_edit_error'])
st.session_state['last_analysis']=(str(folder),raw)

with integration_panel:
    st.subheader('Peak display filters')
    c1,c2,c3=st.columns([1,1,2])
    min_pct=c1.number_input('Minimum area % shown',value=0.,min_value=0.,max_value=100.,step=.1,format='%.3f',key='area_cutoff',help='Cutoff uses original area percentages. Remaining peaks are then renormalized to 100% within each UV channel.')
    min_snr=c2.number_input('Minimum estimated S/N shown',value=0.,min_value=0.,step=1.,key='snr_cutoff')
    show_hidden=c3.checkbox('Show hidden peaks for review',value=False,key='show_hidden')
    with st.expander('Absolute UV area cutoff (optional)'):
        min_uv_area=st.number_input('Minimum UV area shown (mAU*s)',value=0.,min_value=0.,step=1.,key='uv_area_cutoff')
min_ms_area=0.
summary=display_summary(raw,min_pct,min_snr,min_uv_area,min_ms_area)
names_path=store/'peak_names'/(annotation_run_key(summary)+'.json')
peak_names=read_saved(names_path)
summary=apply_peak_names(summary,peak_names)
pending=st.session_state.pop('pending_manual_peak',None)
if pending:
    st.session_state['peak_highlight']={'analysis':str(folder),'channel':pending['channel'],'peak_id':pending['peak_id']}
    peak=next((p for c in summary['uv_channels'] if c['channel']==pending['channel'] for p in c['peaks'] if p['peak_id']==pending['peak_id']),None)
    if peak and not peak['report_included']:st.info(f"Manual integral {peak['peak_id']} was added ({peak['area_response_seconds']:.3f} mAU*s), but is hidden by the current cutoffs. Lower the cutoff or enable Show hidden peaks for review.")
st.session_state['current_view']=summary
scope='Full run' if settings.integration_start_min is None else f'{settings.integration_start_min:g}–{settings.integration_end_min:g} min'
st.caption('Active integration: '+scope+' · UV Area % uses retained peaks · '+('Saved analysis reopened' if hit else 'Analysis saved'))
if summary['errors']:
    st.warning('Some files could not be read. See Files & metadata for details.')

# A single-point external-standard calculation intentionally lives after the
# current run has been integrated, so it uses the exact displayed/edited peak
# area and can never trigger integration of the entire library.
with st.expander('Single-point assay — compare this unknown with one known standard', expanded=False):
    st.caption('Choose one standard injection with a known concentration, then compare one UV peak from that standard with one peak in the current unknown. MS data is not used.')
    assay_rows=[r for r in st.session_state.get('library_rows',[]) if r['path'] != str(chosen)]
    assay_visible=library_picker_filters(assay_rows,'assay-standard-picker-'+runkey)
    assay_lookup={r['path']:r for r in assay_visible}
    with st.form('single-point-assay-'+runkey):
        if assay_lookup:
            assay_standard_path=st.selectbox('Known standard injection',list(assay_lookup),index=None,placeholder='Choose the standard .D injection',format_func=lambda p:assay_lookup[p]['relative']+' | '+assay_lookup[p]['sample'])
        else:
            assay_standard_path=None
        assay_pasted=st.text_input('Or paste a standard .D folder path',placeholder=r'C:\ChemStation\Data\Standard 1.D')
        a,b,c=st.columns(3)
        assay_standard_conc=a.number_input('Standard concentration',min_value=0.,value=1.,format='%.8g',help='Enter the known concentration and its units below.')
        assay_units=b.text_input('Concentration units',value='mg/mL')
        assay_dilution=c.number_input('Unknown dilution factor',min_value=0.,value=1.,format='%.8g',help='1.0 means report the measured unknown solution concentration. Use a larger value to correct for a prior dilution.')
        assay_response=st.number_input('Relative response factor (optional)',min_value=0.,value=1.,format='%.8g',help='Leave at 1.0 unless an established response factor is available.')
        load_assay=st.form_submit_button('Load standard for assay',type='primary')
    if load_assay:
        selected_assay=assay_pasted.strip().strip('"') or assay_standard_path
        if not selected_assay:st.error('Choose or paste the standard injection first.')
        elif assay_standard_conc<=0 or assay_dilution<=0 or assay_response<=0:st.error('Concentration, dilution and response factor must be positive.')
        else:
            st.session_state['single-point-assay-'+runkey]={'path':selected_assay,'standard_concentration':assay_standard_conc,'units':assay_units.strip() or 'concentration units','unknown_dilution_factor':assay_dilution,'response_factor':assay_response}
    assay_config=st.session_state.get('single-point-assay-'+runkey)
    if assay_config:
        try:
            with st.spinner('Opening the selected standard…'):
                assay_folder,assay_raw,_=cached_analysis(assay_config['path'],store,settings)
            assay_uv=sorted(assay_raw['uv_channels'],key=lambda c:(signal_wavelength(c) is None,signal_wavelength(c) or 0,c['channel']))
            current_uv={c['channel']:c for c in summary['uv_channels']}
            if not assay_uv:st.warning('The selected standard has no UV channels available for a single-point assay.')
            else:
                def assay_channel_label(c):
                    wavelength=signal_wavelength(c)
                    return (f'{wavelength:g} nm · {c["channel"]}' if wavelength is not None else c['channel'])
                assay_channel=st.selectbox('Standard UV signal',assay_uv,format_func=assay_channel_label,key='assay-standard-channel-'+runkey)
                matching=next((c for c in summary['uv_channels'] if signal_wavelength(c)==signal_wavelength(assay_channel)),None)
                if matching is None:
                    st.warning('The current unknown does not have a matching UV wavelength.')
                else:
                    standard_peaks=[p for p in assay_channel['peaks'] if p.get('report_included',True)]
                    unknown_peaks=[p for p in matching['peaks'] if p.get('report_included',True)]
                    if not standard_peaks or not unknown_peaks:
                        st.warning('Both the standard and unknown need at least one visible integrated peak at this wavelength.')
                    else:
                        def peak_label(p):
                            return f"{p.get('peak_name',p['peak_id'])} · RT {p['rt_min']:.4f} min · area {p['area_response_seconds']:.3f}"
                        left,right=st.columns(2)
                        standard_peak=left.selectbox('Standard peak',standard_peaks,format_func=peak_label,key='assay-standard-peak-'+runkey)
                        unknown_peak=right.selectbox('Unknown peak in current run',unknown_peaks,format_func=peak_label,key='assay-unknown-peak-'+runkey)
                        result=calculate_single_point(standard_peak['area_response_seconds'],assay_config['standard_concentration'],unknown_peak['area_response_seconds'],assay_config['unknown_dilution_factor'],assay_config['response_factor'])
                        st.success(f"Estimated unknown concentration: {result['unknown_concentration']:.6g} {assay_config['units']}")
                        st.dataframe(pd.DataFrame([{
                            'Standard injection':assay_raw['injection'],'Unknown injection':summary['injection'],
                            'UV signal':assay_channel_label(assay_channel),'Standard peak':peak_label(standard_peak),
                            'Unknown peak':peak_label(unknown_peak),'Standard concentration':f"{assay_config['standard_concentration']:.8g} {assay_config['units']}",
                            'Unknown dilution factor':assay_config['unknown_dilution_factor'],'Response factor':assay_config['response_factor'],
                            'Estimated unknown concentration':f"{result['unknown_concentration']:.8g} {assay_config['units']}"
                        }]),hide_index=True,width='stretch')
                        st.caption('Formula: unknown area ÷ standard area × standard concentration × unknown dilution factor ÷ response factor. Provisional screening result; review peak identity, integration and preparation records before release use.')
        except Exception as exc:
            st.error('Could not open the selected standard: '+str(exc))

alt.data_transformers.disable_max_rows()

def clear_peak_highlight():
    st.session_state['peak_highlight']=None
    st.session_state['peak_selection_epoch']=st.session_state.get('peak_selection_epoch',0)+1


def active_peak(channel):
    highlight=st.session_state.get('peak_highlight') or {}
    if highlight.get('analysis')==str(folder) and highlight.get('channel')==channel['channel']:
        return highlight.get('peak_id')
    return None


def add_manual(channel,bounds):
    if not bounds:return
    lo,hi=bounds
    candidate={'kind':'add','channel':channel['channel'],'start_min':lo,'end_min':hi,'peak_id':'preview'}
    try:
        # Validate before saving so an accidental empty/negative region cannot
        # leave an invalid persistent edit.
        apply_channel_actions(read_array(folder/channel['trace_csv']),channel['peaks'],[candidate],settings)
        action=record_edit(edits_path,{k:v for k,v in candidate.items() if k!='peak_id'})
        st.session_state['pending_manual_peak']=action
        st.session_state.pop('integration_edit_error',None)
    except Exception as exc:st.session_state['integration_edit_error']=str(exc)


def delete_integral(channel,peak_id):
    if not peak_id:return
    try:
        if not any(p['peak_id']==peak_id for p in channel['peaks']):raise ValueError('Choose an existing integral to delete')
        record_edit(edits_path,{'kind':'delete','channel':channel['channel'],'peak_id':peak_id})
        st.session_state.pop('integration_edit_error',None)
        clear_peak_highlight()
    except Exception as exc:st.session_state['integration_edit_error']=str(exc)


def delete_integrals_box(channel,bounds):
    if bounds is None:return
    ids=peaks_in_box(read_array(folder/channel['trace_csv']),shown(channel),bounds)
    if not ids:return
    record_edit(edits_path,{'kind':'delete_many','channel':channel['channel'],'peak_ids':ids,'selection_bounds':list(bounds),'selection_rule':'Peak apex inside rectangle (time and absorbance)'})
    st.session_state.pop('integration_edit_error',None)
    clear_peak_highlight()


def undo_integral(channel,reset=False):
    undo_edits(edits_path,channel['channel'],reset)
    st.session_state.pop('integration_edit_error',None)
    clear_peak_highlight()


def select_peak(channel,peak_id):
    if peak_id and any(p['peak_id']==peak_id for p in shown(channel)):
        st.session_state['peak_highlight']={'analysis':str(folder),'channel':channel['channel'],'peak_id':peak_id}
        st.session_state['open_peak_spectra']={'analysis':str(folder),'channel':channel['channel'],'peak_id':peak_id}
    else:clear_peak_highlight()


outside_click=st.components.v2.component('mironova_peak_outside_click',js=OUTSIDE_CLICK_JS)
active_highlight=st.session_state.get('peak_highlight') or {}
outside_click(data={'active':active_highlight.get('analysis')==str(folder)},key='peak_outside_click',height=0,on_clear_change=clear_peak_highlight)


def shown(channel):return [p for p in channel['peaks'] if show_hidden or p.get('report_included',True)]

def plot_limits(channel,xy,uv):
    prefix='axis-'+runkey+'-'+channel['channel']+'-'
    low=float(min(0,np.min(xy[:,1])))
    high=float(max(np.max(xy[:,1]),1))
    step=float(max(1,10**math.floor(math.log10(high))/10))
    defaults={'xmin':float(xy[0,0]),'xmax':float(xy[-1,0]),'ymin':float(math.floor(low/step)*step),'ymax':float(math.ceil(high*1.1/step)*step)}
    def reset_axes():
        for name,value in defaults.items():st.session_state[prefix+name]=value
        st.session_state[prefix+'view_revision']=st.session_state.get(prefix+'view_revision',0)+1
    def zoom_out():
        current=[st.session_state.get(prefix+k,defaults[k]) for k in ['xmin','xmax','ymin','ymax']]
        for name,value in zip(['xmin','xmax','ymin','ymax'],wider_view(current,(defaults['xmin'],defaults['xmax']))):st.session_state[prefix+name]=value
        st.session_state[prefix+'view_revision']=st.session_state.get(prefix+'view_revision',0)+1
    cols=st.columns([1,1,1,1,.8,.8])
    xmin=cols[0].number_input('X start (min)',value=defaults['xmin'],step=.1,format='%.3f',key=prefix+'xmin')
    xmax=cols[1].number_input('X end (min)',value=defaults['xmax'],step=.1,format='%.3f',key=prefix+'xmax')
    unit='mAU' if uv else 'response'
    ymin=cols[2].number_input('Y minimum ('+unit+')',value=defaults['ymin'],step=step,key=prefix+'ymin')
    ymax=cols[3].number_input('Y maximum ('+unit+')',value=defaults['ymax'],step=step,key=prefix+'ymax')
    cols[4].button('Zoom out',on_click=zoom_out,key=prefix+'zoomout')
    cols[5].button('Reset axes',on_click=reset_axes,key=prefix+'reset')
    if xmin>=xmax or ymin>=ymax:
        st.warning('Each axis maximum must be greater than its minimum.')
        return None
    return alt.Scale(domain=[xmin,xmax],nice=False),alt.Scale(domain=[ymin,ymax],nice=False)


def trace_chart(channel,uv=True):
    xy=read_array(folder/channel['trace_csv'])
    correction=channel.get('blank_correction',{})
    traces=None;trace_names=[]
    main_label='Blank-corrected sample' if correction.get('applied') else 'Original sample'
    overlay_traces,missing_overlays=comparison_traces(channel,uv,comparisons)
    with st.expander('Chart controls — traces, axes and integration',expanded=False):
        if uv:
            traces=uv_display_traces(folder,channel,blank_folder if correction.get('applied') else None,blank_record if correction.get('applied') else None)
            traces.update(overlay_traces)
            trace_names=st.multiselect('Traces to show',list(traces),default=[main_label,'Baseline']+list(overlay_traces),key='traces-'+runkey+'-'+channel['channel']+'-'+str(bool(correction.get('applied')))+'-'+hashlib.sha256('|'.join(overlay_traces).encode()).hexdigest()[:8],help='Select one trace or any combination to overlay. This changes the view only.')
            if not trace_names:
                st.info('Select at least one trace to display.');return
            st.caption('Integration and peak areas use '+main_label.lower()+'. Baseline is the estimated baseline for that signal.')
            if 'Blank' in trace_names:
                st.caption('Blank shows the recorded absorbance, without a response multiplier.'+(' Its time axis is shifted by '+str(correction.get('shift_min',0.))+' min, matching the correction settings.' if correction.get('shift_min',0.) else ''))
        view_xy=np.concatenate([traces[n] for n in trace_names]) if uv else np.concatenate([xy[:,:2],*overlay_traces.values()])
        view_xy=view_xy[np.argsort(view_xy[:,0])]
        limits=plot_limits(channel,view_xy,uv)
        if limits is None:return
        xscale,yscale=limits
        prefix='axis-'+runkey+'-'+channel['channel']+'-'
        analysis_visible=uv and main_label in trace_names
        modes=['Inspect peaks','Pan / drag','Zoom to region','Manual integration','Delete integrals'] if analysis_visible else ['Pan / drag','Zoom to region']
        if st.session_state.get(prefix+'mode') not in modes:st.session_state[prefix+'mode']=modes[0]
        mode=st.radio('Mouse mode',modes,horizontal=True,key=prefix+'mode')
        if uv:
            a,b=st.columns(2)
            has_edits=any(action['channel']==channel['channel'] for action in integration_edits.get('actions',[]))
            a.button('Undo integral edit',on_click=undo_integral,args=(channel,),key=prefix+'undo-integral',disabled=not has_edits)
            b.button('Reset automatic integrals',on_click=undo_integral,args=(channel,True),key=prefix+'reset-integrals',disabled=not has_edits)
            if mode=='Manual integration':st.info('Drag horizontally from the peak start to its end, then release. Area is measured above the displayed baseline. Overlapping integrals are replaced. Bounds must be within the active integration time range.')
            elif mode=='Delete integrals':st.info('Drag a box around the peak markers to remove their integrals. Peaks with their apex inside the box are removed together. Undo restores the whole group.')
        st.caption('Zoom to region: drag a rectangle around the peaks to enlarge it. Mouse wheel zooms; Zoom out widens the view; Reset axes restores the full run. Plot zoom does not change integration.')
    if missing_overlays:st.caption('No matching '+('UV wavelength/settings' if uv else 'MS channel')+' in: '+', '.join(missing_overlays))
    def apply_box_zoom(chartkey):
        bounds=box_zoom_bounds(st.session_state.get(chartkey,{}))
        if bounds:
            for name,value in zip(['xmin','xmax','ymin','ymax'],bounds):st.session_state[prefix+name]=value
            st.session_state[prefix+'view_revision']=st.session_state.get(prefix+'view_revision',0)+1
    def draw_zoomable(chart):
        brush=alt.selection_interval(name='zoom_region',encodings=['x','y'],translate=False,zoom=False,clear=False)
        wheel=alt.selection_interval(name='wheel_zoom',bind='scales',encodings=['x','y'],translate=False,zoom='wheel!')
        chartkey='zoom-'+channel['channel']+'-'+runkey+'-'+str(st.session_state.get(prefix+'view_revision',0))
        st.altair_chart(chart.add_params(brush,wheel),width='stretch',key=chartkey,on_select=lambda:apply_box_zoom(chartkey),selection_mode=['zoom_region'])
    if not uv:
        ms_traces={'Current injection':xy[:,:2],**overlay_traces}
        df=pd.concat([pd.DataFrame({'Time (min)':arr[:,0],'TIC response':arr[:,1],'Trace':name}) for name,arr in ms_traces.items()],ignore_index=True)
        chart=alt.Chart(df).mark_line(strokeWidth=1,clip=True).encode(x=alt.X('Time (min):Q',scale=xscale),y=alt.Y('TIC response:Q',scale=yscale),color='Trace:N',tooltip=['Time (min):Q','TIC response:Q','Trace:N']).properties(height=280)
        if mode=='Zoom to region':draw_zoomable(chart)
        else:st.altair_chart(chart.interactive(),width='stretch',key='trace-'+channel['channel']+'-'+runkey)
        return
    highlighted=active_peak(channel)
    colors={'Original sample':'#23598D','Blank':'#8E6BB5','Baseline':'#B24F39','Blank-corrected sample':'#00897B'}
    overlay_colors=['#CC7A00','#AA3377','#66A61E','#0077BB','#A6761D','#44AA99']
    colors.update({name:overlay_colors[i%len(overlay_colors)] for i,name in enumerate(overlay_traces)})
    trace_colors=[colors[n] for n in trace_names]
    df=pd.concat([pd.DataFrame({'Time (min)':traces[n][:,0],'Response':traces[n][:,1],'Trace':n}) for n in trace_names],ignore_index=True)
    lines=alt.Chart(df).mark_line(strokeWidth=1,clip=True).encode(
        x=alt.X('Time (min):Q',scale=xscale),y=alt.Y('Response:Q',scale=yscale,title='Absorbance (mAU)' if uv else 'MS response'),
        color=alt.Color('Trace:N',scale=alt.Scale(domain=trace_names,range=trace_colors)),tooltip=['Time (min):Q','Response:Q','Trace:N'])
    layers=[];areas=[];markers=[]
    if settings.integration_start_min is not None:
        window=pd.DataFrame({'From':[settings.integration_start_min],'To':[settings.integration_end_min]})
        layers.append(alt.Chart(window).mark_rect(color='#2A8A65',opacity=.08,clip=True).encode(x=alt.X('From:Q',scale=xscale),x2='To:Q'))
        boundaries=pd.DataFrame({'Time (min)':[settings.integration_start_min,settings.integration_end_min]})
        layers.append(alt.Chart(boundaries).mark_rule(color='#2A8A65',strokeDash=[5,4],clip=True).encode(x=alt.X('Time (min):Q',scale=xscale)))
    for p in (shown(channel) if analysis_visible else []):
        mask=(xy[:,0]>=p['start_min'])&(xy[:,0]<=p['end_min'])
        # Add interpolated boundary points, matching the numerical integral.
        times=np.r_[p['start_min'],xy[mask,0],p['end_min']]
        for t in np.unique(times):areas.append({'Time (min)':float(t),'Raw':float(np.interp(t,xy[:,0],xy[:,1])),'Baseline':float(np.interp(t,xy[:,0],xy[:,2])),'peak_key':p['peak_id'],'Peak ID':p['peak_id'],'Peak':p.get('peak_name',p['peak_id']),'Highlight':p['peak_id']==highlighted})
        markers.append({'Time (min)':p['rt_min'],'Response':float(np.interp(p['rt_min'],xy[:,0],xy[:,1])),'peak_key':p['peak_id'],'Peak ID':p['peak_id'],'Peak':p.get('peak_name',p['peak_id']),'Highlight':p['peak_id']==highlighted,'Area %':p.get('area_pct_retained')})
    if areas:
        area_chart=alt.Chart(pd.DataFrame(areas)).mark_area(opacity=.2,color='#23598D',clip=True).encode(x=alt.X('Time (min):Q',scale=xscale),y=alt.Y('Raw:Q',scale=yscale,title='Absorbance (mAU)'),y2='Baseline:Q',detail=['Peak ID:N','peak_key:N'])
        layers.append(area_chart)
        selected_area=[row for row in areas if row['Highlight']]
        if selected_area:
            layers.append(alt.Chart(pd.DataFrame(selected_area)).mark_area(opacity=.6,color='#E77624',clip=True).encode(x=alt.X('Time (min):Q',scale=xscale),y=alt.Y('Raw:Q',scale=yscale,title='Absorbance (mAU)'),y2='Baseline:Q',detail=['Peak ID:N','peak_key:N']))
    if mode in ['Zoom to region','Manual integration','Delete integrals']:
        # Put the region selector on the line layer, whose x/y fields are
        # unambiguous even when integration-window or area layers are present.
        brush=alt.selection_interval(name='zoom_region',encodings=['x','y'],translate=False,zoom=False,clear=False)
        wheel=alt.selection_interval(name='wheel_zoom',bind='scales',encodings=['x','y'],translate=False,zoom='wheel!')
        if mode=='Manual integration':brush=alt.selection_interval(name='manual_region',encodings=['x'],translate=False,zoom=False,clear=False)
        elif mode=='Delete integrals':brush=alt.selection_interval(name='delete_region',encodings=['x','y'],translate=False,zoom=False,clear=False)
        lines=lines.add_params(brush,wheel)
    layers.append(lines)
    if markers:
        points=alt.Chart(pd.DataFrame(markers)).encode(x=alt.X('Time (min):Q',scale=xscale),y=alt.Y('Response:Q',scale=yscale),tooltip=['Peak:N','Peak ID:N',alt.Tooltip('Time (min):Q',format='.4f'),alt.Tooltip('Area %:Q',format='.3f')])
        layers.extend([points.mark_point(clip=True).encode(size=alt.condition(alt.datum.Highlight,alt.value(110),alt.value(28)),color=alt.condition(alt.datum.Highlight,alt.value('#E77624'),alt.value('#18394D'))),points.mark_text(dy=-12,fontSize=10,clip=True).encode(text='Peak:N',color=alt.condition(alt.datum.Highlight,alt.value('#B34D0A'),alt.value('#18394D')))])
    initial={'value':[{'peak_key':highlighted}]} if highlighted and mode!='Delete integrals' else {}
    pick=alt.selection_point(name='peak_pick',fields=['peak_key'],toggle=False,
        clear="click[!event.item || !event.item.datum || !event.item.datum.peak_key]",**initial)
    chart=alt.layer(*layers).properties(height=280)
    edit_revision=hashlib.sha256(json.dumps(integration_edits,sort_keys=True).encode()).hexdigest()[:10]
    if mode=='Zoom to region':
        chartkey='zoom-'+channel['channel']+'-'+runkey+'-'+str(st.session_state.get(prefix+'view_revision',0))
        st.altair_chart(chart,width='stretch',key=chartkey,on_select=lambda:apply_box_zoom(chartkey),selection_mode=['zoom_region'])
    elif mode=='Manual integration':
        chartkey='manual-'+channel['channel']+'-'+runkey+'-'+edit_revision
        st.altair_chart(chart,width='stretch',key=chartkey,on_select=lambda:add_manual(channel,manual_region(st.session_state.get(chartkey,{}))),selection_mode=['manual_region'])
    elif mode=='Delete integrals':
        chartkey='delete-box-'+channel['channel']+'-'+runkey+'-'+edit_revision+'-'+str(st.session_state.get(prefix+'view_revision',0))
        st.altair_chart(chart,width='stretch',key=chartkey,on_select=lambda:delete_integrals_box(channel,deletion_region(st.session_state.get(chartkey,{}))),selection_mode=['delete_region'])
    elif mode=='Pan / drag':
        st.altair_chart(chart.interactive(),width='stretch',key='pan-'+channel['channel']+'-'+runkey)
    else:
        chart=chart.add_params(pick).interactive()
        chartkey='trace-'+channel['channel']+'-'+runkey+'-'+str(st.session_state.get('peak_selection_epoch',0))+'-'+str(st.session_state.get(prefix+'view_revision',0))+'-'+edit_revision
        def clicked_chart():
            peak_id=chart_peak_id(st.session_state.get(chartkey,{}))
            select_peak(channel,peak_id)
        st.altair_chart(chart,width='stretch',key=chartkey,on_select=clicked_chart,selection_mode=['peak_pick'])


def peak_table(channel,uv=True):
    rows=[]
    for p in shown(channel):
        row={'Peak':p.get('peak_name',p['peak_id']),'RT (min)':p['rt_min'],'Area %':p.get('area_pct_retained'),'Area':p['area_response_seconds']}
        if uv:
            links=[a for a in summary['mass_associations'] if a['uv_channel']==channel['channel'] and a['uv_peak_id']==p['peak_id']]
            for a in links:row[a['ms_channel']+' m/z']=a.get('base_peak_mz')
        else:row['Base m/z']=p.get('spectrum',{}).get('base_peak_mz')
        if show_hidden:
            row['Original area %']=p.get('area_pct_all_detected')
            row['Shown in report']=p['report_included']
        rows.append(row)
    if rows:
        df=pd.DataFrame(rows)
        fmt={c:'{:.4f}' if c=='RT (min)' else '{:.3f}' if c in ['Area','Area %','Original area %'] else '{:.2f}' for c in df.select_dtypes(include='number').columns}
        styled=df.style.format(fmt,na_rep='n/e').set_properties(subset=['RT (min)','Area %'],**{'font-weight':'bold','font-size':'16px'})
        tablekey='peak-table-'+runkey+'-'+channel['channel']+'-'+hashlib.sha256((str(folder)+summary['peak_name_revision']+str([p['peak_id'] for p in shown(channel)])).encode()).hexdigest()[:12]+'-'+str(st.session_state.get('peak_selection_epoch',0))
        def clicked_peak():
            index=selected_library_index(st.session_state.get(tablekey,{}),len(rows))
            peak_id=shown(channel)[index]['peak_id'] if index is not None else None
            select_peak(channel,peak_id)
        st.dataframe(styled,hide_index=True,width='stretch',key=tablekey,on_select=clicked_peak,selection_mode='single-cell')
        chosen_peak=next((p for p in shown(channel) if p['peak_id']==active_peak(channel)),None)
        if chosen_peak:
            with st.container(key='peak-name-editor-'+runkey+'-'+channel['channel']):
                with st.form('rename-'+runkey+'-'+channel['channel']+'-'+chosen_peak['peak_id']):
                    new_name=st.text_input('Peak name',value=chosen_peak.get('peak_name',chosen_peak['peak_id']),max_chars=80,help='User label only; this does not confirm chemical identity or alter integration.')
                    st.caption(f"Original ID: {chosen_peak['peak_id']} · RT {chosen_peak['rt_min']:.4f} min")
                    save_name=st.form_submit_button('Save peak name')
                if save_name:
                    cleaned=clean_peak_name(new_name)
                    peak_names[annotation_peak_key(channel['channel'],chosen_peak)]=cleaned
                    write_saved(names_path,peak_names)
                    st.rerun()
                st.button('Open peak spectra',on_click=lambda:st.session_state.update(open_peak_spectra={'analysis':str(folder),'channel':channel['channel'],'peak_id':chosen_peak['peak_id']}),key='spectra-'+channel['channel']+'-'+runkey)
                st.button('Delete selected integral',on_click=delete_integral,args=(channel,chosen_peak['peak_id']),key='delete-'+channel['channel']+'-'+runkey)
                st.button('Clear peak highlight',on_click=clear_peak_highlight,key='clear-'+channel['channel']+'-'+runkey)
    else:st.info('No peaks meet the current display cutoffs.')
    hidden=[p for p in channel['peaks'] if not p.get('report_included',True)]
    st.caption(f"{len(channel['peaks'])} detected · {len(hidden)} hidden ({sum(p.get('area_pct_all_detected') or 0 for p in hidden):.3f}% of original area). Remaining Area % is normalized to 100%; excluded peaks have no share. All measured areas and original percentages remain in exports. Areas: {channel['area_unit']}.")

uvs=sorted(summary['uv_channels'],key=lambda c:(signal_wavelength(c) is None,signal_wavelength(c) or 0,c['channel']))
tabs=st.tabs(['Chromatograms & peaks','Mass spectra','Files & metadata','PDF report'])
with tabs[0]:
    st.caption('Open a wavelength or MS channel to review its chromatogram and table. Close it when you finish.')
    for channel in uvs:
        with st.expander(uv_label(channel),expanded=False,key='uv-section-'+runkey+'-'+channel['channel']):
            correction=channel.get('blank_correction',{})
            if correction.get('applied'):st.caption(f"UV blank corrected using {correction['blank_injection']}; scale {correction['scale']:g}; shift {correction['shift_min']:+g} min; shared coverage {correction['coverage_min'][0]:.3f}–{correction['coverage_min'][1]:.3f} min. Negative residuals are retained.")
            if correction.get('applied'):
                coverage=correction['coverage_min']
                effective_lo=max(coverage[0],settings.integration_start_min if settings.integration_start_min is not None else coverage[0])
                effective_hi=min(coverage[1],settings.integration_end_min if settings.integration_end_min is not None else coverage[1])
                if effective_lo>=effective_hi:st.warning('The selected integration window is outside this channel’s shared sample/blank coverage. No peaks are integrated here. Change the time window to include the coverage shown above.')
                elif settings.integration_start_min is not None and (settings.integration_start_min<coverage[0] or settings.integration_end_min>coverage[1]):st.info(f'Blank-corrected integration uses {effective_lo:.3f}–{effective_hi:.3f} min for this wavelength, limited to the shared sample/blank coverage.')
            if correction.get('review_negative_residual'):st.warning(f"Substantial negative UV residual ({correction['minimum_residual_mau']:.1f} mAU). Review the blank, scale and time alignment before interpreting peak areas. Overlay the original and blank traces for comparison.")
            if summary.get('blank_correction') and not correction.get('applied'):st.caption('This UV channel is uncorrected: '+correction.get('reason','no matching blank'))
            trace_chart(channel);peak_table(channel)
    for channel in summary['ms_channels']:
        with st.expander('Mass chromatogram '+channel['channel'],expanded=False,key='ms-section-'+runkey+'-'+channel['channel']):
            trace_chart(channel,False)
            st.caption('Raw total-ion chromatogram — no baseline correction, peak picking or integration. Mass spectra are linked to the UV peaks.')
with tabs[1]:
    spectra=[]
    for c in uvs:
        ids={p['peak_id'] for p in shown(c)}
        for a in summary['mass_associations']:
            if a['uv_channel']==c['channel'] and a['uv_peak_id'] in ids and a.get('spectrum_csv'):
                spectra.append((f"{uv_label(c)} · {a.get('uv_peak_name',a['uv_peak_id'])} · RT {a['uv_rt_min']:.4f} · {a['ms_channel']}",a['spectrum_csv'],a))
    for c in summary['ms_channels']:
        for p in shown(c):
            if p.get('spectrum_csv'):spectra.append((f"{c['channel']} · {p['peak_id']} · RT {p['rt_min']:.4f}",p['spectrum_csv'],p['spectrum']))
    if spectra:
        index=st.selectbox('Choose a peak / MS scan',range(len(spectra)),format_func=lambda i:spectra[i][0])
        label,csvpath,info=spectra[index];xy=read_array(folder/csvpath);xy=xy[xy[:,1]>0]
        if len(xy):
            df=pd.DataFrame({'m/z':xy[:,0],'Relative intensity %':xy[:,1]/xy[:,1].max()*100})
            chart=alt.Chart(df).mark_rule(color='#23598D').encode(x='m/z:Q',y='Relative intensity %:Q',y2=alt.datum(0),tooltip=[alt.Tooltip('m/z:Q',format='.2f'),alt.Tooltip('Relative intensity %:Q',format='.2f')]).properties(height=320).interactive()
            st.altair_chart(chart,width='stretch');st.dataframe(df.sort_values('Relative intensity %',ascending=False).head(20),hide_index=True,width='stretch')
        st.caption(f"Observed scan RT {info.get('scan_rt_min',0):.4f} min. No background subtraction or confirmed compound identity.")
    else:st.info('No matching mass spectra for the currently shown peaks.')
with tabs[2]:
    st.caption('UV peaks use an estimated baseline and minimum prominence. Broad bumps require analyst review. Area % is detector response, not mass purity; nearby m/z does not assign identity.')
    st.write('Sample metadata');st.json(summary['sample_metadata'])
    st.write('Acquisition headers');st.json({uv_label(c)+' ('+c['channel']+')':c['metadata'] for c in uvs})
    st.write('Files in the injection folder')
    st.dataframe([{'File':p.name,'Bytes':p.stat().st_size} for p in sorted(Path(chosen).iterdir()) if p.is_file()],hide_index=True,width='stretch')
    if summary['errors']:st.json(summary['errors'])
    st.caption('Library search reads filenames and optional SAMPLE.XML names/comments. Acquisition headers and signals are read only when you select an injection.')
with tabs[3]:
    st.write('Create a printable report using the integration range and display cutoffs currently shown above.')
    viewkey=hashlib.sha256((str(folder)+json.dumps(summary['settings'],sort_keys=True)+summary['peak_name_revision']).encode()).hexdigest()
    exports=st.session_state.setdefault('report_exports',{})
    if st.button('Create report',type='primary'):
        try:
            with st.spinner('Creating your PDF…'):report=export_report(folder,summary,store)
            exports[viewkey]=str(report)
        except Exception as exc:st.error(str(exc))
    report=Path(exports[viewkey]) if viewkey in exports else None
    if report and (report/'single_run_report.pdf').is_file():
        st.session_state['single_report']=(str(report),summary)
        pdf=(report/'single_run_report.pdf').read_bytes()
        st.download_button('Open / download PDF',pdf,file_name=summary['injection']+' - report.pdf',mime='application/pdf')
        st.download_button('Download report and all processed data',report_zip(report),file_name=report.name+'.zip',mime='application/zip')
        encoded=base64.b64encode(pdf).decode()
        components.html('<iframe title="Printable PDF report" src="data:application/pdf;base64,'+encoded+'" width="100%" height="650"></iframe>',height=665)
        st.caption('Open the downloaded PDF to print. Your report is saved at '+str(report))


@st.dialog('Peak spectra',width='large')
def peak_spectra_dialog(channel,peak):
    st.subheader(peak.get('peak_name',peak['peak_id'])+' · RT '+f"{peak['rt_min']:.4f} min")
    st.caption('Raw spectra at the peak apex. UV chromatogram blank correction is not applied to these spectra. MS scans use the configured detector time offset; nearby scans do not confirm identity.')
    uvtab,mstab=st.tabs(['UV spectrum','Mass spectra'])
    with uvtab:
        cubes=summary.get('spectral_cubes',[])
        if not cubes:st.info('No UV spectral cube was recorded or decoded for this injection. Fixed-wavelength chromatograms alone cannot provide a UV spectrum.')
        for cube in cubes:
            try:
                scan_rt,spectrum=uv_spectrum(folder,cube,peak['rt_min'])
                if spectrum is None:
                    st.info(cube['file']+': no UV spectrum near this peak RT.');continue
                unit=cube.get('metadata',{}).get('unit','instrument response')
                frame=pd.DataFrame({'Wavelength (nm)':spectrum[:,0],'Absorbance ('+unit+')':spectrum[:,1]})
                y=frame.columns[1]
                chart=alt.Chart(frame).mark_line(color='#23598D').encode(x='Wavelength (nm):Q',y=alt.Y(y+':Q'),tooltip=['Wavelength (nm):Q',alt.Tooltip(y+':Q')]).properties(height=290).interactive()
                st.altair_chart(chart,width='stretch');st.caption(cube['file']+' · recorded scan RT '+f'{scan_rt:.4f} min')
                st.download_button('Download UV spectrum CSV',frame.to_csv(index=False),file_name=peak['peak_id']+'_'+Path(cube['file']).stem+'_UV.csv',mime='text/csv',key='uv-spectrum-download-'+cube['file'])
            except Exception as exc:st.error('Could not open UV spectrum: '+str(exc))
    with mstab:
        links=[a for a in summary['mass_associations'] if a['uv_channel']==channel['channel'] and a['uv_peak_id']==peak['peak_id']]
        available=[a for a in links if a.get('spectrum_csv')]
        if not available:st.info('No nearby mass spectrum is available for this peak.')
        for link in available:
            st.write(link['ms_channel'])
            spectrum=read_array(folder/link['spectrum_csv']);spectrum=spectrum[spectrum[:,1]>0]
            if not len(spectrum):st.info('No positive ion intensities in this scan.');continue
            frame=pd.DataFrame({'m/z':spectrum[:,0],'Relative intensity %':spectrum[:,1]/spectrum[:,1].max()*100})
            chart=alt.Chart(frame).mark_rule(color='#23598D').encode(x='m/z:Q',y='Relative intensity %:Q',y2=alt.datum(0),tooltip=[alt.Tooltip('m/z:Q',format='.2f'),alt.Tooltip('Relative intensity %:Q',format='.2f')]).properties(height=260).interactive()
            st.altair_chart(chart,width='stretch');st.caption('Recorded scan RT '+f"{link['scan_rt_min']:.4f} min")
            st.dataframe(frame.sort_values('Relative intensity %',ascending=False).head(20),hide_index=True,width='stretch')
            st.download_button('Download mass spectrum CSV',frame.to_csv(index=False),file_name=peak['peak_id']+'_'+link['ms_channel']+'_MS.csv',mime='text/csv',key='ms-spectrum-download-'+link['ms_channel'])


spectrum_request=st.session_state.pop('open_peak_spectra',None)
if spectrum_request and spectrum_request.get('analysis')==str(folder):
    requested_channel=next((c for c in uvs if c['channel']==spectrum_request['channel']),None)
    requested_peak=next((p for p in shown(requested_channel) if p['peak_id']==spectrum_request['peak_id']),None) if requested_channel else None
    if requested_peak:peak_spectra_dialog(requested_channel,requested_peak)
