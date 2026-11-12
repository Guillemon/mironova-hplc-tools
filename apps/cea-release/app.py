"""Local CEA sequence review, calibration, assay and report application."""
from pathlib import Path,PurePosixPath
from dataclasses import asdict
import os,json,io,zipfile,stat
import numpy as np
import pandas as pd
import altair as alt
import streamlit as st
import matplotlib.pyplot as plt
from cea_release.raw import IntegrationSettings,import_sequence,save_json,digest
from cea_release.review import default_tables,calculate,fingerprint,clean
from cea_release.calculation import DEFAULT_CRITERIA
from cea_release.reports import export_report,calibration_figure

st.set_page_config(page_title='CEA Release Analysis',page_icon='🧪',layout='wide')
st.title('CEA release analysis')
st.caption('Review a complete Agilent sequence, confirm standards and preparations, calculate CEA assay, and export traceable reports.')
store=Path(os.environ.get('CEA_RELEASE_DATA',str(Path.cwd()/'cea_release_results'))).expanduser().resolve();store.mkdir(parents=True,exist_ok=True)
preference=store/'preferences.json'
try:last=json.loads(preference.read_text())
except (OSError,ValueError):last={}


def upload_sequence(upload):
    data=upload.getvalue();target=store/'uploads'/digest(data)
    if (target/'complete.json').is_file():return target/'sequence'
    root=target/'sequence';root.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(io.BytesIO(data)) as archive:
        if len(archive.infolist())>15000 or sum(v.file_size for v in archive.infolist())>1024**3:raise ValueError('Choose a smaller sequence ZIP (maximum 15,000 files / 1 GB expanded)')
        for member in archive.infolist():
            rel=PurePosixPath(member.filename.replace('\\','/'))
            if rel.is_absolute() or any(p in ['..',''] or ':' in p for p in rel.parts) or stat.S_ISLNK(member.external_attr>>16):raise ValueError('Unsafe ZIP path or link')
            if member.is_dir() or Path(rel.name).suffix.casefold() not in ['.ch','.xml','.txt']:continue
            dest=root.joinpath(*rel.parts);dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(archive.read(member))
    save_json(target/'complete.json',{'sha256':digest(data)});return root


with st.sidebar:
    st.header('Open a sequence')
    source_kind=st.radio('Data source',['Folder path','Sequence ZIP'])
    with st.form('load-sequence'):
        source=st.text_input('Agilent sequence folder',value=last.get('source',''),placeholder=r'C:\ChemStation\1\Data\CEA_Release') if source_kind=='Folder path' else None
        upload=st.file_uploader('ZIP containing the .D injection folders',type=['zip']) if source_kind=='Sequence ZIP' else None
        wave_mode=st.selectbox('UV wavelength',['Automatic (208, then 210 nm)','208 nm','210 nm','Other wavelength'])
        custom=st.number_input('Other wavelength (nm)',value=250.,min_value=100.,step=1.) if wave_mode=='Other wavelength' else None
        with st.expander('Target integration settings',expanded=False):
            rtmin=st.number_input('Target RT start (min)',value=10.3,min_value=0.,step=.1)
            rtmax=st.number_input('Target RT end (min)',value=11.4,min_value=0.,step=.1)
            prominence=st.number_input('Minimum target prominence (mAU)',value=100.,min_value=.01)
            half=st.number_input('Integration half-window (min)',value=.30,min_value=.001,step=.01,format='%.3f')
            guard=st.number_input('Baseline guard half-width (min)',value=.025,min_value=.001,step=.005,format='%.3f')
            minwidth=st.number_input('Minimum target width (min)',value=.04,min_value=.001,format='%.3f')
            maxwidth=st.number_input('Maximum target width (min)',value=.30,min_value=.001,format='%.3f')
        load=st.form_submit_button('Load chromatograms',type='primary')
    if load:
        try:
            if source_kind=='Sequence ZIP':
                if upload is None:raise ValueError('Choose a sequence ZIP first')
                source=str(upload_sequence(upload))
            wave=None if wave_mode.startswith('Automatic') else custom if custom is not None else float(wave_mode.split()[0])
            with st.spinner('Reading sequence and integrating the CEA target…'):
                folder,imported=import_sequence(source.strip().strip('"'),store,IntegrationSettings(rtmin,rtmax,prominence,half,guard,minwidth,maxwidth),wave)
            st.session_state['active_import']=str(folder);st.session_state.pop('calculation',None)
            if source_kind=='Folder path':save_json(preference,{'source':source})
        except Exception as exc:st.error(str(exc))
    st.caption('Raw files are read only. Results remain on this computer. This application does not issue final batch release approval.')

if not st.session_state.get('active_import'):
    st.info('Open the folder containing one sequence of .D injections. Start with the proposed roles and areas, then verify the preparation records before calculating.');st.stop()
folder=Path(st.session_state['active_import']);imported=json.loads((folder/'import.json').read_text(encoding='utf-8'));runid=imported['run_id'];key=runid[:16]
st.subheader(Path(imported['source']).name)
st.caption(f"{len(imported['records'])} injections | {imported['wavelength_nm']:g} nm | Import {key}")
read_errors=[r for r in imported['records'] if r['reader_status']!='EXTRACTED']
if read_errors:st.warning(f'{len(read_errors)} injection(s) could not be read. Resolve them or explicitly exclude with a reason.')
default_assign,default_preps=default_tables(imported)
reviewfile=store/'reviews'/(runid+'.json')
try:saved=json.loads(reviewfile.read_text())
except (OSError,ValueError):saved={}

with st.expander('Raw chromatograms and metadata',expanded=False):
    lookup={r['injection_id']:r for r in imported['records']}
    chosen=st.selectbox('Injection to review',list(lookup));record=lookup[chosen]
    st.caption(record['proposed_kind']+' | '+record.get('channel_file','no UV channel')+' | '+record['integration']['status'])
    if 'trace_npz' in record:
        with np.load(folder/record['trace_npz']) as data:x,y=data['time_min'],data['intensity_mau']
        df=pd.DataFrame({'Time (min)':x,'Raw absorbance (mAU)':y})
        line=alt.Chart(df).mark_line(color='#23598D').encode(x='Time (min):Q',y='Raw absorbance (mAU):Q',tooltip=['Time (min):Q','Raw absorbance (mAU):Q'])
        p=record['integration']
        if p.get('start_min') is not None:
            xx=np.r_[p['start_min'],x[(x>p['start_min'])&(x<p['end_min'])],p['end_min']]
            baseline=np.interp(xx,[p['start_min'],p['end_min']],[p['baseline_start_mau'],p['baseline_end_mau']])
            frame=pd.DataFrame({'Time (min)':xx,'Raw absorbance (mAU)':np.interp(xx,x,y),'Baseline (mAU)':baseline})
            area=alt.Chart(frame).mark_area(opacity=.25,color='#23598D').encode(x='Time (min):Q',y='Raw absorbance (mAU):Q',y2='Baseline (mAU):Q')
            line=alt.layer(area,line)
        st.altair_chart(line.properties(height=300).interactive(),width='stretch')
    st.json({'sample':record['sample_metadata'],'integration':record['integration'],'read_error':record.get('error')})

st.subheader('1. Review injection assignments')
st.caption('Areas come from raw CH signals. Standard mass / volume overrides the nominal concentration when both are entered. An area override or exclusion requires a reason. Blanks remain separate from standards and samples.')
frame=pd.DataFrame(saved.get('assignments',default_assign))
for c in ['area','nominal_concentration_mg_ml','standard_mass_mg','standard_volume_ml']:frame[c]=pd.to_numeric(frame[c],errors='coerce')
assign=st.data_editor(frame,hide_index=True,width='stretch',disabled=['injection_id'],key='assign-'+key,column_config={'kind':st.column_config.SelectboxColumn('Role',options=['standard','sample','solvent_blank','matrix_blank','qc','system_suitability','unassigned'],required=True),'area':st.column_config.NumberColumn('Reviewed area (mAU*s)',format='%.6f'),'nominal_concentration_mg_ml':st.column_config.NumberColumn('Nominal mg/mL',format='%.6f'),'standard_mass_mg':st.column_config.NumberColumn('Standard mass (mg)',format='%.4f'),'standard_volume_ml':st.column_config.NumberColumn('Standard volume (mL)',format='%.4f')})

st.subheader('2. Confirm sample preparations')
st.caption('Mass is formulation mass in mg; volume is final mL. Missing records must be entered. Use the same Sample / batch ID only for preparations of the same batch. Additional dilution defaults to 1 and must be verified.')
prepframe=pd.DataFrame(saved.get('preparations',default_preps),columns=['preparation_id','sample_id','formulation_mass_mg','final_volume_ml','additional_dilution_factor'])
for c in ['formulation_mass_mg','final_volume_ml','additional_dilution_factor']:prepframe[c]=pd.to_numeric(prepframe[c],errors='coerce')
preps=st.data_editor(prepframe,hide_index=True,width='stretch',num_rows='dynamic',key='preps-'+key)

st.subheader('3. Calculate')
method=st.radio('Calculation method',['Single external standard (STD 1 + check standards)','OLS calibration - free intercept'],key='method-'+key)
model='single_standard' if method.startswith('Single') else 'ols'
levels=[]
for r in clean(assign.to_dict('records')):
    if r['include'] and r['kind']=='standard':
        mass,volume=r.get('standard_mass_mg'),r.get('standard_volume_ml');c=mass/volume if mass is not None and volume is not None and volume>0 else r.get('nominal_concentration_mg_ml')
        if c is not None and c>0:levels.append(c)
levels=sorted(set(levels))
reference=st.selectbox('STD 1 nominal concentration (mg/mL)',levels,index=None,placeholder='Select the primary reference standard',key='reference-'+key) if model=='single_standard' else None
allow=st.checkbox('Allow a two-level exploratory OLS fit',value=False,key='two-level-'+key) if model=='ols' else False
purity_mode=st.radio('Reference purity',['Uncorrected estimate - purity not supplied','Apply verified reference purity'],key='purity-mode-'+key)
purity_confirmed=purity_mode.startswith('Apply');purity_percent=st.number_input('Verified purity (%)',value=float(saved.get('purity_fraction',1.)*100),min_value=.0001,max_value=100.,format='%.4f',key='purity-'+key) if purity_confirmed else 100.
with st.expander('Numerical review criteria',expanded=False):
    st.caption('Defaults are screening limits from the recovered engine. Verify against the laboratory method; passing is not release approval.')
    criteria_text=st.text_area('Criteria JSON',value=json.dumps(saved.get('criteria',DEFAULT_CRITERIA),indent=2),height=240,key='criteria-'+key)
analyst=st.text_input('Analyst',value=saved.get('analyst',''),key='analyst-'+key)
verified=st.checkbox('I checked target identity, roles, preparation records, standard concentrations and purity settings',key='verified-'+key)
try:criteria=json.loads(criteria_text)
except ValueError:criteria=None;st.error('Numerical criteria must be valid JSON')
payload=clean({'assignments':assign.to_dict('records'),'preparations':preps.to_dict('records'),'calibration_model':model,'reference_nominal_mg_ml':reference,
               'allow_two_levels':allow,'purity_confirmed':purity_confirmed,'purity_fraction':purity_percent/100,'criteria':criteria,'analyst':analyst,'verified':verified})
if st.button('Calculate results',type='primary',disabled=not verified or not analyst.strip() or criteria is None):
    try:
        analysis=calculate(imported,payload);st.session_state['calculation']=analysis
        reviewfile.parent.mkdir(exist_ok=True);save_json(reviewfile,payload)
    except Exception as exc:st.error(str(exc));st.session_state.pop('calculation',None)
analysis=st.session_state.get('calculation')
if analysis and (analysis['import_run_id']!=runid or analysis['input_fingerprint']!=fingerprint(payload)):
    st.warning('Inputs changed. Recalculate before viewing or exporting updated results.');analysis=None
if analysis:
    result=analysis['result'];st.subheader('CEA assay results')
    for r in result['sample_summaries']:st.metric(r['sample_id']+' - mean assay (% w/w)',f"{r['mean']:.4f}")
    st.caption(analysis['status']+' | '+('Verified reference purity applied once' if purity_confirmed else 'Uncorrected estimate; reference purity not supplied'))
    st.dataframe(result['preparations'],hide_index=True,width='stretch')
    fig=calibration_figure(result);st.pyplot(fig);plt.close(fig)
    st.dataframe(result['standards'],hide_index=True,width='stretch')
    with st.expander('Review findings',expanded=True):st.dataframe(result['flags'],hide_index=True,width='stretch')
    if st.button('Save report package',type='primary'):
        try:
            with st.spinner('Generating PDFs, tables and raw-data audit package…'):report=export_report(folder,imported,analysis,store)
            st.session_state['report_directory']=str(report);st.session_state['report_analysis_id']=analysis['analysis_id']
        except Exception as exc:st.error(str(exc))
    if st.session_state.get('report_analysis_id')==analysis['analysis_id']:
        report=Path(st.session_state['report_directory']);st.success('Report saved: '+str(report))
        st.download_button('Download analytical PDF',(report/'analytical_report.pdf').read_bytes(),file_name='CEA_analytical_report.pdf',mime='application/pdf')
        st.download_button('Download all chromatograms PDF',(report/'chromatographic_review.pdf').read_bytes(),file_name='CEA_chromatographic_review.pdf',mime='application/pdf')
        buffer=io.BytesIO()
        with zipfile.ZipFile(buffer,'w',zipfile.ZIP_DEFLATED) as archive:
            for p in sorted(report.rglob('*')):
                if p.is_file():archive.write(p,p.relative_to(report))
        st.download_button('Download complete report and data ZIP',buffer.getvalue(),file_name=report.name+'.zip',mime='application/zip')
