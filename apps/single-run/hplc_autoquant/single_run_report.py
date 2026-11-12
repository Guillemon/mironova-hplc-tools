"""PDF for one injection; uses measured traces, areas and nearby spectra."""
from pathlib import Path
from xml.sax.saxutils import escape
import json
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet,ParagraphStyle
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,Image,PageBreak,KeepTogether


def signal_wavelength(channel):
    """Read the acquisition wavelength, excluding bandwidth and reference values."""
    signal=str(channel.get('metadata',{}).get('signal',''))
    match=re.search(r'\bSig\s*=\s*(\d+(?:\.\d+)?)',signal,re.IGNORECASE)
    return float(match.group(1)) if match else None


def uv_label(channel):
    wavelength=signal_wavelength(channel)
    return f'UV - {wavelength:g} nm' if wavelength is not None else 'UV - wavelength not recorded'


def create_pdf(run,directory):
    directory=Path(directory);styles=getSampleStyleSheet()
    uv_channels=sorted(run['uv_channels'],key=lambda c:(signal_wavelength(c) is None,signal_wavelength(c) or 0,c['channel']))
    navy=colors.HexColor('#18394D')
    styles.add(ParagraphStyle(name='RunTitle',fontName='Helvetica-Bold',fontSize=20,leading=24,textColor=navy,spaceAfter=12))
    styles.add(ParagraphStyle(name='RunBody',fontName='Helvetica',fontSize=9,leading=12,textColor=navy,spaceAfter=7))
    styles.add(ParagraphStyle(name='RunSmall',fontName='Helvetica',fontSize=7.5,leading=10,textColor=colors.HexColor('#52636E'),spaceAfter=5))
    styles.add(ParagraphStyle(name='RunCell',fontName='Helvetica',fontSize=7.5,leading=10,textColor=navy))
    styles.add(ParagraphStyle(name='RunHead',fontName='Helvetica-Bold',fontSize=7.5,leading=10,textColor=colors.white))
    styles.add(ParagraphStyle(name='RunImportant',fontName='Helvetica-Bold',fontSize=11,leading=14,textColor=navy))
    styles.add(ParagraphStyle(name='RunImportantHead',fontName='Helvetica-Bold',fontSize=10,leading=12,textColor=colors.white))
    story=[]
    def text(t,style='RunBody'):story.append(Paragraph(escape(str(t)),styles[style]))
    def heading(t):story.extend([Spacer(1,8),Paragraph(escape(t),styles['Heading2']),Spacer(1,5)])
    def table(headers,rows,widths,compact=False,emphasis=()):
        cells=[[Paragraph(escape(str(v)),styles['RunImportantHead' if i in emphasis else 'RunHead']) for i,v in enumerate(headers)]]+[[Paragraph(escape(str(v)),styles['RunImportant' if i in emphasis else 'RunCell']) for i,v in enumerate(row)] for row in rows]
        t=Table(cells,colWidths=widths,repeatRows=1,hAlign='LEFT')
        pad=2 if compact else 5
        t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),navy),('VALIGN',(0,0),(-1,-1),'TOP'),('TOPPADDING',(0,0),(-1,-1),pad),('BOTTOMPADDING',(0,0),(-1,-1),pad),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#F1F5F8')])]))
        story.extend([t,Spacer(1,8)])
    def savefig(fig,name,height=180):
        fig.tight_layout();path=directory/'figures'/name;fig.savefig(path,dpi=160);plt.close(fig)
        return Image(str(path),width=516,height=height)
    def peak_label(peak):
        name=peak.get('peak_name',peak['peak_id'])
        return name if name==peak['peak_id'] else name+' ('+peak['peak_id']+')'
    def trace_plot(channel,uv):
        xy=np.loadtxt(directory/channel['trace_csv'],delimiter=',',skiprows=1)
        fig,ax=plt.subplots(figsize=(10,3.5))
        ax.plot(xy[:,0],xy[:,1],lw=.8,color='#23598D',label='Blank-corrected UV' if uv and channel.get('blank_correction',{}).get('applied') else 'Raw signal')
        if uv:ax.plot(xy[:,0],xy[:,2],lw=.9,color='#B24F39',label='AsLS baseline')
        for peak in [p for p in channel['peaks'] if p.get('report_included',True)] if uv else []:
            sel=(xy[:,0]>=peak['start_min'])&(xy[:,0]<=peak['end_min'])
            ax.fill_between(xy[sel,0],xy[sel,1],xy[sel,2],color='#23598D',alpha=.2)
        plotted=sorted([p for p in channel['peaks'] if uv and p.get('report_included',True)],key=lambda p:p['prominence'],reverse=True)[:12]
        for p in plotted:
            ax.annotate(p.get('peak_name',p['peak_id'])[:40],(p['rt_min'],np.interp(p['rt_min'],xy[:,0],xy[:,1])),xytext=(0,7),textcoords='offset points',fontsize=7,ha='center')
        ax.set(xlabel='Time (min)',ylabel='Absorbance (mAU)' if uv else 'MS response (instrument units)',title=(uv_label(channel) if uv else channel['channel'])+(' - UV signal and integrated peaks' if uv else ' - raw TIC, no integration'))
        ax.margins(y=.18);ax.legend(fontsize=8);ax.grid(alpha=.1)
        return savefig(fig,channel['channel']+'_trace.png')
    text('Single-injection UV / MS report','RunTitle')
    text(run['injection'])
    text('PROVISIONAL - integrated peak, area-percentage, estimated S/N and observed m/z review. No assay, mass purity, concentration or compound identity is assigned.')
    if run['errors']:text('PARTIAL DATA: one or more files could not be processed. See the file/error list below.')
    sample=run['sample_metadata']
    if sample:
        table(['Sample record','Value'],[[k,sample.get(k) or 'Not recorded'] for k in ['Name','Description','Dilution','Multiplier']],[110,406])
    heading('Data read from this injection')
    table(['File','Recorded signal / axes','Data points','Detected peaks'],
          [[c['file'],uv_label(c),c['points'],f"{len(c['peaks'])} detected / {sum(p.get('report_included',True) for p in c['peaks'])} shown"] for c in uv_channels]+
          [[c['file'],'MS scans and measured m/z bins',f"{c['scans']} scans / {c['mz_bins']} bins",'Not integrated'] for c in run['ms_channels']]+
          [[c['file'],'UV time x wavelength spectral cube',f"{c['time_points']} x {c['wavelength_points']}",'Spectrum archive'] for c in run['spectral_cubes']],
          [80,240,115,81])
    heading('How to read the peak and mass tables')
    for note in run['review_notes']:text(note,'RunSmall')
    text(f"MS spectra use reader bin width {run['settings']['ms_bin_width_da']:g} Da. This is not a claim of high-resolution mass accuracy. MS-to-UV time offset: {run['settings']['ms_minus_uv_delay_min']:+g} min.",'RunSmall')
    if run['errors']:
        heading('Files requiring review')
        table(['File','Problem'],[[r['file'],r['error']] for r in run['errors']],[120,396])
    maxfig=run['settings']['spectra_plotted_per_channel']
    def metric(value,pattern='.3f'):
        return 'n/e' if value is None else format(value,pattern).replace('e+0','e').replace('e-0','e-')
    def filter_note(channel):
        if run['settings'].get('integration_start_min') is not None or run['settings'].get('integration_end_min') is not None:
            text(f"Integration window: {run['settings'].get('integration_start_min')} to {run['settings'].get('integration_end_min')} min. Measured areas refer to peaks in this window. Area % uses only peaks retained after the cutoffs. Full raw trace is shown.",'RunSmall')
        peaks=channel['peaks'];shown=[p for p in peaks if p.get('report_included',True)];hidden=[p for p in peaks if not p.get('report_included',True)]
        total=sum(p['area_response_seconds'] for p in peaks if p['area_response_seconds']>0)
        if total>0:
            shown_pct=sum(p.get('area_pct_all_detected') or 0 for p in shown);hidden_pct=sum(p.get('area_pct_all_detected') or 0 for p in hidden)
            retained_total=sum(p['area_response_seconds'] for p in shown if p['area_response_seconds']>0)
            if retained_total>0:
                text(f"Area-% denominator after cutoffs: {retained_total:.6g} response*s. The {len(shown)} retained peaks are normalized to 100%. Excluded: {len(hidden)} peaks / {hidden_pct:.4f}% of original detected area. Measured areas are unchanged; original percentages remain in the exports.",'RunSmall')
            else:text('No positive peak area remains after cutoffs; Area % cannot be calculated. Reduce the cutoffs to restore peaks.','RunSmall')
        else:text('Area % not estimated: no positive detected-area denominator. Review the integration.','RunSmall')
        if run['settings'].get('filter_report'):
            text(f"Cutoffs applied to original percentages: area >= {run['settings']['report_min_area_pct']:g}% and estimated S/N >= {run['settings']['report_min_estimated_snr']:g}. Known cutoffs apply independently; unavailable metrics require review. All peaks and reasons remain in CSV/JSON.",'RunSmall')
            cutoff=run['settings'].get('report_min_uv_area' if channel in run['uv_channels'] else 'report_min_ms_area',0)
            if cutoff:text(f"Absolute area display cutoff: {cutoff:g} {channel['area_unit']}.",'RunSmall')
        text('S/N is a screening estimate of peak height / robust noise sigma, not a validated LOD/LOQ. n/e = not estimable. Broad/shared-valley/edge features still require review even with a high S/N. Only retained UV integration regions are shaded on the raw trace.','RunSmall')
    selected_spectra=[]
    for channel in uv_channels:
        story.append(PageBreak());text(uv_label(channel),'RunTitle')
        text('Areas: mAU*s | '+('UV blank subtraction applied' if channel.get('blank_correction',{}).get('applied') else 'No blank subtraction'))
        correction=channel.get('blank_correction',{})
        if correction.get('applied'):text(f"Blank: {correction['blank_injection']} / {correction['blank_file']}; scale {correction['scale']:g}; shift {correction['shift_min']:+g} min. Shared coverage {correction['coverage_min'][0]:.3f}–{correction['coverage_min'][1]:.3f} min. Original sample and scaled blank traces are retained in the data exports.",'RunSmall')
        if correction.get('review_negative_residual'):text(f"REVIEW: substantial negative UV residual ({correction['minimum_residual_mau']:.1f} mAU). Verify blank suitability, scale and time alignment before interpreting areas.",'RunSmall')
        story.append(trace_plot(channel,True))
        rows=[]
        retained_total=sum(p['area_response_seconds'] for p in channel['peaks'] if p.get('report_included',True) and p['area_response_seconds']>0)
        for p in channel['peaks']:
            if not p.get('report_included',True):continue
            assoc=[a for a in run['mass_associations'] if a['uv_channel']==channel['channel'] and a['uv_peak_id']==p['peak_id']]
            labels=[]
            for a in assoc:
                mz=a.get('base_peak_mz')
                labels.append(a['ms_channel']+(': '+f'{mz:.2f}' if mz is not None else ': no nearby scan'))
            rows.append([peak_label(p),f"{p['rt_min']:.4f}",metric(100*p['area_response_seconds']/retained_total if p['area_response_seconds']>0 and retained_total>0 else None),f"{p['area_response_seconds']:.3f}",'; '.join(labels) or 'No MS'])
        heading('Peaks shown in this report')
        table(['Peak','RT (min)','Area %','Area (mAU*s)','Nearby MS base m/z'],rows or [['-','-','-','-','No shown peaks']],[75,75,75,115,176],emphasis=(1,2))
        filter_note(channel)
        text('Base m/z is the strongest bin in the nearby raw MS scan. MS background is not subtracted. It may be an adduct, fragment, matrix or solvent ion. Full spectra and scan time differences are in the data exports.','RunSmall')
        for peak in sorted([p for p in channel['peaks'] if p.get('report_included',True)],key=lambda p:p['prominence'],reverse=True)[:maxfig]:
            for a in run['mass_associations']:
                if a['uv_channel']==channel['channel'] and a['uv_peak_id']==peak['peak_id'] and a.get('spectrum_csv'):
                    selected_spectra.append((a['spectrum_csv'],f"{uv_label(channel)} {peak_label(peak)} at RT {peak['rt_min']:.4f} min / {a['ms_channel']}",a))
    for channel in run['ms_channels']:
        story.append(PageBreak());text(f"Mass chromatogram {channel['channel']}",'RunTitle')
        text('TIC = summed recorded ion intensities at each scan. No background subtraction; no polarity or molecular identity inferred from the filename.')
        story.append(trace_plot(channel,False))
        text('Raw TIC only: no MS baseline correction, peak picking or integration. Observed mass spectra are associated with the UV peaks.')
        text(channel['mass_axis_note']+' Entire decoded time/mass/intensity array is retained in '+channel['matrix_export']+'.','RunSmall')
    # Several UV channels may select the identical MS scan. Plot it once and
    # retain every association in the caption and exported link table.
    grouped={}
    for csvpath,label,info in selected_spectra:
        key=(info['ms_channel'],info['scan_rt_min'])
        if key not in grouped:grouped[key]=[csvpath,[label],info]
        else:grouped[key][1].append(label)
    selected_spectra=[(path,labels[0],dict(info,all_plot_links=labels)) for path,labels,info in grouped.values()]
    # Two spectrum figures per page, with exact scan/association provenance.
    for index,(csvpath,label,info) in enumerate(selected_spectra):
        if index%2==0:story.append(PageBreak());text('Observed mass spectra','RunTitle')
        xy=np.loadtxt(directory/csvpath,delimiter=',',skiprows=1)
        fig,ax=plt.subplots(figsize=(10,3.2))
        maximum=float(np.max(xy[:,1]))
        relative=xy[:,1]/maximum*100 if maximum>0 else xy[:,1]
        nonzero=relative>0
        ax.vlines(xy[nonzero,0],0,relative[nonzero],color='#23598D',lw=.7)
        # Label distinct local maxima, not every neighboring mass bin as a new compound.
        from scipy.signal import find_peaks
        tops,_=find_peaks(relative,prominence=max(2.,float(np.max(relative))*.02))
        labeled=[]
        spacing=max(5.,float(np.ptp(xy[:,0]))*.035)
        for j in sorted(tops,key=lambda j:relative[j],reverse=True):
            if all(abs(xy[j,0]-xy[k,0])>=spacing for k in labeled):
                ax.annotate(f'{xy[j,0]:.2f}',(xy[j,0],relative[j]),xytext=(0,5),textcoords='offset points',fontsize=7,ha='center')
                labeled.append(j)
            if len(labeled)>=6:break
        ax.set(xlabel='Observed m/z',ylabel='Relative intensity (%)',title=label,ylim=(0,120));ax.grid(alpha=.1)
        figure=savefig(fig,f'spectrum_{index+1}.png',165)
        caption=Paragraph(escape(f"Actual MS scan RT {info['scan_rt_min']:.4f} min. Links: {'; '.join(info['all_plot_links'])}. Exact per-link scan offsets remain in the data exports. No background subtraction or compound identification. File: {csvpath}"),styles['RunSmall'])
        story.append(KeepTogether([figure,caption,Spacer(1,9)]))
    for cube in run['spectral_cubes']:
        story.append(PageBreak());text('UV spectral data cube','RunTitle')
        archive=np.load(directory/cube['export']);x=archive['time_min'];wl=archive['wavelength_nm'];y=archive['intensity']
        fig,ax=plt.subplots(figsize=(10,4));low,high=np.percentile(y,[2,99])
        im=ax.imshow(y.T,aspect='auto',origin='lower',extent=[x[0],x[-1],wl[0],wl[-1]],vmin=low,vmax=high,cmap='viridis')
        ax.set(xlabel='Time (min)',ylabel='Wavelength (nm)',title=cube['file']);fig.colorbar(im,ax=ax,label='Recorded absorbance response')
        story.append(savefig(fig,'uv_cube.png',206))
        text('Display color scale clipped at the 2nd and 99th percentiles for readability. The complete decoded time, wavelength and intensity arrays are unchanged in '+cube['export']+'. This heatmap is not a spectral purity test.')
    story.append(PageBreak());text('Processing settings and traceability','RunTitle')
    table(['Setting','Value'],[[k,v] for k,v in run['settings'].items() if k not in ('ms_relative_prominence','report_min_ms_area')],[270,246],compact=True)
    text('Baseline: asymmetric least squares; lambda is scaled by (reference time step / median recorded time step)^4 to preserve the smoothing time scale. Peak detection: SciPy prominence. Boundaries: recorded prominence fraction, with shared-valley partitioning for overlapping bounds. Integrated area is the signed baseline-subtracted trapezoid multiplied by 60 for minutes-to-seconds conversion. This method requires comparison to analyst-reviewed integration before routine use.')
    text(f'At most {maxfig} most-prominent shown UV peaks per UV channel receive spectrum figures. MS TIC traces are not integrated. Every detected peak remains in CSV/JSON; every linked spectrum is exported. All decoded MS and UV spectrum matrices are retained in compressed NPZ files.')
    heading('Consumed original files - SHA-256')
    for source in run['sources']:text(source['file']+': '+source['sha256'],'RunSmall')
    text('Original files were read only and checked unchanged after processing. Source snapshots, CSV tables, spectral arrays, run_summary.json and a SHA-256 output manifest accompany this PDF. Unrecognized sidecars and unsupported proprietary files are not decoded. No calibration, assay or confirmed chemical identity is claimed.')
    def footer(canvas,doc):
        canvas.setFont('Helvetica',7.5);canvas.setFillColor(colors.HexColor('#52636E'))
        canvas.drawString(42,25,'Mironova Labs | Single-injection UV/MS | Provisional review')
        canvas.drawRightString(570,25,f'Page {doc.page}')
    SimpleDocTemplate(str(directory/'single_run_report.pdf'),pagesize=(612,792),leftMargin=48,rightMargin=48,topMargin=42,bottomMargin=45,
                      title='Single injection UV and mass spectrometry report',author='Mironova Labs').build(story,onFirstPage=footer,onLaterPages=footer)
