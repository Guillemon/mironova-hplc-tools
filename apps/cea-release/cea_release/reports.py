"""Portable PDF/CSV report packages tied to reviewed raw-source snapshots."""
from pathlib import Path
from datetime import datetime,timezone
import csv,shutil,uuid,json
from xml.sax.saxutils import escape
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet,ParagraphStyle
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,Image,PageBreak
from .raw import save_json,digest


def write_csv(path,rows):
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with Path(path).open('w',newline='',encoding='utf-8-sig') as f:
        writer=csv.DictWriter(f,fieldnames=keys);writer.writeheader()
        for r in rows:
            safe={k:json.dumps(v) if isinstance(v,(list,dict)) else v for k,v in r.items()}
            safe={k:("'"+v if isinstance(v,str) and v.startswith(('=','+','-','@')) else v) for k,v in safe.items()}
            writer.writerow(safe)


def calibration_figure(result):
    fig,ax=plt.subplots(figsize=(8,3.2));levels=result['standards'];x=np.array([r['concentration_mg_ml'] for r in levels]);y=np.array([r['mean'] for r in levels])
    ref=np.array([r.get('assay_reference',True) for r in levels]);ax.scatter(x[ref],y[ref],label='Calibration/reference means',color='#23598D')
    if np.any(~ref):ax.scatter(x[~ref],y[~ref],marker='s',label='Check standards',color='#D48232')
    line=np.linspace(0 if result['calibration_model']=='single_standard' else min(x)*.95,max(x)*1.05,100)
    ax.plot(line,line*result['slope']+result['intercept'],color='#23598D',lw=1,label='Reference response line' if result['calibration_model']=='single_standard' else 'OLS fitted curve')
    ax.set(xlabel='Effective CEA concentration (mg/mL)',ylabel='Mean area (mAU*s)');ax.grid(alpha=.15);ax.legend(fontsize=8);fig.tight_layout();return fig


def trace_figure(folder,record,settings):
    with np.load(Path(folder)/record['trace_npz']) as data:x,y=data['time_min'],data['intensity_mau']
    fig,axes=plt.subplots(1,2,figsize=(9,2.8),gridspec_kw={'width_ratios':[1,1.4]})
    p=record['integration']
    for ax in axes:
        ax.plot(x,y,lw=.7,color='#23598D');ax.set(xlabel='Time (min)',ylabel='mAU');ax.grid(alpha=.15)
        if p.get('start_min') is not None:
            lo,hi=p['start_min'],p['end_min'];xx=np.r_[lo,x[(x>lo)&(x<hi)],hi];base=np.interp(xx,[lo,hi],[p['baseline_start_mau'],p['baseline_end_mau']])
            ax.plot(xx,base,color='#B24F39',lw=1);ax.fill_between(xx,np.interp(xx,x,y),base,color='#23598D',alpha=.2)
            ax.axvline(lo,color='#B24F39',lw=.6,ls='--');ax.axvline(hi,color='#B24F39',lw=.6,ls='--')
    axes[1].set_xlim(settings['rt_min']-.3,settings['rt_max']+.3)
    mask=(x>=settings['rt_min']-.3)&(x<=settings['rt_max']+.3)
    if np.any(mask):
        low,high=float(np.min(y[mask])),float(np.max(y[mask]));pad=max(5.,(high-low)*.15);axes[1].set_ylim(low-pad,high+pad)
    axes[0].set_title('Full trace',fontsize=9);axes[1].set_title('Target window and local baseline',fontsize=9)
    fig.tight_layout();return fig


def create_pdfs(folder,imported,analysis):
    folder=Path(folder);result=analysis['result'];styles=getSampleStyleSheet();navy=colors.HexColor('#18394D')
    styles.add(ParagraphStyle(name='SmallReport',fontSize=8,leading=11,spaceAfter=6,textColor=navy))
    styles.add(ParagraphStyle(name='TableReport',fontSize=7,leading=9,textColor=navy))
    styles.add(ParagraphStyle(name='TableHead',fontSize=7,leading=9,fontName='Helvetica-Bold',textColor=colors.white))
    def p(value,style='SmallReport'):return Paragraph(escape(str(value)),styles[style])
    def table(headers,rows,widths):
        cells=[[p(v,'TableHead') for v in headers]]+[[p(v,'TableReport') for v in row] for row in rows]
        t=Table(cells,colWidths=widths,repeatRows=1,hAlign='LEFT');t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),navy),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#F1F5F8')]),('VALIGN',(0,0),(-1,-1),'TOP'),('TOPPADDING',(0,0),(-1,-1),5),('BOTTOMPADDING',(0,0),(-1,-1),5)]));return t
    def footer(canvas,doc):
        canvas.saveState();canvas.setFont('Helvetica',7);canvas.setFillColor(navy);canvas.drawString(42,25,'CEA release analysis - provisional; no final batch approval');canvas.drawRightString(570,25,str(doc.page));canvas.restoreState()
    def image(fig,name,height):
        path=folder/'figures'/name;fig.savefig(path,dpi=150);plt.close(fig);return Image(str(path),width=528,height=height)
    fmt=lambda v,n=4:'not estimable' if v is None else f'{v:.{n}f}'
    story=[p('CEA release analysis','Title'),p(analysis['status']),p('Analyst: '+analysis['analyst']+' | '+datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')),p('Sequence: '+imported['source']),p(f"UV wavelength: {imported['wavelength_nm']:g} nm | Analysis ID: {analysis['analysis_id'][:16]}")]
    story+=[p('Assay results','Heading2'),table(['Sample / batch','Preparations','Mean assay (% w/w)','Between-preparation RSD (%)'],[[r['sample_id'],r['n'],fmt(r['mean']),fmt(r['rsd_pct'])] for r in result['sample_summaries']],[160,70,145,153]),Spacer(1,10)]
    purity=analysis['payload'];story.append(p('Reference purity correction: '+(f"{purity['purity_fraction']*100:.4f}%" if purity['purity_confirmed'] else 'not supplied; uncorrected estimate (factor 1.0)')))
    story.extend([p('Calculation method','Heading2'),p(result['model']),p('Slope: '+fmt(result['slope'])+' | Intercept: '+fmt(result['intercept'])+' | R squared: '+('not applicable (single-standard ratio)' if result['r_squared'] is None else fmt(result['r_squared'],8))),p('Assay % w/w = vial concentration x final volume x additional dilution factor / formulation mass x 100. Replicate injections are averaged within preparation; preparations are averaged only within their assigned sample/batch.'),p('Raw peak integration uses the recovered median-guard local linear baseline. Peak identity, coelution, baseline and assay-method suitability remain subject to laboratory review. Numerical checks do not provide a final release decision.')])
    story.extend([PageBreak(),p('Calibration and standards','Title'),image(calibration_figure(result),'calibration.png',212),p('A single-standard response line is not a fitted calibration curve or a linearity validation.' if result['calibration_model']=='single_standard' else 'OLS uses equal weights for concentration means. A two-level exact fit does not establish linearity.'),table(['Effective mg/mL','n','Mean area','Area RSD %','Recovery %','Role'],[[fmt(r['concentration_mg_ml'],6),r['n'],fmt(r['mean'],2),fmt(r['rsd_pct'],3),fmt(r['recovery_pct'],3),'Reference' if r.get('assay_reference') else ('Check' if result['calibration_model']=='single_standard' else 'Curve')] for r in result['standards']],[95,35,110,85,90,113]),Spacer(1,10),p('Preparation results','Heading2'),table(['Preparation','Mass mg','Volume mL','Dilution','Vial mg/mL','Assay % w/w'],[[r['preparation_id'],fmt(r['formulation_mass_mg'],2),fmt(r['final_volume_ml'],2),fmt(r['additional_dilution_factor'],2),fmt(r['vial_concentration_mg_ml'],6),fmt(r['assay_pct_ww'])] for r in result['preparations']],[100,80,80,70,100,98])])
    story.extend([PageBreak(),p('Review findings','Title'),table(['Code','Subject','Finding'],[[r['code'],r['subject'],r['message']] for r in result['flags']],[125,140,263]),Spacer(1,10),p('Exclusions','Heading2')])
    if analysis['exclusions']:story.append(table(['Injection','Reason'],[[r['injection_id'],r['reason']] for r in analysis['exclusions']],[300,228]))
    else:story.append(p('No injections excluded. Solvent and matrix controls are reviewed separately from calibration and samples.'))
    story+=[p('Traceability','Heading2'),p('The package contains raw-source snapshots, every decoded trace, import settings, reviewed assignments, override/exclusion reasons, input fingerprint, calculation JSON and CSV tables. manifest.json contains SHA-256 hashes of all exported files. The chromatographic review shows the raw integrations; any reviewed area overrides are explicitly recorded in assignments.csv. No automatic blank subtraction is applied. Acceptance criteria are screening settings recorded in analysis.json and require method verification.')]
    SimpleDocTemplate(str(folder/'analytical_report.pdf'),pagesize=(612,792),leftMargin=42,rightMargin=42,topMargin=42,bottomMargin=42).build(story,onFirstPage=footer,onLaterPages=footer)
    story=[p('CEA chromatographic review','Title'),p(f"{imported['wavelength_nm']:g} nm - raw signals, local baselines and integration bounds. See analytical_report.pdf for calculations and review findings.")]
    plotted=0
    for i,r in enumerate(imported['records']):
        if 'trace_npz' not in r:continue
        if plotted and plotted%2==0:story.append(PageBreak())
        story.extend([p(r['injection_id'],'Heading3'),p(f"Role proposal: {r['proposed_kind']} | Raw area: {fmt(r['integration'].get('area_mau_s'),3)} mAU*s | {r['integration']['status']}"),image(trace_figure(folder,r,imported['settings']),f'trace-{i:04d}.png',166),Spacer(1,10)]);plotted+=1
    SimpleDocTemplate(str(folder/'chromatographic_review.pdf'),pagesize=(612,792),leftMargin=42,rightMargin=42,topMargin=38,bottomMargin=42).build(story,onFirstPage=footer,onLaterPages=footer)


def export_report(import_folder,imported,analysis,store):
    import_folder=Path(import_folder)
    if analysis['import_run_id']!=imported['run_id']:raise ValueError('The calculation belongs to a different import')
    for s in imported['sources']:
        if digest((import_folder/'source_snapshots'/s['relative_file']).read_bytes())!=s['sha256']:raise ValueError('A saved raw-source snapshot changed; reimport the sequence')
    out=Path(store)/'reports'/('CEA-'+datetime.now(timezone.utc).strftime('%Y%m%d-%H%M%S')+'-'+uuid.uuid4().hex[:8]);out.mkdir(parents=True)
    try:
        shutil.copytree(import_folder,out,dirs_exist_ok=True);(out/'figures').mkdir(exist_ok=True)
        save_json(out/'analysis.json',analysis)
        result=analysis['result']
        for name,rows in [('assignments',analysis['payload']['assignments']),('standards',result['standards']),('preparations',result['preparations']),('injections',result['injections']),('review_findings',result['flags']),('exclusions',analysis['exclusions']),('controls',analysis['controls'])]:write_csv(out/(name+'.csv'),rows)
        create_pdfs(out,imported,analysis)
        save_json(out/'manifest.json',{'algorithm':'SHA-256','files':{p.relative_to(out).as_posix():digest(p.read_bytes()) for p in sorted(out.rglob('*')) if p.is_file() and p.name!='manifest.json'}})
        return out
    except Exception:shutil.rmtree(out,ignore_errors=True);raise
