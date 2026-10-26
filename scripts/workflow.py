"""NMR workflow: solvent survey, preview diagnostics, reviewed quantification.

python workflow.py survey --config sample.json --out survey
python workflow.py inspect --config sample.json --out preview
python workflow.py quantify --config sample.json --review review.json --out results
The review is the agent's evidence record, not a request for user approval.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.signal import find_peaks
from scipy.integrate import cumulative_trapezoid
import engine as e
import solvent
import copy


def config_hash(cfg):
    stable = {k:v for k,v in cfg.items() if k != 'output_dir'}
    return hashlib.sha256(json.dumps(stable,sort_keys=True,ensure_ascii=False).encode()).hexdigest()


def candidate_peaks(x,y,search,noise):
    m=e.mask(x,search)
    xx,yy=x[m],y[m]
    if len(xx)<5: raise ValueError('Search window is empty or too narrow.')
    # Proposals only: no automatic chemical assignment or integration boundary.
    prominence=max(5*noise,.015*float(np.ptp(yy)),1e-12)
    peaks,props=find_peaks(yy,prominence=prominence)
    return [{'ppm':float(xx[i]),'intensity':float(yy[i]),'prominence':float(p)}
            for i,p in zip(peaks,props['prominences'])]


def assert_review(cfg,review,source_hash):
    if review.get('source_sha256')!=source_hash:
        raise ValueError('Review belongs to a different source file; inspect this spectrum first.')
    if review.get('config_sha256')!=config_hash(cfg):
        raise ValueError('Configuration changed since review; regenerate and inspect diagnostics.')
    for key in ['solvent_and_reference_assessment','sample_identity_basis','phase_assessment','baseline_assessment',
                'main_assignment_and_bounds','end_assignment_and_bounds',
                'overlap_assessment','boundary_sensitivity_assessment','formula_applicability']:
        if not isinstance(review.get(key),str) or len(review[key].strip())<12:
            raise ValueError(f'Missing substantive review field: {key}')
    if review.get('decision')!='quantify_provisionally':
        raise ValueError('Unresolved review: do not calculate/report Mn as usable.')
    conflicts=solvent.region_conflicts(cfg)
    if conflicts and len(review.get('solvent_overlap_resolution','').strip())<12:
        raise ValueError('Solvent/water overlap flagged: resolve regions or document why quantification remains meaningful; do not silently cut the peak.')


def survey(cfg,output):
    """Initial view: read solvent, keep native ppm axis, no phase/integrals/Mn."""
    src=Path(cfg['input_jdf']);before=e.sha256(src)
    dic,fid=e.ng.jeol.read(str(src))
    raw=copy.deepcopy(cfg);raw['reference']=None
    x,s,meta=e.transform(dic,fid,raw)
    info=solvent.describe(dic['parameters'],cfg)
    expected=info['suggested_residual_ppm']
    candidates=[]
    fig,axs=e.plt.subplots(2,1,figsize=(12,7),constrained_layout=True)
    axs[0].plot(x,abs(s),lw=.6);axs[0].set_xlim(x[-1],x[0])
    axs[0].set_title('Unreferenced magnitude spectrum: location only, never quantitative')
    if expected is not None:
        bounds=[max(x[0],expected-.3),min(x[-1],expected+.3)]
        if bounds[0]<bounds[1]:
            candidates=candidate_peaks(x,abs(s),bounds,0.)
            m=e.mask(x,bounds);axs[1].plot(x[m],abs(s[m]),lw=.8)
            axs[1].set_xlim(bounds[1],bounds[0]);axs[1].axvline(expected,color='#b64a54',ls=':')
    else:
        axs[1].text(.05,.5,'No fixed solvent reference suggested. Inspect metadata / clarify solvent.',transform=axs[1].transAxes)
    axs[1].set_title('Reference candidates only: verify identity, multiplet centre and overlap')
    for ax in axs:ax.set_xlabel('Original chemical shift / ppm');ax.set_ylabel('Magnitude / a.u.')
    output.mkdir(parents=True,exist_ok=True);fig.savefig(output/'solvent_survey.png',dpi=160);e.plt.close(fig)
    result={'stage':'survey_no_reference_no_Mn','source_sha256':before,
        'solvent':info,'candidate_peaks_original_ppm':candidates,'processing':meta,
        'acquisition':{k:dic['parameters'].get(k) for k in ['x_freq','digital_filter','filter_factor','relaxation_delay','x_angle']}}
    if before!=e.sha256(src):raise RuntimeError('Input changed during survey.')
    (output/'solvent_survey.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    return result


def inspect(cfg,output):
    src=Path(cfg['input_jdf']);before=e.sha256(src)
    dic,fid=e.ng.jeol.read(str(src))
    x,s,meta=e.transform(dic,fid,cfg)
    phased,ph=e.phase_spectrum(x,s,cfg)
    b,bmeta=e.baseline(x,phased.real,cfg['baseline_regions'])
    y=phased.real-b
    output.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(output/'review_arrays.npz',ppm=x,corrected=y,phased=phased,baseline=b)
    regions=cfg.get('review_search_windows')
    if not regions or set(regions)!={'main','end'}:
        raise ValueError('Specify separate broad main/end review_search_windows for this sample.')
    fig,axes=e.plt.subplots(2,2,figsize=(13,8),constrained_layout=True)
    diagnostics={}
    for row,key in enumerate(['main','end']):
        search=regions[key]; bounds=cfg['main_integral_ppm' if key=='main' else 'end_integral_ppm']
        if search[0]>bounds[0]-.05 or search[1]<bounds[1]+.05:
            raise ValueError('Search must extend at least 0.05 ppm beyond both provisional boundaries.')
        peaks=candidate_peaks(x,y,search,bmeta['median_point_noise'])
        area=e.integrate(x,y,bounds)
        # Vary each boundary independently: catches translation as well as width errors.
        offsets=np.linspace(-.05,.05,11)
        sweep=np.full((11,11),np.nan)
        rows=[]
        for i,dl in enumerate(offsets):
            for j,dh in enumerate(offsets):
                a,z=bounds[0]+dl,bounds[1]+dh
                if a>=z:continue
                value=e.integrate(x,y,[a,z]);sweep[i,j]=value
                rows.append([a,z,value])
        np.savetxt(output/f'{key}_boundary_sweep.csv',rows,delimiter=',',
                   header='lower_ppm,upper_ppm,signed_area',comments='')
        diagnostics[key]={'search_ppm':search,'provisional_bounds_ppm':bounds,
            'candidate_peaks':peaks,'provisional_signed_area':area,
            'boundary_sweep_area_min_max':[float(np.nanmin(sweep)),float(np.nanmax(sweep))]}
        ax=axes[row,0];m=e.mask(x,search)
        ax.plot(x[m],phased.real[m],color='#aaa',lw=.7,label='Phased')
        ax.plot(x[m],b[m],color='#bb8b22',lw=.9,label='Baseline')
        ax.plot(x[m],y[m],color='#174a6e',lw=.9,label='Corrected')
        for lo,hi in cfg['baseline_regions']:
            ax.axvspan(lo,hi,color='#aaa',alpha=.12)
        ax.axvspan(*bounds,color='#6ead9a',alpha=.18)
        if peaks:ax.scatter([p['ppm'] for p in peaks],[p['intensity'] for p in peaks],s=7,color='#ac4751')
        ax.set_xlim(search[1],search[0]);ax.set_title(f'{key}: provisional interval; inspect complete multiplet')
        ax.set_xlabel('ppm');ax.axhline(0,color='#777',lw=.5)
        ax.legend(fontsize=8)
        twin=ax.twinx();cum=cumulative_trapezoid(y[m],x[m],initial=0)
        twin.plot(x[m],cum,color='#a34e5a',alpha=.6,lw=.8);twin.set_ylabel('Cumulative signed area')
        ax=axes[row,1]
        relative=(sweep/area-1)*100 if abs(area)>1e-12 else sweep
        im=ax.imshow(relative,origin='lower',extent=[-.05,.05,-.05,.05],aspect='auto',cmap='coolwarm')
        ax.set_xlabel('Upper-bound shift / ppm');ax.set_ylabel('Lower-bound shift / ppm')
        ax.set_title('Area change (%)' if abs(area)>1e-12 else 'Signed area')
        fig.colorbar(im,ax=ax)
    fig.savefig(output/'interval_review.png',dpi=170);e.plt.close(fig)
    after=e.sha256(src)
    if before!=after:raise RuntimeError('Input changed during inspection.')
    result={'stage':'inspect_only_no_Mn','source_sha256':before,'config_sha256':config_hash(cfg),
        'processing':meta,'phase':ph,'baseline':bmeta,'regions':diagnostics,
        'note':'Candidate peaks and stable areas do not establish chemical identity or absence of overlap.'}
    (output/'diagnostics.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    # Always a draft; never manufacture an accepted review from diagnostics.
    review={'source_sha256':before,'config_sha256':config_hash(cfg),'decision':'needs_review',
        **{k:'' for k in ['solvent_and_reference_assessment','sample_identity_basis','phase_assessment','baseline_assessment',
        'main_assignment_and_bounds','end_assignment_and_bounds','overlap_assessment',
        'boundary_sensitivity_assessment','formula_applicability']}}
    (output/'review_draft.json').write_text(json.dumps(review,ensure_ascii=False,indent=2),encoding='utf-8')
    return result


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('stage',choices=['survey','inspect','quantify'])
    ap.add_argument('--config',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--review',type=Path)
    args=ap.parse_args()
    cfg=json.loads(args.config.read_text(encoding='utf-8-sig'))
    if args.stage=='survey':
        r=survey(cfg,args.out)
        print(json.dumps({'stage':r['stage'],'solvent':r['solvent'],'output':str(args.out)},ensure_ascii=False,indent=2))
    elif args.stage=='inspect':
        r=inspect(cfg,args.out)
        print(json.dumps({'stage':r['stage'],'output':str(args.out),'config_sha256':r['config_sha256']},indent=2))
    else:
        if not args.review:ap.error('quantify requires --review')
        review=json.loads(args.review.read_text(encoding='utf-8-sig'))
        assert_review(cfg,review,e.sha256(cfg['input_jdf']))
        result=e.run(args.config,args.out)
        result['interval_review']=review
        (args.out/'result.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
        (args.out/'review_used.json').write_text(json.dumps(review,ensure_ascii=False,indent=2),encoding='utf-8')


if __name__=='__main__':main()
