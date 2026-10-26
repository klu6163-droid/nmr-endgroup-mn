"""Metadata-first solvent handling. Reference positions are suggestions only."""
import re
import numpy as np

SOURCE = 'https://b2b.sigmaaldrich.com/BR/en/technical-documents/technical-article/analytical-chemistry/nuclear-magnetic-resonance/1h-nmr-and-13c-nmr-chemical-shifts-of-impurities-chart'
POSITIONS = {'DMSO-d6':2.50,'CDCl3':7.26,'acetone-d6':2.05,
             'CD3CN':1.94,'CD3OD':3.31,'C6D6':7.16}
ALIASES = {
    'DMSOD6':'DMSO-d6','D6DMSO':'DMSO-d6','DIMETHYLSULFOXIDED6':'DMSO-d6',
    'CDCL3':'CDCl3','CHLOROFORMD':'CDCl3','CHLOROFORMD1':'CDCl3',
    'ACETONED6':'acetone-d6','D6ACETONE':'acetone-d6',
    'CD3CN':'CD3CN','ACETONITRILED3':'CD3CN',
    'CD3OD':'CD3OD','METHANOLD4':'CD3OD','MEOD':'CD3OD',
    'C6D6':'C6D6','BENZENED6':'C6D6','D2O':'D2O','DEUTERIUMOXIDE':'D2O'}


def normalize(value):
    if not isinstance(value,str) or not value.strip():return None
    key=value.translate(str.maketrans('₀₁₂₃₄₅₆₇₈₉','0123456789')).upper()
    key=re.sub(r'[\s_\-()]','',key)
    return ALIASES.get(key)


def describe(parameters,cfg):
    raw=parameters.get('solvent')
    canonical=normalize(raw)
    provenance='file_metadata'
    override=cfg.get('solvent_override')
    if override:
        if len(override.get('basis','').strip())<12:
            raise ValueError('Solvent override needs a documented user/experimental basis.')
        canonical=normalize(override.get('name'));provenance='documented_override'
    raw_text=str(raw or '')
    status='recognized' if canonical else ('missing' if not raw_text.strip() else
        ('mixed_or_ambiguous' if any(t in raw_text for t in ['/', '+', ',', ';', ':']) else 'unrecognized'))
    expected=POSITIONS.get(canonical)
    return {'raw_metadata':raw,'canonical':canonical,'status':status,
        'provenance':provenance,'override':override,
        'suggested_residual_ppm':expected,'reference_table':SOURCE,
        'note':'Metadata describes the acquisition setting, not independently verified sample contents. '
               'Unknown/mixed solvents need clarification; D2O/HOD is not assigned a fixed reference here. '
               'Water signals and analyte shifts must be inspected for this sample.'}


def reference_shift(ppm,spectrum,parameters,cfg,source_hash):
    info=describe(parameters,cfg)
    ref=cfg.get('reference')
    if not ref or ref.get('mode')=='none':
        return 0.,{'status':'original_axis_preserved','solvent':info}
    if ref.get('mode')!='reviewed_peak':
        raise ValueError('Reference requires mode=reviewed_peak and observed_ppm; use survey first. No automatic largest-peak calibration.')
    if ref.get('source_sha256')!=source_hash:
        raise ValueError('Reference review is missing or belongs to another JDF.')
    if len(ref.get('basis','').strip())<12:
        raise ValueError('Document the reference identity/overlap assessment before shifting.')
    if ref.get('kind','solvent')=='solvent':
        if not info['canonical'] or normalize(ref.get('solvent'))!=info['canonical']:
            raise ValueError('Reference solvent is unknown or conflicts with file/confirmed solvent.')
        if info['canonical']=='D2O':
            raise ValueError('Do not fix HOD to a universal shift; use a justified external/internal standard.')
    elif ref.get('kind')!='standard':
        raise ValueError('Reference kind must be solvent or standard.')
    observed=float(ref['observed_ppm']);target=float(ref['target_ppm'])
    if not np.isfinite(observed+target) or not ppm[0]<observed<ppm[-1]:
        raise ValueError('Invalid observed reference position.')
    delta=target-observed
    if abs(delta)>float(ref.get('max_shift_ppm',.25)):
        raise ValueError('Reference shift exceeds configured limit; inspect rather than forcing calibration.')
    return delta,{'status':'reviewed_reference_applied','observed_ppm':observed,
        'target_ppm':target,'basis':ref['basis'],'solvent':info}


def region_conflicts(cfg):
    """Flag, never silently notch out solvent/water or trim analyte integrals."""
    contaminants=list(cfg.get('observed_contaminant_regions',[]))
    ref=cfg.get('reference')
    if ref and ref.get('mode')=='reviewed_peak' and ref.get('kind','solvent')=='solvent':
        c=ref['target_ppm'];half=ref.get('exclusion_half_width_ppm',.05)
        contaminants.append({'label':'reviewed solvent residual','bounds_ppm':[c-half,c+half]})
    warnings=[]
    for item in contaminants:
        a,b=sorted(item['bounds_ppm'])
        for name,regions in [('baseline',cfg.get('baseline_regions',[])),
            ('main_integral',[cfg['main_integral_ppm']] if 'main_integral_ppm' in cfg else []),
            ('end_integral',[cfg['end_integral_ppm']] if 'end_integral_ppm' in cfg else [])]:
            for lo,hi in regions:
                if max(a,lo)<min(b,hi):warnings.append({'region_type':name,'region':[lo,hi],'contaminant':item})
    return warnings
