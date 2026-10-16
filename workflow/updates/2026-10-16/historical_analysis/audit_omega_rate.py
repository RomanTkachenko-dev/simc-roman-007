"""Read-only audit of the existing omega sample's absolute SIMC normalization.

No simulation rerun, no fit/Pythia scale, no invented efficiency correction.
Output is conditional on the archived model, apparatus and generation domain.
"""
from pathlib import Path
import hashlib
import json
import math
import re
import subprocess
import ROOT

BASE = Path(__file__).resolve().parent
RUN = BASE/'runs/omega_Roman_007'
OUT = BASE/'omega_rate_audit_70uA_1h'
COMMIT = 'b6025c462ef33eeb620bdf4b399cd602c92e89bf'
REPO = '/home/severus/eic/ACPOSTUMA_SIMC_GFORTRAN'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    OUT.mkdir(exist_ok=True)
    assert not (OUT/'summary.json').exists(), 'Existing audit will not be overwritten.'
    paths = [RUN/'omega_Roman_007.root', RUN/'omega_Roman_007.hist',
             RUN/'omega_Roman_007.inp', RUN/'simc.log',
             BASE/'meeting_20260916_overlay/summary.json', BASE/'source_b6025c4.tar.gz']
    hashes = {str(p):sha(p) for p in paths}
    hist = paths[1].read_text()
    def number(pattern):
        match = re.search(pattern, hist, flags=re.MULTILINE)
        assert match, pattern
        return float(match.group(1).replace('D','E'))
    norm = number(r'^\s*normfac\s*=\s*([\d.E+\-]+)')
    input_charge = number(r'^\s*charge\s*=\s*([\d.E+\-]+)\s*mC')
    trials = int(number(r'Ntried\s*=\s*(\d+)'))
    saved = int(number(r'Ncontribute\s*=\s*(\d+)'))
    requested = int(number(r'Ngen \(request\)\s*=\s*(\d+)'))
    lum = number(r'^\s*luminosity\s*=\s*([\d.E+\-]+)\s*ub\^-1')
    genvol = number(r'^\s*genvol\s*=\s*([\d.E+\-]+)')
    internal_yield = number(r'MeV: wtcontr=\s*([\d.E+\-]+)')
    meta = json.loads(paths[4].read_text())
    assert requested == saved == 1000
    assert math.isclose(norm, lum*genvol*saved/trials, rel_tol=5e-6)
    original_code = {}
    for file in ['simc.f','event.f','results_write.f','physics_omega.f','init.f']:
        source = subprocess.check_output(['git','-C',REPO,'show',COMMIT+':'+file], text=True)
        original_code[file] = source
    assert 'normfac = luminosity/ntried*nevent' in original_code['simc.f']
    assert 'normfac = normfac * genvol' in original_code['simc.f']
    assert 'wtcontribute/nevent' in original_code['simc.f']
    assert 'ntu(47) = main%weight' in original_code['results_write.f']
    assert 'main%SF_weight*main%jacobian*main%gen_weight*main%sigcc' in original_code['event.f']
    f = ROOT.TFile.Open(str(paths[0]),'READ'); assert f and not f.IsZombie()
    t = f.Get('h10'); assert t.GetEntries() == saved
    formula = ROOT.TTreeFormula('acceptance_cuts',meta['cuts'],t)
    all_w = []; cut_w = []; legacy_w = []; accepted = []
    # Reproduce the archived simc.f internal pass_cuts literally, including
    # its p-delta upper bound referencing SPedge%e%delta%max (10 here).
    # This diagnostic is NOT used as the requested analysis cut.
    legacy = ROOT.TTreeFormula('legacy_internal_cuts',
        'ssdelta>-10 && ssdelta<10 && abs(ssyptar)<0.05 && abs(ssxptar)<0.10'
        ' && hsdelta>-15 && hsdelta<10 && abs(hsyptar)<0.05 && abs(hsxptar)<0.05', t)
    for i in range(t.GetEntries()):
        t.GetEntry(i); w = float(t.Weight)
        assert math.isfinite(w) and w>0
        all_w.append(w)
        if legacy.EvalInstance(): legacy_w.append(w)
        if formula.EvalInstance():
            cut_w.append(w)
            accepted.append({key:float(getattr(t,key)) for key in ['ssdelta','ssdeltai','missmass','Q2','W','t']})
    del formula, legacy; f.Close()
    assert len(cut_w) == meta['channels']['omega']['after_cuts']
    assert math.isclose(math.fsum(cut_w),meta['channels']['omega']['sum_weight_after'], rel_tol=1e-12)
    reconstructed_internal = norm/saved*math.fsum(legacy_w)
    assert math.isclose(reconstructed_internal,internal_yield,rel_tol=5e-6), (reconstructed_internal,internal_yield)
    charge_mC = 70e-6*3600*1000
    scale = norm/saved*charge_mC/input_charge
    expected = scale*math.fsum(cut_w)
    source_refs = {}
    patterns = {
        'simc.f':['targetfac=targ%mass_amu','luminosity = EXPER%charge/targetfac',
                  'normfac = luminosity/ntried*nevent','normfac = normfac * genvol','wtcontribute/nevent'],
        'event.f':['main%sigcc = peep_omega','main%weight = main%SF_weight'],
        'results_write.f':['ntu(47) = main%weight'],
        'physics_omega.f':['sig4 = sig_LT_extrapolation','peep_omega = Jttheta_fx*gamma_T*ntup%sigcm']}
    for file, search in patterns.items():
        source_refs[file] = [{'line':i,'text':line.strip()} for i,line in enumerate(original_code[file].splitlines(),1)
                             if any(p in line for p in search)]
    report = {
        'status':'Conditional model prediction within the existing generated domain, NOT validated recorded-event rate',
        'beam_current_uA':70,'beam_time_s':3600,'requested_charge_mC':charge_mC,
        'source_commit':COMMIT,'source_references':source_refs,
        'normalization':{'run_charge_mC':input_charge,'normfac_printed':norm,'Ntried':trials,
            'Nsuccess_normalization_denominator':saved,'luminosity_per_run_ub_inverse':lum,
            'genvol_printed':genvol,'event_factor_per_1mC':norm/saved,
            'event_factor_for_requested_charge':scale,
            'formula':'N_model = sum(Weight of stored omega events passing cuts) * normfac/Nsuccess * Q_requested/Q_input',
            'fit_or_Pythia_factors_used':False,'additional_acceptance_factor_used':False,
            'external_efficiency_or_live_time_applied':False},
        'cuts':meta['cuts'],'mass_window_added':False,
        'events':{'stored_before_analysis_cuts':saved,'after_cuts':len(cut_w),
            'sum_weight_before':math.fsum(all_w),'sum_weight_after':math.fsum(cut_w),
            'sum_weight_squared_after':math.fsum(w*w for w in cut_w),
            'effective_weighted_sample_size_after':math.fsum(cut_w)**2/math.fsum(w*w for w in cut_w),
            'negative_weights':0},
        'prediction':{'omega_per_mC_after_cuts':norm/saved*math.fsum(cut_w),
            'omega_per_1h_at_70uA_after_cuts':expected,'omega_per_second':expected/3600,
            'omega_per_1h_at_70uA_before_analysis_cuts':scale*math.fsum(all_w)},
        'setup':{'Ebeam_GeV':10.60,'SHMS_electron_p_GeV':6.7,'SHMS_angle_deg':13,
            'HMS_proton_p_GeV':1.1,'HMS_angle_deg':32,
            'target':'hydrogen, A=Z=1','density_g_cm3':.07332,'areal_density_g_cm2':.7332,
            'target_length_cm':10,'using_rad':False,'using_Eloss':False,'mc_smear':False,'doing_decay':False},
        'coverage':{'generated_electron_delta_percent':[-10.5,10.5],
            'requested_SHMS_analysis_delta_percent':[-10,22],
            'accepted_ssdelta_range':[min(v['ssdelta'] for v in accepted),max(v['ssdelta'] for v in accepted)],
            'accepted_ssdeltai_range':[min(v['ssdeltai'] for v in accepted),max(v['ssdeltai'] for v in accepted)],
            'note':'No prediction for ungenerated region above +10.5%; check whether relevant omega kinematics populate it before claiming complete analysis-acceptance coverage.'},
        'normalization_closure':{'hist_internal_wtcontr':internal_yield,
            'recomputed_same_legacy_internal_cut':reconstructed_internal,'stored_events_in_legacy_cut':len(legacy_w),
            'note':'Internal simc.f pass_cuts are different from requested cuts. The archived code uses electron upper delta edge in the proton upper-edge comparison. hard_cuts=0: this does not discard stored events; do not use hist wtcontr as the requested-cut yield.'},
        'remaining_inputs_or_checks':[
            'Bill confirmation of active omega sig_LT_extrapolation model and its absolute normalization in these kinematics.',
            'Confirmation of actual target thickness/density, beam energy and spectrometer settings (existing input says provisional).',
            'Check generation-domain coverage, especially SHMS +10.5 to +22 percent.',
            'Measured tracking/PID/trigger efficiencies and DAQ live time, target-density change at 70 uA if applicable.',
            'Radiation, energy loss and smearing are disabled; quantify before production-rate claim.',
            'Specify whether an omega mass window or additional kinematic cuts should define the counted yield.',
            'No decay branching factor inserted: simulated channel is e p -> e p omega detected via missing mass, not a selected omega decay. Confirm intended observable.'
        ],
        'sources_sha256':hashes,
        'verification':'Original ROOT event sum matched earlier audit; normfac checked against L*genvol*Nsuccess/Ntried; reproduced hist integrated yield using its exact legacy cuts; source hashes unchanged.'}
    for p,digest in hashes.items(): assert sha(Path(p))==digest
    (OUT/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__': main()
