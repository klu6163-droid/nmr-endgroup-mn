"""JEOL 1D FID -> phased, baseline-corrected spectrum -> end-group Mn.

Entry point: python scripts/workflow.py survey/inspect/quantify --config sample.json --out output
Only the configured output directory is written. Input JDF is read-only.
The Mn result is conditional on the supplied assignments (PDLLA-PEG-PDLLA diol).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import tempfile
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "codex-nmr-mpl"))
import numpy as np
import scipy
from scipy.optimize import differential_evolution, least_squares
from scipy.integrate import cumulative_trapezoid, trapezoid
import nmrglue as ng
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import solvent


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def mask(x, bounds):
    a, b = sorted(bounds)
    return (x >= a) & (x <= b)


def integrate(x, y, bounds):
    """Signed integral, including linearly interpolated exact boundaries."""
    a, b = sorted(bounds)
    if a < x[0] or b > x[-1] or a == b:
        raise ValueError(f"Invalid/out-of-range integration interval: {bounds}")
    inside = (x > a) & (x < b)
    if inside.sum() < 3:
        raise ValueError(f"Too few points in {bounds}")
    xx = np.r_[a, x[inside], b]
    yy = np.r_[np.interp(a, x, y), y[inside], np.interp(b, x, y)]
    return float(trapezoid(yy, xx))


def molecular_weight(ia, ie, peg=400.0, repeat=72.0):
    if not np.isfinite(ia + ie) or ia <= 0 or ie <= 0:
        raise ValueError("Both signed integrals must be finite and positive.")
    # Two terminal lactyl CH + four PEG ester-adjacent CH2 protons = six H.
    return (6 * ia / ie + 2) * repeat + peg


def transform(dic, fid, cfg, drop_override=None):
    h, p = dic["header"], dic["parameters"]
    if h["data_dimension_number"] != 1 or fid.ndim != 1:
        raise ValueError("Only 1D JDF is supported.")
    if h["units"][0][2] != "second" or not np.iscomplexobj(fid):
        raise ValueError("This script expects a complex time-domain FID; no implicit FFT of processed data.")
    if h["units"][0][0] != "none" or h["data_offset_start"][0] != 0:
        raise ValueError("Unsupported scaled/nonzero-start time axis; verify acquisition timing manually.")
    sw, obs = float(p["x_sweep"]), float(p["x_freq"]) / 1e6
    dt_header = (h["data_axis_stop"][0] - h["data_axis_start"][0]) / (len(fid) - 1)
    if not np.isclose(dt_header, 1 / sw, rtol=1e-3):
        raise ValueError("Header time step and acquisition sweep disagree.")
    drop = cfg["fid_drop_points"] if drop_override is None else drop_override
    if drop == "auto":
        # Empirical digital-filter onset, not a manufacturer group-delay constant.
        head = np.abs(fid[:min(256, len(fid))])
        drop = int(np.flatnonzero(head >= .2 * head.max())[0])
    drop = int(drop)
    if drop < 0 or drop >= min(256, len(fid) // 4):
        raise ValueError("Unreasonable initial FID truncation.")
    f = fid[drop:].copy()
    f -= np.mean(fid[-max(32, len(fid) // 20):])
    f *= np.exp(-np.pi * cfg["line_broadening_hz"] * np.arange(len(f)) / sw)
    nfft = 1 << int(np.ceil(np.log2(len(f) * cfg["zero_fill_factor"])))
    spectrum = np.fft.fftshift(np.fft.fft(f, n=nfft))
    # nmrglue JEOL reader returns stored real - i*stored imag.
    # With NumPy's negative-exponent FFT, ppm INCREASES with fftfreq.
    ppm = p["x_offset"] + np.fft.fftshift(np.fft.fftfreq(nfft, 1 / sw)) / obs
    source_hash = sha256(cfg['input_jdf']) if cfg.get('input_jdf') else None
    shift, reference_meta = solvent.reference_shift(ppm,spectrum,p,cfg,source_hash)
    ppm += shift
    return ppm, spectrum, {"fid_drop_points": drop, "reference_shift_ppm": shift,
        "nfft": nfft, "sweep_hz": sw, "obs_mhz": obs,
        "delay_method": "empirical onset plus fitted residual first-order phase",
        "reference_review": reference_meta, "region_conflicts":solvent.region_conflicts(cfg)}


def phase_spectrum(x, s, cfg):
    pivot = cfg["phase_pivot_ppm"]
    if cfg.get("phase_anchors") and cfg.get("manual_phase") is None:
        # Complex Lorentzian fits separate absorption/dispersion while allowing
        # a local complex linear baseline. Anchors need not be analyte peaks.
        anchors = []
        for spec in cfg["phase_anchors"]:
            center, half = spec["center_ppm"], spec["half_width_ppm"]
            m = mask(x, [center-half, center+half])
            xx, zz = x[m], s[m]
            if len(xx) < 12: raise ValueError("Too few phase-anchor data points.")
            zz = zz / np.max(abs(zz))
            offsets = np.array(spec.get("offsets_ppm", [0.0]))
            weights = np.array(spec.get("relative_weights", [1.0]))

            def local_fit(par, full=False):
                c, g = par
                line = sum(w*g/(g+1j*(xx-c-o)) for o,w in zip(offsets,weights))
                design = np.column_stack([line, np.ones(len(xx)), xx-center])
                coef = np.linalg.lstsq(design, zz, rcond=None)[0]
                residual = design @ coef - zz
                return (coef, residual) if full else np.r_[residual.real,residual.imag]

            fit = least_squares(local_fit, [center,.001],
                bounds=([center-half/2,.0001],[center+half/2,min(.01,half)]))
            coef, residual = local_fit(fit.x, True)
            rms = float(np.sqrt(np.mean(abs(residual)**2)))
            if not fit.success or rms > .08:
                raise RuntimeError(f"Phase-anchor fit failed at {center} ppm (RMS={rms:.3f}).")
            anchors.append({"center_ppm":float(fit.x[0]),"halfwidth_ppm":float(fit.x[1]),
                "observed_phase_degrees":float(np.angle(coef[0],deg=True)),"normalized_complex_rms":rms})
        anchors.sort(key=lambda a:a["center_ppm"])
        positions = np.array([a["center_ppm"] for a in anchors])
        angles = np.rad2deg(np.unwrap(np.deg2rad([a["observed_phase_degrees"] for a in anchors])))
        design = np.column_stack([np.ones(len(positions)),positions-pivot])
        ph = np.linalg.lstsq(design,-angles,rcond=None)[0]
        residual_deg = angles+design@ph
        if np.max(abs(residual_deg)) > 3:
            raise RuntimeError("Phase anchors disagree by >3 degrees; review anchor assignments/line shapes.")
        y = s*np.exp(1j*np.deg2rad(ph[0]+ph[1]*(x-pivot)))
        return y, {"p0_degrees_at_pivot":float(ph[0]),"p1_degrees_per_ppm":float(ph[1]),
            "pivot_ppm":pivot,"method":"complex line-shape anchors plus linear phase regression",
            "anchors":anchors,"anchor_phase_residual_degrees":residual_deg.tolist()}
    regions = []
    for bounds in cfg["phase_regions"]:
        m = mask(x, bounds)
        if m.sum() < 10:
            raise ValueError("Phase region is too small.")
        regions.append((x[m]-pivot, s[m] / np.max(abs(s[m]))))

    def score(p):
        total = 0.0
        for xx, zz in regions:
            v = (zz * np.exp(1j * np.deg2rad(p[0] + p[1] * xx))).real
            dv = abs(np.diff(v))
            prob = dv / max(dv.sum(), 1e-30)
            entropy = -np.sum(prob * np.log(prob + 1e-30))
            total += entropy + 1000 * np.mean(np.minimum(v, 0)**2)
        return float(total)

    if cfg.get("manual_phase") is not None:
        ph = np.asarray(cfg["manual_phase"], dtype=float)
        success = True
    else:
        fit = differential_evolution(score, [(-180, 180), (-60, 60)],
                                     seed=42, tol=1e-8, maxiter=300, polish=True)
        ph, success = fit.x, bool(fit.success)
    if not success:
        raise RuntimeError("Automatic phasing failed to converge; inspect/manual override required.")
    y = s * np.exp(1j * np.deg2rad(ph[0] + ph[1] * (x-pivot)))
    return y, {"p0_degrees_at_pivot": float(ph[0]),
        "p1_degrees_per_ppm": float(ph[1]), "pivot_ppm": pivot,
        "objective": score(ph), "method": "region-balanced entropy and negative-area penalty"}


def baseline(x, y, intervals, method="linear"):
    """Median noise anchors; interpolation across excluded signal bands."""
    xp, yp, noise = [], [], []
    for bounds in intervals:
        m = mask(x, bounds)
        if m.sum() < 10:
            raise ValueError(f"Insufficient baseline anchor points: {bounds}")
        v = y[m]
        xp.append(float(np.median(x[m])))
        yp.append(float(np.median(v)))
        noise.append(float(np.median(abs(np.diff(v))) / .9538725524))
    order = np.argsort(xp)
    xp, yp = np.asarray(xp)[order], np.asarray(yp)[order]
    if method == "pchip":
        from scipy.interpolate import PchipInterpolator
        b = PchipInterpolator(xp, yp)(np.clip(x, xp[0], xp[-1]))
    else:
        b = np.interp(x, xp, yp)
    return b, {"anchor_ppm": xp.tolist(), "anchor_intensity": yp.tolist(),
               "median_point_noise": float(np.median(noise)), "method": method}


def calculate(x, y, cfg, endpoint=None, main=None):
    aa = cfg["main_integral_ppm"] if main is None else main
    ee = cfg["end_integral_ppm"] if endpoint is None else endpoint
    ia, ie = integrate(x, y, aa), integrate(x, y, ee)
    value = molecular_weight(ia, ie, cfg["peg_mn"], cfg["repeat_mass"])
    return {"main_bounds_ppm": aa, "end_bounds_ppm": ee, "I_main": ia,
        "I_end": ie, "ratio_main_to_end": ia/ie, "Mn_g_mol": value,
        "lactyl_units_total": 6*ia/ie+2}


def sensitivity(dic, fid, x, s, phased, cfg, meta, primary):
    rows = []
    def add(label, xx, yy, **kw):
        try:
            rows.append({"variation": label, **calculate(xx, yy, cfg, **kw)})
        except ValueError as e:
            rows.append({"variation": label, "error": str(e)})
    a, e = cfg["main_integral_ppm"], cfg["end_integral_ppm"]
    for d in [-.02, 0, .02]:
        add(f"end_bounds_expand_{d:+.2f}_ppm", x, primary, endpoint=[e[0]-d, e[1]+d])
    for d in [-.02, .02]:
        add(f"main_bounds_expand_{d:+.2f}_ppm", x, primary, main=[a[0]-d, a[1]+d])
    for degrees in [-2, 2]:
        yy = (phased*np.exp(1j*np.deg2rad(degrees))).real
        b, _ = baseline(x, yy, cfg["baseline_regions"])
        add(f"phase_p0_delta_{degrees:+d}_deg", x, yy-b)
    b, _ = baseline(x, phased.real, cfg["baseline_regions"], "pchip")
    add("baseline_pchip", x, phased.real-b)
    for delta in [-1, 1]:
        drop = meta["fid_drop_points"]+delta
        if drop < 0:
            continue
        xx, ss, _ = transform(dic, fid, cfg, drop_override=drop)
        pp, _ = phase_spectrum(xx, ss, cfg)
        bb, _ = baseline(xx, pp.real, cfg["baseline_regions"])
        add(f"fid_drop_{drop}_refit_phase", xx, pp.real-bb)
    return rows


def plot_result(x, s, phased, b, corrected, cfg, result, output):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, axs = plt.subplots(4, 1, figsize=(12, 12), constrained_layout=True)
    axs[0].plot(x, s.real, color="#adb5bd", lw=.6, label="FFT real, before phase")
    axs[0].plot(x, corrected, color="#174a6e", lw=.7, label="Phase + baseline corrected")
    axs[0].set_xlim(6.2, 0)
    m = mask(x, [0, 6.2]); axs[0].set_ylim(np.min(s.real[m])*1.1, max(corrected[m])*1.1)
    axs[0].legend(loc="upper left", fontsize=9)
    axs[0].set_title(f"PDLLA-PEG-PDLLA | provisional end-group Mn = {result['Mn_g_mol']:,.0f} g/mol")
    for ax, bounds, lim, title in zip(axs[1:3],
        [cfg["main_integral_ppm"], cfg["end_integral_ppm"]],
        [(v[1]+.25,v[0]-.25) for v in [cfg["main_integral_ppm"],cfg["end_integral_ppm"]]], ["Main lactyl CH (assignment assumed)", "Terminal CH + PEG ester-adjacent CH2 (assignment assumed)"]):
        ax.plot(x, phased.real, color="#aaa", lw=.8, label="After phase")
        ax.plot(x, b, color="#c67f12", lw=1, label="Baseline")
        ax.plot(x, corrected, color="#174a6e", lw=1, label="Corrected")
        ax.axvspan(*bounds, color="#70b7a3", alpha=.18)
        for v in bounds: ax.axvline(v, color="#469b83", ls=":", lw=.8)
        mm = mask(x, sorted(lim)); ax.set_xlim(*lim)
        span = max(corrected[mm])-min(corrected[mm]); ax.set_ylim(min(corrected[mm])-.07*span,max(corrected[mm])+.1*span)
        ax.axhline(0,color="#777",lw=.5); ax.set_title(title)
        integ = cumulative_trapezoid(corrected[mask(x,bounds)], x[mask(x,bounds)], initial=0)
        twin = ax.twinx(); twin.plot(x[mask(x,bounds)],integ/result['I_end'],color="#b54853",lw=1)
        twin.set_ylabel("Cumulative integral / I(end)",color="#b54853")
    axs[1].legend(fontsize=8,loc="upper left")
    axs[3].plot(x, corrected, color="#174a6e",lw=.8)
    axs[3].axvspan(*cfg["end_integral_ppm"],alpha=.18,color="#70b7a3",label="Sample-specific endpoint interval")
    
    for lo,hi in cfg.get("comparison_end_intervals", []):
        axs[3].axvspan(lo,hi,alpha=.18,color="#e8a460",label=f"Comparison: {lo}-{hi} ppm")
    for lo, hi in cfg["baseline_regions"]:
        axs[3].axvspan(lo, hi,color="#b9bfc4",alpha=.12)
    axs[3].set_xlim(cfg["main_integral_ppm"][1]+.3,cfg["end_integral_ppm"][0]-.25); mm=mask(x,cfg["end_integral_ppm"]); top=max(corrected[mm])*1.25
    axs[3].set_ylim(-top*.15,top); axs[3].axhline(0,color="#777",lw=.5)
    axs[3].set_title("Weak-signal / baseline check (main peak clipped for display)")
    axs[3].legend(fontsize=8)
    for ax in axs:
        ax.set_xlabel("Chemical shift / ppm"); ax.set_ylabel("Intensity / a.u.")
    fig.savefig(output / "spectrum_review.png",dpi=180)
    fig.savefig(output / "spectrum_review.pdf")
    plt.close(fig)


def run(config_path, output_override=None):
    cfg = json.loads(Path(config_path).read_text(encoding="utf-8-sig"))
    src = Path(cfg["input_jdf"]).resolve()
    out = Path(output_override or cfg["output_dir"])
    if not out.is_absolute(): out = Path(config_path).resolve().parent / out
    out = out.resolve()
    if out == src.parent:
        raise ValueError("Use a separate output folder.")
    out.mkdir(parents=True, exist_ok=True)
    before = sha256(src)
    dic, fid = ng.jeol.read(str(src))
    x, s, meta = transform(dic, fid, cfg)
    phased, phase_meta = phase_spectrum(x, s, cfg)
    b, baseline_meta = baseline(x, phased.real, cfg["baseline_regions"])
    y = phased.real-b
    primary = calculate(x, y, cfg)
    variants = sensitivity(dic, fid, x, s, phased, cfg, meta, y)
    literal = []
    for bounds in cfg.get("comparison_end_intervals", []):
        try:
            literal.append({"status": "comparison only: chemical assignment and coverage require review",
                            **calculate(x,y,cfg,endpoint=bounds)})
        except ValueError as ex:
            literal.append({"end_bounds_ppm":bounds,"error":str(ex)})
    p = dic["parameters"]
    warnings = [
        "PRELIMINARY: assumes the selected end-group region contains exactly six non-exchangeable H per chain; assignments/purity require review.",
        "Integration intervals and assignments must be reviewed for each sample.",
        "FID truncation is empirically estimated; it is not a verified JEOL digital-filter delay constant.",
        "Relaxation delay and pulse angle alone do not establish quantitative response without T1/validation data.",
        "Sensitivity range is parameter sensitivity, NOT a statistical confidence interval.",
        "No statistical confidence interval or absolute-Mn validation is implied."]
    after = sha256(src)
    if before != after: raise RuntimeError("Source hash changed during processing.")
    mn_values = [v["Mn_g_mol"] for v in variants if "Mn_g_mol" in v]
    result = {"status":"provisional_assignment_based_estimate", "source":str(src),
        "source_sha256_before":before,"source_sha256_after":after,
        "formula":"Mn = (6 * I_main / I_end + 2) * repeat_mass + PEG_Mn",
        "sample_identity":cfg.get("sample_identity", "unconfirmed"),
        "acquisition":{k:p.get(k) for k in ["solvent","x_freq","scans","x_points","x_angle","relaxation_delay","repetition_time","digital_filter","filter_factor"]},
        "processing":{**meta, **phase_meta, "baseline":baseline_meta},
        "primary":primary, "literal_image_intervals":literal, "sensitivity":variants,
        "sensitivity_min_max_g_mol":[min(mn_values),max(mn_values)],"warnings":warnings,
        "versions":{"python":sys.version,"numpy":np.__version__,"scipy":scipy.__version__,"nmrglue":ng.__version__}}
    np.savetxt(out/"spectrum.csv",np.column_stack([x,s.real,s.imag,phased.real,phased.imag,b,y]),delimiter=",",
        header="ppm,fft_real,fft_imag,phased_real,phased_imag,baseline,corrected_real",comments="")
    np.savez_compressed(out/"processing_arrays.npz",raw_fid=fid,ppm=x,fft=s,phased=phased,baseline=b,corrected=y)
    (out/"result.json").write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
    (out/"config_used.json").write_text(json.dumps(cfg,ensure_ascii=False,indent=2),encoding="utf-8")
    plot_result(x,s,phased,b,y,cfg,primary,out)
    report = ["# 核磁端基法分子量处理结果", "", "**状态：按暂定峰归属计算的初步估算，不是已验证的绝对分子量。**", "",
        f"- 文件：`{src.name}`；原文件 SHA-256 前后相同。",
        f"- 溶剂：{p.get('solvent')}；频率：{meta['obs_mhz']:.6f} MHz；扫描：{p.get('scans')}。",
        f"- 弛豫延迟：{p.get('relaxation_delay')} s；脉冲角：{p.get('x_angle')}°；没有 T1 数据。",
        f"- FID 起始移除 {meta['fid_drop_points']} 点（经验估计），随后 0/1 阶相位优化、噪声锚点基线扣除。",
        f"- 化学位移整体平移 {meta['reference_shift_ppm']:+.6f} ppm；参照设置：{cfg.get('reference')}。", "",
        "## 计算", "", f"Mn = (6 × I主链 / I端基 + 2) × {cfg['repeat_mass']} + {cfg['peg_mn']}", "",
        "分母按两个末端乳酸 CH（2H）与两侧 PEG 酯相邻 CH2（4H）共 6H 处理；+2 补回末端乳酸单元。",
        f"- 主链积分区间：{cfg['main_integral_ppm']} ppm；积分 {primary['I_main']:.8g}。",
        f"- 端基积分区间：{cfg['end_integral_ppm']} ppm；积分 {primary['I_end']:.8g}。",
        f"- 端基积分归一为 1：主链积分 {primary['ratio_main_to_end']:.6f}。",
        f"- **暂定 Mn = {primary['Mn_g_mol']:.2f} g/mol**。", "",
        "区间是根据本次谱线峰群选定的，可在配置中修改；峰归属仍须结合结构/其他谱学证据确认。", "",
        "## 参数敏感性（不是置信区间）", "", "|变化|Mn / g mol⁻¹|", "|---|---:|"]
    for row in variants: report.append(f"|{row['variation']}|{row.get('Mn_g_mol',float('nan')):.2f}|")
    report += ["", "## 配置中的对照区间", "", "以下只作区间对照，是否适用须结合峰归属和覆盖范围判断。"]
    for row in literal: report.append(f"- {row['end_bounds_ppm']} ppm：{row.get('Mn_g_mol', '积分无效')} g/mol。")
    report += ["", "本次区间选择依据见 review_used.json。", "",
        "## 复核与方法来源", "", "所有积分为有符号实部面积，不取绝对值、不截去负值，不用幅度谱做定量。",
        "查看 spectrum_review.png 中峰脚、基线及端基区；相位参数 manual_phase 为 [枢轴处角度, 每 ppm 角度斜率]。",
        "- [nmrglue JEOL 读取器](https://github.com/jjhelmus/nmrglue/blob/master/nmrglue/fileio/jeol.py)",
        "- [DMSO-d6 2.50 ppm 参照](https://chem.vt.edu/analytical-services/nmr.html)",
        "- [JEOL 定量测量与 T1](https://www.jeol.com/solutions/applications/details/qNMR-issue8.php)", ""]
    (out/"report.md").write_text("\n".join(report),encoding="utf-8")
    print(json.dumps({"output":str(out),"primary":primary,"sensitivity_min_max":result['sensitivity_min_max_g_mol'],"phase":phase_meta,"source_unchanged":before==after},ensure_ascii=False,indent=2))
    return result


if __name__ == "__main__":
    raise SystemExit("Use workflow.py inspect / quantify with a per-sample review record.")
