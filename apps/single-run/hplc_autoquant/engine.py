"""Deterministic calculations on verified, pre-integrated peak areas.

Concentrations are mg/mL; masses mg; volumes mL. RSDs and recoveries
are expressed in percent units (e.g. 2.0 means 2%, not 200%).
"""
from collections import defaultdict
import math
from statistics import mean, stdev

import numpy as np


class ValidationError(ValueError):
    pass


DEFAULT_CRITERIA = {
    "min_r_squared": 0.995,
    "max_standard_rsd_pct": 2.0,
    "max_backcalc_error_pct": 5.0,
    "max_response_factor_deviation_pct": 10.0,
    "max_preparation_rsd_pct": 2.0,
    "max_replicate_rsd_pct": 2.0,
    "qc_recovery_min_pct": 95.0,
    "qc_recovery_max_pct": 105.0,
    "max_blank_area": 0.0,
}


def number(value, field, *, positive=False, nonnegative=False):
    if isinstance(value, bool):
        raise ValidationError(f"{field}: a finite number is required")
    try:
        result = float(value)
    except (TypeError, ValueError):
        raise ValidationError(f"{field}: a finite number is required") from None
    if not math.isfinite(result):
        raise ValidationError(f"{field}: a finite number is required")
    if positive and result <= 0 or nonnegative and result < 0:
        raise ValidationError(f"{field}: invalid value {result}")
    return result


def stats(values):
    values = list(values)
    avg = mean(values)
    sd = stdev(values) if len(values) > 1 else None
    return {"n": len(values), "mean": avg, "sd": sd,
            "rsd_pct": sd / abs(avg) * 100 if sd is not None and avg != 0 else None}


def analyze(injections, preparations, criteria=None, *, allow_two_levels=False,
            calibration_model='ols', reference_concentration_mg_ml=None):
    if calibration_model not in ('ols', 'single_standard'):
        raise ValidationError('Unknown calibration model')
    limits = dict(DEFAULT_CRITERIA)
    if criteria is not None:
        unknown = set(criteria) - set(limits)
        if unknown:
            raise ValidationError(f"Unknown acceptance criteria: {sorted(unknown)}")
        limits.update(criteria)
    for key, value in limits.items():
        limits[key] = number(value, key, nonnegative=True)
    if limits["min_r_squared"] > 1:
        raise ValidationError("min_r_squared must be <= 1")
    if limits["qc_recovery_min_pct"] > limits["qc_recovery_max_pct"]:
        raise ValidationError("QC recovery limits are reversed")

    rows, seen, flags = [], set(), []

    def flag(code, subject, message):
        flags.append({"id": f"{code}:{subject}", "code": code,
                      "subject": subject, "message": message})

    allowed = {"standard", "sample", "blank", "qc", "system_suitability"}
    standards, samples = defaultdict(list), defaultdict(list)
    preps = {}
    for original in preparations:
        p = dict(original)
        pid = str(p.get("preparation_id", "")).strip()
        if not pid or pid in preps:
            raise ValidationError(f"Missing or duplicate preparation_id: {pid}")
        if not str(p.get("sample_id", "")).strip():
            raise ValidationError(f"{pid}: sample_id is required")
        for key in ("formulation_mass_mg", "final_volume_ml", "additional_dilution_factor"):
            p[key] = number(p.get(key), f"{pid}.{key}", positive=True)
        if p["additional_dilution_factor"] < 1:
            raise ValidationError(f"{pid}: additional dilution factor must be >= 1")
        preps[pid] = p
    for original in injections:
        row = dict(original)
        iid = str(row.get("injection_id", "")).strip()
        if not iid or iid in seen:
            raise ValidationError(f"Missing or duplicate injection_id: {iid}")
        seen.add(iid)
        kind = row.get("kind")
        if kind not in allowed:
            raise ValidationError(f"{iid}: unconfirmed/unsupported kind {kind!r}")
        row["area"] = number(row.get("area"), f"{iid}.area", nonnegative=True)
        if kind in {"standard", "qc"}:
            c = number(row.get("nominal_concentration_mg_ml"), f"{iid}.nominal_concentration_mg_ml", positive=True)
            purity = number(row.get("purity_fraction"), f"{iid}.purity_fraction", positive=True)
            if purity > 1:
                raise ValidationError(f"{iid}: purity_fraction must be in (0, 1]")
            row["nominal_concentration_mg_ml"] = c
            row["purity_fraction"] = purity
            row["effective_concentration_mg_ml"] = c * purity
        if kind == "standard":
            standards[row["effective_concentration_mg_ml"]].append(row)
        if kind == "sample":
            pid = row.get("preparation_id")
            if pid not in preps:
                raise ValidationError(f"{iid}: missing preparation metadata for {pid!r}")
            samples[pid].append(row)
        if kind == "blank" and row["area"] > limits["max_blank_area"]:
            flag("BLANK_RESPONSE", iid, "Blank area exceeds configured limit")
        if kind == "system_suitability":
            flag("SST_NOT_EVALUATED", iid, "Phase 1 does not evaluate chromatographic system suitability")
        rows.append(row)
    if not isinstance(allow_two_levels, bool):
        raise ValidationError('allow_two_levels must be an explicit boolean')
    if calibration_model == 'single_standard':
        reference_concentration_mg_ml = number(reference_concentration_mg_ml, 'STD 1 reference concentration', positive=True)
        matches = [c for c in standards if math.isclose(c, reference_concentration_mg_ml, rel_tol=1e-10)]
        if len(matches) != 1:
            raise ValidationError('Select exactly one included effective STD 1 concentration')
        reference_concentration_mg_ml = matches[0]
        flag('SINGLE_POINT_CALIBRATION', 'calibration', 'Single external standard ratio. Other levels are checks only; linearity and a validated concentration range are not established.')
    elif len(standards) == 2 and allow_two_levels:
        flag('TWO_LEVEL_CALIBRATION', 'calibration', 'Explicit two-level exploratory calibration; exact fit does not establish linearity')
    elif len(standards) < 3:
        raise ValidationError("At least three distinct effective standard concentrations are required")
    if not samples:
        raise ValidationError("At least one sample preparation with injections is required")
    unused = set(preps) - set(samples)
    if unused:
        raise ValidationError(f"Preparations have no sample injections: {sorted(unused)}")
    levels = []
    for c, group in sorted(standards.items()):
        s = stats(r["area"] for r in group)
        levels.append({"concentration_mg_ml": c, **s,
                       "injection_ids": [r["injection_id"] for r in group]})
        if s["n"] == 1:
            flag("NO_STANDARD_REPLICATES", str(c), "Standard repeatability cannot be estimated")
        elif s["rsd_pct"] is None or s["rsd_pct"] > limits["max_standard_rsd_pct"]:
            flag("STANDARD_RSD", str(c), "Standard area RSD exceeds limit or is undefined")
    x = np.array([v["concentration_mg_ml"] for v in levels])
    y = np.array([v["mean"] for v in levels])
    # Each concentration mean has equal weight, regardless of replicate count.
    if calibration_model == 'single_standard':
        reference = next(v for v in levels if v['concentration_mg_ml'] == reference_concentration_mg_ml)
        slope, intercept = reference['mean'] / reference_concentration_mg_ml, 0.0
        model = 'Single external standard: mean STD 1 area / effective STD 1 concentration; sample concentration = sample area / mean STD 1 area x effective STD 1 concentration'
        for v in levels:
            v['assay_reference'] = v['concentration_mg_ml'] == reference_concentration_mg_ml
    else:
        coefficients, _, rank, _ = np.linalg.lstsq(
            np.column_stack((x, np.ones_like(x))), y, rcond=None)
        if rank != 2:
            raise ValidationError("Calibration concentrations do not support a numerically independent slope and intercept")
        slope, intercept = (float(v) for v in coefficients)
        model = 'OLS on concentration mean areas, free intercept, unweighted'
    if slope <= 0 or not math.isfinite(slope) or not math.isfinite(intercept):
        raise ValidationError("Calibration must have a finite positive slope")
    residuals = y - (slope * x + intercept)
    total = float(np.sum((y - np.mean(y)) ** 2))
    if total <= 0 and calibration_model == 'ols':
        raise ValidationError("Calibration response has no variance")
    r2 = 1 - float(np.sum(residuals ** 2)) / total if calibration_model == 'ols' else None
    if r2 is not None and r2 < limits["min_r_squared"]:
        flag("CALIBRATION_R2", "calibration", "R squared is below the configured limit")
    avg_rf = mean(v["mean"] / v["concentration_mg_ml"] for v in levels)
    for v, residual in zip(levels, residuals):
        c = v["concentration_mg_ml"]
        back = (v["mean"] - intercept) / slope
        v.update(residual_area=float(residual), response_factor=v["mean"] / c,
                 backcalc_mg_ml=back, recovery_pct=100 * back / c,
                 error_pct=100 * (back - c) / c)
        if abs(v["error_pct"]) > limits["max_backcalc_error_pct"]:
            flag("STANDARD_BACKCALC", str(c), "Mean standard back-calculation error exceeds limit")
        if abs(v["response_factor"] / avg_rf - 1) * 100 > limits["max_response_factor_deviation_pct"]:
            flag("RESPONSE_FACTOR", str(c), "Response factor differs from mean of level response factors")
    for row in rows:
        kind, iid = row["kind"], row["injection_id"]
        if kind in {"sample", "standard", "qc"}:
            c = (row["area"] - intercept) / slope
            row["calculated_concentration_mg_ml"] = c
            if c < 0:
                flag("NEGATIVE_CONCENTRATION", iid, "Negative calculated concentration; no zero clamping applied")
            if calibration_model == 'ols' and kind in {"sample", "qc"} and not min(x) <= c <= max(x):
                flag("EXTRAPOLATION", iid, "Calculated concentration is outside the calibration range")
            if kind in {"standard", "qc"}:
                recovery = c / row["effective_concentration_mg_ml"] * 100
                row["recovery_pct"] = recovery
                row["error_pct"] = recovery - 100
                row["response_factor"] = row["area"] / row["effective_concentration_mg_ml"]
                if kind == "standard" and abs(recovery - 100) > limits["max_backcalc_error_pct"]:
                    flag("INJECTION_BACKCALC", iid, "Individual standard back-calculation error exceeds limit")
                if kind == "qc" and not limits["qc_recovery_min_pct"] <= recovery <= limits["qc_recovery_max_pct"]:
                    flag("QC_RECOVERY", iid, "QC recovery is outside configured limits")
    results = []
    for pid, group in sorted(samples.items()):
        p = preps[pid]
        s = stats(r["area"] for r in group)
        vial = (s["mean"] - intercept) / slope
        mass = vial * p["additional_dilution_factor"] * p["final_volume_ml"]
        results.append({**p, "n": s["n"], "mean_area": s["mean"],
                        "area_sd": s["sd"], "area_rsd_pct": s["rsd_pct"],
                        "vial_concentration_mg_ml": vial,
                        "preparation_concentration_mg_ml": vial * p["additional_dilution_factor"],
                        "cea_mass_mg": mass,
                        "assay_pct_ww": 100 * mass / p["formulation_mass_mg"],
                        "injection_ids": [r["injection_id"] for r in group]})
        if s["n"] == 1:
            flag("NO_SAMPLE_REPLICATES", pid, "Injection repeatability cannot be estimated")
        elif s["rsd_pct"] is None or s["rsd_pct"] > limits["max_replicate_rsd_pct"]:
            flag("SAMPLE_RSD", pid, "Sample injection area RSD exceeds limit or is undefined")
    # Never average unrelated samples into a single assay.
    by_sample = defaultdict(list)
    for p in results:
        by_sample[p["sample_id"]].append(p["assay_pct_ww"])
    summaries = []
    for sid, values in sorted(by_sample.items()):
        s = stats(values)
        summaries.append({"sample_id": sid, **s})
        if s["n"] == 1:
            flag("NO_PREPARATION_REPLICATES", sid, "Between-preparation precision cannot be estimated")
        elif s["rsd_pct"] is None or s["rsd_pct"] > limits["max_preparation_rsd_pct"]:
            flag("PREPARATION_RSD", sid, "Between-preparation assay RSD exceeds limit or is undefined")
    return {"model": model, 'calibration_model': calibration_model,
            'reference_concentration_mg_ml': reference_concentration_mg_ml,
            "slope": slope, "intercept": intercept, "r_squared": r2,
            "calibration_min_mg_ml": float(min(x)), "calibration_max_mg_ml": float(max(x)),
            "criteria": limits, "standards": levels, "injections": rows,
            "preparations": results, "sample_summaries": summaries, "flags": flags}
