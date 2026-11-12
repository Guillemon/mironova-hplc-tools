"""Read-only input snapshots, content-bound review, and audit records."""
import csv
from datetime import datetime, timezone
import hashlib
import io
from importlib.metadata import version
import json
from pathlib import Path
import platform

from . import __version__
from .engine import analyze, ValidationError


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def digest(data):
    return hashlib.sha256(data).hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def json_write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def load_dataset(folder):
    folder = Path(folder).resolve()
    snapshots, sources = {}, []
    for name in ("injections.csv", "preparations.csv", "config.json"):
        path = folder / name
        try:
            raw = path.read_bytes()
        except OSError as exc:
            raise ValidationError(f"Cannot read {path}: {exc}") from exc
        snapshots[name] = raw
        sources.append({"file": name, "source_path": str(path), "sha256": digest(raw), "bytes": len(raw)})
    try:
        config = json.loads(snapshots["config.json"])
    except (ValueError, UnicodeError) as exc:
        raise ValidationError(f"Invalid config.json: {exc}") from exc
    if not isinstance(config, dict):
        raise ValidationError("config.json must be an object")
    allowed = {"dataset_id", "analyte", "area_unit", "concentration_unit", "criteria", "criteria_note", "integration_note"}
    if set(config) - allowed:
        raise ValidationError(f"Unknown config fields: {sorted(set(config) - allowed)}")
    for field in ("dataset_id", "analyte", "area_unit", "criteria_note", "integration_note"):
        if not isinstance(config.get(field), str) or not config[field].strip():
            raise ValidationError(f"config.json: {field} is required")
    if config.get("concentration_unit") != "mg/mL":
        raise ValidationError("Only mg/mL concentration input is supported")
    if not isinstance(config.get("criteria", {}), dict):
        raise ValidationError("criteria must be an object")

    def parse(name, required, optional=()):
        try:
            reader = csv.DictReader(io.StringIO(snapshots[name].decode("utf-8-sig")))
            fields = reader.fieldnames or []
            if len(set(fields)) != len(fields) or not set(required) <= set(fields):
                raise ValidationError(f"{name}: required unique columns: {', '.join(required)}")
            unknown = set(fields) - set(required) - set(optional)
            if unknown:
                raise ValidationError(f"{name}: unrecognized columns {sorted(unknown)}")
            rows = []
            for row in reader:
                if None in row or any(v is None for v in row.values()):
                    raise ValidationError(f"{name}:{reader.line_num}: malformed CSV row")
                row = {k: v.strip() for k, v in row.items()}
                row.update(source_csv=name, source_line=reader.line_num)
                rows.append(row)
            return rows
        except UnicodeError as exc:
            raise ValidationError(f"{name}: expected UTF-8 CSV") from exc

    injections = parse("injections.csv", ["injection_id", "kind", "preparation_id", "area",
                       "nominal_concentration_mg_ml", "purity_fraction"],
                       ["source_raw_file", "acquisition_datetime", "integration_reference"])
    preparations = parse("preparations.csv", ["preparation_id", "sample_id", "formulation_mass_mg",
                           "final_volume_ml", "additional_dilution_factor"])
    result = analyze(injections, preparations, config.get("criteria"))
    # Hash the calculations, full input content, settings, and application version.
    runtime = {"python": platform.python_version(),
               **{name: version(name) for name in ("numpy", "openpyxl", "reportlab")}}
    identity = {"app_version": __version__, "runtime": runtime, "config": config,
                "sources": [{k: s[k] for k in ("file", "sha256")} for s in sources], "result": result}
    return {"analysis_id": digest(canonical(identity)), "app_version": __version__, "runtime": runtime,
            "config": config, "sources": sources, "result": result}, snapshots


def validate_review(analysis, review):
    if not isinstance(review, dict):
        raise ValidationError("Review must be a JSON object")
    if review.get("analysis_id") != analysis["analysis_id"]:
        raise ValidationError("Review does not match these inputs, criteria, calculations, and application version")
    for field in ("analyst", "approved_at_utc", "review_note"):
        if not isinstance(review.get(field), str) or not review[field].strip():
            raise ValidationError(f"Review requires {field}")
    for field in ("records_verified", "assignments_confirmed", "approved"):
        if review.get(field) is not True:
            raise ValidationError(f"Review requires {field}=true")
    expected = {f["id"] for f in analysis["result"]["flags"]}
    reasons = review.get("flag_acknowledgments", {})
    if not isinstance(reasons, dict) or set(reasons) != expected:
        raise ValidationError("Review must acknowledge exactly the current flag IDs")
    if any(not isinstance(v, str) or not v.strip() for v in reasons.values()):
        raise ValidationError("Every flag acknowledgment requires a reason")
    try:
        when = datetime.fromisoformat(review["approved_at_utc"])
        if when.tzinfo is None or when.utcoffset().total_seconds() != 0:
            raise ValueError()
        if when > datetime.now(timezone.utc):
            raise ValueError()
    except (ValueError, TypeError):
        raise ValidationError("Approval time must be a valid non-future UTC ISO timestamp") from None


def review_template(analysis):
    return {"analysis_id": analysis["analysis_id"], "analyst": "", "approved_at_utc": "",
            "records_verified": False, "assignments_confirmed": False, "approved": False,
            "review_note": "", "flag_acknowledgments": {
                f["id"]: "" for f in analysis["result"]["flags"]}}
