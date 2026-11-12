# Mironova HPLC Tools

Local-first tools for analyst review of Agilent HPLC data. This repository contains two companion applications:

- **Single Run UV/MS** (`apps/single-run`) — browse a parent folder of Agilent `.D` injections, open one run, inspect DAD/UV traces and mass spectra, integrate UV peaks, apply analyst-selected blank correction, compare runs, and export a report.
- **CEA Release Analysis** (`apps/cea-release`) — import a CEA sequence, identify standards and samples, apply preparation records/calibration or single-standard calculations, and create release-review reports.

The applications run locally on Windows and macOS. They do not upload instrument files. Use the included Windows `.bat` launchers, or follow the `START_HERE.txt` files in each app folder.

## Important scientific scope

This is analyst-review software, not a validated GMP/GLP release system. It does not replace Agilent OpenLab/MassHunter, an approved laboratory method, second-person review, or your quality system. Verify calculations, integration boundaries, blank choices, sample identity, and report output before using results in a controlled record.

## Data and privacy

Raw Agilent data, cached results, generated reports, and personal computer paths are intentionally excluded from this repository. Do not commit `.D`, `.ch`, `.ms`, `.uv`, `Sample.xml`, confidential batch data, credentials, or client information. For bug reports, use a minimal anonymized fixture or describe the file structure without sharing the data.

## Installation

Each app has its own `pyproject.toml`, Windows requirements file, setup script, and launcher. Python 3.12+ is recommended (Python 3.15 may work, but third-party parser compatibility should be checked).

```text
apps/single-run/Setup_Single_Run_Windows.bat
apps/single-run/Start_Single_Run.bat
apps/cea-release/Setup_Windows.bat
apps/cea-release/Start_CEA_Release.bat
```

The apps use the open-source `rainbow-api` parser for Agilent formats where available. See each app's guide for supported file layouts and troubleshooting.

## Contributing

Issues and pull requests are welcome. Start with [CONTRIBUTING.md](CONTRIBUTING.md), especially the sample-data policy. Suggested contributions include parser compatibility cases, anonymized fixtures, calculation tests, documentation, accessibility improvements, and validation notes.

## License

Original code in this repository is released under the MIT License. Dependencies and bundled third-party components remain under their own licenses. See [LICENSE](LICENSE).

