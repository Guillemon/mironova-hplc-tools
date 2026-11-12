# Contributing

Thank you for helping improve Mironova HPLC Tools.

## Before opening an issue

Search existing issues, include the app and version, operating system, Python version, exact steps, and the smallest useful error message or screenshot. Never attach confidential raw instrument data to a public issue.

## Sample-data policy

Do not upload production `.D` folders, chromatograms, mass spectra, customer information, batch identifiers, or laboratory records. Prefer synthetic data, tiny anonymized traces, metadata-only examples, or a reproducible test. If a private review is necessary, describe the issue first and use your organization's approved secure channel.

## Pull requests

1. Explain the user-visible change and scientific rationale.
2. Add or update tests for parsing and calculations.
3. Update the relevant guide or changelog.
4. Run the app's test command and include the result in the PR description.
5. Keep changes focused; do not commit generated caches, reports, or raw data.

Use clear commit messages such as `fix: handle missing DAD channel` or `test: add single-standard assay case`.

## Review expectations

Every calculation change needs an example with expected values and a note about rounding, units, integration window, and blank handling. Scientific results remain subject to analyst verification and the laboratory quality system.

