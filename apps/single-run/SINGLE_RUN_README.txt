MIRONOVA RUN LIBRARY v0.6.3 - WINDOWS

Browse a parent Agilent data folder, open .D injections inside the app, review
integrated chromatograms and tables, and create printable PDFs when needed.

INSTALL OR UPDATE
1. Close the old app's launcher window before starting this version.
2. Extract this ZIP into your existing application folder and replace app files.
   Keep its single_run_results folder: it contains saved data, reports and state.
   Alternatively, extract to a new folder and COPY single_run_results into it.
3. Use 64-bit Python 3.12 with the Python launcher enabled.
   https://www.python.org/ftp/python/3.12.10/python-3.12.10-amd64.exe
   Python 3.15 can remain installed; setup explicitly selects Python 3.12.
4. Run Setup_Windows.bat once. Internet is needed for prebuilt dependencies.
   Microsoft C++ Build Tools are not required. Setup runs the included tests.
5. Run Start_Single_Run.bat. Keep its window open while working.
   If the browser does not open, use http://127.0.0.1:8502.
If setup fails, share setup-log.txt from the application folder.

OPEN THE DATA LIBRARY
- In the sidebar choose Browse parent folder.
- Enter the local path to the parent folder, for example ChemStation\1\Data.
- Click Open / refresh library to list the .D injections recursively.
  Listing the parent folder does NOT decode UV/MS data or automatically open a run.
  Select an injection or click its table row to load its chromatograms and tables.
- Read sample names and comments reads only the small SAMPLE.XML file in each
  injection. No CH/MS/UV headers, signals, acquisition text or integration are
  read or computed during listing. Uncheck it for filenames-only listing, which
  opens no instrument files or sample-metadata files.
- Filter by sequence/subfolder or day, or search filenames, sample names, lots,
  descriptions/comments and method paths present in SAMPLE.XML.
- Select an injection, use Previous / Next, or expand Browse injections and click
  any cell in the injection table. Chromatograms and tables open automatically;
  no checkbox or export is needed. The injection list is expanded by default.
- The parent path is remembered, but the library is not scanned automatically
  after reopening the browser. Click Open / refresh library when ready.
- Saved-library-index loading and automatic restoration of the selected run,
  filters and settings are disabled in this release. Old index/state
  files are ignored; no deletion is needed.
- Within the current browser session you can navigate and adjust settings normally.
- Library grouping uses a date in the folder name, then folder modified date;
  Date source is shown. Acquisition headers become available after selecting a
  run, rather than being read for every injection during folder listing.
- Library scanning stops at 5,000 runs or 50,000 folders and reports that limit.
  Choose a smaller parent folder if needed. Folder links are not followed.
- You can also open a single .D path or upload one .D folder as before.
Library search covers filenames and optional sample XML, not acquisition text,
binary headers, proprietary signal contents, compound identities or m/z peaks.

VIEW AND ADJUST DATA
- Chromatograms & peaks: UV sections sort by wavelength (e.g. 208, 210, 250 nm),
  followed by both raw, unintegrated MS chromatograms. Retention time and Area %
  are emphasized in UV tables: Peak, RT, Area %, Area, then observed m/z.
  Open Chart controls above a graph for X start/end and Y minimum/maximum controls. For example,
  set Y maximum (mAU) to 5000. Reset axes restores default limits.
  Hover to inspect values, drag to pan, and scroll inside the plot to zoom.
  Zoom only changes the view; it does not change integration.
- Minimum area % shown updates the displayed peaks and tables immediately.
  Set it to zero to restore all detected peaks. Optional absolute UV area
  cutoffs and estimated S/N cutoffs are also available. All cutoffs default to zero.
- Show hidden peaks for review temporarily reveals hidden rows and shading.
  It does not change the saved report cutoffs.
- Open Integration settings: enter start/end times in minutes,
  then click Reintegrate with these settings. There is no full-run checkbox to
  override your entered times. Reintegrate full run resets the entire window.
  The active range is displayed prominently and highlighted on each chart.
  Peaks with their apex inside the selected window are integrated. Integration
  boundaries crossing the window are clipped and flagged in exported data.
  The baseline is fitted to the full trace; full raw traces remain available.
  Percentages are recalculated within that integration window. This differs from
  applying cutoffs, which changes the displayed denominator to retained peaks.
- Each channel uses its acquired overlap with the selected window. A channel with
  no overlap has no integrated peaks; its full decoded trace remains available.
- Mass spectra: choose a UV peak/MS scan; inspect observed m/z and
  the strongest measured bins. Links use nearest scan after the chosen time offset.
- Files & metadata: inspect sample records, acquisition headers, source filenames
  and any reading errors. Unsupported sidecars are listed but not decoded.

SAVED ANALYSIS AND REPORTS
Decoded data and integration results are cached under single_run_results.

SINGLE-POINT ASSAY
- Open the collapsed “Single-point assay” panel after selecting an unknown.
- Choose a standard .D injection, enter its known concentration and units, and
  load it. Select the standard UV signal, its integrated peak, and the matching
  unknown peak in the current injection.
- The result is calculated as unknown area / standard area × standard
  concentration × unknown dilution factor ÷ response factor. The default
  dilution and response factors are 1.0.
- This is a provisional single-point external-standard estimate. It is not a
  calibration curve and does not establish linearity, precision, accuracy,
  purity, identity, or release approval. MS data is not used.
Selecting an unchanged run with the same settings reuses saved analysis results,
including after restarting the app; the selected view itself is not restored. Quick file-stat checks detect changes; changed files
are rehashed and decoded. Reopening unchanged runs does not read raw contents or
parse metadata again. Changing the integration window or peak settings recalculates
from saved decoded arrays, with a separate saved result for each set of settings.
Moving between runs and changing display cutoffs do not generate PDFs.
In PDF report, click Create report to export the current integration and cutoffs.
Download/open the PDF to print, or download its accompanying complete data ZIP.
The website also provides a PDF preview if supported by the browser.
Report tables omit Bounds and Review flags, while exports retain those fields.
Source snapshots, CSV peak/trace tables, spectra, NPZ arrays, processing settings,
source checksums and a SHA-256 output manifest accompany each exported report.

AREA PERCENTAGES AND REVIEW
Original Area % = 100 x positive peak area / sum of all positive detected peak
areas within one UV channel and integration window, before cutoffs.
Displayed Area % = 100 x retained peak area / sum of positive retained peak areas.
If no peaks remain, Area % is not estimated.
Wavelengths and MS channels are never pooled. The denominator excludes features
below detection thresholds and nonpositive integrals. Area % measures detector
response, not mass purity, composition, concentration or the CEA formulation assay.
Small peaks can be real. Compare questionable features with an appropriate blank.
Cutoffs use the original percentages once. Remaining peak percentages are then
recalculated to 100% using retained positive UV areas only. Measured areas stay
unchanged. Excluded area is reported against the original denominator.
Original and retained percentages, excluded peaks and reasons remain in exports.
Known area and S/N cutoffs apply independently. Unavailable S/N does not override
a known area cutoff; nonpositive or unknown area is retained for review.

NOISE AND ANALYTICAL LIMITS
Generic UV prominence default: 5 mAU. MS TICs have no baseline correction, peak
picking, integration or area percentages. UV peaks link to nearby MS scans.
Baseline: asymmetric least squares; detection: prominence; boundaries: prominence
fraction with shared-valley partitioning. These are screening settings, not a
qualified analytical method, and may differ from analyst-reviewed CEA integration.
Estimated S/N = baseline-corrected height / robust noise sigma. Sigma uses MAD of
adjacent peak-free residual differences / 0.67449 / sqrt(2), assuming weakly
correlated noise. It uses local peak-free regions, with a global fallback. A manual
noise region can also be specified. Unknown noise is marked not estimable.
Smoothed/correlated noise can overstate S/N. This is not a validated LOD/LOQ or a
pharmacopoeial S/N calculation. Broad features, shared valleys, run/range edges,
and peak identity require analyst review, even at high estimated S/N.
Observed m/z may be a matrix/solvent ion, adduct or fragment. No background
subtraction, confirmed compound identity, polarity or neutral molecular mass is
assigned. MS bins default to 0.05 Da; this is not a claim of mass accuracy.
No calibration, sample weights, standard concentrations or assay calculation occur.

LOCAL OPERATION
Original instrument files are read only; analysis works offline after setup,
without this chat, an API key or a cloud service. The app listens only on this
computer (127.0.0.1). A writable results folder and completed/stable raw files are
needed. This remains a Python application, not a standalone Windows EXE.
The supplied ChemStation CH, MS and UV formats were tested with actual UV/MS data;
no .dx or HRMS qualification is claimed. Windows setup runs tests on the LC computer.

PEAK NAMES AND TEMPORARY HIGHLIGHTS
Click a UV peak marker/label or a cell in its peak table to highlight it in orange.
Click blank graph space or elsewhere in the page to restore normal colors.
The selected peak opens a Peak name field: edit it and click Save peak name.
Names are saved by injection, signal contents, channel and peak retention time.
Original peak IDs, RTs and measured areas are retained. These are analyst labels,
not automatic compound identifications. New PDFs/exports include saved names.

UV BLANK CORRECTION
Open UV blank correction and explicitly select the actual blank .D injection
from the loaded library, or paste its path. Click Apply UV blank correction.
Only UV/mAU channels with matching wavelength and available bandwidth/reference
settings are subtracted. Unmatched channels stay uncorrected and are labelled.
Correction = sample absorbance minus scale x blank absorbance at aligned times.
Default scale is 1 and time shift is 0. Positive shift moves the blank later.
Choose a suitable blank acquired under the same analytical conditions; automatic
RT alignment or scale fitting is not performed. Corrections are provisional.
The shared time range is used without extrapolation. Negative residuals are kept.
A requested integration window uses the shared sample/blank coverage at each
wavelength. UV baseline, peak picking and areas are recalculated after correction. Minimum Area % then applies
to these corrected UV peaks. Corrected UV percentages are not original assay.
Use Traces to show at each wavelength to select Original sample, Blank, Baseline,
and/or Blank-corrected sample. The displayed blank uses recorded absorbance.
Show uncorrected UV restores the original view. Original sample traces, scaled
blank curves, blank source snapshots/checksums and all settings are exported.
Mass chromatograms and m/z spectra are never subtracted, masked or integrated.

For blank-corrected UV, the integration baseline is the difference of independently
estimated sample and aligned/scaled blank AsLS baselines. It is not refitted
through a negative subtraction trough. A substantial negative residual triggers
a screening warning to review blank suitability, scaling and retention alignment;
it does not automatically discard peaks or claim that a blank is suitable.

CHROMATOGRAM ZOOM
Each UV or MS chromatogram has a Mouse mode selector. Inspect peaks supports
peak clicks/names (UV only). Zoom to region lets you drag a rectangle around
a time/response region; it updates the X and Y limits to enlarge that region.
Mouse wheel zoom remains available. Zoom out widens both axes by a factor of
two (time stays within the recorded run). Reset axes restores full-run limits.
Zoom affects the display only; integration ranges, measured areas and peak
percentages are unchanged. Switch back to Inspect peaks to select/rename peaks.

Manual UV integration and removal
- Choose Pan / drag to move around, or Zoom to region to enlarge a region.
- Choose Manual integration and drag horizontally from peak start to end; release to calculate. Areas use the displayed baseline and interpolated endpoints. A new region replaces overlapping integrals to prevent double counting.
- Choose Delete integrals and drag a rectangle around peak markers. Integrals whose apex is inside the box (time and absorbance) are removed together. Area percentages update immediately; one Undo restores the entire group.
- Undo integral edit restores the previous edit for that UV channel. Reset automatic integrals restores that channel to its automatic integration.
- Manual edits are saved for the current injection, integration settings and UV blank correction. Originals and mass spectra are unchanged. Exported reports include updated integrals and an edit history.

UV trace display
Each wavelength has a Traces to show selector. Choose Original sample, Blank, Baseline, and/or Blank-corrected sample. For a simple overlay select Original sample and Blank; for processed data alone select Blank-corrected sample. Blank displays the recorded absorbance without response scaling; any configured time shift is applied and stated. Baseline refers to the active integration signal (blank-corrected when correction is active). Integrals and editing controls appear only when the integration signal is visible. Changing trace visibility does not change integration, tables or reports. Reset axes fits the selected traces.

Blank-corrected integration windows
Integration time changes remain available with a UV blank applied. Each corrected wavelength integrates the intersection of the requested window and its shared sample/blank coverage; the app displays any limited range. A channel with no overlap has no integrated peaks. Blank signals are never extrapolated beyond recorded coverage. A correction error preserves the selected blank and stops results display until corrected or Show uncorrected UV is chosen.

COMPACT CONTROLS AND COMPARISONS (v0.6.2)
Integration settings, UV blank correction and Overlay chromatograms are collapsed
in one row above the results. Area/SN cutoffs are inside Integration settings.
Each chart has a collapsed Chart controls section for trace selection, axes,
mouse modes, Undo and Reset. Click a section to expand or collapse it.

OVERLAY CHROMATOGRAMS
Open Overlay chromatograms. Select injections from the loaded library and/or
paste .D paths from other folders, one per line. Click Apply overlays, then
collapse the section. Matching UV wavelengths and MS channel names are overlaid
using original recorded signals. There is no automatic RT shift or normalization.
UV overlays can be toggled under Chart controls > Traces to show. Reset axes fits
the selected traces. Clear overlays removes all comparisons. Missing or ambiguous
matches are stated per channel. Only explicitly selected comparison injections
are opened; library browsing remains metadata-only. Integration, manual edits,
blank correction and reports still apply to the current injection, not overlays.

PEAK SPECTRA WINDOW
In Inspect peaks mode, click a shaded integral, peak marker/label or table cell.
A Peak spectra popup shows tabs for the raw UV spectrum and linked mass spectra
at the selected peak apex. Scan times are shown. No UV or MS spectral background
subtraction is performed. UV spectral cubes must actually have been recorded;
fixed-wavelength chromatograms cannot supply a spectrum by themselves. If data
or a nearby scan is unavailable, the popup says so. Close the popup to return to
the chromatogram. Open peak spectra below a selected peak reopens it. Peak names
and manual integration remain available in the main page.

FIND BLANKS AND OVERLAYS (v0.6.2)
Both UV blank correction and Overlay chromatograms include independent search,
Folder / sequence and Day filters above the injection selector. Load the parent
folder in Data library first. Text search covers the catalog's filenames, sample
names and comments. Filters use the already loaded metadata and update immediately;
they do not rescan the library or open every injection. Day uses the same date
source as the main library. Choose the injection, then Apply UV blank correction
or Apply overlays. Applying overlays adds to the existing comparisons, allowing
selections across folders and days; changing filters does not remove them.
Clear overlays removes all comparisons. Pasting external .D paths remains available.

COLLAPSIBLE CHANNELS (v0.6.2)
Each UV wavelength and raw MS chromatogram starts collapsed in Chromatograms &
peaks. Click the wavelength/channel heading to open its plot, chart controls and
peak table; click again to hide it. You can keep every channel hidden, review
one at a time, or open several for comparison. Hiding sections changes only
visibility: integration results, peak edits and report contents are retained.
