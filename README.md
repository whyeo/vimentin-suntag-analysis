# Vimentin SunTag — Tracking pipeline

A compact, documented pipeline to process Vimentin SunTag single-molecule localization data. The package includes an ImageJ/Fiji macro to convert raw microscope files into CSVs and a MATLAB analysis pipeline to track localizations, compute per-track metrics (MSD, displacement, velocity, directional persistence), and generate publication-ready figures.

## Contents

- ImageJ macro: `imagej_script.ijm` — opens raw `.nd2` files, performs bleach correction, background subtraction, creates visualizations, runs ThunderSTORM and exports CSV results.
- MATLAB scripts:
  - `step1_tracking.m` — run per-file tracking (uses `track.m`) and save per-file `.mat` caches.
  - `step2_analyze_per_file.m` — load track caches, compute per-track metrics and per-file figures/tables.
  - `step3_aggregate.m` — combine per-file results and produce condition-level summary figures.
  - `step4_visualize.m` — generate visual summaries (boxplots, bar charts) of condition-level data.
- Helpers / requirements: `track.m` (several implementations may exist in the workspace), plotting helpers such as `formatfig` (optional).

## High-level pipeline

### 1. Image processing & localization (ImageJ/Fiji)

- Run the ImageJ macro (`imagej_script.ijm`) to import microscope stacks, correct background, generate visualizations and run ThunderSTORM to localize single molecules.
- Output: CSV files in `data/<date>/<condition>/output/` (filenames expected to end with `-sub.csv`).

### 2. Tracking (MATLAB) — `step1_tracking.m`

- Reads each CSV, prepares inputs, runs `track.m`, and saves a per-file `*_track.mat` cache.
- Output: `*_track.mat` next to each CSV. The script can also write an optional per-file Excel summary (`*.xlsx`) containing the raw track points and per-track statistics.

### 3. Per-file analysis (MATLAB) — `step2_analyze_per_file.m`

- Loads the track cache, computes per-track statistics (length, displacement, total distance, mean speed), computes MSD and directional persistence, classifies tracks (stationary/short/long), and writes per-file figures and `.mat` analysis files.
- Output: per-file analysis `.mat` (named `*_MSD_DP.mat`) and PNG figure(s). By default these are written to `graphs/<date>/<condition>/`.

### 4. Aggregate / summary (MATLAB) — `step3_aggregate.m`

- Loads per-file analysis `.mat` files, pads and concatenates matrices across files, computes condition-level statistics and creates summary figures ready for publication.
- Output: condition-level figures and a consolidated MAT with assembled results.

### 5. Visualization (MATLAB) — `step4_visualize.m`

- Create condition-level visual summaries from the assembled results produced by previous steps.
- Run `step4_visualize.m` in MATLAB to load the assembled data (expects `graphs/<date>/assemble_track_data.mat`), compute per-condition summaries (fractional distances, speeds, etc.), and generate plots (boxplots and stacked bars) saved under `graphs/<date>/`.
- Before running, edit `experiment_date` and `experiment_conditions` at the top of `step4_visualize.m` to match your dataset and time scaling.

## Required software and plugins

- Fiji / ImageJ
  - Plugins used by the macro: Bio-Formats Importer, Bleach Correction, ThunderSTORM plugin
- MATLAB
  - `track.m` must be available on the MATLAB path (several implementations exist in this workspace). Recommended: place one `track.m` implementation in the repository root (or add `vimentin_suntag/` or `cleaned/` to your MATLAB path). Example candidate files in this workspace include `track.m` in the project root or in subfolders such as `vimentin_suntag/` and `particletracking/` — ensure the correct implementation is on the path before running `step1_tracking.m`.
  - Optional: `formatfig` or other helper functions used to save figures.

## Optional helpers

- `formatfig` (optional): used to export figures with consistent styling and publication settings. If `formatfig` is not available, use MATLAB's `saveas` instead.
- `colors` (optional): returns color palettes used by the plotting code. If missing the scripts use simple fallback color arrays.

## Expected project layout

```
project-root/
  imagej_script.ijm
  data/
    2023-12-17/
      Condition1/
        raw/
          sample1.nd2
          sample2.nd2
        output/
          sample1-sub.csv
          sample1_track.mat
          sample2-sub.csv
          sample2_track.mat
      Condition2/
        ...
  graphs/
    2023-12-17/
      Condition1/
      Condition2/
  step1_tracking.m
  step2_analyze_per_file.m
  step3_aggregate.m
  step4_visualize.m
  track.m (or accessible on MATLAB path)
  README.md
  LICENSE
```

## Quick start

1. Image processing (create CSVs)
   - Open Fiji and run the macro `imagej_script.ijm` (Edit → Macros → Run...).
   - Or run from the command line (replace placeholders):
     ImageJ --ij2 --console --run path/to/imagej_script.ijm "root_path='/path/to/images', input_folder='raw', output_folder='output'"
   - Confirm CSV files appear in `data/<date>/<condition>/output/`.

2. Run MATLAB tracking (`step1_tracking.m`)
   - Edit `step1_tracking.m` to set `experiment_date` and `experiment_conditions`.
   - Run the script in MATLAB. It will create `*_track.mat` files next to each CSV.

3. Run per-file analysis (`step2_analyze_per_file.m`)
   - Run `step2_analyze_per_file.m` in MATLAB. It will read the track caches, compute metrics, and save per-file analysis `.mat` and figures under `graphs/<date>/<condition>/`.

4. Aggregate and plot (`step3_aggregate.m`)
   - Run `step3_aggregate.m` to combine per-file results and generate condition-level figures.

5. Create visual summaries (`step4_visualize.m`)
   - Run `step4_visualize.m` to load assembled data and produce condition-level boxplots and stacked bar charts. Adjust `experiment_date` and `experiment_conditions` at the top of the script before running.

## Notes and troubleshooting

- If you change tracking parameters or the tracker implementation, re-run the ImageJ macro only if localization parameters changed. Otherwise re-run only `step1_tracking.m` or downstream steps as needed.
- Camera parameters (pixel size, frame rate) are used to convert pixels → micrometers and frames → seconds. Verify `pixel_size_nm`, `photon_normalization`, and `frame_rate_hz` in the MATLAB scripts match your setup.
- If a MATLAB script errors on one file, run it interactively for that file to inspect and fix issues; the pipeline is split so you can re-run only the affected stage.

## Packaging for GitHub

- Include the `.ijm` macro and MATLAB scripts in the repo root or a `scripts/` folder.
- Add a clear `LICENSE` file (e.g., MIT) and a short CITATION or USAGE note if you want credit attribution.
- Provide a short example dataset (one small CSV and corresponding track cache) or a link to sample data so users can run a smoke test.

## Related publication

Code developed for the study 'Transport and Organization of Individual Vimentin Filaments Within Dense Networks Revealed by Single Particle Tracking and 3D FIB-SEM' (Renganathan et al., Journal of Cell Biology, Jan 2025). If you use these scripts, please cite: https://doi.org/10.1083/jcb.202406054.
