# Tomato Pollen Phenotyping

R scripts for the manuscript describing an automated tomato pollen phenotyping method. Scripts are run from the repository root and read from `data/`; figures are written to `plots/` (both directories are gitignored due to file size).

## R Scripts

### `R/process_cv_data.R`

Prepares raw computer vision output for downstream analysis. Reads per-frame bounding box inference from the CV pipeline, links wells to accessions via Google Sheets, classifies tracks into biological phenotype classes (ungerminated, germinated, burst, aborted, tube tip), and writes the processed and interpolated tables to `data/processed_data/`.

### `R/plot_classes_over_time.R`

Generates manuscript figures showing pollen class percentages over time. Reads processed and interpolated inference tables from `data/processed_data/` and accession metadata from Google Sheets. Produces plots comparing germination dynamics across temperatures and accessions, including time-to-50%-germination and burst percentage comparisons. Writes figures to `plots/classes_over_time/`.

### `R/run_statistics_on_cv_data.R`

Runs ANCOVA on the processed pollen phenotype data to adjust for time-of-year effects and compare accessions. Reads processed inference and tube length data from `data/processed_data/`, links accessions via Google Sheets, and produces ANCOVA diagnostic and results plots. Writes figures to `plots/ANCOVA/`.

### `R/tube_lengths.R`

Generates manuscript figures for pollen tube length and tube growth speed phenotypes. Reads processed tube length data from `data/processed_data/` and per-frame track data from `data/model_predictions/`. Applies a camera-adjustment factor using the Heinz reference accession and plots tube length distributions, temperature ratios, and per-accession growth speed curves. Writes figures to `plots/tube_lengths/`.

### `R/compare_predictions_with_ground_truth.R`

Validates model accuracy by comparing bounding box predictions against hand-counted ground truth. Reads ground truth counts from `data/model_ground_truth/` and model predictions from `data/model_predictions/`. Computes optimal confidence score thresholds by maximizing r-squared on the training set, then evaluates the selected thresholds on the validation set. Produces confidence threshold plots and ground truth vs. prediction scatter plots. Writes figures to `plots/confidence_threshold_optimization/` and `plots/ground_truth_vs_model_predictions/`.

### `R/model_training_metrics.R`

Generates supplemental figures summarizing model training and selection. Reads TensorBoard training logs from `data/model_training_tensorboard/` and ground truth/prediction files from `data/model_ground_truth/` and `data/model_predictions/`. Produces three figures: (A) a dot plot comparing validation accuracy between models trained on individual cameras vs. both cameras combined, (B) training loss/mAP/recall and combined optimization score over training steps for the combined-camera model, and (C) r-squared vs. confidence threshold for the combined-camera all-classes model on the training set. Writes figures to `plots/model_training_metrics/`.

## Dependencies

All scripts require R with the following packages: `dplyr`, `tidyr`, `ggplot2`, `patchwork`, `stringr`. Additional packages used by specific scripts: `car`, `emmeans` (`run_statistics_on_cv_data.R`), `googlesheets4`, `purrr` (scripts that access accession metadata from Google Sheets).

Scripts that read from Google Sheets (`process_cv_data.R`, `plot_classes_over_time.R`, `run_statistics_on_cv_data.R`, `tube_lengths.R`) require a Google service account key at `~/.credentials/google_sheets_api/service_account.json`.