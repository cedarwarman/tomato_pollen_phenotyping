# Introduction ------------------------------------------------------------
# Manuscript version: tube length and tube growth speed phenotypes.
# Run from repository root. Inputs: data/; figures:
# plots/tube_lengths/.
#
# Tube length plots read from:
#   data/processed_data/processed_tube_lengths.tsv
# Growth speed plots read from:
#   data/model_predictions/2023-06-23_all_tracks_bug_fix.txt
#   + wells_to_accessions from Google Sheets (accession mapping only)

library(dplyr)
library(ggplot2)
library(tidyr)
library(googlesheets4)

pdf(NULL)
set.seed(16)

gs4_auth(path = "~/.credentials/google_sheets_api/service_account.json")

# Temperature colors (consistent with existing tube growth speed plots)
col_26 <- "#00B0F0"
col_34 <- "#C00000"

plot_dir <- file.path(getwd(), "plots", "tube_lengths")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)


# Loading tube length data ------------------------------------------------
tube_lengths <- read.table(
  file = file.path(getwd(), "data", "processed_data", "processed_tube_lengths.tsv"),
  sep = '\t',
  header = TRUE
) %>%
  filter(date > "2021-11-09")


# Per-well mean tube lengths ----------------------------------------------
tube_length_data <- tube_lengths %>%
  mutate(well_name = paste(date, run, well, sep = "_")) %>%
  group_by(camera, temp_target, accession, date, run, well, well_name) %>%
  summarize(
    mean_tube_length = sum(total_length_um) / n(),
    .groups = "drop"
  ) %>%
  mutate(days_since_start = as.numeric(round(difftime(
    as.Date(date), as.Date("2021-11-10"), units = "days"
  ))))

# Camera adjustment using Heinz as reference ------------------------------
# Heinz (CW0000) was run throughout the experiment. The camera-2 mean is
# divided by the camera-1 mean to produce an adjustment factor, then all
# camera-1 lengths are multiplied by it so both cameras are on the same scale.
heinz_mean_camera_1 <- tube_length_data %>%
  filter(accession == "CW0000", camera == 1) %>%
  summarize(mean = mean(mean_tube_length)) %>%
  pull()

heinz_mean_camera_2 <- tube_length_data %>%
  filter(accession == "CW0000", camera == 2) %>%
  summarize(mean = mean(mean_tube_length)) %>%
  pull()

camera_adjustment_factor <- heinz_mean_camera_2 / heinz_mean_camera_1

adjusted_tube_length_data <- tube_length_data %>%
  mutate(
    adjusted_mean_tube_length = ifelse(
      camera == 1,
      mean_tube_length * camera_adjustment_factor,
      mean_tube_length
    )
  )


# Shared theme elements ---------------------------------------------------
tube_length_theme <- function() {
  theme_bw() +
    theme(
      axis.title       = element_text(size = 26, face = "bold"),
      axis.title.y     = element_text(margin = margin(0, 15, 0, 0)),
      axis.text        = element_text(size = 22, face = "bold", color = "black"),
      axis.text.x      = element_blank(),
      plot.title       = element_text(size = 38, face = "bold", margin = margin(0, 0, 10, 0)),
      panel.border     = element_blank(),
      axis.line        = element_line(linewidth = 1, color = "black"),
      axis.ticks       = element_line(linewidth = 1, color = "black"),
      axis.ticks.length = unit(8, "pt"),
      axis.ticks.x     = element_blank(),
      plot.margin      = margin(0.5, 0.5, 0.5, 0.5, "cm"),
      panel.grid       = element_blank(),
      legend.position  = "none"
    )
}


# Plot 1 & 2: Tube length by accession at 26 and 34 ºC -------------------
make_tube_length_plot <- function(temp_string, fill_color, title_string) {
  plot_data <- adjusted_tube_length_data %>%
    filter(temp_target == temp_string) %>%
    group_by(temp_target, accession) %>%
    filter(n() >= 7) %>%
    ungroup()

  ggplot(plot_data,
         aes(x = reorder(accession, adjusted_mean_tube_length, median),
             y = adjusted_mean_tube_length)) +
    geom_boxplot(linewidth = 0.5, color = "black", fill = fill_color,
                 outlier.shape = NA) +
    labs(title = title_string,
         x = "Accession",
         y = "Tube length (µm)") +
    scale_y_continuous(breaks = seq(300, 1500, 200),
                       labels = seq(300, 1500, 200),
                       expand = c(0, 0)) +
    coord_cartesian(ylim = c(300, 1500)) +
    tube_length_theme()

  ggsave(
    filename = file.path(plot_dir, paste0(title_string, ".png")),
    device   = "png",
    width    = 20,
    height   = 10,
    dpi      = 400,
    units    = "in"
  )
}

make_tube_length_plot("26", col_26, "Tube length by accession at 26 \u00baC")
make_tube_length_plot("34", col_34, "Tube length by accession at 34 \u00baC")


# Plot 3: Ratio of mean pollen tube lengths at 26 and 34 ºC --------------
# Day-of-year has no significant effect on tube length ratio (p = 0.805 in
# linear regression; see run_statistics_on_cv_data.R), so unadjusted per-well
# means are used here.
tube_length_ratio_data <- tube_length_data %>%
  group_by(temp_target, accession) %>%
  filter(n() >= 7) %>%
  group_by(accession) %>%
  filter(n_distinct(temp_target) == 2) %>%
  ungroup() %>%
  group_by(temp_target, accession) %>%
  summarize(
    mean_tube_length_val = mean(mean_tube_length),
    .groups = "drop"
  ) %>%
  group_by(accession) %>%
  summarize(
    ratio = mean_tube_length_val[temp_target == 34] /
            mean_tube_length_val[temp_target == 26],
    .groups = "drop"
  )

ggplot(tube_length_ratio_data,
       aes(x = reorder(accession, ratio), y = ratio)) +
  geom_hline(yintercept = 1, linewidth = 1, linetype = 2, color = "gray") +
  geom_point(size = 1.5, color = "black") +
  labs(title = "Ratio of mean pollen tube lengths by accession at 26 and 34 \u00baC",
       x = "Accession",
       y = "Ratio (34/26 \u00baC)") +
  scale_y_continuous(breaks = seq(0.7, 1.2, 0.1),
                     labels = seq(0.7, 1.2, 0.1),
                     limits = c(0.7, 1.2),
                     expand = c(0, 0)) +
  theme_bw() +
  theme(
    axis.title        = element_text(size = 22, face = "bold"),
    axis.title.y      = element_text(margin = margin(0, 11, 0, 0)),
    axis.text         = element_text(size = 18, face = "bold", color = "black"),
    axis.text.x       = element_blank(),
    plot.title        = element_text(size = 26, face = "bold", margin = margin(0, 0, 10, 0)),
    panel.border      = element_blank(),
    axis.line         = element_line(linewidth = 1, color = "black"),
    axis.ticks        = element_line(linewidth = 1, color = "black"),
    axis.ticks.length = unit(8, "pt"),
    axis.ticks.x      = element_blank(),
    plot.margin       = margin(0.5, 0.5, 0.5, 0.5, "cm"),
    panel.grid        = element_blank(),
    legend.position   = "none"
  )

ggsave(
  filename = file.path(plot_dir, "Ratio of mean pollen tube lengths by accession at 26 and 34 \u00baC.png"),
  device   = "png",
  width    = 12,
  height   = 6,
  dpi      = 400,
  units    = "in"
)


# Loading growth speed data -----------------------------------------------
wells_to_accessions <- read_sheet("1yQ5yAKiL6BzwZ-wH-Q44RoUEwMZztTYafzdvVylq6fo")
wells_to_accessions <- wells_to_accessions[, c("date", "run", "well", "accession")]
wells_to_accessions$date <- as.character(wells_to_accessions$date)
wells_to_accessions$run  <- as.character(wells_to_accessions$run)

tracks <- read.table(
  file   = file.path(getwd(), "data", "model_predictions",
                     "2023-06-23_all_tracks_bug_fix.txt"),
  sep    = '\t',
  header = TRUE
)

tracks <- tracks %>%
  filter(object_class == "tube_tip") %>%
  mutate(
    camera = ifelse(as.Date(date) < as.Date("2022-05-27"), 1, 2),
    x      = as.numeric(x),
    y      = as.numeric(y)
  )

tracks <- left_join(tracks, wells_to_accessions, by = c("date", "run", "well"))


# Filtering near-edge tracks ----------------------------------------------
margin_cam_1 <- 35
margin_cam_2 <- 20

tracks <- tracks %>%
  mutate(
    img_x  = ifelse(camera == 1, 2048, 1600),
    img_y  = ifelse(camera == 1, 2048, 1200),
    margin = ifelse(camera == 1, margin_cam_1, margin_cam_2)
  ) %>%
  group_by(date, run, well, track_id) %>%
  filter(all(x > margin & x < img_x - margin)) %>%
  filter(all(y > margin & y < img_y - margin)) %>%
  ungroup()

# Removing tracks that extend beyond the 2-hour endpoint
tracks <- tracks %>%
  filter(ifelse(camera == 1, t <= 82, t <= 65))


# Normalizing track start time --------------------------------------------
tracks <- tracks %>%
  mutate(t = as.numeric(t)) %>%
  group_by(date, run, well, track_id) %>%
  mutate(t = t - min(t)) %>%
  ungroup()


# Calculating speed -------------------------------------------------------
tracks <- tracks %>%
  group_by(date, run, well, track_id) %>%
  mutate(
    dist  = sqrt((x - lag(x))^2 + (y - lag(y))^2),
    dt    = t - lag(t),
    speed = dist / dt
  ) %>%
  ungroup()


# Summarizing to per-well mean speed at each time point -------------------
tracks_summary <- tracks %>%
  group_by(camera, tempc, accession, date, run, well, t) %>%
  summarize(mean_speed = mean(speed, na.rm = TRUE), .groups = "drop") %>%
  filter(ifelse(camera == 1, t <= 82, t <= 65)) %>%
  drop_na()

# Rescale integer time to minutes
tracks_summary <- tracks_summary %>%
  mutate(t_min = ifelse(camera == 1, t * (120 / 82), t * (120 / 65)))

# Convert speed to µm/min
tracks_summary <- tracks_summary %>%
  mutate(
    mean_speed = ifelse(camera == 1,
                        mean_speed * 0.6451613 / 1.463415,
                        mean_speed * 1.190476  / 1.846154)
  )


# Accession list: ≥ 8 reps at both temperatures ---------------------------
accession_list <- tracks_summary %>%
  group_by(tempc, accession, date, run, well) %>%
  summarize(n = n(), .groups = "drop") %>%
  group_by(tempc, accession) %>%
  filter(n() >= 8) %>%
  summarize(n = n(), .groups = "drop") %>%
  group_by(accession) %>%
  filter(n_distinct(tempc) == 2)
accession_list <- unique(accession_list$accession)

tracks_summary_filtered <- tracks_summary %>%
  filter(accession %in% accession_list)


# Shared theme for growth speed plots -------------------------------------
speed_theme <- function() {
  theme_bw() +
    theme(
      axis.title        = element_text(size = 26, face = "bold"),
      axis.title.x      = element_text(margin = margin(15, 0, 0, 0)),
      axis.title.y      = element_text(margin = margin(0, 15, 0, 0)),
      axis.text         = element_text(size = 22, face = "bold", color = "black"),
      plot.title        = element_text(size = 32, face = "bold", margin = margin(0, 0, 10, 0)),
      panel.border      = element_blank(),
      axis.line         = element_line(linewidth = 1, color = "black"),
      axis.ticks        = element_line(linewidth = 1, color = "black"),
      axis.ticks.length = unit(8, "pt"),
      plot.margin       = margin(0.5, 1, 0.5, 0.5, "cm"),
      panel.grid        = element_blank(),
      legend.position   = "none"
    )
}


# Plot 4 & 5: Tube growth speed by accession at 26 and 34 ºC -------------
make_speed_by_accession_plot <- function(temp_string, line_color, title_string) {
  df <- tracks_summary_filtered %>%
    filter(tempc == temp_string)

  ggplot(df, aes(x = t_min, y = mean_speed, group = accession)) +
    geom_line(stat = "smooth", linewidth = 0.5, alpha = 0.35, color = line_color) +
    geom_smooth(linewidth = 2, color = line_color, group = 1) +
    labs(title = title_string,
         x = "Tube growth time (min)",
         y = "Speed (\u00b5m/min)") +
    scale_x_continuous(breaks = seq(0, 120, 20),
                       labels = seq(0, 120, 20),
                       expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, 7)) +
    speed_theme() +
    theme(
      plot.title   = element_text(size = 22, face = "bold", margin = margin(0, 0, 10, 0)),
      axis.title   = element_text(size = 20, face = "bold"),
      axis.title.x = element_text(size = 20, face = "bold", margin = margin(15, 0, 0, 0)),
      axis.title.y = element_text(size = 20, face = "bold", margin = margin(0, 15, 0, 0)),
      axis.text    = element_text(size = 17, face = "bold", color = "black")
    )

  ggsave(
    filename = file.path(plot_dir, paste0(title_string, ".png")),
    device   = "png",
    width    = 7.5,
    height   = 7,
    dpi      = 400,
    units    = "in"
  )
}

make_speed_by_accession_plot("26", col_26, "Tube growth speed by accession at 26 \u00baC")
make_speed_by_accession_plot("34", col_34, "Tube growth speed by accession at 34 \u00baC")


# Plot 6 & 7: Heinz and Tamaulipas growth speed density at 34 ºC ---------
# These use camera 2 data only, which covers the full 2-hour window at
# consistent frame intervals (65 frames over 120 min).
make_single_accession_density_plot <- function(accession_string, temp_string,
                                               line_color, title_string) {
  df <- tracks %>%
    filter(
      accession == accession_string,
      tempc     == temp_string,
      camera    == 2
    ) %>%
    drop_na(speed) %>%
    mutate(
      speed = speed * 1.190476,   # px → µm (camera 2)
      speed = speed / 1.846154,   # integer → min
      t_min = t * (120 / 65)
    )

  ggplot(df, aes(x = t_min, y = speed)) +
    geom_density_2d_filled(alpha = 0.6, bins = 12) +
    scale_fill_grey(start = 0.90, end = 0.00) +
    geom_density_2d(linewidth = 0.25, colour = "black", alpha = 0.2, bins = 12) +
    geom_smooth(linewidth = 2, se = FALSE, color = line_color) +
    labs(title = title_string,
         x = "Tube growth time (min)",
         y = "Speed (\u00b5m/min)") +
    scale_x_continuous(breaks = seq(0, 120, 20),
                       labels = seq(0, 120, 20),
                       expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, 6)) +
    speed_theme()

  ggsave(
    filename = file.path(plot_dir, paste0(title_string, ".png")),
    device   = "png",
    width    = 12,
    height   = 7,
    dpi      = 400,
    units    = "in"
  )
}

make_single_accession_density_plot("CW0000", "34", col_34,
                                   "Heinz tube growth speed at 34 \u00baC")
make_single_accession_density_plot("CW0002", "34", col_34,
                                   "Tamaulipas tube growth speed at 34 \u00baC")
make_single_accession_density_plot("CW0000", "26", col_26,
                                   "Heinz tube growth speed at 26 \u00baC")
make_single_accession_density_plot("CW0002", "26", col_26,
                                   "Tamaulipas tube growth speed at 26 \u00baC")
