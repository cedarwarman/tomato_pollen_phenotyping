# Introduction ------------------------------------------------------------
# Manuscript version: class percentages over time from processed inference.
# Run from repository root. Inputs: data/; figures:
# plots/classes_over_time/.
#
# Processed tables (processed_inference.tsv, interpolated_inference.tsv,
# simplified_inference.tsv) must live under data/processed_data/.
# Either unpack supplemental data there, or set REGENERATE_PROCESSED_DATA <- TRUE
# below to rebuild from raw tracks (slow). After a successful regeneration, set
# REGENERATE_PROCESSED_DATA <- FALSE for normal figure-only runs.

library(car)
library(dplyr)
library(ggplot2)
library(googlesheets4)
library(purrr)
library(stringr)
library(tidyr)

data_root <- file.path(getwd(), "data")
processed_dir <- file.path(data_root, "processed_data")
plot_dir_classes <- file.path(getwd(), "plots", "classes_over_time")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir_classes, recursive = TRUE, showWarnings = FALSE)

TRACKS_FILE <- "2023-06-23_all_tracks_bug_fix.txt"
REGENERATE_PROCESSED_DATA <- TRUE

temp_target_fill <- c("26" = "#00B0F0", "34" = "#C00000")

gs4_auth(path = "~/.credentials/google_sheets_api/service_account.json")

wells_to_accessions <- read_sheet("1yQ5yAKiL6BzwZ-wH-Q44RoUEwMZztTYafzdvVylq6fo")
wells_to_accessions <- wells_to_accessions[, c("date", "run", "well", "temp_target", "accession")]
wells_to_accessions$date <- as.character(wells_to_accessions$date)

if (REGENERATE_PROCESSED_DATA) {
  tracks_path <- file.path(data_root, "model_predictions", TRACKS_FILE)
  if (!file.exists(tracks_path)) {
    stop("Regeneration requested but missing: ", tracks_path)
  }
  inference <- read.table(tracks_path, sep = "\t", header = TRUE)

  process_data <- function(df) {
    df <- df %>%
      mutate(name = paste0(
        date, "_run", run, "_", tempc, "C_", well, "_t", str_pad(t, 3, pad = "0")
      )) %>%
      filter(object_class != "tube_tip", object_class != "object_class") %>%
      group_by(name, object_class) %>%
      summarize(count = n(), .groups = "drop") %>%
      group_by(name) %>%
      mutate(percentage = count / sum(count)) %>%
      ungroup() %>%
      complete(name, object_class) %>%
      mutate(count = replace_na(count, 0)) %>%
      mutate(percentage = replace_na(percentage, 0)) %>%
      mutate(
        time = as.integer(str_sub(name, -3, -1)),
        date = str_sub(name, 1, 10),
        run = as.double(str_sub(name, 15, 15)),
        well = str_sub(name, 21, 22)
      ) %>%
      filter(time <= 82)
    df$camera <- NA
    df$camera[as.Date(df$date) < as.Date("2022-05-27")] <- 1
    df$camera[as.Date(df$date) >= as.Date("2022-05-27")] <- 2
    df
  }

  rescale_by_camera <- function(df) {
    df <- df %>% mutate(name = str_sub(name, 1, -6))
    df1 <- df %>% filter(camera == 1)
    df2 <- df %>% filter(camera == 2)
    rescale_time <- function(time, old_range, new_range) {
      old_min <- old_range[1]
      old_max <- old_range[2]
      new_min <- new_range[1]
      new_max <- new_range[2]
      ((time - old_min) / (old_max - old_min)) * (new_max - new_min) + new_min
    }
    df2 <- df2 %>%
      mutate(time = rescale_time(time, old_range = c(0, 65), new_range = c(0, 82))) %>%
      mutate(time = round(time, 2))
    unique_names <- unique(df2$name)
    process_name <- function(name) {
      df_name <- df2[df2$name == name, ]
      time_df <- expand_grid(
        time = seq(0, max(df_name$time), by = 0.01),
        name = unique(df_name$name),
        object_class = unique(df_name$object_class),
        date = unique(df_name$date),
        run = unique(df_name$run),
        well = unique(df_name$well),
        camera = unique(df_name$camera)
      )
      df_name_processed <- time_df %>%
        left_join(df_name, by = c("time", "name", "object_class", "date", "run", "well", "camera")) %>%
        group_by(name, object_class, date, run, well, camera) %>%
        arrange(time) %>%
        mutate(
          count = approx(time, count, time, rule = 2)$y,
          percentage = approx(time, percentage, time, rule = 2)$y
        ) %>%
        filter(floor(time) == time) %>%
        ungroup()
      df_name_processed
    }
    df_interpolated <- unique_names %>% map(process_name) %>% list_rbind()
    rbind(df1, df_interpolated)
  }

  processed_inference <- process_data(inference)
  interpolated <- rescale_by_camera(processed_inference)
  interpolated <- left_join(interpolated, wells_to_accessions, by = c("date", "run", "well"))
  simplified_df <- interpolated %>%
    group_by(accession, temp_target, time, object_class) %>%
    summarize(mean_percentage = mean(percentage), .groups = "drop") %>%
    mutate(time = as.integer(time)) %>%
    filter(object_class != "aborted")

  write.table(processed_inference,
    file = file.path(processed_dir, "processed_inference.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  write.table(interpolated,
    file = file.path(processed_dir, "interpolated_inference.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  write.table(simplified_df,
    file = file.path(processed_dir, "simplified_inference.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
} else {
  processed_inference <- read.table(
    file = file.path(processed_dir, "processed_inference.tsv"),
    sep = "\t", header = TRUE
  )
  interpolated <- read.table(
    file = file.path(processed_dir, "interpolated_inference.tsv"),
    sep = "\t", header = TRUE
  )
  simplified_df <- read.table(
    file = file.path(processed_dir, "simplified_inference.tsv"),
    sep = "\t", header = TRUE
  )
}


# Making plots with the simplified data -----------------------------------
make_plot <- function(input_df, image_name) {
  color_vec <- c("#FF00FF", # burst
                 "#11e00d", # germinated
                 "#1b74fa", # ungerminated
                 "#FFB000", # unknown_germinated
                 "#787878", # aborted
                 "#ffa6db", # tube_tip_burst
                 "#fffa70", # tube_tip_bulging
                 "#a8ffe1") # tube_tip
  names(color_vec) <- c("burst", 
                        "germinated", 
                        "ungerminated", 
                        "unknown_germinated", 
                        "aborted", 
                        "tube_tip_burst",
                        "tube_tip_bulging",
                        "tube_tip")
  
  ggplot(input_df, aes(x = time, y = mean_percentage, color = object_class, linetype = accession)) +
    geom_line(linewidth = 0.5, alpha = 0.3) +
    scale_color_manual(values = color_vec,
                       name = "Class",
                       breaks = c("ungerminated", 
                                  "germinated", 
                                  "burst", 
                                  "aborted", 
                                  "unknown_germinated", 
                                  "tube_tip", 
                                  "tube_tip_burst", 
                                  "tube_tip_bulging"),
                       labels = c("Ungerminated", 
                                  "Germinated", 
                                  "Burst", 
                                  "Aborted", 
                                  "Unknown germinated", 
                                  "Tube tip", 
                                  "Tube tip burst", 
                                  "Tube tip bulging"),
                       limits = force,
                       guide = guide_legend(override.aes = list(alpha = 1, linewidth = 2))) +
    scale_linetype_manual(values = rep.int(1, 228), guide = "none") +
    scale_x_continuous(breaks = c(0, 20, 40, 60, 80),
                       labels = c(15, 45, 75, 105, 135),
                       # limits = c(0, 82),
                       limits = c(0, 80),
                       expand = c(0, 0)) +
    scale_y_continuous(breaks = c(0, 0.25, .5, .75, 1),
                       labels = c("0%", "25%", "50%", "75%", "100%"),
                       limits = c(0, 1),
                       expand = c(0, 0)) +
    # coord_cartesian(xlim = c(0, 135)) +
    labs(title = image_name,
         x = "Time (minutes)",
         y = "Class percentage") +
    theme_bw() +
    theme(axis.title = element_text(size = 26, face = 'bold'),
          axis.text = element_text(size = 22, face = 'bold', color = 'black'),
          axis.text.x = element_text(size = 26, face = 'bold', color = 'black'),
          plot.title = element_text(size = 28, face = 'bold', margin = margin(0, 0, 10, 0)),
          panel.border = element_blank(),
          axis.line = element_line(linewidth = 1, color = 'black'),
          axis.ticks = element_line(linewidth = 1, color = 'black'),
          axis.ticks.length = unit(8, 'pt'),
          plot.margin = margin(0.5, 0.5, 0.5, 0.5, 'cm'),
          panel.grid = element_blank(),
          legend.position = 'right',
          legend.title = element_text(size = 18, face = 'bold', color = 'black'),
          legend.text = element_text(size = 14, face = 'bold', color = 'black'),
          strip.background = element_blank(),
          strip.placement = "outside")
  
  ggsave(filename = file.path(plot_dir_classes, paste0(image_name, ".png")),
         device = 'png',
         width = 14,
         height = 8,
         dpi = 400,
         units = 'in')
}

make_plot_with_lines <- function(input_df, image_name) {
  color_vec <- c("#FF00FF", # burst
                 "#11e00d", # germinated
                 "#1b74fa", # ungerminated
                 "#FFB000", # unknown_germinated
                 "#787878", # aborted
                 "#ffa6db", # tube_tip_burst
                 "#fffa70", # tube_tip_bulging
                 "#a8ffe1") # tube_tip
  names(color_vec) <- c("burst", 
                        "germinated", 
                        "ungerminated", 
                        "unknown_germinated", 
                        "aborted", 
                        "tube_tip_burst",
                        "tube_tip_bulging",
                        "tube_tip")
  
  ggplot(input_df, aes(x = time, y = mean_percentage, color = object_class)) +
    geom_line(aes(linetype = accession), linewidth = 0.5, alpha = 0.1) +
    geom_smooth(span = 0.4, se = FALSE, linewidth = 2, fullrange = TRUE) +
    scale_color_manual(values = color_vec,
                       name = "Class",
                       breaks = c("ungerminated", 
                                  "germinated", 
                                  "burst", 
                                  "aborted", 
                                  "unknown_germinated", 
                                  "tube_tip", 
                                  "tube_tip_burst", 
                                  "tube_tip_bulging"),
                       labels = c("Ungerminated", 
                                  "Germinated", 
                                  "Burst", 
                                  "Aborted", 
                                  "Unknown germinated", 
                                  "Tube tip", 
                                  "Tube tip burst", 
                                  "Tube tip bulging"),
                       limits = force) +
    scale_linetype_manual(values = rep.int(1, 230), guide = "none") +
    scale_x_continuous(breaks = c(0, 20, 40, 60, 80),
                       labels = c(15, 45, 75, 105, 135),
                       limits = c(0, 80),
                       expand = c(0, 0)) +
    scale_y_continuous(breaks = c(0, 0.25, .5, .75, 1),
                       labels = c("0%", "25%", "50%", "75%", "100%"),
                       limits = c(0, 1),
                       expand = c(0, 0)) +
    labs(title = image_name,
         x = "Time (minutes)",
         y = "Class percentage") +
    theme_bw() +
    theme(axis.title = element_text(size = 26, face = 'bold'),
          axis.text = element_text(size = 22, face = 'bold', color = 'black'),
          axis.text.x = element_text(size = 26, face = 'bold', color = 'black'),
          plot.title = element_text(size = 28, face = 'bold', margin = margin(0, 0, 10, 0)),
          panel.border = element_blank(),
          axis.line = element_line(linewidth = 1, color = 'black'),
          axis.ticks = element_line(linewidth = 1, color = 'black'),
          axis.ticks.length = unit(8, 'pt'),
          plot.margin = margin(0.5, 0.5, 0.5, 0.5, 'cm'),
          panel.grid = element_blank(),
          legend.position = 'right',
          legend.title = element_text(size = 18, face = 'bold', color = 'black'),
          legend.text = element_text(size = 14, face = 'bold', color = 'black'),
          strip.background = element_blank(),
          strip.placement = "outside")
  
  ggsave(filename = file.path(plot_dir_classes, paste0(image_name, "_with_lines.png")),
         device = 'png',
         width = 14,
         height = 8,
         dpi = 400,
         units = 'in')
}

make_plot_nappn <- function(input_df, image_name) {
  color_vec <- c("#FF00FF", # burst
                 "#11e00d", # germinated
                 "#1b74fa", # ungerminated
                 "#FFB000", # unknown_germinated
                 "#787878", # aborted
                 "#ffa6db", # tube_tip_burst
                 "#fffa70", # tube_tip_bulging
                 "#a8ffe1") # tube_tip
  names(color_vec) <- c("burst", 
                        "germinated", 
                        "ungerminated", 
                        "unknown_germinated", 
                        "aborted", 
                        "tube_tip_burst",
                        "tube_tip_bulging",
                        "tube_tip")
  
  ggplot(input_df, aes(x = time, y = mean_percentage, color = object_class)) +
    geom_line(aes(linetype = accession), linewidth = 0.5, alpha = 0.15) +
    geom_smooth(span = 0.4, se = FALSE, linewidth = 2, fullrange = TRUE) +
    scale_color_manual(values = color_vec,
                       name = "Class",
                       breaks = c("ungerminated", 
                                  "germinated", 
                                  "burst", 
                                  "aborted", 
                                  "unknown_germinated", 
                                  "tube_tip", 
                                  "tube_tip_burst", 
                                  "tube_tip_bulging"),
                       labels = c("Ungerminated", 
                                  "Germinated", 
                                  "Burst", 
                                  "Aborted", 
                                  "Unknown germinated", 
                                  "Tube tip", 
                                  "Tube tip burst", 
                                  "Tube tip bulging"),
                       limits = force) +
    scale_linetype_manual(values = rep.int(1, 230), guide = "none") +
    scale_x_continuous(breaks = c(0, 20, 40, 60, 80),
                       labels = c(15, 45, 75, 105, 135),
                       limits = c(0, 80),
                       expand = c(0, 0)) +
    scale_y_continuous(breaks = c(0, 0.25, .5, .75, 1),
                       labels = c("0%", "25%", "50%", "75%", "100%"),
                       limits = c(0, 1),
                       expand = c(0, 0)) +
    labs(title = image_name,
         x = "Time (minutes)",
         y = "Class percentage") +
    theme_bw() +
    theme(axis.title = element_text(size = 26, face = 'bold'),
          axis.text = element_text(size = 22, face = 'bold', color = 'black'),
          axis.text.x = element_text(size = 26, face = 'bold', color = 'black'),
          plot.title = element_blank(),
          # plot.title = element_text(size = 28, face = 'bold', margin = margin(0, 0, 10, 0)),
          panel.border = element_blank(),
          axis.line = element_line(linewidth = 1, color = 'black'),
          axis.ticks = element_line(linewidth = 1, color = 'black'),
          axis.ticks.length = unit(8, 'pt'),
          plot.margin = margin(0.5, 1, 0.5, 0.5, 'cm'),
          panel.grid = element_blank(),
          legend.position = 'none',
          # legend.position = 'bottom',
          # legend.title = element_text(size = 18, face = 'bold', color = 'black'),
          # legend.text = element_text(size = 14, face = 'bold', color = 'black'),
          strip.background = element_blank(),
          strip.placement = "outside")
  
  ggsave(filename = file.path(plot_dir_classes, paste0(image_name, "_nappn.png")),
         device = 'png',
         width = 14,
         height = 7,
         dpi = 400,
         units = 'in')
}

make_plot(simplified_df[simplified_df$temp_target == 26, ], "All accessions at 26 ºC")
make_plot(simplified_df[simplified_df$temp_target == 34, ], "All accessions at 34 ºC")


make_plot_with_lines(simplified_df[simplified_df$temp_target == 26, ], "All accessions at 26 ºC")
make_plot_with_lines(simplified_df[simplified_df$temp_target == 34, ], "All accessions at 34 ºC")

make_plot_nappn(simplified_df[simplified_df$temp_target == 26, ], "All accessions at 26 ºC")
make_plot_nappn(simplified_df[simplified_df$temp_target == 34, ], "All accessions at 34 ºC")



# Making plots with accessions highlighted --------------------------------
make_plot_with_highlight <- function(input_df, image_name, accession_id) {
  color_vec <- c("#FF00FF", # burst
                 "#11e00d", # germinated
                 "#1b74fa", # ungerminated
                 "#FFB000", # unknown_germinated
                 "#787878", # aborted
                 "#ffa6db", # tube_tip_burst
                 "#fffa70", # tube_tip_bulging
                 "#a8ffe1") # tube_tip
  names(color_vec) <- c("burst", 
                        "germinated", 
                        "ungerminated", 
                        "unknown_germinated", 
                        "aborted", 
                        "tube_tip_burst",
                        "tube_tip_bulging",
                        "tube_tip")
  
  ggplot(input_df, aes(x = time, y = mean_percentage, color = object_class)) +
    geom_line(aes(linetype = accession), linewidth = 0.5, alpha = 0.08) +
    geom_line(stat = "smooth", span = 0.4, se = FALSE, linewidth = 1.5, fullrange = TRUE, alpha = 0.5) +
    scale_color_manual(values = color_vec,
                       name = "Class",
                       breaks = c("ungerminated", 
                                  "germinated", 
                                  "burst", 
                                  "aborted", 
                                  "unknown_germinated", 
                                  "tube_tip", 
                                  "tube_tip_burst", 
                                  "tube_tip_bulging"),
                       labels = c("Ungerminated", 
                                  "Germinated", 
                                  "Burst", 
                                  "Aborted", 
                                  "Unknown germinated", 
                                  "Tube tip", 
                                  "Tube tip burst", 
                                  "Tube tip bulging"),
                       limits = force) +
    scale_linetype_manual(values = rep.int(1, 228), guide = "none") +
    geom_line(data = input_df[input_df$accession == {{accession_id}}, ], linetype = '41', linewidth = 2) +
    scale_x_continuous(breaks = c(0, 20, 40, 60, 80),
                       labels = c(15, 45, 75, 105, 135),
                       limits = c(0, 80),
                       expand = c(0, 0)) +
    scale_y_continuous(breaks = c(0, 0.25, .5, .75, 1),
                       labels = c("0%", "25%", "50%", "75%", "100%"),
                       limits = c(0, 1),
                       expand = c(0, 0)) +
    labs(title = image_name,
         x = "Time (minutes)",
         y = "Class percentage") +
    theme_bw() +
    theme(axis.title = element_text(size = 26, face = 'bold'),
          axis.text = element_text(size = 22, face = 'bold', color = 'black'),
          axis.text.x = element_text(size = 26, face = 'bold', color = 'black'),
          plot.title = element_text(size = 28, face = 'bold', margin = margin(0, 0, 10, 0)),
          panel.border = element_blank(),
          axis.line = element_line(linewidth = 1, color = 'black'),
          axis.ticks = element_line(linewidth = 1, color = 'black'),
          axis.ticks.length = unit(8, 'pt'),
          plot.margin = margin(0.5, 0.5, 0.5, 0.5, 'cm'),
          panel.grid = element_blank(),
          legend.position = 'right',
          legend.title = element_text(size = 18, face = 'bold', color = 'black'),
          legend.text = element_text(size = 14, face = 'bold', color = 'black'),
          strip.background = element_blank(),
          strip.placement = "outside")
  
  ggsave(filename = file.path(plot_dir_classes, paste0(image_name, "_with_highlight.png")),
         device = 'png',
         width = 14,
         height = 8,
         dpi = 400,
         units = 'in')
}

make_plot_with_highlight(simplified_df[simplified_df$temp_target == 26, ], "Heinz at 26 °C", "CW0000")
make_plot_with_highlight(simplified_df[simplified_df$temp_target == 34, ], "Heinz at 34 °C", "CW0000")

make_plot_with_highlight(simplified_df[simplified_df$temp_target == 26, ], "Tamaulipas at 26 °C", "CW0002")
make_plot_with_highlight(simplified_df[simplified_df$temp_target == 34, ], "Tamaulipas at 34 °C", "CW0002")

make_plot_with_highlight(simplified_df[simplified_df$temp_target == 26, ], "Nagcarlang at 26 °C", "CW0003")
make_plot_with_highlight(simplified_df[simplified_df$temp_target == 34, ], "Nagcarlang at 34 °C", "CW0003")

make_plot_with_highlight(simplified_df[simplified_df$temp_target == 26, ], "Malintka at 26 °C", "CW0004")
make_plot_with_highlight(simplified_df[simplified_df$temp_target == 34, ], "Malintka at 34 °C", "CW0004")

make_plot_with_highlight(simplified_df[simplified_df$temp_target == 26, ], "Gold Nugget at 26 °C", "CW0005")
make_plot_with_highlight(simplified_df[simplified_df$temp_target == 34, ], "Gold Nugget at 34 °C", "CW0005")

make_plot_with_highlight(simplified_df[simplified_df$temp_target == 26, ], "San Marzano at 26 °C", "CW0006")
make_plot_with_highlight(simplified_df[simplified_df$temp_target == 34, ], "San Marzano at 34 °C", "CW0006")

make_plot_with_highlight(simplified_df[simplified_df$temp_target == 26, ], "Hotset at 26 °C", "CW0007")
make_plot_with_highlight(simplified_df[simplified_df$temp_target == 34, ], "Hotset at 34 °C", "CW0007")

# Figuring out which line is super resistant to heat stress
just_germinated_and_sixty <- simplified_df[simplified_df$time == 60 & simplified_df$object_class == "germinated" & simplified_df$temp_target == 34, ]
just_germinated_and_sixty[which.max(just_germinated_and_sixty$mean_percentage), ]

make_plot_with_highlight(simplified_df[simplified_df$temp_target == 26, ], "CW0169 at 26 °C", "CW0169")
make_plot_with_highlight(simplified_df[simplified_df$temp_target == 34, ], "CW0169 at 34 °C", "CW0169")

# # Some for Ravi
# make_plot_with_highlight(simplified_df[simplified_df$temp_target == 26, ], "CW0079 at 26 °C", "CW0079")
# make_plot_with_highlight(simplified_df[simplified_df$temp_target == 34, ], "CW0079 at 34 °C", "CW0079")
# 
# make_plot_with_highlight(simplified_df[simplified_df$temp_target == 26, ], "CW0119 at 26 °C", "CW0119")
# make_plot_with_highlight(simplified_df[simplified_df$temp_target == 34, ], "CW0119 at 34 °C", "CW0119")
# 
# 
# Figuring out which line is super sensitive to heat stress
just_burst_and_forty <- simplified_df[simplified_df$time == 40 & simplified_df$object_class == "burst" & simplified_df$temp_target == 34, ]
just_burst_and_forty[which.max(just_burst_and_forty$mean_percentage), ]
# 
# # Looks like it's CW0045
make_plot_with_highlight(simplified_df[simplified_df$temp_target == 26, ], "CW0045 at 26 °C", "CW0045")
make_plot_with_highlight(simplified_df[simplified_df$temp_target == 34, ], "CW0045 at 34 °C", "CW0045")


# Plot with just Heinz ----------------------------------------------------
# I'm going to make a 1-off plot with just Heinz for the PAG presentation
make_plot_with_one <- function(input_df, image_name, accession_id) {
  color_vec <- c("#FF00FF", # burst
                 "#11e00d", # germinated
                 "#1b74fa", # ungerminated
                 "#FFB000", # unknown_germinated
                 "#787878", # aborted
                 "#ffa6db", # tube_tip_burst
                 "#fffa70", # tube_tip_bulging
                 "#a8ffe1") # tube_tip
  names(color_vec) <- c("burst", 
                        "germinated", 
                        "ungerminated", 
                        "unknown_germinated", 
                        "aborted", 
                        "tube_tip_burst",
                        "tube_tip_bulging",
                        "tube_tip")
  
  input_df <- input_df[input_df$accession == accession_id, ]
  
  ggplot(input_df, aes(x = time, y = mean_percentage, color = object_class)) +
    geom_line(linewidth = 2) +
    scale_color_manual(values = color_vec,
                       name = "Class",
                       breaks = c("ungerminated", 
                                  "germinated", 
                                  "burst", 
                                  "aborted", 
                                  "unknown_germinated", 
                                  "tube_tip", 
                                  "tube_tip_burst", 
                                  "tube_tip_bulging"),
                       labels = c("Ungerminated", 
                                  "Germinated", 
                                  "Burst", 
                                  "Aborted", 
                                  "Unknown germinated", 
                                  "Tube tip", 
                                  "Tube tip burst", 
                                  "Tube tip bulging"),
                       limits = force) +
    scale_x_continuous(breaks = c(0, 20, 40, 60, 80),
                       labels = c(15, 45, 75, 105, 135),
                       limits = c(0, 80),
                       expand = c(0, 0)) +
    scale_y_continuous(breaks = c(0, 0.25, .5, .75, 1),
                       labels = c("0%", "25%", "50%", "75%", "100%"),
                       limits = c(0, 1),
                       expand = c(0, 0)) +
    labs(title = image_name,
         x = "Time (minutes)",
         y = "Class percentage") +
    theme_bw() +
    theme(axis.title = element_text(size = 26, face = 'bold'),
          axis.text = element_text(size = 22, face = 'bold', color = 'black'),
          axis.text.x = element_text(size = 26, face = 'bold', color = 'black'),
          plot.title = element_text(size = 28, face = 'bold', margin = margin(0, 0, 10, 0)),
          panel.border = element_blank(),
          axis.line = element_line(linewidth = 1, color = 'black'),
          axis.ticks = element_line(linewidth = 1, color = 'black'),
          axis.ticks.length = unit(8, 'pt'),
          plot.margin = margin(0.5, 0.5, 0.5, 0.5, 'cm'),
          panel.grid = element_blank(),
          legend.position = 'right',
          legend.title = element_text(size = 18, face = 'bold', color = 'black'),
          legend.text = element_text(size = 14, face = 'bold', color = 'black'),
          strip.background = element_blank(),
          strip.placement = "outside")
  
  ggsave(filename = file.path(plot_dir_classes, paste0(image_name, "_only.png")),
         device = 'png',
         width = 14,
         height = 8,
         dpi = 400,
         units = 'in')
}

make_plot_with_one(simplified_df[simplified_df$temp_target == 26, ], "Heinz at 26 °C", "CW0000")


























#Percent burst at 80 minutes plot ----------------------------------------
# It would be nice to make a plot that's just the percent burst at 80 minutes, 
# So I can see them all with error bars like I did with the flower measurements.
# I'll use that plot as a starting point.

plot_burst_80_34 <- ggplot(data = simplified_df[simplified_df$object_class == "burst" & simplified_df$temp_target == 34, ], aes(x = reorder(accession, mean_percentage, median), 
                                                                                                     y = mean_percentage)) +
  geom_boxplot(linewidth = 0.5, color = "black") +
  labs(title = "Pollen grain burst at 34 ºC, 80 minutes after adding media",
       y = "Percentage") +
  scale_y_continuous(breaks = c(0, 0.25, .5, .75, 1),
                     labels = c("0%", "25%", "50%", "75%", "100%"),
                     limits = c(0, 1),
                     expand = c(0, 0)) +
  theme_bw() +
  theme(axis.title = element_text(size = 26, face = 'bold'),
        axis.text = element_text(size = 22, face = 'bold', color = 'black'),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1, face = 'bold', color = 'black'),
        plot.title = element_text(size = 28, face = 'bold', margin = margin(0, 0, 10, 0)),
        axis.title.x = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(linewidth = 1, color = 'black'),
        axis.ticks = element_line(linewidth = 1, color = 'black'),
        axis.ticks.length = unit(8, 'pt'),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, 'cm'),
        panel.grid = element_blank(),
        legend.position = 'none')

ggsave(filename = file.path(plot_dir_classes, "burst_at_80_mins_and_34C.png"),
       plot = plot_burst_80_34,
       device = 'png',
       width = 25,
       height = 9,
       dpi = 400,
       units = 'in')


# Getting info about specific accessions ----------------------------------
# I want to pull the vids for specific accessions, so here they are:
heinz_info <- wells_to_accessions[wells_to_accessions$accession == "CW0000", ]
tamaulipas_info <- wells_to_accessions[wells_to_accessions$accession == "CW0002", ]
CW0164_info <- wells_to_accessions[wells_to_accessions$accession == "CW0164", ]

CW0079_info <- wells_to_accessions[wells_to_accessions$accession == "CW0079", ]


# Calculating time to 50% germinated --------------------------------------
# I want to know the median time until 50% germinated for 26 and 34 C, and maybe make a scatter plot.
simplified_df_germ_26 <- simplified_df %>%
  filter(temp_target == 26, object_class == "germinated") %>%
  group_by(accession) %>%
  summarize(
    time_to_50 = ifelse(any(mean_percentage >= 0.50), 
                        min(time[mean_percentage >= 0.50], na.rm = TRUE), 
                        NA),  # Assign NA if never reaches 50%
    reaches_50 = any(mean_percentage >= 0.50)  # Flag if it ever reaches 50%
  ) %>%
  mutate(time_to_50 = ifelse(!reaches_50, 82, time_to_50))
median_time_germ_26 <- median(simplified_df_germ_26$time_to_50, na.rm = TRUE) # 7

simplified_df_germ_34 <- simplified_df %>%
  filter(temp_target == 34, object_class == "germinated") %>%
  group_by(accession) %>%
  summarize(
    time_to_50 = ifelse(any(mean_percentage >= 0.50), 
                        min(time[mean_percentage >= 0.50], na.rm = TRUE), 
                        NA),  # Assign NA if never reaches 50%
    reaches_50 = any(mean_percentage >= 0.50)  # Flag if it ever reaches 50%
  ) %>%
  mutate(time_to_50 = ifelse(!reaches_50, 82, time_to_50))
median_time_germ_34 <- median(simplified_df_germ_34$time_to_50, na.rm = TRUE) # 4

# Doing some statistics:
# Compute time to 50% for both temperatures
df_germ_processed <- simplified_df %>%
  filter(object_class == "germinated") %>%
  group_by(accession, temp_target) %>%
  summarize(
    time_to_50 = ifelse(any(mean_percentage >= 0.50), 
                        min(time[mean_percentage >= 0.50], na.rm = TRUE), 
                        NA),  # Assign NA if never reaches 50%
    reaches_50 = any(mean_percentage >= 0.50)  # Flag if it ever reaches 50%
  ) %>%
  mutate(time_to_50 = ifelse(!reaches_50, 82, time_to_50)) %>%  # Assign max_time if never reaches 50%
  ungroup()

# Shapiro-Wilk errors when x has <3 values or no variance (all identical).
safe_shapiro_test <- function(x, label = deparse(substitute(x))) {
  x <- x[is.finite(x)]
  if (length(x) < 3) {
    return(structure(
      list(
        statistic = c(W = NA_real_),
        p.value = NA_real_,
        method = "Shapiro-Wilk normality test",
        data.name = paste0(label, " (n<3 after filtering)"),
        alternative = "not enough data"
      ),
      class = "htest"
    ))
  }
  if (length(unique(x)) < 2) {
    return(structure(
      list(
        statistic = c(W = NA_real_),
        p.value = NA_real_,
        method = "Shapiro-Wilk normality test",
        data.name = paste0(label, " (all values identical)"),
        alternative = "all values identical"
      ),
      class = "htest"
    ))
  }
  shapiro.test(x)
}

# Check normality. Shapiro-Wilk test. Null hypothesis is that the data is normally distributed.
# So small p-value rejects null hypothesis and suggests data is not normally distributed.
shapiro_26 <- safe_shapiro_test(
  df_germ_processed$time_to_50[df_germ_processed$temp_target == 26],
  label = "time_to_50 at 26C"
)
shapiro_34 <- safe_shapiro_test(
  df_germ_processed$time_to_50[df_germ_processed$temp_target == 34],
  label = "time_to_50 at 34C"
)
shapiro_26
shapiro_34

# Data is non-normal so Wilcoxon rank-sum test
test_result_germ <- wilcox.test(time_to_50 ~ temp_target, data = df_germ_processed)
test_result_germ

# Also comparing variability
df_germination_variability <- df_germ_processed %>%
  group_by(temp_target) %>%
  summarize(
    variance = var(time_to_50, na.rm = TRUE),
    sd = sd(time_to_50, na.rm = TRUE),
    mean = mean(time_to_50, na.rm = TRUE),
    cv = sd / mean  # Coefficient of Variation
  )

# Perform Levene’s test to compare variance between temperatures
levene_time_to_50_germ <- leveneTest(time_to_50 ~ as.factor(temp_target), 
                                     data = df_germ_processed)
levene_time_to_50_germ


plot_time_to_50_germination <- ggplot(df_germ_processed, aes(x = as.factor(temp_target), y = time_to_50, fill = as.factor(temp_target))) +
  geom_boxplot(color = "black", linewidth = 1, alpha = 0.7) +
  scale_fill_manual(
    values = temp_target_fill,
    name = "Temperature (°C)",
    breaks = c("26", "34"),
    labels = c("26 °C", "34 °C")
  ) +
  labs(title = "Time to 50% Germination", x = "Temperature (°C)", y = "Minutes") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 26, face = "bold"),
    axis.title.x = element_text(size = 26, face = "bold", margin = margin(t = 12)),
    axis.text = element_text(size = 22, face = "bold", color = "black"),
    axis.text.x = element_text(size = 22, face = "bold", color = "black"),
    plot.title = element_text(size = 28, face = "bold", margin = margin(0, 0, 10, 0)),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 1, color = "black"),
    axis.ticks = element_line(linewidth = 1, color = "black"),
    axis.ticks.length = unit(8, "pt"),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
    panel.grid = element_blank(),
    legend.position = "none"
  )

ggsave(
  filename = file.path(plot_dir_classes, "time_to_50pct_germination_26_vs_34.png"),
  plot = plot_time_to_50_germination,
  device = "png",
  width = 7,
  height = 8,
  dpi = 400,
  units = "in"
)


# Percent burst at end of experiment --------------------------------------
df_burst_processed <- simplified_df %>%
  filter(object_class == "burst", time == 82) %>%
  group_by(temp_target) %>%
  summarize(median_burst = median(mean_percentage, na.rm = TRUE))
  
# Check normality of burst percentage at time 82 for each temperature
shapiro_burst_26 <- safe_shapiro_test(
  simplified_df$mean_percentage[simplified_df$temp_target == 26 &
    simplified_df$object_class == "burst" &
    simplified_df$time == 82],
  label = "burst mean_percentage at 26C, time 82"
)
shapiro_burst_34 <- safe_shapiro_test(
  simplified_df$mean_percentage[simplified_df$temp_target == 34 &
    simplified_df$object_class == "burst" &
    simplified_df$time == 82],
  label = "burst mean_percentage at 34C, time 82"
)
shapiro_burst_26
shapiro_burst_34

# Data is non-normal so Wilcoxon rank-sum test
test_burst <- wilcox.test(mean_percentage ~ temp_target, 
                          data = simplified_df %>% filter(object_class == "burst", time == 82))
test_burst

# Also comparing variance
df_burst_variability <- simplified_df %>%
  filter(object_class == "burst", time == 82) %>%
  group_by(temp_target) %>%
  summarize(
    variance = var(mean_percentage, na.rm = TRUE),
    sd = sd(mean_percentage, na.rm = TRUE),
    mean = mean(mean_percentage, na.rm = TRUE),
    cv = sd / mean  # Coefficient of Variation
  )

# Perform Levene’s test to compare variance
levene_burst <- leveneTest(mean_percentage ~ as.factor(temp_target), 
                           data = simplified_df %>% filter(object_class == "burst", time == 82))
levene_burst

plot_burst_pct_2h <- ggplot(simplified_df %>% filter(object_class == "burst", time == 82),
       aes(x = as.factor(temp_target), y = mean_percentage, fill = as.factor(temp_target))) +
  geom_boxplot(color = "black", linewidth = 1, alpha = 0.7) +
  scale_fill_manual(
    values = temp_target_fill,
    name = "Temperature (°C)",
    breaks = c("26", "34"),
    labels = c("26 °C", "34 °C")
  ) +
  labs(title = "Burst percentage at 2 hours",
       x = "Temperature (°C)",
       y = "Burst percentage") +
  scale_y_continuous(
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = c("0%", "25%", "50%", "75%", "100%"),
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 26, face = "bold"),
    axis.title.x = element_text(size = 26, face = "bold", margin = margin(t = 12)),
    axis.text = element_text(size = 22, face = "bold", color = "black"),
    axis.text.x = element_text(size = 22, face = "bold", color = "black"),
    plot.title = element_text(size = 28, face = "bold", margin = margin(0, 0, 10, 0)),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 1, color = "black"),
    axis.ticks = element_line(linewidth = 1, color = "black"),
    axis.ticks.length = unit(8, "pt"),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
    panel.grid = element_blank(),
    legend.position = "none"
  )

ggsave(
  filename = file.path(plot_dir_classes, "burst_percentage_2h_26_vs_34.png"),
  plot = plot_burst_pct_2h,
  device = "png",
  width = 7,
  height = 8,
  dpi = 400,
  units = "in"
)



