# Introduction ------------------------------------------------------------
# Comparing model predictions with ground truth hand counts.
#
# This script compares model predictions with ground truth hand counts for
# each type of model training set (all using Centernet Hourglass104):
#
# 1) Both cameras, all classes
# 2) Both cameras, pollen classes
# 3) Both cameras, tube tip class
# 4) First camera, all classes
# 5) First camera, pollen classes
# 6) First camera, tube tip classes
# 7) Second camera, all classes
# 8) Second camera, pollen classes
# 9) Second camera, tube tip classes
#
# The script calculates the optimal confidence score cutoff for each class
# in each model using the inference on the training set, then compares the 
# validation set predictions to the hand-counted ground truth. It visualizes
# these data so that the most accurate model can be determined.
#
# Manuscript version: run from the repository root. Tabular inputs live under
# data/; figures are written to plots/.

library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(patchwork)

data_root <- file.path(getwd(), "data")
plot_root <- file.path(getwd(), "plots")
plot_dir_threshold <- file.path(plot_root, "confidence_threshold_optimization")
plot_dir_lm <- file.path(plot_root, "ground_truth_vs_model_predictions")
plot_dir_processed <- file.path(plot_root, "ground_truth_vs_processed_tracks")
dir.create(plot_dir_threshold, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir_lm, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir_processed, recursive = TRUE, showWarnings = FALSE)

class_colors <- c(
  burst = "#DC267F",
  germinated = "#5fc77b",
  ungerminated = "#2F69FF",
  unknown_germinated = "#a8ffe1",
  aborted = "#787878",
  tube_tip_burst = "#ffa6db",
  tube_tip_bulging = "#fffa70",
  tube_tip = "#FFB000"
)
class_scale_breaks <- c(
  "germinated", "ungerminated", "burst", "aborted", "unknown_germinated",
  "tube_tip", "tube_tip_burst", "tube_tip_bulging"
)
class_scale_labels <- c(
  "Germinated", "Ungerminated", "Burst", "Aborted", "Unknown germinated",
  "Tube tip", "Tube tip burst", "Tube tip bulging"
)
scale_fill_classes <- function() {
  scale_fill_manual(
    values = class_colors,
    name = "Class",
    breaks = class_scale_breaks,
    labels = class_scale_labels,
    limits = force
  )
}
scale_color_classes <- function() {
  scale_color_manual(
    values = class_colors,
    name = "Class",
    breaks = class_scale_breaks,
    labels = class_scale_labels,
    limits = force
  )
}
theme_lm_scatter <- function(title_size = 22) {
  theme_bw() +
    theme(
      axis.title = element_text(size = 20, face = "bold"),
      axis.text = element_text(size = 14, face = "bold", color = "black"),
      plot.title = element_text(size = title_size, face = "bold", margin = margin(0, 0, 10, 0), hjust = 0.5),
      panel.border = element_blank(),
      axis.line = element_line(size = 1, color = "black"),
      axis.ticks = element_line(size = 1, color = "black"),
      axis.ticks.length = unit(8, "pt"),
      plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text = element_text(size = 12, face = "bold", color = "black"),
      legend.position = "none",
      aspect.ratio = 1
    )
}

# Importing ground truth and model predictions ----------------------------
import_ground_truth <- function() {
  gt_dir <- file.path(data_root, "model_ground_truth")
  cameras <- c("combined", "one", "two")
  splits <- c("train", "val")
  pieces <- vector("list", length(cameras) * length(splits))
  k <- 0L
  for (cam in cameras) {
    for (split in splits) {
      k <- k + 1L
      fn <- sprintf("ground_truth_%s_all_%s.tsv", cam, split)
      df <- read.table(file.path(gt_dir, fn), sep = "\t", header = TRUE)
      df$camera <- cam
      df$train_or_val <- split
      pieces[[k]] <- df
    }
  }
  output_df <- bind_rows(pieces)
  output_df[, -1]
}

ground_truth <- import_ground_truth()

import_model_predictions <- function() {
  empty <- data.frame(
    date = character(),
    run = character(),
    well = character(),
    timepoint = numeric(),
    tempc = numeric(),
    class = character(),
    score = numeric(),
    ymin = numeric(),
    xmin = numeric(),
    ymax = numeric(),
    xmax = numeric(),
    train_camera = character(),
    label_type = character(),
    train_or_val = character(),
    inference_camera = character(),
    stringsAsFactors = FALSE
  )
  train_cameras <- c("combined", "one", "two")
  label_types <- c("all", "pollen", "tube_tip")
  train_or_vals <- c("train", "val")
  inference_cameras <- c("one", "two")
  pred_dir <- file.path(data_root, "model_predictions")
  pieces <- list()
  for (train_camera in train_cameras) {
    for (label_type in label_types) {
      for (train_or_val in train_or_vals) {
        for (inference_camera in inference_cameras) {
          filename <- file.path(
            pred_dir,
            paste0(
              train_camera, "_", label_type, "_", train_or_val, "_",
              inference_camera, "_predictions.tsv"
            )
          )
          if (file.exists(filename)) {
            temp_df <- read.table(filename, sep = "\t", header = TRUE)
            temp_df$train_camera <- train_camera
            temp_df$label_type <- label_type
            temp_df$train_or_val <- train_or_val
            temp_df$inference_camera <- inference_camera
            pieces[[length(pieces) + 1L]] <- temp_df
          }
        }
      }
    }
  }
  if (!length(pieces)) {
    return(empty)
  }
  bind_rows(pieces)
}

model_predictions <- import_model_predictions()


# Calculate optimal confidence thresholds ---------------------------------
# Each predicted bounding box has a confidence score. I want to figure out the 
# best confidence score threshold to use for each model/class. I will calculate
# these optimal thresholds using the model inference on the training datasets.
# Optimal confidence score thresholds will be chosen for each class by first 
# calculating the r-squared values at each confidence threshold, then choosing 
# the highest one.

calculate_optimal_confidence_threshold <- function(ground_truth, model_predictions) {
  
  # This function calculates the r-squared values comparing ground truth and 
  # model predictions at set confidence score cutoffs. This function is run 
  # once for each camera (combined, one, two) and class (all, pollen, tube_tip) 
  # combination using the training set inference data. For the combined training 
  # set, r-squared values are calculated for test sets from camera one or two. 
  calculate_r_squared <- function(
    ground_truth, 
    inference, 
    train_camera_string, 
    label_type_string, 
    train_or_val_string,
    inference_camera_string) {
    
    # First subsetting the inference to only include the requested camera and 
    # training class.
    inference <- inference %>%
      filter(train_camera == train_camera_string & label_type == label_type_string & train_or_val == train_or_val_string & inference_camera == inference_camera_string)
    
    # Next subsetting and pre-processing the ground truth
    process_ground_truth <- function(df, inference_camera_string, train_or_val_string) {
      df <- df %>%
        filter(camera == inference_camera_string & train_or_val == train_or_val_string) %>%
        complete(name, class) %>%
        mutate(hand_count = replace_na(size, 0)) %>%
        select(-size)
      return(df)
    }
    
    ground_truth <- process_ground_truth(ground_truth, inference_camera_string, train_or_val_string)
    
    # Sub function that calculates the r-squared for a single threshold and 
    # a single class. 
    calculate_single_r_squared <- function(ground_truth, inference, threshold, class_string) {
      # Processing the inference with confidence score cutoff (threshold)
      process_inference <- function(df, confidence_cutoff){
        df <- df %>%
          mutate(name = paste0(
            date,
            "_run",
            run,
            "_",
            tempc,
            "C_",
            well,
            "_t",
            str_pad(timepoint, 3, pad = "0")
          )) %>%
          filter(score >= confidence_cutoff) %>%
          group_by(name, class) %>%
          summarize(model_count = n(), .groups = "drop") %>%
          complete(name, class) %>%
          mutate(model_count = replace_na(model_count, 0)) 
        return(df)
      }
      
      inference <- process_inference(inference, threshold)
      
      # Combine the two data frames
      df <- full_join(ground_truth, inference, by = c("name", "class"))
      # df <- left_join(inference, ground_truth, by = c("name", "class"))
      
      df <- df %>%
        filter(class == class_string) %>%
        ungroup() %>%
        complete(name, class) %>%
        mutate(model_count = replace_na(model_count, 0))
      
      # # Test print
      # if (threshold == 0.2) {
      #   print("R-squared df at 0.2:")
      #   View(df)
      # }
      
      # Doing the regression
      regression_r <- summary(lm(hand_count ~ model_count, data = df))$adj.r.squared
      output <- data.frame("class" = class_string, 
                           "threshold" = threshold, 
                           "r_squared" = regression_r)
      return(output)
    }
    
    # Going through each class and calculating the r-squared values at intervals of 0.01
    output_df <- data.frame()
    
    for (class_string in unique(inference$class)) {
      thresholds <- seq(0.01, 1, by = 0.01)
      r_squared_df <- bind_rows(lapply(thresholds, function(x) {
        calculate_single_r_squared(ground_truth, inference, x, class_string)
      }))
      output_df <- rbind(output_df, r_squared_df)
    }
    
    output_df <- output_df %>%
      filter(r_squared != 0)
    
    output_df$train_camera <- train_camera_string
    output_df$label_type <- label_type_string
    output_df$train_or_val <- train_or_val_string
    output_df$inference_camera <- inference_camera_string
    
    # Print the max R-squared values for each class
    print(paste0("Max R-squared values for train camera ", 
                 train_camera_string, 
                 " and inference camera ", 
                 inference_camera_string,
                 " and train_or_val_string ",
                 train_or_val_string,
                 " and label type ",
                 label_type_string))
    
    for (class_string in unique(output_df$class)) {
      print_df <- output_df %>%
        filter(class == class_string)
      r_squared_val = as.character(max(print_df$r_squared))
      threshold_at_max_r_squared = as.character(print_df$threshold[which.max(print_df$r_squared)])
      print(paste0(class_string, ": ", r_squared_val, " at threshold ", threshold_at_max_r_squared))
    }
    
    return(output_df)
  } 
  
  # Going through all the different combinations and calculating the r-squareds.
  output_df <- data.frame(
    class = character(),
    threshold = numeric(),
    r_squared = numeric(),
    train_camera = character(),
    label_type = character(),
    train_or_val = character(),
    inference_camera = character()
  )
  
  train_cameras <- c("combined", "one", "two")
  label_types <- c("all", "pollen", "tube_tip")
  train_or_vals <- c("train", "val")
  inference_cameras <- c("one", "two")
  
  # # Test case
  # r_squareds <- calculate_r_squared(ground_truth,
  #                                   model_predictions,
  #                                   "combined",
  #                                   "all",
  #                                   "train",
  #                                   "one")
  # output_df <- rbind(output_df, r_squareds)
  
  for (train_camera in train_cameras) {
    for (label_type in label_types) {
      for (train_or_val in train_or_vals) {
        for (inference_camera in inference_cameras) {
          # Skipping the parts of the loop that don't exist
          if ((train_camera == "one" & inference_camera == "two") |
              (train_camera == "two" & inference_camera == "one")) {
            next
          }

          r_squareds <- calculate_r_squared(ground_truth,
                                            model_predictions,
                                            train_camera,
                                            label_type,
                                            train_or_val,
                                            inference_camera)
          output_df <- rbind(output_df, r_squareds)
        }
      }
    }
  }
  
  return(output_df)
}

r_squared_df <- calculate_optimal_confidence_threshold(ground_truth, model_predictions)

# Summary of highest r-squared values
r_squared_summary <- r_squared_df %>%
  group_by(train_camera, inference_camera, label_type, train_or_val, class) %>%
  summarize(max_r_squared = max(r_squared), threshold_at_max = threshold[which.max(r_squared)])

r_squared_summary_train <- r_squared_summary %>%
  filter(train_or_val == "train")

# Plot r-squared values at all confidence thresholds ----------------------
plot_threshold_facet <- function(df, label_type_string, train_or_val_string) {
  # Subsetting the dataframe for the camera and train/val set
  df <- df %>%
    filter(
      train_or_val == train_or_val_string &
      label_type == label_type_string
    ) %>%
    mutate(name = paste0(train_camera, " training, ", inference_camera, " inference"))
  
  plot_aesthetics <- theme_bw() +
    theme(
      axis.title = element_text(size = 16, face = 'bold'),
      axis.text = element_text(size = 13, face = 'bold', color = 'black'),
      plot.title = element_text(size = 20, face = 'bold', margin = margin(0, 0, 10, 0), hjust = 0.5),
      axis.title.x = element_text(margin = margin(10, 0, 0, 0)),
      panel.border = element_blank(),
      axis.line = element_line(linewidth = 1, color = 'black'),
      axis.ticks = element_line(linewidth = 1, color = 'black'),
      axis.ticks.length = unit(8, 'pt'),
      plot.margin = margin(0.5, 0.5, 0.5, 0.5, 'cm'),
      panel.grid = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 11, face = 'bold'),
      legend.position = 'right',
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text = element_text(size = 14, face = 'bold')
    )
  
  x_axis_settings <- scale_x_continuous(
    breaks = seq(0, 1, 0.2),
    labels = seq(0, 1, 0.2),
    limits = c(0, 1.05),
    expand = c(0, 0)
  ) 
  y_axis_settings <- scale_y_continuous(
    breaks = seq(0, 1, 0.2),
    labels = seq(0, 1, 0.2),
    limits = c(0, 1.05),
    expand = c(0, 0)
  )
  
  if (train_or_val_string == "train") {
    plot_title_string = "Train"
  } else if (train_or_val_string == "val") {
    plot_title_string = "Validation"
  }
  
  plot_all <- ggplot(df, aes(x = threshold, y = r_squared, color = class)) +
    geom_point(size = 1.2) +
    # geom_smooth(linewidth = 1, span = 0.2, se = FALSE) +
    facet_wrap(vars(name), ncol = 1, scales = "free") +
    x_axis_settings +
    y_axis_settings +
    scale_color_classes() +
    labs(title = plot_title_string, 
         x = "Confidence threshold", y = "R-squared") +
    plot_aesthetics
  
  
  if (train_or_val_string == "train") {
    plot_all <- plot_all + theme(legend.position = "none")
  }
  
  return(plot_all)
  
}

all_labels_train_plot <- plot_threshold_facet(r_squared_df, "all", "train")
all_labels_val_plot <- plot_threshold_facet(r_squared_df, "all", "val")

pollen_labels_train_plot <- plot_threshold_facet(r_squared_df, "pollen", "train")
pollen_labels_val_plot <- plot_threshold_facet(r_squared_df, "pollen", "val")

tube_tip_labels_train_plot <- plot_threshold_facet(r_squared_df, "tube_tip", "train")
tube_tip_labels_val_plot <- plot_threshold_facet(r_squared_df, "tube_tip", "val")

all_threshold_plot_title <- plot_annotation(title = "All class confidence thresholds",
  theme = theme(plot.title = element_text(size = 28, face = 'bold', hjust = 0.5, margin = margin(10, 0, 10, 0))))
pollen_threshold_plot_title <- plot_annotation(title = "Pollen class confidence thresholds",
  theme = theme(plot.title = element_text(size = 28, face = 'bold', hjust = 0.5, margin = margin(10, 0, 10, 0))))
tube_tip_threshold_plot_title <- plot_annotation(title = "Tube tip class confidence thresholds",
  theme = theme(plot.title = element_text(size = 28, face = 'bold', hjust = 0.5, margin = margin(10, 0, 10, 0))))

all_threshold_plot <- (all_labels_train_plot | all_labels_val_plot) + all_threshold_plot_title 
pollen_threshold_plot <- (pollen_labels_train_plot | pollen_labels_val_plot) + pollen_threshold_plot_title 
tube_tip_threshold_plot <- (tube_tip_labels_train_plot | tube_tip_labels_val_plot) + tube_tip_threshold_plot_title 

ggsave(all_threshold_plot,
       filename = file.path(plot_dir_threshold, "threshold_plot_all.png"),
       device = 'png',
       width = 11,
       height = 12,
       dpi = 400,
       units = 'in')
ggsave(pollen_threshold_plot,
       filename = file.path(plot_dir_threshold, "threshold_plot_pollen.png"),
       device = 'png',
       width = 11,
       height = 12,
       dpi = 400,
       units = 'in')
ggsave(tube_tip_threshold_plot,
       filename = file.path(plot_dir_threshold, "threshold_plot_tube_tip.png"),
       device = 'png',
       width = 10,
       height = 12,
       dpi = 400,
       units = 'in')


# Make max r-squared threshold bar plot -----------------------------------
# It's a little hard to make sense of the threshold plots from the previous 
# section, so I'll also make a little bar plot that summarizes the max 
# r-squared values from the different models on the training set.
plot_r_squared_cols <- function(df, camera_string, title_string, zoom) {
  df <- df %>%
    filter(inference_camera == camera_string, class != "unknown_germinated") %>%
    mutate(name = paste(train_camera, inference_camera, label_type))
  
  df$class <- factor(df$class, levels = c("ungerminated",
                                          "germinated",
                                          "burst",
                                          "aborted",
                                          "tube_tip"))
  
  y_axis_settings <- scale_y_continuous(
    breaks = seq(0, 1, 0.2),
    labels = seq(0, 1, 0.2),
    limits = c(0, 1),
    expand = c(0, 0)
  )
  
  plot_aesthetics <- theme_bw() +
    theme(
      axis.title = element_text(size = 16, face = 'bold'),
      axis.text = element_text(size = 13, face = 'bold', color = 'black'),
      plot.title = element_text(size = 20, face = 'bold', margin = margin(0, 0, 10, 0), hjust = 0.5),
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust=1),
      panel.border = element_blank(),
      axis.line = element_line(linewidth = 1, color = 'black'),
      axis.ticks = element_line(linewidth = 1, color = 'black'),
      axis.ticks.length = unit(8, 'pt'),
      plot.margin = margin(0.5, 0.5, 0.5, 0.5, 'cm'),
      # panel.grid = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(linewidth = 1, color = "black"),
      panel.grid.minor.y = element_line(linewidth = 0.5, color = "black"),
      legend.title = element_blank(),
      legend.text = element_text(size = 16, face = 'bold'),
      legend.position = 'bottom',
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text = element_text(size = 14, face = 'bold')
    )
  
  output_plot <- ggplot(df, aes(x = name, y = max_r_squared, fill = class)) +
    geom_col(width = 0.8, color = "black", position = position_dodge2(width = 1, preserve = "single")) +
    y_axis_settings +
    scale_fill_classes() +
    labs(title = title_string, 
         x = "Name", y = "r-squared") +
    plot_aesthetics 
  
  if (zoom) {
    output_plot <- output_plot + scale_y_continuous(
        breaks = seq(0, 1, 0.02),
        labels = seq(0, 1, 0.02),
        limits = c(0, 1),
        expand = c(0, 0)
      ) +
      coord_cartesian(ylim = c(0.86, 1))
  }
  
  return(output_plot)
  
}

# Also add the max validation r-squared values at the threshold determined in the training set
r_squared_col_plot_one <- plot_r_squared_cols(df = r_squared_summary_train, camera_string = "one", title_string = "Max training r-squared values camera one", zoom = FALSE)
r_squared_col_plot_one_zoom <- plot_r_squared_cols(df = r_squared_summary_train, camera_string = "one", title_string = "Max training r-squared values camera one", zoom = TRUE)
r_squared_col_plot_two <- plot_r_squared_cols(df = r_squared_summary_train, camera_string = "two", title_string = "Max training r-squared values camera two", zoom = FALSE)
r_squared_col_plot_two_zoom <- plot_r_squared_cols(df = r_squared_summary_train, camera_string = "two", title_string = "Max training r-squared values camera two", zoom = TRUE)

ggsave(r_squared_col_plot_one,
       filename = file.path(plot_dir_threshold, "max_r_squared_plot_camera_one.png"),
       device = 'png',
       width = 10,
       height = 10,
       dpi = 400,
       units = 'in')
ggsave(r_squared_col_plot_one_zoom,
       filename = file.path(plot_dir_threshold, "max_r_squared_plot_camera_one_zoom.png"),
       device = 'png',
       width = 10,
       height = 10,
       dpi = 400,
       units = 'in')
ggsave(r_squared_col_plot_two,
       filename = file.path(plot_dir_threshold, "max_r_squared_plot_camera_two.png"),
       device = 'png',
       width = 10,
       height = 10,
       dpi = 400,
       units = 'in')
ggsave(r_squared_col_plot_two_zoom,
       filename = file.path(plot_dir_threshold, "max_r_squared_plot_camera_two_zoom.png"),
       device = 'png',
       width = 10,
       height = 10,
       dpi = 400,
       units = 'in')


# Make validation set r-squared bar plot ----------------------------------
# I'd also like to visualize the r-squared values from the different models
# in the validation set at the threshold cutoff determined by the max r-squared 
# value in the training set. I'll make the data frame using a join. 

train_join_df <- r_squared_summary %>%
  filter(train_or_val == "train", class != "unknown_germinated") %>%
  mutate(threshold = threshold_at_max, train_r_squared = max_r_squared)

val_join_df <- r_squared_df %>%
  filter(train_or_val == "val", class != "unknown_germinated") %>%
  mutate(val_r_squared = r_squared)

r_squared_val_at_train_r <- left_join(train_join_df, val_join_df, by = c("train_camera", "inference_camera", "label_type", "class", "threshold"))
r_squared_val_at_train_r <- r_squared_val_at_train_r  %>%
  mutate(train_or_val = train_or_val.y, max_r_squared = val_r_squared)

r_squared_col_plot_val_one <- plot_r_squared_cols(df = r_squared_val_at_train_r , camera_string = "one", title_string = "Validation r-squared values camera one", zoom = FALSE)
r_squared_col_plot_val_one_zoom <- plot_r_squared_cols(df = r_squared_val_at_train_r , camera_string = "one", title_string = "Validation r-squared values camera one", zoom = TRUE)
r_squared_col_plot_val_two <- plot_r_squared_cols(df = r_squared_val_at_train_r , camera_string = "two", title_string = "Validation r-squared values camera two", zoom = FALSE)
r_squared_col_plot_val_two_zoom <- plot_r_squared_cols(df = r_squared_val_at_train_r , camera_string = "two", title_string = "Validation r-squared values camera two", zoom = TRUE)

ggsave(r_squared_col_plot_val_one,
       filename = file.path(plot_dir_threshold, "max_r_squared_plot_val_one.png"),
       device = 'png',
       width = 10,
       height = 10,
       dpi = 400,
       units = 'in')
ggsave(r_squared_col_plot_val_one_zoom,
       filename = file.path(plot_dir_threshold, "max_r_squared_plot_val_one_zoom.png"),
       device = 'png',
       width = 10,
       height = 10,
       dpi = 400,
       units = 'in')
ggsave(r_squared_col_plot_val_two,
       filename = file.path(plot_dir_threshold, "max_r_squared_plot_val_two.png"),
       device = 'png',
       width = 10,
       height = 10,
       dpi = 400,
       units = 'in')
ggsave(r_squared_col_plot_val_two_zoom,
       filename = file.path(plot_dir_threshold, "max_r_squared_plot_val_two_zoom.png"),
       device = 'png',
       width = 10,
       height = 10,
       dpi = 400,
       units = 'in')


# Subset validation set predictions with confidence thresholds ------------
# Making different subsets of the model predictions for the validation set 
# based on the max confidence threshold chosen in the training set. 
model_predictions_val <- model_predictions %>%
  filter(train_or_val == "val", class != "unknown_germinated") %>%
  mutate(name = paste0(
    date,
    "_run",
    run,
    "_",
    tempc,
    "C_",
    well,
    "_t",
    str_pad(timepoint, 3, pad = "0")
  )) %>%
  # This strategy uses a join then a filter to only keep the rows that pass the 
  # confidence threshold cutoff determined in the training dataset.
  left_join(r_squared_summary_train, by = c(
    "train_camera", "inference_camera", "label_type", "class"
  )) %>%
  filter(score >= threshold_at_max) 


# Compare ground truth with model predictions -----------------------------
plot_lm_facet <- function(
  ground_truth, 
  inference, 
  train_camera_string, 
  inference_camera_string,
  label_type_string,
  title_string) {
  
  # Processing the ground truth
  process_ground_truth <- function(df){
    df <- df %>%
      complete(name, class) %>%
      mutate(hand_count = replace_na(size, 0)) %>%
      select(-size)
    return(df)
  }
  ground_truth <- process_ground_truth(ground_truth)
  
  # Processing the inference
  process_inference <- function(df, train_camera_string, inference_camera_string, label_type_string){
    df <- df %>%
      filter(
        train_camera == train_camera_string,
        inference_camera == inference_camera_string,
        label_type == label_type_string) %>%
      group_by(name, class) %>%
      summarize(model_count = n(), .groups = "drop") %>%
      complete(name, class) %>%
      mutate(model_count = replace_na(model_count, 0)) 
    
    return(df)
  }
  
  inference <- process_inference(inference, train_camera_string, inference_camera_string, label_type_string)
  
  # Combining the two data frames
  df <- full_join(ground_truth, inference, by = c("name", "class"))
  
  df <- df %>%
    ungroup() %>%
    complete(name, class) %>%
    mutate(model_count = replace_na(model_count, 0)) %>%
    filter(model_count != 0)
  
  # Plotting
  out_plot <- ggplot(df, aes(x = hand_count, y = model_count, fill = class)) +
    geom_abline(intercept = 0, slope = 1, linewidth = 1, linetype = 2) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1, color = "black") +
    geom_point(shape = 21, color = "black", size = 2) +
    scale_fill_classes() +
    facet_wrap(~class, nrow = 1, scales = "free") +
    expand_limits(x = max(c(df$hand_count, df$model_count)), 
                  y = max(c(df$hand_count, df$model_count))) +
    labs(title = title_string, x = "Hand counts", y = "Model predictions") +
    theme_lm_scatter(22)
  
  ggsave(out_plot,
         filename = file.path(plot_dir_lm, paste0(train_camera_string, "_",
                                                  inference_camera_string, "_",
                                                  label_type_string, ".png")),
         device = 'png',
         width = 14,
         height = 5,
         dpi = 400,
         units = 'in')
  
  return(out_plot)
  
  
}

lm_facet_combined_one_all <- plot_lm_facet(ground_truth, model_predictions_val, "combined", "one", "all", "Combined camera train, camera one validation, all classes")
lm_facet_combined_one_pollen <- plot_lm_facet(ground_truth, model_predictions_val, "combined", "one", "pollen", "Combined camera train, camera one validation, pollen classes")
lm_facet_combined_one_tube_tip <- plot_lm_facet(ground_truth, model_predictions_val, "combined", "one", "tube_tip", "Combined camera train, camera one validation, tube tip class")

# lm_facet_combined_one_all / lm_facet_combined_one_pollen / lm_facet_combined_one_tube_tip

lm_facet_combined_two_all <- plot_lm_facet(ground_truth, model_predictions_val, "combined", "two", "all", "Combined camera train, camera two validation, all classes")
lm_facet_combined_two_pollen <- plot_lm_facet(ground_truth, model_predictions_val, "combined", "two", "pollen","Combined camera train, camera two validation, pollen classes")
lm_facet_combined_two_tube_tip <- plot_lm_facet(ground_truth, model_predictions_val, "combined", "two", "tube_tip", "Combined camera train, camera two validation, tube tip class")

lm_facet_one_one_all <- plot_lm_facet(ground_truth, model_predictions_val, "one", "one", "all", "Camera one train and validation, all classes")
lm_facet_one_one_pollen <- plot_lm_facet(ground_truth, model_predictions_val, "one", "one", "pollen", "Camera one train and validation, pollen classes")
lm_facet_one_one_tube_tip <- plot_lm_facet(ground_truth, model_predictions_val, "one", "one", "tube_tip", "Camera one train and validation, tube tip class")

lm_facet_two_two_all <- plot_lm_facet(ground_truth, model_predictions_val, "two", "two", "all", "Camera two train and validation, all classes")
lm_facet_two_two_pollen <- plot_lm_facet(ground_truth, model_predictions_val, "two", "two", "pollen", "Camera two train and validation, pollen classes")
lm_facet_two_two_tube_tip <- plot_lm_facet(ground_truth, model_predictions_val, "two", "two", "tube_tip", "Camera two train and validation, tube tip class")


# Compare ground truth with processed tracks ------------------------------
# These predictions are from the processed tracks, including Bayesian Tracker
# and biological prior-based class inference from track data.
tracks_path <- file.path(data_root, "model_predictions", "2023-06-23_all_tracks_bug_fix.txt")

# Compare processed predictions and ground truth. This one is simplified a bit 
# because based on the previous step we selected the model trained with both 
# cameras and all classes.
plot_lm_facet_processed <- function(
    ground_truth, 
    inference, 
    inference_camera_string,
    title_string) {
  
  # Processing the ground truth
  process_ground_truth <- function(df){
    df <- df %>%
      filter(
        train_or_val == "val", 
        camera == "combined", 
        class != "unknown_germinated"
      ) %>%
      complete(name, class) %>%
      mutate(hand_count = replace_na(size, 0)) %>%
      select(-size)
    return(df)
  }
  ground_truth <- process_ground_truth(ground_truth)
  if (inference_camera_string == "one"){
    ground_truth <- ground_truth[1:500, ]
  }
  if (inference_camera_string == "two"){
    ground_truth <- ground_truth[501:1000, ]
  }
  
  # Processing the inference
  process_inference <- function(df){
    df <- df %>%
      group_by(name, object_class) %>%
      summarize(model_count = n(), .groups = "drop") %>%
      complete(name, object_class) %>%
      mutate(model_count = replace_na(model_count, 0)) %>%
      rename(class = object_class)
    
    return(df)
  }
  if (inference_camera_string == "one"){
    inference <- inference[as.Date(inference$date) < as.Date("2022-05-27"), ]
  }
  if (inference_camera_string == "two"){
    inference <- inference[as.Date(inference$date) >= as.Date("2022-05-27"), ]
  }
  inference <- process_inference(inference)
  
  # Combining the two data frames
  df <- full_join(ground_truth, inference, by = c("name", "class"))
  
  df <- df %>%
    ungroup() %>%
    complete(name, class) %>%
    mutate(model_count = replace_na(model_count, 0)) %>%
    filter(model_count != 0)
  
  # Getting r-squared values
  for (class_string in c("aborted", "ungerminated", "germinated", "burst", "tube_tip")) {
    print(class_string)
    subsetted_df <- df %>%
      filter(class == class_string)
    print(summary(lm(hand_count ~ model_count, data = subsetted_df))$adj.r.squared)
  }
  
  # Plotting
  processed_png <- if (inference_camera_string == "one") {
    "processed_tracks_cam_one.png"
  } else {
    "processed_tracks_cam_two.png"
  }
  out_plot <- ggplot(df, aes(x = hand_count, y = model_count, fill = class)) +
    geom_abline(intercept = 0, slope = 1, linewidth = 1, linetype = 2) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1, color = "black") +
    geom_point(shape = 21, color = "black", size = 2) +
    scale_fill_classes() +
    facet_wrap(~class, nrow = 1, scales = "free") +
    expand_limits(x = max(c(df$hand_count, df$model_count)), 
                  y = max(c(df$hand_count, df$model_count))) +
    labs(title = title_string, x = "Hand counts", y = "Model predictions") +
    theme_lm_scatter(22)
  
  ggsave(out_plot,
         filename = file.path(plot_dir_processed, processed_png),
         device = 'png',
         width = 14,
         height = 5,
         dpi = 400,
         units = 'in')
  
  return(df)
}

# Compare processed predictions and ground truth for individual classes with 
# combined classes (for Plant Bio Savannah talk).
plot_individual_lm_processed <- function(
    ground_truth, 
    inference, 
    class_string,
    title_string) {
  
  process_ground_truth <- function(df){
    df <- df %>%
      filter(
        train_or_val == "val", 
        camera == "combined", 
        class != "unknown_germinated"
      ) %>%
      complete(name, class) %>%
      mutate(hand_count = replace_na(size, 0)) %>%
      select(-size)
    return(df)
  }
  ground_truth <- process_ground_truth(ground_truth)
  
  process_inference <- function(df){
    df <- df %>%
      group_by(name, object_class) %>%
      summarize(model_count = n(), .groups = "drop") %>%
      complete(name, object_class) %>%
      mutate(model_count = replace_na(model_count, 0)) %>%
      rename(class = object_class)
    
    return(df)
  }
  inference <- process_inference(inference)
  
  df <- full_join(ground_truth, inference, by = c("name", "class"))
  
  df <- df %>%
    ungroup() %>%
    complete(name, class) %>%
    mutate(model_count = replace_na(model_count, 0)) %>%
    filter(model_count != 0)
  
  df <- df %>%
    filter(class == class_string)
  
  print(paste0(
    class_string, 
    ": ",
    summary(lm(hand_count ~ model_count, data = df))$adj.r.squared)
  )
  
  axis_limit <- max(c(max(df$hand_count), max(df$model_count))) + 5
  
  out_plot <- ggplot(df, aes(x = hand_count, y = model_count, fill = class)) +
    geom_abline(intercept = 0, slope = 1, linewidth = 1, linetype = 2) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1, color = "black") +
    geom_point(shape = 21, color = "black", size = 3, alpha = 0.7) +
    scale_x_continuous(limits = c(0, axis_limit)) +
    scale_y_continuous(limits = c(0, axis_limit)) +
    scale_fill_classes() +
    labs(title = title_string, x = "Hand counts", y = "Model predictions") +
    theme_lm_scatter(24)
  
  ggsave(out_plot,
         filename = file.path(
           plot_dir_processed,
           paste0("individual_lm_plot_both_cams_", title_string, ".png")),
         device = 'png',
         width = 6,
         height = 6,
         dpi = 400,
         units = 'in')
}

if (!file.exists(tracks_path)) {
  warning("Skipping processed-tracks analysis (file not found): ", tracks_path)
} else {
  ground_truth_images <- unique(ground_truth[ground_truth$camera == "combined" & ground_truth$train_or_val == "val", ]$name)
  tracks_df <- read.table(tracks_path, sep = "\t", header = TRUE)
  tracks_df <- tracks_df %>%
    mutate(name = paste0(
      date,
      "_run",
      run,
      "_",
      tempc,
      "C_",
      well,
      "_t",
      str_pad(t, 3, pad = "0")
    )) %>%
    filter(name %in% ground_truth_images)

lm_facet_processed_cam_one <- plot_lm_facet_processed(
  ground_truth, 
  tracks_df,
  "one",
  "Processed tracks camera one bug fix")

lm_facet_processed_cam_two <- plot_lm_facet_processed(
  ground_truth, 
  tracks_df,
  "two",
  "Processed tracks camera two bug fix")

plot_individual_lm_processed(ground_truth, tracks_df, "ungerminated", "Ungerminated")
plot_individual_lm_processed(ground_truth, tracks_df, "germinated", "Germinated")
plot_individual_lm_processed(ground_truth, tracks_df, "burst", "Burst")
plot_individual_lm_processed(ground_truth, tracks_df, "aborted", "Aborted")
plot_individual_lm_processed(ground_truth, tracks_df, "tube_tip", "Tube tip")

# Pulling out some specific images to check -------------------------------
image_test <- lm_facet_processed_cam_one[lm_facet_processed_cam_one$class == "ungerminated", ]
image_test$diff <- image_test$hand_count - image_test$model_count
}
