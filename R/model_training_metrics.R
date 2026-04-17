# Introduction ------------------------------------------------------------
# Manuscript version: model training metrics supplemental figure.
# Run from repository root. Inputs: data/; figures are written
# to plots/model_training_metrics/.
#
# Data inputs:
#   data/model_training_tensorboard/*.csv  (Plot B)
#   data/model_ground_truth/*.tsv          (Plots A, C)
#   data/model_predictions/*.tsv           (Plots A, C)
#
# Plots produced:
#   A) Individual vs Combined Camera Training.png
#   B) Model Training Metrics.png
#   C) CenterNet HourGlass All Classes.png

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(stringr)

pdf(NULL)
set.seed(16)

data_root <- file.path(getwd(), "data")
plot_dir  <- file.path(getwd(), "plots", "model_training_metrics")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)


# Shared aesthetics -------------------------------------------------------
manuscript_theme <- function(title_size = 22, axis_title_size = 18, axis_text_size = 14) {
  theme_bw() +
    theme(
      axis.title        = element_text(size = axis_title_size, face = "bold"),
      axis.text         = element_text(size = axis_text_size, face = "bold", color = "black"),
      plot.title        = element_text(size = title_size, face = "bold", hjust = 0.5,
                                       margin = margin(0, 0, 10, 0)),
      panel.border      = element_blank(),
      axis.line         = element_line(linewidth = 1, color = "black"),
      axis.ticks        = element_line(linewidth = 1, color = "black"),
      axis.ticks.length = unit(8, "pt"),
      plot.margin       = margin(0.5, 0.5, 0.5, 0.5, "cm"),
      panel.grid        = element_blank()
    )
}

# Class color palette consistent with compare_predictions_with_ground_truth.R
class_colors <- c(
  burst              = "#DC267F",
  germinated         = "#5fc77b",
  ungerminated       = "#2F69FF",
  unknown_germinated = "#a8ffe1",
  aborted            = "#787878",
  tube_tip_burst     = "#ffa6db",
  tube_tip_bulging   = "#fffa70",
  tube_tip           = "#FFB000"
)
class_scale_breaks <- c(
  "germinated", "ungerminated", "burst", "aborted",
  "unknown_germinated", "tube_tip", "tube_tip_burst", "tube_tip_bulging"
)
class_scale_labels <- c(
  "Germinated", "Ungerminated", "Burst", "Aborted",
  "Unknown germinated", "Tube tip", "Tube tip burst", "Tube tip bulging"
)

scale_color_classes <- function() {
  scale_color_manual(
    values = class_colors,
    name   = "Class",
    breaks = class_scale_breaks,
    labels = class_scale_labels,
    limits = force
  )
}


# Loading the data --------------------------------------------------------

# Plot B: tensorboard training metrics data
load_tensorboard_data <- function() {
  df                 <- data.frame()
  metrics            <- c("total_loss", "AR100", "mAP")
  bounding_box_types <- c("all", "pollen", "tube_tip")
  cameras            <- c("combined", "one", "two")
  prefix             <- file.path(data_root, "model_training_tensorboard")

  for (camera in cameras) {
    for (metric in metrics) {
      for (bounding_box_type in bounding_box_types) {
        for (train_or_eval in c("train", "eval")) {
          filename <- file.path(
            prefix,
            paste0(camera, "_", bounding_box_type, "_", metric, "_", train_or_eval, ".csv")
          )
          if (file.exists(filename)) {
            temp <- read.csv(filename, header = TRUE)
            temp$camera            <- camera
            temp$metric            <- metric
            temp$bounding_box_type <- bounding_box_type
            temp$train_or_eval     <- train_or_eval
            df <- rbind(df, temp[, c(
              "camera", "metric", "bounding_box_type",
              "train_or_eval", "Wall.time", "Step", "Value"
            )])
          }
        }
      }
    }
  }
  colnames(df) <- c(
    "camera", "metric", "bounding_box_type", "train_or_eval", "time", "step", "value"
  )
  df$time <- as.POSIXct(df$time, origin = "1970-01-01")
  df
}

# Plots A & C: ground truth hand counts
import_ground_truth <- function() {
  gt_dir  <- file.path(data_root, "model_ground_truth")
  cameras <- c("combined", "one", "two")
  splits  <- c("train", "val")
  pieces  <- vector("list", length(cameras) * length(splits))
  k <- 0L
  for (cam in cameras) {
    for (split in splits) {
      k <- k + 1L
      fn  <- sprintf("ground_truth_%s_all_%s.tsv", cam, split)
      df  <- read.table(file.path(gt_dir, fn), sep = "\t", header = TRUE)
      df$camera       <- cam
      df$train_or_val <- split
      pieces[[k]] <- df
    }
  }
  bind_rows(pieces)[, -1]
}

# Plots A & C: model prediction bounding boxes
import_model_predictions <- function() {
  pred_dir       <- file.path(data_root, "model_predictions")
  train_cameras  <- c("combined", "one", "two")
  label_types    <- c("all", "pollen", "tube_tip")
  train_or_vals  <- c("train", "val")
  inference_cams <- c("one", "two")
  pieces <- list()
  for (tc in train_cameras) {
    for (lt in label_types) {
      for (tv in train_or_vals) {
        for (ic in inference_cams) {
          fn <- file.path(
            pred_dir,
            paste0(tc, "_", lt, "_", tv, "_", ic, "_predictions.tsv")
          )
          if (file.exists(fn)) {
            temp <- read.table(fn, sep = "\t", header = TRUE)
            temp$train_camera     <- tc
            temp$label_type       <- lt
            temp$train_or_val     <- tv
            temp$inference_camera <- ic
            pieces[[length(pieces) + 1L]] <- temp
          }
        }
      }
    }
  }
  bind_rows(pieces)
}

tensorboard_df    <- load_tensorboard_data()
ground_truth      <- import_ground_truth()
model_predictions <- import_model_predictions()


# Plot A: Individual vs Combined Camera Training --------------------------
# For each (train_camera, label_type, inference_camera, class) combination,
# calculate r-squared at each confidence threshold using the training set,
# then retrieve the validation r-squared at the optimal training threshold.
# The combined model (trained on both cameras) is compared against models
# trained on a single camera to show whether pooling training data improves
# prediction accuracy.

calculate_r_squared_df <- function(ground_truth, model_predictions) {

  process_gt <- function(df, ic, tv) {
    df %>%
      filter(camera == ic, train_or_val == tv) %>%
      complete(name, class) %>%
      mutate(hand_count = replace_na(size, 0)) %>%
      select(-size)
  }

  process_inference <- function(df, tc, lt, tv, ic) {
    df %>%
      filter(
        train_camera    == tc,
        label_type      == lt,
        train_or_val    == tv,
        inference_camera == ic
      ) %>%
      mutate(name = paste0(
        date, "_run", run, "_", tempc, "C_", well,
        "_t", str_pad(timepoint, 3, pad = "0")
      ))
  }

  calculate_single_r_squared <- function(gt, inference, threshold, class_string) {
    pred_at_thresh <- inference %>%
      filter(score >= threshold) %>%
      group_by(name, class) %>%
      summarize(model_count = n(), .groups = "drop") %>%
      complete(name, class) %>%
      mutate(model_count = replace_na(model_count, 0))

    df <- full_join(gt, pred_at_thresh, by = c("name", "class")) %>%
      filter(class == class_string) %>%
      complete(name, class) %>%
      mutate(model_count = replace_na(model_count, 0))

    r2 <- summary(lm(hand_count ~ model_count, data = df))$adj.r.squared
    data.frame(class = class_string, threshold = threshold, r_squared = r2)
  }

  train_cameras  <- c("combined", "one", "two")
  label_types    <- c("all", "pollen", "tube_tip")
  inference_cams <- c("one", "two")
  thresholds     <- seq(0.01, 1, by = 0.01)

  out <- list()
  for (tc in train_cameras) {
    for (lt in label_types) {
      for (ic in inference_cams) {
        # Individual-camera models are only evaluated on their own camera
        if ((tc == "one" & ic == "two") | (tc == "two" & ic == "one")) next

        for (tv in c("train", "val")) {
          gt_sub  <- process_gt(ground_truth, ic, tv)
          inf_sub <- process_inference(model_predictions, tc, lt, tv, ic)
          if (nrow(inf_sub) == 0) next

          for (cl in unique(inf_sub$class)) {
            r2_rows <- bind_rows(lapply(thresholds, function(thresh) {
              calculate_single_r_squared(gt_sub, inf_sub, thresh, cl)
            })) %>%
              filter(r_squared != 0)

            r2_rows$train_camera     <- tc
            r2_rows$label_type       <- lt
            r2_rows$train_or_val     <- tv
            r2_rows$inference_camera <- ic
            out[[length(out) + 1L]] <- r2_rows
          }
        }
      }
    }
  }
  bind_rows(out)
}

r_squared_df <- calculate_r_squared_df(ground_truth, model_predictions)

# Optimal threshold from training set
r_squared_summary_train <- r_squared_df %>%
  filter(train_or_val == "train", class != "unknown_germinated") %>%
  group_by(train_camera, inference_camera, label_type, class) %>%
  summarize(
    max_train_r_squared = max(r_squared),
    optimal_threshold   = threshold[which.max(r_squared)],
    .groups = "drop"
  )

# Validation r-squared at the threshold selected from the training set
r_squared_val <- r_squared_df %>%
  filter(train_or_val == "val", class != "unknown_germinated")

camera_comparison_data <- left_join(
  r_squared_summary_train,
  r_squared_val,
  by = c(
    "train_camera", "inference_camera", "label_type", "class",
    "optimal_threshold" = "threshold"
  )
) %>%
  rename(val_r_squared = r_squared)

# For the combined model (tc = "combined"), average r² across both inference
# cameras; individual models contribute a single inference camera each.
camera_comparison_avg <- camera_comparison_data %>%
  mutate(
    training_config = case_when(
      train_camera == "combined" ~ "Combined",
      train_camera == "one"      ~ "Camera 1 only",
      train_camera == "two"      ~ "Camera 2 only"
    ),
    label_type_label = case_when(
      label_type == "all"      ~ "All classes",
      label_type == "pollen"   ~ "Pollen classes",
      label_type == "tube_tip" ~ "Tube tip class"
    ),
    class_label = case_when(
      class == "ungerminated" ~ "Ungerminated",
      class == "germinated"   ~ "Germinated",
      class == "burst"        ~ "Burst",
      class == "aborted"      ~ "Aborted",
      class == "tube_tip"     ~ "Tube tip",
      TRUE                    ~ class
    )
  ) %>%
  group_by(training_config, label_type_label, class_label) %>%
  summarize(mean_val_r_squared = mean(val_r_squared, na.rm = TRUE), .groups = "drop")

class_order <- c("Ungerminated", "Germinated", "Burst", "Aborted", "Tube tip")
camera_comparison_avg <- camera_comparison_avg %>%
  filter(
    class_label      %in% class_order,
    label_type_label == "All classes"
  ) %>%
  mutate(
    class_label     = factor(class_label, levels = rev(class_order)),
    training_config = factor(
      training_config,
      levels = c("Camera 1 only", "Camera 2 only", "Combined")
    )
  )

config_colors <- c(
  "Camera 1 only" = "#56B4E9",
  "Camera 2 only" = "#E69F00",
  "Combined"      = "#009E73"
)
config_shapes <- c(
  "Camera 1 only" = 17,   # filled triangle
  "Camera 2 only" = 15,   # filled square
  "Combined"      = 16    # filled circle
)

camera_comparison_plot <- ggplot(
  camera_comparison_avg,
  aes(x = mean_val_r_squared, y = class_label,
      color = training_config, shape = training_config)
) +
  geom_vline(xintercept = 1, linewidth = 0.5, linetype = "dashed", color = "gray60") +
  geom_point(size = 3.5, position = position_dodge(width = 0.6)) +
  scale_x_continuous(
    breaks = seq(0, 1, 0.2),
    labels = seq(0, 1, 0.2),
    limits = c(0, 1.05),
    expand = c(0, 0)
  ) +
  scale_color_manual(values = config_colors, name = "Training set") +
  scale_shape_manual(values = config_shapes, name = "Training set") +
  labs(
    title = "Individual vs Combined Camera Training",
    x     = "Mean validation R\u00b2",
    y     = NULL
  ) +
  manuscript_theme(title_size = 24, axis_title_size = 18, axis_text_size = 16) +
  theme(
    legend.position    = "bottom",
    legend.title       = element_text(size = 18, face = "bold"),
    legend.text        = element_text(size = 14, face = "bold"),
    panel.grid.major.x = element_line(linewidth = 0.4, color = "gray88"),
    axis.line.y        = element_blank(),
    axis.ticks.y       = element_blank()
  ) +
  guides(
    color = guide_legend(nrow = 3, override.aes = list(size = 5)),
    shape = guide_legend(nrow = 3)
  )

ggsave(
  camera_comparison_plot,
  filename = file.path(plot_dir, "Individual vs Combined Camera Training.png"),
  device   = "png",
  width    = 8,
  height   = 6,
  dpi      = 400,
  units    = "in"
)


# Plot B: Model Training Metrics ------------------------------------------
# Plots total loss, mAP, average recall, and a combined optimization score
# vs training step for the model trained on data from both cameras.
# The red vertical line marks the step selected as the optimal checkpoint.

get_combined_metrics <- function(df) {
  weight_mAP        <- 1
  weight_total_loss <- 0.4

  df %>%
    filter(
      step != 0,
      step %% 1000 == 0,
      step >= 1000 & step <= 30000,
      train_or_eval == "eval"
    ) %>%
    group_by(camera, bounding_box_type, metric) %>%
    mutate(value_scaled = (value - min(value)) / (max(value) - min(value))) %>%
    ungroup() %>%
    group_by(camera, bounding_box_type, step) %>%
    summarize(
      combined_values = (
        value_scaled[metric == "mAP"] * weight_mAP -
        value_scaled[metric == "total_loss"] * weight_total_loss
      ),
      .groups = "drop"
    )
}

get_optimal_epoch <- function(combined_metrics_df) {
  combined_metrics_df %>%
    group_by(camera, bounding_box_type) %>%
    summarize(optimal_epoch = step[which.max(combined_values)], .groups = "drop")
}

combined_metrics_df <- get_combined_metrics(tensorboard_df)
optimal_epoch_df    <- get_optimal_epoch(combined_metrics_df)

make_training_plot <- function(
    camera_input, bounding_box_input,
    source_df, combined_metrics_df, optimal_epoch_df,
    add_x_title) {

  vline_value <- optimal_epoch_df$optimal_epoch[
    optimal_epoch_df$camera == camera_input &
    optimal_epoch_df$bounding_box_type == bounding_box_input
  ]
  label_pos <- ifelse(vline_value > 15000, vline_value - 2000, vline_value + 2000)
  label_hjust <- ifelse(vline_value > 15000, 1, 0)

  point_colors <- list(
    total_loss = "cyan", AR100 = "magenta", mAP = "orange", combined = "green"
  )
  x_title <- if (add_x_title) "Step (x1000)" else " "
  plot_title <- switch(bounding_box_input,
    "all"      = "All classes",
    "pollen"   = "Pollen classes",
    "tube_tip" = "Tube tip class"
  )

  base_theme <- manuscript_theme(
    title_size = 22, axis_title_size = 18, axis_text_size = 16
  ) +
    theme(
      legend.position  = "none",
      strip.background = element_blank(),
      strip.placement  = "outside",
      strip.text       = element_text(size = 24, face = "bold")
    )

  x_axis <- scale_x_continuous(
    breaks = seq(0, 30000, 10000),
    labels = seq(0, 30, 10),
    limits = c(0, 31000),
    expand = c(0, 0)
  )

  filter_metric <- function(df, cam, bb, met, tv) {
    df %>% filter(
      camera == cam, bounding_box_type == bb, metric == met,
      step >= 1000 & step <= 30000, train_or_eval == tv
    )
  }

  total_loss_plot <- filter_metric(source_df, camera_input, bounding_box_input, "total_loss", "eval") %>%
    ggplot(aes(x = step, y = value)) +
    geom_point(size = 3, color = point_colors$total_loss) +
    geom_smooth(linewidth = 1.5, se = FALSE, span = 0.2, color = "black") +
    geom_vline(xintercept = vline_value, linewidth = 1, color = "red") +
    labs(title = plot_title, y = "Total loss") +
    scale_y_continuous(
      breaks = seq(0, 1.5, 0.3), labels = seq(0, 1.5, 0.3),
      limits = c(0, 1.5), expand = c(0, 0)
    ) +
    x_axis + base_theme + theme(axis.title.x = element_blank())

  mAP_plot <- filter_metric(source_df, camera_input, bounding_box_input, "mAP", "eval") %>%
    ggplot(aes(x = step, y = value)) +
    geom_point(size = 3, color = point_colors$mAP) +
    geom_smooth(linewidth = 1.5, se = FALSE, span = 0.2, color = "black") +
    geom_vline(xintercept = vline_value, linewidth = 1, color = "red") +
    labs(y = "mAP") +
    scale_y_continuous(
      breaks = seq(0, 1, 0.2), labels = seq(0, 1, 0.2),
      limits = c(0, 1), expand = c(0, 0)
    ) +
    x_axis + base_theme + theme(axis.title.x = element_blank())

  AR100_plot <- filter_metric(source_df, camera_input, bounding_box_input, "AR100", "eval") %>%
    ggplot(aes(x = step, y = value)) +
    geom_point(size = 3, color = point_colors$AR100) +
    geom_smooth(linewidth = 1.5, se = FALSE, span = 0.2, color = "black") +
    geom_vline(xintercept = vline_value, linewidth = 1, color = "red") +
    labs(y = "Avg. recall") +
    scale_y_continuous(
      breaks = seq(0, 1, 0.2), labels = seq(0, 1, 0.2),
      limits = c(0, 1), expand = c(0, 0)
    ) +
    x_axis + base_theme + theme(axis.title.x = element_blank())

  combined_plot <- combined_metrics_df %>%
    filter(
      camera == camera_input, bounding_box_type == bounding_box_input,
      step >= 1000 & step <= 30000
    ) %>%
    ggplot(aes(x = step, y = combined_values)) +
    geom_point(size = 3, color = point_colors$combined) +
    geom_smooth(linewidth = 1.5, se = FALSE, span = 0.2, color = "black") +
    geom_vline(xintercept = vline_value, linewidth = 1, color = "red") +
    annotate(
      "text",
      x = label_pos, y = 1.1,
      label = vline_value, size = 5, color = "red", fontface = "bold",
      hjust = label_hjust
    ) +
    labs(x = x_title, y = "Combined") +
    scale_y_continuous(
      breaks = seq(0, 1.2, 0.2), labels = seq(0, 1.2, 0.2),
      limits = c(0, 1.2), expand = c(0, 0)
    ) +
    x_axis + base_theme

  total_loss_plot / mAP_plot / AR100_plot / combined_plot
}

all_plot    <- make_training_plot("combined", "all",      tensorboard_df, combined_metrics_df, optimal_epoch_df, FALSE)
pollen_plot <- make_training_plot("combined", "pollen",   tensorboard_df, combined_metrics_df, optimal_epoch_df, TRUE)
tube_plot   <- make_training_plot("combined", "tube_tip", tensorboard_df, combined_metrics_df, optimal_epoch_df, FALSE)

training_metrics_plot <- (all_plot | pollen_plot | tube_plot) +
  plot_annotation(
    title = "Model Training Metrics",
    theme = theme(
      plot.title = element_text(size = 28, face = "bold", hjust = 0.5,
                                margin = margin(10, 0, 0, 0))
    )
  )

ggsave(
  training_metrics_plot,
  filename = file.path(plot_dir, "Model Training Metrics.png"),
  device   = "png",
  width    = 10,
  height   = 12,
  dpi      = 400,
  units    = "in"
)


# Plot C: Model Confidence Threshold Selection ----------------------------
# Shows the r-squared between model predictions and hand counts across all
# confidence score thresholds (0.01–1.00) for the combined-camera all-classes
# model on the training set. R-squared values are averaged across both
# inference cameras (one and two).

threshold_data <- r_squared_df %>%
  filter(
    train_camera == "combined",
    label_type   == "all",
    train_or_val == "train",
    class        != "unknown_germinated"
  ) %>%
  group_by(class, threshold) %>%
  summarize(mean_r_squared = mean(r_squared, na.rm = TRUE), .groups = "drop") %>%
  filter(mean_r_squared > 0)

threshold_plot <- ggplot(threshold_data, aes(x = threshold, y = mean_r_squared, color = class)) +
  geom_point(size = 1.5) +
  scale_x_continuous(
    breaks = seq(0, 1, 0.2), labels = seq(0, 1, 0.2),
    limits = c(0, 1.05), expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.2), labels = seq(0, 1, 0.2),
    limits = c(0, 1.05), expand = c(0, 0)
  ) +
  scale_color_classes() +
  labs(
    title = "Model Confidence Threshold Selection",
    x     = "Confidence threshold",
    y     = "R\u00b2"
  ) +
  manuscript_theme(title_size = 24, axis_title_size = 18, axis_text_size = 16) +
  theme(
    legend.title     = element_blank(),
    legend.text      = element_text(size = 14, face = "bold"),
    legend.position  = "right",
    axis.title.x     = element_text(margin = margin(10, 0, 0, 0))
  )

ggsave(
  threshold_plot,
  filename = file.path(plot_dir, "Model Confidence Threshold Selection.png"),
  device   = "png",
  width    = 8,
  height   = 6,
  dpi      = 400,
  units    = "in"
)
