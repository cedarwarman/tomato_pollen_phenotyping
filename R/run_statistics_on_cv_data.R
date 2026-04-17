# Introduction ------------------------------------------------------------
# In this script ANCOVA is run on raw phenotypic data to adjust values based 
# on time-of-year.

library(dplyr)
library(car)
library(emmeans)
library(ggplot2)
library(googlesheets4)
library(stringr)
library(tidyr)

# Suppress Rplots.pdf: redirect the default graphics device to null so that
# diagnostic plot() calls in batch mode don't produce an unwanted file.
# ggsave() opens its own device and is unaffected by this.
pdf(NULL)

# For reproducible randomness
set.seed(16)

# Adding my Google service account credentials
gs4_auth(path = "~/.credentials/google_sheets_api/service_account.json")

# Some key accessions
# Heinz       - CW0000
# Tamaulipas  - CW0002
# 

# Loading the data --------------------------------------------------------
pollen_inference <- read.table(
    file = file.path(getwd(), "data", "processed_data", "interpolated_inference.tsv"),
    sep = '\t',
    header = TRUE
  ) %>%
  filter(date > "2021-11-09") %>% # Beginning of the experiment
  # Only including viable pollen for this analysis
  select(-percentage) %>%
  filter(object_class != "aborted") %>%
  ungroup() %>%
  group_by(name, time) %>%
  mutate(percentage = count / sum(count)) %>%
  ungroup() %>%
  # Adding the number of days since the start of the experiment
  mutate(days_since_start = as.numeric(round(difftime(date, min(date), units = "days"))))

tube_lengths <- read.table(
    file = file.path(getwd(), "data", "processed_data", "processed_tube_lengths.tsv"),
    sep = '\t',
    header = TRUE
  ) %>%
  filter(date > "2021-11-09") %>%
  mutate(name = paste0(
    date,
    "_run",
    run,
    "_",
    temp_target,
    "C_",
    well
  )) 


# Visualizing time-of-year effect -----------------------------------------
burst_accessions_by_date <- pollen_inference %>%
  filter(object_class == "burst") %>%
  filter((time == 80 & camera == 1) | (time == 65 & camera == 2))

make_burst_date_plot <- function(input_df, image_name) {
  color_vec <- c("#FF00FF", # Camera 1
                 "#1b74fa") # Camera 2
  names(color_vec) <- c("1", 
                        "2") 
  
  ggplot(input_df, aes(x = as.Date(date), 
                       y = percentage, 
                       color = factor(camera))) +
    geom_jitter(width = 0.5, alpha = 0.5) +
    scale_color_manual(values = color_vec,
                       name = "Camera",
                       breaks = c("1", 
                                  "2"), 
                       labels = c("1", 
                                  "2"), 
                       limits = force) +
    geom_smooth(aes(group = 1), method = "lm", color = "black") +
    scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
    scale_y_continuous(breaks = c(0, 0.25, .5, .75, 1),
                       labels = c("0%", "25%", "50%", "75%", "100%"),
                       limits = c(0, 1),
                       expand = c(0, 0)) +
    labs(title = image_name,
         x = "Date",
         y = "Burst % at 80 mins") +
    theme_bw() +
    theme(axis.title = element_text(size = 26, face = 'bold'),
          axis.text = element_text(size = 18, face = 'bold', color = 'black'),
          # axis.text.x = element_text(size = 18, face = 'bold', color = 'black'),
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
          legend.key.width = unit(1.5, "cm"),
          strip.background = element_blank(),
          strip.placement = "outside")
  
  ggsave(filename = file.path(getwd(), "plots", "ANCOVA", paste0(image_name, ".png")),
         device = 'png',
         width = 14,
         height = 8,
         dpi = 400,
         units = 'in')
}

pollen_inference %>%
  filter(
    object_class == "burst",
    time == 80,
    temp_target == 26
  ) %>%
  make_burst_date_plot("Bursting by date at 26 ºC")

pollen_inference %>%
  filter(
    object_class == "burst",
    time == 80,
    temp_target == 34 
  ) %>%
  make_burst_date_plot("Bursting by date at 34 ºC")


# Randomly selecting Heinz reps -------------------------------------------
# Heinz controls were run throughout the experiment. Here I randomly select 
# a subset of the reps for downstream analysis, stratified by temperature so 
# that exactly 8 reps are selected at each temperature.
heinz_reps <- pollen_inference %>%
  filter(accession == "CW0000") %>%
  distinct(temp_target, name) %>%
  group_by(temp_target) %>%
  slice_sample(n = 8) %>%
  pull(name)

# Subsetting the pollen class data with selected Heinz reps.
pollen_inference <- pollen_inference %>%
  filter(
    (accession != "CW0000") | 
    (accession == "CW0000" & name %in% heinz_reps)
  ) 

# Subsetting the tube length data with selected Heinz reps
tube_lengths_all_heinz <- tube_lengths
tube_lengths <- tube_lengths %>%
  filter(
    (accession != "CW0000") | 
    (accession == "CW0000" & name %in% heinz_reps)
  ) 


# ANCOVA: Burst at 2 hours ------------------------------------------------
# Making the models
ancova_burst_26C <- pollen_inference %>%
  filter(
    time == 80,
    object_class == "burst",
    temp_target == 26
  ) %>%
  lm(percentage ~ days_since_start + accession, data = .)
ancova_burst_34C <- pollen_inference %>%
  filter(
    time == 80,
    object_class == "burst",
    temp_target == 34
  ) %>%
  lm(percentage ~ days_since_start + accession, data = .)

# Testing linearity assumption (residuals should not vary systematically with 
# fitted values). There is a diagonal cutoff in this plot from the non-normal
# distribution of the percentages. 
plot(ancova_burst_26C, add.smooth = FALSE, which = 1)
plot(ancova_burst_34C, add.smooth = FALSE, which = 1)

# Testing the normality assumption.
plot(ancova_burst_26C, which = 2)
plot(ancova_burst_34C, which = 2)

# Testing the constant variance assumption.
plot(ancova_burst_26C, add.smooth = FALSE, which = 3)
plot(ancova_burst_34C, add.smooth = FALSE, which = 3)

# Running the ANCOVAs
anova(ancova_burst_26C)
anova(ancova_burst_34C)

# Function to get least squared means.
get_lsmeans <- function(prefix, ancova_output, raw_data) {
  lsmean_values <- as.data.frame(emmeans(ancova_output, "accession")) %>%
    select(accession, adjusted_mean = emmean, adjusted_se = SE) %>%
    rename_with(~ paste0(prefix, .), .cols = c("adjusted_mean", "adjusted_se"))
  
  raw_means <- raw_data %>%
    group_by(accession) %>%
    summarize(unadjusted_mean = mean(percentage),
              unadjusted_se = sd(percentage) / sqrt(n()),
              mean_day = mean(days_since_start)) %>%
    rename_with(~ paste0(prefix, .), .cols = c("unadjusted_mean", "unadjusted_se"))
  
  output_df <- left_join(raw_means, lsmean_values)
  
  return(output_df)
}

ancova_output_burst_2h_26C <- get_lsmeans(
  "burst_2h_26C_", 
  ancova_burst_26C, 
  pollen_inference %>%
    filter(
      time == 80,
      object_class == "burst",
      temp_target == 26 
    )
  )

ancova_output_burst_2h_34C <- get_lsmeans(
  "burst_2h_34C_", 
  ancova_burst_34C, 
  pollen_inference %>%
    filter(
      time == 80,
      object_class == "burst",
      temp_target == 34
    )
  )

# Combine LSMs for both temperatures into a single data frame
lsmeans_burst_2h <- ancova_output_burst_2h_26C %>%
  select(accession, burst_2h_26C_adjusted_mean) %>%
  left_join(
    ancova_output_burst_2h_34C %>% select(accession, burst_2h_34C_adjusted_mean),
    by = "accession"
  )

mean(lsmeans_burst_2h$burst_2h_26C_adjusted_mean)
mean(lsmeans_burst_2h$burst_2h_34C_adjusted_mean)

# Perform the Wilcoxon signed-rank test
wilcox.test(
  lsmeans_burst_2h$burst_2h_26C_adjusted_mean, lsmeans_burst_2h$burst_2h_34C_adjusted_mean, paired = TRUE, mu = 0
)

# Also comparing the variance with Levene's test
# Convert to long format
lsmeans_burst_2h_long <- lsmeans_burst_2h %>%
  pivot_longer(
    cols = starts_with("burst_2h"),
    names_to = "temp_target",
    values_to = "burst_2h_adjusted_mean"
  ) %>%
  mutate(
    temp_target = dplyr::recode(temp_target, 
                         "burst_2h_26C_adjusted_mean" = "26C", 
                         "burst_2h_34C_adjusted_mean" = "34C")
  )

leveneTest(
  burst_2h_adjusted_mean ~ as.factor(temp_target),
  data = lsmeans_burst_2h_long
)

# Visualizing least squared means adjustments by day-of-year
make_ancova_burst_plot <- function(input_df, image_name) {
  save_image_name <- str_replace_all(image_name, "%", "percent")
  
  plot_data <- input_df %>%
    select(-ends_with("_se")) %>%
    pivot_longer(cols = -c(accession, mean_day),
                 names_to = "name",
                 values_to = "mean_percentage") %>%
    mutate(
      mean_type = ifelse(grepl("unadjusted", name), "unadjusted", "adjusted")
    ) %>%
    select(accession, mean_type, mean_day, mean_percentage)
  
  ggplot(plot_data, 
    aes(
      x = mean_day, 
      y = mean_percentage, 
      group = mean_type, 
      shape = mean_type,
      fill = mean_type
    )
  ) +
  geom_line(aes(group = accession), linewidth = 0.6) +
  geom_point(size = 3, shape = 21, stroke = 0.5) +
  scale_fill_manual(values = c("pink", "white")) +
  scale_y_continuous(breaks = seq(0, 1, 0.2),
                     labels = seq(0, 100, 20),
                     limits = c(-0.05, 1.05)) +
  labs(title = image_name,
       y = "Mean percentage", 
       x = "Day",
       fill = "Mean type") +
  theme_bw() +
  theme(axis.title = element_text(size = 20, face = 'bold'),
        axis.text = element_text(size = 16, face = 'bold', color = 'black'),
        plot.title = element_text(size = 22, face = 'bold', margin = margin(0, 0, 10, 0)),
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
  
  ggsave(filename = file.path(getwd(), "plots", "ANCOVA", paste0(save_image_name, ".png")),
         device = 'png',
         width = 14,
         height = 8,
         dpi = 400,
         units = 'in')
}

make_ancova_burst_plot(ancova_output_burst_2h_26C, "Burst % at 2 hours, 26 ºC")
make_ancova_burst_plot(ancova_output_burst_2h_34C, "Burst % at 2 hours, 34 ºC")

# Cleaning up the data for export
export_burst_at_2_hours <- ancova_output_burst_2h_26C %>%
  rename(burst_2h_26C_mean_day = mean_day) %>%
  left_join(ancova_output_burst_2h_34C, by = c("accession")) %>%
  rename(burst_2h_34C_mean_day = mean_day) %>%
  na.omit() %>%
  mutate(
    unadjusted_bursting_increase = burst_2h_34C_unadjusted_mean - burst_2h_26C_unadjusted_mean,
    unadjusted_bursting_increase_se = sqrt(burst_2h_34C_unadjusted_se^2 + burst_2h_26C_unadjusted_se^2),
    adjusted_bursting_increase = burst_2h_34C_adjusted_mean - burst_2h_26C_adjusted_mean,
    adjusted_bursting_increase_se = sqrt(burst_2h_34C_adjusted_se^2 + burst_2h_26C_adjusted_se^2),
  )


# ANCOVA: Germination at 2 hours ------------------------------------------
# This is the percent of viable pollen that has germinated at 2 hours, e.g. it 
# has made some attempt to grow a tube, so this is the sum of percent germinated
# and percent burst, corrected by ANCOVA for day-of-year.
# Making the models
ancova_germ_26C <- pollen_inference %>%
  filter(
    time == 80,
    object_class == "germinated" |
      object_class == "burst",
    temp_target == 26
  ) %>%
  group_by(name, accession, days_since_start) %>%
  summarize(percentage = sum(percentage)) %>%
  lm(percentage ~ days_since_start + accession, data = .)
ancova_germ_34C <- pollen_inference %>%
  filter(
    time == 80,
    object_class == "germinated" |
      object_class == "burst",
    temp_target == 34
  ) %>%
  group_by(name, accession, days_since_start) %>%
  summarize(percentage = sum(percentage)) %>%
  lm(percentage ~ days_since_start + accession, data = .)

# Testing linearity assumption (residuals should not vary systematically with 
# fitted values). There is a diagonal cutoff in this plot from the non-normal
# distribution of the percentages. 
plot(ancova_germ_26C, add.smooth = FALSE, which = 1)
plot(ancova_germ_34C, add.smooth = FALSE, which = 1)

# Testing the normality assumption.
plot(ancova_germ_26C, which = 2)
plot(ancova_germ_34C, which = 2)

# Testing the constant variance assumption.
plot(ancova_germ_26C, add.smooth = FALSE, which = 3)
plot(ancova_germ_34C, add.smooth = FALSE, which = 3)

# Running the ANCOVAs
anova(ancova_germ_26C)
anova(ancova_germ_34C)

# Function to get least squared means.
get_lsmeans_germ <- function(prefix, ancova_output, raw_data) {
  lsmean_values <- as.data.frame(emmeans(ancova_output, "accession")) %>%
    select(accession, adjusted_mean = emmean, adjusted_se = SE) %>%
    rename_with(~ paste0(prefix, .), .cols = c("adjusted_mean", "adjusted_se"))
  
  raw_means <- raw_data %>%
    group_by(accession) %>%
    summarize(unadjusted_mean = mean(percentage),
              unadjusted_se = sd(percentage) / sqrt(n()),
              mean_day = mean(days_since_start)) %>%
    rename_with(~ paste0(prefix, .), .cols = c("unadjusted_mean", "unadjusted_se"))
  
  output_df <- left_join(raw_means, lsmean_values)
  
  return(output_df)
}

ancova_output_germ_2h_26C <- get_lsmeans_germ(
  "germ_2h_26C_", 
  ancova_germ_26C, 
  pollen_inference %>%
    filter(
      time == 80,
      object_class == "germinated" |
        object_class == "burst",
      temp_target == 26
    ) %>%
    group_by(name, accession, days_since_start) %>%
    summarize(percentage = sum(percentage))
)

ancova_output_germ_2h_34C <- get_lsmeans_germ(
  "germ_2h_34C_", 
  ancova_germ_34C, 
  pollen_inference %>%
    filter(
      time == 80,
      object_class == "germinated" |
        object_class == "burst",
      temp_target == 34
    ) %>%
    group_by(name, accession, days_since_start) %>%
    summarize(percentage = sum(percentage))
)

# Visualizing least squared means adjustments by day-of-year
make_ancova_germ_plot <- function(input_df, image_name) {
  save_image_name <- str_replace_all(image_name, "%", "percent")
  
  plot_data <- input_df %>%
    select(-ends_with("_se")) %>%
    pivot_longer(cols = -c(accession, mean_day),
                 names_to = "name",
                 values_to = "mean_percentage") %>%
    mutate(
      mean_type = ifelse(grepl("unadjusted", name), "unadjusted", "adjusted")
    ) %>%
    select(accession, mean_type, mean_day, mean_percentage)
  
  ggplot(plot_data, 
         aes(
           x = mean_day, 
           y = mean_percentage, 
           group = mean_type, 
           shape = mean_type,
           fill = mean_type
         )
  ) +
    geom_line(aes(group = accession), linewidth = 0.6) +
    geom_point(size = 3, shape = 21, stroke = 0.5) +
    scale_fill_manual(values = c("pink", "white")) +
    scale_y_continuous(breaks = seq(0, 1, 0.2),
                       labels = seq(0, 100, 20),
                       limits = c(-0.05, 1.05)) +
    labs(title = image_name,
         y = "Mean percentage", 
         x = "Day",
         fill = "Mean type") +
    theme_bw() +
    theme(axis.title = element_text(size = 20, face = 'bold'),
          axis.text = element_text(size = 16, face = 'bold', color = 'black'),
          plot.title = element_text(size = 22, face = 'bold', margin = margin(0, 0, 10, 0)),
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
  
  ggsave(filename = file.path(getwd(), "plots", "ANCOVA", paste0(save_image_name, ".png")),
         device = 'png',
         width = 14,
         height = 8,
         dpi = 400,
         units = 'in')
}

make_ancova_germ_plot(ancova_output_germ_2h_26C, "Germinated % at 2 hours, 26 ºC")
make_ancova_germ_plot(ancova_output_germ_2h_34C, "Germinated % at 2 hours, 34 ºC")

# Cleaning up the data for export
export_germ_at_2_hours <- ancova_output_germ_2h_26C %>%
  rename(germ_2h_26C_mean_day = mean_day) %>%
  left_join(ancova_output_germ_2h_34C, by = c("accession")) %>%
  rename(germ_2h_34C_mean_day = mean_day) %>%
  na.omit() %>%
  mutate(
    unadjusted_germ_decrease = germ_2h_26C_unadjusted_mean - germ_2h_34C_unadjusted_mean,
    unadjusted_germ_decrease_se = sqrt(germ_2h_26C_unadjusted_se^2 + germ_2h_34C_unadjusted_se^2),
    adjusted_germ_decrease = germ_2h_26C_adjusted_mean - germ_2h_34C_adjusted_mean,
    adjusted_germ_decrease_se = sqrt(germ_2h_26C_adjusted_se^2 + germ_2h_34C_adjusted_se^2)
  )


# ANCOVA: Time until 50% germinated ---------------------------------------
# This is the time until at least half the pollen has germinated (aka in the 
# germinated plus burst category), corrected by ANCOVA for day-of-year.
# Making the models
ancova_time_to_germ_26C <- pollen_inference %>%
  filter(temp_target == 26, object_class %in% c("germinated", "burst")) %>%
  group_by(name, accession, days_since_start, time) %>%
  summarize(total_percentage = sum(percentage), .groups = "drop") %>%
  group_by(name, accession, days_since_start) %>%
  summarize(
    time_to_50 = ifelse(any(total_percentage >= 0.50), 
                        min(time[total_percentage >= 0.50], na.rm = TRUE), 
                        NA),
    reaches_50 = any(total_percentage >= 0.50),
    .groups = "drop"
  ) %>%
  mutate(time_to_50 = ifelse(!reaches_50, 80, time_to_50)) %>%
  lm(time_to_50 ~ days_since_start + accession, data = .)

# Repeat for 34°C
ancova_time_to_germ_34C <- pollen_inference %>%
  filter(temp_target == 34, object_class %in% c("germinated", "burst")) %>%
  group_by(name, accession, days_since_start, time) %>%
  summarize(total_percentage = sum(percentage), .groups = "drop") %>%
  group_by(name, accession, days_since_start) %>%
  summarize(
    time_to_50 = ifelse(any(total_percentage >= 0.50), 
                        min(time[total_percentage >= 0.50], na.rm = TRUE), 
                        NA),
    reaches_50 = any(total_percentage >= 0.50),
    .groups = "drop"
  ) %>%
  mutate(time_to_50 = ifelse(!reaches_50, 80, time_to_50)) %>%
  lm(time_to_50 ~ days_since_start + accession, data = .)


# Testing linearity assumption (residuals should not vary systematically with 
# fitted values). There is a diagonal cutoff in this plot from the non-normal
# distribution of the percentages. 
plot(ancova_time_to_germ_26C, add.smooth = FALSE, which = 1)
plot(ancova_time_to_germ_34C, add.smooth = FALSE, which = 1)

# Testing the normality assumption.
plot(ancova_time_to_germ_26C, which = 2)
plot(ancova_time_to_germ_34C, which = 2)

# Testing the constant variance assumption.
plot(ancova_time_to_germ_26C, add.smooth = FALSE, which = 3)
plot(ancova_time_to_germ_34C, add.smooth = FALSE, which = 3)

# Running the ANCOVAs
anova(ancova_time_to_germ_26C)
anova(ancova_time_to_germ_34C)

# Function to get least squared means.
get_lsmeans_time_to_germ <- function(prefix, ancova_output, raw_data) {
  lsmean_values <- as.data.frame(emmeans(ancova_output, "accession")) %>%
    select(accession, adjusted_mean = emmean, adjusted_se = SE) %>%
    rename_with(~ paste0(prefix, .), .cols = c("adjusted_mean", "adjusted_se"))
  
  raw_means <- raw_data %>%
    group_by(accession) %>%
    summarize(
      unadjusted_mean = mean(time_to_50, na.rm = TRUE),
      unadjusted_se = sd(time_to_50, na.rm = TRUE) / sqrt(n()),
      mean_day = mean(days_since_start, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename_with(~ paste0(prefix, .), .cols = c("unadjusted_mean", "unadjusted_se"))
  
  output_df <- left_join(raw_means, lsmean_values, by = "accession")
  
  return(output_df)
}

ancova_output_time_to_germ_26C <- get_lsmeans_time_to_germ(
  "time_to_germ_26C_", 
  ancova_time_to_germ_26C, 
  pollen_inference %>%
    filter(temp_target == 26, object_class %in% c("germinated", "burst")) %>% # Include burst
    group_by(name, accession, days_since_start, time) %>% 
    summarize(total_percentage = sum(percentage), .groups = "drop") %>% # Sum percentages
    group_by(name, accession, days_since_start) %>%
    summarize(
      time_to_50 = ifelse(any(total_percentage >= 0.50), 
                          min(time[total_percentage >= 0.50], na.rm = TRUE), 
                          NA),
      reaches_50 = any(total_percentage >= 0.50),
      .groups = "drop"
    ) %>%
    mutate(time_to_50 = ifelse(!reaches_50, 80, time_to_50))
)

ancova_output_time_to_germ_34C <- get_lsmeans_time_to_germ(
  "time_to_germ_34C_", 
  ancova_time_to_germ_34C, 
  pollen_inference %>%
    filter(temp_target == 34, object_class %in% c("germinated", "burst")) %>% # Include burst
    group_by(name, accession, days_since_start, time) %>% 
    summarize(total_percentage = sum(percentage), .groups = "drop") %>% # Sum percentages
    group_by(name, accession, days_since_start) %>%
    summarize(
      time_to_50 = ifelse(any(total_percentage >= 0.50), 
                          min(time[total_percentage >= 0.50], na.rm = TRUE), 
                          NA),
      reaches_50 = any(total_percentage >= 0.50),
      .groups = "drop"
    ) %>%
    mutate(time_to_50 = ifelse(!reaches_50, 80, time_to_50))
)

# Combine LSMs for both temperatures into a single data frame
lsmeans_time_to_germ <- ancova_output_time_to_germ_26C %>%
  select(accession, time_to_germ_26C_adjusted_mean) %>%
  left_join(
    ancova_output_time_to_germ_34C %>% select(accession, time_to_germ_34C_adjusted_mean),
    by = "accession"
  )

mean(lsmeans_time_to_germ$time_to_germ_26C_adjusted_mean + 15)
mean(lsmeans_time_to_germ$time_to_germ_34C_adjusted_mean + 15)

# Perform the Wilcoxon signed-rank test
wilcox.test(
  lsmeans_time_to_germ$time_to_germ_26C_adjusted_mean, lsmeans_time_to_germ$time_to_germ_34C_adjusted_mean, paired = TRUE, mu = 0
)

# Also comparing the variance with Levene's test
# Convert to long format
lsmeans_time_to_germ_long <- lsmeans_time_to_germ %>%
  pivot_longer(
    cols = starts_with("time_to_germ"),
    names_to = "temp_target",
    values_to = "time_to_germ_adjusted_mean"
  ) %>%
  mutate(
    temp_target = dplyr::recode(temp_target, 
                         "time_to_germ_26C_adjusted_mean" = "26C", 
                         "time_to_germ_34C_adjusted_mean" = "34C")
  )

leveneTest(
  time_to_germ_adjusted_mean ~ as.factor(temp_target),
  data = lsmeans_time_to_germ_long
)

# Visualizing least squared means adjustments for time until 50% germination
make_ancova_time_to_germ_plot <- function(input_df, image_name) {
  save_image_name <- str_replace_all(image_name, "%", "percent")
  
  plot_data <- input_df %>%
    select(-ends_with("_se")) %>%
    pivot_longer(cols = -c(accession, mean_day),
                 names_to = "name",
                 values_to = "mean_time_to_50") %>%
    mutate(
      mean_type = ifelse(grepl("unadjusted", name), "unadjusted", "adjusted")
    ) %>%
    select(accession, mean_type, mean_day, mean_time_to_50)
  
  ggplot(plot_data, 
         aes(
           x = mean_day, 
           y = mean_time_to_50, 
           group = mean_type, 
           shape = mean_type,
           fill = mean_type
         )
  ) +
    geom_line(aes(group = accession), linewidth = 0.6) +
    geom_point(size = 3, shape = 21, stroke = 0.5) +
    scale_fill_manual(values = c("pink", "white")) +
    scale_y_continuous(breaks = seq(0, 80, 10),
                       limits = c(0, 85)) +
    labs(title = image_name,
         y = "Mean time until 50% germination (minutes)", 
         x = "Day",
         fill = "Mean type") +
    theme_bw() +
    theme(axis.title = element_text(size = 20, face = 'bold'),
          axis.text = element_text(size = 16, face = 'bold', color = 'black'),
          plot.title = element_text(size = 22, face = 'bold', margin = margin(0, 0, 10, 0)),
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
  
  ggsave(filename = file.path(getwd(), "plots", "ANCOVA", paste0(save_image_name, ".png")),
         device = 'png',
         width = 14,
         height = 8,
         dpi = 400,
         units = 'in')
}

# Generate plots for 26°C and 34°C
make_ancova_time_to_germ_plot(ancova_output_time_to_germ_26C, "Time to 50% Germination, 26 ºC")
make_ancova_time_to_germ_plot(ancova_output_time_to_germ_34C, "Time to 50% Germination, 34 ºC")

# Cleaning up the data for export
export_time_to_germ <- ancova_output_time_to_germ_26C %>%
  rename(time_to_germ_26C_mean_day = mean_day) %>%
  left_join(ancova_output_time_to_germ_34C, by = c("accession")) %>%
  rename(time_to_germ_34C_mean_day = mean_day) %>%
  na.omit() %>%
  mutate(
    unadjusted_time_increase = time_to_germ_34C_unadjusted_mean - time_to_germ_26C_unadjusted_mean,
    unadjusted_time_increase_se = sqrt(time_to_germ_26C_unadjusted_se^2 + time_to_germ_34C_unadjusted_se^2),
    adjusted_time_increase = time_to_germ_34C_adjusted_mean - time_to_germ_26C_adjusted_mean,
    adjusted_time_increase_se = sqrt(time_to_germ_26C_adjusted_se^2 + time_to_germ_34C_adjusted_se^2)
  )


# ANCOVA: Numerical integration of bursting percentage --------------------
# An approximation of the numerical integration of the bursting percentage 
# curve will better capture bursting over time. I will use the trapezoidal 
# rule to make the approximation.
burst_integral_data <- pollen_inference %>%
  filter(object_class == "burst") %>%
  select(name, time, date, percentage, camera, accession, temp_target, days_since_start) %>%
  group_by(name, date, camera, accession, temp_target, days_since_start) %>%
  arrange(time) %>%
  mutate(percentage = if_else(is.nan(percentage), 0, percentage)) %>% # Fixes a small number of NaNs that should be 0
  summarize(
    burst_integral = sum((percentage + lag(percentage, default = first(percentage))) / 2 * 
                           (time - lag(time, default = first(time)))),
    .groups = "drop"
  )

ancova_burst_integral_26C <- burst_integral_data %>%
  filter(temp_target == 26) %>%
  lm(burst_integral ~ days_since_start + accession, data = .)
ancova_burst_integral_34C <- burst_integral_data %>%
  filter(temp_target == 34) %>%
  lm(burst_integral ~ days_since_start + accession, data = .)

# Visualizing raw data
ggplot(ancova_burst_integral_26C, aes(x = days_since_start, y = burst_integral, colour = accession)) + 
  geom_point() +
  theme_bw() +
  theme(legend.position = "none")
ggplot(ancova_burst_integral_34C, aes(x = days_since_start, y = burst_integral, colour = accession)) + 
  geom_point() +
  theme_bw() +
  theme(legend.position = "none")

# Testing linearity assumption. 
plot(ancova_burst_integral_26C, add.smooth = FALSE, which = 1)
plot(ancova_burst_integral_34C, add.smooth = FALSE, which = 1)

# Testing the normality assumption.
plot(ancova_burst_integral_26C, which = 2)
plot(ancova_burst_integral_34C, which = 2)

# Testing the constant variance assumption.
plot(ancova_burst_integral_26C, add.smooth = FALSE, which = 3)
plot(ancova_burst_integral_34C, add.smooth = FALSE, which = 3)

# Getting the p-values.
anova(ancova_burst_integral_26C)
anova(ancova_burst_integral_34C)

# Getting the least squared means (and unadjusted means)
get_lsmeans_burst_integral <- function(prefix, ancova_output, raw_data) {
  lsmean_values <- as.data.frame(emmeans(ancova_output, "accession")) %>%
    select(accession, adjusted_mean = emmean, adjusted_se = SE) %>%
    rename_with(~ paste0(prefix, .), .cols = c("adjusted_mean", "adjusted_se"))
  
  raw_means <- raw_data %>%
    group_by(accession) %>%
    summarize(unadjusted_mean = mean(burst_integral),
              unadjusted_se = sd(burst_integral) / sqrt(n()),
              mean_day = mean(days_since_start)) %>%
    rename_with(~ paste0(prefix, .), .cols = c("unadjusted_mean", "unadjusted_se", "mean_day"))
  
  output_df <- left_join(raw_means, lsmean_values)
  
  return(output_df)
}

ancova_output_burst_integral_26C <- get_lsmeans_burst_integral(
  "burst_integral_26C_", 
  ancova_burst_integral_26C, 
  burst_integral_data[burst_integral_data$temp_target == 26, ])
ancova_output_burst_integral_34C <- get_lsmeans_burst_integral(
  "burst_integral_34C_", 
  ancova_burst_integral_34C, 
  burst_integral_data[burst_integral_data$temp_target == 34, ])

# Combine LSMs for both temperatures into a single data frame
lsmeans_burst_integral <- ancova_output_burst_integral_26C %>%
  select(accession, burst_integral_26C_adjusted_mean) %>%
  left_join(
    ancova_output_burst_integral_34C %>% select(accession, burst_integral_34C_adjusted_mean),
    by = "accession"
  )

mean(lsmeans_burst_integral$burst_integral_26C_adjusted_mean)
mean(lsmeans_burst_integral$burst_integral_34C_adjusted_mean)

# Perform the Wilcoxon signed-rank test
wilcox.test(
  lsmeans_burst_integral$burst_integral_26C_adjusted_mean, lsmeans_burst_integral$burst_integral_34C_adjusted_mean, paired = TRUE, mu = 0
)

# Also comparing the variance with Levene's test
# Convert to long format
lsmeans_burst_integral_long <- lsmeans_burst_integral %>%
  pivot_longer(
    cols = starts_with("burst_integral"),
    names_to = "temp_target",
    values_to = "burst_integral_adjusted_mean"
  ) %>%
  mutate(
    temp_target = dplyr::recode(temp_target, 
                         "burst_integral_26C_adjusted_mean" = "26C", 
                         "burst_integral_34C_adjusted_mean" = "34C")
  )

leveneTest(
  burst_integral_adjusted_mean ~ as.factor(temp_target),
  data = lsmeans_burst_integral_long
)

make_ancova_burst_integral_plot <- function(input_df, image_name) {
  names(input_df)[4] <- "mean_day"
  
  plot_data <- input_df %>%
    select(-ends_with("_se")) %>%
    pivot_longer(cols = -c(accession, mean_day),
                 names_to = "name",
                 values_to = "mean_burst_integral") %>%
    mutate(
      mean_type = ifelse(grepl("unadjusted", name), "unadjusted", "adjusted")
    ) %>%
    select(accession, mean_type, mean_day, mean_burst_integral)
  
  ggplot(plot_data, 
         aes(
           x = mean_day, 
           y = mean_burst_integral, 
           group = mean_type, 
           shape = mean_type,
           fill = mean_type
         )
  ) +
    geom_line(aes(group = accession), linewidth = 0.6) +
    geom_point(size = 3, shape = 21, stroke = 0.5) +
    scale_fill_manual(values = c("pink", "white")) +
    scale_y_continuous(breaks = seq(0, 60, 10),
                       labels = seq(0, 60, 10),
                       limits = c(0, 62)) +
    labs(title = image_name,
         y = "Mean burst integral", 
         x = "Day",
         fill = "Mean type") +
    theme_bw() +
    theme(axis.title = element_text(size = 20, face = 'bold'),
          axis.text = element_text(size = 16, face = 'bold', color = 'black'),
          plot.title = element_text(size = 22, face = 'bold', margin = margin(0, 0, 10, 0)),
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
  
  ggsave(filename = file.path(getwd(), "plots", "ANCOVA", paste0(image_name, ".png")),
         device = 'png',
         width = 14,
         height = 8,
         dpi = 400,
         units = 'in')
}

make_ancova_burst_integral_plot(ancova_output_burst_integral_26C, "Burst integral at 26 ºC")
make_ancova_burst_integral_plot(ancova_output_burst_integral_34C, "Burst integral at 34 ºC")

# Formatting data for export, calculating absolute burst integral increase.
export_absolute_burst_integral <- ancova_output_burst_integral_26C %>%
  left_join(ancova_output_burst_integral_34C) %>%
  na.omit() %>%
  mutate(
    unadjusted_integral_increase = burst_integral_34C_unadjusted_mean - burst_integral_26C_unadjusted_mean,
    unadjusted_integral_increase_se = sqrt(burst_integral_34C_unadjusted_se^2 + burst_integral_26C_unadjusted_se^2),
    adjusted_integral_increase = burst_integral_34C_adjusted_mean - burst_integral_26C_adjusted_mean,
    adjusted_integral_increase_se = sqrt(burst_integral_34C_adjusted_se^2 + burst_integral_26C_adjusted_se^2)
  )

# Comparing absolute burst integral increase for Heinz and Tamaulipas to the mean 
mean(export_absolute_burst_integral$adjusted_integral_increase)
export_absolute_burst_integral %>%
  filter(accession == "CW0000") %>%
  select(adjusted_integral_increase)
export_absolute_burst_integral %>%
  filter(accession == "CW0002") %>%
  select(adjusted_integral_increase)

# Plotting adjusted burst integral by accession (sorted by adjusted mean) ---
# Uses ANCOVA lsmeans (adjusted mean ± SE) rather than raw per-rep values, so
# accession means are comparable after removing the day-of-year confound.
make_burst_integral_accession_plot <- function(input_df, image_name) {
  names(input_df)[grep("_adjusted_mean$", names(input_df))] <- "adjusted_mean"
  names(input_df)[grep("_adjusted_se$",   names(input_df))] <- "adjusted_se"

  ggplot(input_df, aes(x = reorder(accession, adjusted_mean), y = adjusted_mean)) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5, linetype = "dashed") +
    geom_errorbar(aes(ymin = adjusted_mean - adjusted_se,
                      ymax = adjusted_mean + adjusted_se),
                  width = 0, color = "black", linewidth = 0.5) +
    geom_point(size = 1.5, color = "black") +
    scale_y_continuous(breaks = seq(-10, 70, 10),
                       limits = c(-10, 70)) +
    labs(title = image_name,
         y = "Burst Integral") +
    theme_bw() +
    theme(axis.title        = element_text(size = 26, face = "bold"),
          axis.title.y      = element_text(size = 26, face = "bold", margin = margin(r = 20)),
          axis.text         = element_text(size = 22, face = "bold", color = "black"),
          axis.title.x      = element_blank(),
          axis.text.x       = element_blank(),
          axis.ticks.x      = element_blank(),
          plot.title        = element_text(size = 28, face = "bold", margin = margin(0, 0, 10, 0)),
          panel.border      = element_blank(),
          axis.line         = element_line(linewidth = 1, color = "black"),
          axis.ticks        = element_line(linewidth = 1, color = "black"),
          axis.ticks.length = unit(8, "pt"),
          plot.margin       = margin(0.5, 0.5, 0.5, 0.5, "cm"),
          panel.grid        = element_blank(),
          legend.position   = "none")

  ggsave(filename = file.path(getwd(), "plots", "ANCOVA", paste0(image_name, ".png")),
         device = "png",
         width = 20,
         height = 7,
         dpi = 400,
         units = "in")
}

# Restrict plots to accessions with >= 7 reps at both temperatures.
# Accessions with fewer reps have disproportionately large SE and represent
# incomplete data relative to the intended experimental design.
min_reps <- 7
accessions_min_reps <- burst_integral_data %>%
  group_by(accession, temp_target) %>%
  summarize(n_reps = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = temp_target, values_from = n_reps,
                     names_prefix = "reps_") %>%
  filter(!is.na(reps_26), !is.na(reps_34),
         reps_26 >= min_reps, reps_34 >= min_reps) %>%
  pull(accession)

make_burst_integral_accession_plot(
  ancova_output_burst_integral_26C %>% filter(accession %in% accessions_min_reps),
  "Adjusted burst integral by accession at 26 \u00baC")
make_burst_integral_accession_plot(
  ancova_output_burst_integral_34C %>% filter(accession %in% accessions_min_reps),
  "Adjusted burst integral by accession at 34 \u00baC")

# Plotting adjusted burst integral increase by accession --------------------
make_adjusted_integral_increase_plot <- function(input_df, image_name) {
  ggplot(input_df, aes(x = reorder(accession, adjusted_integral_increase),
                       y = adjusted_integral_increase)) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5, linetype = "dashed") +
    geom_errorbar(aes(ymin = adjusted_integral_increase - adjusted_integral_increase_se,
                      ymax = adjusted_integral_increase + adjusted_integral_increase_se),
                  width = 0, color = "black", linewidth = 0.5) +
    geom_point(size = 1.5, color = "black") +
    scale_y_continuous(breaks = seq(-10, 70, 10),
                       limits = c(-10, 70)) +
    labs(title = image_name,
         y = "Burst Integral Increase") +
    theme_bw() +
    theme(axis.title        = element_text(size = 26, face = "bold"),
          axis.title.y      = element_text(size = 26, face = "bold", margin = margin(r = 20)),
          axis.text         = element_text(size = 22, face = "bold", color = "black"),
          axis.title.x      = element_blank(),
          axis.text.x       = element_blank(),
          axis.ticks.x      = element_blank(),
          plot.title        = element_text(size = 28, face = "bold", margin = margin(0, 0, 10, 0)),
          panel.border      = element_blank(),
          axis.line         = element_line(linewidth = 1, color = "black"),
          axis.ticks        = element_line(linewidth = 1, color = "black"),
          axis.ticks.length = unit(8, "pt"),
          plot.margin       = margin(0.5, 0.5, 0.5, 0.5, "cm"),
          panel.grid        = element_blank(),
          legend.position   = "none")

  ggsave(filename = file.path(getwd(), "plots", "ANCOVA", paste0(image_name, ".png")),
         device = "png",
         width = 20,
         height = 7,
         dpi = 400,
         units = "in")
}

make_adjusted_integral_increase_plot(
  export_absolute_burst_integral %>% filter(accession %in% accessions_min_reps),
  "Adjusted burst integral increase by accession")


# Tube length -------------------------------------------------------------
# Getting the mean tube length for each rep at 26 and 34.
tube_length_data <- tube_lengths %>%
  group_by(camera, temp_target, accession, name) %>%
  summarize(mean_tube_length = sum(total_length_um) / n()) %>%
  ungroup() %>%
  # Adding the days since the start of the experiment
  mutate(date = as.Date(substr(name, 1, 10)),
         days_since_start = as.numeric(round(difftime(date, min(date), units = "days"))))

tube_length_data_all_heinz <- tube_lengths_all_heinz %>%
  group_by(camera, temp_target, accession, name) %>%
  summarize(mean_tube_length = sum(total_length_um) / n()) %>%
  ungroup() %>%
  # Adding the days since the start of the experiment
  mutate(date = as.Date(substr(name, 1, 10)),
         days_since_start = as.numeric(round(difftime(date, min(date), units = "days"))))

# Camera 2 has a wider field of view than camera 1, so tubes appear longer. 
# However, here we will look at a ratio between the two cameras, which normalizes 
# for the two fields of view. I will also calculate the standard error for the 
# tube length means and propagate the error for the ratios.
tube_length_ratio_data <- tube_length_data %>%
  ungroup() %>%
  group_by(temp_target, accession) %>%
  summarize(
    mean_tube_length_sd = sd(mean_tube_length), 
    mean_tube_length = mean(mean_tube_length), 
    n = n(),
    days_since_start = mean(days_since_start)
  ) %>%
  group_by(accession) %>%
  summarize(
    mean_tube_length_34 = mean(mean_tube_length[temp_target == 34]),
    mean_tube_length_34_se = mean_tube_length_sd[temp_target == 34] / sqrt(n[temp_target == 34]),
    mean_tube_length_26 = mean(mean_tube_length[temp_target == 26]),
    mean_tube_length_26_se = mean_tube_length_sd[temp_target == 26] / sqrt(n[temp_target == 26]),
    days_since_start = mean(days_since_start)
  ) %>%
  mutate(
    ratio = mean_tube_length_34 / mean_tube_length_26,
    ratio_se = ratio * sqrt((mean_tube_length_34_se / mean_tube_length_34)^2 + (mean_tube_length_26_se / mean_tube_length_26)^2)
  ) %>%
  select(accession, ratio, ratio_se, days_since_start)

# In the case of tube length ratio, it doesn't look like there's a trend from day-of-year.
ggplot(tube_length_ratio_data, aes(x = days_since_start, y = ratio)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw()

# Fitting a linear regression to the ratio. day_since_start not significant for 
# tube length ratio (p = 0.805), therefore ANCOVA is not necessary for the tube 
# length ratios and the raw data can be used.
lm_tube_length_ratio <- lm(ratio ~ days_since_start, data = tube_length_ratio_data)
summary(lm_tube_length_ratio)

# Looking at tube lengths by date and camera at 26 and 34 C
make_tube_length_date_plot <- function(input_df, image_name, column_str) {
  color_vec <- c("#FF00FF", # Camera 1
                 "#1b74fa") # Camera 2
  names(color_vec) <- c("1", 
                        "2") 
  
  ggplot(input_df, aes(x = as.Date(date), 
                       y = {{column_str}}, 
                       color = factor(camera))) +
    geom_jitter(size = 3, width = 0.5, alpha = 0.5) +
    scale_color_manual(values = color_vec,
                       name = "Camera",
                       breaks = c("1", 
                                  "2"), 
                       labels = c("1", 
                                  "2"), 
                       limits = force) +
    # geom_smooth(aes(group = 1), method = "lm", color = "black") +
    scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
    scale_y_continuous(breaks = seq(600, 1200, 100),
                       labels = seq(600, 1200, 100),
                       limits = c(500, 1300)) +
    labs(title = image_name,
         x = "Date",
         y = "Tube length at 2 hours (µm)") +
    theme_bw() +
    theme(axis.title = element_text(size = 26, face = 'bold'),
          axis.text = element_text(size = 18, face = 'bold', color = 'black'),
          # axis.text.x = element_text(size = 18, face = 'bold', color = 'black'),
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
          legend.key.width = unit(1.5, "cm"),
          strip.background = element_blank(),
          strip.placement = "outside")
  
  ggsave(filename = file.path(getwd(), "plots", "tube_lengths", paste0(image_name, ".png")),
         device = 'png',
         width = 14,
         height = 8,
         dpi = 400,
         units = 'in')
}

tube_length_data_all_heinz %>%
  filter(
    accession == "CW0000",
    temp_target == 26
  ) %>%
  make_tube_length_date_plot("Heinz tube length by date at 26 ºC", mean_tube_length)

tube_length_data_all_heinz %>%
  filter(
    accession == "CW0000",
    temp_target == 34 
  ) %>%
  make_tube_length_date_plot("Heinz tube length by date at 34 ºC", mean_tube_length)

tube_length_data_all_heinz %>%
  filter(
    temp_target == 26
  ) %>%
  make_tube_length_date_plot("All tube length by date at 26 ºC", mean_tube_length)

tube_length_data_all_heinz %>%
  filter(
    temp_target == 34 
  ) %>%
  make_tube_length_date_plot("All tube length by date at 34 ºC", mean_tube_length)

# To directly compare tube lengths between cameras, they must be adjusted based 
# on a conversion factor calculated from the mean Heinz length in each camera. I 
# will combine the two temperatures for each camera to get the best estimate of 
# the inherent tube length difference between cameras.
heinz_mean_tube_length_camera_1 <- tube_length_data_all_heinz %>% 
  filter(accession == "CW0000", camera == 1) %>%
  summarize(mean = mean(mean_tube_length)) %>%
  pull()
heinz_mean_tube_length_camera_2 <- tube_length_data_all_heinz %>% 
  filter(accession == "CW0000", camera == 2) %>%
  summarize(mean = mean(mean_tube_length)) %>%
  pull()

tube_length_camera_adjustment_factor <- heinz_mean_tube_length_camera_2 / heinz_mean_tube_length_camera_1

# Doing the adjustment
adjusted_tube_length_data <- tube_length_data %>%
  mutate(
    adjusted_mean_tube_length = ifelse(
      camera == 2, mean_tube_length / tube_length_camera_adjustment_factor, 
      mean_tube_length
    )
  )

# Plotting the results
adjusted_tube_length_data %>%
  filter(
    temp_target == 26
  ) %>%
  make_tube_length_date_plot("Adjusted tube lengths by date at 26 ºC", adjusted_mean_tube_length)

adjusted_tube_length_data %>%
  filter(
    temp_target == 34 
  ) %>%
  make_tube_length_date_plot("Adjusted tube lengths by date at 34 ºC", adjusted_mean_tube_length)

# Getting the mean adjusted tube lengths by accession and temperature, 
# calculating standard error.
mean_adjusted_tube_lengths <- adjusted_tube_length_data %>%
  group_by(temp_target, accession) %>%
  summarize(
    mean_adjusted_tube_length_se = sd(adjusted_mean_tube_length) / sqrt(n()), 
    mean_adjusted_tube_length = mean(adjusted_mean_tube_length), 
    n = n(),
    days_since_start = mean(days_since_start)
  )

# Looking at if the tube lengths are different at heat stress:
mean(mean_adjusted_tube_lengths$mean_adjusted_tube_length[mean_adjusted_tube_lengths$temp_target == "26"])
mean(mean_adjusted_tube_lengths$mean_adjusted_tube_length[mean_adjusted_tube_lengths$temp_target == "34"])

wilcox.test(
  mean_adjusted_tube_length ~ temp_target, 
  data = mean_adjusted_tube_lengths, 
  paired = TRUE,
  mu = 0
)

# Counting accessions by their tube length ratios
percentage_ratios <- tube_length_ratio_data %>%
  mutate(category = ifelse(ratio < 1, "<1", ">=1")) %>%
  group_by(category) %>%
  summarize(count = n()) %>%
  mutate(percentage = (count / sum(count)) * 100)  
percentage_ratios

tube_length_ratio_data[tube_length_ratio_data$accession == "CW0000", ]$ratio
tube_length_ratio_data[tube_length_ratio_data$accession == "CW0002", ]$ratio

# # SE quality control, looking at CW1004 to make sure everything is working properly.
# tube_length_data %>%
#   filter(
#     temp_target == 34
#   ) %>%
#   ggplot(aes(x = mean_tube_length)) +
#   geom_histogram(fill = "magenta", binwidth = 20) +
#   theme_minimal()
# tube_lengths %>%
#   filter(
#     temp_target == 34
#   ) %>%
#   ggplot(aes(x = total_length_um)) +
#   geom_histogram(fill = "magenta", binwidth = 20) +
#   theme_minimal()
# 
# tube_length_data %>%
#   filter(
#     accession == "CW0163",
#     temp_target == 34
#   ) %>%
#   ggplot(aes(x = mean_tube_length)) +
#   geom_histogram(fill = "magenta", binwidth = 20) +
#   theme_minimal()
# tube_length_data %>%
#   filter(
#     accession == "CW0163",
#     temp_target == 26 
#   ) %>%
#   ggplot(aes(x = mean_tube_length)) +
#   geom_histogram(fill = "magenta", binwidth = 20) +
#   theme_minimal()
# tube_length_data %>%
#   filter(
#     accession == "CW1004",
#     temp_target == 34
#   ) %>%
#   select(mean_tube_length) %>%
#   summarize(se = sd(mean_tube_length) / sqrt(n()))
# tube_length_data %>%
#   filter(
#     accession == "CW1004",
#     temp_target == 26
#   ) %>%
#   select(mean_tube_length) %>%
#   summarize(se = sd(mean_tube_length) / sqrt(n()))
# 
# tube_length_data %>%
#   filter(
#     accession == "CW0069",
#     temp_target == 34
#   ) %>%
#   ggplot(aes(x = mean_tube_length)) +
#   geom_histogram(fill = "magenta", binwidth = 20) +
#   theme_minimal()
# tube_length_data %>%
#   filter(
#     accession == "CW0069",
#     temp_target == 26 
#   ) %>%
#   ggplot(aes(x = mean_tube_length)) +
#   geom_histogram(fill = "magenta", binwidth = 20) +
#   theme_minimal()
# 
# tube_length_data %>%
#   filter(
#     accession == "CW0058",
#     temp_target == 34
#   ) %>%
#   ggplot(aes(x = mean_tube_length)) +
#   geom_histogram(fill = "magenta", binwidth = 20) +
#   theme_minimal()
# tube_length_data %>%
#   filter(
#     accession == "CW0058",
#     temp_target == 26 
#   ) %>%
#   ggplot(aes(x = mean_tube_length)) +
#   geom_histogram(fill = "magenta", binwidth = 20) +
#   theme_minimal()
# 
# test_df_1 <- tube_lengths %>%
#   filter(name == "2022-03-04_run1_34C_D1")
# tube_lengths %>%
#   filter(name == "2022-03-04_run1_34C_D1") %>%
#   ggplot(aes(x = total_length_um)) +
#   geom_histogram(fill = "magenta", binwidth = 100) +
#   theme_minimal()
# test_df_2 <- tube_lengths %>%
#   filter(name == "2022-03-04_run1_34C_C6")
# tube_lengths %>%
#   filter(name == "2022-03-04_run1_34C_C6") %>%
#   ggplot(aes(x = total_length_um)) +
#   geom_histogram(fill = "magenta", binwidth = 100) +
#   theme_minimal()

# Exporting data ----------------------------------------------------------
# Data output for downstream analysis and visualization.

# Bursting at 2 hours.
# export_burst_at_2_hours

# Adjusted burst integral (at 26 and 34, including absolute burst integral 
# increase, adjusted by day-of-year)
# export_absolute_burst_integral

# Tube length at 2 hours
# For the ratio and the propagated ratio se (mean lengths and se here are unadjusted)
export_tube_length_ratio_data <- tube_length_ratio_data %>%
  select(-days_since_start) %>%
  rename(
    unadjusted_tube_length_ratio = ratio, 
    unadjusted_tube_length_ratio_se = ratio_se
  )

# Formatting the adjusted tube lengths
mean_adjusted_tube_lengths_wide <- mean_adjusted_tube_lengths %>%
  select(-n, -days_since_start) %>%
  pivot_wider(
    id_cols = c(accession),
    names_from = temp_target,
    values_from = c(mean_adjusted_tube_length, mean_adjusted_tube_length_se),
    names_sep = "_"
  ) %>%
  rename(
    tube_length_26C_adjusted_mean = mean_adjusted_tube_length_26,
    tube_length_26C_adjusted_se = mean_adjusted_tube_length_se_26,
    tube_length_34C_adjusted_mean = mean_adjusted_tube_length_34,
    tube_length_34C_adjusted_se = mean_adjusted_tube_length_se_34,
  )

# Combining the data frames
export_df <- export_burst_at_2_hours %>%
  left_join(export_germ_at_2_hours, by = c("accession")) %>%
  left_join(export_absolute_burst_integral, by = c("accession")) %>%
  left_join(mean_adjusted_tube_lengths_wide, by = c("accession")) %>%
  left_join(export_tube_length_ratio_data, by = c("accession"))

# Adding metadata
export_metadata <- read_sheet("1V2kH8G4tfYsYqnYb6bVHpGjwwkjd9le0arGBJ2o4r8s") %>%
  select(accession = name_CW, packet_name_1, name_TGRC, species)

export_df <- left_join(export_df, export_metadata, by = "accession") %>%
  select(accession, 
         name_TGRC, 
         packet_name_1, 
         species, 
         everything())

# Uploading to Google Sheet
# write_sheet(export_df, 
#             ss = "1PXbBhTxpJh_FlhnRDwcegovdqPq8_sYgDFgge0V5Vg0",
#             sheet = "Sheet1")

# Saving for downstream analysis.
export_df %>%
  write.table(
    file = file.path(getwd(), "data", "processed_data", "pollen_measurements.tsv"),
    sep = '\t',
    quote = FALSE,
    row.names = FALSE
  )
