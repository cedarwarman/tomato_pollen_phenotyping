# Introduction ------------------------------------------------------------
# This script prepares data from the computer vision pipeline for downstream 
# analysis. Input data comes from this Python script:
# https://github.com/cedarwarman/pollen_cv/blob/main/python/process_inference.py

library(dplyr)
library(tidyr)
library(googlesheets4)
library(purrr)
library(stringr)

# Adding my Google service account credentials
gs4_auth(path = "~/.credentials/google_sheets_api/service_account.json")


# Importing the computer vision inference ---------------------------------
inference <- read.table(
  file = file.path(getwd(),
                   "data",
                   "model_predictions",
                   "2023-11-01_all_tracks.txt"),
  sep = '\t',
  header = TRUE
)


# Importing the metadata --------------------------------------------------
wells_to_accessions <- read_sheet("1yQ5yAKiL6BzwZ-wH-Q44RoUEwMZztTYafzdvVylq6fo")

# Only keeping the columns we need
wells_to_accessions <- wells_to_accessions %>%
  select(date, run, well, temp_target, accession) %>%
  mutate(date = as.character(date))

# Fixing typos
wells_to_accessions <- wells_to_accessions %>%
  mutate(accession = replace(accession, accession == "CW00173", "CW0173")) %>%
  mutate(accession = replace(accession, accession == "CW01000", "CW1000")) %>%
  mutate(accession = replace(accession, accession == "CW019", "CW0019")) %>%
  mutate(accession = replace(accession, accession == "CW01037", "CW1037")) %>%
  mutate(accession = replace(accession, accession == "CW089", "CW0089")) 


# Cleaning the data for pollen classes ------------------------------------
# Here we will do some cleaning and summarizing of pollen classes. This function
# makes total counts for each class at each time point for each rep.
process_data <- function(df){
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
      str_pad(t, 3, pad = "0")
    )) %>%
    filter(object_class != "tube_tip",
           object_class != "object_class") %>%
    group_by(name, object_class) %>%
    summarize(count = n()) %>%
    mutate(percentage = count / sum(count)) %>%
    ungroup() %>%
    complete(name, object_class) %>%
    mutate(count = replace_na(count, 0)) %>%
    mutate(percentage = replace_na(percentage, 0)) %>%
    mutate(time = as.integer(str_sub(name, -3, -1)),
           date = str_sub(name, 1, 10),
           run = as.double(str_sub(name, 15, 15)),
           well = str_sub(name, 21, 22)) %>%
    filter(time <= 82)

  # Adding a factor for camera
  df$camera <- NA
  df$camera[as.Date(df$date) < as.Date("2022-05-27")] <- 1
  df$camera[as.Date(df$date) >= as.Date("2022-05-27")] <- 2

  return(df)
}

processed_pollen_inference <- process_data(inference)


# Rescaling camera 2 to camera 1 range with linear interpolation ----------
# Because of different capturing speeds, the two cameras have data on different
# scales that both cover 2 hours. Camera one goes from frame 0 to frame 82 while
# camera 2 goes from frame 0 to frame 65. Camera 2 images are rescaled to the
# range of camera 1 images using a linear interpolation.
rescale_by_camera <- function(df, metadata_df) {
  df <- df %>%
    mutate(name = str_sub(name, 1, -6))

  df1 <- df %>% filter(camera == 1)
  df2 <- df %>% filter(camera == 2)

  # Rescale and round time for camera 2
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

  # Doing this one well at a time with Purrr because of memory stuff. Gets a
  # list of unique names then runs it on one name at a time.
  unique_names <- unique(df2$name)

  # Define a function that operates on a single name
  process_name <- function(name) {
    # Filter df for the current name
    df_name <- df2[df2$name == name, ]

    # Create a new data frame with every time point for each combination.
    time_df <- expand_grid(time = seq(0, max(df_name$time), by = 0.01),
                           name = unique(df_name$name),
                           object_class = unique(df_name$object_class),
                           date = unique(df_name$date),
                           run = unique(df_name$run),
                           well = unique(df_name$well),
                           camera = unique(df_name$camera))

    # Do the interpolation from existing values.
    df_name_processed <- time_df %>%
      left_join(df_name, by = c("time", "name", "object_class", "date", "run", "well", "camera")) %>%
      group_by(name, object_class, date, run, well, camera) %>%
      arrange(time) %>%
      mutate(count = approx(time, count, time, rule = 2)$y,
             percentage = approx(time, percentage, time, rule = 2)$y) %>%
      filter(floor(time) == time) %>%
      ungroup()

    # Return the processed dataframe
    return(df_name_processed)
  }

  # Using Purrr
  df_interpolated <- unique_names %>% map(process_name) %>% list_rbind()

  output_df <- rbind(df1, df_interpolated)
  
  # Adding the accession and temp metadata
  output_df <- left_join(output_df, metadata_df, by = c("date", "run", "well"))

  return(output_df)

}

interpolated_pollen_inference <- rescale_by_camera(processed_pollen_inference, wells_to_accessions)


# Cleaning up the data for tube lengths -----------------------------------
# Calculating track lengths. This assumes the input is sorted by time.
calculate_track_length <- function(df, metadata_df){
  # Only keep the tip tips
  df <- df %>%
    filter(object_class == "tube_tip")
  
  # Adding a factor for camera
  df$camera <- NA
  df$camera[as.Date(df$date) < as.Date("2022-05-27")] <- 1
  df$camera[as.Date(df$date) >= as.Date("2022-05-27")] <- 2
  
  df$x <- as.numeric(df$x)
  df$y <- as.numeric(df$y)
  
  # First remove any tracks that grow out of the field of view or into the
  # field of view.
  # x,y image dimensions from camera 1: 2048x2048
  # x,y image dimensions from camera 2: 1600x1200
  percentage = 0.02
  margin_cam_1 = 35
  margin_cam_2 = 20
  df <- df %>%
    mutate(img_x = ifelse(camera == 1, 2048, 1600),
           img_y = ifelse(camera == 1, 2048, 1200),
           margin = ifelse(camera == 1, margin_cam_1, margin_cam_2))

  df <- df %>%
    group_by(date, run, well, track_id) %>%
    filter(all(x > margin & x < img_x - margin)) %>%
    filter(all(y > margin & y < img_y - margin)) %>%
    ungroup()

  # Next calculate the length of each track.
  df <- df %>%
    group_by(track_id) %>%
    mutate(dist = sqrt((x - lag(x))^2 + (y - lag(y))^2)) %>%
    ungroup() %>%
    group_by(camera, date, run, well, track_id) %>%
    summarise(total_length = sum(dist, na.rm = TRUE))

  # Finally scale pixel lengths to actual length in µm depending on the camera.
  # Until I can get a micrometer image, using Heinz pollen diameter on 20 µm.
  # For cam 1 Heinz pollen is 31 px diameter.
  # To get µm from cam 1 px, multiple px by 20/31, 0.6451613
  # For cam 2 Heinze pollen is 16.8 px diameter.
  # To get µm from cam 2 px, multiple px by 20/16.8, 1.190476
  df <- df %>%
    mutate(total_length_um = ifelse(camera == 1,
                                    total_length * 0.6451613,
                                    total_length * 1.190476)) %>%
    select(-total_length)

  df$run <- as.double(df$run)

  df <- left_join(df, metadata_df, by = c("date", "run", "well"))

  return(df)
}

processed_tube_lengths <- calculate_track_length(inference, wells_to_accessions)


# Quality control ---------------------------------------------------------
# # Checking for NAs
# print(processed_tube_lengths[apply(processed_tube_lengths, 1, function(row) any(is.na(row))), ])
# print(interpolated_pollen_inference[apply(interpolated_pollen_inference, 1, function(row) any(is.na(row))), ])
# 
# # Checking to see how many accessions are present in how many reps.
# interpolated_pollen_inference_qc <- interpolated_pollen_inference %>%
#   ungroup() %>%
#   group_by(accession, temp_target, name) %>%
#   summarize(n = n()) %>%
#   ungroup() %>%
#   group_by(accession, temp_target) %>%
#   summarize(n = n()) %>% # Shows the number of reps per accession at each temp
#   filter(n >= 8) %>% # What is happening?
#   ungroup() %>%
#   group_by(accession) %>%
#   summarize(n = n()) %>% # Counting how many rows each accession has. If 2 then they have >= 8 reps at each temp.
#   filter(n == 2)
# 
# # Checking the tube lengths data frame
# processed_tube_lengths_qc <- processed_tube_lengths %>%
#   mutate(name = paste0(
#     date,
#     "_run",
#     run,
#     "_",
#     temp_target,
#     "C_",
#     well
#   )) %>%
#   ungroup() %>%
#   group_by(accession, temp_target, name) %>%
#   summarize(n = n()) %>%
#   ungroup() %>%
#   group_by(accession, temp_target) %>%
#   summarize(n = n()) %>% # Shows the number of reps per accession at each temp
#   filter(n >= 8) %>% # What is happening?
#   ungroup() %>%
#   group_by(accession) %>%
#   summarize(n = n()) %>% # Counting how many rows each accession has. If 2 then they have >= 8 reps at each temp.
#   filter(n == 2)
# 
# setequal(interpolated_pollen_inference_qc$accession, processed_tube_lengths_qc$accession) # TRUE
# 
# # Checking well counts
# pollen_counting_sheet <- read_sheet("10_lG9N0wGvgOmxDGuX5PXILB7QwC7m6CuYXzi78Qe3Q")
# plate_sheet <- read_sheet("1yQ5yAKiL6BzwZ-wH-Q44RoUEwMZztTYafzdvVylq6fo") %>%
#   select(date, run, well, accession, temp_target)
# reps_df <- left_join(plate_sheet, pollen_counting_sheet, by = c("date", "run", "well"))
# reps_df_summary <- reps_df %>%
#   filter(count == "g") %>%
#   mutate(accession = replace(accession, accession == "CW00173", "CW0173")) %>%
#   mutate(accession = replace(accession, accession == "CW01000", "CW1000")) %>%
#   mutate(accession = replace(accession, accession == "CW019", "CW0019")) %>%
#   mutate(accession = replace(accession, accession == "CW01037", "CW1037")) %>%
#   mutate(accession = replace(accession, accession == "CW089", "CW0089")) %>%
#   group_by(accession, temp_target) %>%
#   summarize(n = n()) %>%
#   filter(n >= 8) %>%
#   ungroup() %>%
#   group_by(accession) %>%
#   summarize(n = n()) %>%
#   filter(n == 2)
# 
# # Finding wells that are part of the good wells in the plate worksheets that are no in the inference
# good_wells_plate_worksheets <- reps_df %>%
#   filter(count == "g") %>%
#   mutate(name = paste0(
#     date,
#     "_run",
#     run,
#     "_",
#     temp_target,
#     "C_",
#     well
#   )) %>%
#   select(name, accession)
# 
# wells_pollen_inference <- interpolated_pollen_inference %>%
#   select(name) %>%
#   unique()
# 
# missing_wells <- anti_join(good_wells_plate_worksheets, wells_pollen_inference) %>%
#   arrange(accession)


# Removing duplicated accessions ------------------------------------------
# Some accessions have partial reps from camera 1 and full reps from camera 2.
# Removing reps from camera 1 in this case.
accession_camera <- interpolated_pollen_inference %>%
  filter(accession != "CW0000") %>% # I deal with the Heinz controls later.
  group_by(accession, camera, date, run, well, temp_target) %>%
  summarize(row_count = n(), .groups = "drop") %>%
  group_by(camera, accession, temp_target) %>%
  summarize(row_count = n(), .groups = "drop") %>%
  group_by(accession, camera, temp_target) %>%
  filter(sum(row_count) >= 8) %>%
  group_by(accession, camera) %>%
  filter(n_distinct(temp_target) == 2) %>%
  select(accession, camera)

# Extract camera lists
cam_1_accessions <- accession_camera %>%
  filter(camera == 1) %>%
  pull(accession) %>%
  unique()
cam_2_accessions <- accession_camera %>%
  filter(camera == 2) %>%
  pull(accession) %>%
  unique()

# Checking for overlaps
intersect(cam_1_accessions, cam_2_accessions) # character(0)

# Removing the reps from the wrong camera in the pollen class df.
interpolated_pollen_inference <- interpolated_pollen_inference %>%
  filter((camera == 1 & accession %in% cam_1_accessions) |
         (camera == 2 & accession %in% cam_2_accessions) |
         (accession == "CW0000")) 

# Removing the reps from the wrong camera in the tube df.
processed_tube_lengths <- processed_tube_lengths %>%
  filter((camera == 1 & accession %in% cam_1_accessions) |
         (camera == 2 & accession %in% cam_2_accessions) |
         (accession == "CW0000")) 


# Writing processed data --------------------------------------------------
# Pollen classes
write.table(interpolated_pollen_inference,
  file = file.path(getwd(), "data", "processed_data", "interpolated_inference.tsv"),
  sep = '\t',
  quote = FALSE,
  row.names = FALSE
)

# Tube lengths
write.table(processed_tube_lengths,
  file = file.path(getwd(), "data", "processed_data", "processed_tube_lengths.tsv"),
  sep = '\t',
  quote = FALSE,
  row.names = FALSE
)
