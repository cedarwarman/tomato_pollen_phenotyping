# Introduction ------------------------------------------------------------
# Map accessions by geographic range and species.
# Reads accession metadata from Google Sheets, writes GPS coordinates to
# data/gps_coordinate_output/, and saves the range map to plots/mapping/.

library(dplyr)
library(ggmap)
library(googlesheets4)

pdf(NULL)
set.seed(16)

# Stamen tiles are now served by Stadia Maps (ggmap >= 4.0).
# Register a free API key once at https://client.stadiamaps.com/signup/ then run:
#   register_stadiamaps("YOUR-API-KEY", write = TRUE)
# After that the key is stored in ~/.Renviron and this line is not needed.

gs4_auth(path = "~/.credentials/google_sheets_api/service_account.json")


# Load data ---------------------------------------------------------------
metadata <- read_sheet("1V2kH8G4tfYsYqnYb6bVHpGjwwkjd9le0arGBJ2o4r8s")
metadata <- metadata %>%
  filter(Latitude != "",
         Latitude != "unknown",
         Latitude != "NULL") %>%
  mutate(latitude  = as.numeric(Latitude),
         longitude = as.numeric(Longitude),
         accession = name_CW) %>%
  select(accession, latitude, longitude, species)

dir.create(file.path(getwd(), "data", "gps_coordinate_output"), recursive = TRUE, showWarnings = FALSE)
write.table(
  metadata,
  file = file.path(getwd(), "data", "gps_coordinate_output", "accession_coordinates.tsv"),
  sep       = "\t",
  row.names = FALSE,
  quote     = FALSE
)


# Map ---------------------------------------------------------------------
# get_stadiamap() expects bbox = c(left, bottom, right, top)
bbox <- c(left   = min(metadata$longitude),
          bottom = min(metadata$latitude),
          right  = max(metadata$longitude),
          top    = max(metadata$latitude))

tomato_map <- get_stadiamap(bbox    = bbox,
                            maptype = "stamen_toner",
                            crop    = FALSE,
                            zoom    = 5)

species_colors <- c(lycopersicum    = "#FF00FF",
                    pimpinellifolium = "#1B74FA",
                    cheesmaniae     = "#11E00D",
                    galapagense     = "#FFB000")

ggmap(tomato_map) +
  geom_point(aes(x = longitude, y = latitude, color = species),
             data  = metadata,
             alpha = 0.6,
             size  = 3) +
  scale_color_manual(
    values = species_colors,
    name   = "Species",
    breaks = c("lycopersicum", "pimpinellifolium", "cheesmaniae", "galapagense"),
    labels = c("Lycopersicum", "Pimpinellifolium", "Cheesmaniae", "Galapagense"),
    guide  = guide_legend(ncol = 1)
  ) +
  labs(title = "Diversity panel accessions",
       x     = "Longitude",
       y     = "Latitude") +
  theme(
    axis.title        = element_text(size = 12, face = "bold"),
    axis.text         = element_text(size = 10, face = "bold", color = "black"),
    plot.title        = element_blank(),
    axis.line         = element_blank(),
    axis.ticks        = element_line(linewidth = 1, color = "black"),
    axis.ticks.length = unit(6, "pt"),
    plot.margin       = margin(0.5, 0.5, 0.5, 0.5, "cm"),
    panel.grid        = element_blank(),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 1.5),
    legend.position   = "bottom",
    legend.direction  = "vertical",
    legend.title      = element_text(size = 12, face = "bold", color = "black"),
    legend.text       = element_text(size = 10, face = "bold", color = "black"),
    legend.key        = element_rect(fill = "white"),
    strip.background  = element_blank(),
    strip.placement   = "outside"
  )

dir.create(file.path(getwd(), "plots", "mapping"), recursive = TRUE, showWarnings = FALSE)
ggsave(
  filename = file.path(getwd(), "plots", "mapping", "range_by_species_map.png"),
  device   = "png",
  width    = 10,
  height   = 10,
  dpi      = 400,
  units    = "in"
)
