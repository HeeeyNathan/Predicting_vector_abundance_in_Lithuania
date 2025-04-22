# Panel plot for paper

#==================== Sampling sites ===================
#====== Load packages ======
library(ggplot2)
library(sf)
library(rgeoboundaries)
library(grid)
library(readr)
library(ggspatial)
library(pacman)

#====== Load theme for plotting ======
My_theme <- theme(panel.background = element_blank(),
                  panel.border = element_rect(fill = NA, linewidth = 1.25),
                  strip.background = element_rect(fill = "white",
                                                  color = "white", linewidth = 1.25),
                  legend.position = "bottom",
                  text = element_text(size = 16))

#====== Load data ======
Unique_sites <- read.csv("Outputs/5_unique_sites_for_plotting.csv", h = T, sep = ",", stringsAsFactors = T)
# Load and reproject water bodies from geodatabase
gdb_path <- "Additional data/GeoDatabase/UETK_2024-05-02.gdb"
rivers <- st_read(gdb_path, layer = "upes_l") |>
  st_transform(4326)  # Reproject to WGS84
lakes <- st_read(gdb_path, layer = "ezerai_tvenkiniai") |>
  st_transform(4326)  # Reproject to WGS84

#====== Get Lithuania boundary ======
lithuania_sf <- geoboundaries("Lithuania", adm_lvl = 0)
lithuania_spat <- terra::vect(lithuania_sf)

#====== Get Europe boundary ======
europe <- geoboundaries(c("Lithuania", "Poland", "Latvia", "Estonia", "Belarus",
                          "Russia", "Ukraine", "Germany", "Denmark", "Sweden", "Finland",
                          "Czech Republic", "Slovakia", "Austria", "Hungary", "Romania", "Bulgaria",
                          "Moldova", "Slovenia", "Croatia", "Bosnia and Herzegovina", "Serbia",
                          "Montenegro", "Albania", "North Macedonia", "Kosovo", "Greece", "Turkey",
                          "Cyprus", "Malta", "Italy", "France", "Spain", "Portugal", "Andorra",
                          "Monaco", "San Marino", "Vatican City", "Switzerland", "Liechtenstein",
                          "Luxembourg", "Belgium", "Netherlands", "Ireland", "United Kingdom",
                          "Norway", "Iceland"))
europe <- st_transform(europe, crs = 4326)

# Create main plot
sampling_sites <- ggplot() +
  # Add rivers layer with much thinner lines
  geom_sf(data = rivers,
          color = "lightblue",
          alpha = 0.25,
          size = 0.01) +
  # Add lakes layer
  geom_sf(data = lakes,
          fill = "lightblue",
          color = NA,
          alpha = 0.75,
          size = 0.01) +
  # Base layer - Lithuania boundary
  geom_sf(data = lithuania_sf, color = "black", fill = NA, linewidth = 0.5) +
  # Add sampling points
  geom_point(aes(x = Long,
                 y = Lat,
                 color = waterbody_type,
                 shape = waterbody_type),
             data = Unique_sites,
             size = 1.5) +
  coord_sf(xlim = c(20.9, 26.9), ylim = c(53.85, 56.5)) +
  scale_color_manual(name = "",
                    values = c("#440154", "#47A635")) +
  scale_shape_manual(name = "",
                    values = c(19, 19)) +
  labs(x = "Longitude",
       y = "Latitude") +
  My_theme +
  theme(
    legend.position = "bottom",
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    legend.box.spacing = unit(0.5, "lines"),  # Reduced from 1.5 to 0.5 to bring legend closer
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black")  # Ensuring tick marks are visible
  ) +
  guides(color = guide_legend(override.aes = list(size = 7, shape = 19, fill = c("#440154", "#47A635"))))

# Create inset map with matching style
inset <- ggplot() +
  geom_sf(data = europe, color = "black", fill = "white") +
  geom_sf(data = lithuania_sf, fill = "#440154", alpha = 0.5) +
  theme_void() +
  theme(panel.background = element_rect(fill = "lightblue"),
        panel.border = element_rect(fill = NA, linewidth = 0.5)) +
  coord_sf(xlim = c(10, 39), ylim = c(45, 65))

# Convert inset to grob
inset_grob <- ggplotGrob(inset)

# Add inset to main plot
sampling_sites <- sampling_sites +
  annotation_custom(inset_grob,
                    xmin = 19.99, xmax = 23,
                    ymin = 53.85, ymax = 54.98)

# Add north arrow and scale bar
sampling_sites <- sampling_sites +
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    pad_x = unit(0.25, "cm"),
    pad_y = unit(0.25, "cm"),
    style = north_arrow_fancy_orienteering
  ) +
  annotation_scale(
    location = "bl",
    width_hint = 0.25
  )

ggsave("Plots/Figure2_Sampling_sites_wWater.png", width = 10, height = 8, dpi = 300, bg = "white")

#==================== Average EQC of sampling sites ===================
# Define WFD colors for EQC classes
wfd_colors <- c("High" = "#0000FF",     # Pure blue
                "Good" = "#00AA00",     # Darker green
                "Moderate" = "#FFD700", # Gold
                "Poor" = "#FF6600",     # Darker orange
                "Bad" = "#FF0000")      # Red

# Create the map
average_eqc <- ggplot() +
  # Base layer - Lithuania boundary
  geom_sf(data = lithuania_sf, color = "black", fill = NA, linewidth = 0.5) +
  # Add sampling points - colored by the existing ave_eqc_cat variable
  geom_point(data = Unique_sites,
             aes(x = Long,
                 y = Lat,
                 # shape = waterbody_type,
                 color = avg_eqc_cat),
             size = 1.5) +

  # Set coordinate system and map extent
  coord_sf(xlim = c(20.9, 26.9), ylim = c(53.85, 56.5)) +

  # Set color scale to use WFD colors
  scale_color_manual(name = "Average Ecological Quality Class",
                    values = wfd_colors) +

  # # Set shapes - squares for lakes, circles for rivers
  # scale_shape_manual(name = "Waterbody Type",
  #                   values = c("Lake" = 17, "River" = 19)) +

  # Add labels
  labs(x = "Longitude",
       y = "Latitude") +

  # Apply the theme
  My_theme +

  # Additional theme adjustments
  theme(
    legend.position = "bottom",
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    legend.box.spacing = unit(0.5, "lines"),
    legend.box = "vertical",  # Stack legends vertically
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    legend.title.position = "top",
    legend.title = element_text(hjust = 0.5)
  ) +

  # Customize legend appearance
  guides(
    color = guide_legend(override.aes = list(size = 4), order = 1),
    shape = guide_legend(override.aes = list(size = 4), order = 2)
  )
average_eqc

#==================== Elevation ===================
#====== Load packages ======
library(elevatr)
library(kableExtra)
library(rgeoboundaries)
library(sf)
library(raster)
library(ggplot2)
library(viridis)

#====== Extract elevation data ======
ll_proj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
elev <- elevatr::get_elev_point(pt_df, prj = ll_proj)
elev |>
  kable() |>
  kable_styling(bootstrap_options = c("striped", "hover"))

#====== Extract boundaries ======
lithuania_sf <- rgeoboundaries::geoboundaries("Lithuania", adm_lvl = 0)
elevation_data <- elevatr::get_elev_raster(locations = lithuania_sf, z = 9, clip = "locations")

#====== Link elevation data ======
elevation_data <- as.data.frame(elevation_data, xy = TRUE)
colnames(elevation_data)[3] <- "elevation"
# remove rows of data frame with one or more NA's,using complete.cases
elevation_data <- elevation_data[complete.cases(elevation_data), ]

#====== Plot ======

max_elev <- ceiling(max(elevation_data$elevation, na.rm = TRUE)/100)*100

# contour plot
elev_range <- range(elevation_data$elevation, na.rm = TRUE)

# Create a dummy continuous dataset for the legend
legend_data <- data.frame(
  x = seq(elev_range[1], max_elev, length.out = 100),
  y = rep(1, 100)
)

elevation_plot <- ggplot() +
  # The discrete filled contours for the main plot
  stat_contour_filled(
    data = elevation_data,
    aes(x = x, y = y, z = elevation),
    breaks = seq(floor(elev_range[1]), max_elev, by = 50)
  ) +
  # Add a hidden continuous scale just for the legend
  geom_point(
    data = legend_data,
    aes(x = x, y = y, color = x),
    alpha = 0
  ) +
  geom_sf(data = lithuania_sf, color = "black", fill = NA, linewidth = 0.5) +
  coord_sf(xlim = c(20.9, 26.9), ylim = c(53.85, 56.5)) +
  # Replace with custom color palette including mid-color
  scale_fill_manual(
    values = colorRampPalette(c("#8C7853", "#D2C29D", "#097969"))(length(seq(floor(elev_range[1]), max_elev, by = 50)) - 1),
    guide = "none"
  ) +
  # Add mid-point to the continuous color scale
  scale_color_gradientn(
    colors = c("#8C7853", "#D2C29D", "#097969"),
    limits = c(elev_range[1], max_elev),
    name = "m",
    breaks = c(0, 100, 200, 300),
    guide = guide_colorbar(
      barwidth = 20,
      barheight = 0.5,
      title.position = "top",
      title.hjust = 0.5,
      label.hjust = 0.5
    )
  ) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal() +
  My_theme +
  theme(
    legend.position = "bottom",
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    legend.box.spacing = unit(0.5, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black")
  )

#==================== Temperature ===================
#====== Load packages ======
library(raster)
library(tidyverse)
library(magrittr)
library(rgeoboundaries)
library(sf)
library(terra)


#====== Get Lithuania boundary ======
lithuania_sf <- geoboundaries("Lithuania", adm_lvl = 0)
lithuania_spat <- terra::vect(lithuania_sf)

#====== Load and process TerraClimate data ======
# Create a list of file paths for all years
tmin_files <- paste0("Additional data/TerraClimate/TerraClimate_tmin_", 2013:2022, ".nc")

# Read all files as a raster stack
tmin_stack <- terra::rast(tmin_files)

# Calculate mean across all years
tmin_mean <- terra::mean(tmin_stack)

# Mask to Lithuania
tmin_mean_Lithuania <- terra::mask(tmin_mean, lithuania_spat)

# Clip to Lithuania boundary and convert to dataframe
tmin_mean_Lithuania <- terra::crop(tmin_mean, lithuania_spat, mask = TRUE)
tmin_mean_LT_df <- as.data.frame(tmin_mean_Lithuania, xy = TRUE, na.rm = TRUE)

# Remove any NA values
tmin_mean_LT_df <- tmin_mean_LT_df[complete.cases(tmin_mean_LT_df), ]

# Create a higher resolution template raster
current_res <- res(tmin_mean)
smooth_template <- terra::rast(ext(tmin_mean),
                             resolution = current_res/5,
                             crs = crs(tmin_mean))

# Resample to higher resolution using bilinear interpolation
tmin_smooth <- terra::resample(tmin_mean, smooth_template, method = "cubic")

# Clip to Lithuania boundary and convert to dataframe
tmin_mean_Lithuania <- terra::crop(tmin_smooth, lithuania_spat, mask = TRUE)
tmin_mean_LT_df <- as.data.frame(tmin_mean_Lithuania, xy = TRUE, na.rm = TRUE)

# Remove any NA values
tmin_mean_LT_df <- tmin_mean_LT_df[complete.cases(tmin_mean_LT_df), ]

# Calculate rounded min and max temperatures for nice breaks
max_temp <- ceiling(max(tmin_mean_LT_df$mean, na.rm = TRUE))
min_temp <- floor(min(tmin_mean_LT_df$mean, na.rm = TRUE))

# First get the range of temperatures for proper color mapping
temp_range <- range(tmin_mean_LT_df$mean)
# Create a dummy continuous dataset for the legend
legend_data <- data.frame(
  x = seq(temp_range[1], temp_range[2], length.out = 100),
  y = rep(1, 100)
)
temperature_plot <- ggplot() +
  # The discrete filled contours for the main plot
  stat_contour_filled(
    data = tmin_mean_LT_df,
    aes(x = x, y = y, z = mean),
    breaks = seq(floor(temp_range[1]),
                ceiling(temp_range[2]),
                by = 0.5)  # Changed from 0.5 to 1
  ) +
  # Add a hidden continuous scale just for the legend
  geom_point(
    data = legend_data,
    aes(x = x, y = y, color = x),
    alpha = 0
  ) +
  geom_sf(data = lithuania_sf, color = "black", fill = NA, linewidth = 0.5) +
  coord_sf(xlim = c(20.9, 26.9), ylim = c(53.85, 56.5)) +
  # Replace viridis with yellow-to-red for discrete fill
  scale_fill_manual(
    values = colorRampPalette(c("#FFFF00", "#FF0000"))(length(seq(floor(temp_range[1]), ceiling(temp_range[2]), by = 0.5)) - 1),
    guide = "none"
  ) +
  # Replace viridis with yellow-to-red for continuous color
  scale_color_gradient(
    low = "#FFFF00",  # Yellow
    high = "#FF0000",  # Red
    name = "°C",
    breaks = seq(floor(temp_range[1]),
                ceiling(temp_range[2]),
                by = 0.5),  # Changed to match contour intervals
    guide = guide_colorbar(
      barwidth = 20,
      barheight = 0.5,
      title.position = "top",
      title.hjust = 0.5,
      label.hjust = 0.5
    )
  ) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal() +
  My_theme +
  theme(
    legend.position = "bottom",
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    legend.box.spacing = unit(0.5, "lines"),  # Reduced from 1.5 to 0.5 to bring legend closer
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black")  # Ensuring tick marks are visible
  )

#==================== Precipitation ===================
#====== Load packages ======
library(raster)
library(tidyverse)
library(magrittr)
library(rgeoboundaries)
library(sf)
library(terra)

#====== Get Lithuania boundary ======
lithuania_sf <- geoboundaries("Lithuania", adm_lvl = 0)
lithuania_spat <- terra::vect(lithuania_sf)

#====== Load and process TerraClimate data ======
# Create a list of file paths for all years
ppt_files <- paste0("Additional data/TerraClimate/TerraClimate_ppt_", 2013:2022, ".nc")

# Read all files as a raster stack
ppt_stack <- terra::rast(ppt_files)

# Calculate mean across all years
ppt_mean <- terra::mean(ppt_stack)

# Mask to Lithuania
ppt_mean_Lithuania <- terra::mask(ppt_mean, lithuania_spat)

# Clip to Lithuania boundary and convert to dataframe
ppt_mean_Lithuania <- terra::crop(ppt_mean, lithuania_spat, mask = TRUE)
ppt_mean_LT_df <- as.data.frame(ppt_mean_Lithuania, xy = TRUE, na.rm = TRUE)

# Remove any NA values
ppt_mean_LT_df <- ppt_mean_LT_df[complete.cases(ppt_mean_LT_df), ]

# Create a higher resolution template raster
current_res <- res(ppt_mean)
smooth_template <- terra::rast(ext(ppt_mean),
                             resolution = current_res/5,
                             crs = crs(ppt_mean))

# Resample to higher resolution using bilinear interpolation
ppt_smooth <- terra::resample(ppt_mean, smooth_template, method = "cubic")

# Clip to Lithuania boundary and convert to dataframe
ppt_mean_Lithuania <- terra::crop(ppt_smooth, lithuania_spat, mask=TRUE)
ppt_mean_LT_df <- as.data.frame(ppt_mean_Lithuania, xy = TRUE, na.rm = TRUE)

# Remove any NA values
ppt_mean_LT_df <- ppt_mean_LT_df[complete.cases(ppt_mean_LT_df), ]

# Get the range of precipitation for proper color mapping
ppt_range <- range(ppt_mean_LT_df$mean)

# Expand the range slightly to ensure complete coverage
ppt_min <- floor(ppt_range[1])  # Go 0mm below the minimum
ppt_max <- ceiling(ppt_range[2]) + 1  # Go 1mm above the maximum

# Define breaks with adequate coverage
breaks_seq <- seq(ppt_min, ppt_max, by = 2)

# Create the plot with adjusted breaks
precipitation_plot <- ggplot() +
  # The discrete filled contours for the main plot with expanded breaks
  stat_contour_filled(
    data = ppt_mean_LT_df,
    aes(x = x, y = y, z = mean),
    breaks = breaks_seq
  ) +
  # Add a hidden continuous scale just for the legend
  geom_point(
    data = data.frame(
      x = seq(ppt_range[1], ppt_range[2], length.out = 100),
      y = rep(1, 100)
    ),
    aes(x = x, y = y, color = x),
    alpha = 0
  ) +
  geom_sf(data = lithuania_sf, color = "black", fill = NA, linewidth = 0.5) +
  coord_sf(xlim = c(20.9, 26.9), ylim = c(53.85, 56.5)) +
  # Use more colors to ensure coverage
  scale_fill_manual(
    values = colorRampPalette(c("#ADD8E6", "#0000FF"))(length(breaks_seq) - 1),
    guide = "none"
  ) +
  scale_color_gradient(
    low = "#ADD8E6",  # Light blue
    high = "#0000FF",  # Blue
    name = "mm",
    limits = c(ppt_min, ppt_max),  # Match the limits with the breaks
    breaks = seq(ppt_min, ppt_max, by = 6),  # Fewer breaks for a cleaner legend
    guide = guide_colorbar(
      barwidth = 20,
      barheight = 0.5,
      title.position = "top",
      title.hjust = 0.5,
      label.hjust = 0.5
    )
  ) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal() +
  My_theme +
  theme(
    legend.position = "bottom",
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    legend.box.spacing = unit(0.5, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black")
  )

#==================== Wind Speed ===================
#====== Load packages ======
library(raster)
library(tidyverse)
library(magrittr)
library(rgeoboundaries)
library(sf)
library(terra)

#====== Get Lithuania boundary ======
lithuania_sf <- geoboundaries("Lithuania", adm_lvl = 0)
lithuania_spat <- terra::vect(lithuania_sf)

#====== Load and process TerraClimate data ======
# Create a list of file paths for all years
ws_files <- paste0("Additional data/TerraClimate/TerraClimate_ws_", 2013:2022, ".nc")

# Read all files as a raster stack
ws_stack <- terra::rast(ws_files)

# Calculate mean across all years
ws_mean <- terra::mean(ws_stack)

# Mask to Lithuania
ws_mean_Lithuania <- terra::mask(ws_mean, lithuania_spat)

# Clip to Lithuania boundary and convert to dataframe
ws_mean_Lithuania <- terra::crop(ws_mean, lithuania_spat, mask = TRUE)
ws_mean_LT_df <- as.data.frame(ws_mean_Lithuania, xy = TRUE, na.rm = TRUE)

# Remove any NA values
ws_mean_LT_df <- ws_mean_LT_df[complete.cases(ws_mean_LT_df), ]

# Create a higher resolution template raster
current_res <- res(ws_mean)
smooth_template <- terra::rast(ext(ws_mean),
                             resolution = current_res/5,
                             crs = crs(ws_mean))

# Resample to higher resolution using bilinear interpolation
ws_smooth <- terra::resample(ws_mean, smooth_template, method = "cubic")

# Clip to Lithuania boundary and convert to dataframe
ws_mean_Lithuania <- terra::crop(ws_smooth, lithuania_spat, mask = TRUE)
ws_mean_LT_df <- as.data.frame(ws_mean_Lithuania, xy = TRUE, na.rm = TRUE)

# Remove any NA values
ws_mean_LT_df <- ws_mean_LT_df[complete.cases(ws_mean_LT_df), ]

# Get the range of wind speeds for proper color mapping
ws_range <- range(ws_mean_LT_df$mean)

# Create a dummy continuous dataset for the legend
legend_data <- data.frame(
  x = seq(ws_range[1], ws_range[2], length.out = 100),
  y = rep(1, 100)
)

# Get the range of wind speeds for proper color mapping
ws_range <- range(ws_mean_LT_df$mean)

# Create a dummy continuous dataset for the legend
legend_data <- data.frame(
  x = seq(ws_range[1], ws_range[2], length.out = 100),
  y = rep(1, 100)
)

wind_speed_plot <- ggplot() +
  # The discrete filled contours for the main plot
  stat_contour_filled(
    data = ws_mean_LT_df,
    aes(x = x, y = y, z = mean),
    breaks = seq(floor(ws_range[1]),
                ceiling(ws_range[2]),
                by = 0.25)  # Adjusted interval for wind speed
  ) +
  # Add a hidden continuous scale just for the legend
  geom_point(
    data = legend_data,
    aes(x = x, y = y, color = x),
    alpha = 0
  ) +
  geom_sf(data = lithuania_sf, color = "black", fill = NA, linewidth = 0.5) +
  coord_sf(xlim = c(20.9, 26.9), ylim = c(53.85, 56.5)) +
  # Replace viridis with three-color green gradient for discrete fill
  scale_fill_manual(
    values = colorRampPalette(c("#C5E8B7", "#64A866", "#0B6623"))(length(seq(floor(ws_range[1]), ceiling(ws_range[2]), by = 0.25)) - 1),
    guide = "none"
  ) +
  # Replace viridis with three-color green gradient for continuous color
  scale_color_gradientn(
    colors = c("#C5E8B7", "#64A866", "#0B6623"),  # Light, medium, and dark green
    name = "m/s",  # Changed units to m/s for wind speed
    breaks = seq(floor(ws_range[1]),
                ceiling(ws_range[2]),
                by = 0.5),  # Adjusted to match contour intervals
    guide = guide_colorbar(
      barwidth = 20,
      barheight = 0.5,
      title.position = "top",
      title.hjust = 0.5,
      label.hjust = 0.5
    )
  ) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal() +
  My_theme +
  theme(
    legend.position = "bottom",
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    legend.box.spacing = unit(0.5, "lines"),  # Reduced from 1.5 to 0.5 to bring legend closer
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black")  # Ensuring tick marks are visible
  )

#==================== CORINE Land Cover ===================
#====== Load packages ======
library(raster)
library(tidyverse)
library(rgeoboundaries)
library(sf)
library(terra)
library(tidyterra)

# !!!!!! IMPORTANT: YOU WILL NEED TO DOWNLOAD THE CORINE LANDCOVER DATA YOURSELF AND LINK IT TO YOUR OWN WORKING DIRECTORY!!!!!!!!!!!!!!!!!!!

#====== Get Lithuania boundary ======
lithuania_sf <- geoboundaries("Lithuania", adm_lvl = 0)
lithuania_spat_3035 <- terra::project(terra::vect(lithuania_sf), "EPSG:3035")

#====== Load and process CORINE data ======
# Read CORINE raster
corine_raster <- terra::rast("C:/Users/natha/OneDrive/University of Johannesburg/Data/Lithuanian data/Landcover data/LandUse/Corine2018/u2018_clc2018_v2020_20u1_raster100m/DATA/U2018_CLC2018_V2020_20u1.tif")
print(cats(corine_raster))

# Load legends
LegendCorine_2012 <- as_tibble(read.table("Additional data/Corine Landcover/clc_legend.csv", h = T, sep = ";")) |>
  select(CLC_CODE, LABEL1, LABEL2)

LegendCorine_2018 <- as_tibble(read.csv("Additional data/Corine Landcover/CLC2018_CLC2018_V2018_20_QGIS.txt", h = FALSE, sep = ",")) |>
  rename(CLC_CODE = V1, RED = V2, GREEN = V3, BLUE = V4, ALPHA = V5, LABEL3 = V6) |>
  mutate(GRID_CODE = row_number()) |>
  left_join(LegendCorine_2012, by = "CLC_CODE") |>
  select(GRID_CODE, CLC_CODE, RED, GREEN, BLUE, ALPHA, LABEL1, LABEL2, LABEL3)

# Crop to Lithuania
corine_lt <- terra::crop(corine_raster, lithuania_spat_3035)
corine_lt <- terra::mask(corine_lt, lithuania_spat_3035)

# Project to WGS84 (EPSG:4326)
corine_lt_4326 <- terra::project(corine_lt, "EPSG:4326")

# Create a levels dataframe for LABEL1
levels_df <- data.frame(
  ID = 1:nrow(LegendCorine_2018),
  LABEL1 = LegendCorine_2018$LABEL1
)

# Set the new levels
levels(corine_lt_4326) <- levels_df

# Create color scheme based on LABEL1
unique_labels <- unique(LegendCorine_2018$LABEL1)
label1_colors <- c(
  "Artificial surfaces" = "#ff0000",            # Reddish color for urban/built-up areas
  "Agricultural areas" = "#ffffa8",             # Wheat/crop color for agriculture
  "Forest and semi natural areas" = "#00a600",  # Dark green for forests
  "Wetlands" = "#a6a6ff",                       # Light green-blue for wetlands
  "Water bodies" = "#00ccf2",                   # Dark blue for water
  "NODATA" = "#ffffff"                         # White for no data
)

landcover_plot <- ggplot() +
  geom_spatraster(data = corine_lt_4326) +
  geom_sf(data = lithuania_sf, color = "black", fill = NA, linewidth = 0.5) +
  coord_sf(xlim = c(20.9, 26.9), ylim = c(53.85, 56.5)) +
  scale_fill_manual(
    name = "Land Cover Type",
    values = label1_colors,
    breaks = unique_labels,
    na.value = "white"
  ) +
  labs(
    title = "CORINE Land Cover in Lithuania",
    caption = "corine 2018",
    x = "Longitude",
    y = "Latitude"
  ) +
  My_theme +
  theme(
    legend.position = "bottom",
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    legend.box.spacing = unit(0.5, "lines"),  # Reduced from 1.5 to 0.5 to bring legend closer
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    panel.background = element_rect(fill = "white"),
    panel.grid = element_blank()
  )
landcover_plot

#==================== Combining plots ===================
#====== Load packages ======
library(patchwork)

# Function to modify legend appearance and remove titles
modify_legend <- function(plot) {
  plot +
    theme(
      legend.position = "bottom",
      legend.key.size = unit(0.5, "cm"),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.margin = margin(t = 3, r = 0, b = 0, l = 0),
      legend.box.margin = margin(t = 3, r = 0, b = 0, l = 0),
      # Remove plot title
      plot.title = element_blank()
    )
}

# Modify the landcover plot legend with two columns
landcover_modified <- landcover_plot +
  guides(
    fill = guide_legend(
      ncol = 2,
      byrow = TRUE
    )
  )

#====== Create add_plot_label function with top-right positioning ======
add_plot_label <- function(plot, label) {
  plot +
    annotate(
      "text",
      x = Inf,  # Changed from -Inf to Inf for right side
      y = Inf,
      label = label,
      hjust = 1.5,  # Changed to move label inward from right
      vjust = 1.5,
      size = 7,     # Increased from 5 to 7
      fontface = "bold"
    )
}

plots_list <- list(
  # # Plot A (top-left): Show y-axis, hide x-axis and x-label
  # add_plot_label(modify_legend(sampling_sites +
  #   labs(x = NULL, y = NULL, caption = NULL) +
  #   theme(
  #     axis.text.x = element_blank(),
  #     axis.title.x = element_blank())), "A"),

  # Plot A (top-left): Show y-axis, hide x-axis and x-label
  add_plot_label(modify_legend(average_eqc +
    labs(x = NULL, y = NULL, caption = NULL) +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank())), "A"),

  # Plot B (top-right): Hide both axes and labels
  add_plot_label(modify_legend(elevation_plot +
    labs(x = NULL, y = NULL, caption = NULL) +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())), "B"),

  # Plot C (middle-left): Show y-axis, hide x-axis and x-label
  add_plot_label(modify_legend(temperature_plot +
    labs(x = NULL, y = NULL, caption = NULL) +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank())), "C"),

  # Plot D (middle-right): Hide both axes and labels
  add_plot_label(modify_legend(precipitation_plot +
    labs(x = NULL, y = NULL, caption = NULL) +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())), "D"),

  # Plot E (bottom-left): Show both axes
  add_plot_label(modify_legend(wind_speed_plot +
    labs(x = NULL, y = NULL, caption = NULL)), "E"),

  # Plot F (bottom-right): Show x-axis, hide y-axis and y-label
  add_plot_label(modify_legend(landcover_modified +
    labs(x = NULL, y = NULL, caption = NULL) +
    theme(
      axis.text.y = element_blank(),
      axis.title.y = element_blank())), "F")
)

# Create layout using patchwork
combined_plot <- (
  plots_list[[1]] + plots_list[[2]] +
  plots_list[[3]] + plots_list[[4]] +
  plots_list[[5]] + plots_list[[6]]
) +
  plot_layout(
    ncol = 2,
    guides = "keep"
  ) +
  plot_annotation(
    theme = theme(
      axis.title.x = element_text(size = 12, margin = margin(t = 20)),
      axis.title.y = element_text(size = 12, margin = margin(r = 20)),
      # Ensure no title in combined plot
      plot.title = element_blank()
    )
  ) &
  xlab("Longitude") &
  ylab("Latitude")

ggsave("Plots/Figure3_covariate_panel_plot.png", plot = combined_plot, width = 10.5, height = 14, dpi = 300)

###############################################################################################################
# CLEAN UP WORKSPACE
rm(list = ls())       # Remove all objects from environment
gc()                  # Frees up unused memory
p_unload(all)         # Unload all loaded packages
graphics.off()        # Close all graphical devices
cat("\014")           # Clear the console
# Clear mind :)
