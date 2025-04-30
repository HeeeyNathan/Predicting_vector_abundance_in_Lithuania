#' Load packages
library(googleway)
library(httr)
library(tidyverse)
library(rgeoboundaries)
library(sf)
library(janitor)
library(pacman)

# !!!!!!IMPORTANT, THIS ANALYSIS WILL REQUIRE THAT YOU HAVE YOUR OWN GOOGLE API KEY!!!!!!!!

#' Load the dipteran vector data.
# Import the data
df <- read.csv("Outputs/3_diptera_taxonomic_indices_wCorine2018_TerraClimate.csv", h = T, sep = ",", stringsAsFactors = FALSE, check.names = FALSE) |>
  clean_names() |>
  arrange(desc(waterbody_type), site_id, year) |>  # order the data.frame
  as.data.frame() # Convert tibble to dataframe because some older code does not recognize tibble

# create dataframe with unique sites and coordinates
xy <- df |>
  dplyr::select(site_id, latitude, longitude) |>
  distinct(site_id, .keep_all = TRUE) # Keeps the first occurrence of each site_id

# Split the dataset into batches (e.g., 50 coordinates per batch)
batch_size <- 50
batches <- split(xy, ceiling(seq_along(xy$latitude) / batch_size))

# Function to fetch elevation for a batch
get_elevation_batch <- function(batch) {
  google_elevation(data.frame(lat = batch$latitude, lon = batch$longitude), key = "YOUR KEY HERE") # <----- add your own google api key here
}

# Apply the function to each batch
elevation_results <- lapply(batches, get_elevation_batch)

# Combine all results into a single dataframe
# Extract the 'results' and reset row names for each batch
elevation_data <- bind_rows(
  lapply(seq_along(elevation_results), function(i) {
    # Extract 'results' and add corresponding site_id
    batch <- elevation_results[[i]]$results
    batch$site_id <- xy$site_id[(i - 1) * length(batch$elevation) + seq_along(batch$elevation)]
    as.data.frame(batch)
  })
)

# Merge elevation data back with the original dataframe
xy$elevation <- elevation_data$elevation

# rename columns for ease
xy <- xy |>
  rename(
    x = longitude,
    y = latitude,
  )

# Get Boundry for Lithuania
lithuania.shp <- geoboundaries("Lithuania")
lithuania <- st_transform(lithuania.shp, crs = 4326)

# Plot the map
plot(st_geometry(lithuania),
     col = "white",
     border = "black",
     xlim = c(21, 27),
     ylim = c(53, 57),
     main = "Lithuania Map with Points")
# Overlay points on the map
points(x = df$longitude, y = df$latitude, col = "blue", pch = 1, cex = 1)

# create sf object
xy_sf <- st_as_sf(xy, coords = c("x", "y"), crs = 4326)

# Set finer resolution for interpolation
resolution <- 0.05  # Adjust this value for finer grids (smaller = finer)

# Generate a grid covering Lithuania's bounding box
lith_bbox <- st_bbox(lithuania)
grid <- expand.grid(
  x = seq(lith_bbox["xmin"], lith_bbox["xmax"], by = resolution),
  y = seq(lith_bbox["ymin"], lith_bbox["ymax"], by = resolution)
)

# Interpolate elevation data
interp <- with(xy, akima::interp(
  x = x, y = y, z = elevation,
  xo = unique(grid$x),
  yo = unique(grid$y),
  duplicate = "mean"
))

# Convert interpolation to a data frame for plotting
interp_df <- as.data.frame(expand.grid(x = interp$x, y = interp$y))
interp_df$elevation <- as.vector(interp$z)

# Convert the interpolated data to an sf object
interp_sf <- st_as_sf(interp_df, coords = c("x", "y"), crs = st_crs(lithuania))

# Clip the interpolated raster to Lithuania's boundary
interp_clipped <- st_intersection(st_as_sf(interp_df, coords = c("x", "y"), crs = st_crs(lithuania)), lithuania)

# Extract coordinates and elevation for raster plotting
interp_clipped_df <- interp_clipped |>
  st_as_sf() |>
  st_coordinates() |>
  as.data.frame() |>
  cbind(elevation = interp_clipped$elevation)

ggplot() +
  geom_sf(data = lithuania, fill = "lightgrey", color = "black") + # Lithuania map
  geom_raster(data = interp_clipped_df, aes(x = X, y = Y, fill = elevation)) + # Shaded elevation
  scale_fill_viridis_c(option = "viridis") + # Elevation colour scale
  geom_sf(data = xy_sf, aes(), color = "red", size = 1) + # Overlay points
  labs(
    title = "Elevation Data in Lithuania",
    x = "Longitude",
    y = "Latitude",
    fill = "Elevation"
  ) +
  theme_minimal()

# Join elevation data to main dataset
df <- df |>
  left_join(xy, by = "site_id") |>
  dplyr::select(-x, -y)

# save output
write.csv(df, "Outputs/4_diptera_taxonomic_indices_wCorine2018_TerraClimate_elevation.csv", row.names = F)


#### Predcition data #####

# Data for prediction using our model. We will use the coordinates from a site 16.5km from Rybachy, very near the Lithuanian border
# Rybachy on the Curonian Lagoon: 55.153730, 20.857780 Decimal degree; 490937.22, 6111907.79 UTM
# Our point is at: 55.285238, 20.970616 Decimal degree; 498133.71, 6126533.56 UTM
# We use km instead of meters so we need to divide by 1000. The euclidean distance between the sites is 16.41 km

predict_data <- read.csv("Outputs/6_prediction_data.csv") |>
  mutate(site_id        = as.factor(site_id),
         sample_id      = as.factor(sample_id),
         waterbody_type = as.factor(waterbody_type),
         date           = as.Date(date, format = "%Y-%m-%d"))

# create dataframe with unique sites and coordinates
xy <- predict_data |>
  dplyr::select(site_id, latitude, longitude) |>
  distinct(site_id, .keep_all = TRUE) # Keeps the first occurrence of each site_id

# Split the dataset into batches (e.g., 50 coordinates per batch)
batch_size <- 2
batches <- split(xy, ceiling(seq_along(xy$latitude) / batch_size))

# Function to fetch elevation for a batch
get_elevation_batch <- function(batch) {
  google_elevation(data.frame(lat = batch$latitude, lon = batch$longitude), key = "YOUR KEY HERE") # <----- add your own google api key here
}

# Apply the function to each batch
elevation_results <- lapply(batches, get_elevation_batch)

# Combine all results into a single dataframe
# Extract the 'results' and reset row names for each batch
elevation_data <- bind_rows(
  lapply(seq_along(elevation_results), function(i) {
    # Extract 'results' and add corresponding site_id
    batch <- elevation_results[[i]]$results
    batch$site_id <- xy$site_id[(i - 1) * length(batch$elevation) + seq_along(batch$elevation)]
    as.data.frame(batch)
  })
)

# Merge elevation data back with the original dataframe
xy$elevation <- elevation_data$elevation

# rename columns for ease
xy <- xy |>
  rename(
    x = longitude,
    y = latitude,
  )

# Join elevation data to main dataset
predict_data <- predict_data |>
  left_join(xy, by = "site_id") |>
  dplyr::select(-x, -y)

# save output
write.csv(predict_data, "Outputs/6_prediction_data.csv", row.names = F)

###############################################################################################################
# CLEAN UP WORKSPACE
rm(list = ls())       # Remove all objects from environment
gc()                  # Frees up unused memory
p_unload(all)         # Unload all loaded packages
graphics.off()        # Close all graphical devices
cat("\014")           # Clear the console
# Clear mind :)
