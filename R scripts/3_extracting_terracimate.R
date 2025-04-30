# Load the required packages
library(tidyverse) # General processing and data formating
library(QBMS) # for extracting TerraClimate data
library(lubridate) # for formating dates and timeseries data
library(purrr) # for programing funcitons
library(pacman) # For cleaning workspace

# Load data
Data <- as_tibble(read.csv("Outputs/2_diptera_taxonomic_indices_wCorine2018.csv", sep = ",", h = T)) |>
    mutate(sample_id = as.factor(sample_id)) |>
    mutate(site_id = as.factor(site_id)) |>
    mutate(site_code = as.factor(site_code)) |>
    mutate(date = as.Date(date, format = "%d/%m/%Y")) |>
    mutate(waterbody_type = as.factor(waterbody_type)) |>
    mutate(latitude = as.double(latitude)) |>
    mutate(longitude = as.double(longitude)) |>
    arrange(desc(waterbody_type), site_id, year)

# Calculate the observation period and number of repeated samplings for each site
sampling_info <- Data |>
  group_by(site_id) |>
  summarize(
    observation_period = max(year) - min(year) + 1,
    sampling_events = n_distinct(year))

# Join data
Data <- Data |>
  left_join(sampling_info, by = "site_id") |>
  arrange(desc(waterbody_type), site_id, year)

# Remove duplicated rows
Data_filtered <- Data |>
  dplyr::select(sample_id, site_code, site_id, date, day, month, year, waterbody_type, latitude, longitude) |>
  distinct()

# Extracting the data
lon <- Data_filtered$longitude
lat <- Data_filtered$latitude

# Download offline versions of the TerraClimate data to aid in speeding up extraction
# ini_terraclimate('2000-01-01', '2024-12-31', c('ppt', 'q', 'tmax', 'tmin', 'ws'), data_path = "Additional data/TerraClimate/") # downloads offline datasets for the years and variables specified

# Get TerraClimate data for a given coordinate(s)
# terraclimate_data <- get_terraclimate(lat, lon, '2009-01-01', '2024-12-31', c('ppt', 'q', 'tmax', 'tmin', 'ws'), offline = TRUE, data_path = "Additional data/TerraClimate/") # links downloaded data to sets of coordinates and years provided

# saveRDS(terraclimate_data, "Additional data/TerraClimate/linked_terraclimate_data.RDS") # save extracted data to save time

terraclimate_data <- readRDS("Additional data/TerraClimate/linked_terraclimate_data.RDS") # saves time to always extract these data
str(terraclimate_data$climate) # check structure of linked climate data
print(terraclimate_data$climate[[1]]) # check climate data for site LTR1_4.10.2010

# Calculate means for an individual sample_id to test if the code below works as intended
## means calculated from the 12 months preceeding biological sampling e.g., if sampling for LTR1_2010 was conducted on the 4.10.2010, then calculate mean environmental data from 01.10.2009 to 30.09.2010

# extract climate data for the first row of the diptera_filtered dataset (i.e., LTR1_4.10.2010)
LTR1_4.10.2010 <- as_tibble(terraclimate_data$climate[[1]])
# Specific date and the preceding 12 months period
specific_date <- ymd(Data_filtered$date[1]) # date of biological sampling
start_date <- (specific_date %m-% months(12)) |>
  floor_date(unit = "month") # identifies the start date for mean calculations (12 months prior to biological sampling)
end_date <- (specific_date %m-% months(1)) |>
  floor_date(unit = "month") |>
  ceiling_date(unit = "month") %m-% days(1) # identifies the end date for mean calculations (the month prior to biological sampling)
# Create a Date column in your data
LTR1_4.10.2010 <- LTR1_4.10.2010 |>
  mutate(date = make_date(year, month, 1))
# Filter data for the 12 months preceding the specific date
LTR1_4.10.2010_filtered <- LTR1_4.10.2010 |>
  filter(date >= start_date & date <= end_date)
# Calculate means
LTR1_4.10.2010_means <- LTR1_4.10.2010_filtered |>
  summarise(across(c(ppt, q, tmax, tmin, ws), ~ round(mean(.x, na.rm = TRUE), 3))) # means calculated from the 12 months prior to biological sampling
print(LTR1_4.10.2010_means) # inspect calculated means

# streamlining the process
process_climate_data <- function(climate_data, sample_date) {
  specific_date <- ymd(sample_date)
  start_date <- (specific_date %m-% months(12)) |> floor_date(unit = "month")
  end_date <- (specific_date %m-% months(1)) |> floor_date(unit = "month") |> ceiling_date(unit = "month") %m-% days(1)

  climate_data |>
    mutate(date = make_date(year, month, 1)) |>
    filter(date >= start_date & date <= end_date) |>
    summarise(across(c(ppt, q, tmax, tmin, ws), ~ round(mean(.x, na.rm = TRUE), 3)))
}

# Calculate means with the function defined above
calc_means <- map2_df(terraclimate_data$climate, Data_filtered$date, process_climate_data)
# Add sample_id to the results
calc_means <- bind_cols(Data_filtered |> dplyr::select(sample_id), calc_means)
# Print the final dataframe
print(calc_means)

# Do the same test but on a random row to test if the code is correct
random_index <- sample(nrow(Data_filtered), 1)
# extract climate data for the first row of the diptera_filtered dataset
random_sample <- as_tibble(terraclimate_data$climate[[random_index]])
# Specific date and the preceding 12 months period
specific_date <- ymd(Data_filtered$date[random_index]) # date of biological sampling
start_date <- (specific_date %m-% months(12)) |>
  floor_date(unit = "month") # identifies the start date for mean calculations (12 months prior to biological sampling)
end_date <- (specific_date %m-% months(1)) |>
  floor_date(unit = "month") |>
  ceiling_date(unit = "month") %m-% days(1) # identifies the end date for mean calculations (the month prior to biological sampling)
# Create a Date column in your data
random_sample <- random_sample |>
  mutate(date = make_date(year, month, 1))
# Filter data for the 12 months preceding the specific date
random_sample_filtered <- random_sample |>
  filter(date >= start_date & date <= end_date)
# Calculate means
random_sample_means <- random_sample_filtered |>
  summarise(across(c(ppt, q, tmax, tmin, ws), ~ round(mean(.x, na.rm = TRUE), 3))) # means calculated from the 12 months prior to biological sampling
print(random_sample_means) # inspect calculated means
print(calc_means[random_index,]) # compare to loop calculated means

# Join env data to diptera data
Data_env <- Data |>
  left_join(calc_means, by = "sample_id")

# save output
write.csv(Data_env, "Outputs/3_diptera_taxonomic_indices_wCorine2018_TerraClimate.csv", row.names = F)

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

pred_lon <- predict_data$longitude
pred_lat <- predict_data$latitude

# Get TerraClimate data for a given coordinate(s)
# pred_terraclimate_data <- get_terraclimate(pred_lat, pred_lon, '2000-01-01', '2024-12-31', c('ppt', 'q', 'tmax', 'tmin', 'ws'), offline = TRUE, data_path = "Additional data/TerraClimate/") # links downloaded data to sets of coordinates and years provided
# saveRDS(pred_terraclimate_data, "Additional data/TerraClimate/linked_pred_terraclimate_data.RDS") # save extracted data to save time

pred_terraclimate_data <- readRDS("Additional data/TerraClimate/linked_pred_terraclimate_data.RDS") # saves time to always extract these data
str(pred_terraclimate_data$climate) # check structure of linked climate data
print(pred_terraclimate_data$climate[[1]]) # check climate data for site Rybachy_17.7.2003

# Calculate means for an individual sample_id to test if the code below works as intended
## means calculated from the 12 months preceeding biological sampling e.g., if sampling for LTR1_2010 was conducted on the 4.10.2010, then calculate mean environmental data from 01.10.2009 to 30.09.2010

# extract climate data for the first row of the diptera_filtered dataset (i.e., Rybachy_17.7.2003)
Rybachy_17.7.2003 <- as_tibble(pred_terraclimate_data$climate[[1]])
# Specific date and the preceding 12 months period
specific_date <- ymd(predict_data$date[1]) # date of biological sampling
start_date <- (specific_date %m-% months(12)) |>
  floor_date(unit = "month") # identifies the start date for mean calculations (12 months prior to biological sampling)
end_date <- (specific_date %m-% months(1)) |>
  floor_date(unit = "month") |>
  ceiling_date(unit = "month") %m-% days(1) # identifies the end date for mean calculations (the month prior to biological sampling)
# Create a Date column in your data
Rybachy_17.7.2003 <- Rybachy_17.7.2003 |>
  mutate(date = make_date(year, month, 1))
# Filter data for the 12 months preceding the specific date
Rybachy_17.7.2003_filtered <- Rybachy_17.7.2003 |>
  filter(date >= start_date & date <= end_date)
# Calculate means
Rybachy_17.7.2003_means <- Rybachy_17.7.2003_filtered |>
  summarise(across(c(ppt, q, tmax, tmin, ws), ~ round(mean(.x, na.rm = TRUE), 3))) # means calculated from the 12 months prior to biological sampling
print(Rybachy_17.7.2003_means) # inspect calculated means

# Calculate means with the function defined above
pred_calc_means <- map2_df(pred_terraclimate_data$climate, predict_data$date, process_climate_data)
# Add sample_id to the results
pred_calc_means <- bind_cols(predict_data |> dplyr::select(sample_id), pred_calc_means)

# Do the same test but on a random row to test if the code is correct
random_index <- sample(nrow(predict_data), 1)
# extract climate data for the first row of the diptera_filtered dataset
random_sample <- as_tibble(pred_terraclimate_data$climate[[random_index]])
# Specific date and the preceding 12 months period
specific_date <- ymd(predict_data$date[random_index]) # date of biological sampling
start_date <- (specific_date %m-% months(12)) |>
  floor_date(unit = "month") # identifies the start date for mean calculations (12 months prior to biological sampling)
end_date <- (specific_date %m-% months(1)) |>
  floor_date(unit = "month") |>
  ceiling_date(unit = "month") %m-% days(1) # identifies the end date for mean calculations (the month prior to biological sampling)
# Create a Date column in your data
random_sample <- random_sample |>
  mutate(date = make_date(year, month, 1))
# Filter data for the 12 months preceding the specific date
random_sample_filtered <- random_sample |>
  filter(date >= start_date & date <= end_date)
# Calculate means
random_sample_means <- random_sample_filtered |>
  summarise(across(c(ppt, q, tmax, tmin, ws), ~ round(mean(.x, na.rm = TRUE), 3))) # means calculated from the 12 months prior to biological sampling
  print(random_sample_means) # inspect calculated means
print(pred_calc_means[random_index,]) # compare to loop calculated means

# Join env data to diptera data
pred_env <- predict_data |>
  left_join(pred_calc_means, by = "sample_id")

# save output
write.csv(pred_env, "Outputs/6_prediction_data.csv", row.names = F)


###############################################################################################################
# CLEAN UP WORKSPACE
rm(list = ls())       # Remove all objects from environment
gc()                  # Frees up unused memory
p_unload(all)         # Unload all loaded packages
graphics.off()        # Close all graphical devices
cat("\014")           # Clear the console
# Clear mind :)
