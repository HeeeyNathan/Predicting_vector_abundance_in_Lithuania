# Section 1: Data description-----

#* Subsection 1.1: Data source----

#'  The data used in this paper were provided by the Lithuanian
#'  Environmental Protection Agency (LEPA). If you want to use these data
#'  for any other purpose, please contact LEPA directly. Their website address
#'  is: https://aaa.lrv.lt/lt/

#'  A subset of this data were also published in:
#'  Baker, N.J., Baker, N.J., Pilotto, F., Welti, E.A.R., Osadčaja, D. &
#'  Palinauskas, V. (2024). Recovery or reorganisation? Long-term increases
#'  in riverine taxonomic and functional diversity are confounded by compositional
#'  dynamics. Hydrobiologia. https://link.springer.com/article/10.1007/s10750-024-05665-5

#'  Sinclair, J.S., Welti, E.A.R., Altermatt, F., Álvarez-Cabria, M., Aroviita, J.,
#'  Baker, N.J., Barešová, L., Barquín, J., Bonacina, L., Bonada, N., Cañedo-Argüelles,
#'  M., Csabai, Z., de Eyto, E., Dohet, A., Dörflinger, G., Eriksen, T.E., Evtimova, V.,
#'  Feio, M.J., Ferréol, M., … Haase, P. 2024. Multi-decadal improvements in the ecological
#'  quality of European rivers are not consistently reflected in biodiversity metrics.
#'  Nature Ecology & Evolution, 8, 430-441. https://doi.org/10.1038/s41559-023-02305-4

#'  Cano-Barbacil, C., Sinclair, J.S., Welti, E.A.R., & Haase, P. 2025. Recovery and
#'  Degradation Drive Changes in the Dispersal Capacity of Stream Macroinvertebrate
#'  Communities. Global Change Biology, 31:e70054. https://doi.org/10.1111/gcb.70054


#* Subsection 1.2: Funding source----

#'  The collation and processing of this dataset was supported by the Diana Osadcaja and
#'  the Lithuanian Environmental Protection Agency (https://aaa.lrv.lt/) who collected
#'  the data. The project was funded by the Lithuanian Research Council (project number
#'  S-PD-22-72).


#* Subsection 1.3: Background information----

#' This dataset contains count data of 962 macroinvertebrate taxa (taxa is used due to
#' mixed identification levels) collected from freshwater bodies across lithuania.
#' The initial biomonitoring dataset was created as part of Lithuania's national monitoring
#' scheme, in line with the European Union's Water Framwork Directive. The dataset spans
#' a period from 2010 to 2023, covering both lentic (rivers) and lotic (lakes) freshwater
#' bodies. The dataset includes count data of macroinvertebrates from 953 sampling sites; 615
#' riverine sites and 338 lake sites. LEPA employs a staggered sampling approach, whereby
#' most sites are resampled in 3-year increments, with only a few sites (~20%) being sampled
#' more consistently. Rivers were sampled from 2010 to 2022. Lakes were sampled from 2013 to
#' 2023.
#'
#' For river sites (n = 615), the average number of sites sampled per year is 151, with
#' distances between sites averaging 132 km and remaining relatively consistent through
#' time. On average, sites were resampled 3 times within the 14-year observation period
#' (2010 - 2023), with the majority of sites being sampled twice (mode). Around 85% (521)
#' of sites had time series shorter than 5 years in length, while 15% (94 sites) had time
#' series greater than or equal to 5 years. The dataset contains count data of 761
#' unique taxa, with the total abundance of these taxa being 1 156 915 individuals. The
#' identification of these taxa was mostly at the species (61%) and genus (29%) levels,
#' with the remaining 10% of taxa being identified to the family (9%) or sub-family-levels
#' (1%).
#'
#' For lake sites (n = 338), the average number of sites sampled per year was 60, with the
#' distances between sites averaging 127 km and remianing relatively consistent through time.
#' On average, sites were resampled twice within the 11-year observation period (2013 - 2023),
#' with most sites being resampled twice (mode). Around 86% (290 sites) of sites had timeseries
#' shorter than 3 years in length, while 13% (38 sites) had time series longer than or equal
#' to 3 years. The dataset contains count data of 675 unique taxa, with the total abundance
#' of these taxa being 320 837 individuals. The identification of these taxa was mostly at
#' the species (66%) and genus (27%) levels, with the remaining 8% of taxa being identified
#' to the family (7%) or sub-family-levels (1%). There majority of lakes are located in the
#' north-east part of the country.
#'
#' In Baker et al. (2024), a subset of this dataset was used to determine temporal changes
#' in macroinvertebrate communities between 2010 and 2022. This investigation is a regional
#' follow-up to Haase et al. (2023) - Nature, and Sinclair et al. (2024) - Nat. Ecol. Evol.
#' which also included parts of the larger Lithuanian dataset.
#'
#' Lithuania is in the boreal ecoregion, with mean annual precipitation and temperature
#' (from 1991 to 2020) being 679 mm and 7.38 °C, respectively. Situated at the edge of
#' the East European Plain, its low-lying post-glaciated landscape has a maximum elevation
#' of 297 m a.s.l., with 758 streams greater than 10 km long transecting the country. The
#' largest and most economically important river is the Nemunas, which originates in Belarus,
#' forms the catchment area for 72% of Lithuania’s 64,800 km2 territory, and drains into
#' the Baltic Sea. Other noteworthy catchments include the basins of the Venta, Lielupė,
#' and Dauguva rivers. All Lithuanian rivers are calcareous and lowland (< 200 m a.s.l.),
#' with their river typology being based on catchment size and slope. Despite heavy
#' deforestation during the soviet regime, the coverage of postsoviet forests in Lithuania
#' comprised mostly of pine, spruce, and birch—has recovered from 31.3% in 1992 to 35.5% in
#' 2023. Lithuania’s built-up areas (3.73% coverage) have increased since 1971, while the
#' coverage of arable land (45.99%), meadows/pastures (5.55%), wetlands (5.62%), and other
#' land uses (5.78%) has declined. Lithuanian agricultural land has gradually shifted from
#' less extensive to more intensive agricultutral practices, with high nitrogen levels being
#' attributed to these modern agricultural practices, as well as the legacy of soviet-era
#' fertilization practices.
#'
#' Most of Lithuania’s arable land is restricted to a central strip that runs north to south,
#' from above the city of Šiauliai to below the city of Kaunas. This arable region transects
#' the country into two distinct climatic regions across an east-west gradient. The east of the
#' country is categorized by higher average temperatures, higher precipitation, more sunshine
#' hours and higher wind speeds, while the east is cooler, drier, and calmer. These regions can
#' also be delimited by lowlands in the west and highlands in the east. More information is
#' available here: https://www.meteo.lt/en/climate/lithuanian-climate/air-temperature/

#' Lithuania is geographically located along an important migratory pathway of birds migrating
#' northward in the spring and southward in the fall. The worlds second oldest ornithological
#' station, Ventas Ragas, is located in the western part of the country, on the Curonian lagoon.
#' Here, and at a nearby ornithological station called Ribachy located on the Curonian Spit,
#' birds are screened for avian malarial parasites of the Haemosporidian genera Leucocytozoon,
#' Haemoproteus, and Plasmodium. These parasites are transmitted by insect vectors, namely members
#' of the Simulidae (black flies), the Ceratopogoniidae (biting midges), the Culicidae (mosquitos),
#' during feeding on blood and have been increasing in prevalence among birds since observations
#' began (see below).

print(readRDS("Plots/Figure1_parasite_prevalence_dynamics.RDS"))


#' Evidence suggests that Diptera (i.e., true flies) are increasing in both abundance and richness
#' throughout Lithuania (Baker et al., 2024), potentially explaining the increase in prevalence of
#' malarial haemosporidian parasite prevalence in birds.

print(readRDS("Plots/Baker_et.al._2024_trends.rds"))

#' However, the driving factors behind these increases and what factors control their presence and
#' abundance are still not well understood. Furthermore, initial data exploration appears to show
#' that vector presence and abundance is highest in the eastern Lithuania, which is comparatively
#' undersampled, both in terms of bird malaria and vector dynamics, compared to western Lithuania
#' near the curonian lagoon.

#' The aim of this study is to determine the driving forces behind vector presence and abundance, in
#' an attempt to better understand vector and malarial dynamics, while simultaneously identifying
#' underrepresented and under sampled regions that may yeild better more sucess in finding bird
#' malarial parasites, thereby shedding light on the processes shaping malarial spread amongst birds
#' and thus potentially informing human-malaria relationships.

#' The species investigated in the paper are:
#'   -Freshwater macroinvertebrates, which are categorized as:
#'    insects in their nymph and larval stages, snails, worms,
#'    crayfish, and clams that spend at least part of their
#'    lives in water and large enough to see without the aid
#'    of a microscope.
#'   -The initial data set contained counts for all the macro-
#'    invertebrate species collected  at a site, however, only
#'    dipteran data were retained.

#' Below are some points describing the data collection.
#'  1. River and lake sites form part of Lithuania's national biomonitoring
#'     scheme

#'  2. At each site, macroinvertebrate sampling included both quantitative
#'     (standard multihabitat kick-sampling) and qualitative (hand-picking
#'     organisms from underwater objects) components (Šidagytė-Copilas &
#'     Arbačiauskas, 2022).

#'  3. The quantitative samples were collected using the multihabitat method.
#'     Within a site, 10 subsamples (later to be pooled into one sample) were
#'     collected from all available microhabitats proportionally to their
#'     estimated distribution. Each subsample was collected using the kick
#'     method by disturbing the 40 cm length area in front of the standard
#'     handnet (25 × 25 cm opening, 0.5 mm mesh size), resulting in 0.1 m2 of
#'     the bottom area sampled per subsample, and 1 m2 per the whole sample.

#'  4. Additional qualitative samples were collected by searching for
#'     macroinvertebrates attached to underwater objects (roots, stones, plants,
#'     etc.). In cases where qualitative sampling led to the inclusion of a
#'     taxon not identified during the quantitative sampling, taxa were
#'     included with an assigned abundance of one.

#'  5. Most sites are resampled in 3-year increments, with few sites having
#'     more regular sampling.

#' Below are some points describing the initial data processing collection.
#'  1. To ensure our biological dataset conformed to contemporary scientific
#'     nomenclature (Grenié et al., 2021, 2022), taxon names were individually
#'     checked by Lithuanian macroinvertebrate experts and verified against the
#'     2023 Global Biodiversity Information Facility (GBIF) taxonomic backbone.
#'     All non-freshwater macroinvertebrates were removed, including
#'     microinvertebrates (e.g. hydrachnidia, copepoda) and non-freshwater
#'     invertebrates (e.g. terrestrial, semi-terrestrial taxa).

#'  2. For each site-year, the following taxonomic diversity indices were
#'     calculated: total diptera abundance, diptera richness, vector abundance,
#'     and vector richness. The vegan R package (Oksanen et al., 2022) was used
#'     to calculate all other taxonomic indices.

#'  3. For each site and year, an ecological quality ratio was calculated according
#'     to the methods outlined in (Šidagytė-Copilas & Arbačiauskas, 2022). The EQR
#'     is in accordance with the EU's water framework directive and describes the
#'     ecological quality of a site compared to a comparable reference site. The
#'     EQR is bound between zero and one. From the EQR, an ecological quality class
#'     (EQC) is assisned. This metric is the primary covariate in our analysis.

#'  4. Environmental data linked to the study were accessible for a subset of sites.
#'     The environmental dataset comprised monthly measurements of diverse variables,
#'     encompassing climatic parameters including water temperature (°C) and discharge
#'     (m3 s−1), and water quality parameters such as alkalinity (mmol l−1), dissolved
#'     oxygen (mg l−1), pH, electrical conductivity (μS cm−1), ammonium (mg l−1),
#'     nitrate (mg l− 1), nitrite (mg l−1), mineralised nitrogen (mg l−1), total
#'     nitrogen (mg l−1), ortho-phosphate (mg l−1), total phosphorus (mg l−1), and
#'     suspended solids (mg l−1). Additionally, worldclim data was used to extract
#'     climatic variables from each site-year combination. The following parameters
#'     were extracted for each site and year: precipitation (ppt), discharge (q),
#'     maximum temperature (tmax), minimum temperature (tmin), and wind speed (ws).
#'
#'     For each biological site-year, the specific month of the biological sampling
#'     was determined. Subsequently, mean annual values of the predictor variables
#'     were computed using data from the preceding 12 months.
#'
#'  5. Elevation data for each unique site-year combination was extracted using
#'     Google's Maps Elevation API.
#'
#'  6. Landuse data for each unique site-year combination was extracted using
#'     the 2018 Corine landcover dataset. The coarsest level of naming was used for
#'     the landcover (LABEL 1) to avoid too many collinear covaraites. We used 1000m
#'     buffer zones around the sites to account for relatively small-scale effects
#'     of landuse on dipteran communities.

#' The basic statistical unit of the analysis is the sum of abundances of vectors
#' species counted within sample from a biomonitoring survey within a given year.

#' Since many site-year observations had vector abundances of zero, the data are
#' likely zero inflated with a large range (more on this later).


#* Subsection 1.4: The variables----
#' Response variable in this exercise:
#'   Counts:
#'    The number of vector individuals (black flies, biting midges, and mosquitos)
#'    at a sampling site in a specific year


#' Covariates:
#'   -EQR (0 - 1)
#'   -precipitation (ppt)
#'   -discharge (q)
#'   -maximum temperature (tmax)
#'   -minimum temperature (tmin)
#'   -wind speed (ws)
#'   -elevation (m.a.s.l.)
#'   -Agricultural areas (%)
#'   -Artificial.surfaces (%)
#'   -Forest.and.semi.natural.areas (%)
#'   -Water.bodies (%)
#'   -Wetlands (%)
#'   -Year
#'   -Day of year (doy)
#'   -Waterbody_type (Lake vs River)
#'   -Latitude (y coordinate of sampling site)
#'   -Longitude (x coordinate of sampling site)


#' Dependency:
#'  - Spatial correlation.
#'  - Temporal correlation.
#'  - Spatial-Temporal correlation.


#* Subsection 1.5: Underlying questions----

#' Aim of the analysis:
#'  - Determine factors influencing vector abundance (and presence).
#'  - Determine areas in which sampling should be conducted to maximise chances of catching malarial vectors.



# Section 2: Import the data----
#* Subsection 2.1: Load the packages----

#' Load all packages and our support file.
library(pacman)
library(lattice)
library(ggplot2)
library(mgcv)
library(sf)
library(gratia)
library(plyr)
library(tidyr)
library(DHARMa)
library(rgeoboundaries)
library(gstat)
library(INLA)
library(fmesher)
library(MASS)
library(animation)
library(scales)
library(httr)
library(googleway)
library(tidyverse)
library(janitor)
library(GGally)
library(ggmap)
library(cowplot)
library(reshape)
library(grid)
library(gridExtra)
library(kableExtra)
source(file = "Additional functions/HighstatLibV15.R") # <---- the use of these functions requires a citation: Zuur, A.F., Ieno, E.N., Walker, N., Saveliev, A.A., Smith, G.M., 2009. Mixed Effects Models and Extensions in Ecology with R. Springer, New York https://doi.org/10.1007/978-0-387-87458-6.


#* Subsection 2.2: Load the data----

#' Load the dipteran vector data.
# Import the data
df <- read.csv("Outputs/4_diptera_taxonomic_indices_wCorine2018_TerraClimate_elevation.csv", h = T, sep = ",", stringsAsFactors = FALSE, check.names = FALSE) |>
  clean_names() |>
  mutate(dipt_abund_pa    = ifelse(dipt_abund == 0, 0, 1), # create vector with presence-absence data
         dipt_abund_pos   = ifelse(dipt_abund > 0, dipt_abund, NA), # create vector with only positive data (i.e., not zero)
         vec_abund_pa     = ifelse(vec_abund == 0, 0, 1), # create vector with presence-absence data
         vec_abund_pos    = ifelse(vec_abund > 0, vec_abund, NA), # create vector with only positive data (i.e., not zero)
         fyear            = factor(year, levels = c(2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023)), # make year a factor
         cyear            = year - median(year), # centered year - helps model convergence to center variables for the model
         iyear            = year - min(year) + 1, # year as an index starting from 1
         year             = as.numeric(year), # year as numeric
         month            = factor(month, levels = c(4, 5, 9, 10, 11),
                                          labels = c("April", "May", "September", "October", "November"),
                                          ordered = F), # make month a factor
         state            = factor(state, levels = c("A", "HM", "N"),
                                          labels = c("Artificial", "Heavily_modified", "Natural"),
                                          ordered = F), # make state a factor
         eqc              = factor(eqc, levels = c("Bad", "Poor", "Moderate", "Good", "High"),
                                        labels = c("Bad", "Poor", "Moderate", "Good", "High"),
                                        ordered = F), # make EQC a factor
         waterbody_name   = factor(waterbody_name), # make waterbody_name a factor
         waterbody_type   = factor(waterbody_type, levels = c("lake", "river"),
                                                   labels = c("Lake", "River"),
                                                   ordered = T), # make waterbody_type a factor
         date             = as.Date(date, format = "%Y-%m-%d"), # make the dates dates
         doy              = yday(date)) |> # calculate sampling day of year (doy)
  dplyr::rename(agriculture      = agricultural_areas,
         artificial       = artificial_surfaces,
         natural          = forest_and_semi_natural_areas,
         water            = water_bodies,
         wetlands         = wetlands) |>
  dplyr::select(-c(observation_period, sampling_events)) |> # removes unnecessary columns
  arrange(desc(waterbody_type), site_id, year) |>  # order the data.frame
  as.data.frame() # Convert tibble to dataframe because some older code does not recognize tibble

#' What do we have?
names(df)
glimpse(df)
head(df)


# Section 3: Data coding----

#' We will relabel some of the names to make them more intuitive.

#' This is the response variable.
df$Counts <- df$vec_abund

#' And the covariates:
df$Lat    <- df$latitude
df$Long   <- df$longitude
df$Year   <- df$year

#' For some graphs, we need year as a categorical variable.
df$fyear


#' Convert longitude and latitude to UTM in zone 20. The UTM_Transform
#' function is on our support file.
XY.utm <- UTM_Transform(x = df$Long,
                        y = df$Lat,
                        zone = 34,
                        Hemisphere = "north")

df$Xkm <- XY.utm[,"X"] / 1000
df$Ykm <- XY.utm[,"Y"] / 1000

#' We define some common settings for INLA.
MyControlCompute  <- list(config = TRUE,    #' Allow for posterior simulation
                          dic = TRUE,       #' Calculate AIC
                          waic = TRUE,      #' Calculate WAIC
                          residuals = TRUE) #' Get residuals (see below)
MyControlPredictor  <- list(compute = TRUE, #' Calculate fitted values
                            link = 1)       #' Predict on the scale of
                                            #' the response variable.

# Define theme for plotting via ggplot
My_theme <- theme(panel.background = element_blank(),
                  panel.border = element_rect(fill = NA, linewidth = 1.25),
                  strip.background = element_rect(fill = "white",
                                                  color = "white", linewidth = 1.25),
                  # legend.position = "right",
                  # axis.text.x = element_text(angle = 45, hjust = 1),
                  text = element_text(size = 16))

# Section 4: House keeping----

#* Subsection 4.1: Number of observations by site and year----

#'  How many observations do we have per year?
#'  Calculate how many observations we have per year?
table(df$year)

#'  Likely need to drop 2023 because it has far fewer sites than the other years.
#'  Also, 2010, 2011, and 2012 are only available for rivers, so might as well drop
#'  those too.

df <- df |>
  filter(year >= 2013 & year <= 2022) # Remove years less than 2013 or greater than 2022

#'  How many observations do we have per location per year?
#'  Calculate the observation period and number of repeated samplings for each site
df <- df |>
  group_by(site_id) |>
  mutate(
    observation_period = max(year) - min(year) + 1,
    sampling_events = n_distinct(year)
  ) |>
  ungroup() |>
  arrange(desc(waterbody_type), site_id, year)

df[, c("site_id", "observation_period", "sampling_events")] |> unique()

#'  After removing the years 2010, 2011, 2012, and 2023, some sites were not sampled in
#'  the remaining years, leaving them with a value of zero. Further, many sites were only
#'  sampled once. These site should be remove since they are not technically time series.
#'  Lets keep these for now

#'  Keep sites with less than 1 sampling events throughout the observation period
df <- df |>
  filter(sampling_events >= 1) |>
  arrange(desc(waterbody_type), site_id, year)

#* Subsection 4.2: Missing values----

#'  Missing values?
colSums(is.na(df))
100 * colSums(is.na(df)) / nrow(df)

#'  There are many NA values, including some ~1% in our focal covariate: EQR. Let's remove
#'  the rows with missing EQR values and then remove columns containing many NAs which are
#'  unnecessary for this analysis.

#'  Remove rows where EQR is NA
df <- df |>
  filter(!is.na(eqr))

#'  Remove columns containing NA values
df <- df |>
  dplyr::select(where(~ all(!is.na(.))))

#* Subsection 4.3: Check values after cleaning----

#'  How many sites do we have in total?
NROW(unique(df$site_id))

#'  We still have 934 sites

#'  How many lake sites vs rivers sites do we have?
df |>
  group_by(waterbody_type) |>
  summarise(unique_sites = n_distinct(site_id)) |>
  arrange(desc(unique_sites))

#'  We have more than double the amount of lake sites compared to river sites

#'  How sites were sampled per year?
table(df$year)

# That is better, but 2018 and 2019 still have fewer sites sampled. Lets keep them in for now

#'  Average number of sites sampled per year (overall)
df |>
  group_by(year) |>
  summarise(num_sites = n_distinct(site_id)) |>
  summarise(avg_sites = mean(num_sites))

# 209 sites sampled on average per year

#'  Average number of river sites sampled per year
df |>
  filter(waterbody_type == "River") |>
  group_by(year) |>
  summarise(num_river_sites = n_distinct(site_id)) |>
  summarise(avg_river_sites = mean(num_river_sites))

# 152 river sites sampled on average per year

#'  Average number of lake sites sampled per year
df |>
  filter(waterbody_type == "Lake") |>
  group_by(year) |>
  summarise(num_river_sites = n_distinct(site_id)) |>
  summarise(avg_river_sites = mean(num_river_sites))

# 58 river sites sampled on average per year (rounded up from)

#'  Calculate average distance between sites by year and waterbody

# Function to calculate mean pairwise distance within a group of sites
calc_mean_distance <- function(x, y) {
  if (length(x) <= 1) {
    return(NA)  # Need at least 2 sites to calculate distance
  }

  site_pairs <- combn(length(x), 2)
  distances <- numeric(ncol(site_pairs))

  for (i in 1:ncol(site_pairs)) {
    idx1 <- site_pairs[1, i]
    idx2 <- site_pairs[2, i]
    distances[i] <- sqrt((x[idx1] - x[idx2])^2 + (y[idx1] - y[idx2])^2)
  }

  return(mean(distances))
}

#' Average distance between sites per year
df |>
  dplyr::select(year, site_id, Xkm, Ykm) |>
  distinct() |>
  group_by(year) |>
  summarise(
    num_sites = n_distinct(site_id),
    avg_distance = calc_mean_distance(Xkm, Ykm)) |>
  summarise(avg_dist_sites = mean(avg_distance))

# Average of 138 km between sites in each year

#' Average distance between sites per year by waterbody type
df |>
  dplyr::select(year, site_id, waterbody_type, Xkm, Ykm) |>
  distinct() |>
  group_by(year, waterbody_type) |>
  dplyr::summarise(
    num_sites = n_distinct(site_id),
    avg_distance = calc_mean_distance(Xkm, Ykm)) |>
  ungroup() |>
  group_by(waterbody_type) |>
  dplyr::summarise(
    avg_site_dist_wb = mean(avg_distance, na.rm = TRUE)
  )

# Average of 127 km between lake sites in each year
# Average of 133 km between river sites in each year

#' Get a dataframe containing the unique sampling sites used in the analysis (not duplicated)
unique_sites <- df |>
  mutate(eqc_numeric = as.numeric(eqc)) |>
  dplyr::select(site_id:eqc, eqc_numeric, everything()) |>
  group_by(site_id) |>
  dplyr::summarize(
    latitude = first(latitude),
    longitude = first(longitude),
    Lat = first(latitude),
    Long = first(longitude),
    Xkm = first(Xkm),
    Ykm = first(Ykm),
    waterbody_type = first(waterbody_type),
    avg_eqr = mean(eqr, na.rm = TRUE),
    mode_eqc = get_mode(eqc),
    avg_eqc_num = round(mean(eqc_numeric, na.rm = TRUE), 0),
    observation_count = n()
  ) |>
  ungroup() |>
  mutate(avg_eqc_cat = factor(avg_eqc_num,
                              levels = c(1, 2, 3, 4, 5),
                              labels = c("Bad", "Poor", "Moderate", "Good", "High"))) |>
  dplyr::arrange(desc(waterbody_type), site_id)

write_csv(unique_sites, "Outputs/5_unique_sites_for_plotting.csv")

# For each waterbody type in the unique_sites dataset, calculate the percentage of sites with observation periods less than 5 years, and greater than or equal to 5 years

unique_sites |>
  group_by(waterbody_type) |>
  dplyr::summarise(
    less_than_5_years = sum(observation_count < 5),
    greater_than_or_equal_to_5_years = sum(observation_count >= 5),
    total_sites = n()
  ) |>
  mutate(
    less_than_5_years_percent = round((less_than_5_years / total_sites) * 100, 2),
    greater_than_or_equal_to_5_years_percent = round((greater_than_or_equal_to_5_years / total_sites) * 100, 0)
  )

# for each waterbody type, calculate summaries of the observation period (e.g., mean, min, max, mode, range, etc.

unique_sites |>
  group_by(waterbody_type) |>
  dplyr::summarise(
    mean_observation_period = mean(observation_count, na.rm = TRUE),
    min_observation_period = min(observation_count, na.rm = TRUE),
    max_observation_period = max(observation_count, na.rm = TRUE),
    mode_observation_period = get_mode(observation_count),
    range_observation_period = max(observation_count, na.rm = TRUE) - min(observation_count, na.rm = TRUE)
  )

# Section 5: Data exploration----

# Start by defining a preferred figure format, called 'My_theme'
My_theme <- theme(panel.background = element_blank(),
                  panel.border = element_rect(fill = NA, linewidth = 1),
                  strip.background = element_rect(fill = "white",
                                                  color = "white", linewidth = 1),
                  text = element_text(size = 14),
                  panel.grid.major = element_line(colour = "white", linewidth = 0.1),
                  panel.grid.minor = element_line(colour = "white", linewidth = 0.1),
                  legend.position = "none")

#* Subsection 5.1: Outliers----

#' Make a Cleveland dotplot of the response variable 'vector richness' and the
#' continuous covariates.
MyVar <- c("Counts", "eqr", "ppt", "tmax", "tmin", "q", "ws", "elevation",
           "agriculture", "artificial", "natural",
           "year")
Mydotplot(df[,MyVar])
#' All ok.


#* Subsection 5.2: Collinearity----

df |>
  ggpairs(columns = MyVar,
          aes(alpha = 0.8), lower = list(continuous = "smooth_loess",
          combo = wrap("facethist", binwidth = 5))) + My_theme

df |>
  dplyr::select(all_of(MyVar)) |>
  corvif()

# Perhaps tmin and q can cause some trouble. Let's remove them and see

MyVar_red <- c("eqr", "ppt", "tmin", "ws", "elevation",
               "agriculture", "artificial", "natural",
               "year")

df |>
  ggpairs(columns = MyVar_red,
          aes(alpha = 0.8), lower = list(continuous = "smooth_loess",
          combo = wrap("facethist", binwidth = 5))) + My_theme

df |>
  dplyr::select(all_of(MyVar_red)) |>
  corvif()

# seems okay now, but maybe we can drop some terms later?


#* Subsection 5.3: Relationships----

# Then plot
p1 <- df |>
  ggplot(aes(y = Counts, x = eqr)) +
  geom_smooth(method = "gam") +
  labs(y = "Vector abundance",
       x = "Ecological quality") +
  geom_jitter(shape = 19, size = 3.5, height = 0.5,
              width = 0.5, alpha = 0.5) +
  # xlim(round(range(df$eqr), 0))  + ylim(range(df$Counts)) +
  xlim(round(range(df$eqr), 0))  + ylim(range(df$Counts[df$Counts <= 300])) +
  My_theme
p1
# a slightly positive relationship
p1 + facet_grid(~waterbody_type)
# a similar relationship between lakes and rivers


p2 <- df |>
  ggplot(aes(y = Counts, x = ppt)) +
  geom_smooth(method = "gam") +
  labs(y = "Vector counts",
       x = "Precipitation") +
  geom_jitter(shape = 19, size = 3.5, height = 0.5,
              width = 0.5, alpha = 0.5) +
  # xlim(round(range(df$ppt), 0))  + ylim(range(df$Counts)) +
  xlim(round(range(df$ppt), 0))  + ylim(range(df$Counts[df$Counts <= 300])) +
  My_theme
p2
# not much of a relationship
p2 + facet_grid(~waterbody_type)
# maybe some difference between lakes and rivers... interaction?

p3 <- df |>
  ggplot(aes(y = Counts, x = tmin)) +
  geom_smooth(method = "gam") +
  labs(y = "Vector counts",
       x = "Minimum temperature") +
  geom_jitter(shape = 19, size = 3.5, height = 0.5,
              width = 0.5, alpha = 0.5) +
  # xlim(round(range(df$tmin), 0))  + ylim(range(df$Counts)) +
  xlim(round(range(df$tmin), 0))  + ylim(range(df$Counts[df$Counts <= 300])) +
  My_theme
p3
# maybe some slightly negative relationship
p3 + facet_grid(~waterbody_type)
# differences between lakes and rivers, maybe some non-linearity there?

p4 <- df |>
  ggplot(aes(y = vec_abund, x = ws)) +
  geom_smooth(method = "gam") +
  labs(y = "Vector counts",
       x = "Wind speed") +
  geom_jitter(shape = 19, size = 3.5, height = 0.5,
              width = 0.5, alpha = 0.5) +
  # xlim(round(range(df$ws), 0))  + ylim(range(df$Counts)) +
  xlim(round(range(df$ws), 0))  + ylim(range(df$Counts[df$Counts <= 300])) +
  My_theme
p4
# not much of a relationship
p4 + facet_grid(~waterbody_type)
# not much there.

p5 <- df |>
  ggplot(aes(y = Counts, x = elevation)) +
  geom_smooth(method = "gam") +
  labs(y = "Vector counts",
       x = "Elevation") +
  geom_jitter(shape = 19, size = 3.5, height = 0.5,
              width = 0.5, alpha = 0.5) +
  # xlim(round(range(df$elevation), 0))  + ylim(range(df$Counts)) +
  xlim(round(range(df$elevation), 0))  + ylim(range(df$Counts[df$Counts <= 300])) +
  My_theme
p5
# a postive relationship
p5 + facet_grid(~waterbody_type)
# not much difference between lakes and rivers

p6 <- df |>
  ggplot(aes(y = Counts, x = agriculture)) +
  geom_smooth(method = "gam") +
  labs(y = "Vector counts",
       x = "Agricultural areas (%)") +
  geom_jitter(shape = 19, size = 3.5, height = 0.5,
              width = 0.5, alpha = 0.5) +
  # xlim(round(range(df$agriculture), 0))  + ylim(range(df$Counts)) +
  xlim(round(range(df$agriculture), 0))  + ylim(range(df$Counts[df$Counts <= 300])) +
  My_theme
p6
# not much there
p6 + facet_grid(~waterbody_type)
# maybe a slight negative relationship for rivers, but not clear


p7 <- df |>
  ggplot(aes(y = Counts, x = artificial)) +
  geom_smooth(method = "gam") +
  labs(y = "Vector counts",
       x = "Artificial surfaces (%)") +
  geom_jitter(shape = 19, size = 3.5, height = 0.5,
              width = 0.5, alpha = 0.5) +
  # xlim(round(range(df$artificial), 0))  + ylim(range(df$Counts)) +
  xlim(round(range(df$artificial), 0))  + ylim(range(df$Counts[df$Counts <= 300])) +
  My_theme
p7
# not much there
p7 + facet_grid(~waterbody_type)
# not much there

p8 <- df |>
  ggplot(aes(y = Counts, x = natural)) +
  geom_smooth(method = "gam") +
  labs(y = "Vector counts",
       x = "Forest & semi natural areas (%)") +
  geom_jitter(shape = 19, size = 3.5, height = 0.5,
              width = 0.5, alpha = 0.5) +
  # xlim(round(range(df$natural), 0))  + ylim(range(df$Counts)) +
  xlim(round(range(df$natural), 0))  + ylim(range(df$Counts[df$Counts <= 300])) +
  My_theme
p8
# not much there
p8 + facet_grid(~waterbody_type)
# maybe a slight postive relationship for rivers, but not much there


p9 <- df |>
  ggplot(aes(y = Counts, x = year)) +
  geom_smooth(method = "gam") +
  labs(y = "Vector counts",
       x = "Time") +
  geom_jitter(shape = 19, size = 3.5, height = 0.5,
              width = 0.5, alpha = 0.5) +
  scale_x_continuous(breaks = 2013:2022, limits = c(2013, 2022)) +
  # ylim(range(df$Counts)) +
  ylim(range(df$Counts[df$Counts <= 300])) +
  My_theme
p9
# not much there
p9 + facet_grid(~waterbody_type)
# some non-linearity in rivers perhaps?

#' What about using using year as a factor? A factor allows for
#' more sudden changes as compared to a smoother.
p10 <- df |>
  ggplot(aes(y = Counts, x = fyear)) +
  geom_boxplot(size = 0.5) +
  labs(y = "Vector counts",
       x = "Year") +
  ylim(range(df$Counts[df$Counts <= 300])) +
  My_theme
p10
# not much there
p10 + facet_grid(~waterbody_type)
# no clear patterns


p11 <- df |>
  ggplot(aes(y = Counts, x = waterbody_type)) +
  geom_boxplot(size = 0.5) +
  labs(y = "Vector counts",
       x = "Waterbody_type") +
  ylim(range(df$Counts[df$Counts <= 300])) +
  My_theme
p11
# the range for rivers is much larger

#' What about plotting the time series for each sampling site?
p12 <- df |>
  ggplot(aes(x = year, y = Counts,
           group = site_id)) +
  scale_x_continuous(breaks = 2013:2022, limits = c(2013, 2022)) +
  xlab("Year") + ylab("Counts") +
  geom_line()
p12
table(df$site_id, df$year)


#* Subsection 5.4: Normality and zero inflation----

# Frequency plot
p1 <- df |>
  ggplot(aes(Counts)) +
  geom_freqpoly(bins = 15) +
  labs(x = "Vector counts",
       y = "Frequency") +
  My_theme
p1
# High number of zeros and count positively skewed

# How many zeros do we have?
df |>
  dplyr::summarise(percentage_zero = sum(Counts == 0) / n() * 100)
# 64% zeros - too many!!!


#* Subsection 5.5: Dependency----

#' Get shape files
surroundings.shp <- geoboundaries(c("Latvia", "Estonia", "Poland", "Belarus", "Russia"))
surroundings <- st_transform(surroundings.shp, crs = 4326)
surroundings <- st_make_valid(surroundings)

lithuania.shp <- geoboundaries("Lithuania")
lithuania <- st_transform(lithuania.shp, crs = 4326)

#' Plot Lithuania using ggplot, and limit the plotting to a certain xlim and ylim.
ggplot(data = lithuania) +
    geom_sf(fill = "transparent") +
    theme_minimal() +
    xlim(-21, -27) + ylim(53, 57)

ggplot(data = surroundings) +
    geom_sf(fill = "transparent") +
    theme_minimal() +
    xlim(-21, -27) + ylim(53, 57)

#' Once we present the results of the spatial models, we would
#' like to remove that part of the the spatial term that covers the baltic sea.
#' First we want to get the study area as we see it, without having to add
#' to xlim and ylim commands. To do this, define a bounding box and
#' crop the bounding box and the Lithuania map.

#' Define bounding box (xmin, xmax, ymin, ymax).
bb <- st_bbox(c(xmin = 20.9, xmax = 26.9, ymin = 53.85, ymax = 56.5), crs = st_crs(lithuania))

#' Crop Lithuania to the bounding box.
CroppedLithuania <- st_crop(lithuania, bb)

#' This is what we now have.
MyCex <- 3 * sqrt(df$Counts + 1) / 10
p <- ggplot(data = CroppedLithuania) +
       geom_sf(fill = "white") +
       theme_minimal() +
       # labs(title = "Study area") +
       geom_point(aes(x = Long,
                      y = Lat,
                      color = waterbody_type,
                      shape = waterbody_type),
                      size = MyCex,
                      data = df) +
       scale_color_manual(name = "Waterbody type",
                          values = c("black", "red"),
                          guide = guide_legend(override.aes = list(size = 4))) +
       scale_shape_manual(name = "Waterbody type", values = c(19, 18)) +
       xlab("Longitude") + ylab("Latitude") +
       My_theme + theme(legend.position = "bottom")
p

# facet_wrap the plot by year
p + facet_wrap(~fyear) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

CroppedSurroundings <- st_crop(surroundings, bb)

#'  Define bounding box as a polygon.
bb_coords <- matrix(
  c(20.9, 53.85,
    26.9, 53.85,
    26.9, 56.5,
    20.9, 56.5,
    20.9, 53.85), # Closing the loop
  byrow = TRUE, ncol = 2)

bb_coords

#' Create a polygon from the bb_coords coordinates.
bb_polygon <- st_polygon(list(bb_coords))

#' Convert this polygon to an sf object. Give it the same crs as CroppedLithuania
bb_sf <- st_sfc(bb_polygon, crs = st_crs(CroppedLithuania))

#' Subtract the cropped land area (CroppedLithuania) from the bounding box to create the "surroundings" polygon
SeaPolygon <- st_difference(bb_sf, st_geometry(CroppedSurroundings))

#' Convert CroppedLithuania, Surorunds, and SeaPolygon to UTM coordinates. Transform to UTM zone 34
CroppedLithuania_UTM <- st_transform(x = CroppedLithuania,
                                     crs = "+proj=utm +zone=34 +north +ellps=WGS84 +units=km +datum=WGS84")

SeaPolygon_UTM <- st_transform(x = SeaPolygon,
                               crs = "+proj=utm +zone=34 +north +ellps=WGS84 +units=km +datum=WGS84")

CroppedSurroundings_UTM <- st_transform(x = CroppedSurroundings,
                                        crs = "+proj=utm +zone=34 +north +ellps=WGS84 +units=km +datum=WGS84")

#' Plot the result (in UTM).
ggplot() +
  geom_sf(data = SeaPolygon_UTM,
          fill = "lightblue",
          color = "black") + # Land
  geom_sf(data = st_sf(CroppedSurroundings_UTM),
          fill = "grey") + # surrounds
  geom_sf(data = CroppedLithuania_UTM,
          fill = "white",
          color = "black") + # Land
  theme_minimal() +
  labs(title = "Study area") +
  geom_point(data = df,
             aes(x = Xkm, y = Ykm),
             size = 0.5,
             alpha = 0.5)

#* Subsection 5.6: Conclusions----

#' 1. We have count data with plenty of zeros. The range of the data is large.
#'    I anticipate a negative binomial distribution or worst case a ZIP or ZINB
#'    (unless I can explain them with a time component or spatial component).
#' 2. We have spatial-temporal data.
#' 3. There appears to be a positive relationship between EQR, precipiation, minimum
#'    temperature and elevation, and there seems to be some variation with waterbody
#'    type. Interaction term?
#' 4. Seasonality should be captured by waterbody type: Rivers always sampled in autumn,
#'    lakes always sampled in spring.



# Section 6: Data standardization ----

#* Subsection 6.1: Standardize covariates----

df <- df |>
  mutate(
    eqr.std         = MyStd(eqr),
    ppt.std         = MyStd(ppt),
    q.std           = MyStd(q),
    tmin.std        = MyStd(tmin),
    tmax.std        = MyStd(tmax),
    ws.std          = MyStd(ws),
    elevation.std   = MyStd(elevation),
    agriculture.std = MyStd(agriculture),
    artificial.std  = MyStd(artificial),
    natural.std     = MyStd(natural),
    water.std       = MyStd(water),
    wetlands.std    = MyStd(wetlands),
    year.std        = MyStd(year),
    doy.std         = MyStd(doy)
)



# Section 7: NB GLM with spatial dependency ----

#' We will execute the following model.
#' Model 2. log(mu_ij) =  Intercept + covariates + u_i


#' We will implement the following 8 steps:
#'  1. Make a mesh.
#'  2. Define the weighting factors a_ik (i.e. the projector matrix A).
#'  3. Define the SPDE.
#'  4. Define the spatial field.
#'  5. Make a stack. It tells INLA at which points on the mesh we sampled
#'     the response variable and the covariates.
#'  6. Specify the model formula in terms of the response variable,
#'     covariates and the spatial correlated term.
#'  7. Run the spatial model in INLA.
#'  8. Inspect the results.


#* Subsection 7.1: Make a mesh ----

#' We first need to get a sense what the distances are between the
#' sampling locations.
Loc <- cbind(df$Xkm, df$Ykm)
head(Loc)

#' This is in km. Avoid using coordinates with large values. Use km instead of meters!

#' Distances between sites (i.e. trees).
D <- dist(Loc)
png("Plots/FigureS1_distances_between_sites.png", width = 10, height = 10, units = "in", bg = "white", res = 300)
par(mfrow = c(1,1), mar = c(5,5,2,2))
hist(D,
     freq = TRUE,
     main = "",
     xlab = "Distance between sites (km)",
     ylab = "Frequency")
abline(v = 75, col = "red", lwd = 2, lty = 2)
dev.off()

#' Small scale distances for these data is anything between 0 and 100-ish km (maybe around 50km????)


#' Next we make the mesh.
#' Select a value for range guess. After trying various values, and
#' running the entire code in this file I settled on the following value.
RangeGuess <- 50
MaxEdge    <- RangeGuess / 5

#' One of the reasons for not using a smaller value for RangeGuess
#' is computing time for the spatial-temporal models. We could use
#' RangeGuess <- 50, and let the computer run for a full night.

#' Make the mesh.
mesh1 <- fm_mesh_2d(loc = Loc,
                    max.edge = c(1, 5) * MaxEdge,
                    cutoff = MaxEdge / 5)

#' Here is the mesh we just created.
par(mfrow = c(1, 1), mar = c(0, 0, 0, 0))
plot(mesh1, asp = 1, main = "")
points(Loc, col = 2, pch = 16, cex = 1)


#' This is the size of the mesh:
mesh1$n
#' That is a smallish mesh size. But for the spatial-temporal
#' models with 10 years, we will end up with 10 * mesh1$n nodes!
#' That is:
10 * mesh1$n   #' That is a lot, but managable


#' We can also use ggplot to plot the mesh.
ggplot() +
  theme_minimal() +
  labs(title = "Border in UTM") +
  geom_fm(data = mesh1) +
  geom_point(data = df,
             aes(x = Xkm,
                 y = Ykm),
             alpha = 0.5)
#' Nice....but we are consuming a fair amount of vertices
#' in the southern part of the study area. Can we improve this?
#' Those are also in the sea!


#* Subsection 7.2: Make another (convex) mesh ----

#' Make another mesh, with few triangles in the southern part.

#' fm_nonconvex_hull makes a non-convex area around the sampling
#' locations. You can control how close the blue line is to the
#' sites with the 'convex' argument. Useful if you have isolated
#' patches of sampling locations.
NonConvexHull <- fm_nonconvex_hull(Loc, convex = -0.08)
plot(NonConvexHull)
points(Loc)


RangeGuess <- 50  #' Again, we could do with a smaller value,
                  #' but that would increase computing time for the
                  #' spatial-temporal models. It was initially 150

#' Recommended settings
MaxEdge <- RangeGuess / 5
mesh2 <- fm_mesh_2d(boundary = NonConvexHull,
                    max.edge = c(1, 5) * MaxEdge,
                    cutoff   = MaxEdge / 5)
#' max.edge: Maximum allowed triangle edge lengths in
#'           the inner domain and in the outer extension
#' cutoff:   Minimum allowed distance between points. Points
#'           at a closer distance than the supplied value are
#'           replaced by a single vertex


#' Use ggplot to plot mesh2b.
ggplot() +
  theme_minimal() +
  labs(title = "Border in UTM") +
  geom_fm(data = mesh2) +
  geom_point(data = df,
             aes(x = Xkm,
                 y = Ykm),
             alpha = 0.5)
#' Note that the blue line now follows the spatial locations of
#' the sampling locations.

mesh2$n #' Fine for the moment, but we need to multiply this with
        #' 10 for the spatial-temporal models.
#' That is:
10 * mesh2$n   #' That is a lot, but again, manageable


#* Subsection 7.3: Define the projector matrix A----

#' We now define the projector matrix. This is used to
#' calculate:  u = A * w.
A1 <- inla.spde.make.A(mesh1, loc = Loc)
A2 <- inla.spde.make.A(mesh2, loc = Loc)
dim(A1)  #' 2089 sites and 2960 nodes in the mesh.
dim(A2)  #' 2089 sites and 2871 nodes in the mesh.


#* Subsection 7.4: Define the SPDE----

#' We need to specify priors for the two Matern correlation
#' parameters Range and sigma. This is not the RangeGuess that we
#' defined earlier.

#' In short: The user needs to select range0, sigma0 and alpha in:
#'   P(Range < range0) = alpha
#'   P(sigma > sigma0) = alpha

#' These are defined in INLA as:
#'    prior.range = c(range0, alpha)
#'    prior.sigma = c(sigma0, alpha)

#' Prior for the range.
#'  I have no idea about the dispersal capacity of dipteran vectors
#'  Also, their dispersal capacity may differ by species / group.
#'  Thus, this is a bit a blind guess, but I will use:
#'  P(Range < 75) = 0.5
#'  This states that it we do not know what the range is and there is an equal probability that it is smaller or larger than 75 km.

#' Prior for sigma.
#'  Here is a quick and dirty way to get an impression of
#'  likely values for sigma. This trick only works for models
#'  with a log-link function.

#'  1. Log-transform the Y data:
df$LogCounts <- log(df$Counts + 1)

#' 2. Apply a lm on LogVec_abund with only an intercept.
Test <- lm(LogCounts ~ 1, data = df)
summary(Test)

#' This output can be written as:
#'  LogCounts_i = 0.89291 + Residuals_i
#'  Residual_i ~ N(0, 1.419 ^2)

#' And this can be written as:
#'   LogCounts_i = exp(0.89291 + Residuals_i)
#'
#' It is actually a little bit more complicated (look up a
#' lognormal distribution). But this gives an indication that
#' the u_i ~ N(0, 1.419^2) would be a decent starting point.
#' Maybe we should be conservative, and use:
#'  u_i ~ N(0, 2^2)

#' This gives the following prior:  P(sigma_u > 2) = 0.05
#' This states that it is unlikely that sigma_u is larger than 2.


#' Summarising, we will use the following PC priors:

#'  P(Range < 75 km) = 0.5
#'  This states that it we do not know what the range is and there is an equal probability that it is smaller or larger than 75 km.
#'  75km is small scale for our study

#'  P(sigma > 2) = 0.05
#'  This states that it is unlikely that sigma is larger than 2.


#' In INLA coding, this is:
spde1 <- inla.spde2.pcmatern(mesh1,
                             prior.range = c(75, 0.5),
                             prior.sigma = c(2, 0.05))

spde2 <- inla.spde2.pcmatern(mesh2,
                             prior.range = c(75, 0.5),
                             prior.sigma = c(2, 0.05))


#* Subsection 7.5: Define the spatial field----

#' Next, we define the spatial random intercept u.
#' For mesh1 we use:
w1.index <- inla.spde.make.index(name = 'w',
                                 n.spde = spde1$n.spde,
                                 n.group = 1,
                                 n.repl = 1)

#' For mesh2 we use:
w2.index <- inla.spde.make.index(name = 'w',
                                 n.spde = spde2$n.spde,
                                 n.group = 1,
                                 n.repl = 1)



#* Subsection 7.6: Make a stack----

#' Make the X matrix using model.matrix()
X <- model.matrix(~ eqr.std + ppt.std + tmin.std + ws.std + elevation.std +
                    agriculture.std + artificial.std + natural.std +
                    year.std + waterbody_type, data = df)
X <- as.data.frame(X) #' Avoids an error message in INLA
colnames(X)


#' This is not relevant for this specific model:
#' In INLA, you can't use ':' in the variable names. Note that the
#' interaction term uses : in its name. Replace the : by '_'.
#OldColnames  <- colnames(X)
#NewColNames  <- gsub(pattern = ":", replacement = "_", x = OldColnames)
#colnames(X)  <- NewColNames
#head(X)


#' Sample size
N <- nrow(df)


#' We now define the stack for mesh2a.
Stack1 <- inla.stack(
  tag = "Fit",                       #' Name of stack
  data = list(Counts = df$Counts), #' Response variable
  A = list(1, 1, A1),               #' The order matches the order in effects.
  effects = list(
    Intercept = rep(1, N),  #' Use our own intercept
    X         = X[,-1],     #' Dump the default intercept from the model.matrix
    w         = w1.index)) #' Spatial random field


#' And this the stack based on mesh2b.
Stack2 <- inla.stack(
  tag = "Fit",                       #' Name of stack
  data = list(Counts = df$Counts), #' Response variable
  A = list(1, 1, A2),               #' The order matches the order in effects.
  effects = list(
    Intercept = rep(1, N),       #' Use our own intercept
    X         = X[,-1],          #' Dump the default intercept from the model.matrix
    w         = w2.index))      #' Spatial random field


#* Subsection 7.7: Specify the model formula----

#' Specify the model formula in terms of the response variable,
#' covariates and the spatial correlated term. Having the colnames
#' is handy at this stage:
colnames(X)


#' This is a model without spatial dependency.
f0 <-  Counts ~ eqr.std + ppt.std + tmin.std + ws.std + elevation.std +
                agriculture.std + artificial.std + natural.std +
                year.std + waterbody_type


#' This is a model with spatial dependency, based on mesh 1.
f1 <-  Counts ~ -1 + Intercept + eqr.std + ppt.std + tmin.std + ws.std + elevation.std +
                     agriculture.std + artificial.std + natural.std +
                     year.std + waterbody_type.L +
                     f(w, model = spde1)


#' This is a model with spatial dependency, based on mesh2.
f2 <-  Counts ~ -1 + Intercept + eqr.std + ppt.std + tmin.std + ws.std + elevation.std +
                     agriculture.std + artificial.std + natural.std +
                     year.std + waterbody_type.L +
                     f(w, model = spde2)


#* Subsection 7.8: Execute the INLA models----

#' This is the NB GLM without spatial dependency.
I0 <- inla(f0,
           family = "nbinomial",
           data = df,
           control.compute = MyControlCompute)


#' This is the spatial NB GLM based on mesh1.
I1 <- inla(f1,
            family = "nbinomial",
            data = inla.stack.data(Stack1),
            control.compute = MyControlCompute,
            control.predictor = list(A = inla.stack.A(Stack1)))


#' This is the spatial NB GLM based on mesh2.
I2 <- inla(f2,
           family = "nbinomial",
           data = inla.stack.data(Stack2),
           control.compute = MyControlCompute,
           control.predictor = list(A = inla.stack.A(Stack2)))


#* Subsection 7.9: Compare the INLA models----

#' Compare DIC and WAIC.
DICs <- c(I0$dic$dic, I1$dic$dic, I2$dic$dic)
WAICs <- c(I0$waic$waic, I1$waic$waic, I2$waic$waic)
Results <- data.frame(Models = c("NB GLM",
                                 "NB GLM + SRF with mesh1",
                                 "NB GLM + SRF with mesh2"),
                      DIC = DICs,
                      WAIC = WAICs)
Results

#' Conclusions:
#'  - Adding spatial dependency improves the models.
#'  - Mesh1 is better than mesh2, but the difference is within 10

#' Therefore, we will use mesh1 in the remaining part of this analysis.

summary(I1)



# Section 8: Imposed spatial dependency----

#' Here is some code to extract the posterior mean values of the
#' spatial parameters.
SpatialParams <- MySpatialParams(Model = I1,
                                 ThisSpde = spde1)

Kappa   <- SpatialParams[1]  #' Parameter kappa for the Mattern correlation function
Sigma.u <- SpatialParams[2]  #' Sigma of the u expressed in metres.
Range   <- SpatialParams[3]  #' Range.

#' This is the important part:
Sigma.u   #' u_i ~ N(0, 1.085252^2 * Spatial correlation)
Range     #' Distance at which the correlation diminishes.
          #' 27.51841. That is not as large as we thought!


#' Visualise the correlation structure.
#' First we obtain the locations of each point of the mesh.
LocMesh <- mesh1$loc[,1:2]

#' And then we calculate the distance between each vertex.
D <- as.matrix(dist(LocMesh))

#' Using the estimated parameters from the model (see above)
#' we can calculate the imposed Matern correlation values.
d.vec <- seq(0, max(D), length = 100)
Cor.M <- (Kappa * d.vec) * besselK(Kappa * d.vec, 1)
Cor.M[1] <- 1

png("Plots/FigureS6_imposed_matern_correlation.png", width = 10, height = 10, units = "in", bg = "white", res = 300)
#' Which we plot here:
par(mfrow=c(1,1), mar = c(5,5,2,2))
plot(x = d.vec,
     y = Cor.M,
     pch = 16,
     type = "l",
     cex.lab = 1.5,
     xlab = "Distance (km)",
     ylab = "Correlation",
     xlim = c(0, 100))
abline(h = 0.8, lty = 2, col = 2)
abline(h = 0.5, lty = 2, col = 2)
abline(h = 0.1, lty = 2, col = 2)
abline(v = 5.5, lty = 2, col = 4)
abline(v = 12.5, lty = 2, col = 4)
abline(v = 32, lty = 2, col = 4)
dev.off()

#' Define strong correlation as correlation between 1 - 0.8.
#' Define moderate correlation as correlation between 0.8 - 0.5.
#' Define weak correlation as correlation between 0.5 - 0.1.
#' Define diminishing correlation as correlation smaller than 0.1.

#' In this case, we have:
#'  - Strong correlation between sites separated between 0 - 5.5 km.
#'  - Moderate correlation between sites separated between 5.5 - 12.5 km.
#'  - Weak correlation between sites separated between 12 - 32 km.
#'  - Diminishing correlation between sites separated by more than 32 km



# Section 9: Plotting the spatial dependency----

#' Next we present the spatial component, the wks.
#' Their posterior mean values can be obtained via
w.pm <- I1$summary.random$w$mean
length(w.pm)


#' This is a vector of length 1674 by 1. Each value in w.pm belongs
#' to a specific vertex on mesh 1. It is tempting to plot these,
#' but that won't work as geom_raster() wants to have data on a
#' regular grid. INLA has nice tools to predict the w.pm values on a
#' grid.

#' More recent INLA versions use fm_evaluator() and fm_evaluate().
#' We first need to call fm_evaluator() followed by a call to
#' fm_evaluate(). See below.
#'  - The dims = c(400, 400) specifies a 400-by-400 grid within
#'    the ranges of the mesh.
#'  - Proj contains the A matrix for this 400-by-400 grid.
#'  - fm_evaluate() gives the values of the SRF on a regular grid
#'    defined by fm_evaluator().
Proj <- fm_evaluator(mesh = mesh1,
                     dims = c(400, 400))
#' Proj$x contains the 400 values along the x-axis.
#' Proj$y contains the 400 values along the y-axis.
#' Proj$proj$A contains the projector matrix.

SRF.Proj <- fm_evaluate(projector = Proj,
                        field  = w.pm)
#' SRF.Proj contains SRF values on the grid that was just defined.

#' Extract the relevant information so that we can use ggplot
#' code from previous exercises to plot the results
MyData      <- expand.grid(Xkm = Proj$x, Ykm = Proj$y)
MyData$w.pm <- as.vector(SRF.Proj)



#' We will now plot SRF, but we want to omit those parts
#' of the SRF that are outside the study area (i.e. in other countries surrounding Lithuania).
#' To do this, we carry out the following steps.
#'  1. Convert MyData to a sf object.
#'  2. Determine which points in inside the study. Use st_contains() for this.
#'  3. Set the SRF values of the sites outside the study area to NA.


#' 1. Convert MyData to a sf object.
MyData.sp <- st_as_sf(MyData,
                      coords = c("Xkm", "Ykm"),
                      crs = "+proj=utm +zone=34 +north +ellps=WGS84 +units=km  +datum=WGS84")

#' Determine which points in MyData are contained within the study area
Contained <- st_contains(CroppedLithuania_UTM, MyData.sp)

#' Extract indices of the points that are on land.
contained_indices <- unlist(Contained)

#' Create a logical vector where TRUE indicates the point is in Lithuania
OnLand <- rep(FALSE, nrow(MyData))
OnLand[contained_indices] <- TRUE
OnLand

#' 3. Set the values of the SRF that are not on land to NA.
MyData$w.pm[!OnLand]  <- NA

#' Here is a method to get the map in the same range as the panel plot:

# Create a data frame with the corner points of your desired plot area
desired_corners <- data.frame(
  lon = c(20.9, 26.9, 26.9, 20.9),
  lat = c(53.85, 53.85, 56.5, 56.5)
)

# Convert to sf object with WGS84 CRS
desired_corners_sf <- st_as_sf(desired_corners,
                              coords = c("lon", "lat"),
                              crs = 4326) # WGS84

# Transform to UTM Zone 34 North (same as your data)
desired_corners_utm <- st_transform(desired_corners_sf,
                                   "+proj=utm +zone=34 +north +ellps=WGS84 +units=km +datum=WGS84")

# Extract the coordinates
utm_coords <- st_coordinates(desired_corners_utm)
utm_xlim <- range(utm_coords[,"X"])
utm_ylim <- range(utm_coords[,"Y"])

range(MyData$w.pm, na.rm = TRUE)

# convert prosterior means of w's to percentage
MyData <- MyData |>
  mutate(w.pm_converted = ((exp(w.pm) - 1) * 100))

range(MyData$w.pm_converted, na.rm = TRUE)


# Now use these limits in your plot
ggplot() +

  # geom_sf(data = st_sf(SeaPolygon_UTM),
  #         color = "lightblue",
  #         fill = "lightblue",
  #         linewidth = 0.5) +
  #
  # geom_sf(data = st_sf(CroppedSurroundings_UTM),
  #         color = "black",
  #         fill = "white",
  #         linewidth = 0.5) +

  geom_raster(data = MyData,
              aes(x = Xkm,
                  y = Ykm,
                  fill = w.pm)) +

  geom_sf(data = CroppedLithuania_UTM,
          fill = NA,
          color = "black",
          linewidth = 0.5) +

  scale_fill_gradient2(
    name = "Posterior means (w)",
    limits = c(min(MyData), max(w.pm)),
    midpoint = 0,
    low = "#21918c",
    mid = "white",
    high = "#440154",
    na.value = NA,
    guide = guide_colorbar(
      barwidth = 20,
      barheight = 0.5,
      title.position = "top",
      title.hjust = 0.5,
      label.hjust = 0.5
    )
  ) +
  # Set the UTM coordinate limits that correspond to your desired lat/long bounds
  coord_sf(xlim = utm_xlim, ylim = utm_ylim) +
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

ggsave("Plots/Figure6_SRF_spatial_dependency.png", width = 10, height = 8, units = "in", bg = "white", dpi = 300)
#' This is quite clear spatial correlation.


# Section 10: Results of spatial-temporal GLM with mesh 1----

#' Numerical output for the regression parameters.
BetasI1 <- I1$summary.fixed[, c("mean", "0.025quant", "0.975quant")]
print(BetasI1, digits = 2)


#' Posterior mean of theta (as in: mu + mu^2 / theta)
theta.pd <- I1$marginals.hyperpar$`size for the nbinomial observations (1/overdispersion)`
theta.pm <- inla.emarginal(function(x) x, theta.pd)
theta.pm

#' Visualise the betas using ggplot2.
#' 1. Create a data frame.
BetasI1_df <- data.frame(Covariate = rownames(BetasI1),
                         Mean      = BetasI1[,"mean"],
                         Lower     = BetasI1[,"0.025quant"],
                         Upper     = BetasI1[, "0.975quant"])

#' 2. Explicitly set the order of the covariates to match their original order
BetasI1_df$Covariate <- factor(BetasI1_df$Covariate, levels = rownames(BetasI1))

#' 3. Plot the info in BetasI1_df.
covariate_labels <- c(
  "Intercept" = "Intercept (Rivers)",
  "eqr.std" = "Ecological Quality Ratio (EQR)",
  "ppt.std" = "Precipitation (mm)",
  "tmin.std" = "Minimum Temperature (°C)",
  "ws.std" = "Wind Speed (m/s)",
  "elevation.std" = "Elevation (m.a.s.l.)",
  "agriculture.std" = "Agricultural Land (%)",
  "artificial.std" = "Artificial Surface (%)",
  "natural.std" = "Natural Area (%)",
  "year.std" = "Year",
  "waterbody_type.L" = "Lakes"
)

# First, let's add a column to identify significant effects
BetasI1_df <- BetasI1_df|>
  mutate(is_significant = !(Lower <= 0 & Upper >= 0))

# Create a function to format the labels based on significance
format_labels <- function(variable, is_significant, labels_map) {
  result <- labels_map[variable]
  result[is_significant] <- sprintf("<b>%s</b>", result[is_significant])
  return(result)
}

# Join significance information with covariate names for labeling
label_data <- BetasI1_df |>
  dplyr::select(Covariate, is_significant) |>
  mutate(label = format_labels(Covariate, is_significant, covariate_labels))

# Modified plot with bold labels for significant covariates
ggplot(BetasI1_df, aes(x = Mean, y = Covariate)) +
  # Add interval lines with varying thickness
  geom_linerange(aes(xmin = Lower, xmax = Upper,
                     linewidth = is_significant,
                     colour = is_significant),
                     show.legend = FALSE) +
  geom_pointrange(aes(xmin = Lower, xmax = Upper,
                      color = is_significant),
                      fatten = 0.5,
                      size = 10,
                      show.legend = FALSE) +
  geom_vline(aes(xintercept = 0),
             linetype = "dashed",
             color = "gray40",
             linewidth = 0.8) +
  labs(x = "Standardized Effect Size",
       y = NULL) +
  theme_minimal() +
  My_theme +
  # Use scale_y_discrete with richtext for HTML formatting
  scale_y_discrete(limits = BetasI1_df$Covariate,
                   labels = setNames(label_data$label, label_data$Covariate)) +
  scale_linewidth_manual(values = c("TRUE" = 1.5, "FALSE" = 1)) +
  scale_color_manual(values = c("gray75", "#21918c")) +
  guides(linewidth = "none") +
  # Add this theme element to render HTML in text
  theme(axis.text.y = ggtext::element_markdown())

ggsave("Plots/Figure4_fixed_effects.png", width = 10, height = 10, bg = "white", dpi = 300)

# Create a simplified table of model output
result_table <- create_simple_inla_table(I1)
result_table

save_kable(result_table, "Plots/TableS2_model_results.html") # <---- save table as png manually (go to 'export' within the viewer window)

# Section 11: Model validation for the INLA NB GLM I1 ----

#* Subsection 11.1: Getting scaled quantile residuals ----

#' We would like to simulate 1000 data sets, because:
#'  1. We can use these to assess whether the model is overdispersed.
#'  2. In case of zero-inflation, we can use the 1000 simulated data
#'     sets to assess whether the model can cope with the zero
#'     inflation.
#'  3. We can use them to obtain scaled quantile residuals.


#' Do the posterior simulation of regression parameters 1000 times.
SimData <- inla.posterior.sample(n = 1000, I1)


#' Extract the 1000 sets of regression parameters
MyVar <- rownames(I2$summary.fixed)
MyVar    #' Names of the regression parameters
Betas1000 <- inla.posterior.sample.eval(MyVar,
                                        SimData)
Betas1000[,1] #' First set of simulated betas
Betas1000[,2] #' Second set of simulated betas
Betas1000[,3] #' Third set of simulated betas


#' Now the spatial correlation random effects.
#' Determine the names of the random effects.
#' In the code below, we are using:  ^w:\d+$
#' This will only match strings that start with "w:", followed
#' exclusively by one or more digits, and then end.
names.SimData <- attributes(SimData[[1]]$latent)$dimnames[[1]]
names.SimData
wNames <- names.SimData[grep("^w:\\d+$", names.SimData)]

#' This is probably easier: wNames <- paste("w",1:mesh3$n, sep = ":")
#' Determine on which rows the random effects are.
MyID  <- function(x){ which(rownames(SimData[[1]]$latent) == x) }
RowsW <- lapply(wNames, MyID)
RowsW <- as.numeric(RowsW)
RowsW
NROW(RowsW)
mesh1$n


#' Calculate 1000 times the fitted values via: mu = exp(X * beta + A1 * w).
#' Get the design matrix X.
X <- model.matrix(~ eqr.std + ppt.std + tmin.std + ws.std + elevation.std +
                    agriculture.std + artificial.std + natural.std +
                    year.std + waterbody_type, data = df)
X <- as.matrix(X)


#' Get the theta from the selected NB GLM:
theta.pd <- I1$marginals.hyperpar$`size for the nbinomial observations (1/overdispersion)`
theta.pm <- inla.emarginal(function(x) x, theta.pd)
theta.pm


#' Start a loop. In each iteration:
#'   1. Calculate the fitted values mu.
#'   2. Simulate NB data with the mean equal to mu.
N    <- nrow(df)                        #' Sample size
Ysim <- matrix(nrow = N, ncol = 1000)   #' Create space

#' Start the loop (can be done more efficient with lapply).
for (i in 1:1000){
  Betas <- Betas1000[,i]
  wk    <- SimData[[i]]$latent[RowsW]
  FixedPart   <- X %*% Betas
  SpatialPart <- A1 %*% wk    #' Note the A2a. This is based on mesh2a.
  mu  <- exp(FixedPart + SpatialPart)[,1]
  Ysim[,i] <- rnegbin(n = N,                   #' Simulated NB Owl count data
                      mu = mu,
                      theta = theta.pm)
}

#' We now have 1000 simulated data sets from the model.
Ysim[, 1] #' First simulated data set.
Ysim[, 2] #' Second simulated data set.
Ysim[, 3] #' Third simulated data set.
#' Etc.

#' Or:
par(mfrow = c(2,2))
hist(df$Counts, main ="Observed vector abundance")
hist(Ysim[,1], main = "First simulated data set")
hist(Ysim[,2], main = "Second simulated data set")
hist(Ysim[,3], main = "Third simulated data set")
par(mfrow = c(1,1))
#' Quite similar! Can we improve? Still, no idea whether the observed data is similar from the simulated data.



#' Now we have 1000 simulated data sets. Give them all to
#' DHARMa, and it will calculate scaled quantile residuals.
N <- nrow(df)
df$Fit2 <- I1$summary.fitted.values[1:N, "mean"]  #' Fitted values.

E2.sqr <- createDHARMa(simulatedResponse = Ysim,
                       observedResponse = df$Counts,
                       fittedPredictedResponse = df$Fit2,
                       integerResponse = TRUE)
#' Now we have scaled quantile residuals.



#* Subsection 11.2: Check for overdispersion----

#' We will use the scaled quantile residuals to assess for
#' overdispersion.
par(mfrow = c(1,1), mar = c(5,5,2,2))
testDispersion(E2.sqr)
#' Fine!


#* Subsection 11.3: Check for homogeneity of variance----

#' Plot the scaled quantile residuals versus (ranked) fitted values.
plotResiduals(E2.sqr, quantreg = TRUE, smoothScatter = FALSE)
#' Some trouble, but tolerable :)


#* Subsection 11.4: Check for uniformity of the residuals----

#' In DHARMa, we verify whether the scaled quantile residuals are
#' uniform distributed.
par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(E2.sqr, testUniformity = TRUE,
           testOutliers = TRUE, testDispersion = FALSE)
#' Fine!


#* Subsection 11.5: Plot residuals versus the covariates----

#' Plot the scaled quantile residuals versus each covariate
#' in the model.
plotResiduals(E2.sqr, form = df$eqr, quantreg = TRUE)            #' Tolerable
plotResiduals(E2.sqr, form = df$ppt, quantreg = TRUE)            #' Fine
plotResiduals(E2.sqr, form = df$tmin, quantreg = TRUE)           #' Fine
plotResiduals(E2.sqr, form = df$ws, quantreg = TRUE)             #' Fine
plotResiduals(E2.sqr, form = df$elevation, quantreg = TRUE)      #' Tolerable
plotResiduals(E2.sqr, form = df$agriculture, quantreg = TRUE)    #' Fine
plotResiduals(E2.sqr, form = df$artificial, quantreg = TRUE)     #' Fine
plotResiduals(E2.sqr, form = df$natural, quantreg = TRUE)        #' Fine
plotResiduals(E2.sqr, form = df$year, quantreg = TRUE)           #' Tolerable
plotResiduals(E2.sqr, form = df$year, asFactor = TRUE, quantreg = TRUE)           #' Tolerable
plotResiduals(E2.sqr, form = df$waterbody_type, asFactor = TRUE, quantreg = TRUE) #' Fine


#* Subsection 11.6: Check for zero inflation----
testZeroInflation(E2.sqr)
#' Fine!


#* Subsection 11.7: Check for spatial dependency----

#' Option 1: Plot the residuals vs spatial locations. Look for patterns.
#' Option 2: Apply Moran's I test on the residuals.
#' Option 3: Make a variogram of the residuals


#' Option 1:  Plot the residuals vs spatial locations
#' Use scaled quantile residuals:
df$E2 <- residuals(E2.sqr) # <-- better if you have zero inflation

#' #' Or use Pearson residuals.
#' df$E2 <- (df$Counts - df$Fit2) / sqrt(df$Fit2 + df$Fit2^2 / theta.pm)

#' Note that scaled quantile residuals are around 0.5, whereas
#' Pearson residuals are around 0.
#' Define colour and point size.
df$MyCol  <- ifelse(df$E2 >= 0.5, "red", "blue") # <- if you use pearson, then it should be 0, not 0.5
df$MySize <- rescale(abs(df$E2), to = c(0, 3))

p <- ggplot() +
  geom_point(data = df,
             aes(x = Xkm,
                 y = Ykm,
                 col = MyCol,
                 size = MySize),
             alpha = 0.5)  +
  scale_color_identity() +
  scale_size_continuous(range = c(1, 3)) +
  theme_minimal() +
  theme(legend.position = "none") +
  guides(fill = guide_legend(title = NULL)) +
  labs(title = "Residuals")
p
p + facet_grid(~waterbody_type)
#' Fine! No clear pattern.


#' Option 2: Apply Moran's I test on the residuals
#' DHARMa can test for spatial correlation using Moran's I
#' test.
# testSpatialAutocorrelation(df$E2, #E2.sqr,
#                            x = df$Xkm,
#                            y = df$Ykm,
#                            plot = FALSE)

# test spatial autocorrelation
groupLocations = aggregate(df[, c("Xkm", "Ykm")], list(df$site_id), mean)
res_space = recalculateResiduals(E2.sqr, group = df$site_id)
testSpatialAutocorrelation(res_space, groupLocations$Xkm, groupLocations$Ykm)
#' Some residual spatial dependency, but not too much.

#' Option 3: Make a variogram.
#' Use scaled quantile residuals:
df$E2 <- residuals(E2.sqr) # <-- better if you have zero inflation
#' Make a sample variogram of the residuals.
MyData <- data.frame(E2  = df$E2,
                     Xkm = df$Xkm,
                     Ykm = df$Ykm)

#' Convert to sf object.
MyData_sf <- st_as_sf(x = MyData,
                      coords = c("Xkm", "Ykm"),
                      crs = NA)  #' Non-Cartesian coordinates.

#' Apply the variogram function from gstat.
V1 <- variogram(E2 ~ 1,
                data = MyData_sf,
                # cutoff = 150,
                cressie = TRUE)

#' Plot the variogram
p1 <- ggplot(data = V1, aes(x = dist, y = gamma)) +
         geom_point() +
         geom_smooth(se = FALSE) +
         labs(x = "Distance (in km)", y = "Semi-variogram") +
         theme(text = element_text(size=15),
               legend.position="none")
p1
#' Conclusion:
#'   -That is a strange pattern
#'   -But, observe the y-axis scale. The values are very small.
#    -Next is to do autoregressive spatial model temporal with AR1 structure


#* Subsection 11.9: Conclusions model validation----

#' We have no clear patterns in the residuals of all the covariates.
#' We no longer have zero inflation. Negative binomial distribution can handle this.
#' The spatial autocorrelation is a a bit weird, but if we look at it in perspective of the y-axis, it is not that bad.



# Section 12: Compare models I0, I1 and I2 ----

#' Let's plot the results of the model, without and with the
#' spatial correlation side by side.
Out1 <- I0$summary.fixed[,c("mean", "0.025quant", "0.975quant")]
Out2 <- I1$summary.fixed[,c("mean", "0.025quant", "0.975quant")]
Out3 <- I2$summary.fixed[,c("mean", "0.025quant", "0.975quant")]
rownames(Out1) <- rownames(Out2) <- rownames(Out3) <- I1$names.fixed

#' Names for the models:
MyNames <- c("GLM",
             "Spatial GLM mesh 1",
             "Spatial GLM mesh 2")

#' Compare results of the two models using MyCompareBetasofModels
#' (which is in our support file).
MyCompareBetasofModels(AllModels = list(Out1, Out2, Out3),
                       ModelNames = MyNames)

#' There are some differences between the model without, and with
#' spatial dependency. But, the models with spatial dependency are
#' not so different, though they have slightly larger credible intervals
#' compared to the model without spatial dependency structures.



# Section 13: NB GLM with replicate spatial-temporal term----

#* Subsection 13.1: Define blocking structure----

#' Define a vector that identifies the blocking structure.
df$YearNum <- as.numeric(as.factor(df$year))
Repl <- df$YearNum
Repl
#' This is a vector with values 1 1 1 2 2 2 3 3 3 3... etc.
#' All observations in the same block will be spatially correlated.
#' The correlation itself can differ per block, but all blocks
#' share the same range and sigma_u.


#' How many blocks do we have? Should be 10.
NRepl <- length(unique(Repl)) #' Number of groups.
NRepl


#* Subsection 13.2: Define the projector matrix A----

#' Define the projector matrix. We can use
#' mesh2a for this. No need to change the mesh.
A3Repl <- inla.spde.make.A(mesh1,
                           loc = Loc,
                           repl = Repl)  #' 18 blocks.
dim(A3Repl)
#' The rows refer to the sites.
#' The columns are weighting factors.
#' Weighting factors 10 blocks.
#' recall: u = A * w



#* Subsection 13.3: Define spde----

#' Define the SPDE. We use the same priors for the range and sigma_u
#' as in the previous section.
spde3Repl <- inla.spde2.pcmatern(mesh1,
                                 prior.range = c(75, 0.5),
                                 prior.sigma = c(2, 0.05))


#* Subsection 13.4: Define SRF w----

#' Define the SRF
w3Repl <- inla.spde.make.index('w',
                               n.spde = mesh1$n,
                               n.repl = NRepl)


#* Subsection 13.5: Stack for the replicate spatial-temporal GLM----

#' Make the X matrix using model.matrix()
X <- model.matrix(~ eqr.std + ppt.std + tmin.std + ws.std + elevation.std +
                    agriculture.std + artificial.std + natural.std +
                    year.std + waterbody_type, data = df)
X <- as.data.frame(X) #' Avoids an error message in INLA
colnames(X)


#' This is not relevant for this specific model:
#' In INLA, you can't use ':' in the variable names. Note that the
#' interaction term uses : in its name. Replace the : by '_'.
#OldColnames  <- colnames(X)
#NewColNames  <- gsub(pattern = ":", replacement = "_", x = OldColnames)
#colnames(X)  <- NewColNames
#head(X)


#' We now define the stack. This is INLA's way to combine data that
#' has been sampled at different spatial resolutions. .

#' Sample size
N <- nrow(df)

Stack3SpatTemp <- inla.stack(
  tag = "Fit",                       #' Name of stack
  data = list(Counts = df$Counts), #' Response variable
  A = list(1, 1, A3Repl),            #' The order matches the order in effects.
  effects = list(
    Intercept = rep(1, N),           #' Use our own intercept
    X         = X[,-1],              #' Dump the intercept from the model.matrix
    w         = w3Repl))             #' Spatial random field

#' This is the stack for the NB GLM with replicate spatial-temporal correlation.


#* Subsection 13.6: Specify the model formula----

#' Specify the model formula in terms of the response variable,
#' covariates and the spatial correlated term. Having the colnames
#' is handy at this stage:
colnames(X)


#' This is the model with spatial dependency, based on mesh 2a.
f3 <-  Counts ~ -1 + Intercept + eqr.std + ppt.std + tmin.std + ws.std + elevation.std +
                     agriculture.std + artificial.std + natural.std +
                     year.std + waterbody_type.L +
                 f(w, model = spde3Repl, replicate = w.repl) # <---  w.repl is correct, should not be w3Repl



#* Subsection 13.7: Execute the replicate spatial-temporal GLM----

#' This is the NB GLM with replicate spatial dependency.
I3 <- inla(f3,
           family = "nbinomial",
           data = inla.stack.data(Stack3SpatTemp),
           control.compute = MyControlCompute,
           control.predictor = list(A = inla.stack.A(Stack3SpatTemp)))


#* Subsection 13.8: Compare the models----

#' Compare DIC and WAIC.
DICs <- c(I0$dic$dic, I1$dic$dic, I2$dic$dic, I3$dic$dic)
WAICs <- c(I0$waic$waic, I1$waic$waic, I2$waic$waic, I3$waic$waic)
Results <- data.frame(Models = c("NB GLM",
                                 "NB GLM + SRF with mesh1",
                                 "NB GLM + SRF with mesh2",
                                 "NB GLM + replicate SRF with mesh2"),
                      DIC = DICs,
                      WAIC = WAICs)
Results

#' Conclusion:
#'  - Adding replicate spatial-temporal dependency does not
#'    improve the model.


#* Subsection 13.9: Present the replicate SRF----

#' The posterior mean values of the replicate SRF can be obtained
#' via:
w.pm <- I3$summary.random$w$mean
length(w.pm)

#' For each year we have a set of w's. Hence, we can make 10 pictures;
#' one for each year. We want to have 10 panels in 1 ggplot graph.
#' Each panel should show the SRF for that specific year.

#' First we need the levels of fYear.
NamesYear <- levels(as.factor(df$year))
NamesYear

#' Create an object in which we can store the 10 sets of w.pm on a
#' 100-by-100 grid.
MyData.All <- NULL

#' Start a loop.
#'  1. Extract the SRF (i.e. the w's) for year i.
#'  2. Project the SRF on a 100-by-100 grid for that year.
#'  3. Store the results in MyData$w.pm.
#'  4. Determine which of the 100 * 100 rows in MyData
#'     are not on land. Set these MyData$w.pm to NA for plotting
#'     purposes.
#'  5. Store the MyData$w.pm in MyData.All

#' We are using a 100-by-100 grid to speed up the ggplot2 plotting.

for (i in 1:NRepl){
  #' 1. Extract the SRF for year i.
  wpm.i <- w.pm[w3Repl$w.repl == i]
  #wsd.i <- wsd[wRepl$w.repl == i]

  #' 2. Project the SRF on a 100-by-100 grid for year i.
  Proj     <- fm_evaluator(mesh = mesh1, dims = c(100, 100))
  SRF.Proj <- fm_evaluate(projector = Proj,  field  = wpm.i)

  #' 3. Store the results in MyData$w.pm.
  MyData      <- expand.grid(Xkm = Proj$x, Ykm = Proj$y)
  MyData$w.pm <- as.vector(SRF.Proj)

  #'  4. Determine which of the 100 * 100 rows in MyData
  #'     are not on land. Set these MyData$w.pm to NA for plotting
  #'     purposes.

  #'     Convert MyData to a sf object.
  MyData.sp <- st_as_sf(MyData,
                        coords = c("Xkm", "Ykm"),
                        crs = "+proj=utm +zone=34 +north +ellps=WGS84 +units=km  +datum=WGS84")

  #'     Determine which points in MyData are contained within the
  #'     study area
  Contained <- st_contains(CroppedLithuania_UTM, MyData.sp)

  #'     Extract indices of the points that are on land.
  contained_indices <- unlist(Contained)

  #'     Create a logical vector where TRUE indicates the point is
  #'     on land.
  OnLand <- rep(FALSE, nrow(MyData))
  OnLand[contained_indices] <- TRUE

  #'     Not on land: set to NA
  MyData$w.pm[!OnLand] <- NA


  #'  5. Store the MyData$w.pm in MyData.All
  MyData$Year <- NamesYear[i]
  MyData.All  <- rbind(MyData.All,MyData)
}

#' Convert Year to a factor for facetting in ggplot.
MyData.All$fYear <- factor(MyData.All$Year)


#' And some ggplot2 coding to plot the results.
ggplot() +
  geom_sf(data = st_sf(SeaPolygon_UTM),
          alpha = 0.5,
          fill = "grey80") + # Sea
  geom_sf(data = CroppedLithuania_UTM,
          fill = "transparent",
          color = "black") + # Land
  geom_raster(data = MyData.All,
              aes(x = Xkm,
                  y = Ykm,
                  fill = w.pm)) +
  scale_fill_gradient2(name = "Spatial smoother",
                       midpoint = 0,
                       low = "#21918c",
                       mid = "white",
                       high = "#440154",
                       na.value = NA) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  labs(title = "Study area") +
  geom_point(data = df,
             aes(x = Xkm, y = Ykm),
             size = 0.1,
             alpha = 0.5) +
  facet_wrap(~ fYear, ncol = 5)
#' From year to year, the spatial random field changes.



# Section 14: AR1 NB GLM----

#' We will now allow for spatial correlation that changes over time
#' according to an auto-regressive process.

#' u_ij = phi * u_{i,j-1} + noise_ij
#' j is the year subscript.


#* Subsection 14.1: Define blocking structure----

#' The ar1 coding is very similar to that of the replicate
#' correlation. First we define a blocking structure.
#' All observations in the same block are spatially correlated.
#' From block to block the changes over time follow the ar1 structure.

#' Blocking structure, not necessarily ordered.
df$YearSeq <- as.numeric(as.factor(df$year))
df$YearSeq

#' We have 10 blocks
NYears <- length(unique(df$YearSeq))
NYears


#' If you run the GLM with spatial-temporal ar1 correlation,
#' then computing time will be about 5 - 10 minutes on a fast
#' computer.

#' We will define time knots. And at each knot, we will
#' calculate a spatial random field. From knot to knot,
#' we will impose the ar1 correlation. As an example:

Knots <- seq(1, 10, by = 1)
Knots

#' The model will now calculate a spatial random field for year = 1,
#' year = 4, year = 7, and year 10.

#' To get the spatial random effects u_ij in year 2, the model will take
#' a weighted average of the SRFs of knots 1 and 4.
#' The weights are given by the time difference (year 2 is closer to
#' the first knot, so its weighting factor is larger that of the second
#' knot).
#'
#' Knots: 1        4        7        10
#'        1  2  3  4  5  6  7  8  9  10

#' u_ij in year 3 =  0.66 * SRF-at-knot-2  +  0.33 * SRF-at-knot-5


#' This is the mesh for the knots.
mesh.t <- inla.mesh.1d(Knots, degree = 1)

#' If you want to see the weighting factors, here they are:
Atest <- inla.spde.make.A(mesh = mesh.t,
                          loc = df$YearSeq)
Atest

#' Summarising: The random effects u_ij in years between the knots
#'              are weighted averages of the two SRF at the
#'              nearest 2 knots.


#' The total numer of knots (and sptial random fields):
NGroups <- length(Knots)
NGroups


#' If you use:
#' Knots <- seq(1, 10, by = 1)
#' Then you will get a SRF in each year. Takes about 15 minutes on
#' a fast computer.


#* Subsection 14.2: Define projector matrix A----

#' Define the projector matrix. We can use
#' mesh2a for this. No need to change the mesh.
A4.ar1 <- inla.spde.make.A(mesh1,
                           loc = Loc,
                           group = df$YearSeq,
                           n.group = NGroups,  #' <--- number of knots
                           group.mesh = mesh.t)
dim(A4.ar1) #' As before, 2089 observations, and 10 * mesh2a$n nodes.





#* Subsection 14.3: Define spde----

#' Define the SPDE. We use the same priors for the range and sigma_u
#' as in the previous section.
spde4ar1 <- inla.spde2.pcmatern(mesh1,
                                prior.range = c(75, 0.5),
                                prior.sigma = c(2, 0.05))



#* Subsection 14.4: Define SRF w----

#' Define the SRF
#' We now need to define the number of groups.
w4ar1 <- inla.spde.make.index('w',
                               n.spde  = spde4ar1$n.spde,
                               n.group = NGroups)    #' <- Number of knots



#* Subsection 14.5: Stack for the replicate GAM----

#' Sample size
N <- nrow(df)


#' Define the stack.
Stack4.ar1 <- inla.stack(
  tag = "Fit",                       #' Name of stack
  data = list(Counts = df$Counts),   #' Response variable
  A = list(1, 1, A4.ar1),            #' The order matches the order in effects.
  effects = list(
    Intercept = rep(1, N),           #' Use our own intercept
    X         = X[,-1],              #' Dump the intercept from the model.matrix
    w         = w4ar1))              #' Spatial random field


#* Subsection 14.6: Specify the model formula----

#' This is the model with spatial-temporal auto-regression
#' correlation, based on mesh2a.
f4 <-  Counts ~ -1 + Intercept + eqr.std + ppt.std + tmin.std + ws.std + elevation.std +
                     agriculture.std + artificial.std + natural.std +
                     year.std + waterbody_type.L +
                     f(w,
                       model = spde4ar1,
                       group = w.group,
                       control.group = list(model='ar1'))


#* Subsection 14.7: Execute the ar1 spatial-temporal GLM----

#' This is the NB GLM with ar1 spatial-temporal dependency.
#' We are going to use some numerical approximations for faster
#' calculations.
#' If you use this, then it is an option to use it for all earlier
#' models as well.
#

#' Computing time is about 5 minutes. You can also load our workspace;
#' then you don't have to execute the following model. See the course
#' website, or the code in Subsection 2.1.
I4 <- inla(f4,
           family = "nbinomial",
           data = inla.stack.data(Stack4.ar1),
           control.compute = MyControlCompute,
           control.predictor = list(A = inla.stack.A(Stack4.ar1))) #,
#'         Faster calculation:
#'         control.inla = list(strategy = 'gaussian',
#'                             int.strategy = 'eb')) # <--- remove these lines for more precision

# comment that faster calculation part out if we want a more precise model.

summary(I4)
#' Note that the rho is 0.795  . That is a medium value and indicates changes over time. But we already
#' knew that. This comes from "GroupRho for w" in the hyperparameters part


#* Subsection 14.9: Results ar1 spatial-temporal GLM----


#' Numerical output for the regression parameters.
BetasI4 <- I4$summary.fixed[, c("mean", "0.025quant", "0.975quant")]
print(BetasI4, digits = 2)


#' Posterior mean of theta (as in: mu + mu^2 / theta)
theta.pd <- I4$marginals.hyperpar$`size for the nbinomial observations (1/overdispersion)`
theta.pm <- inla.emarginal(function(x) x, theta.pd)
theta.pm


#' Visualise the betas using ggplot2.
#' 1. Create a data frame.
BetasI4_df <- data.frame(Covariate = rownames(BetasI4),
                         Mean      = BetasI4[,"mean"],
                         Lower     = BetasI4[,"0.025quant"],
                         Upper     = BetasI4[, "0.975quant"])

#' 2. Explicitly set the order of the covariates to match their original order
BetasI4_df$Covariate <- factor(BetasI4_df$Covariate, levels = rownames(BetasI4))

#' 3. Plot the info in BetasI1_df.
covariate_labels <- c(
  "Intercept" = "Intercept (Rivers)",
  "eqr.std" = "Ecological Quality Ratio (EQR)",
  "ppt.std" = "Precipitation (mm)",
  "tmin.std" = "Minimum Temperature (°C)",
  "ws.std" = "Wind Speed (m/s)",
  "elevation.std" = "Elevation (m.a.s.l.)",
  "agriculture.std" = "Agricultural Land (%)",
  "artificial.std" = "Artificial Surface (%)",
  "natural.std" = "Natural Area (%)",
  "year.std" = "Year",
  "waterbody_type.L" = "Lakes"
)

# First, let's add a column to identify significant effects
BetasI4_df <- BetasI4_df|>
  mutate(is_significant = !(Lower <= 0 & Upper >= 0))

# Create a function to format the labels based on significance
format_labels <- function(variable, is_significant, labels_map) {
  result <- labels_map[variable]
  result[is_significant] <- sprintf("<b>%s</b>", result[is_significant])
  return(result)
}

# Join significance information with covariate names for labeling
label_data <- BetasI4_df |>
  dplyr::select(Covariate, is_significant) |>
  mutate(label = format_labels(Covariate, is_significant, covariate_labels))

# Modified plot with bold labels for significant covariates
ggplot(BetasI4_df, aes(x = Mean, y = Covariate)) +
  # Add interval lines with varying thickness
  geom_linerange(aes(xmin = Lower, xmax = Upper,
                     linewidth = is_significant,
                     colour = is_significant),
                     show.legend = FALSE) +
  geom_pointrange(aes(xmin = Lower, xmax = Upper,
                      color = is_significant),
                      fatten = 0.5,
                      size = 10,
                      show.legend = FALSE) +
  geom_vline(aes(xintercept = 0),
             linetype = "dashed",
             color = "gray40",
             linewidth = 0.8) +
  labs(x = "Standardized Effect Size",
       y = NULL) +
  theme_minimal() +
  My_theme +
  # Use scale_y_discrete with richtext for HTML formatting
  scale_y_discrete(limits = BetasI4_df$Covariate,
                   labels = setNames(label_data$label, label_data$Covariate)) +
  scale_linewidth_manual(values = c("TRUE" = 1.5, "FALSE" = 1)) +
  scale_color_manual(values = c("gray75", "#21918c")) +
  guides(linewidth = "none") +
  # Add this theme element to render HTML in text
  theme(axis.text.y = ggtext::element_markdown())


#' We want to have 10 panels in 1 ggplot graph. Each panel should
#' show the SRF for that specific year.

#' The posterior mean values can be obtained via:
w.pm <- I4$summary.random$w$mean
length(w.pm)


#' Start a loop.
#'  1. Extract the SRF for year i.
#'  2. Project the SRF on a 100-by-100 grid for that year.
#'  3. Store the results in MyData$w.pm.
#'  4. Determine which of the 100 * 100 rows in MyData
#'     are not on land. Set these MyData$w.pm to NA for plotting
#'     purposes.
#'  5. Store the MyData$w.pm in MyData.All

#' These are the years.
NamesYear  <- levels(as.factor(df$year))
NamesKnots <- Knots


#' Create an object in which we can store the w.pm on a
#' 100-by-100 grid.
MyData.All <- NULL

#' Start a loop to extract the SRF for each knot.
MyData.All <- NULL
for (i in 1:NGroups){
  #' 1. Extract the SRF for year i.
  wpm.i <- w.pm[w4ar1$w.group == i]
  #wsd.i <- wsd[wRepl$w.group == i]

  #' 2. Project the SRF on a 100-by-100 grid for year i.
  Proj     <- fm_evaluator(mesh = mesh1, dims = c(400, 400))
  SRF.Proj <- fm_evaluate(projector = Proj,  field  = wpm.i)

  #'  3. Store the results in MyData$w.pm.
  MyData      <- expand.grid(Xkm = Proj$x, Ykm = Proj$y)
  MyData$w.pm <- as.vector(SRF.Proj)

  #'  4. Determine which of the 100 * 100 rows in MyData
  #'     are not on land. Set these MyData$w.pm to NA for plotting
  #'     purposes.
  MyData.sp <- st_as_sf(MyData,
                        coords = c("Xkm", "Ykm"),
                        crs = "+proj=utm +zone=34 +north +ellps=WGS84 +units=km  +datum=WGS84")

  #' Determine which points in MyData are contained within the study area
  Contained <- st_contains(CroppedLithuania_UTM, MyData.sp)

  #' Extract indices of the points that are on land.
  contained_indices <- unlist(Contained)

  #' Create a logical vector where TRUE indicates the point is on Lithuania
  OnLand <- rep(FALSE, nrow(MyData))
  OnLand[contained_indices] <- TRUE

  #' Not on land: set to NA
  MyData$w.pm[!OnLand] <- NA

  #'  5. Store the MyData$w.pm in MyData.All
  MyData$Year <- NamesYear[i]
  MyData$Knots <- NamesKnots[i]
  MyData.All  <- rbind(MyData.All,
                       MyData)
}

#' Add year as a factor for facetting in ggplot.
MyData.All$fKnots <- factor(MyData.All$Knots)
MyData.All$fYear <- factor(MyData.All$Year)

# add fYear column to match MyData.All
df <- df |>
  mutate(fYear = factor(year))

# Get the range of the spatial smoother values
w_pm_range <- range(MyData.All$w.pm, na.rm = TRUE)

# Create a dummy dataset for the legend (similar to temperature plot)
legend_data <- data.frame(
  x = seq(w_pm_range[1], w_pm_range[2], length.out = 100),
  y = rep(1, 100)
)

# make bounding box in UTM to match panel plot figure
bb_utm <- st_bbox(st_transform(st_as_sfc(bb), crs = "+proj=utm +zone=34 +north +ellps=WGS84 +units=km +datum=WGS84"))

# Create a dataset for year labels
# Using the known bounding box and placing labels in the upper right
# Adjust these values to fine-tune the position
years_data <- MyData.All |>
  group_by(fYear) |>
  summarise(
    # Calculate panel extents (in case they vary by year)
    x_min = min(Xkm, na.rm = TRUE),
    x_max = max(Xkm, na.rm = TRUE),
    y_min = min(Ykm, na.rm = TRUE),
    y_max = max(Ykm, na.rm = TRUE)
  ) |>
  mutate(
    # Position in the upper right area
    # Placing at approximately 85% across and 85% up from bottom
    x_pos = x_min + (x_max - x_min) * 0.8,
    y_pos = y_min + (y_max - y_min) * 0.8
  )

# And plot everything in one ggplot with internal year labels
ggplot() +
  geom_sf(data = CroppedLithuania_UTM,
          fill = NA,
          color = "black",
          linewidth = 0.5) +
  # Set the extent using coord_sf with the UTM bounding box
  coord_sf(xlim = c(bb_utm["xmin"], bb_utm["xmax"]),
           ylim = c(bb_utm["ymin"], bb_utm["ymax"])) +
  # W's
  geom_raster(data = MyData.All,
              aes(x = Xkm,
                  y = Ykm,
                  fill = w.pm)) +
  # Add a hidden continuous scale just for the legend
  geom_point(
    data = legend_data,
    aes(x = x, y = y, color = x),
    alpha = 0
  ) +
  # Add year labels inside each panel - positioned lower right
  geom_text(
    data = years_data,
    aes(x = x_pos, y = y_pos, label = fYear),
    size = 5,
    fontface = "bold",
    hjust = 0.5,  # Center horizontally
    vjust = 0.5,  # Center vertically
    color = "black",
    # Add a white background/outline to make text more readable
  ) +
  # Color scheme for the raster fill
  scale_fill_gradient2(
    name = "Posterior mean values (w)",
    midpoint = 0,
    low = "#21918c",
    mid = "white",
    high = "#440154",
    na.value = NA,
    guide = "none"  # Hide this legend as we'll use the color legend
  ) +
  # Color scale for the legend, using the same Color scheme for the raster fill
  scale_color_gradient2(
    name = "Posterior mean values (w)",
    midpoint = 0,
    low = "#21918c",
    mid = "white",
    high = "#440154",
    # Coarser breaks with increments of 1
    breaks = seq(floor(w_pm_range[1]), ceiling(w_pm_range[2]), by = 1),
    guide = guide_colorbar(
      barwidth = 20,
      barheight = 0.5,
      title.position = "top",
      title.hjust = 0.5,
      label.hjust = 0.5
    )
  ) +
  # Labels
  labs(x = "Longitude", y = "Latitude",
       caption = "Mesh size (w) = 2871 • Prior range = c(75, 0.5) • Sigma_u = c(2, 0.05)") +
  # Themes
  theme_minimal() +
  My_theme +
  # Additional theme adjustments
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    legend.box.spacing = unit(0.5, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.3),  # Smaller tick marks
    axis.ticks.length = unit(0.1, "cm"),  # Shorter tick length
    strip.text = element_blank(),  # Remove facet strip text since we're adding internal labels
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # Tick text size 10
    axis.text.y = element_text(size = 10),  # Tick text size 10
    axis.title.x = element_text(size = 16),  # Axis label size 16
    axis.title.y = element_text(size = 16),   # Axis label size 16
    # Instead of legend.position = "bottom", use:
    legend.position = c(0.69, 0.1),  # x and y coordinates (adjust as needed)
    legend.justification = c(0.5, 0.5),  # Center the legend at the position
    legend.direction = "horizontal",
    legend.box = "horizontal",
    plot.caption = element_text(size = 6, hjust = 1, margin = margin(t = 10))
  ) +
  # Faceting without labels (since we're adding them inside)
  facet_wrap(~ fYear, ncol = 3)

#' From year-to-year the SRF is changing slowly. That is because the rho is relatively high.
summary(I4)



# Section 15: Quick comparison of all models----

#' Compare DIC and WAIC.
DICs <- c(I1$dic$dic, I1$dic$dic, I2$dic$dic, I3$dic$dic, I4$dic$dic)
WAICs <- c(I1$waic$waic, I1$waic$waic, I2$waic$waic, I3$waic$waic, I4$waic$waic)
Results <- data.frame(Models = c("NB GLM",
                                 "NB GLM + SRF with mesh2a",
                                 "NB GLM + SRF with mesh2b",
                                 "NB GLM + replicate SRF",
                                 "NB GLM + ar1 SRF"),
                      DIC = DICs,
                      WAIC = WAICs)
Results
#' The NB GLM is the best!



# Section 16: Model validation ar1 spatial-temporal GLM----

#* Subsection 16.1: Getting scaled quantile residuals----


#' Do the posterior simulation of regression parameters 1000 times.
SimData <- inla.posterior.sample(n = 1000, I4)

#' Here is some fancy code to grab the positions of specific variables.
#' Custom function.
MyGrep <- function(x, SimulatedData){
  # SimulatedData is the name of the object containing the simulation results
  # x is an element of BetasInModel
  names.SimData <- attributes(SimulatedData[[1]]$latent)$dimnames[[1]]
  names.SimData[grep(x, names.SimData)]
}


#' What are the names of the regression parameters?
BetasInModel <- rownames(I4$summary.fixed)

#' Get their location in the object with simulation results.
BetaNames    <- unlist(lapply(BetasInModel,
                              FUN = MyGrep,
                              SimulatedData = SimData))
MyID     <- function(x){ which(rownames(SimData[[1]]$latent) == x) }
RowsBeta <- lapply(BetaNames, MyID)
RowsBeta <- as.numeric(RowsBeta)
RowsBeta



#' Now the spatial correlation random effects.
#' Determine the names of the random effects.
#' In the code below, we are using:  ^w:\d+$
#' This will only match strings that start with "w:", followed
#' exclusively by one or more digits, and then end.
names.SimData <- attributes(SimData[[1]]$latent)$dimnames[[1]]
names.SimData
wNames <- names.SimData[grep("^w:\\d+$", names.SimData)]

#' This is probably easier: wNames <- paste("w",1:mesh3$n, sep = ":")
#' Determine on which rows the random effects are.
RowsW <- lapply(wNames, MyID)
RowsW <- as.numeric(RowsW)
RowsW



#' Calculate 1000 times the fitted values via: mu = exp(X * beta + A * w).
#' Get the design matrix X.
X <- model.matrix(~ eqr.std + ppt.std + tmin.std + ws.std + elevation.std +
                    agriculture.std + artificial.std + natural.std +
                    year.std + waterbody_type, data = df)
X <- as.matrix(X)


#' Start a loop. In each iteration:
#'   1. Calculate the fitted values mu.
#'   2. Simulate NB data with the mean equal to mu.
N    <- nrow(df)                      #' Sample size.
Ysim <- matrix(nrow = N, ncol = 1000) #' Create space.

#' Start the loop (can be done more efficient with lapply).
for (i in 1:1000){
  Betas       <- SimData[[i]]$latent[RowsBeta]
  wk          <- SimData[[i]]$latent[RowsW]
  FixedPart   <- X %*% Betas
  SpatialPart <- A4.ar1 %*% wk            #' Note the A4.ar1.                         #' Will change later.
  mu  <- exp(FixedPart + SpatialPart)[,1]
  Ysim[,i] <- rnegbin(n = N,               #' Simulated count data
                      mu = mu,
                      theta = theta.pm)
}

par(mfrow = c(2,2))
hist(df$Counts, main ="Observed Vector abundance")
hist(Ysim[,1],  main = "First simulated data set")
hist(Ysim[,2],  main = "Second simulated data set")
hist(Ysim[,3],  main = "Third simulated data set")
par(mfrow = c(1,1))
#' The observed data seems to be slightly different?


#' Now we have 1000 simulated data sets. Give them all to
#' DHARMa, and it will calculate scaled quantile residuals.
N <- nrow(df)
df$Fit <- I4$summary.fitted.values[1:N, "mean"]  #' Fitted values.

E5.sqr <- createDHARMa(simulatedResponse = Ysim,
                       observedResponse = df$Counts,
                       fittedPredictedResponse = df$Fit,
                       integerResponse = TRUE)
#' Now we have scaled quantile residuals.

plot(df$Counts, df$Fit)

#* Subsection 16.2: Check for overdispersion and zero-inflation----

testDispersion(E5.sqr)
#' A bit of underdispersion?
#' Try the Poisson version of the model.

testZeroInflation(E5.sqr)
#' Fine

#* Subsection 16.3: Check for homogeneity of variance----

#' Plot the scaled quantile residuals versus (ranked) fitted values.
plotResiduals(E5.sqr, quantreg = TRUE, smoothScatter = FALSE)
#' Tolerable



#* Subsection 16.4: Check for uniformity of the residuals----

#' In DHARMa, we verify whether the scaled quantile residuals are
#' uniform distributed.
par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(E5.sqr, testUniformity = TRUE,
           testOutliers = TRUE, testDispersion = FALSE)
#' OK!



#* Subsection 16.5: Plot residuals versus the covariates----

#' Plot the scaled quantile residuals versus each covariate
#' in the model.
plotResiduals(E5.sqr, form = df$eqr, quantreg = TRUE)            #' Tolerable
plotResiduals(E5.sqr, form = df$ppt, quantreg = TRUE)            #' Fine
plotResiduals(E5.sqr, form = df$tmin, quantreg = TRUE)           #' Fine
plotResiduals(E5.sqr, form = df$ws, quantreg = TRUE)             #' Fine
plotResiduals(E5.sqr, form = df$elevation, quantreg = TRUE)      #' Tolerable
plotResiduals(E5.sqr, form = df$agriculture, quantreg = TRUE)    #' Fine
plotResiduals(E5.sqr, form = df$artificial, quantreg = TRUE)     #' Fine
plotResiduals(E5.sqr, form = df$natural, quantreg = TRUE)        #' Fine
plotResiduals(E5.sqr, form = df$year, quantreg = TRUE)           #' Fine
plotResiduals(E5.sqr, form = df$year, asFactor = TRUE, quantreg = TRUE)           #' Fine
plotResiduals(E5.sqr, form = df$waterbody_type, asFactor = TRUE, quantreg = TRUE) #' Fine

E5.sqr$scaledResiduals
# actual residuals that we can use for the gam trick
# grab it, then make numeric, then do the rest.
# loop that runs the plot residuals for each covariate and extracted the residuals and the p-value from the test
# store everything under each other and then plot it via ggplot
# instead of 11 pictures, have everything in one ggplot and colour everything by p-value
# can present this as an online supplement but only if the referee asks / to anticipate what a reviewer will ask for



#* Subsection 16.6: Check spatial dependency----

#' Combine residuals and coordinates. We also added year to this
#' data frame. We will explain in a moment why.
Data4Vario <- data.frame(Res   = E5.sqr$scaledResiduals,
                         Xkm   = df$Xkm,
                         Ykm   = df$Ykm,
                         fYear = as.factor(df$year))

#' Convert to sf object.
MyData_sf <- st_as_sf(x = Data4Vario,
                      coords = c("Xkm", "Ykm"),
                      crs = NA)  #' Non-Cartesian coordinates.

#' Calculate the variogram.
V4 <- variogram(Res ~ 1,
                data = MyData_sf,
                cutoff = 100,
                cressie = TRUE)

#' Plot the variogram.
ggplot() +
  geom_point(data = V4,
             aes(x = dist, y = gamma)) +
  geom_smooth(data = V4,
              aes(x = dist, y = gamma),
              se = TRUE,
              span = 1,
              alpha = 0.3) +
  labs(x = "Distance (in km)", y = "Semi-variogram") +
  theme(text = element_text(size = 15))
#' There is still some dependency in the residuals. However, if we consider the y-axis scale, then the dependency is not strong.
#' We can conclude that the model is fine.



#' Well...not entire correct. What is this variogram telling
#' is? We gave it residuals from all years in one long
#' vector plus the spatial coordinates. It does not know
#' that the residuals come from 10 blocks. We should also
#' calculate a variogram for the residuals of each year.



#' That means 10 variograms.
#' Instead of doing this manually, we wrote some fancy code to do this
#' automatically. There is no need to fully understand this, but we
#' did annotate the R code a little bit.

# Step 1: Split the MyData_sf object by fYear.
SplitData <- split(MyData_sf, MyData_sf$fYear)
SplitData

#' Step 2: lapply()
#' Use lapply() to execute the function CalculateVariogram() for each
#' element (i.e. residuals from a specific year) in SplitData.
CalculateVariogram <- function(year) {
  v <- variogram(Res ~ 1,
                 data <- SplitData[[year]],
                 cutoff = 100,
                 cressie = TRUE)
  v$fYear <- as.numeric(year)
  return(v)}

#' For each element in SplitData (i.e. for each year), calculate the variogram.
VarioList <- lapply(names(SplitData), CalculateVariogram)
VarioList
#' np:    number of site combinations within a distance band.
#' dist:  Distance band for which the variogram is being calculated.
#' gamma: Variogram values.
#' fYear: Year identifier.
#' The rest is not relevant.

#' Step 3: unlist
#' VarioList is a list with results. We want to have it as a data frame.
V_all <- do.call(what = rbind,
                 args = VarioList)

#' This function is stacking the elements of VarioList on top of
#' each other to create a single data frame.
#' The function do.call calls the rbind function with as argument VarioList.
#' Do something (what) with something (args).


#' Now we can plot the results.
ggplot() +
  geom_point(data = V_all,
             aes(x = dist, y = gamma)) +
  geom_smooth(data = V_all,
              aes(x = dist, y = gamma),
              se = TRUE,
              span = 1,
              alpha = 0.3) +
  labs(x = "Distance (in km)", y = "Semi-variogram") +
  facet_wrap(~ fYear, scales = "fixed") +
  theme(text = element_text(size = 15))
#' The variograms of the residuals by year do certainly not show spatial
#' dependency.

#' We conclude that there is no strong spatial dependency in the residuals.
#' To be more precise:
#'  -The graph on the right shows that there is no spatial dependency within
#'   a year.
#'  -The previous graph showed that if we look at all residuals (irrespective
#'   of year) there is neither spatial dependency.

# 2015, 2018, 2020, 2021, 2022
# try multiple likelihood
# Or use the AR1 term per year, not with knots.


# Section 17: Detailed comparison of all models----

#' Let's plot the results of the model, without and with the
#' spatial correlation side by side.
Out2                     <- I0$summary.fixed[,c("mean", "0.025quant", "0.975quant")]
Out2.mesh1               <- I1$summary.fixed[,c("mean", "0.025quant", "0.975quant")]
Out2.mesh2               <- I2$summary.fixed[,c("mean", "0.025quant", "0.975quant")]
Out2.mesh2repl           <- I3$summary.fixed[,c("mean", "0.025quant", "0.975quant")]
Out2.mesh2ar1            <- I4$summary.fixed[,c("mean", "0.025quant", "0.975quant")]
rownames(Out2)           <- I0$names.fixed
rownames(Out2.mesh1)     <- I1$names.fixed
rownames(Out2.mesh2)     <- I2$names.fixed
rownames(Out2.mesh2repl) <- I3$names.fixed
rownames(Out2.mesh2ar1)  <- I4$names.fixed

#' Names for the models:
MyNames <- c("NB GLM",
             "Spatial NB GLM mesh 1",
             "Spatial NB GLM mesh 2",
             "Spatial NB GLM mesh 2 replicates",
             "Spatial NB GLM mesh 2 ar1")

#' Compare results of the two models using MyCompareBetasofModels
#' (which is in our support file).
MyCompareBetasofModels(AllModels = list(Out2, Out2.mesh1, Out2.mesh2, Out2.mesh2repl, Out2.mesh2ar1),
                       ModelNames = MyNames)

ggsave("Plots/FigureS7_fixed _effect_model_comparisons.png", width = 10, height = 10, units = "in", bg = "white", dpi = 300)

#' There are some differences between the model without and with spatial dependency.
#' There are no major differences between the spatial models.


#' if the theta is ultra small, then you should be worried
I0$summary.hyperpar[1, "mean"]
I1$summary.hyperpar[1, "mean"]
I2$summary.hyperpar[1, "mean"]
I3$summary.hyperpar[1, "mean"]
I4$summary.hyperpar[1, "mean"]



# Section 18: Sensitivity analysis a la Dambly et al. (2023)----

#' We may, or we may not discuss this section during the course.
#' That depends on the available time and questions.
#' Important:
#'    Do not focus on the R code in this section!
#'    It contains a lot of loops. Try to understand
#'    the graphs!


#' We wrote this section after the publication of this paper.
#' https://onlinelibrary.wiley.com/doi/full/10.1111/ecog.06391

#' Before everyone working with INLA starts to panic, we are not too
#' impressed by this paper. Its conclusions are that mesh configuration
#' can be rather influential when using INLA. That is partly true,
#' but the paper uses a set of Mickey Mouse meshes. And there may be
#' some data specific issues as well.

#' But it does not harm to try a variety of meshes and see whether
#' that influences the results. Then you can write in your paper
#' that you did this (no need to show all these graphs).



#' In this section, we will:
#'  - Define 6  meshes with different resolution.
#'  - Use 5 different PC priors for the Range.

#' This gives 6 * 5 = 30 combinations, and therefore 30 INLA models.
#' Computing time is about 10 minutes.


#' In general, you can save a lot of computing time by using a
#' non-convex boundary, particularly when your study area is
#' non-square or non-rectangular. Therefore, we will use a
#' non-convex boundary for all models.

NonConBoundary <- fm_nonconvex_hull(Loc, convex = -0.08)


#* Subsection 18.1: Defining 6 meshes with different resolution----

#' We will consider 6 meshes with different resolutions. Recall
#' that the resolution of the mesh was defined via the construction:
#'  RangeGuess <- 50 <---- based on variogram of non-spatial GLM
#'  MaxEdge    <- RangeGuess / 5

#' Choose a couple values below and above the range guess values.
#' Instead of the 50, we will now use:
#'  Mesh 1: RangeGuess = 10
#'  Mesh 2: RangeGuess = 20
#'  Mesh 3: RangeGuess = 30
#'  Mesh 4: RangeGuess = 50 <--- orignal value
#'  Mesh 5: RangeGuess = 70
#'  Mesh 6: RangeGuess = 100


#' Define the vector of RangeGuess values:
RangeGuessValues <- c(10, 20, 30, 50, 70, 100)
NumberOfMeshes   <- length(RangeGuessValues)
NumberOfMeshes




#* Subsection 18.2: Defining 5 priors for the range----

#' In the previous subsection, we defined 6 meshes with different
#' resolutions. Now we will define 5 different priors for the range.


#' Recall that we used:
#'   P(Range < 75) = 0.5.



#' Now, we will use:
#'  P(Range < 10) = 0.5  #' Unlikely that the range is < than 10 km. <--- allow for small
#'  P(Range < 25) = 0.5  #' Unlikely that the range is < than 25 km.
#'  P(Range < 50) = 0.5  #' Unlikely that the range is < than 50 km.
#'  P(Range < 75) = 0.5  #' Unlikely that the range is < than 75 km. <--- original value
#'  P(Range < 100) = 0.5  #' Unlikely that the range is < than 100 km. <--- allow for big

RangePriorValues  <- c(10, 25, 50, 75, 100)
NumOfPriors4Range <- length(RangePriorValues)

#' VERY IMPORTANT
#' DONT FORGET TO TRY DIFFERENT SIGMA Us - 3 different sigmas
#' VERY IMPORTANT

#* Subsection 18.3: Create the 6 meshes----

#' We make 6 meshes and put them in a list.

# recall our original mesh
# mesh1 <- fm_mesh_2d(loc = Loc,
#                     max.edge = c(1, 5) * MaxEdge,
#                     cutoff = MaxEdge / 5)

#' Initialize an empty list to store the meshes.
MeshList <- list()

#' Loop over each RangeGuessValues value.
for (i in 1:NumberOfMeshes) {
  RangeGuess <- RangeGuessValues[i]
  MaxEdge    <- RangeGuess / 5

  #' Define mesh.
  Mesh <- fm_mesh_2d(loc = Loc,
                     max.edge = c(1, 5) * MaxEdge,
                     cutoff = MaxEdge / 5)

  #' Store the mesh in the list.
  MeshList[[i]] <- Mesh
}
#' MeshList now contains all the meshes.



#' How many w's (nodes) have these meshes?
#' Loop through each mesh in the list.
MeshSizes <- vector(length = NumberOfMeshes)
Counter <- 1
for (mesh in MeshList) {
  MeshSizes[Counter] <- mesh$n
  Counter <- Counter + 1
}
MeshSizes
#' The first mesh has 66555 w's, the 6th has 1196 w's.



#' Plot the meshes.
plotList <- list()
for (i in 1:NumberOfMeshes) {
  mesh <- MeshList[[i]]
  fm_crs(mesh) <- st_crs(CroppedLithuania_UTM)

  p <- ggplot() +
    theme_minimal() +
    labs(title = paste("Mesh ", i, " (",MeshSizes[i],")", sep = "")) +
    geom_fm(data = mesh) +
    #geom_point(data = df, aes(x = Xkm, y = Ykm), alpha = 0.5) +
    geom_sf(data = CroppedLithuania_UTM, fill = "transparent", col = "red")

  plotList[[i]] <- p
}

#' And plot all meshes.
CombinedPlot <- plot_grid(plotlist = plotList,
                          ncol = 3)
CombinedPlot


ggsave(CombinedPlot,
       filename = "Sensitivity/FigureS2_sensisitivity_mesh_comparison.png",
       width = 14,
       height = 7,
       bg = "white",
       dpi = 300)

#* Subsection 18.4: Define the SPDE----

#' We need to specify priors for the two Matern correlation
#' parameters Range and sigma.

#' In the analysis above, we used:
#'  P(Range < 75) = 0.5
#'  P(sigma > 2) = 0.05


#' In this sensitivity analysis, we will use:
#'  P(Range < 10) = 0.5   #' likely that the range is < than 10
#'  P(Range < 25) = 0.5   #' likely that the range is < than 25
#'  P(Range < 50) = 0.5   #' likely that the range is < than 50
#'  P(Range < 75) = 0.5   #' likely that the range is < than 75
#'  P(Range < 100) = 0.5  #' likely that the range is < than 100

#' We will keep the prior for the sigma as it is. Feel
#' free to change that too!


#' The spde depends on mesh and prior. Hence, there will be
#' 6 * 5 unique spde's.

#' Initialize an empty list to store the SPDE models
spdeList <- list()

#' Loop over each mesh and each value in RangePriorValues.
#' Put each SPDE in a list.
Counter <- 1
for (i in 1:NumberOfMeshes) {
  mesh <- MeshList[[i]]
  for (j in 1:NumOfPriors4Range) {
    rangePrior <- RangePriorValues[j]

    #' Create the SPDE model using the current
    #' mesh and prior for the range.
    spde <- inla.spde2.pcmatern(mesh,
                                prior.range = c(rangePrior, 0.5),
                                prior.sigma = c(2, 0.05))

    #' Store the SPDE model in the list.
    spdeList[[Counter]] <- spde
    Counter <- Counter + 1
  }}

#' spdeList now contains the 30 SPDE models.


#* Subsection 18.5: Execute INLA 30x----

#' Make the X matrix using model.matrix()
X <- model.matrix(~ eqr.std + ppt.std + tmin.std + ws.std + elevation.std +
                    agriculture.std + artificial.std + natural.std +
                    year.std + waterbody_type, data = df)
X <- as.data.frame(X) #' Avoids an error message in INLA

#' In INLA, you can't use ':' in the variable names. Note that the
#' interaction term uses : in its name. Replace the : by '_'.
# OldColnames  <- colnames(X)
# NewColNames  <- gsub(pattern = ":", replacement = "_", x = OldColnames)
# colnames(X)  <- NewColNames


#' Define sample size.
N <- nrow(df)


#' Initialize a list to store the INLA results.
InlaResultsList <- list()

#' Loop over each mesh (6x).
Counter <- 1
for (i in 1:NumberOfMeshes) {
  print(i)
  mesh <- MeshList[[i]]
  A    <- inla.spde.make.A(mesh, loc = Loc)

  #' Loop over each PC prior range value (5x).
  for (j in 1:NumOfPriors4Range) {
    spde       <- spdeList[[Counter]]
    w.index <- inla.spde.make.index(name = 'w', n.spde = spde$n.spde)

    #' Make a stack.
    Stack <- inla.stack(tag = "Fit",
                        data = list(vec_abund = df$vec_abund),
                        A = list(1, 1, A),
                        effects = list(Intercept = rep(1, N),
                                       X = X[,-1],
                                       w = w.index))

    #' Specify the model formula
    f.mesh <- vec_abund ~ -1 + Intercept + eqr.std + ppt.std + tmin.std + ws.std + elevation.std +
                     agriculture.std + artificial.std + natural.std +
                     year.std + waterbody_type.L + f(w, model = spde)

    #' Execute the INLA model.
    Inla.Results <- inla(f.mesh,
                         family = "nbinomial",
                         data = inla.stack.data(Stack),
                         control.compute = list(dic = TRUE, waic = TRUE),
                         control.predictor = list(A = inla.stack.A(Stack)))

    #' Store the model output in the list.
    InlaResultsList[[Counter]] <- Inla.Results
    Counter <- Counter + 1
  }
}
#' Now we can apply a sensitivity analysis and see how
#' much "things" are changing for the different meshes and
#' priors.




#* Subsection 18.6: Results; Compare DICs----

#' We will plot the DICs of all models in a heatmap.


#' Initialize a vector to store the DIC values
dicValues <- vector(length = NumberOfMeshes * NumOfPriors4Range)

#' Loop over the INLA results and extract DIC values.
for (i in 1:length(InlaResultsList)) {
  modelResult <- InlaResultsList[[i]]

  #' Extract the DIC value and store it.
  dicValues[i] <- modelResult$dic$dic
}
#' dicValues now contains the DIC values from each model.


#' Make a heatmap of these DIC values.
dicMatrix <- matrix(dicValues,
                    nrow = NumberOfMeshes,
                    ncol = NumOfPriors4Range)

rownames(dicMatrix) <- paste("w", MeshSizes, sep = "=")
colnames(dicMatrix) <- paste(" P(Range<", RangePriorValues, ")=0.5", sep = "")
dicMatrix


#' w=2871:
#'  -This is the mesh with 2871 ws.
#'  -We obtained it by using RangeGuess = 75

#' P(Range<75)
#'  - We used P(Range < 75) = 0.05
#'  - It is unlikely that the range is smaller than 75.


#' In other words:
#'   Top rows have lots of ws
#'   Left columns: We allow the range to be small...if it wants to.
#'   Right columns: We do not allow the range to be small.
#'   Bottom right: Coarse mesh and range is most likely large


#' Visualise the results.
#' Convert the matrix to a long format.
dicLong <- melt(dicMatrix)

#' Rename the columns for clarity.
colnames(dicLong) <- c("Mesh", "Prior", "DIC")

ggplot(dicLong,
       aes(x = Prior,
           y = Mesh,
           fill = DIC)) +
  geom_tile() +
  scale_fill_gradient(low = "#21918c", high = "#440154") +
  theme_minimal() +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1)) +
  labs(title = "DIC Values Heatmap",
       x = "Prior value for range",
       y = "Mesh resolution",
       fill = "DIC",
       caption = "Sigma_u = c(2, 0.05)")

ggsave("Sensitivity/FigureS3.1_sensitivity_DIC_values.png", width = 10, height = 10, dpi = 300, bg = "white")

#' This is what we used in the actual analysis:
mesh1$n
#' Conclusions:
#'   -As long as we use a sensible mesh configuration, and a
#'    sensible prior, results do not seem to be too different.
#'   -If we choose odd mesh configurations or priors, then
#'    we have higher DICs.
#'  -Within the same mesh, prior choice matters...just pick a sensible one.
#'  -For a fine mesh, choice of prior seems to be less important.
#'  -DIC vlue is very similar regardless of prior choice or mesh resolution.


#* Subsection 18.7: Compare R^2 values----

#' Not sure if this makes sense for a GLM. R^2 is not the best tool to
#' compare observed and fitted data for a NB model.


#* Subsection 18.8: Compare imposed dependency----

#' Show how the posterior mean values of the range change
#' for the different meshes and priors.

#' Extract the estimated ranges and Kappa values for each model.
#' Initialize a matrix to store the DIC values.
EstimatedRanges <- matrix(nrow = NumberOfMeshes, ncol = NumOfPriors4Range)
EstimatedKappas <- matrix(nrow = NumberOfMeshes, ncol = NumOfPriors4Range)

#' Loop over each mesh and prior.
k <- 1
for (i in 1:NumberOfMeshes) {
  for (j in 1:NumOfPriors4Range) {
    modelResult <- InlaResultsList[[k]]
    SpatialParams <- MySpatialParams(Model = modelResult,
                                     ThisSpde = spdeList[[k]])
    EstimatedRanges[i,j] <- SpatialParams[3]
    EstimatedKappas[i,j] <- SpatialParams[1]
    k <- k + 1
  }}

rownames(EstimatedRanges) <- paste("w=", MeshSizes, sep = "")
colnames(EstimatedRanges) <- paste("P(Range<", RangePriorValues, ")=0.5", sep = "")
EstimatedRanges


#' Visualise the results
#' Convert the matrix to a long format.
EstimatedRangesLong <- melt(EstimatedRanges)


#' Rename the columns for clarity.
colnames(EstimatedRangesLong) <- c("Mesh", "Prior", "Ranges")


#' Plot the results.
ggplot(EstimatedRangesLong,
       aes(x = Prior,
           y = Mesh,
           fill = Ranges)) +
  geom_tile() +
  scale_fill_gradient(low = "#21918c", high = "#440154") +
  theme_minimal() +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1)) +
  labs(title = "Posterior mean values of the range",
       x = "Prior Range Value",
       y = "Mesh Resolution",
       fill = "P(Range < ..) = 0.5",
       caption = "Sigma_u = c(2, 0.05)")

ggsave("Sensitivity/FigureS3.2_sensitivity_range_estimates.png", width = 10, height = 10, dpi = 300, bg = "white")

#' Conclusions:
#' -For any mesh resolution, choice of prior is influential, no matter the coarseness of the mesh.
#' -If we choose a reasonable prior and a reasonable mesh, the range should be reasonable.



#' Visualise the Matern correlation for the different meshes
#' and prior range values.

#' Calculate the distance between sites.
D <- as.matrix(dist(Loc))
d.vec <- seq(0, max(D), length = 100)

#' Initialize an empty data frame to store the
#' results for plotting
plotData <- data.frame(Distance = d.vec)

#' Add nice column and row labels.
rownames(EstimatedKappas) <- paste("w=", MeshSizes, sep = "")
colnames(EstimatedKappas) <- paste("P(Range<", RangePriorValues, ")", sep = "")


#' Loop over each kappa and calculate the Matern correlation.
MyData <- NULL
for (i in 1:nrow(EstimatedKappas)) {
  for (j in 1:ncol(EstimatedKappas)) {
    kappa <- EstimatedKappas[i, j]

    Cor.M <- kappa * d.vec * besselK(kappa * d.vec, 1)
    Cor.M[1] <- 1

    # Store the results in the data frame
    Data.i <- data.frame(Distance   = d.vec,
                         MaternC    = Cor.M,
                         PriorRange = paste("P(Range<", RangePriorValues[j],")", sep = ""),
                         MeshName   = paste("w=",MeshSizes[i], sep =""))
    MyData <- rbind(MyData,
                    Data.i)
  }
}

#' For ggplot2 plotting purposes:
MyData$PriorRange <- factor(MyData$PriorRange,levels = unique(MyData$PriorRange))
MyData$MeshName   <- factor(MyData$MeshName,levels = unique(MyData$MeshName))


#' Plot results with facet_grid.
ggplot(MyData,
       aes(x = Distance,
           y = MaternC,
           group = PriorRange,
           col = PriorRange)) +
  geom_line() +
  facet_wrap( ~ MeshName, ncol = 3) +
  theme_minimal() +
  xlim(0,150)+
  labs(title = "Matern correlation",
       x = "Distance",
       y = "Correlation",
       color = "Prior",
       caption = "Sigma_u = c(2, 0.05)")

ggsave("Sensitivity/FigureS3.3_sensitivity_MaternCorrelation1.png", width = 10, height = 10, dpi = 300, bg = "white")

#' This is in principle the same information as we discussed for the posterior mean values of the range. Just go back one graph.
#'  -For finer meshes, the ranges differences are not so different, and the lines are smoother. But not by much.
#'  -For coarser meshes, the range is larger, but the convergence of the lines is better. But again, not by much.



#' And flip it.
ggplot(MyData,
       aes(x = Distance,
           y = MaternC,
           group = MeshName,
           col = MeshName)) +
  geom_line() +
  facet_wrap( ~ PriorRange, ncol = 3) +
  theme_minimal() +
  xlim(0,150)+
  labs(title = "Matern correlation",
       x = "Distance",
       y = "Correlation",
       color = "mesh resolution",
       caption = "Sigma_u = c(2, 0.05)")

ggsave("Sensitivity/FigureS3.4_sensitivity_MaternCorrelation2.png", width = 10, height = 10, dpi = 300, bg = "white")

#' Same as in the previous two pictures.
#'  -No so much difference in the Matern correlation for different meshes.
#'  -Bigger differences if we select a Micky-Mouse prior.


#* Subsection 18.9: Regression parameters----

#' We now show how/if the regression parameters are changing
#' for different mesh configurations and priors for the range.

#' Pick one of the following covariates:
ThisX <- "eqr.std"
# ThisX <- "ppt.std"
# ThisX <- "tmin.std"
# ThisX <- "ws.std"
# ThisX <- "elevation.std"
# ThisX <- "agriculture.std"
# ThisX <- "artificial.std"
# ThisX <- "natural.std"
# ThisX <- "year.std"
# ThisX <- "waterbody_type.L"

#' Grab the posterior mean value and the 95% credible interval
#' for the selected covariate, for each mesh and prior combi.
MyData <- NULL
k <- 1
for (i in 1:NumberOfMeshes) {
  for (j in 1:NumOfPriors4Range) {
    modelResult <- InlaResultsList[[k]]

    Betas <- modelResult$summary.fixed[, c("mean",
                                           "0.025quant",
                                           "0.975quant")]

    # Store the results in the data frame
    Data.i <- data.frame(beta.pm    = Betas[ThisX, "mean"],
                         SeUp       = Betas[ThisX, "0.975quant"],
                         SeLo       = Betas[ThisX, "0.025quant"],
                         PriorRange = RangePriorValues[j],
                         MeshName   = MeshSizes[i])
    MyData <- rbind(MyData,
                    Data.i)
    k <- k + 1
  }}


#' Add some nice labels for ggplot.
MyData$PriorRange <- factor(MyData$PriorRange,
                            levels = RangePriorValues,
                            labels = paste("P(Range<",
                                           RangePriorValues, ")",
                                           sep = ""))

MyData$MeshName   <- factor(MyData$MeshName,
                            levels = MeshSizes,
                            labels = paste("w=",
                                           MeshSizes,
                                           sep = ""))


#'  Plot the results.
ggplot(data = MyData,
       aes(x = MeshName,
           y = beta.pm,
           ymax = SeUp,
           ymin = SeLo)) +
  geom_point() +
  geom_errorbar(alpha = 0.5) +
  labs(x = "Mesh size",
       y = ThisX,
       caption = "Sigma_u = c(2, 0.05)") +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1)) +
  facet_wrap(~ PriorRange, ncol = 3)


eqr
ppt
tmin
elevation

# Combine the plots in one row with three columns
combined_plot_sensitivity <- eqr + elevation + ppt + tmin +
  patchwork::plot_layout(ncol = 2) +
  patchwork::plot_annotation(tag_levels = 'A')  # This adds a, b, c labels

# Display the combined plot
combined_plot_sensitivity


ggsave("Sensitivity/FigureS5_sensitivity_RegressionParameters.png", width = 20, height = 20, dpi = 300, bg = "white")

#' Conclusions:
#'  - it seems range prior does not have much influence on the regression parameters. The mesh size does not seem to have much influence either.
#'  - The regression parameters are very similar for all models. This shows that our model is doing a good job in estimating the regression parameters.


#* Subsection 18.10: Compare the spatial random fields----

#' We now plot the 30 spatial random fields side by side.
#' Given what we have seen so far, we do not expect much
#' differences.



#' Loop over each mesh and prior range values.
#' Calculate the spatial random field for each combi, remove
#' the points in the sea and in Northern Ireland, and
#' plot the results.

MyData <- NULL
k <- 1
for (i in 1:NumberOfMeshes) {
  #' Grab mesh number i.
  mesh <- MeshList[[i]]

  for (j in 1:NumOfPriors4Range) {
    #' Grab the posterior mean values of the w.
    modelResult <- InlaResultsList[[k]]
    w.pm <- modelResult$summary.random$w$mean

    #' Project the w on a 400-by-400 grid.
    Proj <- fm_evaluator(mesh = mesh,  dims = c(400, 400))
    SRF.Proj <- fm_evaluate(projector = Proj, field  = w.pm)

    #' Extract and combine the relevant information.
    Data.i      <- expand.grid(Xkm = Proj$x, Ykm = Proj$y)
    Data.i$w.pm <- as.vector(SRF.Proj)

    #' For all sites not inside Ireland, set the w.pm to NA.
    #' Create an sf object from the coordinates
    sputm <- st_as_sf(Data.i,
                      coords = c("Xkm", "Ykm"),
                      crs = "+proj=utm +zone=34 +north +ellps=WGS84 +units=km  +datum=WGS84")

    #' Determine which points in Data.i are contained within Ireland_UTM
    InLithuania <- st_contains(CroppedLithuania_UTM, sputm)

    #' Get the indices of the points contained within Ireland
    contained_indices <- unlist(st_contains(CroppedLithuania_UTM, sputm))

    #' Create a logical vector where TRUE indicates the point is in Ireland
    in_Lithuania <- rep(FALSE, nrow(Data.i))
    in_Lithuania[contained_indices] <- TRUE

    #' Update the 'w.pm' column of Data.i.
    Data.i$w.pm[!in_Lithuania] <- NA

    #' Add mesh and prior info for plotting.
    Data.i$PriorRange <- RangePriorValues[j]
    Data.i$MeshName   <- MeshSizes[i]

    #' Combine the w.pm with those of other meshes and
    #' ranges.
    MyData <- rbind(MyData,
                    Data.i)

    k <- k + 1
  }}


#' Add some nice labels for ggplot.
MyData$PriorRange <- factor(MyData$PriorRange,
                            levels = RangePriorValues,
                            labels = paste("P(Range<",
                                           RangePriorValues, ")",
                                           sep = ""))

MyData$MeshName   <- factor(MyData$MeshName,
                            levels = MeshSizes,
                            labels = paste("w=",
                                           MeshSizes,
                                           sep = ""))


#' Plot the w.pm values for each mesh and prior.
ggplot(data = CroppedLithuania_UTM) +
  coord_fixed() +
  geom_sf(fill = "transparent") +
  geom_raster(data = MyData,
              aes(x = Xkm,
                  y = Ykm,
                  fill = w.pm )) +
  scale_fill_gradient2(midpoint = 0,
                       low = "#21918c",
                       mid = "white",
                       high = "#440154",
                       na.value = NA) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())    +
  guides(fill = guide_legend(title = NULL)) +
  labs(title = "Spatial random field",
       caption = "Sigma_u = c(2, 0.05)")  +
  facet_grid(PriorRange ~ MeshName)

ggsave("Sensitivity/FigureS4_sensitivity_SpatialRandomField.png", width = 10, height = 8, dpi = 300, bg = "white")

#' Interpretations
#'   -All panels look very similar. This is good news.
#'   -It means that our model is robust to different mesh configurations and priors.
#'   -There is still an underlying pattern that is not being captured by the covariates in the model,
#'    but given that this pattern aligns with precipitation, windspeed, temeprature, and elevation,
#'    we should account for this, albeit in a different way. There could of course be some other
#'    covariates that we have not considered, like agriculture, nutrients, etc. but we are happy for now.


# Potential next steps:
#'   -Plot fitted values versus observed data for this setting.
#'   -Check for underdispersion.
#'   -Check for residual spatial dependency.


#* Subsection 18.11: Conclusions----

#' For this data set and the model that we used:
#'  -Don't use a sensible mesh resolution + sensible priors (likely around 50km).


#* Subsection 18.12: PC Prior density distributions----

# Function to calculate PC prior densities
calculate_pc_prior <- function(x, param_type, theta, alpha) {
  # param_type: "range" or "sigma"
  # theta: rate parameter of the exponential distribution
  # alpha: probability level for the constraint

  if (param_type == "range") {
    # PC prior for range parameter
    density <- theta * exp(-theta * x)
  } else if (param_type == "sigma") {
    # PC prior for standard deviation parameter
    # This follows an exponential distribution on the distance scale d = -log(sigma/sigma0)
    # where sigma0 is a reference value (often set to 1)
    sigma0 <- 1
    d <- -log(x/sigma0)
    density <- theta * exp(-theta * d) * (1/x)
  }

  return(density)
}

# Derive theta parameters from the constraints:
# P(range < 75) = 0.5
# P(sigma < 2) = 0.05

# For range: P(range < 75) = 0.5
# In exponential distribution: P(X < x) = 1 - exp(-theta * x)
# Therefore: 0.5 = 1 - exp(-theta_range * 75)
# Solving for theta_range:
theta_range <- -log(1 - 0.5) / 75

# For sigma: P(sigma < 2) = 0.05
# This is more complex because PC prior for sigma is not directly exponential
# We need to solve: 0.05 = P(sigma < 2)
# For a PC prior on precision (1/sigma^2), we have P(sigma > sigma0) = exp(-theta * (-log(sigma/sigma0)))
# So P(sigma < 2) = 1 - exp(-theta * (-log(2)))
theta_sigma <- -log(1 - 0.05) / (-log(2))

# Create plots
par(mfrow=c(2,1), mar=c(4,4,2,1))

# Plot for range PC prior
range_values <- seq(0.1, 200, length.out=1000)
range_densities <- calculate_pc_prior(range_values, "range", theta_range, 0.5)

plot(range_values, range_densities, type="l",
     xlab="Range", ylab="Density",
     main="PC Prior for Range Parameter with P(range < 75) = 0.5",
     col="blue", lwd=2)
abline(v=75, col="red", lty=2)
abline(v=29.5785, col="green3", lty=2, lwd=2)
text(80, max(range_densities)*0.8, "P(range < 75) = 0.5", pos=4, col="red")
text(30, max(range_densities)*0.9, "Estimated Range: 29.5785", pos=4, col="green3")

# Plot for sigma PC prior
sigma_values <- seq(0.01, 5, length.out=1000)
sigma_densities <- calculate_pc_prior(sigma_values, "sigma", theta_sigma, 0.05)

plot(sigma_values, sigma_densities, type="l",
     xlab="Sigma", ylab="Density",
     main="PC Prior for Sigma Parameter with P(sigma < 2) = 0.05",
     col="blue", lwd=2)
abline(v=2, col="red", lty=2)
abline(v=1.079655, col="green3", lty=2, lwd=2)
text(2.1, -2, "P(sigma < 2) = 0.05", pos=4, col="red")
text(1.1, -1, "Estimated Sigma: 1.079655", pos=4, col="green3")

# Range PC prior with ggplot2
range_df <- data.frame(range = range_values, density = range_densities)
p1 <- ggplot(range_df, aes(x = range, y = density)) +
  geom_line(color = "blue", size = 1) +
  geom_vline(xintercept = 75, linetype = "dashed", color = "red") +
  geom_vline(xintercept = I1$summary.hyperpar$mean[[2]], linetype = "dashed", color = "green3", size = 1) +
  annotate("text", x = 75+1, y = max(range_densities)*0.8, label = "Prior Range: P(range < 75) = 0.5", hjust = 0, color = "red") +
  annotate("text", x = I1$summary.hyperpar$mean[[2]]+1, y = max(range_densities)*0.9, label = paste("Estimated Range: ", round(I1$summary.hyperpar$mean[[2]], 3), " km", sep = ""), hjust = 0, color = "green3") +
  labs(title = "PC Prior for Range Parameter",
       x = "Range", y = "Density") +
  theme_minimal() +
  My_theme
p1

# Sigma PC prior with ggplot2
sigma_df <- data.frame(sigma = sigma_values, density = sigma_densities)
p2 <- ggplot(sigma_df, aes(x = sigma, y = density)) +
  geom_line(color = "blue", size = 1) +
  geom_vline(xintercept = 2, linetype = "dashed", color = "red") +
  geom_vline(xintercept = I1$summary.hyperpar$mean[[3]], linetype = "dashed", color = "green3", size = 1) +
  annotate("text", x = 2.05, y = 0.5, label = "Prior sigma: P(sigma < 2) = 0.05", hjust = 0, color = "red") +
  annotate("text", x = I1$summary.hyperpar$mean[[3]]+0.05, y = -0.5, label = paste("Estimated sigma: ", round(I1$summary.hyperpar$mean[[3]], 3), sep = ""), hjust = 0, color = "green3") +
  labs(title = "PC Prior for Sigma Parameter",
       x = "Sigma", y = "Density") +
  theme_minimal() +
  My_theme
p2

# Combine the plots in one row with three columns
combined_plot <- p1 + p2 +
  patchwork::plot_layout(nrow = 2) +
  patchwork::plot_annotation(tag_levels = 'A')  # This adds a, b, c labels

# Display the combined plot
combined_plot

ggsave("Sensitivity/FigureS8_DensityDistribution_PCPriors.png", width = 10, height = 10, dpi = 300, bg = "white")

#' Section 19: Clean up----
library(pacman)

###############################################################################################################
# CLEAN UP WORKSPACE
rm(list = ls())       # Remove all objects from environment
gc()                  # Frees up unused memory
p_unload(all)         # Unload all loaded packages
graphics.off()        # Close all graphical devices
cat("\014")           # Clear the console
# Clear mind :)

