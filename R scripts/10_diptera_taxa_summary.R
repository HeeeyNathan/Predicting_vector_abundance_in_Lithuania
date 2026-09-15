# Summary of dipteran vector taxa used in the analysis and reproduction of
# the descriptive statistics reported in the manuscript paragraph on dipteran
# diversity. Restricted to the three vector families (Ceratopogonidae,
# Simuliidae, Culicidae) and to the study period 2013-2022, matching the
# downstream model in 5_model_implementation.R.

library(tidyverse)
library(readxl)
library(writexl)
library(janitor)
library(pacman)

###############################################################################################################
# Load raw data and restrict to the site-codes used in the main analysis
site_list <- read_xlsx("Outputs/7_site_list.xlsx") |> pull(1)

all_lf <- read.csv("Data/1_raw_macroinvertebrate_data_long.csv", header = TRUE, sep = ",",
                   stringsAsFactors = FALSE, check.names = FALSE) |>
  as_tibble() |>
  filter(!grepl("delete", notes)) |>
  filter(site_code %in% site_list)

vector_families <- c("Ceratopogonidae", "Simuliidae", "Culicidae")

###############################################################################################################
# Sampling-event frame: every (site, date) the team visited in 2013-2022,
# regardless of whether dipteran vectors were caught. Used as the denominator
# for prevalence statistics.
sampling_events <- all_lf |>
  distinct(sample_id, site_id, site_code, date, year, waterbody_type)

###############################################################################################################
# Vector-only long table
vec_lf <- all_lf |>
  filter(order == "Diptera", family %in% vector_families)

###############################################################################################################
# 1. Taxon list with higher taxonomy and habitat occurrence
taxa_summary <- vec_lf |>
  group_by(phylum, class, order, family, taxon_name) |>
  summarise(
    n_sampling_events = n_distinct(sample_id),
    total_abundance   = sum(abundance, na.rm = TRUE),
    waterbody_type    = paste(sort(unique(tolower(waterbody_type))), collapse = ", "),
    .groups = "drop"
  ) |>
  mutate(
    waterbody_type = case_when(
      waterbody_type == "lake"        ~ "Lake only",
      waterbody_type == "river"       ~ "River only",
      waterbody_type == "lake, river" ~ "Both",
      TRUE                            ~ waterbody_type
    ),
    id_level = case_when(
      str_detect(taxon_name, regex("Gen[._ ]+sp\\.?$", ignore_case = TRUE)) ~ "Family",
      str_detect(taxon_name, regex("[_ ]sp\\.?$",      ignore_case = TRUE)) ~ "Genus",
      TRUE                                                                   ~ "Species"
    )
  ) |>
  arrange(family, taxon_name)

# Percentage of taxa identified at each level (also weighted by abundance)
id_level_summary <- taxa_summary |>
  group_by(id_level) |>
  summarise(
    n_taxa            = n(),
    total_abundance   = sum(total_abundance),
    .groups = "drop"
  ) |>
  mutate(
    pct_taxa      = round(100 * n_taxa / sum(n_taxa), 1),
    pct_abundance = round(100 * total_abundance / sum(total_abundance), 1)
  )

###############################################################################################################
# 2. Paragraph statistics

# Per-sampling-event vector counts (zero-filled for events with no vectors)
event_counts <- sampling_events |>
  left_join(
    vec_lf |>
      group_by(sample_id) |>
      summarise(
        vec_abund = sum(abundance, na.rm = TRUE),
        n_taxa    = n_distinct(taxon_name),
        .groups = "drop"
      ),
    by = "sample_id"
  ) |>
  mutate(
    vec_abund = replace_na(vec_abund, 0),
    n_taxa    = replace_na(n_taxa, 0)
  )

# Per-family per-sampling-event counts
family_event_counts <- sampling_events |>
  dplyr::select(sample_id) |>
  expand_grid(family = vector_families) |>
  left_join(
    vec_lf |>
      group_by(sample_id, family) |>
      summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop"),
    by = c("sample_id", "family")
  ) |>
  mutate(abundance = replace_na(abundance, 0))

# Headline numbers
overall_stats <- tibble(
  total_sampling_events = nrow(event_counts),
  events_with_zero      = sum(event_counts$vec_abund == 0),
  events_with_vectors   = sum(event_counts$vec_abund > 0),
  pct_zero              = mean(event_counts$vec_abund == 0) * 100,
  pct_with_vectors      = mean(event_counts$vec_abund > 0)  * 100,
  total_individuals     = sum(event_counts$vec_abund),
  mean_abund_nonzero    = mean(event_counts$vec_abund[event_counts$vec_abund > 0]),
  events_lt_100_when_present = sum(event_counts$vec_abund > 0 & event_counts$vec_abund < 100),
  min_taxa_when_present = min(event_counts$n_taxa[event_counts$n_taxa > 0]),
  max_taxa_when_present = max(event_counts$n_taxa[event_counts$n_taxa > 0])
)

# Per-family occurrence and abundance (denominator = all sampling events)
family_stats <- family_event_counts |>
  group_by(family) |>
  summarise(
    n_events_present  = sum(abundance > 0),
    total_individuals = sum(abundance),
    .groups = "drop"
  ) |>
  arrange(desc(n_events_present))

# Mean abundance per waterbody type, across all sampling events
waterbody_stats <- event_counts |>
  group_by(waterbody_type) |>
  summarise(
    n_events       = n(),
    mean_abund_all = mean(vec_abund),
    mean_abund_nonzero = mean(vec_abund[vec_abund > 0]),
    .groups = "drop"
  )

###############################################################################################################
# 3. Sampling design statistics by waterbody (sites, communities, season)

# Sites and communities (= sampling events) per waterbody
sampling_design <- sampling_events |>
  group_by(waterbody_type) |>
  summarise(
    n_sites       = n_distinct(site_id),
    n_communities = n_distinct(sample_id),
    .groups = "drop"
  ) |>
  mutate(mean_visits_per_site = n_communities / n_sites)

# Sampling month distribution per waterbody
sampling_months <- sampling_events |>
  mutate(month = lubridate::month(as.Date(date, format = "%d/%m/%Y"))) |>
  count(waterbody_type, month) |>
  arrange(waterbody_type, month)

###############################################################################################################
# Write summary statistics to disk so the manuscript paragraph can be reproduced
write_xlsx(taxa_summary,    "Outputs/8.1_diptera_taxa_summary.xlsx")
write_xlsx(overall_stats,   "Outputs/8.2_diptera_overall_summary.xlsx")
write_xlsx(family_stats,    "Outputs/8.3_diptera_family_summary.xlsx")
write_xlsx(waterbody_stats, "Outputs/8.4_diptera_waterbody_summary.xlsx")
write_xlsx(sampling_design, "Outputs/8.5_sampling_design_summary.xlsx")
write_xlsx(sampling_months, "Outputs/8.6_sampling_months_summary.xlsx")
write_xlsx(id_level_summary,"Outputs/8.7_diptera_id_level_summary.xlsx")

###############################################################################################################
# CLEAN UP WORKSPACE
rm(list = ls())
gc()
p_unload(all)
graphics.off()
cat("\014")
