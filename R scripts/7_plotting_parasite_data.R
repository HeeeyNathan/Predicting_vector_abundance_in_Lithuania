# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(rphylopic)
library(scales)
library(stringr)
library(rphylopic)
library(vegan)
library(patchwork)
library(pacman)

My_theme <- theme(
  # Panel settings
  panel.background = element_blank(),
  panel.border = element_rect(fill = NA, linewidth = 1.25),

  # Strip settings for facets
  strip.background = element_rect(fill = "white", color = "white", linewidth = 1.25),

  # Legend settings
  legend.position = "bottom",        # Default legend position
  # legend.position = "right",       # TOGGLE: Uncomment for right legend
  # legend.position = "top",         # TOGGLE: Uncomment for top legend
  # legend.position = "left",        # TOGGLE: Uncomment for left legend
  # legend.position = "none",        # TOGGLE: Uncomment to hide legend
  # legend.direction = "horizontal", # TOGGLE: Uncomment for horizontal legend (default)
  # legend.direction = "vertical",   # TOGGLE: Uncomment for vertical legend
  # legend.box = "horizontal",       # TOGGLE: Uncomment for side-by-side legends
  # legend.box = "vertical",         # TOGGLE: Uncomment for stacked legends
  # legend.title = element_text(face = "bold"), # TOGGLE: Uncomment for bold legend title
  legend.title = element_blank(),  # TOGGLE: Uncomment to remove legend title
  # legend.text = element_text(size = 10), # TOGGLE: Uncomment for smaller legend text
  # legend.key.size = unit(0.8, "cm"), # TOGGLE: Uncomment for smaller legend keys
  # legend.margin = margin(t = 0, r = 0, b = 0, l = 0), # TOGGLE: Uncomment for tighter legend margins
  # legend.box.spacing = unit(0.2, "cm"), # TOGGLE: Uncomment to bring legend closer to plot

  # Grid settings
  # panel.grid.major = element_blank(), # TOGGLE: Comment out for major grid lines
  # panel.grid.minor = element_blank(), # TOGGLE: Comment out for minor grid lines
  panel.grid.major = element_line(color = "gray90", linewidth = 0.2), # TOGGLE: Uncomment for light major grid lines
  panel.grid.minor = element_line(color = "gray90", linewidth = 0.2), # TOGGLE: Uncomment for light minor grid lines
  # panel.grid.major.x = element_blank(), # TOGGLE: Uncomment to remove only horizontal grid lines
  # panel.grid.major.y = element_blank(), # TOGGLE: Uncomment to remove only vertical grid lines

  # Axis settings - enhanced ticks and consistent text size
  axis.ticks = element_line(color = "black", linewidth = 0.8), # Bolder ticks
  axis.ticks.length = unit(0.25, "cm"),                       # Longer ticks
  axis.text = element_text(size = 14),                        # Text at tick marks (12pt)
  axis.title = element_text(size = 16, face = "bold"),        # Axis labels (14pt, bold)
  # axis.text.x = element_text(angle = 45, hjust = 1),        # TOGGLE: Uncomment for angled x labels
  # axis.title.x = element_blank(),                           # TOGGLE: Uncomment to hide x axis title
  # axis.title.y = element_blank(),                           # TOGGLE: Uncomment to hide y axis title

  # Text settings
  text = element_text(size = 16),
  # text = element_text(size = 12),                           # TOGGLE: Uncomment for smaller text

  # Caption with smaller text
  plot.caption = element_text(size = 10, hjust = 1)           # Smaller caption text, right-aligned
)

pastel_palette <- c(
  "#440154",  # Soft teal
  # "#BEBADA",  # Lavender
  # "#FB8072",  # Salmon pink
  # "#80B1D3",  # Sky blue
  # "#B3DE69",  # Lime green
  # "#FCCDE5",  # Pink
  # "#D9D9D9",  # Light gray
  # "#BC80BD",  # Light purple
  # "#CCEBC5",   # Mint green
  "#21918c"  # Peach
)

# Load your dataset
haemosporidian <- read.csv("Data/3_haemosporidian_parasite_data_long.csv")

# Summarize infection prevalence by period, calculating percentages
infection_summary <- haemosporidian|>
  group_by(period)|>
  summarise(Infected = sum(infection_binary),
            Total = n())|>
  mutate(Not_Infected = Total - Infected,
         Infected_Percentage = (Infected / Total) * 100,
         Not_Infected_Percentage = (Not_Infected / Total) * 100)

# Melt data into long format for easy plotting
infection_long <- infection_summary|>
  dplyr::select(period, Infected_Percentage, Not_Infected_Percentage)|>
  gather(key = "Infection_Status", value = "Percentage", Infected_Percentage, Not_Infected_Percentage)

# Rename infection statuses for clarity
infection_long$Infection_Status <- factor(infection_long$Infection_Status,
                                          levels = c("Not_Infected_Percentage", "Infected_Percentage"),
                                          labels = c("Not Infected", "Infected"))

# Update the period labels for the y-axis
infection_long$period <- factor(infection_long$period,
                                levels = c("old", "new"),
                                labels = c("2003/2004", "2018/2019"))
                                # labels = c("2003/2004\n(n = 112)", "2018/2019\n(n = 95)"))

# Create contingency table for Chi-square test
contingency_table <- table(haemosporidian$period, haemosporidian$infection_binary)

# Perform Chi-square test
chisq_test <- chisq.test(contingency_table)
p_value <- chisq_test$p.value
chi_square_stat <- chisq_test$statistic

# Format the p-value for the caption as requested
formatted_p_value <- if(p_value < 0.001) {
  "p < 0.001"
} else {
  paste0("p = ", format(round(p_value, 3), nsmall = 3))
}

# Plot overall infection prevalence as percentages across the two periods with labels
p1 <- ggplot(infection_long, aes(x = period, y = Percentage, fill = Infection_Status)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(round(Percentage, 1), "%"), color = Infection_Status),
            position = position_stack(vjust = 0.5), size = 5) +
  labs(
       # caption = paste0("Chi-Squared Test: χ² = ", round(chi_square_stat, 2), ", ", formatted_p_value),
       x = "Sampling period",
       y = "Prevalence of infection (%)") +
  scale_fill_manual(
    values = c("Not Infected" = "#440154", "Infected" = "#21918c"),
    name = "Infection Status",
    labels = c("Not Infected", "Infected")
  ) +
  scale_color_manual(values = c("Not Infected" = "white", "Infected" = "black")) +
  guides(
    color = "none",
    fill = guide_legend(
      ncol = 2,
      byrow = TRUE
    )
  ) +
  theme_minimal() +
  My_theme +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.5, "cm"),
    legend.margin = margin(t = 3, r = 0, b = 0, l = 0),
    legend.box.margin = margin(t = 3, r = 0, b = 0, l = 0),
    legend.box.spacing = unit(0.5, "lines"),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
    axis.ticks = element_line(color = "black"),
    plot.caption = element_text(hjust = 0, size = 9, face = "italic"),
    plot.margin = unit(c(1, 1, 1, 1), "lines")
  )
p1

# ggsave("Plots/Figure1_parasite_prevalence_old_vs_new.png", plot = p1, width = 10, height = 10, units = "in", dpi = 300, device = "png", bg = "white")

##### Lets do it by genus now #####
# Create contingency tables for each parasite genus
haem_table <- table(haemosporidian$period, !is.na(haemosporidian$Haemoproteus))
leuco_table <- table(haemosporidian$period, !is.na(haemosporidian$Leucocytozoon))
plasmo_table <- table(haemosporidian$period, !is.na(haemosporidian$Plasmodium))

# Perform Chi-Squared test for Haemoproteus
haem_chi <- chisq.test(haem_table)

# Perform Fisher's Exact Test for Leucocytozoon and Plasmodium
leuco_fisher <- fisher.test(leuco_table)
plasmo_fisher <- fisher.test(plasmo_table)

# Extract p-values from the tests
haem_p_value <- haem_chi$p.value
leuco_p_value <- leuco_fisher$p.value
plasmo_p_value <- plasmo_fisher$p.value

# Simplify p-values for subtitle display
haem_p_simplified <- ifelse(haem_p_value <= 0.05, "p ≤ 0.05", "p > 0.05")
leuco_p_simplified <- ifelse(leuco_p_value <= 0.05, "p ≤ 0.05", "p > 0.05")
plasmo_p_simplified <- ifelse(plasmo_p_value <= 0.05, "p ≤ 0.05", "p > 0.05")

# Calculate the number of birds with mixed infections using the "infection_status" column
mixed_infection_count <- haemosporidian|>
  filter(infection_status == "mixed")|>
  nrow()

# Summarize infection prevalence by parasite genus and period
parasite_summary <- haemosporidian|>
  group_by(period)|>
  summarise(
    Haemoproteus = sum(!is.na(Haemoproteus)) / n() * 100,
    Leucocytozoon = sum(!is.na(Leucocytozoon)) / n() * 100,
    Plasmodium = sum(!is.na(Plasmodium)) / n() * 100
  )|>
  gather(key = "Parasite", value = "Prevalence", Haemoproteus, Leucocytozoon, Plasmodium)

# Update the period labels
parasite_summary$period <- factor(parasite_summary$period,
                                  levels = c("old", "new"),
                                  labels = c("2003/2004", "2018/2019"))


# Create a data frame for significant asterisk labels with test types
significance_labels <- data.frame(
  Parasite = c("Haemoproteus", "Plasmodium"),
  period = rep("2018/2019", 2),
  Prevalence = c(
    max(parasite_summary$Prevalence[parasite_summary$Parasite == "Haemoproteus" & parasite_summary$period == "2018/2019"]) + 2,
    max(parasite_summary$Prevalence[parasite_summary$Parasite == "Plasmodium" & parasite_summary$period == "2018/2019"]) + 2
  ),
  label = c("* Chi-Squared: χ²", "* Fisher's exact")
)

# Plot the prevalence of each parasite genus by period (position = "dodge")
# Plot the prevalence of each parasite genus by period (position = "dodge")
p2 <- ggplot(parasite_summary, aes(x = period, y = Prevalence, fill = Parasite)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = paste0(round(Prevalence, 1), "%")),
            position = position_dodge(width = 0.9), vjust = -0.3, size = 5) +
  # Add the Chi-Squared label for Haemoproteus in #f2cc84
  geom_text(data = significance_labels[significance_labels$Parasite == "Haemoproteus",],
            aes(x = period, y = Prevalence, label = label),
            position = position_dodge(width = 1), size = 6, vjust = -0.1, hjust = 1,
            color = "#440154", fontface = "bold") +
  # Add the Fisher's exact label for Plasmodium in #84a9f2
  geom_text(data = significance_labels[significance_labels$Parasite == "Plasmodium",],
            aes(x = period, y = Prevalence, label = label),
            position = position_dodge(width = 1), size = 6, vjust = -0.1, hjust = 0.0,
            color = "#FB8072", fontface = "bold") +
  labs(
       # title = "Prevalence of haemosporidian genera in Parus caeruleus \nLocation: Curonian Lagoon (Rybachy), juveniles sampled",
       # subtitle = paste0("Mixed infections: ", mixed_infection_count,
       #                   " | Haemoproteus: ", haem_p_simplified,
       #                   " | Leucocytozoon: ", leuco_p_simplified,
       #                   " | Plasmodium: ", plasmo_p_simplified),
       x = "Sampling period",
       y = "Prevalence of infection (%)") +
  scale_fill_manual(values = c("Haemoproteus" = "#440154",
                               "Leucocytozoon" = "#21918c",
                               "Plasmodium" = "#FB8072")) +
  theme_minimal() +
  My_theme +  # Apply the My_theme from the first plot
  theme(
    # Visual styling only - matching the first plot's aesthetics
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.5, "cm"),
    legend.margin = margin(t = 3, r = 0, b = 0, l = 0),
    legend.box.margin = margin(t = 3, r = 0, b = 0, l = 0),
    legend.box.spacing = unit(0.5, "lines"),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
    axis.ticks = element_line(color = "black"),
    plot.margin = unit(c(1, 1, 1, 1), "lines"),
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )
p2

# ggsave("Plots/FigureS8_parasite_prevalence_b_genus_old_vs_new.png", plot = p2, width = 10, height = 10, units = "in", dpi = 300, device = "png", bg = "white")

# combine the plots using pathcwork, add A and B
combined_plot <- p1 + p2 +
  # plot_annotation(
  #   title = "Haemosporidian parasite prevalence in Parus caeruleus",
  #   subtitle = "Location: Curonian Lagoon (Rybachy), juveniles sampled",
  #   caption = "Data source: Rybachy Ornithological Station, 2003-2019"
  # ) +
  plot_layout(nrow = 1) +
  plot_annotation(tag_levels = 'A') # Add A and B labels
combined_plot

ggsave("Plots/Figure6_parasite_prevalence_dynamics.png", plot = combined_plot, width = 20, height = 10, units = "in", dpi = 300, device = "png", bg = "white")
saveRDS(combined_plot, "Plots/Figure6_parasite_prevalence_dynamics.RDS")

###############################################################################################################
# CLEAN UP WORKSPACE
rm(list = ls())       # Remove all objects from environment
gc()                  # Frees up unused memory
p_unload(all)         # Unload all loaded packages
graphics.off()        # Close all graphical devices
cat("\014")           # Clear the console
# Clear mind :)
