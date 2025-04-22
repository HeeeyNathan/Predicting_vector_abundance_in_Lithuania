# Visualising the fitted model----

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
library(patchwork)
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
par(mfrow = c(1,1), mar = c(5,5,2,2))
hist(D,
     freq = TRUE,
     main = "",
     xlab = "Distance between sites (km)",
     ylab = "Frequency")
abline(v = 75, col = "red", lwd = 2, lty = 2)

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

#* Subsection 7.3: Define the projector matrix A----

#' We now define the projector matrix. This is used to
#' calculate:  u = A * w.
A1 <- inla.spde.make.A(mesh1, loc = Loc)
dim(A1)  #' 2089 sites and 2960 nodes in the mesh.


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

#* Subsection 7.5: Define the spatial field----

#' Next, we define the spatial random intercept u.
#' For mesh1 we use:
w1.index <- inla.spde.make.index(name = 'w',
                                 n.spde = spde1$n.spde,
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


#* Subsection 7.7: Specify the model formula----

#' Specify the model formula in terms of the response variable,
#' covariates and the spatial correlated term. Having the colnames
#' is handy at this stage:
colnames(X)

#' This is a model with spatial dependency, based on mesh 1.
f1 <-  Counts ~ -1 + Intercept + eqr.std + ppt.std + tmin.std + ws.std + elevation.std +
                     agriculture.std + artificial.std + natural.std +
                     year.std + waterbody_type.L +
                     f(w, model = spde1)

#* Subsection 7.8: Execute the INLA models----

#' This is the spatial NB GLM based on mesh1.
I1 <- inla(f1,
            family = "nbinomial",
            data = inla.stack.data(Stack1),
            control.compute = MyControlCompute,
            control.predictor = list(A = inla.stack.A(Stack1)))


#' To visualise the model fit, we will plot the LogAlt effect for
#' each of the two Forested values, and for average SDI. We therefore
#' create a MyData object with 50 LogAlt value for each Forested level.


# Create combined data frames for each covariate
# We only need one set of predictions since we're not splitting by waterbody type
MyDataCombined1 <- data.frame(eqr.std = seq(min(df$eqr.std), max(df$eqr.std), length = 100))
MyDataCombined2 <- data.frame(ppt.std = seq(min(df$ppt.std), max(df$ppt.std), length = 100))
MyDataCombined3 <- data.frame(tmin.std = seq(min(df$tmin.std), max(df$tmin.std), length = 100))
MyDataCombined4 <- data.frame(elevation.std = seq(min(df$elevation.std), max(df$elevation.std), length = 100))

# Set all other covariates to their mean values
# EQR
MyDataCombined1$ppt.std         <- mean(df$ppt.std)
MyDataCombined1$tmin.std        <- mean(df$tmin.std)
MyDataCombined1$ws.std          <- mean(df$ws.std)
MyDataCombined1$elevation.std   <- mean(df$elevation.std)
MyDataCombined1$agriculture.std <- mean(df$agriculture.std)
MyDataCombined1$artificial.std  <- mean(df$artificial.std)
MyDataCombined1$natural.std     <- mean(df$natural.std)
MyDataCombined1$year.std        <- (max(df$year) - mean(df$year)) / sd(df$year)
MyDataCombined1$waterbody_type  <- factor(levels(df$waterbody_type)[1], levels = levels(df$waterbody_type)) # For waterbody_type, create a proper factor with the same levels as in the original data

# PPT
MyDataCombined2$eqr.std         <- mean(df$eqr.std)
MyDataCombined2$tmin.std        <- mean(df$tmin.std)
MyDataCombined2$ws.std          <- mean(df$ws.std)
MyDataCombined2$elevation.std   <- mean(df$elevation.std)
MyDataCombined2$agriculture.std <- mean(df$agriculture.std)
MyDataCombined2$artificial.std  <- mean(df$artificial.std)
MyDataCombined2$natural.std     <- mean(df$natural.std)
MyDataCombined2$year.std        <- (max(df$year) - mean(df$year)) / sd(df$year)
MyDataCombined2$waterbody_type  <- factor(levels(df$waterbody_type)[1], levels = levels(df$waterbody_type))

# TMIN
MyDataCombined3$eqr.std         <- mean(df$eqr.std)
MyDataCombined3$ppt.std         <- mean(df$ppt.std)
MyDataCombined3$ws.std          <- mean(df$ws.std)
MyDataCombined3$elevation.std   <- mean(df$elevation.std)
MyDataCombined3$agriculture.std <- mean(df$agriculture.std)
MyDataCombined3$artificial.std  <- mean(df$artificial.std)
MyDataCombined3$natural.std     <- mean(df$natural.std)
MyDataCombined3$year.std        <- (max(df$year) - mean(df$year)) / sd(df$year)
MyDataCombined3$waterbody_type  <- factor(levels(df$waterbody_type)[1], levels = levels(df$waterbody_type))

# ELEVATION
MyDataCombined4$eqr.std         <- mean(df$eqr.std)
MyDataCombined4$ppt.std         <- mean(df$ppt.std)
MyDataCombined4$tmin.std        <- mean(df$tmin.std)
MyDataCombined4$ws.std          <- mean(df$ws.std)
MyDataCombined4$agriculture.std <- mean(df$agriculture.std)
MyDataCombined4$artificial.std  <- mean(df$artificial.std)
MyDataCombined4$natural.std     <- mean(df$natural.std)
MyDataCombined4$year.std        <- (max(df$year) - mean(df$year)) / sd(df$year)
MyDataCombined4$waterbody_type  <- factor(levels(df$waterbody_type)[1], levels = levels(df$waterbody_type))

# Make the X matrix using model.matrix()
XpCombined1 <- model.matrix(~ eqr.std + ppt.std + tmin.std + ws.std + elevation.std +
                     agriculture.std + artificial.std + natural.std +
                     year.std + waterbody_type,
                     data = MyDataCombined1)
XpCombined1 <- as.data.frame(XpCombined1)

XpCombined2 <- model.matrix(~ eqr.std + ppt.std + tmin.std + ws.std + elevation.std +
                     agriculture.std + artificial.std + natural.std +
                     year.std + waterbody_type,
                     data = MyDataCombined2)
XpCombined2 <- as.data.frame(XpCombined2)

XpCombined3 <- model.matrix(~ eqr.std + ppt.std + tmin.std + ws.std + elevation.std +
                     agriculture.std + artificial.std + natural.std +
                     year.std + waterbody_type,
                     data = MyDataCombined3)
XpCombined3 <- as.data.frame(XpCombined3)

XpCombined4 <- model.matrix(~ eqr.std + ppt.std + tmin.std + ws.std + elevation.std +
                     agriculture.std + artificial.std + natural.std +
                     year.std + waterbody_type,
                     data = MyDataCombined4)
XpCombined4 <- as.data.frame(XpCombined4)

# Create prediction stacks
StackPredCombined1 <- inla.stack(
               tag = "Predict",
               data = list(Counts = NA),
               A = list(1, 1),
               effects = list(
                  Intercept = rep(1, nrow(XpCombined1)),
                  XpCombined1 = XpCombined1[,-1]))

StackPredCombined2 <- inla.stack(
               tag = "Predict",
               data = list(Counts = NA),
               A = list(1, 1),
               effects = list(
                  Intercept = rep(1, nrow(XpCombined2)),
                  XpCombined2 = XpCombined2[,-1]))

StackPredCombined3 <- inla.stack(
               tag = "Predict",
               data = list(Counts = NA),
               A = list(1, 1),
               effects = list(
                  Intercept = rep(1, nrow(XpCombined3)),
                  XpCombined3 = XpCombined3[,-1]))

StackPredCombined4 <- inla.stack(
               tag = "Predict",
               data = list(Counts = NA),
               A = list(1, 1),
               effects = list(
                  Intercept = rep(1, nrow(XpCombined4)),
                  XpCombined4 = XpCombined4[,-1]))

# Combine with the original stack
All.stacksCombined1 <- inla.stack(Stack1, StackPredCombined1)
All.stacksCombined2 <- inla.stack(Stack1, StackPredCombined2)
All.stacksCombined3 <- inla.stack(Stack1, StackPredCombined3)
All.stacksCombined4 <- inla.stack(Stack1, StackPredCombined4)

# Run INLA with the combined stacks
I1.PredCombined1 <- inla(f1,
                 family = "nbinomial",
                 data = inla.stack.data(All.stacksCombined1),
                 control.compute = list(dic = TRUE,
                                        waic = TRUE),
                 control.predictor = list(link = 1,
                                 A = inla.stack.A(All.stacksCombined1)))

I1.PredCombined2 <- inla(f1,
                 family = "nbinomial",
                 data = inla.stack.data(All.stacksCombined2),
                 control.compute = list(dic = TRUE,
                                        waic = TRUE),
                 control.predictor = list(link = 1,
                                 A = inla.stack.A(All.stacksCombined2)))

I1.PredCombined3 <- inla(f1,
                 family = "nbinomial",
                 data = inla.stack.data(All.stacksCombined3),
                 control.compute = list(dic = TRUE,
                                        waic = TRUE),
                 control.predictor = list(link = 1,
                                 A = inla.stack.A(All.stacksCombined3)))

I1.PredCombined4 <- inla(f1,
                 family = "nbinomial",
                 data = inla.stack.data(All.stacksCombined4),
                 control.compute = list(dic = TRUE,
                                        waic = TRUE),
                 control.predictor = list(link = 1,
                                 A = inla.stack.A(All.stacksCombined4)))

# Extract the relevant indices
index.FitCombined1  <- inla.stack.index(All.stacksCombined1, tag = "Fit")$data
index.PredCombined1 <- inla.stack.index(All.stacksCombined1, tag = "Predict")$data

index.FitCombined2  <- inla.stack.index(All.stacksCombined2, tag = "Fit")$data
index.PredCombined2 <- inla.stack.index(All.stacksCombined2, tag = "Predict")$data

index.FitCombined3  <- inla.stack.index(All.stacksCombined3, tag = "Fit")$data
index.PredCombined3 <- inla.stack.index(All.stacksCombined3, tag = "Predict")$data

index.FitCombined4  <- inla.stack.index(All.stacksCombined4, tag = "Fit")$data
index.PredCombined4 <- inla.stack.index(All.stacksCombined4, tag = "Predict")$data

# Extract the predicted values
Fit2.fitCombined1  <- I1.PredCombined1$summary.fitted.values[index.FitCombined1, c(1,3,5)]
Fit2.predCombined1 <- I1.PredCombined1$summary.fitted.values[index.PredCombined1, c(1,3,5)]

Fit2.fitCombined2  <- I1.PredCombined2$summary.fitted.values[index.FitCombined2, c(1,3,5)]
Fit2.predCombined2 <- I1.PredCombined2$summary.fitted.values[index.PredCombined2, c(1,3,5)]

Fit2.fitCombined3  <- I1.PredCombined3$summary.fitted.values[index.FitCombined3, c(1,3,5)]
Fit2.predCombined3 <- I1.PredCombined3$summary.fitted.values[index.PredCombined3, c(1,3,5)]

Fit2.fitCombined4  <- I1.PredCombined4$summary.fitted.values[index.FitCombined4, c(1,3,5)]
Fit2.predCombined4 <- I1.PredCombined4$summary.fitted.values[index.PredCombined4, c(1,3,5)]

# Add the predicted values to the MyData objects
MyDataCombined1$mu    <- Fit2.predCombined1[, "mean"]
MyDataCombined1$selow <- Fit2.predCombined1[, "0.025quant"]
MyDataCombined1$seup  <- Fit2.predCombined1[, "0.975quant"]

MyDataCombined2$mu    <- Fit2.predCombined2[, "mean"]
MyDataCombined2$selow <- Fit2.predCombined2[, "0.025quant"]
MyDataCombined2$seup  <- Fit2.predCombined2[, "0.975quant"]

MyDataCombined3$mu    <- Fit2.predCombined3[, "mean"]
MyDataCombined3$selow <- Fit2.predCombined3[, "0.025quant"]
MyDataCombined3$seup  <- Fit2.predCombined3[, "0.975quant"]

MyDataCombined4$mu    <- Fit2.predCombined4[, "mean"]
MyDataCombined4$selow <- Fit2.predCombined4[, "0.025quant"]
MyDataCombined4$seup  <- Fit2.predCombined4[, "0.975quant"]

# Convert standardized EQR back to original scale for plotting
# We need the mean and sd of original covariates
eqr_mean <- mean(df$eqr)
eqr_sd <- sd(df$eqr)

ppt_mean <- mean(df$ppt)
ppt_sd <- sd(df$ppt)

tmin_mean <- mean(df$tmin)
tmin_sd <- sd(df$tmin)

elevation_mean <- mean(df$elevation)
elevation_sd <- sd(df$elevation)

# Convert standardized values back to original scales
MyDataCombined1$eqr       <- MyDataCombined1$eqr.std * eqr_sd + eqr_mean
MyDataCombined2$ppt       <- MyDataCombined2$ppt.std * ppt_sd + ppt_mean
MyDataCombined3$tmin      <- MyDataCombined3$tmin.std * tmin_sd + tmin_mean
MyDataCombined4$elevation <- MyDataCombined4$elevation.std * elevation_sd + elevation_mean


# Define the plotting theme
My_theme <- theme(panel.background = element_blank(),
                  panel.border = element_rect(fill = NA, linewidth = 1.25),
                  strip.background = element_rect(fill = "white",
                                                 color = "white", linewidth = 1.25),
                  text = element_text(size = 16))

#' EQR PLOT WITH COMBINED DATA
# Define ecological quality classes and their colors with explicit ordering
eco_classes <- data.frame(
  class = factor(c("Bad", "Poor", "Moderate", "Good", "High"),
                levels = c("Bad", "Poor", "Moderate", "Good", "High")),
  start = c(0.0, 0.3, 0.4, 0.6, 0.8),
  end = c(0.3, 0.4, 0.6, 0.8, 1.0),
  color = c("#FF0000", "#FF6600", "#FFD700", "#00AA00", "#0000FF")
)

# Create the EQR plot with combined data
P1_combined <- ggplot() +
  # Add ecological quality class bands
  geom_rect(data = eco_classes,
            aes(xmin = start,
                xmax = end,
                ymin = 0,
                ymax = Inf,
                fill = class),
            alpha = 0.4) +

  # Set colors for ecological quality classes
  scale_fill_manual(values = setNames(eco_classes$color, eco_classes$class),
                  name = "EQC",
                  breaks = levels(eco_classes$class),
                  guide = guide_legend(override.aes = list(alpha = 0.5))) +

  # Add observed data points colored by waterbody type
  geom_point(data = df,
             aes(x = eqr,
                 y = Counts,
                 color = waterbody_type),
             alpha = 0.6,
             size = 2) +

  # Color scale for points
  scale_color_manual(values = c("LAKE" = "darkblue", "RIVER" = "darkgreen"),
                     name = "Waterbody Type") +

  # Add the fitted line from the combined model
  geom_line(data = MyDataCombined1,
            aes(x = eqr,
                y = mu),
            size = 1.2,
            color = "black") +

  # Add the credible intervals as ribbons
  geom_ribbon(data = MyDataCombined1,
              aes(x = eqr,
                  ymax = seup,
                  ymin = selow),
              alpha = 0.2,
              fill = "black") +

  # Limit y-axis range
  ylim(range(df$Counts[df$Counts <= 250])) +

  # Set x-axis breaks at specified intervals
  scale_x_continuous(breaks = c(0.0, 0.3, 0.4, 0.6, 0.8, 1.0)) +

  # Labels and theme
  labs(x = "Ecological Quality Ratio",
       y = "Vector Abundance") +

  # Custom theming
  theme_minimal() +
  My_theme +
  theme(
    legend.title.position = "top",
    legend.title.align = 0.5,
    legend.position = "bottom",
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black")
  )

#' PPT PLOT WITH COMBINED DATA
# Create the PPT plot with combined data
# Define the alpha value for gradient
rect_alpha <- 0.7
# Create a function to blend colors with the white background
blend_with_alpha <- function(color, alpha, bg = "white") {
  rgb_col <- col2rgb(color) / 255
  rgb_bg <- col2rgb(bg) / 255
  blended <- rgb_bg * (1 - alpha) + rgb_col * alpha
  rgb(blended[1,], blended[2,], blended[3,])
}
# Blended colors for legend
light_blue_alpha <- blend_with_alpha("lightblue", rect_alpha)
dark_blue_alpha <- blend_with_alpha("#0000FF", rect_alpha)
# Determine the precipitation range
ppt_range <- range(df$ppt, na.rm = TRUE)
# Create a hidden dataset for the legend
legend_data <- data.frame(
  x = seq(ppt_range[1], ppt_range[2], length.out = 100),
  y = rep(0, 100)
)
# Create segments for the gradient
n_segments <- 100
ppt_segments <- seq(ppt_range[1], ppt_range[2], length.out = n_segments + 1)
# Create gradient dataframe
gradient_df <- data.frame(
  xmin = ppt_segments[-length(ppt_segments)],
  xmax = ppt_segments[-1],
  ymin = 0,
  ymax = Inf,
  precipitation = (ppt_segments[-length(ppt_segments)] + ppt_segments[-1]) / 2
)
# Create the PPT plot
P2_combined <- ggplot() +
  # Add the gradient band
  geom_rect(data = gradient_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = precipitation),
            alpha = rect_alpha) +
  # Configure the gradient fill
  scale_fill_gradient(
    low = "lightblue",
    high = "#0000FF",
    guide = "none"
  ) +
  # Add observed data points colored by waterbody type
  geom_point(data = df,
             aes(x = ppt,
                 y = Counts,
                 color = waterbody_type),
             alpha = 0.6,
             size = 2) +
  # Add the fitted line
  geom_line(data = MyDataCombined2,
            aes(x = ppt,
                y = mu),
            size = 1.2,
            color = "black") +
  # Add the credible intervals
  geom_ribbon(data = MyDataCombined2,
              aes(x = ppt,
                  ymax = seup,
                  ymin = selow),
              alpha = 0.2,
              fill = "black") +
  # Color scale for points
  scale_color_manual(values = c("LAKE" = "darkblue", "RIVER" = "darkgreen"),
                     name = "Waterbody Type") +
  # Add hidden points for precipitation color bar (using a new aesthetic)
  geom_point(
    data = legend_data,
    aes(x = x, y = y, fill = x),
    alpha = 0
  ) +
  # Configure the color scale for the hidden points
  scale_fill_gradient(
    low = light_blue_alpha,
    high = dark_blue_alpha,
    name = "mm",
    guide = guide_colorbar(
      barwidth = 20,
      barheight = 0.5,
      title.position = "top",
      title.hjust = 0.5,
      label.hjust = 0.5
    )
  ) +
  # Limit y-axis range
  ylim(range(df$Counts[df$Counts <= 250])) +
  # Labels and theme
  labs(
       # caption = "All other covariates held at their mean, year set to 2022, and waterbody type not considered",
       x = "Precipitation (mm)",
       y = "Vector Abundance") +
  # Custom theming
  theme_minimal() +
  My_theme +
  theme(
    legend.position = "bottom",
    legend.title.align = 0.5,
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black")
  )


#' TMIN PLOT WITH COMBINED DATA
# Determine temperature range
tmin_range <- range(df$tmin, na.rm = TRUE)
# Blended colors for legend
yellow_alpha <- blend_with_alpha("yellow", rect_alpha)
red_alpha <- blend_with_alpha("red", rect_alpha)
# Create a hidden dataset for the legend
legend_data <- data.frame(
  x = seq(tmin_range[1], tmin_range[2], length.out = 100),
  y = rep(0, 100)
)
# Create segments for the gradient
n_segments <- 100
tmin_segments <- seq(tmin_range[1], tmin_range[2], length.out = n_segments + 1)
# Create gradient dataframe
gradient_df <- data.frame(
  xmin = tmin_segments[-length(tmin_segments)],
  xmax = tmin_segments[-1],
  ymin = 0,
  ymax = Inf,
  temperature = (tmin_segments[-length(tmin_segments)] + tmin_segments[-1]) / 2
)
# Create the TMIN plot
P3_combined <- ggplot() +
  # Add the gradient band
  geom_rect(data = gradient_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = temperature),
            alpha = rect_alpha) +
  # Configure the gradient fill
  scale_fill_gradient(
    low = "yellow",
    high = "red",
    guide = "none"
  ) +
  # Add observed data points colored by waterbody type
  geom_point(data = df,
             aes(x = tmin,
                 y = Counts,
                 color = waterbody_type),
             alpha = 0.6,
             size = 2) +
  # Add the fitted line
  geom_line(data = MyDataCombined3,
            aes(x = tmin,
                y = mu),
            size = 1.2,
            color = "black") +
  # Add the credible intervals
  geom_ribbon(data = MyDataCombined3,
              aes(x = tmin,
                  ymax = seup,
                  ymin = selow),
              alpha = 0.2,
              fill = "black") +
  # Color scale for points
  scale_color_manual(values = c("LAKE" = "darkblue", "RIVER" = "darkgreen"),
                     name = "Waterbody Type") +
  # Add hidden points for temperature color bar (using a different aesthetic)
  geom_point(
    data = legend_data,
    aes(x = x, y = y, fill = x),  # Changed from color to fill
    alpha = 0
  ) +
  # Configure color scale for the temperature gradient hidden points
  scale_fill_gradient(
    low = yellow_alpha,
    high = red_alpha,
    name = "°C",
    guide = guide_colorbar(
      barwidth = 20,
      barheight = 0.5,
      title.position = "top",
      title.hjust = 0.5,
      label.hjust = 0.5
    )
  ) +
  # Limit y-axis range
  ylim(range(df$Counts[df$Counts <= 250])) +
  # Labels and theme
  labs(x = "Minimum Temperature (°C)",
       y = "Vector Abundance") +
  # Custom theming
  theme_minimal() +
  My_theme +
  theme(
    legend.position = "bottom",
    legend.title.align = 0.5,
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black")
  )

#' ELEVATION PLOT WITH COMBINED DATA
# Determine elevation range
elevation_range <- range(df$elevation, na.rm = TRUE)
# Blended colors for legend
brown_alpha <- blend_with_alpha("#8C7853", rect_alpha)
khaki_alpha <- blend_with_alpha("#D2C29D", rect_alpha)
teal_alpha <- blend_with_alpha("#097969", rect_alpha)
# Create a hidden dataset for the legend
legend_data <- data.frame(
  x = seq(elevation_range[1], elevation_range[2], length.out = 100),
  y = rep(0, 100)
)
# Create segments for the gradient
n_segments <- 100
elevation_segments <- seq(elevation_range[1], elevation_range[2], length.out = n_segments + 1)
# Create gradient dataframe
gradient_df <- data.frame(
  xmin = elevation_segments[-length(elevation_segments)],
  xmax = elevation_segments[-1],
  ymin = 0,
  ymax = Inf,
  elevation = (elevation_segments[-length(elevation_segments)] + elevation_segments[-1]) / 2
)
# Create the ELEVATION plot
P4_combined <- ggplot() +
  # Add the gradient band
  geom_rect(data = gradient_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = elevation),
            alpha = rect_alpha) +
  # Configure the gradient fill
  scale_fill_gradientn(
    colors = c("#8C7853", "#D2C29D", "#097969"),
    guide = "none"
  ) +
  # Add observed data points colored by waterbody type
  geom_point(data = df,
             aes(x = elevation,
                 y = Counts,
                 color = waterbody_type),
             alpha = 0.6,
             size = 2) +
  # Add the fitted line
  geom_line(data = MyDataCombined4,
            aes(x = elevation,
                y = mu),
            size = 1.2,
            color = "black") +
  # Add the credible intervals
  geom_ribbon(data = MyDataCombined4,
              aes(x = elevation,
                  ymax = seup,
                  ymin = selow),
              alpha = 0.2,
              fill = "black") +
  # Color scale for waterbody type points
  scale_color_manual(values = c("LAKE" = "darkblue", "RIVER" = "darkgreen"),
                     name = "Waterbody Type") +
  # Add hidden points for elevation color bar (using fill)
  geom_point(
    data = legend_data,
    aes(x = x, y = y, fill = x),  # Using fill instead of color
    alpha = 0
  ) +
  # Configure the fill scale for the hidden points
  scale_fill_gradientn(
    colors = c(brown_alpha, khaki_alpha, teal_alpha),
    name = "m a.s.l.",
    guide = guide_colorbar(
      barwidth = 20,
      barheight = 0.5,
      title.position = "top",
      title.hjust = 0.5,
      label.hjust = 0.5
    )
  ) +
  # Limit y-axis range
  ylim(range(df$Counts[df$Counts <= 250])) +
  # Labels and theme
  labs(x = "Elevation (m a.s.l.)",
       y = "Vector Abundance") +
  # Custom theming
  theme_minimal() +
  My_theme +
  theme(
    legend.position = "bottom",
    legend.title.align = 0.5,
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black")
  )

P6_combined <- (P1_combined + P4_combined + plot_layout(axes = "collect_y")) /
               (P3_combined + P2_combined + plot_layout(axes = "collect_y")) +
               plot_annotation(tag_levels = 'A')  # This adds a, b, c labels

# Save the plot
ggsave("Plots/Figure5_predicted_fixed_effects_without_spatial.png", width = 10, height = 12, bg = "white", dpi = 300)

###############################################################################################################
# CLEAN UP WORKSPACE
rm(list = ls())       # Remove all objects from environment
gc()                  # Frees up unused memory
p_unload(all)         # Unload all loaded packages
graphics.off()        # Close all graphical devices
cat("\014")           # Clear the console
# Clear mind :)
