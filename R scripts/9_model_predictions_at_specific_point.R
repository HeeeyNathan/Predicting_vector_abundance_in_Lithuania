# Predicting using our fitted model----

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
library(multcompView) # compact letter display (CLD) for significant differences
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

summary(I1)

#  Now, let's make some predictions
#'  Can we predict for a specific location? The answer is yes.

#' Let us first write down the fitted model.
#' Here are the results (again).

#' Posterior mean values and 95% CI for the regression parameters:
BetaFinal <- I1$summary.fixed[,c("mean", "0.025quant", "0.975quant")]
print(BetaFinal, digits = 2)

#' Visualise the betas using ggplot2.
#' 1. Create a data frame.
BetaFinal_df <- data.frame(Covariate = rownames(BetaFinal),
                           Mean      = round(BetaFinal[,"mean"], 3),
                           Lower     = round(BetaFinal[,"0.025quant"], 3),
                           Upper     = round(BetaFinal[,"0.975quant"], 3))

#' Explicitly set the order of the covariates to match their original order
BetaFinal_df$Covariate <- factor(BetaFinal_df$Covariate, levels = rownames(BetaFinal))

#' 2. Plot the info in Betas2_df
ggplot(BetaFinal_df, aes(x = Mean, y = Covariate)) +
  geom_pointrange(aes(xmin = Lower, xmax = Upper),
                  fatten = 1.5,
                  color = "blue") +
  geom_vline(aes(xintercept = 0),
             linetype = "dotted",
             col = "red") +
  labs(title = "Posterior means and 95% Credible Intervals",
       x = "Value",
       y = "Covariate") +
  theme_minimal()
# If a 95% CI crosses the dotted red line, then the covariate is important.
# Intercept is the baseline where everything is at the average. Lakes have higher dipteran counts, but this is not statistically important.
# If you use a spatial-temporal model, you may very well find that the confidence intervals will get a bit smaller, making the differences between lakes and rivers statistically important.

#' Here is the posterior mean of the theta:
theta.pd <- I1$marginals.hyperpar$`size for the nbinomial observations (1/overdispersion)`
theta.pm <- inla.emarginal(function(x) x, theta.pd)
round(theta.pm, digits = 2)
BetaFinal_df

#' The model is:
#' Counts_i ~ NB(mu_i, theta = 0.14)
#' E[Counts_i] = mu_i
#' var[Counts_i] = mu_i + mu_i^2 / 0.14

#' mu_i = exp(1.415 +
#'            0.753 * eqr.std +
#'            0.472 * ppt.std +
#'            0.228 * tmin.std +
#'            0.120 * ws.std +
#'            0.636 * elevation.std +
#'            0.206 * agriculture.std +
#'            0.171 * artificial.std +
#'            0.171 * natural.std +
#'           -0.022 * year.std +
#'            0.296 * waterbody_type.L + u_i)

#' get shape files
lithuania.shp <- geoboundaries("Lithuania")
lithuania <- st_transform(lithuania.shp, crs = 4326)
bb <- st_bbox(c(xmin = 20.9, xmax = 26.9, ymin = 53.85, ymax = 56.5), crs = st_crs(lithuania))
CroppedLithuania <- st_crop(lithuania, bb)
CroppedLithuania_UTM <- st_transform(x = CroppedLithuania,
                                     crs = "+proj=utm +zone=34 +north +ellps=WGS84 +units=km +datum=WGS84")

#' Recall that this is the mesh:
p1 <- ggplot() +
         geom_sf(data = CroppedLithuania_UTM, fill = "grey80") +
         geom_point(data = df,
                    aes(x = Xkm, y = Ykm),
                    col = grey(0.2),
                    size = 1,
                    alpha = 0.5) +
         theme_minimal()
p1

#' Suppose that we want to predict for the red point (it is a new point)
#' We are predicting for a point in Lithuania very close to where the birds were sampled (Rybachy on the Curonian Lagoon: 55.285238, 20.970616 Decimal degree; 498133.71, 6126533.56 UTM)

# We have already prepared the covariate data for this point (see related R-scripts)

MyData2 <- read.csv("Outputs/6_prediction_data.csv") |>
  clean_names() |>
  filter(site_id == "predicted_site") |>
  mutate(site_id          = as.factor(site_id),
         sample_id        = as.factor(sample_id),
         date             = as.Date(date, format = "%Y-%m-%d"),
         waterbody_type   = factor(waterbody_type, levels = c("Lake", "River"),
                                                   labels = c("Lake", "River"),
                                                   ordered = T)) |>
  dplyr::rename(natural  = forest_and_semi_natural_areas,
                water    = water_bodies,
                wetlands = wetlands,
                Xkm      = xkm,
                Ykm      = ykm)

#* Plot sampling point ----
p1 + geom_point(data = MyData2,
                aes(x = Xkm,
                    y = Ykm),
                  col = "red",
                  cex = 5)

#' That is this point.
unique(MyData2[, c("latitude", "longitude")])

#' In order to predict the the number of dipteran vectors, we also need covariate values.
#' Because the red dot is not a sampling location, we need to obtain relevant covariate values for it.

MyData2$eqr.std              <- mean(df$eqr.std)                                   # <- we do not know the water quality at the red dot, so we will use the mean from all the years
# MyData2$ppt.std              <- mean(df$ppt.std)                                   # <- we will extract this via the terraclimate data
# MyData2$tmin.std             <- mean(df$tmin.std)                                  # <- we will extract this via the terraclimate data
# MyData2$ws.std               <- mean(df$ws.std)                                    # <- we will extract this via the terraclimate data
# MyData2$elevation.std        <- mean(df$elevation.std)                             # <- we will extract this via google elevation API
MyData2$agriculture           <- 0                                                 # <- we will extract this via Copernicus land cover data, this was NA so we changed it to 0 since the sum was already 1
MyData2$artificial            <- 0                                                 # <- we will extract this via Copernicus land cover data, this was NA so we changed it to 0 since the sum was already 1
# MyData2$natural.std          <- mean(df$natural.std)                               # <- we will extract this via Copernicus land cover data
# MyData2$year.std             <- (2023 - mean(df$year)) / sd(df$year)               # <- we will use the dates that are in MyData2, i.e., the dates of the actual sampling
# MyData2$waterbody_type.L     <- factor("Lake", levels = levels(df$waterbody_type)) # <- The Curonian Spit does not have any rivers, therefore we will use Lakes. We use the already set factor.

#' In other words, we pretend that at the red site we have the following
#' covariate conditions:
MyData2

#* Standardize covariates----
MyData2 <- MyData2 |>
  mutate(
    ppt.std         = (ppt - mean(df$ppt)) / sd(df$ppt),
    q.std           = (q - mean(df$q)) / sd(df$q),
    tmin.std        = (tmin - mean(df$tmin)) / sd(df$tmin),
    tmax.std        = (tmax - mean(df$tmax)) / sd(df$tmax),
    ws.std          = (ws - mean(df$ws)) / sd(df$ws),
    elevation.std   = (elevation - mean(df$elevation)) / sd(df$elevation),
    agriculture.std = (agriculture - mean(df$agriculture)) / sd(df$agriculture),
    artificial.std  = (artificial - mean(df$artificial)) / sd(df$artificial),
    natural.std     = (natural - mean(df$natural)) / sd(df$natural),
    water.std       = (water - mean(df$water)) / sd(df$water),
    wetlands.std    = (wetlands - mean(df$wetlands)) / sd(df$wetlands),
    year.std        = (year - mean(df$year)) / sd(df$year)
  )


#' Because of the stack we need to do the prediction slightly different:
#'  1. Use the stack for the original data.
#'  2. Make a stack for the data for which predictions are needed.
#'  3. Combine them.
#'  4. Run INLA
#'  5. Extract the relevant pieces and plot them.



#' 1. See above for the stack of the observed data.

#' 2. Make a stack for the data for which predictions are needed.
#' For these covariate conditions, we make an X matrix for prediction:
Xp <- model.matrix(~ eqr.std + ppt.std + tmin.std + ws.std + elevation.std +
                     agriculture.std + artificial.std + natural.std +
                     year.std + waterbody_type,
                     data = MyData2)
Xp <- as.data.frame(Xp)
Xp

#' If we predict based on only the covariates, then we can use the
#' stack below. It will calculate:
#'    mu = exp(Intercept + Xp * beta)

StackPred <- inla.stack(
               tag = "Predict",
               data = list(Counts = NA),
               A = list(1, 1),
               effects = list(
               Intercept = rep(1, nrow(Xp)),  #' Intercept
               Xp = Xp[,-1]))                 #' Covariates for prediction..

#' But we also want to include the effect of the spatial term. In other words, we want to predict based on:
#' mu = exp(Intercept + Xp * beta + A * w)

#' This means we also need to make an A matrix for the spatial term.
#' We do that as follows. We first grab the coordinates for which we
#' want to make predictions. In this case, that is only one place
#' (you can easily add more).
LocPred <- MyData2[, c("Xkm", "Ykm")] #' Give coordinates to inla.spde.make.A()
LocPred <- as.matrix(LocPred)         #' Has to be a matrix to avoid error.
A3.Pred <- inla.spde.make.A(mesh1, loc = LocPred)
dim(A3.Pred)  #' 65 rows and 2960 vertices


#' This is our stack for prediction:  exp(Intercept + Xp * beta + A3.Pred * w)
StackPred <- inla.stack(
  tag = "Predict",
  data = list(Counts = NA),
  A = list(1, 1, A3.Pred),
  effects = list(
    Intercept = rep(1, nrow(Xp)),  #' Intercept.
    Xp        = Xp[,-1],           #' Covariates for prediction.
    w         = w1.index))         #' w's from mesh 1.



#' 3. Combine the two stacks.
All.stacks <- inla.stack(Stack1, StackPred)


#' 4. Run INLA with the combined stack.
#' This is the spatial NB GLM based on mesh 3.
Pred2.mesh1 <- inla(f1, # formula from mesh2a model
                    family = "nbinomial",
                    data = inla.stack.data(All.stacks),
                    control.compute = MyControlCompute,
                    control.predictor = list(link = 1,
                                    A = inla.stack.A(All.stacks)))

#' The link = 1 will ensure that the exp() is applied for the
#' predicted values.



#' 5. Extract the relevant pieces and plot it.

#' It is now a little bit a headache to extract the predicted values.
#' This is the crucial part for extracting the correct rows.
index.Fit  <- inla.stack.index(All.stacks, tag = "Fit")$data
index.Pred <- inla.stack.index(All.stacks, tag = "Predict")$data

#' The `tag` option gives an index for the rows belonging to the
#' 'Fit' stack, and an index for the rows of the 'Predict' stack.


#' These allow us to extract the correct rows:
Fit2.fit  <- Pred2.mesh1$summary.fitted.values[index.Fit, c(1,3,5)]   #' 245  by 3
Fit2.pred <- Pred2.mesh1$summary.fitted.values[index.Pred, c(1,3,5)]  #' 100 by 3


#' It is the second one we need as these are for the covariate values in
#' MyData. Extract the relevant information and add this to the
#' MyData object.
MyData2$mu    <- Fit2.pred[, "mean"]
MyData2$selow <- Fit2.pred[, "0.025quant"]
MyData2$seup  <- Fit2.pred[, "0.975quant"]
MyData2

#' Hence, at the red dot, and for the given (made up) covariate values, the expected number of dipteran vectors and its 95% credible intervals.


#' Let's explore the differences between the sampling periods
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


# First, update the period labels in your data
MyData2 <- MyData2|>
  mutate(period_label = case_when(
    period == "old" ~ "2003/2004",
    period == "new" ~ "2018/2019"
  ))

# Make sure the factor levels are in chronological order
MyData2$period_label <- factor(MyData2$period_label,
                              levels = c("2003/2004", "2018/2019"))

# Create a summary of vector abundance by period
abundance_summary <- MyData2|>
  group_by(period_label)|>
  dplyr::summarise(
    mean_abundance = mean(mu),
    median_abundance = median(mu),
    ci_low = mean(selow),
    ci_up = mean(seup),
    n_samples = n()
  )

# Create a year-specific summary for the individual years flipped barplot
year_summary <- MyData2|>
  group_by(year)|>
  dplyr::summarise(
    mean_abundance = mean(mu),
    ci_low = mean(selow),
    ci_up = mean(seup),
    n_samples = n()
  )|>
  # Convert year to factor and ensure chronological order
  mutate(year = factor(year, levels = sort(unique(year))))

# Conduct a statistical test (t-test) to compare means between periods
t_test_result <- t.test(mu ~ period_label, data = MyData2)
t_value <- t_test_result$statistic
p_value <- t_test_result$p.value

# Format the p-value and prepare significance indicator
formatted_p_value <- if(p_value < 0.001) {
  "p < 0.001"
} else {
  paste0("p = ", format(round(p_value, 3), nsmall = 3))
}

is_significant <- p_value < 0.05
sig_label <- if(is_significant) "*" else ""

# Create the boxplot with raw data points (Plot 1)
p1 <- ggplot(MyData2, aes(x = period_label, y = mu, fill = period_label)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "black") +
  {if(is_significant) geom_text(data = data.frame(x = "2018/2019",
                                                  y = max(MyData2$mu[MyData2$period_label == "2018/2019"]) + 0.5),
                               aes(x = x, y = y, label = "*"),
                               size = 12, inherit.aes = FALSE) else NULL} +
  scale_fill_manual(values = c("2003/2004" = "#440154", "2018/2019" = "#21918c")) +
  labs(
    x = "Sampling period",
    y = "Vector abundance (predicted)"
  ) +
  theme_minimal() +
  My_theme +
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
p1

# Create vertical bar plot for combined periods (was horizontal p2)
p2 <- ggplot(abundance_summary, aes(x = period_label, y = mean_abundance, fill = period_label)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  # geom_text(aes(label = paste0(round(mean_abundance, 1))),
  #           position = position_dodge(width = 0.9), vjust = -0.5, size = 5) +
  # Add significance asterisk if significant
  # {if(is_significant) geom_text(data = subset(abundance_summary, period_label == "2018/2019"),
  #                              aes(label = "*", y = ci_up + 0.5),
  #                              size = 10, vjust = 0) else NULL} +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_up),
                width = 0.2, position = position_dodge(0.9)) +
  scale_fill_manual(values = c("2003/2004" = "#440154", "2018/2019" = "#21918c")) +
  labs(
    x = "Sampling period",
    y = "Mean vector abundance (predicted)"
  ) +
  theme_minimal() +
  My_theme +
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

# Perform ANOVA to test for differences between years
year_anova <- aov(mu ~ factor(year), data = MyData2)
summary(year_anova)

# Perform Tukey's HSD test for post-hoc comparisons
tukey_result <- TukeyHSD(year_anova)
print(tukey_result)

# Generate the letters that show which groups differ significantly
tukey_groups <- multcompLetters4(year_anova, tukey_result)

# Extract the letter assignments
letter_assignments <- tukey_groups$`factor(year)`$Letters

# Create a data frame for adding letters to the plot
letter_df <- data.frame(
  year = names(letter_assignments),
  letter = unname(letter_assignments),
  stringsAsFactors = FALSE
)

# Match the letter_df with the year_summary data
letter_df$mean_abundance <- year_summary$mean_abundance[match(letter_df$year, year_summary$year)]
letter_df$y_pos <- year_summary$ci_up[match(letter_df$year, year_summary$year)] + 0.5  # Position letters above the error bars

# Create vertical bar plot for individual years
p3 <- ggplot(year_summary, aes(x = year, y = mean_abundance, fill = as.factor(year))) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_up),
                width = 0.2, position = position_dodge(0.9)) +
  # Add the letters showing significant differences
  geom_text(data = letter_df,
            aes(x = year, y = y_pos + 2, label = letter),
            size = 6, fontface = "bold") +
  scale_fill_manual(values = c("2003" = "#FB8072", "2004" = "#21918c", "2018" = "#414487", "2019" = "#440154")) +
  labs(
    x = "Sampling year",
    y = "Mean vector abundance (predicted)",
    caption = paste0("In A, * indicates Welch Two Sample t-test results: t = ", round(t_value, 2), ", ", formatted_p_value, "\n",
                     "In B, letters indicate significant differences (p < 0.05) based on Tukey's HSD test")
  ) +
  theme_minimal() +
  My_theme +
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
p3

# Combine all four plots using patchwork
# Arrange in a 2x2 grid
combined_plot <- (p1 + p3) +
  plot_layout(axes = "collect_y") +
  plot_annotation(tag_levels = 'A')
combined_plot

ggsave("Plots/Figure7_predicted_vector_abundance.png", plot = combined_plot, width = 20, height = 10, units = "in", dpi = 300, device = "png", bg = "white")
