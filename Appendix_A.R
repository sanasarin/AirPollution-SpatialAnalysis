#Load libraries
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(ggsn)
library(sp)
library(sf)
library(phylin)
library(dplyr)
library(automap)
library(gstat)
library(readr)
library(raster)

# set working directory with setwd(dir)
setwd("/Users/FROG/ggr376/a3/GGR376-A3")

########## Data Munging

# Read the data
air_data <- read_csv("data/annual_conc_by_monitor_2021.csv", show_col_types = FALSE)

# List all attributes
names(air_data)

#And, let's load a few rows.
air_data

#The relevant attributes for our case are 
#- `Latitude` and `Longitude` -- probably needed later for spatial analysis
#- `State Name` -- filter by states
#- `City Name` and `County Name` may be useful later
#- `Parameter Name` -- filter the three air pollutants (NO2, Ground Ozone, PM2.5)
#- `Year` -- Filter the year 
#- `Arithmetic Mean` -- useful for air pollution concentration
#- `Units of Measure` -- ensure consistency in the pollutant concentration units

#Let's view all values for `Parameter Name` that may contain nitrogen, Ozone or PM2.5
unique(air_data$`Parameter Name`)

#Since that is a lot going on, let's try to filter by substrings. 
unique(air_data$`Parameter Name`[grepl("Nitrogen", air_data$`Parameter Name`)])
#We may filter by `Nitrogen dioxide (NO2)`
unique(air_data$`Parameter Name`[grepl("Ozone", air_data$`Parameter Name`)])
unique(air_data$`Parameter Name`[grepl("O3", air_data$`Parameter Name`)])
unique(air_data$`Parameter Name`[grepl("PM2.5", air_data$`Parameter Name`)])

#Filter by:
#    Ozone - Daily maximum of 8-hour running average	- Pollutant Standard: Ozone 8-hour 2015	
#    PM2.5 - Local Conditions - Daily Mean	- Pollutant Standard: PM25 24-hour 2012
#    Nitrogen dioxide (NO2) - Daily Maximum 1-hour average	- Pollutant Standard: NO2 Annual 1971
# Ensure each monitor only has one measurement by using distinct()

pollutants <- c("Ozone 8-hour 2015", "PM25 24-hour 2012", "NO2 Annual 1971")

air_data_pollutants <- air_data %>%
  filter(`Pollutant Standard` %in% pollutants) %>%
  filter(`Year` == 2021) %>%
  distinct(`State Code`,`County Code`,`Site Num`, `Pollutant Standard`, .keep_all = TRUE)

# Sanity check
unique(air_data_pollutants$`Parameter Code`)
unique(air_data_pollutants$`Pollutant Standard`)

#Let's see all state names 
unique(air_data$`State Name`)

#Our selected states are: New York, Pennsylvania, New Jersey, Delaware, Maryland
#Let's filter by these states so we have a smaller csv file.

selected_states <- c("New York", "Pennsylvania", "New Jersey", "Delaware", "Maryland")

air_data_states_pollutants <- air_data_pollutants %>%
  filter(`State Name` %in% selected_states) %>%
  filter(`Year` == 2021)

# Sanity check
unique(air_data_states_pollutants$`State Name`)

# Let's export this! 
# (This script was originally split into multiple files, 
# thus exporting and importing was necessary)
write_csv(air_data_states_pollutants, "data/air_data_states_pollutants.csv")

#Map the points by pollutant
no2_data <- air_data_states_pollutants %>%
  filter(`Parameter Name` == "Nitrogen dioxide (NO2)")
ozone_data <- air_data_states_pollutants %>%
  filter(`Parameter Name` == "Ozone")
pm_data <- air_data_states_pollutants %>%
  filter(`Parameter Name` == "PM2.5 - Local Conditions")

# convert to sf object and transform to EPSG 102004 (NAD 1983 Lambert contiguous USA)
crs = st_crs("EPSG:4326")
new_crs = st_crs("ESRI:102004")

states <- sf::st_read("data/states.shp") %>% st_transform(new_crs)
state_names <- c("Delaware", "Maryland", "New Jersey", "New York", "Pennsylvania")
states_filtered <- states %>% filter(NAME %in% selected_states)

no2_df_sf <- st_as_sf(no2_data, coords = c("Longitude", "Latitude"), crs = crs)%>%
  st_transform(new_crs)
ozone_df_sf <- st_as_sf(ozone_data, coords = c("Longitude", "Latitude"), crs = crs)%>%
  st_transform(new_crs)
pm_df_sf <- st_as_sf(pm_data, coords = c("Longitude", "Latitude"), crs = crs)%>%
  st_transform(new_crs)

state_centroids <- st_centroid(states, of_largest_polygon = TRUE) %>% 
  cbind(st_coordinates(.)) %>%
  as.data.frame() %>%
  dplyr::select(NAME, X, Y)
# Error appears but the centroids were still created

no2_df <- no2_df_sf %>% cbind(st_coordinates(.)) %>% as.data.frame()
ozone_df <- ozone_df_sf %>% cbind(st_coordinates(.)) %>% as.data.frame()
pm_df <- pm_df_sf %>% cbind(st_coordinates(.)) %>% as.data.frame()

# position of text labels
state_centroids <- state_centroids %>%
  mutate(nudge_x = case_when(
           NAME == "Delaware" ~ 0.05,
           NAME == "Maryland" ~ 0.05,
           NAME == "New Jersey" ~ 0.05,
           NAME == "New York" ~ 0.05,
           NAME == "Pennsylvania" ~ 0.05
         ),
         nudge_y = case_when(
           NAME == "Delaware" ~ 0.05,
           NAME == "Maryland" ~ -0.05,
           NAME == "New Jersey" ~ 0.05,
           NAME == "New York" ~ 0.05,
           NAME == "Pennsylvania" ~ -0.05
         ))

ggplot() +
  geom_sf(data = states_filtered, fill = "#eeeeee", color = "black", size = 0.1) +
  geom_text_repel(data = state_centroids, aes(x = X, y = Y, label = NAME, nudge_x =nudge_x, nudge_y=nudge_y), size = 2.5, fontface = "bold", force=3)+
  geom_point(data=no2_df, aes(x = X, y= Y, colour="Nitrogen dioxide (NO2)"), size=1.6) +
  geom_point(data=ozone_df, aes(x = X, y= Y, colour="Ozone"), size=1) +
  geom_point(data=pm_df, aes(x = X, y= Y, colour="PM2.5 - Local Conditions"), size=0.3) +
  labs(title="Monitors by Pollutants in Study Area", y="Latitude (degrees)", x="Longitude (degrees)", colour="Pollutant") +
  scale_color_manual(values = c("red", "turquoise", "blue") )+
  theme_minimal() +
  north(states_filtered)+
  ggsn::scalebar(states_filtered, dist = 100, dist_unit = "km",transform = FALSE, model = 'NAD83', st.size = 2, border.size= 0.5, location="topleft")

# The following is exports for clustering
# Lets export for specific pollutants -- Ozone 
air_data_states_pollutants_ozone <- air_data_states_pollutants %>%
  filter(`Parameter Name` == 'Ozone') %>%
  dplyr::select(Latitude, Longitude, `Arithmetic Mean`, `Parameter Name`)

write_csv(air_data_states_pollutants_ozone, "data/air_data_states_pollutants_ozone_v2.csv")

# NO2
air_data_states_pollutants_no2 <- air_data_states_pollutants %>%
  filter(`Parameter Name` == 'Nitrogen dioxide (NO2)') %>%
  dplyr::select(Latitude, Longitude, `Arithmetic Mean`, `Parameter Name`)

write_csv(air_data_states_pollutants_no2, "data/air_data_states_pollutants_no2.csv")

#PM2.5 
air_data_states_pollutants_pm <- air_data_states_pollutants %>%
  filter(`Parameter Name` == 'PM2.5 - Local Conditions') %>%
  dplyr::select(Latitude, Longitude, `Arithmetic Mean`, `Parameter Name`)

write_csv(air_data_states_pollutants_pm, "data/air_data_states_pollutants_pm.csv")




########## NO2 Interpolation
### NO2 IDW

air_data <- read_csv("data/air_data_states_pollutants.csv", show_col_types = FALSE, guess_max = 9500)
no2_data <- air_data%>%
  filter(`Parameter Name` == "Nitrogen dioxide (NO2)")

# plot no2 data
no2_data %>%
  as.data.frame() %>%
  ggplot(aes(Longitude, Latitude)) + 
  geom_point(color="blue", alpha=3/4) + 
  ggtitle("Nitrogen dioxide (NO2) Concentration") + 
  coord_equal() + 
  theme_bw()

# convert to sf object and transform to EPSG 102004 (NAD 1983 Lambert contiguous USA)
crs = st_crs("EPSG:4326")
new_crs = st_crs("ESRI:102004")
no2_df_sf <- st_as_sf(no2_data, coords = c("Longitude", "Latitude"), crs = crs)%>%
  st_transform(new_crs)

states <- sf::st_read("data/states.shp") %>% st_transform(new_crs)
statesOutline <- fortify(states, region="Name")

#make grid
grd_5000_sf <- states %>% 
  st_make_grid(
    cellsize = c(5000, 5000), # 5000m pixel size
    what = "corners"
  ) %>%
  st_as_sf()
grd_5000_sf  <- grd_5000_sf %>%
  cbind(. , st_coordinates(.))

# get the points within the state boundaries
grd_sp <- as(grd_5000_sf, "Spatial")
states_sp <- as(states, "Spatial")

grid_in_poly <- sp::over(grd_sp, states_sp)
grid_predict <- grd_sp[complete.cases(grid_in_poly),]
plot(grid_predict)

# plot points
ggplot() + 
  geom_sf(data=no2_df_sf) +
  geom_sf(data = statesOutline, color="black", fill=NA) +
  coord_sf(crs=new_crs)

## get coordinates and value to interpolate from sf object
data <- no2_df_sf %>% 
  cbind(. , st_coordinates(.)) %>% 
  as.data.frame() %>% 
  dplyr::select(X,Y, `Arithmetic.Mean`)

## get the coordinates of the grid
grd <- grid_predict %>% 
  as.data.frame() %>% 
  dplyr::select(X, Y)

## perform the IDW and join back the coordinates of the grid
no2_i <- phylin::idw(data[,3], 
                 data[,1:2],
                 grd) %>% 
  cbind(grd)

# Cross validation (Leave one out)

# convert sf object to sp object
no2_spdf <- as(no2_df_sf, "Spatial")

# This will give a plot of the RMSE values for each K value and the value 
# of K with the lowest RMSE will be indicated with a red line (run the code). 
# We can use the value of K with the lowest RMSE for our final model.

# Let's find the best k/p value
# create an empty data frame to store p values and RMSE values
results <- data.frame(p = integer(), RMSE = numeric())

test_val <- c(0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

# loop through different p values and perform cross-validation
for (p in test_val) {
  LOOCV <- krige.cv(formula=Arithmetic.Mean~1, locations=no2_spdf, nfold=nrow(no2_spdf), set = list(idp = p))
  
  # calculate RMSE
  RMSE <- sqrt(mean(LOOCV@data$residual^2))
  # add p and RMSE values to the results data frame
  results <- rbind(results, data.frame(p = p, RMSE = RMSE))
}

# find the p value with the lowest RMSE
best_idw <- results[which.min(results$RMSE),]
cat("The best p value is", best_idw$p, "with an RMSE of", best_idw$RMSE)
# The best p value is 1 with an RMSE of 3.82308

ggplot(results, aes(x = p, y = RMSE)) + 
  geom_line() + 
  geom_point() +
  scale_x_continuous(breaks = seq(1, 20, by = 1)) + 
  labs(x = "p value", y = "RMSE",  title="RMSE for test p values for NO2")

## perform the IDW again
no2_i <- phylin::idw(data[,3], 
                 data[,1:2],
                 grd, p=1) %>% 
  cbind(grd)

# rename Z to Mean
no2_i <- no2_i %>% mutate(Mean = Z) %>% dplyr::select(-Z)

#get summary
summary(no2_i)

# tiled (filled) interpolation, displayed in EPSG 102004 coordinate system
ggplot() + 
  geom_tile(data=no2_i, aes(x=X, y=Y, fill=Mean)) +
  geom_sf(data = statesOutline, color="black", fill=NA) +
  coord_sf(crs=new_crs) +
  labs(title="IDW Interpolated NO2 Values", 
       subtitle="For New York, Pennsylvania, New Jersey, Delaware, Maryland in 2021",
       fill="Mean (ppb)") +
  xlab("Longitude")+
  ylab("Latitude")+
  north(statesOutline)+
  ggsn::scalebar(statesOutline, dist = 100, dist_unit = "km",transform = FALSE, model = 'NAD83', st.size = 2, border.size= 0.5, location="topleft")




### NO2 Kriging

# assumption 1: data is normal
hist(no2_data$`Arithmetic Mean`)
shapiro.test(no2_data$`Arithmetic Mean`)
# p-value = 0.2284, so it is normal

# assumption 2: mean and variance are equal

# Find extent values
# Using functions in the raster package
area <- extent(no2_spdf@bbox)
area

# Create a raster of 5 x 5 covering the area
raster_area <- raster(nrows = 4, ncols = 4, ext = area,
                      vals = -999)
# Convert Raster to Polygons
polygons_area <- rasterToPolygons(raster_area)
crs(polygons_area) <- crs(no2_spdf)

# Removing all columns except arithmetic mean
no2_spdf_mean <- no2_spdf[,26]
# Number of obserations per polygon
assessment <- over(polygons_area, no2_spdf_mean)
# Mean of values
assessment$mean <- over(polygons_area, no2_spdf_mean,
                        fn = mean, na.rm = T)[,1]
# Variance of Values
assessment$var <- over(polygons_area, no2_spdf_mean,
                       fn = var, na.rm = T)[,1]
# Change Names
names(assessment)[1] <- c("counts")
# Replace Polygon Data
polygons_area@data <- assessment
summary(polygons_area@data)
spplot(polygons_area, c("counts", "mean", "var"), names.attr = c("Count","Mean", "Variance"),
       colorkey=list(space="bottom"), scales = list(draw = TRUE),
       main = "Stationarity Test",
       as.table = TRUE)

summary(polygons_area)

# assumption 3: isotrophy
# Look at IDW plot again to see if there's any spatial trend
ggplot() + 
  geom_tile(data=no2_i, aes(x=X, y=Y, fill=Mean)) +
  geom_sf(data = statesOutline, color="white", fill=NA) +
  coord_sf(crs=new_crs)
# seems fine

# filtering variogram

no2_vgm <- variogram(Arithmetic.Mean~1, data = no2_spdf)

# see all types of variograms 
show.vgms()

no2_fit <- fit.variogram(no2_vgm, model = vgm(model="Sph"))
plot(no2_vgm, no2_fit) # plot the sample values, along with the fit model

no2_fit
#   model    psill    range
#   Sph      20.33064 72403.14

# Performing kriging
# for the points grid, we already made that in the previous step for IDW
# Let's look at it:
plot(grid_predict)
plot(states_sp)

#krige
no2_ord_krige <- krige(Arithmetic.Mean~1, no2_spdf, grid_predict, model=no2_fit)

#view 
spplot(no2_ord_krige, "var1.pred")

# plot the variance - estimate of estimation error
spplot(no2_ord_krige, "var1.var")

# cross validation

# create array of models to test
test_models = c("Nug", "Exp", "Sph", "Gau", "Exc", "Mat", "Ste", "Cir",
                "Lin", "Bes", "Pen", "Per", "Wav", "Hol", "Log", "Pow", "Spl")

# create array to store results
results_krige <- data.frame()

# loop through different models and perform cross-validation
for (m in test_models) {
  no2_fit <- fit.variogram(no2_vgm, model = vgm(model=m))
  LOOCV <- krige.cv(Arithmetic.Mean~1, no2_spdf, model=no2_fit)
  
  # calculate RMSE
  RMSE <- sqrt(mean(LOOCV@data$residual^2))
  # add p and RMSE values to the results data frame
  results_krige <- rbind(results_krige, data.frame(model=m, RMSE = RMSE, 
                                                   psill=no2_fit$psill, range=no2_fit$range))
}

# find the model with the lowest RMSE
best_krige <- results_krige[which.min(results_krige$RMSE),]
cat("The best model is", best_krige$m, "with an RMSE of", best_krige$RMSE,
    "psill =", best_krige$psill, "range =", best_krige$range)
results_krige
# model     RMSE    psill    range
# Lin    3.141394  19.8527585    48389.63386

# retry kriging with this model

#krige
no2_fit <- fit.variogram(no2_vgm, model = vgm(model="Lin"))
plot(no2_vgm, no2_fit, main = "NO2 Variogram (Linear)", 
     sub = "RMSE = 3.14  psill = 19.85  range = 48389",
     xlab = "Distance (m)",
     ylab = "Semivariance")

no2_ord_krige <- krige(Arithmetic.Mean~1, no2_spdf, grid_predict, model=no2_fit)

# check sum of squared error
attr(no2_fit, "SSErr")
# 0.000004474448, very low

#view 
spplot(no2_ord_krige, "var1.pred")

# plot the variance - estimate of estimation error
spplot(no2_ord_krige, "var1.var")

#plot the same way as the IDW 
no2_krige_sf <- st_as_sf(no2_ord_krige) %>% 
  cbind(. , st_coordinates(.))
# tiled (filled) interpolation, displayed in EPSG 102004 coordinate system
ggplot() + 
  geom_tile(data=no2_krige_sf, aes(x=X, y=Y, fill=var1.pred)) +
  geom_sf(data = statesOutline, color="black", fill=NA) +
  coord_sf(crs=new_crs) + 
  labs(title="Kriging Interpolated NO2 Values", 
       subtitle="For New York, Pennsylvania, New Jersey, Delaware, Maryland in 2021",
       fill="Mean (ppb)") +
  xlab("Longitude")+
  ylab("Latitude")+
  north(statesOutline)+
  ggsn::scalebar(statesOutline, dist = 100, dist_unit = "km",transform = FALSE, model = 'NAD83', st.size = 2, border.size= 0.5, location="topleft")

# get summary
summary(no2_krige_sf$var1.pred)

########## Ozone Interpolation
### Ozone IDW
ozone_data <- air_data%>%
  filter(`Parameter Name` == "Ozone")

# plot ozone data
ozone_data %>%
  as.data.frame() %>%
  ggplot(aes(Longitude, Latitude)) + 
  geom_point(color="blue", alpha=3/4) + 
  ggtitle("Ozone Concentration") + 
  coord_equal() + 
  theme_bw()


# convert to sf object and transform to EPSG 102004 (NAD 1983 Lambert contiguous USA)

crs = st_crs("EPSG:4326")
new_crs = st_crs("ESRI:102004")
ozone_df_sf <- st_as_sf(ozone_data, coords = c("Longitude", "Latitude"), crs = crs)%>%
  st_transform(new_crs)

# plot points
ggplot() + 
  geom_sf(data=ozone_df_sf) +
  geom_sf(data = statesOutline, color="black", fill=NA) +
  coord_sf(crs=new_crs)

# grid already exists in this file
plot(grid_predict)

## get coordinates and value to interpolate from sf object
data <- ozone_df_sf %>% 
  cbind(. , st_coordinates(.)) %>% 
  as.data.frame() %>% 
  dplyr::select(X,Y, `Arithmetic.Mean`)

## perform the IDW and join back the coordinates of the grid
ozone_i <- phylin::idw(data[,3], 
                      data[,1:2],
                      grd) %>% 
  cbind(grd)


# Cross validation (Leave one out)

# convert sf object to sp object
ozone_spdf <- as(ozone_df_sf, "Spatial")

# This will give a plot of the RMSE values for each K value and the value 
# of K with the lowest RMSE will be indicated with a red line (run the code). 
# We can use the value of K with the lowest RMSE for our final model.

# Let's find the best k/p value
# create an empty data frame to store p values and RMSE values
r <- data.frame(p = integer(), RMSE = numeric())

test_val <- c(0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

# loop through different p values and perform cross-validation
for (p in test_val) {
  LOOCV <- krige.cv(formula=Arithmetic.Mean~1, locations=ozone_spdf, nfold=nrow(ozone_spdf), set = list(idp = p))
  
  # calculate RMSE
  RMSE <- sqrt(mean(LOOCV@data$residual^2))
  # add p and RMSE values to the results data frame
  r <- rbind(r, data.frame(p = p, RMSE = RMSE))
}

# find the p value with the lowest RMSE
best_idw <- r[which.min(r$RMSE),]
cat("The best p value is", best_idw$p, "with an RMSE of", best_idw$RMSE)
# The best p value is 1 with an RMSE of 0.002599975

ggplot(r, aes(x = p, y = RMSE)) + 
  geom_line() + 
  geom_point() +
  scale_x_continuous(breaks = seq(1, 20, by = 1)) + 
  labs(x = "p value", y = "RMSE", title="RMSE for test p values for Ozone")

# do idw with p = 1
ozone_i <- phylin::idw(data[,3], 
                 data[,1:2],
                 grd, p=1) %>% 
  cbind(grd)

# rename Z to Mean
ozone_i <- ozone_i %>% mutate(Mean = Z) %>% dplyr::select(-Z)

#get summary
summary(ozone_i)

# tiled (filled) interpolation, displayed in EPSG 102004 coordinate system
ggplot() + 
  geom_tile(data=ozone_i, aes(x=X, y=Y, fill=Mean)) +
  geom_sf(data = statesOutline, color="black", fill=NA) +
  coord_sf(crs=new_crs) +
  labs(title="IDW Interpolated Ozone Values", 
       subtitle="For New York, Pennsylvania, New Jersey, Delaware, Maryland in 2021",
       fill="Mean (ppm)") +
  xlab("Longitude")+
  ylab("Latitude")+
  north(statesOutline)+
  ggsn::scalebar(statesOutline, dist = 100, dist_unit = "km",transform = FALSE, model = 'NAD83', st.size = 2, border.size= 0.5, location="topleft")




### Ozone Kriging

# assumption 1: data is normal
hist(ozone_data$`Arithmetic Mean`)
shapiro.test(ozone_data$`Arithmetic Mean`)
# It looks normal and the p value > 0.05

# assumption 2: mean and variance are equal

# Find extent values
# Using functions in the raster package

area <- extent(ozone_spdf@bbox)
area

# Create a raster of 5 x 5 covering the area
ozone_raster <- raster(nrows = 4, ncols = 4, ext = area,
                       vals = -999)
# Convert Raster to Polygons
polygons_area <- rasterToPolygons(ozone_raster)
crs(polygons_area) <- crs(ozone_spdf)

# Removing all columns except arithmetic mean
ozone_spdf_mean <- ozone_spdf[,26]
# Number of obserations per polygon
assessment <- over(polygons_area, ozone_spdf_mean)
# Mean of values
assessment$mean <- over(polygons_area, ozone_spdf_mean,
                        fn = mean, na.rm = T)[,1]
# Variance of Values
assessment$var <- over(polygons_area, ozone_spdf_mean,
                       fn = var, na.rm = T)[,1]
# Change Names
names(assessment)[1] <- c("counts")
# Replace Polygon Data
polygons_area@data <- assessment
summary(polygons_area@data)
spplot(polygons_area, c("counts", "mean", "var"), names.attr = c("Count","Mean", "Variance"),
       colorkey=list(space="bottom"), scales = list(draw = TRUE),
       main = "Stationarity Test",
       as.table = TRUE)

summary(polygons_area)

# assumption 3: isotrophy
# Look at IDW plot again to see if there's any spatial trend
ggplot() + 
  geom_tile(data=ozone_i, aes(x=X, y=Y, fill=Mean)) +
  geom_sf(data = statesOutline, color="white", fill=NA) +
  coord_sf(crs=new_crs)

ozone_vgm <- variogram(Arithmetic.Mean~1, data = ozone_spdf)

# see all types of variograms 
show.vgms()

ozone_fit <- fit.variogram(ozone_vgm, model = vgm(model="Sph"))
plot(ozone_vgm, ozone_fit) # plot the sample values, along with the fit model

ozone_fit

# psill: 0.000005988684
# range: 15404.1

# Performing kriging
# for the points grid, we already made that in the previous step for IDW

#krige
ozone_ord_krige <- krige(Arithmetic.Mean~1, ozone_spdf, grid_predict, model=ozone_fit)

#view 
spplot(ozone_ord_krige, "var1.pred")

# plot the variance - estimate of estimation error
spplot(ozone_ord_krige, "var1.var")

# cross validation
# create array of models to test
test_models = c("Nug", "Exp", "Sph", "Gau", "Exc", "Mat", "Ste", "Cir",
                "Lin", "Bes", "Pen", "Per", "Wav", "Hol", "Log", "Pow", "Spl")

# create array to store results
results_krige <- data.frame()

# loop through different models and perform cross-validation
for (m in test_models) {
  ozone_fit <- fit.variogram(ozone_vgm, model = vgm(model=m))
  LOOCV <- krige.cv(Arithmetic.Mean~1, ozone_spdf, model=ozone_fit)
  
  # calculate RMSE
  RMSE <- sqrt(mean(LOOCV@data$residual^2))
  # add p and RMSE values to the results data frame
  results_krige <- rbind(results_krige, data.frame(model=m, RMSE = RMSE, 
                                                   psill=ozone_fit$psill, range=ozone_fit$range))
}

# find the model with the lowest RMSE
best_krige <- results_krige[which.min(results_krige$RMSE),]
cat("The best model is", best_krige$m, "with an RMSE of", best_krige$RMSE,
    "psill =", best_krige$psill, "range =", best_krige$range)
results_krige
# model     RMSE    psill    range
# Exc 0.002712457 0.0000061333824    2234.319
# but looking at the data, it seems like linear would fit the best

#krige
ozone_fit <- vgm(model="Lin", psill=0.000004, range=300000, nugget = 0.000005)
plot(ozone_vgm, ozone_fit, main = "Ozone Variogram (Linear)", 
     sub = "RMSE = 0.002  psill = 0.000004  range = 300000  nugget = 0.000005",
     xlab = "Distance (m)",
     ylab = "Semivariance") # plot the sample values, along with the fit model
LOOCV <- krige.cv(Arithmetic.Mean~1, ozone_spdf, model=ozone_fit)
RMSE <- sqrt(mean(LOOCV@data$residual^2))
# RMSE 0.002389295800 is lower than Exc
ozone_ord_krige <- krige(Arithmetic.Mean~1, ozone_spdf, grid_predict, model=ozone_fit)

# check sum of squared error
attr(ozone_fit, "SSErr")
# NULL

#view 
spplot(ozone_ord_krige, "var1.pred")

# plot the variance - estimate of estimation error
spplot(ozone_ord_krige, "var1.var")

#plot the same way as the IDW 
ozone_krige_sf <- st_as_sf(ozone_ord_krige) %>% 
  cbind(. , st_coordinates(.))
# tiled (filled) interpolation, displayed in EPSG 102004 coordinate system
ggplot() + 
  geom_tile(data=ozone_krige_sf, aes(x=X, y=Y, fill=var1.pred)) +
  geom_sf(data = statesOutline, color="black", fill=NA) +
  coord_sf(crs=new_crs) + 
  labs(title="Kriging Interpolated Ozone Values", 
       subtitle="For New York, Pennsylvania, New Jersey, Delaware, Maryland in 2021",
       fill="Mean (ppm)") +
  xlab("Longitude")+
  ylab("Latitude")+
  north(statesOutline)+
  ggsn::scalebar(statesOutline, dist = 100, dist_unit = "km",transform = FALSE, model = 'NAD83', st.size = 2, border.size= 0.5, location="topleft")

# get summary
summary(ozone_krige_sf$var1.pred)


########## PM2.5 Interpolation
### PM2.5 IDW
pm_data <- air_data%>%
  filter(`Parameter Name` == "PM2.5 - Local Conditions")

# plot pm data
pm_data %>%
  as.data.frame() %>%
  ggplot(aes(Longitude, Latitude)) + 
  geom_point(color="blue", alpha=3/4) + 
  ggtitle("PM2.5 - Local Conditions Concentration") + 
  coord_equal() + 
  theme_bw()

# convert to sf object and transform to EPSG 102004 (NAD 1983 Lambert contiguous USA)

crs = st_crs("EPSG:4326")
new_crs = st_crs("ESRI:102004")
pm_df_sf <- st_as_sf(pm_data, coords = c("Longitude", "Latitude"), crs = crs)%>%
  st_transform(new_crs)

# plot points
ggplot() + 
  geom_sf(data=pm_df_sf) +
  geom_sf(data = statesOutline, color="black", fill=NA) +
  coord_sf(crs=new_crs)

# grid already exists in this file
plot(grid_predict)

## get coordinates and value to interpolate from sf object
data <- pm_df_sf %>% 
  cbind(. , st_coordinates(.)) %>% 
  as.data.frame() %>% 
  dplyr::select(X,Y, `Arithmetic.Mean`)

## perform the IDW and join back the coordinates of the grid
pm_i <- phylin::idw(data[,3], 
                      data[,1:2],
                      grd) %>% 
  cbind(grd)

# Cross validation (Leave one out)

# convert sf object to sp object
pm_spdf <- as(pm_df_sf, "Spatial")


# This will give a plot of the RMSE values for each K value and the value 
# of K with the lowest RMSE will be indicated with a red line (run the code). 
# We can use the value of K with the lowest RMSE for our final model.

# Let's find the best k/p value
# create an empty data frame to store p values and RMSE values
r <- data.frame(p = integer(), RMSE = numeric())

test_val <- c(0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

# loop through different p values and perform cross-validation
for (p in test_val) {
  LOOCV <- krige.cv(formula=Arithmetic.Mean~1, locations=pm_spdf, nfold=nrow(pm_spdf), set = list(idp = p))
  
  # calculate RMSE
  RMSE <- sqrt(mean(LOOCV@data$residual^2))
  # add p and RMSE values to the results data frame
  r <- rbind(r, data.frame(p = p, RMSE = RMSE))
}

# find the p value with the lowest RMSE
best_idw <- r[which.min(r$RMSE),]
cat("The best p value is", best_idw$p, "with an RMSE of", best_idw$RMSE)
# The best p value is 2 with an RMSE of 1.101714
# result is 2, so we keep the previous result

ggplot(r, aes(x = p, y = RMSE)) + 
  geom_line() + 
  geom_point() +
  scale_x_continuous(breaks = seq(1, 20, by = 1)) + 
  labs(x = "p value", y = "RMSE", title="RMSE for test p values for PM2.5")

# rename Z to Mean
pm_i <- pm_i %>% mutate(Mean = Z) %>% dplyr::select(-Z)

#get summary
summary(pm_i)

# tiled (filled) interpolation, displayed in EPSG 102004 coordinate system
ggplot() + 
  geom_tile(data=pm_i, aes(x=X, y=Y, fill=Mean)) +
  geom_sf(data = statesOutline, color="black", fill=NA) +
  coord_sf(crs=new_crs) +
  labs(title="IDW Interpolated PM2.5 Values",
       subtitle="For New York, Pennsylvania, New Jersey, Delaware, Maryland in 2021",
       fill="Mean (μg/cubic meter)") +
  xlab("Longitude")+
  ylab("Latitude")+
  north(statesOutline)+
  ggsn::scalebar(statesOutline, dist = 100, dist_unit = "km",transform = FALSE, model = 'NAD83', st.size = 2, border.size= 0.5, location="topleft")



### PM2.5 Kriging

# assumption 1: data is normal
hist(pm_data$`Arithmetic Mean`)
shapiro.test(pm_data$`Arithmetic Mean`)
# It looks normal and the p value > 0.05

# assumption 2: mean and variance are equal

# Find extent values
# Using functions in the raster package

area <- extent(pm_spdf@bbox)
area

# Create a raster of 5 x 5 covering the area
pm_raster <- raster(nrows = 4, ncols = 4, ext = area,
                    vals = -999)
# Convert Raster to Polygons
polygons_area <- rasterToPolygons(pm_raster)
crs(polygons_area) <- crs(pm_spdf)

# Removing all columns except arithmetic mean
pm_spdf_mean <- pm_spdf[,26]
# Number of obserations per polygon
assessment <- over(polygons_area, pm_spdf_mean)
# Mean of values
assessment$mean <- over(polygons_area, pm_spdf_mean,
                        fn = mean, na.rm = T)[,1]
# Variance of Values
assessment$var <- over(polygons_area, pm_spdf_mean,
                       fn = var, na.rm = T)[,1]
# Change Names
names(assessment)[1] <- c("counts")
# Replace Polygon Data
polygons_area@data <- assessment
summary(polygons_area@data)
spplot(polygons_area, c("counts", "mean", "var"), names.attr = c("Count","Mean", "Variance"),
       colorkey=list(space="bottom"), scales = list(draw = TRUE),
       main = "Stationarity Test",
       as.table = TRUE)

summary(polygons_area)

# assumption 3: isotrophy
# Look at IDW plot again to see if there's any spatial trend
ggplot() + 
  geom_tile(data=pm_i, aes(x=X, y=Y, fill=Mean)) +
  geom_sf(data = statesOutline, color="white", fill=NA) +
  coord_sf(crs=new_crs)

pm_vgm <- variogram(Arithmetic.Mean~1, data = pm_spdf)

# see all types of variograms 
show.vgms()

pm_fit <- fit.variogram(pm_vgm, model = vgm(model="Sph"))
plot(pm_vgm, pm_fit) # plot the sample values, along with the fit model

pm_fit

# psill: 1.227854
# range: 29670.48

# Performing kriging
# for the points grid, we already made that in the previous step for IDW

#krige
pm_ord_krige <- krige(Arithmetic.Mean~1, pm_spdf, grid_predict, model=pm_fit)

#view 
spplot(pm_ord_krige, "var1.pred")

# plot the variance - estimate of estimation error
spplot(pm_ord_krige, "var1.var")

# cross validation
# create array of models to test
test_models = c("Nug", "Exp", "Sph", "Gau", "Exc", "Mat", "Ste", "Cir",
                "Lin", "Bes", "Pen", "Per", "Wav", "Hol", "Log", "Pow", "Spl")

# create array to store results
results_krige <- data.frame()

# loop through different models and perform cross-validation
for (m in test_models) {
  pm_fit <- fit.variogram(pm_vgm, model = vgm(model=m))
  LOOCV <- krige.cv(Arithmetic.Mean~1, pm_spdf, model=pm_fit)
  
  # calculate RMSE
  RMSE <- sqrt(mean(LOOCV@data$residual^2))
  # add p and RMSE values to the results data frame
  results_krige <- rbind(results_krige, data.frame(model=m, RMSE = RMSE, 
                                                   psill=pm_fit$psill, range=pm_fit$range))
}

# find the model with the lowest RMSE
best_krige <- results_krige[which.min(results_krige$RMSE),]
cat("The best model is", best_krige$m, "with an RMSE of", best_krige$RMSE,
    "psill =", best_krige$psill, "range =", best_krige$range)
results_krige
# model     RMSE    psill    range
# Exc  0.9940905 1.95995239  59588.793135

# retry kriging with this model

#krige
# first refit the model manually to try to decrease the RMSE
pm_fit <- vgm(model="Exc", psill = 2, range = 150000, nugget = 0.3)
plot(pm_vgm, pm_fit, main = "PM2.5 Variogram (Exponential class)", 
     sub = "RMSE = 0.979  psill = 2  range = 150000  nugget = 0.3",
     xlab = "Distance (m)",
     ylab = "Semivariance") # plot the sample values, along with the fit model
LOOCV <- krige.cv(Arithmetic.Mean~1, pm_spdf, model=pm_fit)
RMSE <- sqrt(mean(LOOCV@data$residual^2))
#RMSE 0.978716099

pm_ord_krige <- krige(Arithmetic.Mean~1, pm_spdf, grid_predict, model=pm_fit)

# check sum of squared error
attr(pm_fit, "SSErr")
# NULL

#view 
spplot(pm_ord_krige, "var1.pred")

# plot the variance - estimate of estimation error
spplot(pm_ord_krige, "var1.var")

#plot the same way as the IDW 
pm_krige_sf <- st_as_sf(pm_ord_krige) %>% 
  cbind(. , st_coordinates(.))
# tiled (filled) interpolation, displayed in EPSG 102004 coordinate system
ggplot() + 
  geom_tile(data=pm_krige_sf, aes(x=X, y=Y, fill=var1.pred)) +
  geom_sf(data = statesOutline, color="black", fill=NA) +
  coord_sf(crs=new_crs) + 
  labs(title="Kriging Interpolated PM2.5 Values", 
       subtitle="For New York, Pennsylvania, New Jersey, Delaware, Maryland in 2021",
       fill="Mean (μg/cubic meter)") +
  xlab("Longitude")+
  ylab("Latitude")+
  north(statesOutline)+
  ggsn::scalebar(statesOutline, dist = 100, dist_unit = "km",transform = FALSE, model = 'NAD83', st.size = 2, border.size= 0.5, location="topleft")

# get summary
summary(pm_krige_sf$var1.pred)




########## NO2 Clustering
# Clustering was performed in GeoDa. This is to visualize the result
# Reading in csv exported from GeoDa
no2_clusters <- read_csv("data/no2_clusters_v3.csv", show_col_types = FALSE)
states <- sf::st_read("data/states.shp")

state_names <- c("Delaware", "Maryland", "New Jersey", "New York", "Pennsylvania")
states_filtered <- states %>% filter(NAME %in% state_names)
state_centroids <- st_centroid(states_filtered, of_largest_polygon = TRUE) %>% 
  cbind(st_coordinates(.)) %>%
  as.data.frame() %>%
  dplyr::select(NAME, X, Y)

# Count the number of points per cluster
cluster_counts <- no2_clusters %>%
  group_by(CL) %>%
  summarize(count = n()) %>%
  mutate(cluster_label = paste("Cluster", CL, "(", count, ")", sep = " "))

# Define individual nudges for each state
state_centroids <- state_centroids %>%
  mutate(nudge_x = case_when(
    NAME == "Delaware" ~ 0.05,
    NAME == "Maryland" ~ 0.05,
    NAME == "New Jersey" ~ -0.05,
    NAME == "New York" ~ -0.05,
    NAME == "Pennsylvania" ~ 0.05
  ),
  nudge_y = case_when(
    NAME == "Delaware" ~ 0.05,
    NAME == "Maryland" ~ -0.05,
    NAME == "New Jersey" ~ 0.05,
    NAME == "New York" ~ 0.05,
    NAME == "Pennsylvania" ~ -0.05
  ))


# Calculate min and max arithmetic mean for each cluster
cluster_ranges <- no2_clusters %>%
  group_by(CL) %>%
  summarize(min_mean = min(`Arithmetic Mean`),
            max_mean = max(`Arithmetic Mean`)) %>%
  mutate(cluster_label = paste("Cluster", CL, ": ", round(min_mean, 3), " - ", round(max_mean, 3), sep = " "))


custom_palette <- c(brewer.pal(8, "Set2"), "#FDB813", "#B2DF8A")

# Plot the data using ggplot2 with the states background and state names
ggplot() +
  geom_sf(data = states_filtered, fill = "#eeeeee", color = "black", size = 0.1) +
  #geom_text_repel(data = state_centroids, aes(x = X, y = Y, label = NAME, nudge_x = nudge_x, nudge_y = nudge_y), size = 2.45, fontface = "bold") +
  geom_point(data = no2_clusters, aes(x = Longitude, y = Latitude, color = factor(CL))) +
  labs(title = "Spatial Clustering of NO2 Monitoring Stations",
       subtitle = paste("States:", paste(state_names, collapse = ", ")),
       x = "Longitude",
       y = "Latitude",
       color = "Cluster") +
  theme_minimal() +
  scale_color_manual(values = custom_palette, name = "Clusters with Pollution Range (PPB)", labels = cluster_ranges$cluster_label)

# Histogram
cluster_counts <- no2_clusters %>%
  group_by(CL) %>%
  summarize(count = n()) %>%
  mutate(cluster_label = factor(CL))

# Plot the bar plot using ggplot2
ggplot(cluster_counts, aes(x = cluster_label, y = count, fill = cluster_label)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.9) +
  labs(title = "Number of NO2 Monitors per Cluster",
       x = "Clusters",
       y = "Number of NO2 Monitors",
       fill = "Cluster") +
  theme_minimal() +
  scale_fill_manual(values = custom_palette, name = "Clusters with Pollution Range (PPB)", labels = cluster_ranges$cluster_label)

# Descriptive statistics
no2_clusters

no2_clusters %>%
  filter(CL == 1) %>%
  dplyr::select(`Arithmetic Mean`) %>%
  summary()

no2_clusters %>%
  filter(CL == 2) %>%
  dplyr::select(`Arithmetic Mean`) %>%
  summary()


no2_clusters %>%
  filter(CL == 3) %>%
  dplyr::select(`Arithmetic Mean`) %>%
  summary()

no2_clusters %>%
  filter(CL == 4) %>%
  dplyr::select(`Arithmetic Mean`) %>%
  summary()

########## Ozone Clustering
# Clustering was performed in GeoDa. This is to visualize the result
# Reading in csv exported from GeoDa
ozone_clusters <- read_csv("data/ozone_clusters_v2.csv", show_col_types = FALSE)

# Include the number of clusters

# Count the number of points per cluster
cluster_counts <- ozone_clusters %>%
  group_by(CL) %>%
  summarize(count = n()) %>%
  mutate(cluster_label = paste("Cluster", CL, "(", count, ")", sep = " "))

# State centroids already defined in this file

# Calculate min and max arithmetic mean for each cluster
cluster_ranges <- ozone_clusters %>%
  group_by(CL) %>%
  summarize(min_mean = min(`Arithmetic Mean`),
            max_mean = max(`Arithmetic Mean`)) %>%
  mutate(cluster_label = paste("Cluster", CL, ": ", round(min_mean, 3), " - ", round(max_mean, 3), sep = " "))

custom_palette <- c(brewer.pal(8, "Set2"), "#FDB813", "#B2DF8A")

# Plot the data using ggplot2 with the states background and state names
ggplot() +
  geom_sf(data = states_filtered, fill = "#eeeeee", color = "black", size = 0.1) +
  #geom_text_repel(data = state_centroids, aes(x = X, y = Y, label = NAME, nudge_x = nudge_x, nudge_y = nudge_y), size = 2.45, fontface = "bold") +
  geom_point(data = ozone_clusters, aes(x = Longitude, y = Latitude, color = factor(CL))) +
  labs(title = "Spatial Clustering of Ozone Monitoring Stations",
       subtitle = paste("States:", paste(state_names, collapse = ", ")),
       x = "Longitude",
       y = "Latitude",
       color = "Cluster") +
  theme_minimal() +
  scale_color_manual(values = custom_palette, name = "Clusters with Pollution Range (PPM)", labels = cluster_ranges$cluster_label)

# Create a data frame of the count of observations in each cluster
cluster_counts <- ozone_clusters %>%
  group_by(CL) %>%
  summarize(count = n()) %>%
  mutate(cluster_label = factor(CL))

# Plot the bar plot using ggplot2
ggplot(cluster_counts, aes(x = cluster_label, y = count, fill = cluster_label)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  labs(title = "Number of Ozone Monitors per Cluster",
       x = "Clusters",
       y = "Number of Ozone Monitors",
       fill = "Cluster") +
  theme_minimal() +
  scale_fill_manual(values = custom_palette, name = "Clusters with Pollution Range (PPM)", labels = cluster_ranges$cluster_label)


# Lets generate some statistics 
ozone_clusters

# For each cluster 
ozone_clusters %>%
  filter(CL == 1) %>%
  dplyr::select(`Arithmetic Mean`) %>%
  summary()

ozone_clusters %>%
  filter(CL == 2) %>%
  dplyr::select(`Arithmetic Mean`) %>%
  summary()

ozone_clusters %>%
  filter(CL == 3) %>%
  dplyr::select(`Arithmetic Mean`) %>%
  summary()

ozone_clusters %>%
  filter(CL == 4) %>%
  dplyr::select(`Arithmetic Mean`) %>%
  summary()

ozone_clusters %>%
  filter(CL == 5) %>%
  dplyr::select(`Arithmetic Mean`) %>%
  summary()

ozone_clusters %>%
  filter(CL == 6) %>%
  dplyr::select(`Arithmetic Mean`) %>%
  summary()

ozone_clusters %>%
  filter(CL == 7) %>%
  dplyr::select(`Arithmetic Mean`) %>%
  summary()

ozone_clusters %>%
  filter(CL == 8) %>%
  dplyr::select(`Arithmetic Mean`) %>%
  summary()

ozone_clusters %>%
  filter(CL == 9) %>%
  dplyr::select(`Arithmetic Mean`) %>%
  summary()

ozone_clusters %>%
  filter(CL == 10) %>%
  dplyr::select(`Arithmetic Mean`) %>%
  summary()

########## PM2.5 Clustering
# Clustering was performed in GeoDa. This is to visualize the result
# Reading in csv exported from GeoDa
pm2_clusters <- read_csv("data/pm2_clusters.csv", show_col_types = FALSE)

# Count the number of points per cluster
cluster_counts <- pm2_clusters %>%
  group_by(CL) %>%
  summarize(count = n()) %>%
  mutate(cluster_label = paste("Cluster", CL, "(", count, ")", sep = " "))

# State centroids already defined in this file

# Calculate min and max arithmetic mean for each cluster
cluster_ranges <- pm2_clusters %>%
  group_by(CL) %>%
  summarize(min_mean = min(`Arithmetic Mean`),
            max_mean = max(`Arithmetic Mean`)) %>%
  mutate(cluster_label = paste("Cluster", CL, ": ", round(min_mean, 3), " - ", round(max_mean, 3), sep = " "))

custom_palette <- c(brewer.pal(8, "Set2"), "#FDB813", "#B2DF8A")

# Plot the data using ggplot2 with the states background and state names
ggplot() +
  geom_sf(data = states_filtered, fill = "#eeeeee", color = "black", size = 0.1) +
  #geom_text_repel(data = state_centroids, aes(x = X, y = Y, label = NAME, nudge_x = nudge_x, nudge_y = nudge_y), size = 2.45, fontface = "bold") +
  geom_point(data = pm2_clusters, aes(x = Longitude, y = Latitude, color = factor(CL))) +
  labs(title = "Spatial Clustering of PM2.5 Monitoring Stations",
       subtitle = paste("States:", paste(state_names, collapse = ", ")),
       x = "Longitude",
       y = "Latitude",
       color = "Cluster") +
  theme_minimal() +
  scale_color_manual(values = custom_palette, name = "Clusters with Pollution Range (μg/m3)", labels = cluster_ranges$cluster_label)

# Histogram

cluster_counts <- pm2_clusters %>%
  group_by(CL) %>%
  summarize(count = n()) %>%
  mutate(cluster_label = factor(CL))

# Plot the bar plot using ggplot2
ggplot(cluster_counts, aes(x = cluster_label, y = count, fill = cluster_label)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.9) +
  labs(title = "Number of PM2.5 Monitors per Cluster",
       x = "Clusters",
       y = "Number of PM2.5 Monitors",
       fill = "Cluster") +
  theme_minimal() +
  scale_fill_manual(values = custom_palette, name = "Clusters with Pollution Range (μg/m3)", labels = cluster_ranges$cluster_label)


# Descriptive statistic for each cluster

pm2_clusters

pm2_clusters %>%
  filter(CL == 1) %>%
  dplyr::select(`Arithmetic Mean`) %>%
  summary()

pm2_clusters %>%
  filter(CL == 2) %>%
  dplyr::select(`Arithmetic Mean`) %>%
  summary()

pm2_clusters %>%
  filter(CL == 3) %>%
  dplyr::select(`Arithmetic Mean`) %>%
  summary()

pm2_clusters %>%
  filter(CL == 4) %>%
  dplyr::select(`Arithmetic Mean`) %>%
  summary()

pm2_clusters %>%
  filter(CL == 5) %>%
  dplyr::select(`Arithmetic Mean`) %>%
  summary()

pm2_clusters %>%
  filter(CL == 6) %>%
  dplyr::select(`Arithmetic Mean`) %>%
  summary()

