
library(tidyverse)
library(sf)
library(phylin)
library(dplyr)
library(automap)
library(sp)
library(gstat)
library(ggsn)

# set working directory with setwd(dir)

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

# plot points
ggplot() + 
  geom_sf(data=no2_df_sf) +
  geom_sf(data = statesOutline, color="black", fill=NA) +
  coord_sf(crs=new_crs)

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
i <- phylin::idw(data[,3], 
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
i <- phylin::idw(data[,3], 
                 data[,1:2],
                 grd, p=1) %>% 
  cbind(grd)


# rename Z to Mean
i <- i %>% mutate(Mean = Z) %>% dplyr::select(-Z)


#get summary
summary(i)

# tiled (filled) interpolation, displayed in EPSG 102004 coordinate system
ggplot() + 
  geom_tile(data=i, aes(x=X, y=Y, fill=Mean)) +
  geom_sf(data = statesOutline, color="black", fill=NA) +
  coord_sf(crs=new_crs) +
  labs(title="IDW Interpolated NO2 Values", 
       subtitle="For New York, Pennsylvania, New Jersey, Delaware, Maryland in 2021",
       fill="Mean (ppb)") +
  xlab("Longitude")+
  ylab("Latitude")+
  north(statesOutline)+
  ggsn::scalebar(statesOutline, dist = 100, dist_unit = "km",transform = FALSE, model = 'NAD83', st.size = 2, border.size= 0.5, location="topleft")






# Kriging

# assumption 1: data is normal
hist(no2_data$`Arithmetic Mean`)
shapiro.test(no2_data$`Arithmetic Mean`)
# p-value = 0.2284, so it is normal

# assumption 2: mean and variance are equal

# Find extent values
# Using functions in the raster package
library(raster)
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
# mean is even


# assumption 3: isotrophy
# Look at IDW plot again to see if there's any spatial trend
ggplot() + 
  geom_tile(data=i, aes(x=X, y=Y, fill=Mean)) +
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
plot(grd_sp)
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
  scalebar(statesOutline, dist = 100, dist_unit = "km",transform = FALSE, model = 'NAD83', st.size = 2, border.size= 0.5, location="topleft")


# get summary
summary(no2_krige_sf$var1.pred)

