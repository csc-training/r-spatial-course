#####################################
##### Session 1 - Raster basics #####
#####################################
require(raster)
require(viridis) # Optional, this library includes good-looking color pallettes

# Create a raster from scracth, no crs, no geographic extent!
?raster
x <- 1:3
y <- 4:6
z <- 7:9
r1 <- raster(rbind(x,y,z))

# Create another raster, this time from a matrix
m <- matrix(1:25, ncol=5, nrow=5, byrow=T) # Create a 5 x 5 matrix
r2 <- raster(m) 
plot(r2, col=viridis(25)) 

# Create a raster from scracth, CRS=WGS84 and whole globe as default extent
# cell size 1*1
r3 <- raster() # only raster geometry, no values!
hasValues(r3) # Logical test, whether RasterLayer is associated with values
values(r3) <- 1:ncell(r3) # "Fill" raster with values 
plot(r3, col=inferno(50))

# Create a raster with cell size of 50
r4 <- raster(res=50)
values(r4) <- runif(ncell(r4)) # Fill RasterLayer with random values from uniform distribution
plot(r4, col=magma(40))

# Create a raster from scracth, defined CRS and extent 
m2 <- matrix(runif(200), ncol=20, nrow=10, byrow=T) # 10 x 20 matrix with random values
r5 <- raster(m2, xmn=-180, xmx=180, ymx=90, ymn=30,
             crs=CRS("+init=epsg:4326"))
plot(r5, col=rainbow(50), ylab="Latitude", xlab="Longitude")

# Reading in raster files
r6 <- raster("D:\\Opetus\\CSC2018\\Data\\Tmean_Monthly_1961-2014\\Tmon_19610101.tif")
plot(r6, main="Raster from file")
summary(r6) # Basic summary
extent(r6) # Print geographical extent of the raster
bbox(r6) # Same thing, slightly different presentation
cellStats(r6, stat = mean, na.rm=T) # Obtain a mean over the raster
freq(r6) # Frequency table
rvals <- getValues(r6) # Extract pixel values and stored them as a vector
hist(rvals, main="Raster values") # Plot a histogram of the raster values
df <- as.data.frame(r6, xy=T) # Convert RasterLayer to a Data Frame
plot(df$x, df$y, pch=19, cex=0.4)
df <- df[complete.cases(df), ] # Extract rows with values
points(df$x, df$y, pch=19, cex=0.6, col="red")

### Creating a raster stack ###
# Import two rasters
r7 <- raster("D:\\Opetus\\CSC2018\\Data\\Tmean_Monthly_1961-2014\\Tmon_19610201.tif")
r8 <- raster("D:\\Opetus\\CSC2018\\Data\\Tmean_Monthly_1961-2014\\Tmon_19610301.tif")
st <- stack(r6, r7, r8) # Stack three rasters
plot(st) # Plot all layers in the stack
summary(st) # Print summaries of the layers
names(st) <- c("Tmon1", "Tmon2", "Tmon3") # Change the layers' names

### Creating raster stack in a loop ###
st2 <- stack() # Create an empty raster stack
setwd("D:\\Opetus\\CSC2018\\Data\\Tmean_Monthly_1961-2014\\") # Define working directory
files <- list.files()[1:12] # List first 12 rasters in the folder

# Loop across the files and add to the stack
for (f in files) {
  print(f)
  st2 <- stack(st2, f)
}

# Or use lapply to import the rasters at once as a list, and then create a stack
st3 <- do.call(stack, lapply(files, FUN=raster))

# Form a new raster stack
new_st <- st2[[c(1:4, 6)]] # Extract layers 1-4 & 6

# Write raster file out 
?writeRaster
writeRaster(r6, filename = "F:\\Testi\\Rasteri.tif", format="GTiff")
?dataType

##########################
##### Assignment 1.1 #####
##########################
# a) Create a raster object covering roughly the spatial extent of Finland (in lat lon)
# at the spatial resolution of 0.1
# b) Assign any values to the raster and plot it
# c) Write the raster out as a geotiff

##########################
##### Assignment 1.2 #####
##########################
# a) Read in the raster "Tmon_20100101.tif" and plot it
# b) Plot a histogram of the raster values, 
# c) include a red vertical line to the histogram indicating the mean of distribution

##########################
##### Assignment 1.3 #####
##########################
# a) read in three raster layers ("Tmon_20100101.tif", "Tmon_20100201.tif" and "Tmon_20100301.tif")
# and create a raster stack
# b) calculate mean values of each raster
# c) Plot first and third raster in the stack

###########################################
##### Session 2 - Raster manipulation #####
###########################################
rm(list=ls()) # Clear the workspace, DELETES EVERY OBJECT!

### Reprojection ###
# www.spatialreference.org
r <- raster("D:\\Opetus\\CSC2018\\Data\\Tmean_Monthly_1961-2014\\Tmon_19610201.tif")
proj4string(r) # Prints out the CRS string (here Euref-Fin)

# Reproject for Euref-Fin to ETRS-LAEA using bilinear resampling
r2 <- projectRaster(from=r, crs = CRS("+init=epsg:3035"), res=10000, method="bilinear")
proj4string(r2)

# Reproject for Euref-Fin to WGS84 using bilinear resampling (note resolution)
r3 <- projectRaster(from=r, crs = CRS("+init=epsg:4326"), res=0.1, method="bilinear")
proj4string(r3)

# Plot the three rasters with different projections
par(mfrow=c(1,3))
plot(r, main="Euref-Fin", col=magma(50))
plot(r2, main="ETRS-LAEA", col=plasma(50))
plot(r3, main="WGS84", col=viridis(50))

### Making spatial subsets i.e. cropping ###
?crop
par(mfrow=c(1,1)) # Return to normal plotting settings
plot(r3)
e <- drawExtent() # draw extent manually
cr <- crop(r3, e) # Perform the crop
plot(cr)

# Manual entry
e <- extent(20, 28, 66, 70) # order: xmin, xmax, ymin, ymax
plot(r3); plot(e, add=T, lwd=2)
plot(r3, ext=e)
cr <- crop(r3, e) # perform crop
plot(cr2, breaks=seq(-14, -5, length.out = 5), col=inferno(5),
     main="Lapland")

### Merge rasters ###
plot(r3); plot(e, add=T, lwd=2)
e2 <- extent(22, 30, 62, 66) # Define an another area to be cropped
plot(e2, add=T, lwd=2)

cr3 <- crop(r3, e2) # perform crop
m <- merge(cr, cr3) # merge the two cropped rasters as one
plot(m) # Plot the merged rasters

### Spatial aggregation ### 
?aggregate
plot(r, main="Resolution 10 x 10 km")
a <- aggregate(r, fact=2, fun=mean) # Increase the cell size by factor of two, based on mean function
plot(a, main="Resolution 20 x 20 km")
a2 <- aggregate(r, fact=5, fun=max) # Increase the cell size by factor of two, based on max function
plot(a2, main="Resolution 50 x 50 km")

### Moving window statistics ###
?focal
plot(cr)
f <- focal(cr, w=matrix(1,5,5), fun=mean, na.rm=T) # Focal mean filter of 5 x 5 pixels
plot(f, main="Focal mean 5 x 5 pixels")
f2 <- focal(cr, w=matrix(1,5,5), fun=min, na.rm=T) # Focal minimum filter of 5 x 5 pixels
plot(f2, main="Focal minimum 5 x 5 pixels")

# Compare air temperature at a cell to its surroundings 
plot(cr-f) # with focal mean
plot(cr-f2) # with focal minimum

### Reclassification ###
# Accessing cell values
plot(r)
values(r)[values(r) < -5] <- -5  # Change values less than -5 to -5
plot(r)

plot(r2)
values(r2)[values(r2) < 0 & values(r2) > -5] <- NA # Change values less than 0 AND greater than -5 to NA  
plot(r2)

?reclassify
r <- raster("D:\\Opetus\\CSC2018\\Data\\Tmean_Monthly_1961-2014\\Tmon_19610201.tif")
plot(r)

# All values < 0 are recoded to 0, all values > 0 are recode to 1
rcl <- reclassify(r, c(-Inf, 0, 0,  0, Inf, 1)) 
plot(rcl, col=c("blue", "red"))

# All values < -5 are recoded to 1, 
# all values from -5 to 0 are recoded to 2,
# all values > 0 are recoded to 3
rcl2 <- reclassify(r, c(-Inf,-5, 1,  
                        -5, 0, 2,  
                        0, Inf, 3))
plot(rcl2, col=viridis(3))

##########################
##### Assignment 2.1 #####
##########################
# a) Work with the previous reclassified air temperature object ("rcl2"); 
# re-project the raster to Finland Uniform Coordinate System (Zone 3, a.k.a YKJ) using
# nearest-neighbourhood interpolation at 20 km spatial resolution (see ?projectRaster)

# b) Read in the June 2000 air temperature raster "Tmon_20000601.tif" 
# and reproject it to WGS84 (spatial resolution 0.1 degrees).
# Crop the raster to cover areas north of Arctic Circle
# c) For the cropped raster in B, apply a moving average of 3x3 pixels

# BONUS d) For the cropped raster in b, write a function that calculates local 95 % range 
# on variation (i.e. 97.5 percentile - 2.5 percentile) inside a 5x5 neighbourhood

###################################
##### Session 3 - Map algebra #####
###################################
### Raster calculations 
r <- raster("D:\\Opetus\\CSC2018\\Data\\Tmean_Monthly_1961-2014\\Tmon_19610201.tif")
plot(r)
r2 <- r * 10 # Multiply all cell values by 10
r3 <- round(r2/10) # Divide all cell values by 10, and then round to integers
r3 <- calc(r2, fun=function(x) {x/10}) # Same thing, but using calc() -function

### Calculate annual mean air temperature for the year 2000 ###
# Define directory to the file location
setwd("D:\\Opetus\\CSC2018\\Data\\Tmean_Monthly_1961-2014\\")

# Filter data to include only year 2000
files <- list.files(pattern="2000")

# Create a raster stack from the files
st <- stack(files)

# Summarize the variation of each rasters
cellStats(st, "mean") # Mean 
cellStats(st, "range") # Range of variation

# Calculate mean over all the layers in the stack
avg <- mean(st)
avg <- calc(st, fun = mean)

### Raster overlay ###
# Calculate amount of snow precipitation in a month
# Read in an air temperature raster
temp <- raster("D:\\Opetus\\CSC2018\\Data\\Tmean_Monthly_1961-2014\\Tmon_19610401.tif")

# Read in a precipitation raster (Ideally the same year and month...)
rr <- raster("D:\\Opetus\\CSC2018\\Data\\RRmon_Monthly_1961-2014\\RRmon_19610401.tif")

# Reclassify the air temperature raster to identify cells with subzero temperatures
rcl <- reclassify(temp, c(-Inf,0, 1, 0.0001, Inf, 0)) # Air T < 0 is recoded to 1

# Multiply reclassified air temperature rasters with the precipitation raster
rr2 <- overlay(rcl, rr, fun=function(x,y){(x*y)} )
plot(rr2, main="Snow precipitation", col=rev(viridis(20)))

### DEMO ###
### Loop throuh the whole year and calculate snow precipitation sum ###
# List all air temperature rasters in the folder
files_t <- list.files(path="D:\\Opetus\\CSC2018\\Data\\Tmean_Monthly_1961-2014\\", 
                      full.names = T)

# List all precipitation rasters in the folder
files_rr <- list.files(path="D:\\Opetus\\CSC2018\\Data\\RRmon_Monthly_1961-2014\\",
                       full.names = T)

# Apply a file mask to find files for the year 2010
year <- "2010"

# Pattern matching ... 
t_mask <- files_t[substr(files_t, 
                         nchar(files_t)-11, nchar(files_t)-8)==year]
rr_mask <- files_rr[substr(files_rr, 
                         nchar(files_rr)-11, nchar(files_rr)-8)==year]

# Create an empty raster stack for the monthly snow precipitation rasters
st <- stack()

# Loop through the months (1-12) and calculate snow precipitation
for (i in 1:12) {
  print(paste("Month =", i))
  
  # Read in the n:th raster
  temp <- raster(t_mask[i])
  prec <- raster(rr_mask[i])
  
  # Reclassify air temperature raster
  rcl <- reclassify(temp, c(-Inf,0, 1, 0.0001, Inf, 0))
  
  # Overlay with precipitation raster
  ovrl <- overlay(rcl, prec, fun=function(x,y){(x*y)} )
  
  # Add calculated snow precipitation to the stack
  st <- stack(st, ovrl)
}

# Calculate pixel-wise sum across the whole year
snow <- calc(st, fun=sum)
plot(snow, main="Snow precipitation 2010", col=rev(viridis(20)))

### Pixel-wise regression ####
# List all January air temperature files in the data folder
setwd("D:\\Opetus\\CSC2018\\Data\\Tmean_Monthly_1961-2014\\")
files <- list.files()

# Extract only files representing March air temperatures
file_sub <- files[substr(files, 10, 11) == "03"]

# Create a raster stack to indlude all March air temperature rasters
st <- do.call(stack, lapply(file_sub, FUN=raster))

# Create a variable time (i.e. sequence from 1 to 54) indicating the number of years
time <- 1:nlayers(st) 

### Define functions to fit a model and extract relevant information ###
# For obtaining estimated model coefficients
fun_model_coeffs <- function(x) { if (is.na(x[1])){ NA } 
                  else { m <- lm(x ~ time); return(summary(m)$coefficients[2]) }}

# For obtaining p-values
fun_p_values <- function(x) { if (is.na(x[1])){ NA } 
                  else { m <- lm(x ~ time); return(summary(m)$coefficients[8]) }}

### Apply functions over the raster stack to fit pixel-wise regression ###
model_coeffs <- calc(st, fun_model_coeffs)# Returns estimated regression slopes
model_p_values <- calc(st, fun_p_values) # Return p-values

# Plot the results
plot(model_coeffs, main="Estimated regression slope", col=plasma(20))
plot(model_p_values, breaks = c(0, 0.05, 1), 
    col=c("pink3", "lightblue"), main="Significance of the trend, red areas p<0.05")

# Convert p-value raster to basic data frame
pvals <- as.data.frame(model_p_values, xy=T)

# Get rid of missing values
pvals <- pvals[complete.cases(pvals), ]

# Subset p-value data frame to hold values < 0.05 indicating statistical significance
sub <- subset(pvals, layer <= 0.05)

# Overlay points with p<0.05 on the top of the regression slope raster
plot(model_coeffs, main="Estimated regression slope", col=plasma(20))
points(sub$x, sub$y, pch=19, cex=.3)

##########################
##### Assignment 3.1 #####
##########################
# Convert an air temperature raster (no matter which of those) from celsius to kelvins 
# and write the new raster out as a geotiff to a new directory

##########################
##### Assignment 3.2 #####
##########################
# Can you do the same for the first 10 air temperature rasters in your directory? 

##########################
##### Assignment 3.3 #####
##########################
# a) Create a RasterStack that holds all precipitation rasters for 1991.   
# b) Now calculate precipitation sum for summer months (Jun-Aug)

##########################
##### Assignment 3.4 #####
##########################
# a) Investigate, whether there has been a significant trend (time period of 1961-2014)
# in May water precipitation. Plot the estimated coefficients and p-values. 
# b) Summarize the variation in regression coefficients over the whole finland using 
# (i) mean and (ii) 95 % range of variation.

###########################################################
##### SESSION 4 - Raster package in spatial modelling #####
###########################################################
# Draw random pixel values from a raster
# sampleRandom()
r <- raster("D:\\Opetus\\CSC2018\\Data\\Tmean_Monthly_1961-2014\\Tmon_19950901.tif")
plot(r)

pts <- sampleRandom(r, size = 200, xy=T, sp=T) # Draw 200 random pixels
plot(pts, add=T, pch=19, cex=.5)

### Extract raster values to points ###
# First, import some additional environmental data 
elev <- raster("D:\\Opetus\\CSC2018\\Data\\Elevation_10km.tif") # Digital elevation model
lake <- raster("D:\\Opetus\\CSC2018\\Data\\Lake_10km.tif") # Lake cover
envs <- stack(elev, lake) # Stack the two rasters

ext <- extract(envs, pts, sp=T) # Extract raster values to the randomly chosen points

# Rename the variables in both point and RasterStack
# It's critical, that the variable names match exactly!
names(envs) <- c("elev", "lake")
names(ext) <- c("x", "y", "tmon", "elev", "lake")

# Create a multiple regression model to predict air temperatures
mod <- lm(tmon~x+y+elev, data=ext) # The notation means, that the variation in air temperatures 
                                   # are being explained with geographical location and elevation  
summary(mod) # Print regression model results (e.g. estimated coefficients, R-squared, etc.)

# Now based on the multiple regression model, make a prediction to the raster stack
pred <- interpolate(object = envs, model=mod, xyOnly=FALSE,  
                      xyNames=c("x", "y"), fun=predict)
plot(pred, main="Linear model prediction") # Plot the prediction

# DEMO: Fit a generalized additive model (GAM) and make a prediction
require(mgcv)

# This is a syntax for GAM model, with k-arguments defining the maximun
# smoothing function (here, roughly the same idea as 3rd order polynomials ...)
mod2 <- gam(tmon ~ s(elev, k=3) + te(x, y, k=3), data=ext, family = gaussian)
summary(mod2) # Print GAM model results

# Make a prediction based on the fitted GAM model
pred2 <- interpolate(object = envs, model=mod2, xyOnly=FALSE,  
                    xyNames=c("x", "y"), fun=predict.gam)
plot(pred2, main="GAM prediction")

### Going parallel ###
require(snow)

# Set up a CPU cluster
beginCluster(2)

# Or detect number of cores available
require(parallel)
ncores <- detectCores()

# Make a GAM prediction to the raster stack using multicore calculations
pred <- clusterR(envs, raster::interpolate,  
                 args=list(model=mod2, xyOnly=FALSE,  
                           xyNames=c("x", "y"), fun=predict.gam, type="response"))

# End the local cluster
endCluster()

##########################
##### Assignment 4.1 #####
##########################
# Update the previous linear regression model (tmon ~ x + y + elev)
# by adding variable "lake" as a predictor. 
# Make a prediction, plot the results and write out the prediction as a geotiff




