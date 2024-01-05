#packages -------------------------------------------------------------------------------------
library(tidyverse)
library(stats)
library(amt)
library(sf)
library(stars)
library(mapview)
library(terra)
library(raster)
library(lubridate)
library(adehabitatHR)
library(rgdal)
library(leaflet)
library(leafem)
library(leafpop)
library(RColorBrewer)
library(leafsync)
library(ggplot2)
library(ggmap)
library(ggspatial)
library(data.table)
library(nasapower)
library(MODISTools)
library(geosphere)
library(move)
library(utils)
library(lme4)
library(GGally)
library(xtable)
library(viridis)
library(sp)
library(ggpubr)
library(leafsync)
library(coxme)

#Functions ---------------------------------------------------------------------------
landclass <- function(x) {
  if (x==1){
    "Deciduous woodland" 
  } else if (x== 2){
    "Coniferous woodland"
  } else if (x==3) {
    "Arable"
  } else if (x==4) {
    "Improve grassland" 
  } else if (x==5) {
    "Neutral grassland"
  } else if (x==6) {
    "Calcareous grassland"
  } else if (x==7) {
    "Acid grassland"
  } else if (x==8) {
    "Fen"
  } else if (x==9) {
    "Heather"
  } else if (x==10) {
    "Heather grassland"
  } else if (x==11) {
    "Bog"
  } else if (x==12) {
    "Inland rock"
  } else if (x==13) {
    "Saltwater"
  } else if (x==14) {
    "Freshwater"
  } else if (x==15) {
    "Superlittoral rock"
  } else if (x==16) {
    "Superlittoral sediment"
  } else if (x==17) {
    "Lutteroal rock"
  } else if (x==18) {
    "Littoral sediment"
  } else if (x==19) {
    "Saltmarsh"
  } else if (x==20) {
    "Urban"
  } else if (x==21) {
    "Suburban"
  }
}

landclass2 <- function(x) {
  ifelse(x==1,"Deciduous woodland",
  ifelse(x== 2, "Coniferous woodland",
  ifelse(x==3, "Arable",
  ifelse(x==4, "Improve grassland",
  ifelse(x==5, "Neutral grassland",
  ifelse(x==6, "Calcareous grassland",
  ifelse(x==7, "Acid grassland",
  ifelse(x==8, "Fen",
  ifelse(x==9, "Heather",
  ifelse(x==10, "Heather grassland",
  ifelse(x==11, "Bog",
  ifelse(x==12, "Inland rock",
  ifelse(x==13, "Saltwater",
  ifelse(x==14, "Freshwater",
  ifelse(x==15, "Superlittoral rock",
  ifelse(x==16, "Superlittoral sediment",
  ifelse(x==17, "Lutteroal rock",
  ifelse(x==18, "Littoral sediment",
  ifelse(x==19, "Saltmarsh",
  ifelse(x==20, "Urban",
  ifelse(x==21,"Suburban","Uknown")
  ))))))))))))))))))))
}


#Import data ---------------------------------------------------------------------------------
meta <- read.csv("Deer metadata and schedules2.csv")
  meta$Collar.ID <- as.factor(meta$Collar.ID)
deer <- read.csv("deer_data_all_20230905v2.csv")
land <- terra::rast("UKCEH_Landcover2021/data/LCM.tif")
road <- read_sf("Hotos_grb_scotland_road/hotosm_gbr_scotland_roads_lines.shp")
tempjan <-terra::rast("Copernicus_data/temperature_january.nc") #
  names(tempjan) <- time(tempjan)
tempfeb <- terra::rast("Copernicus_data/temperature_febuary.nc") #
  names(tempfeb) <- time(tempfeb)
tempmar <- terra::rast("Copernicus_data/temperature_march.nc") #
  names(tempmar) <- time(tempmar)
tempapr <- terra::rast("Copernicus_data/temperature_april.nc") #
  names(tempapr) <- time(tempapr)
tempmay <- terra::rast("Copernicus_data/temperature_may.nc") #
  names(tempmay) <- time(tempmay)
tempjun <- terra::rast("Copernicus_data/temperature_june.nc") #
  names(tempjun) <- time(tempjun)
tempjul <- terra::rast("Copernicus_data/temperature_july.nc") #
  names(tempjul) <- time(tempjul)
tempaug <- terra::rast("Copernicus_data/temperature_august.nc") #
  names(tempaug) <- names(tempaug)
tempsep <- terra::rast("Copernicus_data/temperature_september.nc") #
  names(tempsep) <- time(tempsep)
  temp <- terra::merge(sprc(c(tempjan, tempfeb, tempmar, tempapr, tempmay, tempjun, tempjul, tempaug, tempsep)))# combine datespe
  time(temp) <- as.POSIXct(time(temp), origin= "1970-01-01 00:00:00 UTC")
  names(temp) <- as.POSIXct(time(temp), origin= "1970-01-01 00:00:00 UTC")
perc <- terra::rast("Copernicus_data/percepitation.nc")
  time(perc) <- as.POSIXct(time(perc), origin= "1970-01-01 00:00:00 UTC")
  names(perc) <- as.POSIXct(time(perc), origin= "1970-01-01 00:00:00 UTC")
elev1 <- terra::rast("USGS_EROS_digital_elevation/n57_w003_1arc_v3.tif")
elev2 <- terra::rast("USGS_EROS_digital_elevation/n56_w003_1arc_v3.tif")
  elev <- terra::merge(sprc(elev1, elev2)) #combining elevation
rm(elev1) #delete
rm(elev2) #delete
rm(tempjan, tempfeb, tempmar, tempapr, tempmay, tempjun, tempjul, tempaug, tempsep) #delete



#Data extraction -------------------------------------------------------------------------------
st_bbox(st_transform(deerbuf, st_crs(4326)))# limites
    #MODITools
  mt_products() #List of poducts
  mt_bands("MOD13Q1") #bands of data for NDVI 
  
  NDVI_AREA <- mt_subset(product = "MOD13Q1", #NDVI extract
                         lat =  57.09111,
                         lon = -2.623463,
                         band = "250m_16_days_NDVI",
                         start = "2023-01-21",
                         end = "2023-09-05",
                         km_lr =25, #+- from center of lat for download
                         km_ab = 25, #+- from center of long for download
                         internal = TRUE,
                         progress = TRUE)
    #turn into raster
  NDVI_AREA$calendar_date <- as.POSIXct(NDVI_AREA$calendar_date, format="%Y-%m-%d", tz="GMT")
  dates<-as.POSIXct(unique(NDVI_AREA$calendar_date), format="%Y-%m-%d", tz="GMT") #list dates
  
  NDVI.list<-list()
  NDVI_r.list<-list()
  
  for (i in 1:length(dates)){ #for each date
    NDVI.list[[i]]<-subset(NDVI_AREA, calendar_date==dates[i]) #subset by date
    NDVI_r.list[[i]] <- MODISTools::mt_to_terra(df = NDVI.list[[i]], reproject = TRUE) #convert to spatraster #convert to spatraster
    names(NDVI_r.list[[i]]) <- as.POSIXlt(dates[i], format="%Y-%m-%d", tz="GMT", origin= "1970-01-01 00:00:00 UTC") #rename
    crs(NDVI_r.list[[i]]) <- crs("EPSG:4326") #set crs
    NDVI_r.list[[i]] <- terra::project(NDVI_r.list[[i]], "EPSG:27700") #reproject to UK
    terra::time(NDVI_r.list[[i]]) <- dates[i] #set date
    NDVI.r <- terra::rast(NDVI_r.list) #combine
    svMisc::progress( i , length(dates), progress.bar = TRUE) #show progress rate  
  }#From subset take value and turn into raster
  
  
  rm(NDVI_r.list) #delete
  rm(NDVI.list) #delete
 
  
  
#Coordinate ref -------------------------------------------------------------
WGS84 <- st_crs(4326) #epsg for WGS84
OSGB36 <- st_crs(27700) #epsg for OSGB36/British National Grid

#Project to UK (OSGB36) ----------------------------------------------------------------
deersf <- st_as_sf(deer, coords = c("x","y"), na.fail = T)
st_crs(deersf) <- WGS84
deersf <- st_transform(deersf, crs = OSGB36) #re-project for UK
road <- st_transform(road, crs = OSGB36) #re-project roads to UK
elev <- terra::project(elev, "EPSG:27700") #re-preojct raster eleve to UK
NDVI.r <- terra::project(NDVI.r, "EPSG:27700") #re-preojct stacked raster NDVI to UK
temp <- terra::project(temp, "EPSG:27700") #re-preojct stacked raster temperature to UK
perc <- terra::project(perc, "EPSG:27700") #re-preojct stacked raster percipitation to UK


#Topography factors ------------------------------------------------------------
slope <- terrain(elev, "slope")
aspect <- terrain(elev, "aspect") #Direction of face (NESW),In radians
  Northness<-cos(aspect*3.14/180) #1=North -1=South
  Eastness<-sin(aspect*3.14/180) #1=East -1=West
names(elev) <- "elevation" #rename
names(Northness) <- "northness"#rename
names(Eastness) <- "eastness" #rename
topo <-rast(list(elev, slope, Northness, Eastness)) #Stack topography


#Deer data csv -> shape ------------------------------------------------------
#convert time to Posix
deer$Date <- as.POSIXct(deer$Date, format = "%Y-%m-%d %H:%M:%S", origin= "1970-01-01 00:00:00 UTC")

# convert coordinates to numeric values
deer$y  <- as.numeric(deer$Latitude)
deer$x <- as.numeric(deer$Longitude)

# check that there are no spatial outliers in x-y directions
plot(deer$x,deer$y)
plot(deer$Date,deer$x)
plot(deer$Date,deer$y)

#looking for duplicates
sum(duplicated(paste(deer$Date, deer$ID))) #number of ID and datetime duplicated
sum(duplicated(paste(deer$Date, deer$ID, deer$x, deer$y))) #number of ID, datetime and location duplicated

#Remove duplicates from Tom's advice
deer <- deer[!duplicated(paste(deer$Name,deer$Date)),]

deer <- deer[!is.na(deer$ID),] #remove points with no ID

#remove points with less then 2 satellites
deer <- deer[!deer$nSats<=2,]


#Make trajectory -------------------------------------------------------------------------
# convert points into a trajectory using 'adehabitat' package 
deerdf <- as.data.frame(deer) #must be dataframe
deertraj <- as.ltraj(cbind(as.numeric(deerdf$x),
                           as.numeric(deerdf$y)),
                     date = deerdf$Date,
                     id = deerdf$Name, 
                     typeII = TRUE) #makeing trajectory (Turnde into type II since only fixes are "real" data)

qplot(deertraj[[1]]$rel.angle)

  #trajectory into spacial lines
deerline <- ltraj2sldf(deertraj)
deerline <- st_as_sf(deerline)
st_crs(deerline) <- WGS84
deerline <- st_transform(deerline, crs = OSGB36) #reproject for UK

#Buffer ----------------------------------------------------------------------------------------
deerbuf <- deerline %>% 
              sf::st_transform(5683) %>%  # transform to a metric CRS
              sf::st_buffer(1000) # buffer by 1km
  deerbuf <- sf::st_transform(deerbuf, 27700) #re tranform to UK crs
deerbufrast <- terra::rast(deerbuf) #make raster
  st_bbox(st_transform(deerbuf, st_crs(4326))) #Coordinate boundry

##subset using buffer -------------------------------------------------------------------------

  #land
land <- land %>% 
  crop(extent(deerbuf)) %>% #crop
  mask(deerbuf) #Mask
  #roads
road <- st_intersection(road, deerbuf) #crop
  #elevation
elev <- elev %>% 
  crop(extent(deerbuf)) %>%  #Crop
  mask(deerbuf) #Mask
  #NDVI
NDVI.r <- NDVI.r %>% 
  crop(extent(deerbuf)) %>%  #Crop
  mask(deerbuf) #Mask

##Functions-----------------------------------------------------------------------------------------
#Distance matrix -------------------------------------------------------------------------
#Urban landcover matrix
uland <- ifel(land$LCM_1>19, land$LCM_1, NA) #subset only to only urban areas
uland.d <- distance(uland)
    #subset
uland.d <- uland.d %>% 
  crop(extent(deerbuf)) %>%  #Crop
  mask(deerbuf) #Mask
names(uland.d) <- "urban_dist" #rename for ease

#human landcover matrix
hland <- ifel(land$LCM_1 >=20, land$LCM_1, NA) #subset only to only urban areas
hland <- raster::merge(hland$lyr1, hland$lyr2, hland$lyr3) #combine layers
hland.d <- distance(hland)
  #subset
hland.d <- hland.d %>% 
  crop(hland.d) %>%  #Crop
  mask(deerbuf) #Mask
names(hland.d) <- "human_dist" #rename for ease

#road matrix
  #road shp -> raster
road.r<-st_rasterize(road) #shp -> SpatRast
plot(road.r)
road.r<-rast(road.r) # SpatRast -> rast
  #Distance matrix
road.rd<-distance(road.r)
    #subset
road.rd <- road.rd %>% 
  crop(extent(deerbuf)) %>%  #Crop
  mask(deerbuf) #Mask
names(road.rd) <- "road_dist" #rename for ease



#visualize all -----------------------------------------------------------------------------
mapview(deerline, zcol="id",lwd=2) +
mapview(road, zcol="highway") +
mapview(st_as_stars(land)) +
mapview(st_as_stars(elev))

#Aggregate trajectory ---------------------------------------------------------------------
  #Aggregate daily step lengths by individual and month: bc step lengths are not normally distributed, let's just take the med values
deertrajdf <- as.data.frame(ld(deertraj)) #convert daily ltraj object into spatial lines object
deertrajdf$id <- as.factor(deertrajdf$id)
deertrajdf <- na.omit(deertrajdf)
deertrajdf$month <- month(deertrajdf$date) #extract month number
deertrajdf$week <- week(deertrajdf$date) #extract week number
deertrajdf$tday <- as.numeric(format(deertrajdf$date, "%H")) + #extracting time of day
as.numeric(format(deertrajdf$date, "%M"))/60
deertrajdf$day <- day(deertrajdf$date) #extract day
deertrajdf$x2 <- deertrajdf$x+deertrajdf$dx #next point 
deertrajdf$y2 <- deertrajdf$y+deertrajdf$dy #next point


deertrajdf$distm1 <- distGeo(deertrajdf[,c("x", "y")],   #first point in step
                             deertrajdf[,c("x2", "y2")], #second point in step
                             a=6377563.396,              #major (equatorial) radius of ellipsoid value, set for OSGB36
                             f=1/299.3249646)           #ellipsoid flattening value, set for OSGB36
deertrajdf$speed <- (deertrajdf$distm/1000)/(deertrajdf$dt/(60*60)) #speed in km/h


#Outlier removal ----------------------------------------------------------------------------
deertrajdfall <- deertrajdf #save dataframe with outliers

deer.list.outlier <- list()
for(i in 1:length(ID)){ 
  deer.list.outlier[[i]] <- subset(deertrajdf, id == ID[i]) #subset by individual
  deer.list.outlier[[i]]$outlierpos <- deer.list.outlier[[i]]$speed > 
                                    stats::quantile(deer.list.outlier[[i]]$speed, prob=0.95, na.rm=TRUE) #note data with values over 1km/h
  deer.list.outlier[[i]] <- deer.list.outlier[[i]][deer.list.outlier[[i]]$outlierpos==FALSE,]
  deer.list.outlier[[i]]$id <- ID[[i]]
  deertrajdfc <- rbindlist(deer.list.outlier)
  svMisc::progress( i , length(dates), progress.bar = TRUE) #show progress rate  
}

(nrow(deertrajdfc)/nrow(deertrajdfall))*100 # % of points removed 


ggplot(deertrajdfc, aes(x=id, y=speed, fill=id)) + 
  geom_boxplot()
    
  deertrajc <- as.ltraj(cbind(as.numeric(deertrajdfc$x),
                               as.numeric(deertrajdfc$y)),
                         date = as.POSIXct(deertrajdfc$date, format = "%Y-%m-%d %H:%M:%S", origin= "1970-01-01 00:00:00 UTC"),
                         id = deertrajdfc$id, 
                         typeII = TRUE) #makeing trajectory (Turnde into type II since only fixes are "real" data)
    
    #trajectory into spacial lines
    deerlinec <- ltraj2sldf(deertrajc)
    deerlinec <- st_as_sf(deerlinec)
    st_crs(deerlinec) <- WGS84
    deerlinec <- st_transform(deerlinec, crs = OSGB36) #reproject for UK
    deerlinec <- left_join(deerlinec, meta, join_by(id == Collar.ID))
    mapview(deerlinec, zcol="id",lwd=2) + mapview(deerline, zcol="id",lwd=2)
  


##deer square net displacement visualization over time
ggplot(deertrajdf, aes(date,R2n)) + 
  geom_point() + 
  facet_wrap(vars(id)) + 
  theme_bw()

#AMT format trajectory ---------------------------------------------------------------------------

#list of individuals

  deer.track <- amt::mk_track(tbl = as.data.frame(na.omit(deertrajdfc)),
                                .x = x,.y = y,.t =date, all_cols=TRUE,crs=4326)#make track
  deer.track <- amt::transform_coords(deer.track, crs_from=4326, crs_to=27700) #reproject
  deer.track <- amt::track_resample(deer.track,
                                              rate = hours(3),tolerance = minutes(15)) #sub-sample for 3 hour rate w/ 15min variation
  deer.track <- amt::filter_min_n_burst(deer.track,min_n = 3) # only burst with more then 3 points
  deer.track <- amt::steps_by_burst(deer.track, all_cols = TRUE, keep_cols = NULL)  


deer.list <- list()
deer.track.list <- list()
ID <- unique(deer$Name)

#list of individuals
for (i in 1:length(ID)){ #for each individual
  deer.list[[i]] <- subset(deertrajdfc, id==ID[i]) #subset by individual
  deer.track.list[[i]] <- amt::mk_track(tbl = as.data.frame(deer.list[[i]]),
                              .x = x,.y = y,.t =date, all_cols=TRUE,crs=4326)#make track
  deer.track.list[[i]] <- amt::transform_coords(deer.track.list[[i]], crs_from=4326, crs_to=27700) #reproject
  deer.track.list[[i]] <- amt::track_resample(deer.track.list[[i]],
                              rate = hours(3),tolerance = minutes(15)) #sub-sample for 3 hour rate w/ 15min variation
  deer.track.list[[i]] <- amt::filter_min_n_burst(deer.track.list[[i]],min_n = 3) # only burst with more then 3 points
  deer.track.list[[i]] <- amt::steps_by_burst(deer.track.list[[i]], all_cols = TRUE) 
  svMisc::progress( i , max.value = length(ID)) #show progress rate  
}


#ISSF ------------------------------------------------------------------------------------
#extract covariates
  #list method
  issf.list <- list()
  issf.list.df <- list()
  
 for(i in 1:length(ID)){
     issf.list[[i]] <- deer.track.list[[i]] %>%
       random_steps(n_control = 10) %>% #
       extract_covariates(land, where = "both") %>% #"Both" argument to sample start and end of step
       extract_covariates(topo, where = "both") %>% 
       extract_covariates(road.rd, where = "both") %>% 
       extract_covariates(uland.d, where = "both") %>% #max_time refers to maximum time difference between the point and the raster. I've set it to 10 years because kind of makes the function redundent
       extract_covariates_var_time(temp, where="both",
                                   max_time=hours(4), name_covar = "temperature") %>%
       extract_covariates_var_time(perc, where="both",
                                   max_time=hours(4), name_covar = "precipitation") %>%
       extract_covariates_var_time(NDVI.r, where="both",
                                   max_time=days(11), name_covar = "NDVI") %>% #8 days due to 16day frequancy
       mutate(land_start =  LCM_1_start,  #scale and center covariates
              land_end = LCM_1_end,                             #\/ 
              elev_starts = scale(elevation_start),               #
              elev_ends = scale(elevation_end),                   #
              slope_starts = scale(slope_start),                  #
              slope_ends = scale(slope_end),                      #
              northness_starts = scale(northness_start),          #
              northness_ends = scale(northness_end),              #
              eastness_starts = scale(eastness_start),            #
              eastness_ends = scale(eastness_end),                #
              temperature_starts = scale(temperature_start),      #
              temperature_ends = scale(temperature_end),          #
              distroad_starts = scale(road_dist_start),           #
              distroad_ends = scale(road_dist_end),               #
              disturb_starts = scale(urban_dist_start),           #
              disturb_ends = scale(urban_dist_end),               #
              precipitation_starts = scale(precipitation_start),  #
              precipitation_ends = scale(precipitation_end),      #
              NDVI_starts = scale(NDVI_start),                    #
              NDVI_ends = scale(NDVI_end),                        #
              cos_ta_ = cos(ta_),                                 #
              log_sl_ = log(sl_))                              #/\
     issf.list[[i]]$land_start <- as.factor(issf.list[[i]]$land_start) #factor for land type
     issf.list[[i]]$land_end <- as.factor(issf.list[[i]]$land_end)
     issf.list.df[[i]] <- as.data.frame(issf.list[[i]]) #turn into df
     issf.list.df[[i]]$id <- as.factor(ID[i]) #add id column to df
     issf.df <- rbindlist(issf.list.df, fill=TRUE) #combine dataframes
     issf.df$step_id_ <- as.factor(issf.df$step_id_)
     issf.df$stratum <- as.factor(paste(issf.df$step_id_, issf.df$id)) #make stratum with id and step id
      svMisc::progress( i, length(ID)) #progress.bar = TRUE
 }
issf.df <- left_join(issf.df, meta, join_by(id == Collar.ID))
issf.df$deploy_location <- as.factor(issf.df$deploy_location)
issf.df$t1p <- as.numeric(hour(issf.df$t1_))
issf.df$t2p <- as.numeric(hour(issf.df$t2_))

issf.dfs <-issf.df #make backup 
issf.lists <- issf.list #make backup
  
saveRDS(issf.list,file="Extracted_covariates_list_n=10.Rda") #save to save loading times
saveRDS(issf.df,file="Extracted_covariates_n=10.Rda") #save to save loading times
issf.list <- readRDS("Extracted_covariates_list_n=10.Rda") #import covariate list
issf.df <- readRDS("Extracted_covariates_n=10.Rda") #import covariates




ggpairs(subset(issf.df, select=c(elevation_end, NDVI_end, #habit co-variates
                               temperature_end, precipitation_end)))

cor <- cor(na.omit(subset(issf.df, select=c(elevation_end, NDVI_end, #habit co-variates
                                 temperature_end, precipitation_end, road_dist_end, urban_dist_end)))) #enviromental factors
max(cor[cor!=max(cor)]) #highest correlation value                                


##Site and ID diagnostics
    #road distance
sum.stats.id <- as.data.frame(issf.df[issf.df$case_==T,] %>% 
                group_by(id) %>% 
                  summarize(mean(road_dist_start), mean(urban_dist_start)))
sum.stats.id <- sum.stats.id %>% 
                  rename(
                    road_dist.id = "mean(road_dist_start)",
                    urban_dist.id =  "mean(urban_dist_start)"
                  )


sum.stats.site <- as.data.frame(issf.df[issf.df$case_==T,] %>% 
                                group_by(deploy_location) %>% 
                                summarize(mean(road_dist_start), mean(urban_dist_start)))
sum.stats.site <- sum.stats.site %>% 
                  rename(
                    road_dist.site = "mean(road_dist_start)",
                    urban_dist.site =  "mean(urban_dist_start)"
                  )
                


#INTEGRATED SSF 
  #everything model
m1i <- issf.df %>% 
  fit_clogit(case_ ~ elevation_end + land_end + road_dist_end + urban_dist_end + road_dist_end:urban_dist_end + NDVI_end + #habit co-variates
             temperature_end + precipitation_end  +#enviromental factors
            log_sl_:cos_ta_ + sl_ + cos_ta_ + #movement metrics
             strata(stratum),#stratum to group ie. burst 
           model = TRUE)#to store the mode

#everything model
m1ii.omit <- issf.df[!issf.df$id == "T3HS-7359",] %>% 
  fit_clogit(case_ ~ elevation_end + land_end + road_dist_end + urban_dist_end + road_dist_end:urban_dist_end  + NDVI_end + #habit co-variates
               temperature_end + precipitation_end  +#enviromental factors
               log_sl_:cos_ta_ + sl_ + cos_ta_ + #movement metrics
               strata(stratum),#stratum to group ie. burst 
             model = TRUE)#to store the mode

#everything model
m1i.omit <- issf.df[!issf.df$id == "T3HS-7359",] %>% 
  fit_clogit(case_ ~ elevation_end + land_end + urban_dist_end + road_dist_end + road_dist_end:urban_dist_end + NDVI_end  + #habit co-variates
               temperature_end + precipitation_end  +#enviromental factors
               log_sl_:cos_ta_ + sl_ + cos_ta_ + #movement metrics
               strata(stratum),#stratum to group ie. burst 
             model = TRUE)#to store the mode
  
  #Integrate SSF

  summary(m1i) #model with strata
  summary(m1i.omit) #model with strata
  summary(m1ii.omit)
 
  #mode per individual
  m1i.list <- list()
  for(i in 1:length(ID)){
  m1i.list[[i]] <-  fit_clogit(issf.list[[i]], case_ ~ elevation_end + land_end + slope_end + NDVI_end +
                                 road_dist_end + urban_dist_end +road_dist_end:urban_dist_end +  #habit co-variates
                 temperature_end + precipitation_end  + #enviromental factors
                 log_sl_:cos_ta_ + sl_ + cos_ta_ + #movement metrics
                 strata(step_id_), #stratum to group ie. burst
                model = TRUE)#to store the mode
  names(m1i.list) <- ID[i] #give ID values to individuals in list
  svMisc::progress( i, max.value = (length(ID)))
  }
  summary(m1i.list[[3]])


# (RSS) Relative Selection Strength -------------------------------------------------------------
  ###Individual s1, s2 -------------------
   # Make a new data.frame for s1 #distanc form road = 1
    s1road <- data.frame(
      elevation_end = 0,
      slope_end = 0,
      land_end = factor(levels(issf.df$land_end)[1], levels = levels(issf.df$land_end)),
      road_dist_end = seq(from = -1, to = 4, length.out = 1000), #Create number sequencedisturb_end = 0,
      urban_dist_end = 0,
      temperature_end = 0,
      precipitation_end = 0,
      humidity_end = 0,
      NDVI_end = 0,
      sl_ = 100,
      log_sl_ = log(100),
      cos_ta_ = 1,
      id = as.vector(levels(issf.df$id)[1]),
      step_id_ = as.vector(levels(issf.df$step_id_)[1]),
      stratum = as.vector(levels(issf.df$stratum)[1]))
    # data.frame for s2 #distanc form road = 0
    s2road <- data.frame(
      elevation_end = 0,
      slope_end = 0, 
      land_end = factor(levels(issf.df$land_end)[1], levels = levels(issf.df$land_end)),
      distroad_end = 0,
      urban_dist_end = 0,
      road_dist_end = 1.5,
      temperature_end = 0,
      precipitation_end = 0,
      humidity_end = 0,
      NDVI_end = 0,
      sl_ = 100,
      log_sl_ = log(100),
      cos_ta_ = 1,
      id = as.vector(levels(issf.df$id)[1]),
      step_id_ = as.vector(levels(issf.df$step_id_)[1]),
      stratum = as.vector(levels(issf.df$stratum)[1]))

    
    # Make a new data.frame for s1 #distanc form road = 1
    s1urban <- data.frame(
      elevation_end = 0,
      slope_end = 0,
      land_end = factor(levels(issf.df$land_end)[1], levels = levels(issf.df$land_end)),
      distroad_end = 0,
      urban_dist_end = seq(from = 0, to = 1000, length.out = 1000), #Create number sequence, Scaled to be 0
      road_dist_end = 0,
      temperature_end = 0,
      precipitation_end = 0,
      humidity_end =0,
      NDVI_end = 0,
      sl_ = 100,
      log_sl_ = log(100),
      cos_ta_ = 1,
      id = as.vector(levels(issf.df$id)[1]),
      step_id_ = as.vector(levels(issf.df$step_id_)[1]),
      stratum = as.vector(levels(issf.df$stratum)[1]))
    # data.frame for s2 #distanc form road = 0
    s2urban <- data.frame(
      elevation_end = 0,
      slope_end = 0, 
      land_end = factor(levels(issf.df$land_end)[1], levels = levels(issf.df$land_end)),
      distroad_end = 0,
      urban_dist_end = 0,
      road_dist_end = 0,
      temperature_end = 0,
      precipitation_end = 0,
      humidity_end = 0,
      NDVI_end = 0,
      sl_ = 100,
      log_sl_ = log(100),
      cos_ta_ = 1,
      id = as.vector(levels(issf.df$id)[1]),
      step_id_ = as.vector(levels(issf.df$step_id_)[1]),
      stratum = as.vector(levels(issf.df$stratum)[1]))
    
    # Make a new data.frame for s1 #distanc form road = 1
    s1land <- data.frame(
      elevation_end = 0,
      slope_end = 0,
      land_end =rep(unique(issf.df$land_end), each=10),
      distroad_end = 0,
      urban_dist_end = 0, #Create number sequence, Scaled to be 0
      road_dist_end = 0,
      temperature_end = 0,
      precipitation_end = 0,
      humidity_end =0,
      NDVI_end = 0,
      sl_ = 100,
      log_sl_ = log(100),
      cos_ta_ = 1,
      id = as.vector(levels(issf.df$id)[1]),
      step_id_ = as.vector(levels(issf.df$step_id_)[1]),
      stratum = as.vector(levels(issf.df$stratum)[1]))
    # data.frame for s2 #distanc form road = 0
    s2land <- data.frame(
      elevation_end = 0,
      slope_end = 0, 
      land_end = factor(levels(issf.df$land_end)[1], levels = levels(issf.df$land_end)),
      distroad_end = 0,
      urban_dist_end = 0,
      road_dist_end = 0,
      temperature_end = 0,
      precipitation_end = 0,
      humidity_end = 0,
      NDVI_end = 0,
      sl_ = 100,
      log_sl_ = log(100),
      cos_ta_ = 1,
      id = as.vector(levels(issf.df$id)[1]),
      step_id_ = as.vector(levels(issf.df$step_id_)[1]),
      stratum = as.vector(levels(issf.df$stratum)[1]))
    
    # Make a new data.frame for s1 #distanc form road = 1
    s1NDVI <- data.frame(
      elevation_end = 0,
      slope_end = 0,
      land_end = factor(levels(issf.df$land_end)[1], levels = levels(issf.df$land_end)),
      distroad_end = 0,
      urban_dist_end = 0, #Create number sequence, Scaled to be 0
      road_dist_end = 0,
      temperature_end = 0,
      precipitation_end = 0,
      humidity_end =0,
      NDVI_end = seq(from = 0, to = 0.9294498, length.out = 1000),
      sl_ = 100,
      log_sl_ = log(100),
      cos_ta_ = 1,
      id = as.vector(levels(issf.df$id)[1]),
      step_id_ = as.vector(levels(issf.df$step_id_)[1]),
      stratum = as.vector(levels(issf.df$stratum)[1]))
    # data.frame for s2 #distanc form road = 0
    s2NDVI <- data.frame(
      elevation_end = 0,
      slope_end = 0, 
      land_end = factor(levels(issf.df$land_end)[1], levels = levels(issf.df$land_end)),
      distroad_end = 0,
      urban_dist_end = 0,
      road_dist_end = 0,
      temperature_end = 0,
      precipitation_end = 0,
      humidity_end = 0,
      NDVI_end = 0,
      sl_ = 100,
      log_sl_ = log(100),
      cos_ta_ = 1,
      id = as.vector(levels(issf.df$id)[1]),
      step_id_ = as.vector(levels(issf.df$step_id_)[1]),
      stratum = as.vector(levels(issf.df$stratum)[1]))
    
    # Make a new data.frame for s1 #distanc form road = 1
    s1elev <- data.frame(
      elevation_end = seq(from = 0, to = 400, length.out = 1000),
      slope_end = 0,
      land_end = factor(levels(issf.df$land_end)[1], levels = levels(issf.df$land_end)),
      distroad_end = 0,
      urban_dist_end = 0, #Create number sequence, Scaled to be 0
      road_dist_end = 0,
      temperature_end = 0,
      precipitation_end = 0,
      humidity_end =0,
      NDVI_end = 0,
      sl_ = 100,
      log_sl_ = log(100),
      cos_ta_ = 1,
      id = as.vector(levels(issf.df$id)[1]),
      step_id_ = as.vector(levels(issf.df$step_id_)[1]),
      stratum = as.vector(levels(issf.df$stratum)[1]))
    # data.frame for s2 #distanc form road = 0
    s2elev <- data.frame(
      elevation_end = 0,
      slope_end = 0, 
      land_end = factor(levels(issf.df$land_end)[1], levels = levels(issf.df$land_end)),
      distroad_end = 0,
      urban_dist_end = 0,
      road_dist_end = 0,
      temperature_end = 0,
      precipitation_end = 0,
      humidity_end = 0,
      NDVI_end = 0,
      sl_ = 100,
      log_sl_ = log(100),
      cos_ta_ = 1,
      id = as.vector(levels(issf.df$id)[1]),
      step_id_ = as.vector(levels(issf.df$step_id_)[1]),
      stratum = as.vector(levels(issf.df$stratum)[1]))

    
    # Make a new data.frame for s1 #distanc form road = 1
    s1slope <- data.frame(
      elevation_end = 0,
      slope_end = seq(from = 0, to = 22, length.out = 1000),
      land_end = factor(levels(issf.df$land_end)[1], levels = levels(issf.df$land_end)),
      distroad_end = 0,
      urban_dist_end = 0, #Create number sequence, Scaled to be 0
      road_dist_end = 0,
      temperature_end = 0,
      precipitation_end = 0,
      humidity_end =0,
      NDVI_end = 0,
      sl_ = 100,
      log_sl_ = log(100),
      cos_ta_ = 1,
      id = as.vector(levels(issf.df$id)[1]),
      step_id_ = as.vector(levels(issf.df$step_id_)[1]),
      stratum = as.vector(levels(issf.df$stratum)[1]))
    # data.frame for s2 #distanc form road = 0
    s2slope <- data.frame(
      elevation_end = 0,
      slope_end = 0, 
      land_end = factor(levels(issf.df$land_end)[1], levels = levels(issf.df$land_end)),
      distroad_end = 0,
      urban_dist_end = 0,
      road_dist_end = 0,
      temperature_end = 0,
      precipitation_end = 0,
      humidity_end = 0,
      NDVI_end = 0,
      sl_ = 100,
      log_sl_ = log(100),
      cos_ta_ = 1,
      id = as.vector(levels(issf.df$id)[1]),
      step_id_ = as.vector(levels(issf.df$step_id_)[1]),
      stratum = as.vector(levels(issf.df$stratum)[1]))
    
    
    # Make a new data.frame for s1 #distanc form road = 1
    s1temp <- data.frame(
      elevation_end = 0,
      slope_end = 0,
      land_end = factor(levels(issf.df$land_end)[1], levels = levels(issf.df$land_end)),
      distroad_end = 0,
      urban_dist_end = 0, #Create number sequence, Scaled to be 0
      road_dist_end = 0,
      temperature_end = seq(from = 262., to = 284, length.out = 1000),
      precipitation_end = 0,
      humidity_end =0,
      NDVI_end = 0,
      sl_ = 100,
      log_sl_ = log(100),
      cos_ta_ = 1,
      id = as.vector(levels(issf.df$id)[1]),
      step_id_ = as.vector(levels(issf.df$step_id_)[1]),
      stratum = as.vector(levels(issf.df$stratum)[1]))
    # data.frame for s2 #distanc form road = 0
    s2temp <- data.frame(
      elevation_end = 0,
      slope_end = 0, 
      land_end = factor(levels(issf.df$land_end)[1], levels = levels(issf.df$land_end)),
      distroad_end = 0,
      urban_dist_end = 0,
      road_dist_end = 0,
      temperature_end = 273.15,
      precipitation_end = 0,
      humidity_end = 0,
      NDVI_end = 0,
      sl_ = 100,
      log_sl_ = log(100),
      cos_ta_ = 1,
      id = as.vector(levels(issf.df$id)[1]),
      step_id_ = as.vector(levels(issf.df$step_id_)[1]),
      stratum = as.vector(levels(issf.df$stratum)[1]))
    
    
    # Make a new data.frame for s1 #distanc form road = 1
    s1perc <- data.frame(
      elevation_end = 0,
      slope_end = 0,
      land_end = factor(levels(issf.df$land_end)[1], levels = levels(issf.df$land_end)),
      distroad_end = 0,
      urban_dist_end = 0, #Create number sequence, Scaled to be 0
      road_dist_end = 0,
      temperature_end = 0,
      precipitation_end = seq(from = 0, to = 0.003, length.out = 1000),
      humidity_end =0,
      NDVI_end = 0,
      sl_ = 100,
      log_sl_ = log(100),
      cos_ta_ = 1,
      id = as.vector(levels(issf.df$id)[1]),
      step_id_ = as.vector(levels(issf.df$step_id_)[1]),
      stratum = as.vector(levels(issf.df$stratum)[1]))
    # data.frame for s2 #distanc form road = 0
    s2perc <- data.frame(
      elevation_end = 0,
      slope_end = 0, 
      land_end = factor(levels(issf.df$land_end)[1], levels = levels(issf.df$land_end)),
      distroad_end = 0,
      urban_dist_end = 0,
      road_dist_end = 0,
      temperature_end = 0,
      precipitation_end = 0,
      humidity_end = 0,
      NDVI_end = 0,
      sl_ = 100,
      log_sl_ = log(100),
      cos_ta_ = 1,
      id = as.vector(levels(issf.df$id)[1]),
      step_id_ = as.vector(levels(issf.df$step_id_)[1]),
      stratum = as.vector(levels(issf.df$stratum)[1]))
    
      #RSS ----------------
  lrssroad <- log_rss(m1i, s1road, s2road, ci= "se") #calculate reserource selection strength
  lrssroad.omit <- log_rss(m1i.omit, s1road, s2road, ci= "se") #calculate reserource selection strength
    lrssroad.list <- list()
    lrssroad.land <- list()
    urban.dist <- seq(from =0 , to = 6000, length.out = 6)
    urban.dist.omit <- seq(from =0 , to = 400, length.out = 6)
    land.types <- levels(issf.df$land_end)
        #road per landcover type
        for(i in 1:length(land.types)){    
          lrssroad.list[[i]] <- log_rss(m1i, 
                                    data.frame(s1road[, -which(names(s1road)=="land_end")], #remove general land_end
                                               land_end = factor(levels(issf.df$land_end)[i], 
                                                                 levels = levels(issf.df$land_end))), #add specific to individual land_end
                                    data.frame(s2road[, -which(names(s1road)=="land_end")],
                                               land_end =  factor(levels(issf.df$land_end)[i], 
                                                                  levels = levels(issf.df$land_end))),
                                    ci= "se") #calculate reserource selection strength
          lrssroad.list[[i]] <- lrssroad.list[[i]]$df #safe df
          lrssroad.list[[i]]$land_type <- as.factor(land.types[i]) #ID collumn
          lrssroad.land <- rbindlist(lrssroad.list, fill=T) #turn into 1 dataframe
          lrssroad.urban <- left_join(lrssroad.land, meta, join_by(id_x1 == Collar.ID)) #add meta data
          svMisc::progress(i, max.value = (length(land.types)))
        }
    
    
  
        #road per distance proximity to urban
    lrss.list.urb <- list()
  for(i in 1:length(urban.dist)){    
    lrssroad.list[[i]] <- log_rss(m1i, 
                                  data.frame(s1road[, -which(names(s1road)=="urban_dist_end")], #remove general land_end
                                             urban_dist_end = urban.dist[i]), #add specific to individual land_end
                                  data.frame(s2road[, -which(names(s2road)=="urban_dist_end")],
                                             urban_dist_end = urban.dist[i]),
                                  ci= "se") #calculate reserource selection strength
    lrss.list.urb[[i]] <- lrssroad.list[[i]]$df#safe df
    lrss.list.urb[[i]]$urb.dist <- as.factor(urban.dist[i]) #ID collumn
    lrssroad.disturb <- rbindlist(lrss.list.urb, fill=TRUE) #turn into 1 dataframe
    svMisc::progress(i, max.value = (length(urban.dist)))
  }
    
    lrss.list.urb.omit <- list()
    for(i in 1:length(urban.dist.omit)){    
      lrss.list.urb.omit[[i]] <- log_rss(m1i.omit, 
                                    data.frame(s1road[, -which(names(s1road)=="urban_dist_end")], #remove general land_end
                                               urban_dist_end = urban.dist.omit[i]), #add specific to individual land_end
                                    data.frame(s2road[, -which(names(s2road)=="urban_dist_end")],
                                               urban_dist_end = urban.dist.omit[i]),
                                    ci= "se") #calculate reserource selection strength
      lrss.list.urb.omit[[i]] <- lrss.list.urb.omit[[i]]$df#safe df
      lrss.list.urb.omit[[i]]$urb.dist <- as.factor(urban.dist.omit[i]) #ID collumn
      lrssroad.disturb.omit <- rbindlist(lrss.list.urb.omit, fill=TRUE) #turn into 1 dataframe
      svMisc::progress(i, max.value = (length(urban.dist.omit)))
    }
    
    
    
  lrssurb <- log_rss(m1i, s1urban, s2urban, ci= "se") #calculate reserource selection strength
  lrssurb.omit <- log_rss(m1i.omit, s1urban, s2urban, ci= "se") #calculate reserource selection strength
  
  
  
  
    #RSS per indvidual
    lrss.list.road <- list()
    lrss.list <- list()
      #distance road per individual
  for(i in 1:length(ID)){    
    lrss.list[[i]] <- log_rss(m1i.list[[i]], 
                              data.frame(s1road[, -which(names(s1road)=="land_end")], #remove general land_end
                                         land_end = factor(levels(issf.list[[i]]$land_end)[1], 
                                                           levels = levels(issf.list[[i]]$land_end))), #add specific to individual land_end
                              data.frame(s2road[, -which(names(s1road)=="land_end")],
                                         land_end =  factor(levels(issf.list[[i]]$land_end)[1], 
                                                            levels = levels(issf.list[[i]]$land_end))),
                              ci= "se") #calculate reserource selection strength
    lrss.list.road[[i]] <- lrss.list[[i]]$df #safe df
    lrss.list.road[[i]]$ID <- as.factor(ID[i]) #ID collumn
    lrss.idroad <- rbindlist(lrss.list.road, fill=T) #turn into 1 dataframe
    lrss.idroad <- left_join(lrss.idroad, meta, join_by(ID == Collar.ID)) #add meta data
    svMisc::progress(i, max.value = (length(ID)))
  }
    
   #distance urban per individual
  for(i in 1:length(ID)){    
    lrss.list.urb[[i]] <- log_rss(m1i.list[[i]], 
                                 data.frame(s1urban[, -which(names(s1road)=="land_end")], #remove general land_end
                                            land_end = factor(levels(issf.list[[i]]$land_end)[1], 
                                                              levels = levels(issf.list[[i]]$land_end))), #add specific to individual land_end
                                 data.frame(s2urban[, -which(names(s1road)=="land_end")],
                                            land_end =  factor(levels(issf.list[[i]]$land_end)[1], 
                                                               levels = levels(issf.list[[i]]$land_end))),
                                 ci= "se") #calculate reserource selection strength
    lrss.list.urb[[i]] <- lrss.list.urb[[i]]$df #safe df
    lrss.list.urb[[i]]$ID <- as.factor(ID[i]) #ID collumn
    lrss.idurb <- rbindlist(lrss.list.urb, fill=T) #turn into 1 dataframe
    lrss.idurb <- left_join(lrss.idurb, meta, join_by(ID == Collar.ID)) #add meta data
    svMisc::progress(i, max.value = (length(ID)))
  }
    
    #RSS road for proximity to urban area per individual
    lrss.list.road.id <- list()
    lrss.list.road.urb <- list()
  for(j in 1:length(urban.dist)){
    for(i in 1:length(ID)){ 
      lrss.list.road.id[[i]] <- list()
      lrss.list.road.id[[i]] <- log_rss(m1i.list[[i]], 
                                    data.frame(s1road[, -which(names(s1road)%in%c("land_end", "urban_dist_end"))], #remove general land_end
                                               land_end = factor(levels(issf.list[[i]]$land_end)[1], 
                                                                 levels = levels(issf.list[[i]]$land_end)),
                                               urban_dist_end = urban.dist[j]), #add specific to individual urban distance categories
                                    data.frame(s2road[, -which(names(s1road)%in%c("land_end", "urban_dist_end"))],
                                               land_end =  factor(levels(issf.list[[i]]$land_end)[1],
                                                                  levels = levels(issf.list[[i]]$land_end)),
                                               urban_dist_end = urban.dist[j]),
                                    ci= "se") #calculate reserource selection strength
      lrss.list.road.id[[i]] <- lrss.list.road.id[[i]]$df #safe df
      lrss.list.road.id[[i]]$id <- ID[i] #ID collumn
      lrss.list.road.id[[i]]$urban.dist <- as.factor(urban.dist[j]) #ID collumn
      lrss.list.road.urb[[j]] <- rbindlist(lrss.list.road.id, fill=T) #combine list into individual df 
      svMisc::progress(i, max.value = (length(ID)), progress.bar=TRUE)
    }
    lrss.idroad.urb <- rbindlist(lrss.list.road.urb, fill=T) #turn into 1 dataframe
    lrss.idroad.urb$id <- as.factor(lrss.idroad.urb$id)
    lrss.idroad.urb <- left_join(lrss.idroad.urb, meta, join_by(id == Collar.ID)) #add meta data
    svMisc::progress(j, max.value = (length(urban.dist)))
  }
    
    #land cover RSS
    lrssland <- log_rss(m1i, s1land, s2land, ci= "se") #calculate reserource selection strength
    
    lrssland.list <- list()
    lrssland.list.df <- list()
    for(i in 1:length(land.types)){    
      lrssland.list[[i]] <- log_rss(m1i, 
                                    data.frame(s1land), 
                                    data.frame(s2land[, -which(names(s2land)=="land_end")],#remove general land_end
                                               land_end = factor(levels(issf.df$land_end)[i], 
                                                                  levels = levels(issf.df$land_end))),#add specific to individual land_end
                                    ci= "se") #calculate reserource selection strength
      lrssland.list.df[[i]] <- lrssland.list[[i]]$df #safe df
      lrssland.list.df[[i]]$land_type <- landclass(unique(issf.df$land_end)[i]) #ID collumn
      lrsslandvland <- rbindlist(lrssland.list.df, fill=T) #turn into 1 dataframe 
      lrsslandvland$land_type <- landclass(lrsslandvland$land_type) #convert to usable name
      lrsslandvland <- left_join(lrsslandvland, meta, join_by(id_x1 == Collar.ID)) #add meta data
      svMisc::progress(i, max.value = length(land.types))
    }
    
    
    #lrss NDVI
    lrssNDVI <- log_rss(m1i.omit, s1NDVI, s2NDVI, ci= "se")
    
    #lrss slope
    lrssslope <- log_rss(m1i.omit, s1slope, s2slope, ci= "se")
    
    #lrss elevation
    lrsselev <- log_rss(m1i.omit, s1elev, s2elev, ci= "se")
    
    #lrss temperature
    lrsstemp <- log_rss(m1i.omit, s1temp, s2temp, ci="se")
    
    #lrss precipitation
    lrssperc <- log_rss(m1i.omit, s1perc, s2perc, ci="se")
                
#Visualization -------------------------------------------------------------------------------
    ##rss-------------------------------------------------------------------------------------
    ggplot(lrssroad$df, aes(x = (road_dist_end_x1*attr(issf.df$distroad_end, "scaled:scale")
                             + attr(issf.df$distroad_end, "scaled:center")),
                        y = log_rss))+
      geom_ribbon(aes(ymin = lwr, ymax= upr), alpha = 0.3, fill="skyblue")+
      geom_line(size = 1.5, color="darkorange") +
      geom_hline(yintercept=0, linetype="dashed") +
      xlab("Road distance (m)") +
      ylab("log-RSS") +
      theme_classic()
    
    road.plot <- ggplot(lrssroad.omit$df, aes(x = (road_dist_end_x1*attr(issf.df$distroad_end, "scaled:scale")
                                 + attr(issf.df$distroad_end, "scaled:center")),
                            y = log_rss))+
      geom_ribbon(aes(ymin = lwr, ymax= upr), alpha = 1, fill="beige")+
      geom_line(size = 1.5, color="darkgreen") +
      geom_hline(yintercept=0, linetype="dashed") +
      xlab("Road distance (m)") +
      ylab("log-RSS") +
      theme_classic() +
      theme(text = element_text(size=20))
    
        #rss road plot per land cover
        
        ggplot(lrssroad.land, aes(x = (road_dist_end_x1*attr(issf.df$distroad_end, "scaled:scale")
                                       +attr(issf.df$distroad_end, "scaled:center")),
                                  y = log_rss)) +
          geom_ribbon(aes(ymin = lwr, ymax= upr, fill = land_type), alpha = 0.3) +
          geom_line(size=1.5, aes(group = land_type, color = land_type)) +
          geom_hline(yintercept=0, size=1, linetype="dashed") +
          xlab("Distance from urban area (m)") +
          ylab("log-RSS") +
          facet_wrap(~land_type) +
          theme_classic()
        
        #rss road plot per urban proximity
        
        ggplot(lrssroad.disturb, aes(x = (road_dist_end_x1*attr(issf.df$distroad_end, "scaled:scale")
                                       +attr(issf.df$distroad_end, "scaled:center")),
                                  y = log_rss)) +
          geom_line(size=1.5, aes(group = urb.dist, color = urb.dist)) +
          geom_hline(yintercept=0, size=1, linetype="dashed") +
          xlab("Distance from road (m)") +
          ylab("log-RSS") +
          labs(color = "Distance from urban area (m)") +
          theme_classic() + 
          scale_color_viridis(discrete = TRUE, option = "B")+
          scale_fill_viridis(discrete = TRUE)
        
        ggplot(lrssroad.disturb.omit, aes(x = (road_dist_end_x1*attr(issf.df$distroad_end, "scaled:scale")
                                          +attr(issf.df$distroad_end, "scaled:center")),
                                     y = log_rss)) +
          geom_line(size=1.5, aes(group = urb.dist, color = urb.dist)) +
          geom_hline(yintercept=0, size=1, linetype="dashed") +
          xlab("Distance from road (m)") +
          ylab("log-RSS") +
          labs(color = "Distance from urban area (m) ") +
          theme_classic() + 
          theme(legend.position="bottom", text = element_text(size=30)) +
          scale_color_viridis(discrete = TRUE, option = "B")+
          scale_fill_viridis(discrete = TRUE)
    
    #urban distance all
      ggplot(lrssurb$df, aes(x = urban_dist_end_x1,
                            y = log_rss))+
      geom_ribbon(aes(ymin = lwr, ymax= upr), alpha = 0.3, fill="skyblue") +
      geom_line(size = 1.5, color="orange") +
      geom_hline(yintercept=0, linetype="dashed") +
      xlab("Distance urban (m)") +
      ylab("log-RSS") +
      theme_classic()
    
      urban.plot <- ggplot(lrssurb.omit$df, aes(x = urban_dist_end_x1,
                             y = log_rss))+
        geom_ribbon(aes(ymin = lwr, ymax= upr), alpha = 0.3, fill="skyblue") +
        geom_line(size = 1.5, color="orange") +
        geom_hline(yintercept=0, linetype="dashed") +
        xlab("Distance urban (m)") +
        ylab("log-RSS") +
        theme_classic() +
        theme(text = element_text(size=20))
      
  
  ggplot(lrss.idroad, aes(x = (road_dist_end_x1*attr(issf.df$distroad_end, "scaled:scale")
                             +attr(issf.df$distroad_end, "scaled:center")),
                        y = exp(log_rss))) +
      geom_ribbon(aes(ymin = lwr, ymax= upr, fill = ID), alpha = 0.3) +
      geom_line(size=1.5, aes(group = ID, color = ID)) +
      geom_hline(yintercept=0, size=1, linetype="dashed") +
      xlab("Distance from road (m)") +
      ylab("log-RSS") +
      facet_wrap(~deploy_location) +
      theme_classic()
    
 ggplot(lrss.idurb, aes(x = urban_dist_end_x1,
                        y = log_rss)) +
      geom_ribbon(aes(ymin = lwr, ymax= upr, fill = ID), alpha = 0.3) +
      geom_line(size=1.5, aes(group = ID, color = ID)) +
      geom_hline(yintercept=0, size=1, linetype="dashed") +
      xlab("Distance from urban area (m)") +
      ylab("log-RSS") +
      facet_wrap(~deploy_location) +
      theme_classic()
 
 #rss NDVI plot
 NDVIplot <- ggplot(lrssNDVI$df, aes(x = NDVI_end_x1,
                         y = log_rss))+
   geom_ribbon(aes(ymin = lwr, ymax= upr), alpha = 0.3, fill="skyblue")+
   geom_line(size = 1.5, color="darkorange") +
   geom_hline(yintercept=0, linetype="dashed") +
   xlab("Normalized difference vegetation index (NDVI)") +
   ylab("log-RSS") +
   theme_classic() +
   theme(text = element_text(size=20))
 
 #rss elevtaion plot
 elevationplot <- ggplot(lrsselev$df, aes(x = elevation_end_x1,
                         y = log_rss))+
   geom_ribbon(aes(ymin = lwr, ymax= upr), alpha = 0.3, fill="navyblue")+
   geom_line(size = 1.5, color="yellow") +
   geom_hline(yintercept=0, linetype="dashed") +
   xlab("Elevation (m)") +
   ylab("log-RSS") +
   theme_classic() +
   theme(text = element_text(size=20))
 
  #rss slope plot
 slopeplot <- ggplot(lrssslope$df, aes(x = slope_end_x1,
                         y = log_rss))+
   geom_ribbon(aes(ymin = lwr, ymax= upr), alpha = 0.3, fill="gold")+
   geom_line(size = 1.5, color="red") +
   geom_hline(yintercept=0, linetype="dashed") +
   xlab("Slope (°)") +
   ylab("log-RSS") +
   theme_classic() +
   theme(text = element_text(size=20))
 
 #rss temperature plot
 temperatureplot <- ggplot(lrsstemp$df, aes(x = (temperature_end_x1-273.15),
                          y = log_rss))+
   geom_ribbon(aes(ymin = lwr, ymax= upr), alpha = 0.3, fill="purple")+
   geom_line(size = 1.5, color="brown") +
   geom_hline(yintercept=0, linetype="dashed") +
   xlab("Temperature (c°)") +
   ylab("log-RSS") +
   theme_classic() +
   theme(text = element_text(size=20))
 
 #rss percipitation plot
 percipitationplot <- ggplot(lrssperc$df, aes(x = precipitation_end_x1,
                                            y = log_rss))+
   geom_ribbon(aes(ymin = lwr, ymax= upr), alpha = 0.3, fill="lightgreen")+
   geom_line(size = 1.5, color="khaki") +
   geom_hline(yintercept=0, linetype="dashed") +
   xlab("Precipitation (M)") +
   ylab("log-RSS") +
   theme_classic()
 
 #Combined enviromental plot
 ggarrange(ggarrange(NDVIplot, temperatureplot, ncol = 2, labels = c("A", "B")), 
           elevationplot,
           labels = c(" ", "C"),
           nrow=2,
           ncol=1)
 
 ggarrange(urban.plot, road.plot,
           labels = c("A", "B"),
           ncol=2)
 
    #RSS for road by proximity to urban area per individual
 lrss.idroad.urb <- left_join(lrss.idroad.urb, sum.stats.id, join_by(id == id)) #adding summery date per id
 lrss.idroad.urb <- left_join(lrss.idroad.urb, sum.stats.site, join_by(deploy_location == deploy_location )) #adding summery data per site
                              
 ggplot(lrss.idroad.urb, aes(x = (road_dist_end_x1*attr(issf.df$distroad_end, "scaled:scale")
                                  +attr(issf.df$distroad_end, "scaled:center")),
                        y = log_rss)) +
   geom_line(size=1.5, aes(group = urban_dist_end_x1, color = urban_dist_end_x1)) +
   geom_hline(yintercept=0, size=1, linetype="dashed") +
   xlab("Distance from road area (m)") +
   ylab("log-RSS") +
   labs(color = "Distance from urban area (m)") +
   facet_wrap(deploy_location~id, nrow=2) +
   theme(legend.position="bottom", legend.box = "horizontal", legend.text = element_text(angle = -45)) +
   scale_color_viridis(discrete=FALSE, option="inferno")

 #Combine plots
 ggarrange(urban.plot, road.plot,
          labels = c("A", "B"),
          nrow=2)
 
 
 #land type all
 ggplot(lrssland$df[!lrssland$df$land_end_x1=="5" & !lrssland$df$land_end_x1=="20",], aes(x = landclass2(land_end_x1),
                        y = log_rss)) +
   geom_hline(yintercept=0, size=0.5, linetype="dashed") +
   geom_point(size=4, aes(color = landclass2(land_end_x1))) +
   geom_errorbar(size=1.5, aes(ymin = lwr, ymax= upr, color = landclass2(land_end_x1))) +
   xlab("Land cover type") +
   ylab("log-RSS") +
   theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust=1),legend.position = "none",text = element_text(size=20))
 
    #landvoer type per type
 ggplot(lrsslandvland, aes(x = land_end_x1,
                         y = log_rss)) +
   geom_point(aes(group = land_type, color = land_type)) +
   geom_hline(yintercept=0, linetype="dashed") +
   xlab("Land cover type") +
   ylab("log-RSS") +
   theme_classic()

 
 
 
 
 

  ##Histogram -----------------------------------------------------------------------------------
 all.hist <- ggplot(issf.df[!issf.df$case_ == FALSE,], aes(x=urban_dist_end, fill=deploy_location)) +
   geom_histogram(bins=100)+
   xlab("Urban distance (m)") +
   ylab("Frequancy") +
   labs(fill = "Site") +
   theme_classic() +
   theme(legend.position="bottom", text = element_text(size=20))
 
 ggplot(issf.df[!issf.df$case_ == FALSE,], aes(x=urban_dist_end, fill=deploy_location)) +
   geom_boxplot()+
   xlab("Site") +
   ylab("Frequancy") +
   labs(fill = "Land cover type") +
   theme_classic() +
   theme(legend.position="bottom")
 
 omit.hist <- ggplot(issf.df[!issf.df$case_ == FALSE & !issf.df$id == "T3HS-7359",], aes(x=urban_dist_end, fill=id)) +
   geom_histogram()+
   xlab("Urban distance (m)") +
   ylab("Frequancy") +
   theme(legend.position="none", theme(text = element_text(size=20))) +
   theme_classic()
 
 #anova  
  #equal variance test
 var.test(na.omit(issf.df[!issf.df$case_ == FALSE & !issf.df$id == "T3HS-7359",]$urban_dist_end),
          na.omit(issf.df[!issf.df$case_ == FALSE & issf.df$id == "T3HS-7359",]$urban_dist_end),
          alternative = "two.sided")
  #normality
 qqPlot(na.omit(issf.df[!issf.df$case_ == FALSE & !issf.df$id == "T3HS-7359",]$urban_dist_end))
 qqPlot(na.omit(issf.df[!issf.df$case_ == FALSE & issf.df$id == "T3HS-7359",]$urban_dist_end))
 ks.list <- list()
 for (i in 1:length(unique(issf.df$deploy_location))){
   ks.list[[i]] <- ks.test(na.omit(issf.df[!issf.df$case_ == FALSE & issf.df$deploy_location == unique(issf.df$deploy_location)[i],]$urban_dist_end), 
                           'pnorm')
 }
 
 
 
  site.anova <- aov(urban_dist_end ~ deploy_location,
     data = issf.df[issf.df$case_ == TRUE]
 )
  summary(site.anova)
  TukeyHSD(site.anova)
  plot(TukeyHSD(site.anova))
  
  
  
    ## Maps -------------------------------------------------------------------------------------------
sites <- unique(deerlinec$deploy_location)
deer.map <- list()
  id.col <- colorFactor("Paired", 
                       domain = deerlinec$id)

#leaflet

  
  all.map <- leaflet() %>%
      addTiles(group= "OSM background") %>%
      addProviderTiles("Esri.WorldImagery", 
                       group= "Satellite view") %>%
      addLayersControl(baseGroups = "Satellite view", 
                       position = "bottomright")  %>%
      addPolylines(data=deerlinec,
                   opacity=0.8,
                   weight = 2,
                   color=~id.col(deerlinec$id),
                   label=deerlinec$id) %>%
      addMiniMap(zoomLevelOffset = -7) %>%
      addScaleBar() %>%
      addLegend("topleft", pal = id.col, values =deerlinec$id,
              title = "Individual",
              opacity = 1)

 for(i in 1:length(unique(deerlinec$deploy_location))){
 deer.map[[i]] <- leaflet() %>%
   addTiles(group= "OSM background") %>%
   addProviderTiles("Esri.WorldImagery", 
                    group= "Satellite view") %>%
   addLayersControl(baseGroups = "Satellite view", 
                    position = "bottomright")  %>%
     addPolylines(data=deerlinec[deerlinec$deploy_location==sites[i],],
                opacity=0.8,
                weight = 2,
                color=~id.col(deerlinec[deerlinec$deploy_location==sites[i],]$id),
                label=deerlinec$id) %>%
   addScaleBar() %>%
   addMiniMap(zoomLevelOffset = -7) %>%
   addLegend("topleft", pal = id.col, values =deerlinec[deerlinec$deploy_location==sites[i],]$id,
             title = "Individual",
             opacity = 1)
 svMisc::progress(i, max.value = (length(sites)))
 }
  
  latticeview(deer.map[[1]],deer.map[[2]], deer.map[[3]], deer.map[[4]], deer.map[[5]], all.map,
              ncol=2)
  
  ggpubr::ggarrange(deer.map[[1]], deer.map[[2]],
            labels = c("A", "B"),
            nrow=2)
deer.map[[3]]

#ggmaps
  # get the map info
Aberdeenshire <- c(lon=-2.632163, lat=57.09297)
background <- get_stadiamap(bbox = c(left=-2.961833, bottom=56.949590, right=-2.302493, bottom=57.236358), 
                            zoom = 8, maptype = "stamen_terrain")

center=c(lon=-2.632163 , lat=57.09297)

ggmap(background) 

  geom_polygon(data=deerlinec, fill=id)







