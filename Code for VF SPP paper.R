library("readxl")
library("terra")
library("sp")
library("spatstat")
library("sf")
library("RColorBrewer")
library("tidyverse")
library("dplyr")
library("purrr")
library("scales")
library("spatstat.model")


vfdf<-read.csv("\\VF_Simulated_Data.csv") #Forthcoming: Will contain randomized marks based on age-sex category
vf.ys<-list()
for (i in 1:22){
  vf.ys[[i]]<-subset.ppp(vfdf,marks(vfdf)$Year_sem==i)
} #There are 22 semesters starting at 2013 semester 1 and 2023 semester 2

AZ_Counties<-st_read('//counties_-6648357430706609459 (1)')
AZ_Counties<-st_transform(AZ_Counties, crs=st_crs("EPSG:3857"))
NDVI.dat<-list() #NOTE: we are including 2012 NDVI for lag effect
for (i in 5){
  NDVI.dat[[i]]<-rast(paste0("//knb-lter-cap//NDVI_multiseason_CAPLTER_",2011+i,".tif"))
  NDVI.dat[[i]]<-project(NDVI.dat[[i]], crs(AZ_Counties), method= "cubicspline")
  # NDVI.dat[[i]]<-aggregate(NDVI.dat[[i]], fact=1000/res(NDVI.dat[[i]])[1], fun="mean", na.rm=TRUE)
}

NDVI.pts<-which(!is.na(values(NDVI.dat[[1]]))) # we choose the first element of the raster list simply because all of the rasters in the list have the same irregular shape
NDVI.cds<-xyFromCell(NDVI.dat[[1]],NDVI.pts)
points<-vect(NDVI.cds, type="points", crs=crs(NDVI.dat[[1]]))
cvh<-convHull(points)

#The above three steps simply allow us to get the actual convex hull of the 'tilted' rectangular NDVI raster shape (where values are not NA)

city.cent.cds<-data.frame(Longitude=c(-112.0740,-111.8315,-111.9400,-111.9217,-111.8413),
                          Latitude=c(33.4484,33.4152,33.4255,33.4949,33.3062),
                          city.names=c("Phoenix", "Mesa", "Tempe", "Scottsdale", "Chandler"))

city.cent.cds<-st_as_sf(city.cent.cds,coords=c("Longitude","Latitude"), crs="EPSG:4326")

city.cent.cds<-st_transform(city.cent.cds,crs=st_crs("EPSG:3857"))

city.cds<-st_coordinates(city.cent.cds)
#Normalized difference between 2013 and 2023 NDVI

avg.NDVI.13.23.diff<-(NDVI.dat[[12]]$`4_fall`-NDVI.dat[[1]]$`4_fall`)/(NDVI.dat[[12]]$`4_fall`+NDVI.dat[[1]]$`4_fall`)
plot(avg.NDVI.13.23.diff, range=c(-1,1) )
points(city.cds[,1],city.cds[,2], col="red2", pch=16)
text(city.cds[,1],city.cds[,2],labels=city.cent.cds$city.names,
     pos=c(2,4,1,2,3), col="red2", cex=1, font=2)
plot(AZ_Counties[1], add=TRUE,
     col="transparent", border="red2", lwd=2)


# Read in Phoenix rasters
phx_osm<-rast("H:\\Guest\\jginos\\OSM_Phx_L_image.tif")
phx_mask_resc<-mask(phx_osm, cvh) #create masked phx_osm for plotting
ext(phx_mask_resc)[1:4]<-ext(phx_mask_resc)[1:4]/1000 # rescale it for compatibility with data

phx_osm.resc<-phx_osm #unmasked phx_osm
ext(phx_osm.resc)[1:4]<-ext(phx_osm.resc)[1:4]/1000 #rescaled (for possible alternative when plotting: Not recommended for plotting, however)

phx_W_sf<-st_as_sf(cvh)
phx_W<-as.owin(phx_W_sf)

#5km bandwidth kernel density estimates
par(mar=c(0,0,0,0))
par(mfrow=c(2,2))
vf.mt.den<-list()
for (i in 1:11){
  vf.mt.den[[i]]<-density(vf.ys[[i]], sigma=5, dimyx=250) #The default is edge=TRUE as well, for edge correction
  vf.mt.den[[i]]$v<-(ifelse(vf.mt.den[[2*i-1]]$v>0.025,vf.mt.den[[i]]$v,NA))
  plot(vf.mt.den[[i]],
       col=alpha(reds,0.7),
       main="", zlim=c(0,0.57))
  plot(phx_mask_resc, add=TRUE)
  plot(vf.mt.den[[i]],
       col=alpha(reds,0.7),
       add=TRUE, zlim=c(0,0.57))
}

#Transform crs of sf point object to make it compatible with the phx_W_sf (sf object).
# (projected coordinate system: WGS 84/EPSG 3857). 
vf.sf<-st_transform(vf.sf, crs = st_crs("EPSG:3857"))

#Create point pattern from the sf object. 
vf.pp<-ppp(st_coordinates(vf.sf)[,1], y=st_coordinates(vf.sf)[,2],
           window = phx_W, marks = st_drop_geometry(vf.sf))

#Rescale pattern and window. After, we subset the point pattern by the window to get ride of rejected points.
phx_W.resc<-rescale.owin(phx_W, s=1000, unitname = "Kilometers")


### Plotting the VF case non-parametric adaptive density
## Code forthcoming (currently R files are in secure environment)

#######################################################################
# Read in raster data for modeling point pattern as a function of covariates
# Next: Convert rasters to pixel images for use with Spatstat functions
# We need the following function for converting rasters to pixel images.
# This code (as.im.SpatRaster2 function) was obtained from https://stackoverflow.com/questions/77912041/convert-raster-terra-to-im-object-spatstat
# Author: Robert Hijmans and edited by Adrian Baddeley
as.im.SpatRaster2 <- function(X) {
  X <- X[[1]]
  g <- as.list(X, geom=TRUE)
  
  isfact <- is.factor(X)
  if (isfact) {
    v <- matrix(as.data.frame(X)[, 1], nrow=g$nrows, ncol=g$ncols, byrow=TRUE)
  } else {
    v <- as.matrix(X, wide=TRUE)
  }
  vtype <- if(isfact) "factor" else typeof(v)
  if(vtype == "double") vtype <- "real"
  tv <- v[g$nrows:1, ]
  if(isfact) tv <- factor(tv, levels=levels(X))
  out <- list(
    v = tv,
    dim = c(g$nrows, g$ncols),
    xrange = c(g$xmin, g$xmax),
    yrange = c(g$ymin, g$ymax),
    xstep = g$xres[1],
    ystep = g$yres[1],
    xcol = g$xmin + (1:g$ncols) * g$xres[1] + 0.5 * g$xres,
    yrow = g$ymax - (g$nrows:1) * g$yres[1] + 0.5 * g$yres,
    type = vtype,
    units  = list(singular=g$units, plural=g$units, multiplier=1)
  )
  attr(out$units, "class") <- "unitname"
  attr(out, "class") <- "im"
  out
}

####Create hsi pixel image by cropping original HSI raster to phx_osm image
hsi.raster<-rast("\\HSI_Raster_10m_0_to_1.tif")
hsi.raster.c<-crop(hsi.raster,cvh)
hsi.ras.msk<-mask(hsi.raster.c,cvh)
values(hsi.ras.msk)<-ifelse(values(hsi.ras.msk)>1,NA,values(hsi.ras.msk))
#hsi.ras.ag<-aggregate(hsi.ras.msk, fact= 100/res(hsi.ras.msk)[1], 
#                      fun="max",na.rm=TRUE)

## We set Certain unpublished areas to NA that were erroneously recorded as 0 in the original data.
## Also, we set values above 1, which include areas that are mostly likely water, equal to NA 
hsi.test<-hsi.ras.msk

#Topright major subregion
bbox<-ext(-12473000,-12411979,3915242,4015610)
subreg<-crop(hsi.test,bbox)
values(subreg)<-ifelse(values(subreg)==0,NA,values(subreg))
#plot(subreg)
subreg.coords<-xyFromCell(subreg, 1:ncell(subreg))
new_cells<-cellFromXY(hsi.test, subreg.coords)
hsi.test[new_cells]<-values(subreg)

#bottom left subregion
bbox<-ext(-12560836,-12540000,3915242,3925500)
subreg<-crop(hsi.test,bbox)
values(subreg)<-ifelse(values(subreg)==0,NA,values(subreg))
#plot(subreg)
subreg.coords<-xyFromCell(subreg, 1:ncell(subreg))
new_cells<-cellFromXY(hsi.test, subreg.coords)
hsi.test[new_cells]<-values(subreg)
plot(hsi.test)

hsi.im<-as.im.SpatRaster2(hsi.test)

#standardized HSI covariate

hsi.stan<-hsi.im


hsi.stan<-(hsi.stan-mean(hsi.stan,na.rm=T))/sd(hsi.stan, na.rm=TRUE)

#Rescale to Kilometers
hsi.resc.stn<-rescale.im(hsi.stan, s=1000, unitname = "Kilometers")
plot(hsi.resc.stn) #HSI covariate used


#NDVI.ws without NA's replaced
NDVI.ws.im.wna<-list()
for (i in 5){
  avg.NDVI<- (NDVI.dat[[i]]$`1_winter`+NDVI.dat[[i]]$`2_spring`)/2
  NDVI.ws.mask<-mask(avg.NDVI,cvh) 
  NDVI.ws.im.wna[[i]]<-as.im.SpatRaster2(NDVI.ws.mask)
}
#### standardize NDVI.ws.im.stan
NDVI.ws.im.stan<-NDVI.ws.im.wna
for (i in 1:length(NDVI.ws.im.stan)){
  NDVI.ws.im.stan[[i]]$v<-(NDVI.ws.im.stan[[i]]$v-mean(NDVI.ws.im.stan[[i]]$v,na.rm=T))/sd(NDVI.ws.im.stan[[i]]$v, na.rm=T)
}
# Rescale to Kilometers
NDVI.ws.rsc.stn<-list()
for(i in 1:12){
  NDVI.ws.rsc.stn[[i]]<-rescale.im(NDVI.ws.im.stan[[i]], s=1000, unitname = "Kilometers")
}


###### NDVI.sf with NA
NDVI.sf.im.wna<-list()
for (i in 1:length(NDVI.dat)){
  avg.NDVI<- (NDVI.dat[[i]]$`3_summer`+NDVI.dat[[i]]$`4_fall`)/2
  NDVI.sf.mask<-mask(avg.NDVI,cvh) 
  NDVI.sf.im.wna[[i]]<-as.im.SpatRaster2(NDVI.sf.mask)
}

#### standardize NDVI.sf.im.stan
NDVI.sf.im.stan<-NDVI.sf.im.wna
for (i in 1:length(NDVI.sf.im.stan)){
  NDVI.sf.im.stan[[i]]<-(NDVI.sf.im.stan[[i]]-mean(NDVI.sf.im.stan[[i]],na.rm=T))/sd(NDVI.sf.im.stan[[i]], na.rm=T)
}

#Rescale to Kilometers
NDVI.sf.rsc.stn<-list()
for(i in 1:12){
  NDVI.sf.rsc.stn[[i]]<-rescale.im(NDVI.sf.im.stan[[i]], s=1000, unitname = "Kilometers")
}

NDVI.diff.im<-list()
for (i in 1:11){
  NDVI.diff.im[[2*i-1]]<-(NDVI.ws.im.wna[[i+1]]-NDVI.sf.im.wna[[i]])
  NDVI.diff.im[[2*i]]<-(NDVI.sf.im.wna[[i+1]]-NDVI.ws.im.wna[[i+1]])
}
# differenced NDIV standardized covariate

NDVI.diff.stan<-NDVI.diff.im
for (i in 1:length(NDVI.diff.stan)){
  NDVI.diff.stan[[i]]<-(NDVI.diff.stan[[i]]-mean(NDVI.diff.stan[[i]]$v,na.rm=T))/sd(NDVI.diff.stan[[i]], na.rm=T)
}
### Rescale to Kilometers
NDVI.d.rsc.stn<-list()
for(i in 1:22){
  NDVI.d.rsc.stn[[i]]<-rescale.im(NDVI.diff.stan[[i]], s=1000, unitname = "Kilometers")
}

#Land cover change data prep
#Projecting first and then aggregating did not produce a large difference compared to aggregating fist and then projecting the raster (in the above case for 2013 and 2014) 
# Moreover, it seems intuitively reasonable to project first so that aggregation with respect to resolution takes place in the projected coordinate system used for analysis

lnd.cov.ch.Dev<-list()
lc.Dev.categ<-c(5121:5124,5221:5224)

for (i in 1:11){
  lnd.cov.ch.Dev[[i]]<-rast(paste0("\\NLCD_g7SX3oX8v6FFj9gn0jdV\\Annual_NLCD_LndChg_",2012+i,"_CU_C1V0_g7SX3oX8v6FFj9gn0jdV.tiff"))
  vals<-as.character(values(lnd.cov.ch.Dev[[i]]))
  values(lnd.cov.ch.Dev[[i]])<-ifelse((values(lnd.cov.ch.Dev[[i]])%in%lc.Dev.categ) ,1,0)
  lnd.cov.ch.Dev[[i]]<-project(lnd.cov.ch.Dev[[i]],crs(AZ_Cities), method="near")
  lnd.cov.ch.Dev[[i]]<-aggregate(lnd.cov.ch.Dev[[i]], fact=1000/res(lnd.cov.ch.Dev[[i]])[1], fun="any", na.rm=TRUE)
}


lnd.ch.Dev.im<-list()
for (i in 1:length(lnd.cov.ch.Dev)){
  lch.crop<-crop(lnd.cov.ch.Dev[[i]],cvh)
  lch.mask<-mask(lch.crop,cvh)
  lnd.ch.Dev.im[[i]]<-as.im.SpatRaster2(lch.mask)
  lnd.ch.Dev.im[[i]]<-as.im(lnd.ch.Dev.im[[i]])
  lnd.ch.Dev.im[[i]]<-cut(lnd.ch.Dev.im[[i]], 
                          breaks=2, 
                          labels=c("No_Change", "Change"))
}

lch.resc<-list()
for(i in 1:11){
  lch.resc[[i]]<-rescale.im(lnd.ch.Dev.im[[i]], s=1000, unitname = "Kilometers")
}
# Land cover data prep
lnd.cov<-list()
for (i in 1:11){
  lnd.cov[[i]]<-rast(paste0("\\NLCD_g7SX3oX8v6FFj9gn0jdV\\Annual_NLCD_LndCov_",2012+i,"_CU_C1V0_g7SX3oX8v6FFj9gn0jdV.tiff"))
  values(lnd.cov[[i]])<-ifelse(values(lnd.cov[[i]])==21,"Developed_Open",
                               ifelse(values(lnd.cov[[i]])==22,"Developed_LI",
                                      ifelse(values(lnd.cov[[i]])==23,"Developed_MI",
                                             ifelse(values(lnd.cov[[i]])==24,"Developed_HI"
                                                    ,"Other"))))
  
  lnd.cov[[i]]<-project(lnd.cov[[i]],crs(AZ_Cities), method = "mode")
  #lnd.cov[[i]]<-aggregate(lnd.cov[[i]], fact=500/res(lnd.cov[[i]])[1], fun="modal")
  
}

lnd.cov.Dev.im<-list()
for (i in 1:length(lnd.cov)){
  lcov.crop<-crop(lnd.cov[[i]],cvh)
  lcov.mask<-mask(lcov.crop,cvh)
  lnd.cov.Dev.im[[i]]<-as.im.SpatRaster2(lcov.mask)
  lnd.cov.Dev.im[[i]]<-cut(lnd.cov.Dev.im[[i]],breaks=5, 
                           labels=c("Developed_Open",
                                    "Developed_LI",
                                    "Developed_MI",
                                    "Developed_HI",
                                    "Other"))
}
#Rescale to Kilometers
lcov.resc<-list()
for(i in 1:11){
  lcov.resc[[i]]<-rescale.im(lnd.cov.Dev.im[[i]], s=1000, unitname = "Kilometers")
}



#Relevel Categories so We can see model outputs in terms of change in intensity with change in Land Use Category
for (i in 1:length(lcov.resc)){
  f_vals<-relevel(as.factor(lcov.resc[[i]]$v), ref="Other")
  lcov.resc[[i]]<-im(f_vals,
                     xcol = lcov.resc[[i]]$xcol ,
                     yrow = lcov.resc[[i]]$yrow )
}

#####################################################################
# Read in and prepare census data for use in the medicaid population density covariate

#Read in medicaid ACTUAL data by county/month from AZHCCCS
Med.dat.actual<-read_xlsx("//AHCCCS By County 2013-2024.xlsx")
Med.dat.actual<-as.data.frame(Med.dat.actual[1:(nrow(Med.dat.actual)-1),1:(ncol(Med.dat.actual)-2)])
#Transpose the data to make it easier to work with (rows are ordered by month-year and columns are ordered by county)
Med.dat.actual<-as.data.frame(t(Med.dat.actual))
colnames(Med.dat.actual)<-Med.dat.actual[1,] #create column names

#Get rid of first column which we have turned into the column names
Med.dat.actual<-Med.dat.actual[2:nrow(Med.dat.actual),]

#make the actual dates the row names
dates.col<-seq(as.Date("2013-01-01"),as.Date("2024-12-01"),by="months")
rownames(Med.dat.actual)<-as.character(dates.col)

#Add the columns with months and years for quickly accessing subsets of the data
months<-c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
Month.col<-rep(months,12)
Year.col<-rep(seq(2013,2024,1), each=12)

Med.dat.actual$Month<-Month.col
Med.dat.actual$Year<-Year.col


#################################################################################
# Compile Medicaid by county from block group numbers for normalization by actual AZHCCCS numbers


#Set columns we want from medicaid data ACS CSV files
Medicaid.bg.estimates<-c("GEO_ID","NAME", "B27010_007E", "B27010_007M", "B27010_013E", "B27010_013M", "B27010_023E",
                         "B27010_023M", "B27010_029E", "B27010_029M", "B27010_039E", "B27010_039M", "B27010_046E",
                         "B27010_046M", "B27010_062E", "B27010_062M")
#Give useful labels for the above columns we take from medicaid ACS data by block group
Medicaid.labels<-c("Medicaid_Only_U19","MOE_Medicaid_Only_U19","Medicaid&Medicare_U19", "MOE_Medicaid&Medicare_U19",
                   "Medicaid_Only_19_to_34","MOE_Medicaid_Only_19_to_34","Medicaid&Medicare_19_to_34", "MOE_Medicaid&Medicare_19_to_34",
                   "Medicaid_Only_35_to_64","MOE_Medicaid_Only_35_to_64","Medicaid&Medicare_35_to_64", "MOE_Medicaid&Medicare_35_to_64",
                   "Medicaid&Medicare_65&Ab", "MOE_Medicaid&Medicare_65&Ab","All_Medicaid", "All_MOE")

# Create space for importing csv's with medicaid by bg info
med.dat<-list()
#Organize data into data frames and store them in the med.dat list object
for (i in 1:11){
  med.dat[[i]]<-as.data.frame(read.csv(paste0("//ACSDT5Y",2012+i,".B27010-Data.csv")))
  med.dat[[i]]<-med.dat[[i]][2:nrow(med.dat[[i]]),Medicaid.bg.estimates]
  for (k in 3:16){
    med.dat[[i]][,k]<-as.numeric(med.dat[[i]][,k])
  }
  med.dat[[i]]$All_Medicaid<-med.dat[[i]][,3]+med.dat[[i]][,5]+med.dat[[i]][,7]+med.dat[[i]][,9]+med.dat[[i]][,11]
  +med.dat[[i]][,13]+med.dat[[i]][,15]
  
  med.dat[[i]]$All_MOE<-med.dat[[i]][,4]+med.dat[[i]][,6]+med.dat[[i]][,8]+med.dat[[i]][,10]+med.dat[[i]][,12]
  +med.dat[[i]][,14]+med.dat[[i]][,16]
  for (j in 1:16){
    names(med.dat[[i]])[j+2]<-paste(2012+i,Medicaid.labels[j], sep=" ")
  }
}

#Read in shape files for the 2010 block groups and the 2020 block groups
bg.shp.2010<-st_read(dsn = "//tl_2019_04_bg")
bg.shp.2020<-st_read(dsn = "//AZ_2022_Block_Group_Shapefile") #this shouldbe the 2016 Tiger/Lines shapefile (this is currently being fixed as it is not the 2020 block groups)
AZ_Counties<-st_read(dsn='//AZ_Counties_SHP')

#To standardize the county names across data sets, it is helpful to make the changes below to the AZ_Counties shape file
Counties<-c("Apache","Cochise","Coconino","Gila","Graham","Greenlee","La Paz", 
            "Maricopa","Mohave", "Navajo", "Pima","Pinal","Santa Cruz","Yavapai","Yuma")

for (i in 1:length(AZ_Counties$NAME)){
  AZ_Counties$NAME<- gsub(AZ_Counties$NAME[i],Counties[i],AZ_Counties$NAME)
}

#Change CRS for block group shape files to that of AZ_Counties
bg.shp.2010<-st_transform(bg.shp.2010, crs = st_crs(AZ_Counties))
bg.shp.2020<-st_transform(bg.shp.2020, crs = st_crs(AZ_Counties))
#The following standardizes the name of the column with GEOID so that the subsequent lines will run and we can connect data from block group shape files and medicaid bg information
bg.shp.2020$GEOID<- bg.shp.2020$GEOID20
#bg.shp.2010$GEOID<- paste0(0,bg.shp.2010$GEOID)

#Change bg GEOID format of the medicaid bg information to match that of the bg shape files for 2010 and 2020
for (i in 1:11){
  med.dat[[i]]$GEO_ID<-sub("1500000US","",med.dat[[i]]$GEO_ID, fixed=TRUE)
  med.dat[[i]]$GEOID<- med.dat[[i]]$GEO_ID
}

#Reduce data frames in list to one single data frame and join it to the block group shape files for 2010 and 2020

med.dat.bg.2010<-med.dat[c(seq(1,7,1))] %>% reduce(left_join,by='GEOID')
med.dat.bg.2020<-med.dat[c(seq(8,11,1))] %>% reduce(left_join,by='GEOID')

bg.medicaid.2010<-merge((bg.shp.2010),med.dat.bg.2010, by="GEOID", all=TRUE)
bg.medicaid.2020<-merge((bg.shp.2020),med.dat.bg.2020, by="GEOID", all=TRUE)

sum(is.na(st_drop_geometry(bg.shp.2010)))
sum(is.na(st_drop_geometry(bg.shp.2020)))

#Count medicaid in bg's and compare to counties
#Store the results by county in a matrix with appropriate column and row names for years and counties, respectively
Medicaid.by.County<-as.data.frame(matrix(data = rep(NA,15*11),ncol=11, nrow=15))
rownames(Medicaid.by.County)<-Counties
colnames(Medicaid.by.County)<-seq(2013,2023,1)

for (i in 1:length(AZ_Counties$NAME)){
  bg.Counties<-st_intersection(AZ_Counties[AZ_Counties$NAME==AZ_Counties$NAME[i],],bg.medicaid.2010 )
  for (j in 1:7){
    Medicaid.by.County[i,j]<-  sum(bg.Counties[,paste0("X",2012+j,".All_Medicaid")][[1]]) #+ 0.5*sum(bg.Counties[,paste0("X",2012+j,".All_MOE")][[1]])
  }
}


for (i in 1:length(AZ_Counties$NAME)){
  bg.Counties<-st_intersection(AZ_Counties[AZ_Counties$NAME==AZ_Counties$NAME[i],],bg.medicaid.2020)
  for (j in 8:11){
    Medicaid.by.County[i,j]<- sum(bg.Counties[,paste0("X",2012+j,".All_Medicaid")][[1]]) 
    #+0.5*sum(bg.Counties[,paste0("X",2012+j,".All_MOE")][[1]])
  }
}

Medicaid.by.County

#Obtain the sum of number of people on medicaid by county by decomposing the 'NAME' attribute in the census bg information
############################################################### Read-in and incorporate population data by bg

for (j in 1:7){
  for (i in colnames(Med.dat.actual)[1:15]){
    Avg_Act_Medicaid<-mean(as.numeric(Med.dat.actual[Med.dat.actual$Year==(2012+j),which(colnames(Med.dat.actual)==i)]))
    
    ACS_Medicaid_sum_bg<-Medicaid.by.County[which(rownames(Medicaid.by.County)==i),j] #See Medicaid.by.County above
    
    bg_in_county_i<-st_intersection(AZ_Counties[AZ_Counties$NAME==i,],bg.medicaid.2010 )["GEOID"][[1]]
    
    bg.subset.all.med<- bg.medicaid.2010[bg.medicaid.2010$GEOID %in% bg_in_county_i,paste0(2012+j," All_Medicaid")][[1]]
    print(sum(bg.subset.all.med == 0))
    bg.subset.all.moe<-bg.medicaid.2010[bg.medicaid.2010$GEOID %in% bg_in_county_i,paste0(2012+j," All_MOE")][[1]]
    
    bg.medicaid.2010[bg.medicaid.2010$GEOID %in%bg_in_county_i,paste0(2012+j," Normalized_Medicaid")]<-
      ifelse( bg.subset.all.med==0,
              (0.25)*bg.subset.all.moe*Avg_Act_Medicaid/ACS_Medicaid_sum_bg,
              bg.subset.all.med*Avg_Act_Medicaid/ACS_Medicaid_sum_bg )
  }
}

for (j in 8:11){
  for (i in colnames(Med.dat.actual)[1:15]){
    Avg_Act_Medicaid<-mean(as.numeric(Med.dat.actual[Med.dat.actual$Year==(2012+j),which(colnames(Med.dat.actual)==i)]))
    
    ACS_Medicaid_sum_bg<-Medicaid.by.County[which(rownames(Medicaid.by.County)==i),j] #See Medicaid.by.County above
    
    bg_in_county_i<-st_intersection(AZ_Counties[AZ_Counties$NAME==i,],bg.medicaid.2020 )["GEOID"][[1]]
    
    bg.subset.all.med<- bg.medicaid.2020[bg.medicaid.2020$GEOID %in% bg_in_county_i,paste0(2012+j," All_Medicaid")][[1]]
    
    bg.subset.all.moe<-bg.medicaid.2020[bg.medicaid.2020$GEOID %in% bg_in_county_i,paste0(2012+j," All_MOE")][[1]]
    
    bg.medicaid.2020[bg.medicaid.2020$GEOID %in%bg_in_county_i,paste0(2012+j," Normalized_Medicaid")]<-
      ifelse( bg.subset.all.med==0,
              (0.25)*bg.subset.all.moe*Avg_Act_Medicaid/ACS_Medicaid_sum_bg,
              bg.subset.all.med*Avg_Act_Medicaid/ACS_Medicaid_sum_bg )
  }
}
# Compare normalized sums with actual county medicaid values from AZHCCCS: They should be equal mdf mn nmm 
Normalized.bg<-as.data.frame(matrix(data = rep(NA,15*11),ncol=11, nrow=15))
rownames(Normalized.bg)<-Counties
colnames(Normalized.bg)<-seq(2013,2023,1)

for (i in 1:length(rownames(Normalized.bg))){
  bg.Counties<-st_intersection(AZ_Counties[AZ_Counties$NAME==AZ_Counties$NAME[i],],bg.medicaid.2010 )
  for (j in 1:7){
    Normalized.bg[i,j]<-sum(bg.Counties[,paste0("X",2012+j,".Normalized_Medicaid")][[1]])
  }
}

for (i in 1:length(rownames(Normalized.bg))){
  bg.Counties<-st_intersection(AZ_Counties[AZ_Counties$NAME==AZ_Counties$NAME[i],],bg.medicaid.2020)
  for (j in 8:11){
    Normalized.bg[i,j]<-sum(bg.Counties[,paste0("X",2012+j,".Normalized_Medicaid")][[1]]) 
  }
}
Normalized.bg
#Obtain sum of squared difference between normalized values and actual values
#First make a similar matrix by year and county for actual values (using the mean of the monthly values for each county)

Med.actual.county<-as.data.frame(matrix(data = rep(NA,15*11),ncol=11, nrow=15))
rownames(Med.actual.county)<-Counties
colnames(Med.actual.county)<-seq(2013,2023,1)
for (i in colnames(Med.dat.actual)){
  for (j in 2013:2023){
    Med.actual.county[which(rownames(Med.actual.county)==i),j-2012]<-
      mean(as.numeric(Med.dat.actual[Med.dat.actual$Year==j,which(colnames(Med.dat.actual)==i)]))
  }
}
sum(abs(Normalized.bg-Med.actual.county)) #Total difference in entries is 161449.9 if we include years 2013-2023
#This number is currently 152017 since we added a portion of the margin of error of the medicaid population estimate in block groups suspected to have zero individuals on medicaid

AZ_Metro.bg.2010<-st_intersection(caplter.proj, bg.medicaid.2010) #combining cities shapefile with 2010 block groups
AZ_Metro.bg.2020<-st_intersection(caplter.proj, bg.medicaid.2020) #combining cities shapefile with 2020 block groups

#### No difference when re running code related to creation of AZ_metro.bg files
# Create medicaid population density covariate that is a (rasterized, and then converted to the) pixel-image version of the shapefile 
# with pixels at the 1000-meter resoluation
for (i in 1:7){
  AZ_Metro.bg.2010[,paste0("Med_Pop_Den",2012+i)]<-
    AZ_Metro.bg.2010[,paste0("X",2012+i,".Normalized_Medicaid")]/(AZ_Metro.bg.2010$ALAND/(2.59*(10^6)))
}
for (i in 8:11){
  AZ_Metro.bg.2020[,paste0("Med_Pop_Den",2012+i)]<-
    AZ_Metro.bg.2020[,paste0("X",2012+i,".Normalized_Medicaid")][[1]]/(AZ_Metro.bg.2020$Shape__Are/(2.59*(10^6)))   
}

med.pop.den.im<-list()

v<-vect(AZ_Metro.bg.2010)
y<-rast(v,res=1000)

for(i in 1:7){
  r<-rasterize(x=v,y=y, field=paste0("Med_Pop_Den",2012+i), fun="mean")
  med.pop.den.im[[i]]<-as.im.SpatRaster2(r) #Creates 'im' object with medicaid by bg information
}
# no issues rerunning the rasterize command either when it comes to reproducibility
v2<-vect(AZ_Metro.bg.2020)
y2<-rast(v,res=1000)

for(i in 8:11){
  r<-rasterize(x=v2,y=y2, field=paste0("Med_Pop_Den",2012+i), fun="mean")
  med.pop.den.im[[i]]<-as.im.SpatRaster2(r) #Creates 'im' object with medicaid by bg information
}

#Smoothing the medicaid population density raster
blr.med.den<-list()
for (i in 1:length(med.pop.den.im)){
  blr.med.den[[i]]<-blur(med.pop.den.im[[i]], sigma=2000, normalise = TRUE)
}

#### Create an offset that is the pixel-sum-normalized (divide each pixel by the pixel sum for the whole image) medicaid population density.
blr.mpd.resc.ofst.stan<-list()
for(i in 1:11){
  blr.mpd.resc.ofst.stan[[i]]<-rescale.im(blr.med.den[[i]], s=1000, unitname = "Kilometers")
  blr.mpd.resc.ofst.stan[[i]]<-blr.mpd.resc.ofst.stan[[i]]/sum(blr.mpd.resc.ofst.stan[[i]], na.rm=TRUE)
}

############################################################################################################
####################### Offset only models

ofst.only.s1<-list()
ofst.only.s2<-list()

#Loop through all 6 covariates and develop a quadscheme that avoids areas where these covariates are NA 
for (i in 1:11){
  seed<-2003+(2*i-1)
  set.seed(seed)
  covs.s1<-list(lcov.resc[[i]],
                lch.resc[[i]],
                NDVI.d.rsc.stn[[2*i-1]],
                NDVI.ws.rsc.stn[[i+1]],
                hsi.resc.stn,
                blr.mpd.resc.ofst.stan[[i]])
  
  
  Q.s1<-quadscheme(vf.ys[[2*i-1]],nd=250)
  X<-union.quad(Q.s1)
  cov.lkup1<-lookup.im(covs.s1[[1]], X$x, X$y)
  cov.lkup.na<-(is.na(cov.lkup1))
  for (j in 2:6){
    cov.lkup<-lookup.im(covs.s1[[j]], X$x, X$y)
    cov.lkup.na<-cov.lkup.na+(is.na(cov.lkup))
  }
  Q.s1$data<-as.ppp(X[(is.data(Q.s1))&(cov.lkup.na==0)])
  Q.s1$dummy<-X[(!is.data(Q.s1))&(cov.lkup.na==0)]
  Q.s1$w<-Q.s1$w[(cov.lkup.na==0)]
  
  ofst.only.s1[[i]]<-ppm(Q.s1~marks
                         +offset(log(sm.mpd.ofst.stan))
                         # +lnd.cov
                         # +lnd.ch
                         # +NDVI_cur.stan+I(NDVI_cur.stan^2)
                         # +hsi.stan+I(hsi.stan^2)
                         # +NDVI.diff.stan+I(NDVI.diff.stan^2)
                         ,
                         covariates=list(lnd.cov=lcov.resc[[i]],
                                         lnd.ch=lch.resc[[i]],
                                         NDVI.diff.stan=NDVI.d.rsc.stn[[2*i-1]],
                                         NDVI_cur.stan=NDVI.ws.rsc.stn[[i+1]],
                                         hsi.stan=hsi.resc.stn,
                                         sm.mpd.ofst.stan= blr.mpd.resc.ofst.stan[[i]])
                         
                         , improve.type="ho", nsim=100)
  attr(ofst.only.s1[[i]],"seed")<-seed
  
  seed<-2003+(2*i)
  set.seed(seed)
  
  covs.s2<-list(lcov.resc[[i]],
                lch.resc[[i]],
                NDVI.d.rsc.stn[[2*i]],
                NDVI.sf.rsc.stn[[i+1]],
                hsi.resc.stn,
                blr.mpd.resc.ofst.stan[[i]])
  
  
  Q.s2<-quadscheme(vf.ys[[2*i]],nd=250)
  X<-union.quad(Q.s2)
  cov.lkup1<-lookup.im(covs.s2[[1]], X$x, X$y)
  cov.lkup.na<-(is.na(cov.lkup1))
  for (j in 2:6){
    cov.lkup<-lookup.im(covs.s2[[j]], X$x, X$y)
    cov.lkup.na<-cov.lkup.na+(is.na(cov.lkup))
  }
  Q.s2$data<-as.ppp(X[(is.data(Q.s2))&(cov.lkup.na==0)])
  Q.s2$dummy<-X[(!is.data(Q.s2))&(cov.lkup.na==0)]
  Q.s2$w<-Q.s2$w[(cov.lkup.na==0)]
  
  ofst.only.s2[[i]]<-ppm(Q.s2~marks
                         +offset(log(sm.mpd.ofst.stan))
                         #+lnd.cov
                         #+lnd.ch
                         #+NDVI_cur.stan+I(NDVI_cur.stan^2)
                         #+hsi.stan+I(hsi.stan^2)
                         #+NDVI.diff.stan+I(NDVI.diff.stan^2)
                         ,
                         covariates=list(lnd.cov=lcov.resc[[i]],
                                         lnd.ch=lch.resc[[i]],
                                         NDVI.diff.stan=NDVI.d.rsc.stn[[2*i]],
                                         NDVI_cur.stan=NDVI.sf.rsc.stn[[i+1]],
                                         hsi.stan=hsi.resc.stn,
                                         sm.mpd.ofst.stan= blr.mpd.resc.ofst.stan[[i]])
                         
                         , improve.type="ho", nsim=100)
  attr(ofst.only.s2[[i]], "seed")<-seed
}  

#### View results
# Create the table of model results for the offset only models (model 1)
cm.ofst.only<-matrix(nrow= 22, 
                     ncol = 2*length(names(ofst.only.s2[[1]]$coef))+3)

years.vf<-rep(seq(2013,2023,1), each=2)
sem.rep<-rep(c(1,2),11)
mat.names<-paste0(years.vf," ",
                  sem.rep)
rownames(cm.ofst.only)<-mat.names
colnames(cm.ofst.only)<-c(rep(names(ofst.only.s2[[1]]$coef),each=2),"AIC","M_LogLik", "Converged")


for (i in 1:(nrow(cm.ofst.only)/2)){
  for (j in 1:((ncol(cm.ofst.only)-3)/2)){
    cm.ofst.only[2*i-1,2*j-1]<-coef(ofst.only.s1[[i]])[[j]]
    cm.ofst.only[2*i,2*j-1]<-coef(ofst.only.s2[[i]])[[j]]
    
    smry1<-summary(ofst.only.s1[[i]])
    smry2<-summary(ofst.only.s2[[i]])
    cm.ofst.only[2*i-1,2*j]<-as.character(smry1$coefs.SE.CI$Ztest[[j]])
    cm.ofst.only[2*i,2*j]<-as.character(smry2$coefs.SE.CI$Ztest[[j]])
  }
}
for (i in 1:(nrow(cm.ofst.only)/2)){
  cm.ofst.only[2*i-1,"AIC"]<-
    AIC.ppm(ofst.only.s1[[i]])
  cm.ofst.only[2*i,"AIC"]<-
    AIC.ppm(ofst.only.s2[[i]])
}
for (i in 1:(nrow(cm.ofst.only)/2)){
  cm.ofst.only[2*i-1,"M_LogLik"]<-
    logLik.ppm(ofst.only.s1[[i]])
  cm.ofst.only[2*i,"M_LogLik"]<-
    logLik.ppm(ofst.only.s2[[i]])
}

for (i in 1:11){
  summs1<-summary(ofst.only.s1[[i]])
  cm.ofst.only[2*i-1,"Converged"]<-summs1$converged
  
  summs2<-summary(ofst.only.s2[[i]])
  cm.ofst.only[2*i,"Converged"]<- summs2$converged
}


##################################################
##################################################
#### Offset with land cover covariates. 

lc.mods.s1<-list()
lc.mods.s2<-list()


#Loop through all 6 covariates and develop a quadscheme that avoids areas where these covariates are NA 
for (i in 1:11){
  seed<-3003+(2*i-1)
  set.seed(seed)
  
  covs.s1<-list(lcov.resc[[i]],
                lch.resc[[i]],
                NDVI.d.rsc.stn[[2*i-1]],
                NDVI.ws.rsc.stn[[i+1]],
                hsi.resc.stn,
                blr.mpd.resc.ofst.stan[[i]])
  
  
  Q.s1<-quadscheme(vf.ys[[2*i-1]],nd=250)
  X<-union.quad(Q.s1)
  cov.lkup1<-lookup.im(covs.s1[[1]], X$x, X$y)
  cov.lkup.na<-(is.na(cov.lkup1))
  for (j in 2:6){
    cov.lkup<-lookup.im(covs.s1[[j]], X$x, X$y)
    cov.lkup.na<-cov.lkup.na+(is.na(cov.lkup))
  }
  Q.s1$data<-as.ppp(X[(is.data(Q.s1))&(cov.lkup.na==0)])
  Q.s1$dummy<-X[(!is.data(Q.s1))&(cov.lkup.na==0)]
  Q.s1$w<-Q.s1$w[(cov.lkup.na==0)]
  
  lc.mods.s1[[i]]<-ppm(Q.s1~marks
                           +offset(log(sm.mpd.ofst.stan))
                           +lnd.cov
                           +lnd.ch
                           # +NDVI_cur.stan+I(NDVI_cur.stan^2)
                           # +hsi.stan+I(hsi.stan^2)
                           # +NDVI.diff.stan+I(NDVI.diff.stan^2)
                           ,
                           covariates=list(lnd.cov=lcov.resc[[i]],
                                           lnd.ch=lch.resc[[i]],
                                           NDVI.diff.stan=NDVI.d.rsc.stn[[2*i-1]],
                                           NDVI_cur.stan=NDVI.ws.rsc.stn[[i+1]],
                                           hsi.stan=hsi.resc.stn,
                                           sm.mpd.ofst.stan= blr.mpd.resc.ofst.stan[[i]])
                           
                           , improve.type="ho", nsim=100)
  
  attr(lc.mods.s1[[i]],"seed")<-seed
  
  seed<-3003+(2*i)
  set.seed(seed)
  
  covs.s2<-list(lcov.resc[[i]],
                lch.resc[[i]],
                NDVI.d.rsc.stn[[2*i]],
                NDVI.sf.rsc.stn[[i+1]],
                hsi.resc.stn,
                blr.mpd.resc.ofst.stan[[i]])
  
  
  Q.s2<-quadscheme(vf.ys[[2*i]],nd=250)
  X<-union.quad(Q.s2)
  cov.lkup1<-lookup.im(covs.s2[[1]], X$x, X$y)
  cov.lkup.na<-(is.na(cov.lkup1))
  for (j in 2:6){
    cov.lkup<-lookup.im(covs.s2[[j]], X$x, X$y)
    cov.lkup.na<-cov.lkup.na+(is.na(cov.lkup))
  }
  Q.s2$data<-as.ppp(X[(is.data(Q.s2))&(cov.lkup.na==0)])
  Q.s2$dummy<-X[(!is.data(Q.s2))&(cov.lkup.na==0)]
  Q.s2$w<-Q.s2$w[(cov.lkup.na==0)]
  
  lc.mods.s2[[i]]<-ppm(Q.s2~marks
                           +offset(log(sm.mpd.ofst.stan))
                           +lnd.cov
                           +lnd.ch
                           #+NDVI_cur.stan+I(NDVI_cur.stan^2)
                           #+hsi.stan+I(hsi.stan^2)
                           #+NDVI.diff.stan+I(NDVI.diff.stan^2)
                           ,
                           covariates=list(lnd.cov=lcov.resc[[i]],
                                           lnd.ch=lch.resc[[i]],
                                           NDVI.diff.stan=NDVI.d.rsc.stn[[2*i]],
                                           NDVI_cur.stan=NDVI.sf.rsc.stn[[i+1]],
                                           hsi.stan=hsi.resc.stn,
                                           sm.mpd.ofst.stan= blr.mpd.resc.ofst.stan[[i]])
                           
                           , improve.type="ho", nsim=100)
  attr(lc.mods.s2[[i]],"seed")<-seed
}  


#### View results
# create the table of model results for the land cover models (model 2)

cm.lc.mods<-matrix(nrow= 22, 
                       ncol = 2*length(names(lc.mods.s2[[1]]$coef))+3)

years.vf<-rep(seq(2013,2023,1), each=2)
sem.rep<-rep(c(1,2),11)
mat.names<-paste0(years.vf," ",
                  sem.rep)
rownames(cm.lc.mods)<-mat.names
colnames(cm.lc.mods)<-c(rep(names(lc.mods.s2[[1]]$coef),each=2),"AIC","M_LogLik", "Converged")


for (i in 1:(nrow(cm.lc.mods)/2)){
  for (j in 1:((ncol(cm.lc.mods)-3)/2)){
    cm.lc.mods[2*i-1,2*j-1]<-coef(lc.mods.s1[[i]])[[j]]
    cm.lc.mods[2*i,2*j-1]<-coef(lc.mods.s2[[i]])[[j]]
    
    smry1<-summary(lc.mods.s1[[i]])
    smry2<-summary(lc.mods.s2[[i]])
    cm.lc.mods[2*i-1,2*j]<-as.character(smry1$coefs.SE.CI$Ztest[[j]])
    cm.lc.mods[2*i,2*j]<-as.character(smry2$coefs.SE.CI$Ztest[[j]])
  }
}
for (i in 1:(nrow(cm.lc.mods)/2)){
  cm.lc.mods[2*i-1,"AIC"]<-
    AIC.ppm(lc.mods.s1[[i]])
  cm.lc.mods[2*i,"AIC"]<-
    AIC.ppm(lc.mods.s2[[i]])
}
for (i in 1:(nrow(cm.lc.mods)/2)){
  cm.lc.mods[2*i-1,"M_LogLik"]<-
    logLik.ppm(lc.mods.s1[[i]])
  cm.lc.mods[2*i,"M_LogLik"]<-
    logLik.ppm(lc.mods.s2[[i]])
}

for (i in 1:11){
  summs1<-summary(lc.mods.s1[[i]])
  cm.lc.mods[2*i-1,"Converged"]<-summs1$converged
  
  summs2<-summary(lc.mods.s2[[i]])
  cm.lc.mods[2*i,"Converged"]<- summs2$converged
}
View(cm.lc.mods)

################################## Fitting base models with linear terms only (we initially called this the base model, hence the name)

base.mod.no.dup.s1<-list()
base.mod.no.dup.s2<-list()

#Loop through all 6 covariates and develop a quadscheme that avoids areas where these covariates are NA 

for (i in 1:11){
  seed<-4003+(2*i-1)
  set.seed(seed)
  covs.s1<-list(lcov.resc[[i]],
                lch.resc[[i]],
                NDVI.d.rsc.stn[[2*i-1]],
                NDVI.ws.rsc.stn[[i+1]],
                hsi.resc.stn,
                blr.mpd.resc.ofst.stan[[i]])
  
  
  Q.s1<-quadscheme(vf.ys[[2*i-1]],nd=250)
  X<-union.quad(Q.s1)
  cov.lkup1<-lookup.im(covs.s1[[1]], X$x, X$y)
  cov.lkup.na<-(is.na(cov.lkup1))
  for (j in 2:6){
    cov.lkup<-lookup.im(covs.s1[[j]], X$x, X$y)
    cov.lkup.na<-cov.lkup.na+(is.na(cov.lkup))
  }
  Q.s1$data<-as.ppp(X[(is.data(Q.s1))&(cov.lkup.na==0)])
  Q.s1$dummy<-X[(!is.data(Q.s1))&(cov.lkup.na==0)]
  Q.s1$w<-Q.s1$w[(cov.lkup.na==0)]
  
  base.mod.no.dup.s1[[i]]<-ppm(Q.s1~marks
                               +offset(log(sm.mpd.ofst.stan))
                               +lnd.cov
                               +lnd.ch
                               +NDVI_cur.stan
                               +hsi.stan
                               +NDVI.diff.stan,
                               
                               covariates=list(lnd.cov=lcov.resc[[i]],
                                               lnd.ch=lch.resc[[i]],
                                               NDVI.diff.stan=NDVI.d.rsc.stn[[2*i-1]],
                                               NDVI_cur.stan=NDVI.ws.rsc.stn[[i+1]],
                                               hsi.stan=hsi.resc.stn,
                                               sm.mpd.ofst.stan= blr.mpd.resc.ofst.stan[[i]])
                               
                               ,improve.type="ho", nsim=100)
  attr(base.mod.no.dup.s1[[i]],"seed")<-seed
  
  seed<-4003+(2*i)
  set.seed(seed)
  
  covs.s2<-list(lcov.resc[[i]],
                lch.resc[[i]],
                NDVI.d.rsc.stn[[2*i]],
                NDVI.sf.rsc.stn[[i+1]],
                hsi.resc.stn,
                blr.mpd.resc.ofst.stan[[i]])
  
  
  Q.s2<-quadscheme(vf.ys[[2*i]],nd=250)
  X<-union.quad(Q.s2)
  cov.lkup1<-lookup.im(covs.s2[[1]], X$x, X$y)
  cov.lkup.na<-(is.na(cov.lkup1))
  for (j in 2:6){
    cov.lkup<-lookup.im(covs.s2[[j]], X$x, X$y)
    cov.lkup.na<-cov.lkup.na+(is.na(cov.lkup))
  }
  Q.s2$data<-as.ppp(X[(is.data(Q.s2))&(cov.lkup.na==0)])
  Q.s2$dummy<-X[(!is.data(Q.s2))&(cov.lkup.na==0)]
  Q.s2$w<-Q.s2$w[(cov.lkup.na==0)]
  
  base.mod.no.dup.s2[[i]]<-ppm(Q.s2~marks
                               +offset(log(sm.mpd.ofst.stan))
                               +lnd.cov
                               +lnd.ch
                               +NDVI_cur.stan
                               +hsi.stan
                               +NDVI.diff.stan,
                               
                               covariates=list(lnd.cov=lcov.resc[[i]],
                                               lnd.ch=lch.resc[[i]],
                                               NDVI.diff.stan=NDVI.d.rsc.stn[[2*i]],
                                               NDVI_cur.stan=NDVI.sf.rsc.stn[[i+1]],
                                               hsi.stan=hsi.resc.stn,
                                               sm.mpd.ofst.stan= blr.mpd.resc.ofst.stan[[i]])
                               
                               ,improve.type="ho", nsim=100)
  attr(base.mod.no.dup.s2[[i]],"seed")<-seed
  
}
# create the table of model results for the linear model (model 3)

cm.base.no.dup<-matrix(nrow= 22, 
                       ncol = 2*length(names(base.mod.no.dup.s1[[1]]$coef))+3)

years.vf<-rep(seq(2013,2023,1), each=2)
sem.rep<-rep(c(1,2),11)
mat.names<-paste0(years.vf," ",
                  sem.rep)
rownames(cm.base.no.dup)<-mat.names
colnames(cm.base.no.dup)<-c(rep(names(base.mod.no.dup.s1[[1]]$coef),each=2),"AIC","M_LogLik", "Converged")


for (i in 1:(nrow(cm.base.no.dup)/2)){
  for (j in 1:((ncol(cm.base.no.dup)-3)/2)){
    cm.base.no.dup[2*i-1,2*j-1]<-coef(base.mod.no.dup.s1[[i]])[[j]]
    cm.base.no.dup[2*i,2*j-1]<-coef(base.mod.no.dup.s2[[i]])[[j]]
    
    smry1<-summary(base.mod.no.dup.s1[[i]])
    smry2<-summary(base.mod.no.dup.s2[[i]])
    cm.base.no.dup[2*i-1,2*j]<-as.character(smry1$coefs.SE.CI$Ztest[[j]])
    cm.base.no.dup[2*i,2*j]<-as.character(smry2$coefs.SE.CI$Ztest[[j]])
    
    
  }
}
for (i in 1:(nrow(cm.base.no.dup)/2)){
  cm.base.no.dup[2*i-1,"AIC"]<-
    AIC.ppm(base.mod.no.dup.s1[[i]])
  cm.base.no.dup[2*i,"AIC"]<-
    AIC.ppm(base.mod.no.dup.s2[[i]])
}
for (i in 1:(nrow(cm.base.no.dup)/2)){
  cm.base.no.dup[2*i-1,"M_LogLik"]<-
    logLik.ppm(base.mod.no.dup.s1[[i]])
  cm.base.no.dup[2*i,"M_LogLik"]<-
    logLik.ppm(base.mod.no.dup.s2[[i]])
}

for (i in 1:11){
  summs1<-summary(base.mod.no.dup.s1[[i]])
  cm.base.no.dup[2*i-1,"Converged"]<-summs1$converged
  
  summs2<-summary(base.mod.no.dup.s2[[i]])
  cm.base.no.dup[2*i,"Converged"]<- summs2$converged
}
View(cm.base.no.dup)


############################## Fitting full or quadratic models that have quadratic environmental terms

#fl.env.qd.no.dup.s1<-list()
#fl.env.qd.no.dup.s2<-list()

#Loop through all 6 covariates and develop a quadscheme that avoids areas where these covariates are NA 
for (i in 1:11){
  seed<-5003+(2*i-1)
  set.seed(seed)
  covs.s1<-list(lcov.resc[[i]],
                lch.resc[[i]],
                NDVI.d.rsc.stn[[2*i-1]],
                NDVI.ws.rsc.stn[[i+1]],
                hsi.resc.stn,
                blr.mpd.resc.ofst.stan[[i]])
  
  
  Q.s1<-quadscheme(vf.ys[[2*i-1]],nd=250)
  X<-union.quad(Q.s1)
  cov.lkup1<-lookup.im(covs.s1[[1]], X$x, X$y)
  cov.lkup.na<-(is.na(cov.lkup1))
  for (j in 2:6){
    cov.lkup<-lookup.im(covs.s1[[j]], X$x, X$y)
    cov.lkup.na<-cov.lkup.na+(is.na(cov.lkup))
  }
  Q.s1$data<-as.ppp(X[(is.data(Q.s1))&(cov.lkup.na==0)])
  Q.s1$dummy<-X[(!is.data(Q.s1))&(cov.lkup.na==0)]
  Q.s1$w<-Q.s1$w[(cov.lkup.na==0)]
  
  test<-ppm(Q.s1~marks
            +offset(log(sm.mpd.ofst.stan))
            +lnd.cov
            +lnd.ch
            +NDVI_cur.stan+I(NDVI_cur.stan^2)
            +hsi.stan+I(hsi.stan^2)
            +NDVI.diff.stan+I(NDVI.diff.stan^2)
            ,
            covariates=list(lnd.cov=lcov.resc[[i]],
                            lnd.ch=lch.resc[[i]],
                            NDVI.diff.stan=NDVI.d.rsc.stn[[2*i-1]],
                            NDVI_cur.stan=NDVI.ws.rsc.stn[[i+1]],
                            hsi.stan=hsi.resc.stn,
                            sm.mpd.ofst.stan= blr.mpd.resc.ofst.stan[[i]])
            
  )
  attr(fl.env.qd.no.dup.s1[[i]],"seed")<-seed
  
  seed<-5003+(2*i)
  set.seed(seed)
  covs.s2<-list(lcov.resc[[i]],
                lch.resc[[i]],
                NDVI.d.rsc.stn[[2*i]],
                NDVI.sf.rsc.stn[[i+1]],
                hsi.resc.stn,
                blr.mpd.resc.ofst.stan[[i]])
  
  
  Q.s2<-quadscheme(vf.ys[[2*i]],nd=250)
  X<-union.quad(Q.s2)
  cov.lkup1<-lookup.im(covs.s2[[1]], X$x, X$y)
  cov.lkup.na<-(is.na(cov.lkup1))
  for (j in 2:6){
    cov.lkup<-lookup.im(covs.s2[[j]], X$x, X$y)
    cov.lkup.na<-cov.lkup.na+(is.na(cov.lkup))
  }
  Q.s2$data<-as.ppp(X[(is.data(Q.s2))&(cov.lkup.na==0)])
  Q.s2$dummy<-X[(!is.data(Q.s2))&(cov.lkup.na==0)]
  Q.s2$w<-Q.s2$w[(cov.lkup.na==0)]
  
  fl.env.qd.no.dup.s2[[i]]<-ppm(Q.s2~marks
                                +offset(log(sm.mpd.ofst.stan))
                                +lnd.cov
                                +lnd.ch
                                +NDVI_cur.stan+I(NDVI_cur.stan^2)
                                +hsi.stan+I(hsi.stan^2)
                                +NDVI.diff.stan+I(NDVI.diff.stan^2)
                                ,
                                covariates=list(lnd.cov=lcov.resc[[i]],
                                                lnd.ch=lch.resc[[i]],
                                                NDVI.diff.stan=NDVI.d.rsc.stn[[2*i]],
                                                NDVI_cur.stan=NDVI.sf.rsc.stn[[i+1]],
                                                hsi.stan=hsi.resc.stn,
                                                sm.mpd.ofst.stan= blr.mpd.resc.ofst.stan[[i]])
                                
                                , improve.type="ho", nsim=100)
  attr(fl.env.qd.no.dup.s2[[i]],"seed")<-seed
}
# create the table of model results for the quadratic model (model 4)
cm.fl.env.qd.no.dup<-matrix(nrow= 22, 
                            ncol = 2*length(names(fl.env.qd.no.dup.s1[[1]]$coef))+3)

years.vf<-rep(seq(2013,2023,1), each=2)
sem.rep<-rep(c(1,2),11)
mat.names<-paste0(years.vf," ",
                  sem.rep)
rownames(cm.fl.env.qd.no.dup)<-mat.names
colnames(cm.fl.env.qd.no.dup)<-c(rep(names(fl.env.qd.no.dup.s1[[1]]$coef),each=2),"AIC", "Converged", "M_LogLik")


for (i in 1:(nrow(cm.fl.env.qd.no.dup)/2)){
  for (j in 1:((ncol(cm.fl.env.qd.no.dup)-3)/2)){
    cm.fl.env.qd.no.dup[2*i-1,2*j-1]<-coef(fl.env.qd.no.dup.s1[[i]])[[j]]
    cm.fl.env.qd.no.dup[2*i,2*j-1]<-coef(fl.env.qd.no.dup.s2[[i]])[[j]]
    
    smry1<-summary(fl.env.qd.no.dup.s1[[i]])
    smry2<-summary(fl.env.qd.no.dup.s2[[i]])
    cm.fl.env.qd.no.dup[2*i-1,2*j]<-as.character(smry1$coefs.SE.CI$Ztest[[j]])
    cm.fl.env.qd.no.dup[2*i,2*j]<-as.character(smry2$coefs.SE.CI$Ztest[[j]])
    
  }
} 
for (i in 1:(nrow(cm.fl.env.qd.no.dup)/2)){
  cm.fl.env.qd.no.dup[2*i-1,"AIC"]<-
    AIC.ppm(fl.env.qd.no.dup.s1[[i]])
  cm.fl.env.qd.no.dup[2*i,"AIC"]<-
    AIC.ppm(fl.env.qd.no.dup.s2[[i]])
}

for (i in 1:11){
  summs1<-summary(fl.env.qd.no.dup.s1[[i]])
  cm.fl.env.qd.no.dup[2*i-1,"Converged"]<-summs1$converged
  
  summs2<-summary(fl.env.qd.no.dup.s2[[i]])
  cm.fl.env.qd.no.dup[2*i,"Converged"]<- summs2$converged
}

for (i in 1:(nrow(cm.fl.env.qd.no.dup)/2)){
  cm.fl.env.qd.no.dup[2*i-1,"M_LogLik"]<-
    logLik.ppm(fl.env.qd.no.dup.s1[[i]])
  cm.fl.env.qd.no.dup[2*i,"M_LogLik"]<-
    logLik.ppm(fl.env.qd.no.dup.s2[[i]])
}

View(cm.fl.env.qd.no.dup)


### NDVI coefficient analysis:
### Define and Plot Avg NDVI and median NDVI
NDVI.avg.sem<-vector()
for(i in 1:length(NDVI.sf.im.wna)){
  NDVI.avg.sem[2*i-1]<-mean(NDVI.ws.im.wna[[i]], na.rm=TRUE)
  NDVI.avg.sem[2*i]<-mean(NDVI.sf.im.wna[[i]], na.rm=TRUE)
}

NDVI.med.sem<-vector()
for(i in 1:length(NDVI.sf.im)){
  NDVI.med.sem[2*i-1]<-quantile.im(NDVI.ws.im[[i]])[[3]]
  NDVI.med.sem[2*i]<-quantile.im(NDVI.sf.im[[i]])[[3]]
}


### NDVI, HSI, and NDVI diff time series of coefficients
par(mar=c(5,4,4,4)+0.25)
par(mfrow=c(2,3))
plot(cm.fl.env.qd.no.dup[,"NDVI_cur.stan"], type="l",col="red",
     main="NDVI", ylab="Coefficient",
     xlab="", xaxt='n', ylim=c(-1.2,1.2))
pts.c<-data.frame(matrix(nrow=nrow(cm.fl.env.qd.no.dup), ncol=2))
coef.cols<-which(colnames(cm.fl.env.qd.no.dup)=="NDVI_cur.stan")
coef.vals<-data.frame(coefs=as.numeric(cm.fl.env.qd.no.dup[,coef.cols[1]]),
                      signif_cds=cm.fl.env.qd.no.dup[,coef.cols[2]])
for (i in 1:nrow(cm.fl.env.qd.no.dup)){
  if(coef.vals[i,2] %in% c("*","**","***")){
    pts.c[i,]<-c(coef.vals[i,1], i)
  }
}

pts.c<-na.omit(pts.c)
axis(1, at=seq(1,22,1), labels=rep(seq(2013,2023,1),each=2))
points(y=coef.vals[,1],x=seq(1,22,1), col="red", pch=16, cex=1)
points(y=pts.c[,1],x=pts.c[,2], col="darkred", pch=16, cex=1)

par(new=TRUE)
plot(NDVI.avg.sem[3:length(NDVI.avg.sem)],type="l", col="green2",
     main="",
     ylim=c(-0.25,0.25),
     axes=FALSE,
     xlab="",
     ylab="", bty='n',
     lwd=2)
axis(4)

#### Now a plot of NDVI diff over time

plot(cm.fl.env.qd.no.dup[,"NDVI.diff.stan"], type="l",col="red",
     main="NDVI Difference",
     ylab="",
     xlab="", xaxt='n', ylim=c(-1.2,1.2))
pts<-data.frame(matrix(nrow=nrow(cm.fl.env.qd.no.dup), ncol=2))
coef.cols<-which(colnames(cm.fl.env.qd.no.dup)=="NDVI.diff.stan")
coef.vals<-data.frame(coefs=as.numeric(cm.fl.env.qd.no.dup[,coef.cols[1]]),
                      signif_cds=cm.fl.env.qd.no.dup[,coef.cols[2]])
for (i in 1:nrow(cm.fl.env.qd.no.dup)){
  if(coef.vals[i,2] %in% c("*","**","***")){
    pts[i,]<-c(coef.vals[i,1], i)
  }
}

pts<-na.omit(pts)
axis(1, at=seq(1,22,1), labels=rep(seq(2013,2023,1),each=2))
points(y=coef.vals[,1],x=seq(1,22,1), col="red", pch=16, cex=1)
points(y=pts[,1],x=pts[,2], col="darkred", pch=16, cex=1)


par(new=TRUE)
plot(NDVI.avg.sem[-(1:2)],
     type="l",
     col="green2",
     main="",
     ylim=c(-0.25,0.25),
     axes=FALSE,
     xlab="",
     ylab="", bty='n',
     lwd=2)
axis(4)

#################################### Plotting HSI coefficients over time 

plot(cm.fl.env.qd.no.dup[,"hsi.stan"], type="l",col="red",
     main="HSI",
     ylab="",
     xlab="", xaxt='n',ylim=c(-1.2,1.2))
pts<-data.frame(matrix(nrow=nrow(cm.fl.env.qd.no.dup), ncol=2))
coef.cols<-which(colnames(cm.fl.env.qd.no.dup)=="hsi.stan")
coef.vals<-data.frame(coefs=as.numeric(cm.fl.env.qd.no.dup[,coef.cols[1]]),
                      signif_cds=cm.fl.env.qd.no.dup[,coef.cols[2]])
for (i in 1:nrow(cm.fl.env.qd.no.dup)){
  if(coef.vals[i,2] %in% c("*","**","***")){
    pts[i,]<-c(coef.vals[i,1], i)
  }
}

pts<-na.omit(pts)
axis(1, at=seq(1,22,1), labels=rep(seq(2013,2023,1),each=2))
points(y=coef.vals[,1],x=seq(1,22,1), col="red", pch=16, cex=1)

points(y=pts[,1],x=pts[,2], col="darkred", pch=16, cex=1)

par(new=TRUE)
plot(NDVI.avg.sem[-(1:2)],
     type="l",
     col="green2",
     main="",
     ylim=c(-0.25,0.25),
     axes=FALSE,
     xlab="",
     ylab="", bty='n',
     lwd=2)
axis(4)


#################################################### Plotting time series of quadratic terms of the model
####################################################
####################################################

plot(cm.fl.env.qd.no.dup[,"I(NDVI_cur.stan^2)"], type="l",col="red",
     main="NDVI quadratic", ylab="Coefficient",
     xlab="", xaxt='n', ylim=c(-1.2,1.2))
pts.c<-data.frame(matrix(nrow=nrow(cm.fl.env.qd.no.dup), ncol=2))
coef.cols<-which(colnames(cm.fl.env.qd.no.dup)=="I(NDVI_cur.stan^2)")
coef.vals<-data.frame(coefs=as.numeric(cm.fl.env.qd.no.dup[,coef.cols[1]]),
                      signif_cds=cm.fl.env.qd.no.dup[,coef.cols[2]])
for (i in 1:nrow(cm.fl.env.qd.no.dup)){
  if(coef.vals[i,2] %in% c("*","**","***")){
    pts.c[i,]<-c(coef.vals[i,1], i)
  }
}

pts.c<-na.omit(pts.c)
axis(1, at=seq(1,22,1), labels=rep(seq(2013,2023,1),each=2))
points(y=coef.vals[,1],x=seq(1,22,1), col="red", pch=16, cex=1)

points(y=pts.c[,1],x=pts.c[,2], col="darkred", pch=16, cex=1) 

par(new=TRUE)
plot(NDVI.avg.sem[-(1:2)],type="l", col="green2",
     main="",
     ylim=c(-0.25,0.25),
     axes=FALSE,
     xlab="",
     ylab="", bty='n',
     lwd=2)
axis(4)

#### Now a plot of NDVI diff over time

plot(cm.fl.env.qd.no.dup[,"I(NDVI.diff.stan^2)"], type="l",col="red",
     main="NDVI Difference quadratic",
     ylab="",
     xlab="", xaxt='n', ylim=c(-1.2,1.2))
pts<-data.frame(matrix(nrow=nrow(cm.fl.env.qd.no.dup), ncol=2))
coef.cols<-which(colnames(cm.fl.env.qd.no.dup)=="I(NDVI.diff.stan^2)")
coef.vals<-data.frame(coefs=as.numeric(cm.fl.env.qd.no.dup[,coef.cols[1]]),
                      signif_cds=cm.fl.env.qd.no.dup[,coef.cols[2]])
for (i in 1:nrow(cm.fl.env.qd.no.dup)){
  if(coef.vals[i,2] %in% c("*","**","***")){
    pts[i,]<-c(coef.vals[i,1], i)
  }
}

pts<-na.omit(pts)
axis(1, at=seq(1,22,1), labels=rep(seq(2013,2023,1),each=2))
points(y=coef.vals[,1],x=seq(1,22,1), col="red", pch=16, cex=1)

points(y=pts[,1],x=pts[,2], col="darkred", pch=16, cex=1)


par(new=TRUE)
plot(NDVI.avg.sem[-(1:2)],type="l", col="green2",
     main="",
     ylim=c(-0.25,0.25),
     axes=FALSE,
     xlab="",
     ylab="", bty='n',
     lwd=2)
axis(4)

######################################### HSI quadratic term coefficient time series 
plot(cm.fl.env.qd.no.dup[,"I(hsi.stan^2)"], type="l",col="red",
     main="HSI quadratic",
     ylab="",
     xlab="", xaxt='n', ylim=c(-1.2,1.2))
pts<-data.frame(matrix(nrow=nrow(cm.fl.env.qd.no.dup), ncol=2))
coef.cols<-which(colnames(cm.fl.env.qd.no.dup)=="I(hsi.stan^2)")
coef.vals<-data.frame(coefs=as.numeric(cm.fl.env.qd.no.dup[,coef.cols[1]]),
                      signif_cds=cm.fl.env.qd.no.dup[,coef.cols[2]])
for (i in 1:nrow(cm.fl.env.qd.no.dup)){
  if(coef.vals[i,2] %in% c("*","**","***")){
    pts[i,]<-c(coef.vals[i,1], i)
  }
}

pts<-na.omit(pts)
axis(1, at=seq(1,22,1), labels=rep(seq(2013,2023,1),each=2))
points(y=coef.vals[,1],x=seq(1,22,1), col="red", pch=16, cex=1)

points(y=pts[,1],x=pts[,2], col="darkred", pch=16, cex=1)

par(new=TRUE)
plot(NDVI.avg.sem[-(1:2)],
     type="l",
     col="green2",
     main="",
     ylim=c(-0.25,0.25),
     axes=FALSE,
     xlab="",
     ylab="", bty='n',
     lwd=2)
axis(4)
dev.off()
###############################################
###############################################

#### LRT comparing offset only model to 
# quadratic model 
lrt.qd.no_ofst<-data.frame("Deviance"=rep(NA,11), "Pr.Chi.grt"=rep(NA,11),"less.thn.0.001"=rep(NA,11))
for (i in 1:11){
  anva.lrt<-anova.ppm(fl.env.qd.no.dup.s1[[i]], ofst.only.s1[[i]],fine=TRUE, test="LRT")
  lrt.qd.no_ofst[2*i-1,1]<-anva.lrt$Deviance[2]
  lrt.qd.no_ofst[2*i-1,2]<-anva.lrt$`Pr(>Chi)`[2]
  lrt.qd.no_ofst[2*i-1,3]<-ifelse(anva.lrt$`Pr(>Chi)`[2]<0.001,"***",
                                  ifelse((anva.lrt$`Pr(>Chi)`[2]<0.01)&(anva.lrt$`Pr(>Chi)`[2]>0.001),"**",
                                         ifelse((anva.lrt$`Pr(>Chi)`[2]<0.05)&(anva.lrt$`Pr(>Chi)`[2]>0.01),"*","")))
  
  anva.lrt<-anova.ppm(fl.env.qd.no.dup.s2[[i]],ofst.only.s2[[i]],fine=TRUE, test="LRT")
  lrt.qd.no_ofst[2*i,1]<-anva.lrt$Deviance[2]
  lrt.qd.no_ofst[2*i,2]<-anva.lrt$`Pr(>Chi)`[2]
  lrt.qd.no_ofst[2*i,3]<-ifelse(anva.lrt$`Pr(>Chi)`[2]<0.001,"***",
                                ifelse((anva.lrt$`Pr(>Chi)`[2]<0.01)&(anva.lrt$`Pr(>Chi)`[2]>0.001),"**",
                                       ifelse((anva.lrt$`Pr(>Chi)`[2]<0.05)&(anva.lrt$`Pr(>Chi)`[2]>0.01),"*","")))
}
View(lrt.qd.no_ofst)

##################################################################################
# Code for generating Relative intensity plots is forthcoming 
# It is currently being exported from a secure environment 
##################################################################################

##################################################################################
# Marginal effect analysis (and plots)
min.hsi.stn<-min(hsi.resc.stn)
max.hsi.stn<-max(hsi.resc.stn)
hsi.vals<-seq(min.hsi.stn,max.hsi.stn,0.5)*sd(hsi.im)+mean(hsi.im)#For HSI ME plotting 
NDVI.vals<-seq(-0.5,1,0.1)


ymax=0.002
ef.fun<-list()
ef.fun.ndvi.dif<-list()
ef.fun.hsi<-list()

par(mfrow=c(2,3))
for (i in 1:11){
  ### NDVIME
  ef.fun[[2*i-1]]<-effectfun(fl.env.qd.no.dup.s1[[i]], 
                             covname = "NDVI_cur.stan",
                             NDVI.diff.stan= quantile.im(NDVI.d.rsc.stn[[2*i-1]])[[4]],
                             hsi.stan=quantile.im(hsi.resc.stn)[[4]],
                             sm.mpd.ofst.stan = quantile.im(blr.mpd.resc.ofst.stan[[i]])[[4]],
                             lnd.ch = lch.resc[[i]]$v[which(lch.resc[[i]]$v=="No_Change")[1]],
                             lnd.cov=lcov.resc[[i]]$v[which(lcov.resc[[i]]$v=="Other")[1]],
                             marks="under_45_Male", se.fit = TRUE, ylim=c(0,ymax))
  plot(ef.fun[[2*i-1]], xaxt='n',main=""
       ,ylab="Predicted intensity", 
       xlab="NDVI: Winter/Spring", ylim=c(0,ymax),
       legend=FALSE,xlim=c(-5,5), col="red",shadecol=alpha("pink",0.4))
  axis(side=1, at=seq(-5,5,length.out=16), labels=round(NDVI.vals,2))
  
  ### NDVI Diff MED
  ef.fun.ndvi.dif[[2*i-1]]<-effectfun(fl.env.qd.no.dup.s1[[i]], 
                                      covname = "NDVI.diff.stan",
                                      NDVI_cur.stan= quantile.im(NDVI.ws.rsc.stn[[i+1]])[[4]],
                                      hsi.stan=quantile.im(hsi.resc.stn)[[4]],
                                      sm.mpd.ofst.stan = quantile.im(blr.mpd.resc.ofst.stan[[i]])[[4]],
                                      lnd.ch = lch.resc[[i]]$v[which(lch.resc[[i]]$v=="No_Change")[1]],
                                      lnd.cov=lcov.resc[[i]]$v[which(lcov.resc[[i]]$v=="Other")[1]],
                                      marks="under_45_Male", se.fit = TRUE)
  
  
  plot(ef.fun.ndvi.dif[[2*i-1]], xaxt='n',main="",
       ylab="Predicted intensity", 
       xlab="NDVI Diff.",  ylim=c(0,ymax),
       legend=FALSE, xlim=c(-5,5),col="blue",shadecol=alpha("lightblue",0.4))
  
  axis(side=1, at=seq(-5,5,length.out=16), labels=round(NDVI.vals,2))
  
  ##### HSI ME
  ef.fun.hsi[[2*i-1]]<-effectfun(fl.env.qd.no.dup.s1[[i]], 
                                 covname = "hsi.stan",
                                 NDVI.diff.stan= quantile.im(NDVI.d.rsc.stn[[2*i-1]])[[4]],
                                 NDVI_cur.stan= quantile.im(NDVI.ws.rsc.stn[[i+1]])[[4]],
                                 sm.mpd.ofst.stan = quantile.im(blr.mpd.resc.ofst.stan[[i]])[[4]],
                                 lnd.ch = lch.resc[[i]]$v[which(lch.resc[[i]]$v=="No_Change")[1]],
                                 lnd.cov=lcov.resc[[i]]$v[which(lcov.resc[[i]]$v=="Other")[1]],
                                 marks="under_45_Male", se.fit = TRUE)
  
  plot(ef.fun.hsi[[2*i-1]], main="",xaxt='n',
       ylab="Predicted intensity", ylim=c(0,ymax),
       xlab="HSI", 
       legend=FALSE, xlim=c(min.hsi.stn,max.hsi.stn), col="forestgreen",shadecol=alpha("lightgreen",0.4))
  
  axis(side=1, at=seq(min.hsi.stn,max.hsi.stn,0.5), labels=round(hsi.vals,2))
  
  
  ###########################################
  ####Sem 2 same covs as above and same order
  ### NDVI ME
  ef.fun[[2*i]]<-effectfun(fl.env.qd.no.dup.s2[[i]], 
                           covname = "NDVI_cur.stan",
                           NDVI.diff.stan= quantile.im(NDVI.d.rsc.stn[[2*i]])[[4]],
                           hsi.stan=quantile.im(hsi.resc.stn)[[4]],
                           sm.mpd.ofst.stan = quantile.im(blr.mpd.resc.ofst.stan[[i]])[[4]],
                           lnd.ch = lch.resc[[i]]$v[which(lch.resc[[i]]$v=="No_Change")[1]],
                           lnd.cov=lcov.resc[[i]]$v[which(lcov.resc[[i]]$v=="Other")[1]],
                           marks="under_45_Male", se.fit = TRUE)
  plot(ef.fun[[2*i]], xaxt='n',main=""
       ,ylab="Predicted intensity", 
       xlab="NDVI: Summer/Fall",  ylim=c(0,ymax),
       legend=FALSE,xlim=c(-5,5), col="red",shadecol=alpha("pink",0.4))
  axis(side=1, at=seq(-5,5,length.out=16), labels=round(NDVI.vals,2))
  
  ### NDVI Diff ME
  ef.fun.ndvi.dif[[2*i]]<-effectfun(fl.env.qd.no.dup.s2[[i]], 
                                    covname = "NDVI.diff.stan",
                                    NDVI_cur.stan= quantile.im(NDVI.sf.rsc.stn[[i+1]])[[4]],
                                    hsi.stan=quantile.im(hsi.resc.stn)[[4]],
                                    sm.mpd.ofst.stan = quantile.im(blr.mpd.resc.ofst.stan[[i]])[[4]],
                                    lnd.ch = lch.resc[[i]]$v[which(lch.resc[[i]]$v=="No_Change")[1]],
                                    lnd.cov=lcov.resc[[i]]$v[which(lcov.resc[[i]]$v=="Other")[1]],
                                    marks="under_45_Male", se.fit = TRUE)
  
  
  plot(ef.fun.ndvi.dif[[2*i]], xaxt='n',main="",
       ylab="Predicted intensity", 
       xlab="NDVI Diff.", ylim=c(0,ymax),
       legend=FALSE, xlim=c(-5,5),col="blue",shadecol=alpha("lightblue",0.4))
  
  axis(side=1, at=seq(-5,5,length.out=16), labels=round(NDVI.vals,2))
  
  ###### hsi ME
  ef.fun.hsi[[2*i]]<-effectfun(fl.env.qd.no.dup.s2[[i]], 
                               covname = "hsi.stan",
                               NDVI.diff.stan= quantile.im(NDVI.d.rsc.stn[[2*i]])[[4]],
                               NDVI_cur.stan= quantile.im(NDVI.sf.rsc.stn[[i+1]])[[4]],
                               sm.mpd.ofst.stan = quantile.im(blr.mpd.resc.ofst.stan[[i]])[[4]],
                               lnd.ch = lch.resc[[i]]$v[which(lch.resc[[i]]$v=="No_Change")[1]],
                               lnd.cov=lcov.resc[[i]]$v[which(lcov.resc[[i]]$v=="Other")[1]],
                               marks="under_45_Male", se.fit = TRUE)
  
  plot(ef.fun.hsi[[2*i]], xaxt='n',main="",
       ylab="Predicted intensity", 
       xlab="HSI", ylim=c(0,ymax),
       legend=FALSE, xlim=c(min.hsi.stn,max.hsi.stn), col="forestgreen",shadecol=alpha("lightgreen",0.4))
  
  axis(side=1, at=seq(min.hsi.stn,max.hsi.stn,by=0.5), labels=round(hsi.vals,2))
  
}
################# Relative intensity analysis

test.den.qd.mod<-list()
test.den.ofst.mod<-list()
test.den.lc.mod<-list()
test.den.lin.mod<-list()
pct.removed<-data.frame(matrix(ncol=4,nrow=22))
names(pct.removed)<-c("quadratic", "linear","land cover","offset")

vals<-seq(0,20, length.out=256)
white_pos<-which.min(abs(vals-1.5))
n1<-white_pos
cols1<-colorRampPalette(c("blue","white"))(n1)
n2<-256-n1
cols2<-colorRampPalette(c("white","darkred"))(n2)
my_colors<-c(cols1,cols2)

gd.pixel.pct<-data.frame(matrix(ncol=4,nrow=22)) #Create a data frame to store the pixels between 0.5 and 2
names(gd.pixel.pct)<-c("quadratic", "linear","land cover","offset") 
pct.removed<-data.frame(matrix(ncol=4,nrow=22)) #Create a data frame to store the percentage of points removed as outliers
names(pct.removed)<-c("quadratic", "linear","land cover","offset")

par(mfrow=c(2,2))
par(mar=c(0,0,0,1))
epsilon=0.005 #Threshold for removing data with fitted intensities that are outliers.
for (i in 1:11){
  ##################################
  #### quadratic full model Semester 1
  pred.int<-predict.ppm(fl.env.qd.no.dup.s1[[i]], locations=fl.env.qd.no.dup.s1[[i]]$Q$data)#locations=subset.ppp(test.pp,as.owin(test.mpd.im)))
  pct.removed[2*i-1,1]<-length(pred.int[pred.int<epsilon])/fl.env.qd.no.dup.s1[[i]]$Q$data$n
  pred.int<-ifelse(pred.int<epsilon,min(pred.int[pred.int>epsilon|pred.int==epsilon], na.rm=TRUE),pred.int)
  test.den.qd.mod[[2*i-1]]<-adaptive.density(fl.env.qd.no.dup.s1[[i]]$Q$data, method="kernel", weights=1/pred.int, ho=5, dimyx=250, edge=TRUE)

  adap.den<- test.den.qd.mod[[2*i-1]]
  gd.pixel.pct[2*i-1,1]<-sum((adap.den$v<2)&(adap.den$v>0.5), na.rm = T)/length(adap.den$v)

  adap.den$v<-ifelse((adap.den$v<2)&(adap.den$v>0.5),NA,adap.den$v)
  plot(adap.den, col=alpha(my_colors,1), main="",zlim=c(0,20))
  plot(phx_mask_resc, add=TRUE)
  plot(adap.den, col=alpha(my_colors,0.7), add=TRUE,zlim=c(0,20))
  ################################# Linear model sem 1
  pred.int<-predict.ppm(base.mod.no.dup.s1[[i]], locations=base.mod.no.dup.s1[[i]]$Q$data)#locations=subset.ppp(test.pp,as.owin(test.mpd.im)))
  pct.removed[2*i-1,2]<-length(pred.int[pred.int<epsilon])/base.mod.no.dup.s1[[i]]$Q$data$n
  pred.int<-ifelse(pred.int<epsilon,min(pred.int[pred.int>epsilon|pred.int==epsilon], na.rm=TRUE),pred.int)
  test.den.lin.mod[[2*i-1]]<-adaptive.density(base.mod.no.dup.s1[[i]]$Q$data, method="kernel", weights=1/pred.int, ho=5, dimyx=250, edge=TRUE)

  adap.den<- test.den.lin.mod[[2*i-1]]
  gd.pixel.pct[2*i-1,2]<-sum((adap.den$v<2)&(adap.den$v>0.5), na.rm = T)/length(adap.den$v)

  adap.den$v<-ifelse((adap.den$v<2)&(adap.den$v>0.5),NA,adap.den$v)
  plot(adap.den, col=alpha(my_colors,1), main="",zlim=c(0,20))
  plot(phx_mask_resc, add=TRUE)
  plot(adap.den, col=alpha(my_colors,0.7), add=TRUE,zlim=c(0,20))

  ################################# Land cover model sem 1
  pred.int<-predict.ppm(ofst.NEcovs.s1[[i]], locations=ofst.NEcovs.s1[[i]]$Q$data)#locations=subset.ppp(test.pp,as.owin(test.mpd.im)))
  pct.removed[2*i-1,3]<-length(pred.int[pred.int<epsilon])/ofst.NEcovs.s1[[i]]$Q$data$n
  pred.int<-ifelse(pred.int<epsilon,min(pred.int[pred.int>epsilon|pred.int==epsilon], na.rm=TRUE),pred.int)
  test.den.lc.mod[[2*i-1]]<-adaptive.density(ofst.NEcovs.s1[[i]]$Q$data, method="kernel", weights=1/pred.int, ho=5, dimyx=250, edge=TRUE)

  adap.den<- test.den.lc.mod[[2*i-1]]
  gd.pixel.pct[2*i-1,3]<-sum((adap.den$v<2)&(adap.den$v>0.5), na.rm = T)/length(adap.den$v)

  adap.den$v<-ifelse((adap.den$v<2)&(adap.den$v>0.5),NA,adap.den$v)
  plot(adap.den, col=alpha(my_colors,1), main="",zlim=c(0,20))
  plot(phx_mask_resc, add=TRUE)
  plot(adap.den, col=alpha(my_colors,0.7), add=TRUE,zlim=c(0,20))

  ##################################
  #### Offset Only Model semester 1
  pred.int<-predict.ppm(ofst.only.s1[[i]], locations=ofst.only.s1[[i]]$Q$data)
  pct.removed[2*i-1,4]<-length(pred.int[pred.int<epsilon])/ofst.only.s1[[i]]$Q$data$n
  pred.int<-ifelse(pred.int<epsilon,min(pred.int[pred.int>epsilon|pred.int==epsilon], na.rm=TRUE),pred.int)

  test.den.ofst.mod[[2*i-1]]<-adaptive.density(ofst.only.s1[[i]]$Q$data, method="kernel", weights=1/pred.int, ho=5, dimyx=250, edge=TRUE)

  adap.den<- test.den.ofst.mod[[2*i-1]]
  gd.pixel.pct[2*i-1,4]<-sum((adap.den$v<2)&(adap.den$v>0.5), na.rm = T)/length(adap.den$v)

  adap.den$v<-ifelse((adap.den$v<2)&(adap.den$v>0.5),NA,adap.den$v)
  plot(adap.den, col=alpha(my_colors,1), main="",zlim=c(0,20))
  plot(phx_mask_resc, add=TRUE)
  plot(adap.den, col=alpha(my_colors,0.7), add=TRUE,zlim=c(0,20))

  ###########################
  ###########################
  ###########################
  ##########################
  ##################################
  #### quadratic full model Semester 2
  pred.int<-predict.ppm(fl.env.qd.no.dup.s2[[i]], locations=fl.env.qd.no.dup.s2[[i]]$Q$data)#locations=subset.ppp(test.pp,as.owin(test.mpd.im)))
  pct.removed[2*i,1]<-length(pred.int[pred.int<epsilon])/fl.env.qd.no.dup.s2[[i]]$Q$data$n
  pred.int<-ifelse(pred.int<epsilon,min(pred.int[pred.int>epsilon|pred.int==epsilon], na.rm=TRUE),pred.int)
  test.den.qd.mod[[2*i]]<-adaptive.density(fl.env.qd.no.dup.s2[[i]]$Q$data, method="kernel", weights=1/pred.int, ho=5, dimyx=250, edge=TRUE)

  adap.den<- test.den.qd.mod[[2*i]]
  gd.pixel.pct[2*i,1]<-sum((adap.den$v<2)&(adap.den$v>0.5), na.rm = T)/length(adap.den$v)
  adap.den$v<-ifelse((adap.den$v<2)&(adap.den$v>0.5),NA,adap.den$v)
  plot(adap.den, col=alpha(my_colors,1), main="",zlim=c(0,20))
  plot(phx_mask_resc, add=TRUE)
  plot(adap.den, col=alpha(my_colors,0.7), add=TRUE,zlim=c(0,20))

  ################################# Linear model sem 2
  pred.int<-predict.ppm(base.mod.no.dup.s2[[i]], locations=base.mod.no.dup.s2[[i]]$Q$data)#locations=subset.ppp(test.pp,as.owin(test.mpd.im)))
  pct.removed[2*i,2]<-length(pred.int[pred.int<epsilon])/base.mod.no.dup.s2[[i]]$Q$data$n
  pred.int<-ifelse(pred.int<epsilon,min(pred.int[pred.int>epsilon|pred.int==epsilon], na.rm=TRUE),pred.int)
  test.den.lin.mod[[2*i]]<-adaptive.density(base.mod.no.dup.s2[[i]]$Q$data, method="kernel", weights=1/pred.int, ho=5, dimyx=250, edge=TRUE)

  adap.den<- test.den.lin.mod[[2*i]]
  gd.pixel.pct[2*i,2]<-sum((adap.den$v<2)&(adap.den$v>0.5), na.rm = T)/length(adap.den$v)

  adap.den$v<-ifelse((adap.den$v<2)&(adap.den$v>0.5),NA,adap.den$v)
  plot(adap.den, col=alpha(my_colors,1), main="",zlim=c(0,20))
  plot(phx_mask_resc, add=TRUE)
  plot(adap.den, col=alpha(my_colors,0.7), add=TRUE,zlim=c(0,20))
  # ################################# Land cover model sem 2
  pred.int<-predict.ppm(ofst.NEcovs.s2[[i]], locations=ofst.NEcovs.s2[[i]]$Q$data)
  pct.removed[2*i,3]<-length(pred.int[pred.int<epsilon])/ofst.NEcovs.s2[[i]]$Q$data$n
  pred.int<-ifelse(pred.int<epsilon,min(pred.int[pred.int>epsilon|pred.int==epsilon], na.rm=TRUE),pred.int)
  test.den.lc.mod[[2*i]]<-adaptive.density(ofst.NEcovs.s2[[i]]$Q$data, method="kernel", weights=1/pred.int, ho=5, dimyx=250, edge=TRUE)

  adap.den<- test.den.lc.mod[[2*i]]
  gd.pixel.pct[2*i,3]<-sum((adap.den$v<2)&(adap.den$v>0.5), na.rm = T)/length(adap.den$v)

  adap.den$v<-ifelse((adap.den$v<2)&(adap.den$v>0.5),NA,adap.den$v)
  plot(adap.den, col=alpha(my_colors,1), main="",zlim=c(0,20))
  plot(phx_mask_resc, add=TRUE)
  plot(adap.den, col=alpha(my_colors,0.7), add=TRUE,zlim=c(0,20))

  ##################################
  #### Offset Only Model semester 2
  pred.int<-predict.ppm(ofst.only.s2[[i]], locations=ofst.only.s2[[i]]$Q$data)
  pct.removed[2*i,4]<-length(pred.int[pred.int<epsilon])/ofst.only.s2[[i]]$Q$data$n
  pred.int<-ifelse(pred.int<epsilon,min(pred.int[pred.int>epsilon|pred.int==epsilon], na.rm=TRUE),pred.int)
  test.den.ofst.mod[[2*i]]<-adaptive.density(ofst.only.s2[[i]]$Q$data, method="kernel", weights=1/pred.int, ho=5, dimyx=250,edge=TRUE)

  adap.den<- test.den.ofst.mod[[2*i]]
  gd.pixel.pct[2*i,4]<-sum((adap.den$v<2)&(adap.den$v>0.5), na.rm = T)/length(adap.den$v)
  
  adap.den$v<-ifelse((adap.den$v<2)&(adap.den$v>0.5),NA,adap.den$v)
  plot(adap.den, col=alpha(my_colors,1), main="",zlim=c(0,20))
  plot(phx_mask_resc, add=TRUE)
  plot(adap.den, col=alpha(my_colors,0.7), add=TRUE,zlim=c(0,20))
  
}