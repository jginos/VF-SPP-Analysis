library("readxl")
library("terra")
library("sp")
library("spatstat")
library("sf")
library("RColorBrewer")
library("purrr")
library("scales")
library("spatstat.model")

load("/vfdf.Rdata") #Load simulated data
#vfdf is the simulated data set for each semester over the 22 semesters and in the window created from the window phx_W_sf
# Keep in mind that this is a simulated data set and may contain points in the study window where observed VF cases are not possible, these should be handled accordingly.

# Read in Phoenix raster and convert it to a spatstat owin object
phx_osm<-rast("\\OSM_Phx_L_image.tif")
phx_W_sf<-st_as_sf(st_as_sfc(st_bbox(phx_osm)))
phx_W<-as.owin(phx_W_sf)


#######################################################################
# Read in raster data for modeling point pattern as a function of covariates
# Next: Convert rasters to pixel images for use with Spatstat functions
# We need the following function for converting rasters to pixel images
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
hsi.raster<-rast("//HSI_Raster_10m_0_to_1.tif")
hsi.raster<-crop(hsi.raster,phx_W_sf)
values(hsi.raster)<-ifelse(values(hsi.raster)>1,NA,values(hsi.raster))
hsi.test<-hsi.raster
## We set Certain unpublished areas to NA that were erroneously recorded as 0 in the original data.
## Also, we set values above 1, which include areas that are mostly likely water, equal to NA 
#Topright major subregion
bbox<-ext(-12473000,-12410005,3913118,4015610)

subreg<-crop(hsi.test,bbox)
values(subreg)<-ifelse(values(subreg)==0,NA,values(subreg))
#plot(subreg)
subreg.coords<-xyFromCell(subreg, 1:ncell(subreg))
new_cells<-cellFromXY(hsi.test, subreg.coords)
hsi.test[new_cells]<-values(subreg)

#bottom left subregion
bbox<-ext(-12560836,-12540000,3913118,3925500)

subreg<-crop(hsi.test,bbox)
values(subreg)<-ifelse(values(subreg)==0,NA,values(subreg))
#plot(subreg)
subreg.coords<-xyFromCell(subreg, 1:ncell(subreg))
new_cells<-cellFromXY(hsi.test, subreg.coords)
hsi.test[new_cells]<-values(subreg)


hsi.im<-as.im.SpatRaster2(hsi.test)

#standardized HSI covariate
hsi.stan<-hsi.im

hsi.stan<-(hsi.stan-mean(hsi.stan,na.rm=T))/sd(hsi.stan, na.rm=TRUE)

#Rescale to Kilometers
hsi.resc.stn_1km<-rescale.im(hsi.stan, s=1000, unitname = "Kilometers")


hsi.resc.stn_1km<- blur(hsi.resc.stn, sigma=1, bleed=FALSE, normalise = TRUE)




#Create pixel images of NDVI, land cover, and land cover change covariates
# read in NDVI data and projecting it 
NDVI.dat<-list() #NOTE: we are including 2012 NDVI for lag effect
for (i in 1:12){
  NDVI.dat[[i]]<-rast(paste0("//knb-lter-cap//NDVI_multiseason_CAPLTER_",2011+i,".tif"))
  NDVI.dat[[i]]<-project(NDVI.dat[[i]], "epsg:3857", method= "cubicspline")
}

NDVI.ws.im.wna<-list()
for (i in 1:length(NDVI.dat)){
  avg.NDVI<- (NDVI.dat[[i]]$`1_winter`+NDVI.dat[[i]]$`2_spring`)/2
  NDVI.ws.im.wna[[i]]<-as.im.SpatRaster2(avg.NDVI)
}

#### standardize NDVI.ws.im.stan
NDVI.ws.im.stan<-NDVI.ws.im.wna
for (i in 1:length(NDVI.ws.im.stan)){
  NDVI.ws.im.stan[[i]]$v<-(NDVI.ws.im.stan[[i]]$v-mean(NDVI.ws.im.stan[[i]]$v,na.rm=T))/sd(NDVI.ws.im.stan[[i]]$v, na.rm=T)
}

# Rescale to Kilometers

NDVI.ws.rsc.stn_1km<-list()

for(i in 1:length(NDVI.ws.im.stan)){
  NDVI.ws.rsc.stn_1km[[i]]<-rescale.im(NDVI.ws.im.stan[[i]], s=1000, unitname = "Kilometers")
  NDVI.ws.rsc.stn_1km[[i]]<-blur(NDVI.ws.rsc.stn_1km[[i]], sigma=1,bleed=FALSE,normalise=TRUE)
  
}

###### NDVI.sf with NA

NDVI.sf.im.wna<-list()
for (i in 1:length(NDVI.dat)){
  avg.NDVI<- (NDVI.dat[[i]]$`3_summer`+NDVI.dat[[i]]$`4_fall`)/2
  NDVI.sf.im.wna[[i]]<-as.im.SpatRaster2(avg.NDVI)
}

#### standardize NDVI.sf.im.stan
NDVI.sf.im.stan<-NDVI.sf.im.wna
for (i in 1:length(NDVI.sf.im.stan)){
  NDVI.sf.im.stan[[i]]<-(NDVI.sf.im.stan[[i]]-mean(NDVI.sf.im.stan[[i]],na.rm=T))/sd(NDVI.sf.im.stan[[i]], na.rm=T)
}

#Rescale to Kilometers
NDVI.sf.rsc.stn_1km<-list()
for(i in 1:length(NDVI.sf.im.stan)){
  NDVI.sf.rsc.stn_1km[[i]]<-rescale.im(NDVI.sf.im.stan[[i]], s=1000, unitname = "Kilometers")
  NDVI.sf.rsc.stn_1km[[i]]<-blur(NDVI.sf.rsc.stn_1km[[i]], sigma=1,bleed=FALSE,normalise=TRUE)
  
}

# Estimate intensity of VF cases as a function of a particular covariate (mnarginally)
rho_hat_NDVI_s1_1km<-list()
rho_hat_NDVI_s1_lag1_1km<-list()
rho_hat_NDVI_s1_lag2_1km<-list()
rho_hat_NDVI_s2_1km<-list()
rho_hat_NDVI_s2_lag1_1km<-list()
rho_hat_hsi_s1_1km<-list()
rho_hat_hsi_s2_1km<-list()
rho_hat_lc_HI_s1_1km<-list()
rho_hat_lc_HI_s2_1km<-list()
rho_hat_lc_MI_s1_1km<-list()
rho_hat_lc_MI_s2_1km<-list()
rho_hat_lc_LI_s1_1km<-list()
rho_hat_lc_LI_s2_1km<-list()
rho_hat_lc_open_s1_1km<-list()
rho_hat_lc_open_s2_1km<-list()
rho_hat_lchange_s1_1km<-list()
rho_hat_lchange_s2_1km<-list()

for (i in 1:11){
  # semester 1
  rho_hat_NDVI_s1_1km[[i]]<-rhohat(unmark(vfdf[[2*i-1]]),
                          as.im(NDVI.sf.rsc.stn_1km[[i+1]], dimyx=2500))

  rho_hat_NDVI_s1_lag1_1km[[i]]<-rhohat(unmark(vfdf[[2*i-1]]),
                               as.im(NDVI.ws.rsc.stn_1km[[i+1]], dimyx=2500)) #lag 1
  
  rho_hat_NDVI_s1_lag2_1km[[i]]<-rhohat(unmark(vfdf[[2*i-1]]),
                               as.im(NDVI.sf.rsc.stn_1km[[i]], dimyx=2500)) #lag 2
  # semester 2
  rho_hat_hsi_s1_1km[[i]]<-rhohat(unmark(vfdf[[2*i-1]]),
                              as.im(hsi.resc.stn_1km, dimyx=2500))
  # Lcover High
  rho_hat_lc_HI_s1_1km[[i]]<-rhohat(unmark(vfdf[[2*i-1]]),
                                    as.im(lc.dev.HI_1km_smth[[i]], dimyx=2500))
  
  rho_hat_lc_MI_s1_1km[[i]]<-rhohat(unmark(vfdf[[2*i-1]]),
                                    as.im(lc.dev.MI_1km_smth[[i]], dimyx=2500))
  
  rho_hat_lc_LI_s1_1km[[i]]<-rhohat(unmark(vfdf[[2*i-1]]),
                                    as.im(lc.dev.LI_1km_smth[[i]], dimyx=2500))
  
  rho_hat_lc_open_s1_1km[[i]]<-rhohat(unmark(vfdf[[2*i-1]]),
                                    as.im(lc.dev.open_1km_smth[[i]], dimyx=2500))
  
  rho_hat_lchange_s1_1km[[i]]<-rhohat(unmark(vfdf[[2*i-1]]),
                                      as.im(lch.resc_1km_smth[[i]], dimyx=2500))
  
  rho_hat_NDVI_s2_1km[[i]]<-rhohat(unmark(vfdf[[2*i]]),
                          as.im(NDVI.ws.rsc.stn_1km[[i+1]], dimyx=2500))
  
  rho_hat_NDVI_s2_lag1_1km[[i]]<-rhohat(unmark(vfdf[[2*i]]),
                               as.im(NDVI.sf.rsc.stn_1km[[i]], dimyx=2500)) #lag 1

  
  rho_hat_NDVI_s2_lag2_1km[[i]]<-rhohat(unmark(vfdf[[2*i]]),
                               as.im(NDVI.ws.rsc.stn_1km[[i]], dimyx=2500)) #lag 2

  rho_hat_hsi_s2_1km[[i]]<-rhohat(unmark(vfdf[[2*i]]),
                              as.im(hsi.resc.stn_1km, dimyx=2500))
  
  rho_hat_lc_HI_s2_1km[[i]]<-rhohat(unmark(vfdf[[2*i]]),
                                    as.im(lc.dev.HI_1km_smth[[i]], dimyx=2500))
  
  rho_hat_lc_MI_s2_1km[[i]]<-rhohat(unmark(vfdf[[2*i]]),
                                    as.im(lc.dev.MI_1km_smth[[i]], dimyx=2500))
  
  rho_hat_lc_LI_s2_1km[[i]]<-rhohat(unmark(vfdf[[2*i]]),
                                    as.im(lc.dev.LI_1km_smth[[i]], dimyx=2500))
  
  rho_hat_lc_open_s2_1km[[i]]<-rhohat(unmark(vfdf[[2*i]]),
                                      as.im(lc.dev.open_1km_smth[[i]], dimyx=2500))
  
  rho_hat_lchange_s2_1km[[i]]<-rhohat(unmark(vfdf[[2*i]]),
                                      as.im(lch.resc_1km_smth[[i]], dimyx=2500))
}



lnd.cov.ch.Dev<-list()
lc.Dev.categ<-c(5121:5124,5221:5224) # the values here are

for (i in 1:4){
  lnd.cov.ch.Dev[[i]]<-rast(paste0("\\NLCD_g7SX3oX8v6FFj9gn0jdV\\Annual_NLCD_LndChg_",2012+i,"_CU_C1V0_g7SX3oX8v6FFj9gn0jdV.tiff"))
  vals<-as.character(values(lnd.cov.ch.Dev[[i]]))
  values(lnd.cov.ch.Dev[[i]])<-ifelse((values(lnd.cov.ch.Dev[[i]])%in%lc.Dev.categ) ,1,0)
  lnd.cov.ch.Dev[[i]]<-project(lnd.cov.ch.Dev[[i]],crs(AZ_Cities), method="near")
  #lnd.cov.ch.Dev[[i]]<-aggregate(lnd.cov.ch.Dev[[i]], fact=1000/res(lnd.cov.ch.Dev[[i]])[1], fun="any", na.rm=TRUE)
}


# We are going to smooth the changed land cover as a covariate that is continuos rather than binary
lnd.ch.Dev.im<-list()

for (i in 1:length(lnd.cov.ch.Dev)){
  lnd.ch.Dev.im[[i]]<-as.im.SpatRaster2(lnd.cov.ch.Dev[[i]])
  
}


#lch.resc_1km_smth_1km_smth<-list()

# Sigma should be either 1 or 0.5; for consistentcy with last draft, we keep it at 1
for(i in 1:4){
  lch.resc_1km_smth[[i]]<-rescale.im(lnd.ch.Dev.im[[i]], s=1000, unitname = "Kilometers")
  lch.resc_1km_smth[[i]]<-blur(lch.resc_1km_smth_1km_smth[[i]], sigma=1,bleed=FALSE,normalise = TRUE)
}

for (i in 1:11){
  lch.resc_1km_smth[[i]]<-lch.resc_1km_smth[[i]][phx_W.resc, drop=FALSE]
}

plot(test)
plot(NDVI.ws.rsc.stn_1km[[1]],add=TRUE)
test<-subset(test,lc.dev.LI_1km_smth_1km_smth[[1]])
#Projecting first and then aggregating did not produce a large difference compared to aggregating fist and then projecting the raster (in the above case for 2013 and 2014) 
# Moreover, it seems intuitively reasonable to project first so that aggregation with respect to resolution takes place in the projected coordinate system used for analysis

lc.dev.open_1km_smth_1km<-list()
lc.dev.LI_1km_smth_1km<-list()
lc.dev.MI_1km_smth_1km<-list()
lc.dev.HI_1km_smth_1km<-list()

for (i in 1:11){
  lnd.cov<-rast(paste0("\\NLCD_g7SX3oX8v6FFj9gn0jdV\\Annual_NLCD_LndCov_",2012+i,"_CU_C1V0_g7SX3oX8v6FFj9gn0jdV.tiff"))
  lnd.cov<-project(lnd.cov,crs(AZ_Cities), method="near")
  
  lc.dev.open_1km_smth_1km[[i]]<-lnd.cov
  values(lc.dev.open_1km_smth_1km[[i]]) <-ifelse(values(lnd.cov)==21,1,0)
  
  lc.dev.LI_1km_smth_1km[[i]]<-lnd.cov
  values(lc.dev.LI_1km_smth_1km[[i]]) <-ifelse(values(lnd.cov)==22,1,0)
  
  lc.dev.MI_1km_smth_1km_smth[[i]]<-lnd.cov
  values(lc.dev.MI_1km_smth_1km_smth[[i]]) <-ifelse(values(lnd.cov)==23,1,0)
  
  lc.dev.HI_1km_smth_1km_smth[[i]]<-lnd.cov
  values(lc.dev.HI_1km_smth_1km_smth[[i]]) <-ifelse(values(lnd.cov)==24,1,0)
  # 
  # lc.dev.other[[i]]<-lnd.cov
  # values(lc.dev.other[[i]]) <-ifelse(values(lnd.cov) %in% seq(21,24,1),0,1)
}


sigma=1
for (i in 1:11){
  lc.dev.open_1km_smth[[i]]<-as.im.SpatRaster2(lc.dev.open_1km_smth[[i]])
  lc.dev.open_1km_smth[[i]]<-rescale.im(lc.dev.open_1km_smth[[i]], s=1000, unitname = "Kilometers")
  lc.dev.open_1km_smth[[i]]<-blur(lc.dev.open_1km_smth[[i]], sigma=sigma,bleed=FALSE,normalise = TRUE)
  
  lc.dev.LI_1km_smth[[i]]<-as.im.SpatRaster2(lc.dev.LI_1km_smth[[i]])
  lc.dev.LI_1km_smth[[i]]<-rescale.im(lc.dev.LI_1km_smth[[i]], s=1000, unitname = "Kilometers")
  lc.dev.LI_1km_smth[[i]]<-blur(lc.dev.LI_1km_smth[[i]], sigma=sigma,bleed=FALSE,normalise = TRUE)
  
  lc.dev.MI_1km_smth[[i]]<-as.im.SpatRaster2(lc.dev.MI_1km_smth[[i]])
  lc.dev.MI_1km_smth_1km_smth[[i]]<-rescale.im(lc.dev.MI_1km_smth[[i]], s=1000, unitname = "Kilometers")
  lc.dev.MI_1km_smth[[i]]<-blur(lc.dev.MI_1km_smth[[i]], sigma=sigma,bleed=FALSE,normalise = TRUE)
  
  lc.dev.HI_1km_smth[[i]]<-as.im.SpatRaster2(lc.dev.HI_1km_smth[[i]])
  lc.dev.HI_1km_smth[[i]]<-rescale.im(lc.dev.HI_1km_smth[[i]], s=1000, unitname = "Kilometers")
  lc.dev.HI_1km_smth[[i]]<-blur(lc.dev.HI_1km_smth[[i]], sigma=sigma,bleed=FALSE,normalise = TRUE)
  
}


]#####################################################################
# Read in and prepare census data for use in the medicaid population density covariate

#Read in medicaid ACTUAL data by county/month from AZHCCCS
Med.dat.actual<-read_xlsx("\\AHCCCS By County 2013-2024.xlsx", sheet = "Sheet1")
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
  med.dat[[i]]<-as.data.frame(read.csv(paste0("\\ACSDT5Y",2012+i,".B27010-Data.csv")))
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
bg.shp.2010<-st_read(dsn = "\\tl_2019_04_bg")
bg.shp.2020<-st_read(dsn = "\\AZ_2022_Block_Group_Shapefile")
AZ_Counties<-st_read(dsn="\\AZ_Counties_SHP")
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
bg.shp.2010$GEOID<- paste0(0,bg.shp.2010$GEOID)

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

AZ_Metro.bg.2010<-st_intersection(phx_W_sf, bg.medicaid.2010) #combining cities shapefile with 2010 block groups
AZ_Metro.bg.2020<-st_intersection(phx_W_sf, bg.medicaid.2020) #combining cities shapefile with 2020 block groups
# Create medicaid population density covariate that is a (rasterized, and then converted to the) pixel-image version of the shapefile 
# with pixels at the 1000-meter resolution
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

v2<-vect(AZ_Metro.bg.2020)
y2<-rast(v,res=1000)

for(i in 8:11){
  r<-rasterize(x=v2,y=y2, field=paste0("Med_Pop_Den",2012+i), fun="mean")
  med.pop.den.im[[i]]<-as.im.SpatRaster2(r) #Creates 'im' object with medicaid by bg information
}

#Smoothing the medicaid population density raster

#### Create an offset that is the pixel-sum-standardized (divide each pixel by the pixel sum for the whole image) medicaid population density.
blr.mpd.resc.ofst.stan<-list()
plot(blr.mpd.resc.ofst.stan[[1]])
for(i in 1:11){
  blr.mpd.resc.ofst.stan[[i]]<-rescale.im(med.pop.den.im[[i]], s=1000, unitname = "Kilometers")
  blr.mpd.resc.ofst.stan[[i]]<-blur(blr.mpd.resc.ofst.stan[[i]], bleed=FALSE, normalise = TRUE,sigma=2)
  blr.mpd.resc.ofst.stan[[i]]<-blr.mpd.resc.ofst.stan[[i]]/sum(blr.mpd.resc.ofst.stan[[i]], na.rm=TRUE)
}

############################################################################################################
####################### Offset only models
# a function for cleaning up NA's in covaraites before creating quadscheme for fitting:

trim_quad_na_any<-function(Q,covs){ # this function was generated by AI (ChatGPT )
  X<-union.quad(Q)
  idat<-is.data(Q)
  
  keep<- !is.na(lookup.im(covs[[1]],X$x,X$y))
  if (length(covs)>1L) {
    for (k in 2:length(covs)) {
      keep<- keep & !is.na(lookup.im(covs[[k]],X$x,X$y))
    }
  }
  Q$data<-as.ppp(X[idat & keep ])
  Q$dummy<-X[(!idat) & keep]
  Q$w<-Q$w[keep] #[!idat]
  
  Q
  
}
ofst.only.s1[[i]]$Q$data

ofst.only.s1<-list()
ofst.only.s2<-list()
i=1
#Loop through all covariates and develop a quadscheme that avoids areas where these covariates are NA 
for (i in 1:11){
  seed<-2003+(2*i-1)
  set.seed(seed)
  covs.s1<-list(lc.dev.HI_1km_smth[[i]],
                lc.dev.MI_1km_smth[[i]],
                lc.dev.LI_1km_smth[[i]],
                lc.dev.open_1km_smth[[i]],
                lch.resc_1km_smth[[i]],
                NDVI.ws.rsc.stn_1km[[i+1]],
                NDVI.sf.rsc.stn_1km[[i]],
                hsi.resc.stn_1km,
                blr.mpd.resc.ofst.stan[[i]]
  )
  
  
  Q.s1<-quadscheme(vfdf[[2*i-1]],nd=250)
  Q.s1<-trim_quad_na_any(Q.s1,covs.s1)
  
  ofst.only.s1[[i]]<-ppm(Q.s1~marks
                         +offset(log(sm.mpd.ofst.stan))
                         
                         ,
                         covariates=list(lnd.cov.HI=lc.dev.HI_1km_smth[[i]],
                                         lnd.cov.MI=lc.dev.MI_1km_smth[[i]],
                                         lnd.cov.LI=lc.dev.LI_1km_smth[[i]],
                                         lnd.cov.open=lc.dev.open_1km_smth[[i]],
                                         lnd.ch=lch.resc_1km_smth[[i]],
                                         NDVI_cur.stan=NDVI.ws.rsc.stn_1km[[i+1]],
                                         NDVI_lag1=NDVI.sf.rsc.stn_1km[[i]],
                                         NDVI_lag2=NDVI.ws.rsc.stn_1km[[i]],
                                         hsi.stan=hsi.resc.stn_1km,
                                         sm.mpd.ofst.stan= blr.mpd.resc.ofst.stan[[i]])
                         
                         , improve.type="ho", nsim=100)
  attr(ofst.only.s1[[i]],"seed")<-seed
  
  seed<-2003+(2*i)
  set.seed(seed)
  
  covs.s2<-list(lc.dev.HI_1km_smth[[i]],
                lc.dev.MI_1km_smth[[i]],
                lc.dev.LI_1km_smth[[i]],
                lc.dev.open_1km_smth[[i]],
                lch.resc_1km_smth[[i]],
                #NDVI.sf.rsc.stn_1km[[i]],
                NDVI.sf.rsc.stn_1km[[i+1]],
                NDVI.ws.rsc.stn_1km[[i+1]],
                hsi.resc.stn_1km,
                blr.mpd.resc.ofst.stan[[i]])
  
  
  Q.s2<-quadscheme(vfdf[[2*i]],nd=250)
  Q.s2<-trim_quad_na_any(Q.s2,covs.s2)
  
  ofst.only.s2[[i]]<-ppm(Q.s2~marks
                         +offset(log(sm.mpd.ofst.stan))
                         ,
                         covariates=list(lnd.cov.HI=lc.dev.HI_1km_smth[[i]],
                                         lnd.cov.MI=lc.dev.MI_1km_smth[[i]],
                                         lnd.cov.LI=lc.dev.LI_1km_smth[[i]],
                                         lnd.cov.open=lc.dev.open_1km_smth[[i]],
                                         lnd.ch=lch.resc_1km_smth[[i]],
                                         NDVI_cur.stan=NDVI.sf.rsc.stn_1km[[i+1]],
                                         NDVI_lag1=NDVI.ws.rsc.stn_1km[[i+1]],
                                         NDVI_lag2=NDVI.sf.rsc.stn_1km[[i]],
                                         
                                         hsi.stan=hsi.resc.stn_1km,
                                         sm.mpd.ofst.stan= blr.mpd.resc.ofst.stan[[i]])
                         
                         , improve.type="ho", nsim=100)
  attr(ofst.only.s2[[i]], "seed")<-seed
}  

# Tabulate and View offset-only model results
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
View(cm.ofst.only)


################################## Fitting linear models


linear_model_s1<-list()
linear_model_s2<-list()

summary(linear_model_s1[[i]])
for (i in 1:11){
  seed<-4003+(2*i-1) #add 200 to make the second run seed different
  set.seed(seed)
  covs.s1<-list(lc.dev.HI_1km_smth[[i]],
                lc.dev.MI_1km_smth[[i]],
                lc.dev.LI_1km_smth[[i]],
                lc.dev.open_1km_smth[[i]],
                lch.resc_1km_smth[[i]],
                NDVI.ws.rsc.stn_1km[[i+1]],
                NDVI.sf.rsc.stn_1km[[i]],
                hsi.resc.stn_1km,
                blr.mpd.resc.ofst.stan[[i]])
  
  Q.s1<-quadscheme(vfdf[[2*i-1]],nd=250)
  Q.s1<-trim_quad_na_any(Q.s1,covs.s1)
  
  linear_model_s1[[i]]<-ppm(Q.s1~marks
                                     +offset(log(sm.mpd.ofst.stan))
                                     +lnd.cov.HI+lnd.cov.MI+lnd.cov.LI+lnd.cov.open+
                                       +lnd.ch
                                     +NDVI_cur.stan
                                     +I(NDVI_cur.stan-NDVI_lag1)
                                     +hsi.stan,
                                     
                                     
                                     
                                     covariates=list(
                                       lnd.cov.HI= lc.dev.HI_1km_smth[[i]],
                                       lnd.cov.MI=lc.dev.MI_1km_smth[[i]],
                                       lnd.cov.LI=lc.dev.LI_1km_smth[[i]],
                                       lnd.cov.open=lc.dev.open_1km_smth[[i]],
                                       lnd.ch= lch.resc_1km_smth[[i]],
                                       NDVI_cur.stan= NDVI.ws.rsc.stn_1km[[i+1]],
                                       NDVI_lag1= NDVI.sf.rsc.stn_1km[[i]],
                                       hsi.stan= hsi.resc.stn_1km,
                                       sm.mpd.ofst.stan=  blr.mpd.resc.ofst.stan[[i]])
                                     
                                     
                                     ,improve.type="ho", nsim=100)
  attr(linear_model_s1[[i]],"seed")<-seed
  
  seed<-4003+(2*i)
  set.seed(seed)
  
  covs.s2<-list(lc.dev.HI_1km_smth[[i]],
                lc.dev.MI_1km_smth[[i]],
                lc.dev.LI_1km_smth[[i]],
                lc.dev.open_1km_smth[[i]],
                lch.resc_1km_smth[[i]],
                #NDVI.sf.rsc.stn_1km[[i]],
                NDVI.sf.rsc.stn_1km[[i+1]],
                NDVI.ws.rsc.stn_1km[[i+1]],
                hsi.resc.stn_1km,
                blr.mpd.resc.ofst.stan[[i]])
  
  Q.s2<-quadscheme(vfdf[[2*i]],nd=250)
  Q.s2<-trim_quad_na_any(Q.s2,covs.s2)
  
  linear_model_s2[[i]]<-ppm(Q.s2~marks
                                     +offset(log(sm.mpd.ofst.stan))
                                     +lnd.cov.HI+lnd.cov.MI+lnd.cov.LI+lnd.cov.open+
                                       +lnd.ch
                                     +NDVI_cur.stan
                                     +I(NDVI_cur.stan-NDVI_lag1)
                                     +hsi.stan,,
                                     
                                     covariates=list(lnd.cov.HI=lc.dev.HI_1km_smth[[i]],
                                                     lnd.cov.MI=lc.dev.MI_1km_smth[[i]],
                                                     lnd.cov.LI=lc.dev.LI_1km_smth[[i]],
                                                     lnd.cov.open=lc.dev.open_1km_smth[[i]],
                                                     lnd.ch=lch.resc_1km_smth[[i]],
                                                     NDVI_cur.stan=NDVI.sf.rsc.stn_1km[[i+1]],
                                                     NDVI_lag1=NDVI.ws.rsc.stn_1km[[i+1]],
                                                     hsi.stan=hsi.resc.stn_1km,
                                                     sm.mpd.ofst.stan= blr.mpd.resc.ofst.stan[[i]])
                                     
                                     ,improve.type="ho", nsim=100)
  attr(linear_model_s2[[i]],"seed")<-seed
  
}

# Tabulating and viewing results for linear model
cm_linear_mod<-matrix(nrow= 22, 
                       ncol = 2*length(names(linear_model_s1[[1]]$coef))+3)

years.vf<-rep(seq(2013,2023,1), each=2)
sem.rep<-rep(c(1,2),11)
mat.names<-paste0(years.vf," ",
                  sem.rep)
rownames(cm_linear_mod)<-mat.names
colnames(cm_linear_mod)<-c(rep(names(linear_model_s1[[1]]$coef),each=2),"AIC","M_LogLik", "Converged")


for (i in 1:(nrow(cm_linear_mod)/2)){
  for (j in 1:((ncol(cm_linear_mod)-3)/2)){
    cm_linear_mod[2*i-1,2*j-1]<-coef(linear_model_s1[[i]])[[j]]
    cm_linear_mod[2*i,2*j-1]<-coef(linear_model_s2[[i]])[[j]]
    
    smry1<-summary(linear_model_s1[[i]])
    smry2<-summary(linear_model_s2[[i]])
    cm_linear_mod[2*i-1,2*j]<-as.character(smry1$coefs.SE.CI$Ztest[[j]])
    cm_linear_mod[2*i,2*j]<-as.character(smry2$coefs.SE.CI$Ztest[[j]])
    
    
  }
}
for (i in 1:(nrow(cm_linear_mod)/2)){
  cm_linear_mod[2*i-1,"AIC"]<-
    AIC.ppm(linear_model_s1[[i]])
  cm_linear_mod[2*i,"AIC"]<-
    AIC.ppm(linear_model_s2[[i]])
}
for (i in 1:(nrow(cm_linear_mod)/2)){
  cm_linear_mod[2*i-1,"M_LogLik"]<-
    logLik.ppm(linear_model_s1[[i]])
  cm_linear_mod[2*i,"M_LogLik"]<-
    logLik.ppm(linear_model_s2[[i]])
}

for (i in 1:11){
  summs1<-summary(linear_model_s1[[i]])
  cm_linear_mod[2*i-1,"Converged"]<-summs1$converged
  
  summs2<-summary(linear_model_s2[[i]])
  cm_linear_mod[2*i,"Converged"]<- summs2$converged
}
View(cm_linear_mod)  


predicted.ints.linear.mod<-list()
for (i in 1:length(linear_model_s1)){
  predicted.ints.linear.mod[[2*i-1]]<-predict.ppm(linear_model_s1[[i]])
  predicted.ints.linear.mod[[2*i]]<-predict.ppm(linear_model_s2[[i]])
  
}

############################## Fitting full or quadratic models 
quadratic_model_s1<-list()
quadratic_model_s2<-list()
#Loop through all covariates and develop a quadscheme that avoids areas where these covariates are NA 
for (i in 1:11){
  seed<-5003+(2*i-1) # making seeds that save for each iteration
  set.seed(seed)
  covs.s1<-list(lc.dev.HI_1km_smth[[i]],
                lc.dev.MI_1km_smth[[i]],
                lc.dev.LI_1km_smth[[i]],
                lc.dev.open_1km_smth[[i]],
                lch.resc_1km_smth[[i]],
                #NDVI.ws.rsc.stn[[i]],
                NDVI.ws.rsc.stn_1km[[i+1]],
                NDVI.sf.rsc.stn_1km[[i]],
                hsi.resc.stn_1km,
                blr.mpd.resc.ofst.stan[[i]]
  )
  
  Q.s1<-quadscheme(vfdf[[2*i-1]],nd=250)
  Q.s1<-trim_quad_na_any(Q.s1,covs.s1)
  
  quadratic_model_s1[[i]]<-ppm(Q.s1~marks
                                +offset(log(sm.mpd.ofst.stan))
                                +lnd.cov.HI+lnd.cov.MI+lnd.cov.LI+lnd.cov.open+
                                  +lnd.ch
                                +NDVI_cur.stan+I(NDVI_cur.stan^2)
                                +I(NDVI_cur.stan-NDVI_lag1)+I((NDVI_cur.stan-NDVI_lag1)^2)
                                +hsi.stan+I(hsi.stan^2),
                                
                                
                                
                                covariates=list(
                                  lnd.cov.HI= lc.dev.HI_1km_smth[[i]],
                                  lnd.cov.MI=lc.dev.MI_1km_smth[[i]],
                                  lnd.cov.LI=lc.dev.LI_1km_smth[[i]],
                                  lnd.cov.open=lc.dev.open_1km_smth[[i]],
                                  lnd.ch= lch.resc_1km_smth[[i]],
                                  NDVI_cur.stan= NDVI.ws.rsc.stn_1km[[i+1]],
                                  NDVI_lag1= NDVI.sf.rsc.stn_1km[[i]],
                                  hsi.stan= hsi.resc.stn_1km,
                                  sm.mpd.ofst.stan=  blr.mpd.resc.ofst.stan[[i]]),
                                
                                
                                ,improve.type="ho", nsim=100)
  attr(quadratic_model_s1[[i]],"seed")<-seed
  
  seed<-5003+(2*i) 
  set.seed(seed)
  covs.s2<-list(lc.dev.HI_1km_smth[[i]],
                lc.dev.MI_1km_smth[[i]],
                lc.dev.LI_1km_smth[[i]],
                lc.dev.open_1km_smth[[i]],
                lch.resc_1km_smth[[i]],
                NDVI.sf.rsc.stn_1km[[i+1]],
                NDVI.ws.rsc.stn_1km[[i+1]],
                hsi.resc.stn_1km,
                blr.mpd.resc.ofst.stan[[i]]
  )
  
  Q.s2<-quadscheme(vfdf[[2*i]],nd=250)
  Q.s2<-trim_quad_na_any(Q.s2,covs.s2)
  quadratic_model_s2[[i]]<-ppm(Q.s2~marks
                                +offset(log(sm.mpd.ofst.stan))
                                +lnd.cov.HI+lnd.cov.MI+lnd.cov.LI+lnd.cov.open+
                                  +lnd.ch
                                +NDVI_cur.stan+I(NDVI_cur.stan^2)
                                +I(NDVI_cur.stan-NDVI_lag1)+I((NDVI_cur.stan-NDVI_lag1)^2)
                                +hsi.stan+I(hsi.stan^2),
                                
                                
                                covariates=list(
                                  lnd.cov.HI= lc.dev.HI_1km_smth[[i]],
                                  lnd.cov.MI=lc.dev.MI_1km_smth[[i]],
                                  lnd.cov.LI=lc.dev.LI_1km_smth[[i]],
                                  lnd.cov.open=lc.dev.open_1km_smth[[i]],
                                  lnd.ch= lch.resc_1km_smth[[i]],
                                  NDVI_cur.stan= NDVI.sf.rsc.stn_1km[[i+1]],
                                  NDVI_lag1= NDVI.ws.rsc.stn_1km[[i+1]],
                                  hsi.stan= hsi.resc.stn_1km,
                                  sm.mpd.ofst.stan=  blr.mpd.resc.ofst.stan[[i]]),
                                
                                
                                ,improve.type="ho", nsim=100)
  attr(quadratic_model_s2[[i]],"seed")<-seed
}

# Tabulate and view quadratic model results
cm_quadratic_mod<-matrix(nrow= 22, 
                            ncol = 2*length(names(quadratic_model_s1[[1]]$coef))+3)

years.vf<-rep(seq(2013,2023,1), each=2)
sem.rep<-rep(c(1,2),11)
mat.names<-paste0(years.vf," ",
                  sem.rep)
rownames(cm_quadratic_mod)<-mat.names
colnames(cm_quadratic_mod)<-c(rep(names(quadratic_model_s1[[1]]$coef),each=2),"AIC", "Converged", "M_LogLik")


for (i in 1:(nrow(cm_quadratic_mod)/2)){
  for (j in 1:((ncol(cm_quadratic_mod)-3)/2)){
    cm_quadratic_mod[2*i-1,2*j-1]<-coef(quadratic_model_s1[[i]])[[j]]
    cm_quadratic_mod[2*i,2*j-1]<-coef(quadratic_model_s2[[i]])[[j]]
    
    smry1<-summary(quadratic_model_s1[[i]])
    smry2<-summary(quadratic_model_s2[[i]])
    cm_quadratic_mod[2*i-1,2*j]<-as.character(smry1$coefs.SE.CI$Ztest[[j]])
    cm_quadratic_mod[2*i,2*j]<-as.character(smry2$coefs.SE.CI$Ztest[[j]])
    
  }
} 
for (i in 1:(nrow(cm_quadratic_mod)/2)){
  cm_quadratic_mod[2*i-1,"AIC"]<-
    AIC.ppm(quadratic_model_s1[[i]])
  cm_quadratic_mod[2*i,"AIC"]<-
    AIC.ppm(quadratic_model_s2[[i]])
}

for (i in 1:11){
  summs1<-summary(quadratic_model_s1[[i]])
  cm_quadratic_mod[2*i-1,"Converged"]<-summs1$converged
  
  summs2<-summary(quadratic_model_s2[[i]])
  cm_quadratic_mod[2*i,"Converged"]<- summs2$converged
}

for (i in 1:(nrow(cm_quadratic_mod)/2)){
  cm_quadratic_mod[2*i-1,"M_LogLik"]<-
    logLik.ppm(quadratic_model_s1[[i]])
  cm_quadratic_mod[2*i,"M_LogLik"]<-
    logLik.ppm(quadratic_model_s2[[i]])
}

View(cm_quadratic_mod)

# Predicted intensities for quadratic model
predicted.ints.quad.mod<-list()
for (i in 1:length(quadratic_model_s1)){
  predicted.ints.quad.mod[[2*i-1]]<-predict.ppm(quadratic_model_s1[[i]])
  predicted.ints.quad.mod[[2*i]]<-predict.ppm(quadratic_model_s2[[i]])
  
}

## Relative Intensity plots
rel.den.qd.mod_ep_0.001<-list()
rel.den.ofst.mod_ep_0.001<-list()
rel.den.lin.mod_ep_0.001<-list()
pred.quad.mod<-list()
pred.ofst.mod<-list()
pred.lin.mod<-list()


vals<-seq(0,20, length.out=256)
white_pos<-which.min(abs(vals-1.5))
n1<-white_pos
cols1<-colorRampPalette(c("blue","grey"))(n1)
n2<-256-n1
cols2<-colorRampPalette(c("grey","darkred"))(n2)
my_colors<-c(cols1,cols2)
epsilon=0.001

ho=pilot.bw=15 # used for rel int plos

for (i in 1:11){
  ##################################
  #### quadratic full model Semester 1
  ### Obtain prediction 'window' that will eliminate areas where covariates in this semester/year are NA
  M3_quad_model_s1<-quadratic_model_s1[[i]]
  
  preds<-predict(M3_quad_model_s1, dimyx=250)
  pred.quad.mod[[2*i-1]]<-preds[[1]]+preds[[2]]+preds[[3]]+preds[[4]]
  Window(M3_quad_model_s1$Q$data)<-as.owin(preds[[1]])
  
  pred.int<-predict.ppm(M3_quad_model_s1, locations=M3_quad_model_s1$Q$data)#locations=subset.ppp(test.pp,as.owin(test.mpd.im)))
  pred.int<-ifelse(pred.int<epsilon,min(pred.int[pred.int>epsilon|pred.int==epsilon], na.rm=TRUE),pred.int)
  rel.den.qd.mod_ep_0.001[[2*i-1]]<-adaptive.density(M3_quad_model_s1$Q$data, method="kernel", weights=(1/pred.int), ho=pilot.bw, dimyx=250, edge=TRUE)
  
  adap.den<- rel.den.qd.mod_ep_0.001[[2*i-1]]

  adap.den$v<-ifelse((adap.den$v<2)&(adap.den$v>0.5),NA,adap.den$v)
  plot(adap.den, col=alpha(my_colors,1), main=paste0("Quad. Model ", 2012+i," sem. 1"),zlim=c(0,20))
  #plot(phx_osm, add=TRUE)
  plot(adap.den, col=alpha(my_colors,0.7), add=TRUE,zlim=c(0,20))
    # 
  ################################# Linear model sem 1
  M2_lin_model_s1<-linear_model_s1[[i]]
  
  preds<-predict(M2_lin_model_s1, dimyx=250)
  pred.lin.mod[[2*i-1]]<-preds[[1]]+preds[[2]]+preds[[3]]+preds[[4]]
  Window(M2_lin_model_s1$Q$data)<-as.owin(preds[[1]])
  
  pred.int<-predict.ppm(M2_lin_model_s1, locations=M2_lin_model_s1$Q$data)#locations=subset.ppp(test.pp,as.owin(test.mpd.im)))
  pred.int<-ifelse(pred.int<epsilon,min(pred.int[pred.int>epsilon|pred.int==epsilon], na.rm=TRUE),pred.int)
  rel.den.lin.mod_ep_0.001[[2*i-1]]<-adaptive.density(M2_lin_model_s1$Q$data, method="kernel", weights=1/pred.int, ho=pilot.bw, dimyx=250, edge=TRUE)

  adap.den<- rel.den.lin.mod_ep_0.001[[2*i-1]]

  adap.den$v<-ifelse((adap.den$v<2)&(adap.den$v>0.5),NA,adap.den$v)
  plot(adap.den, col=alpha(my_colors,1), main=paste0("Linear Model ", 2012+i," sem. 1"),zlim=c(0,20))
  plot(phx_osm, add=TRUE)
  plot(adap.den, col=alpha(my_colors,0.7), add=TRUE,zlim=c(0,20))
  
  ################################## ofst.only.s1
  #### Offset Only Model semester 1
  
  M1_ofst_s1<-ofst.only.s1[[i]]
  Window(M1_ofst_s1$Q$data)<-as.owin(pred.quad.mod[[2*i-1]]) # note the change of the window to that of the quadratic model's predicted intensity: this is to make sure that areas
  																#removed due to NA covariate values do not get assigned relative intensity values at these locations and the plots
  																# for all models match (this could have been from either the linear or quadratic model)
 
   preds<-predict(M1_ofst_s1, dimyx=250)
  pred.ofst.mod[[2*i-1]]<-preds[[1]]+preds[[2]]+preds[[3]]+preds[[4]]
  pred.int<-predict.ppm(M1_ofst_s1, locations=M1_ofst_s1$Q$data)
  pred.int<-ifelse(pred.int<epsilon,min(pred.int[pred.int>epsilon|pred.int==epsilon], na.rm=TRUE),pred.int)

  rel.den.ofst.mod_ep_0.001[[2*i-1]]<-adaptive.density(M1_ofst_s1$Q$data, method="kernel", weights=1/pred.int, ho=pilot.bw, dimyx=250, edge=TRUE)

  adap.den<- rel.den.ofst.mod_ep_0.001[[2*i-1]]

  adap.den$v<-ifelse((adap.den$v<2)&(adap.den$v>0.5),NA,adap.den$v)
  plot(adap.den, col=alpha(my_colors,1), main=paste0("Offset Model ", 2012+i," sem. 1"),zlim=c(0,20))
  plot(phx_osm, add=TRUE)
  plot(adap.den, col=alpha(my_colors,0.7), add=TRUE,zlim=c(0,20))
 ###########################
  ###########################
  ###########################
  ##########################
  ##################################
  plot.new()
  #### quadratic full model Semester 2
  M3_quad_model_s2<-quadratic_model_s2[[i]]
  
  preds<-predict(M3_quad_model_s2, dimyx=250)
  pred.quad.mod[[2*i]]<-preds[[1]]+preds[[2]]+preds[[3]]+preds[[4]]
  Window(M3_quad_model_s2$Q$data)<-as.owin(preds[[1]])
  
  pred.int<-predict.ppm(M3_quad_model_s2, locations=M3_quad_model_s2$Q$data)#locations=subset.ppp(test.pp,as.owin(test.mpd.im)))
  pred.int<-ifelse(pred.int<epsilon,min(pred.int[pred.int>epsilon|pred.int==epsilon], na.rm=TRUE),pred.int)
  rel.den.qd.mod_ep_0.001[[2*i]]<-adaptive.density(M3_quad_model_s2$Q$data, method="kernel", weights=1/pred.int, ho=pilot.bw, dimyx=250, edge=TRUE)
  adap.den<- rel.den.qd.mod_ep_0.001[[2*i]]

  adap.den$v<-ifelse((adap.den$v<2)&(adap.den$v>0.5),NA,adap.den$v)
  plot(adap.den, col=alpha(my_colors,1), main=paste0("Quad. Model ", 2012+i," sem. 2"),zlim=c(0,20))
  plot(phx_osm, add=TRUE)
  plot(adap.den, col=alpha(my_colors,0.7), add=TRUE,zlim=c(0,20))
  
  ################################# Linear model sem 2
  M2_lin_model_s2<-linear_model_s2[[i]]
  
  preds<-predict(M2_lin_model_s2, dimyx=250)
  pred.lin.mod[[2*i]]<-preds[[1]]+preds[[2]]+preds[[3]]+preds[[4]]
  Window(M2_lin_model_s2$Q$data)<-as.owin(preds[[1]])
  
  pred.int<-predict.ppm(M2_lin_model_s2, locations=M2_lin_model_s2$Q$data)#locations=subset.ppp(test.pp,as.owin(test.mpd.im)))
  pred.int<-ifelse(pred.int<epsilon,min(pred.int[pred.int>epsilon|pred.int==epsilon], na.rm=TRUE),pred.int)
  rel.den.lin.mod_ep_0.001[[2*i]]<-adaptive.density(M2_lin_model_s2$Q$data, method="kernel", weights=1/pred.int, ho=pilot.bw, dimyx=250, edge=TRUE)
  
  adap.den<- rel.den.lin.mod_ep_0.001[[2*i]]
  
  adap.den$v<-ifelse((adap.den$v<2)&(adap.den$v>0.5),NA,adap.den$v)
  plot(adap.den, col=alpha(my_colors,1), main=paste0("Linear Model ", 2012+i," sem. 2"),zlim=c(0,20))
  plot(phx_osm, add=TRUE)
  plot(adap.den, col=alpha(my_colors,0.7), add=TRUE,zlim=c(0,20))   
 
  #### Offset Only Model semester 2
  M1_ofst_s2<-ofst.only.s2[[i]]
  
  Window(M1_ofst_s2$Q$data)<-as.owin(pred.quad.mod[[2*i]])
  preds<-predict(M1_ofst_s2, dimyx=250)
  pred.ofst.mod[[2*i]]<-preds[[1]]+preds[[2]]+preds[[3]]+preds[[4]]
  
  pred.int<-predict.ppm(M1_ofst_s2, locations=M1_ofst_s2$Q$data)
  pred.int<-ifelse(pred.int<epsilon,min(pred.int[pred.int>epsilon|pred.int==epsilon], na.rm=TRUE),pred.int)
  rel.den.ofst.mod_ep_0.001[[2*i]]<-adaptive.density(M1_ofst_s2$Q$data, method="kernel", weights=1/pred.int, ho=pilot.bw, dimyx=250,edge=TRUE)
  adap.den<- rel.den.ofst.mod_ep_0.001[[2*i]]

  adap.den$v<-ifelse((adap.den$v<2)&(adap.den$v>0.5),NA,adap.den$v)
  plot(adap.den, col=alpha(my_colors,1), main=paste0("Offset Model ", 2012+i," sem. 2"),zlim=c(0,20))
  plot(phx_osm, add=TRUE)
  plot(adap.den, col=alpha(my_colors,0.7), add=TRUE,zlim=c(0,20))

}



