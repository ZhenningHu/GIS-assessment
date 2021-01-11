# Load Packages
library(dplyr)
library(ggplot2)
library(leaflet)
library(rgdal)
library(raster)
library(deldir)
library(mgcv)
library(rgeos)
library(mapproj)
library(pacman) #p_load(mapproj)
library(sf)
library(spatstat)
library(here)
library(sp)
library(rgeos)
library(maptools)
library(GISTools)
library(tmap)
library(geojson)
library(geojsonio)
library(tmaptools)
library(tidyverse)
library(stringr)
library(spdep)

# Data Loading
# Set working directory
getwd()
setwd("C:/Users/Arios/Documents/GIS/input")

# Starbucks Data 
starbucks <- read.csv("../input/Starbucks_london.csv", header = TRUE, sep = ",")

# Shapefiles with london borough
ldn<- st_read(here::here("GIS",
                         "input",
                         "london.shp"))%>%
  st_set_crs(., 27700) %>% st_transform(.,27700)

# UK pubs
pub <- read.csv("../input/open_pubs.csv", header = TRUE, sep = ",")

# Data Cleaning
# Filter by local authority
AuthoritiesLondon = c("Barking and Dagenham","Barnet","Bexley","Brent","Bromley",
                      "Camden","City of London","Croydon","Ealing","Enfield","Greenwich",
                      "Hackney","Hammersmith and Fulham","Haringey","Harrow","Havering",
                      "Hillingdon","Hounslow","Islington","Kensington and Chelsea","Kingston upon Thames",
                      "Lambeth","Lewisham","Merton","Newham","Redbridge",
                      "Richmond upon Thames","Southwark","Sutton","Tower Hamlets","Waltham Forest",
                      "Wandsworth","Westminster")
pub %>% filter(local_authority %in% AuthoritiesLondon) -> publ
publ %>% group_by(local_authority) %>% summarize(pub = length(name)) -> publondon
# Rename column
names(starbucks)[names(starbucks) == 'local_authority'] <- 'NAME'
names(ldn)[names(ldn) == 'name'] <- 'borough'
names(publondon)[names(publondon) == 'local_authority'] <- 'NAME'

#Join the he summarized data into the London borough
ldn <- ldn %>%
  left_join(.,
            publondon,
            by = "NAME")%>%
  left_join(.,
            starbucks,
            by = "NAME")


# Data Visualization
# Starbucks in boroughs of London
tm1 <- tm_shape(ldn) +
  tm_polygons("starbucks",
              palette="Reds",
              midpoint=NA,
              title="Starbucks")+
  tm_layout(frame=FALSE,
            legend.position = c("right","bottom"), 
            legend.text.size=0.6, 
            legend.title.size = 1)+
  tm_scale_bar(position=c(0,0), text.size=0.4)+
  tm_compass(type = "8star",size=3,position=c(0,0.15))+
  tm_credits("(c) Greater London Authority",position=c(0,0),size = 0.5)
tm1

# Pubs in boroughs of London
tm2 <- tm_shape(ldn) +
  tm_polygons("pub",
              style="jenks",
              palette=brewer.pal(5, "Reds"),
              midpoint=NA,
              title="Pub")+
  tm_layout(frame=FALSE,
            legend.position= c("right","bottom"),
            legend.text.size=0.6,
            legend.title.size = 1)+
  tm_scale_bar(position=c(0,0), text.size=0.4)+
  tm_compass(type = "8star",size=3,position=c(0,0.15))+
  tm_credits("(c) Greater London Authority",position=c(0,0),size = 0.5)
tm2

# Spatial Auto-correlation Analysis
# calculate the centroids of all Wards in London
coordsW <- ldn %>%
  st_centroid()%>%
  st_geometry()
#plot(coordsW,axes=TRUE)

#generate a spatial weight matrix of nearest k-nearest neighbours (k=4)
knn_wards <-coordsW %>%
  knearneigh(., k=4) 
ldn_knn <- knn_wards %>%
  knn2nb() %>%
  nb2listw(., style="C")

# Global Moran's I Test
# For Starbucks’ density
I_global_starbucks <- ldn %>%
  pull(starbucks) %>%
  as.vector()%>%
  moran.test(., ldn_knn)
I_global_starbucks

# For pubs’ density
I_global_pub <- ldn %>%
  pull(pub) %>%
  as.vector()%>%
  moran.test(., ldn_knn)
I_global_pub

# Global Geary's C Test
# For Starbucks’ density
C_global_starbucks <-ldn %>%
  pull(starbucks) %>%
  as.vector()%>%
  geary.test(., ldn_knn)
C_global_starbucks

# For pubs’ density
C_global_pub <-ldn %>%
  pull(pub) %>%
  as.vector()%>%
  geary.test(., ldn_knn)
C_global_pub

# Global Gentis' G Test
# For Starbucks’ density
G_global_starbucks <-
  ldn %>%
  pull(starbucks) %>%
  as.vector() %>%
  globalG.test(., ldn_knn)
G_global_starbucks

# For pubs’ density
G_global_pub <-
  ldn %>%
  pull(pub) %>%
  as.vector()%>%
  globalG.test(., ldn_knn)
G_global_pub

# Local Getis' G z-score
# calculate local versions of the Getis' G statistic 
G_local_starbucks <- ldn %>%
  pull(starbucks) %>%
  as.vector()%>%
  localG(., ldn_knn)

G_local_pub <- ldn %>%
  pull(pub) %>%
  as.vector()%>%
  localG(., ldn_knn)

# Covert into a dataframe and append into borough data frame
G_local_pub_df <- data.frame(matrix(unlist(G_local_pub), nrow=33, byrow=T))
ldn_pub <- ldn %>%
  mutate(z_score = as.numeric(G_local_pub_df$matrix.unlist.G_local_pub...nrow...33..byrow...T.))

G_local_starbucks_df <- data.frame(matrix(unlist(G_local_starbucks), nrow=33, byrow=T))
ldn_starbucks <- ldn %>%
  mutate(z_score = as.numeric(G_local_starbucks_df$matrix.unlist.G_local_starbucks...nrow...33..byrow...T.))

# Set the breaks and color bar manually
breaks1<-c(-1000.00,-2.58,-1.96,-1.65,1.65,1.96,2.58,1000.00)

# create a new diverging colour brewer palette
GIColours<- rev(brewer.pal(8, "RdBu"))
tmap_mode("plot")

# Starbucks’ Local Getis-Ord's G z-score map
GI_starbucks<-tm_shape(ldn_starbucks) +
  tm_polygons("z_score",
              style="fixed",
              breaks=breaks1,
              palette=GIColours,
              midpoint=NA,
              title="Gi* Starbucks")+
  tm_layout(frame=FALSE,
            legend.position = c("right","bottom"),
            legend.text.size=0.6,
            legend.title.size = 1)+
  tm_scale_bar(position=c(0,0), text.size=0.4)+
  tm_compass(type = "8star",size=3,position=c(0,0.15))+
  tm_credits("(c) Greater London Authority",position=c(0,0),size = 0.5)
GI_starbucks

# Pub’s Local Getis-Ord's G z-score map
GI_pub<-tm_shape(ldn_pub) +
  tm_polygons("z_score",
              style="fixed",
              breaks=breaks1,
              palette=GIColours,
              midpoint=NA,
              title="Gi* Pub")+
  tm_layout(frame=FALSE,
            legend.position = c("right","bottom"),
            legend.text.size=0.6,
            legend.title.size = 1)+
  tm_scale_bar(position=c(0,0), text.size=0.4)+
  tm_compass(type = "8star",size=3,position=c(0,0.15))+
  tm_credits("(c) Greater London Authority",position=c(0,0),size = 0.5)
GI_pub

