# Mesh exploration
# Off Aqua
# Tim Szewczyk

# This script is just to get a handle on the high resolution mesh


library(tidyverse); library(ncdf4); library(sf); library(WeStCOMS)


hiRes.dir <- "D:/hydroOut/linnhe7/"
hiRes_1h.dir <- paste0(hiRes.dir, "linnhe7_tides_met_tsobc/netcdf_2021/")
hiRes_5min.dir <- paste0(hiRes.dir, "linnhe7_tides_met_tsobc_mf_5min/")




load_WeStCOMS <- function(type="node", nc_file, make_sf=FALSE,
                          n_hours=24, n_sigLay=10, n_sigLev=11) {
  library(ncdf4)
  ncin <- nc_open(nc_file)
  
  isNode <- type=="node"
  lon <- c(ncvar_get(ncin, ifelse(isNode, "lon", "lonc")))
  lat <- c(ncvar_get(ncin, ifelse(isNode, "lat", "latc")))
  n_elements <- length(lon)
  sigLayVals <- ncvar_get(ncin, "siglay")[1,]
  n_sigLay <- length(sigLayVals)
  sigLevVals <- ncvar_get(ncin, "siglev")[1,]
  n_sigLev <- length(sigLevVals)
  
  # load mesh data
  mesh <- tibble(i=1:n_elements, lon, lat,
                 h=c(ncvar_get(ncin, ifelse(isNode, "h", "h_center"))))
  if(isNode) {
    mesh <- mesh %>% add_column(area=c(ncvar_get(ncin,"art1")))
  }
  
  
  # load hourly data: [elem/node][hour]
  hour <- tibble(i=rep(1:n_elements, n_hours), 
                 hour=rep(0:(n_hours-1), each=n_elements)) %>%
    full_join(., mesh, by="i")
  if(isNode) {
    hour <- hour %>% 
      add_column(zeta=c(ncvar_get(ncin,"zeta")),
                 short_wave=c(ncvar_get(ncin,"short_wave")),
                 net_heat_flux=c(ncvar_get(ncin,"net_heat_flux")),
                 precip=c(ncvar_get(ncin,"precip")),
                 evap=c(ncvar_get(ncin,"evap")))
  } else {
    hour <- hour %>%
      add_column(uwind_speed=c(ncvar_get(ncin,"uwind_speed")),
                 vwind_speed=c(ncvar_get(ncin,"vwind_speed")),
                 tauc=c(ncvar_get(ncin,"tauc")),
                 ua=c(ncvar_get(ncin,"ua")),
                 va=c(ncvar_get(ncin,"va")),
                 velocity2D=sqrt(ua^2 + va^2),
                 direction2D=atan2(va, ua))
  }
  
  
  # load sigma level data [elem/node][sigma][hour]
  sigmaLevel <- tibble(i=rep(1:n_elements, n_hours*n_sigLev), 
                       hour=rep(rep(0:(n_hours-1), each=n_elements*n_sigLev)),
                       sigLevel=rep(rep(1:n_sigLev, each=n_elements), n_hours)) %>%
    full_join(., mesh, by="i") %>%
    mutate(levelDepth=-h*(sigLevVals[sigLevel]))
  if(isNode) {
    sigmaLevel <- sigmaLevel %>%
      add_column(omega=c(ncvar_get(ncin, "omega")))
  } else {
    # check variables... not sure anything is measured here
  }
  
  
  # load sigma layer data
  sigmaLayer <- tibble(i=rep(1:n_elements, n_hours*n_sigLay), 
                       hour=rep(rep(0:(n_hours-1), each=n_elements*n_sigLay)),
                       sigLayer=rep(rep(1:n_sigLay, each=n_elements), n_hours)) %>%
    full_join(., mesh, by="i") %>%
    mutate(layerDepth=-h*(sigLayVals[sigLayer]))
  if(isNode) {
    sigmaLayer <- sigmaLayer %>%
      add_column(temp=c(ncvar_get(ncin, "temp")),
                 salinity=c(ncvar_get(ncin, "salinity")))
  } else {
    sigmaLayer <- sigmaLayer %>%
      add_column(u=c(ncvar_get(ncin, "u")),
                 v=c(ncvar_get(ncin, "v")),
                 ww=c(ncvar_get(ncin, "ww")),
                 velocity2D=sqrt(u^2 + v^2),
                 velocity3D=sqrt(u^2 + v^2 + ww^2),
                 direction2D=atan2(v, u))
  }
  
  
  # generate sf if desired
  if(make_sf) {
    library(sf)
    mesh <- mesh %>% st_as_sf(coords=c("lon", "lat"))
    hour <- hour %>% st_as_sf(coords=c("lon", "lat"))
    sigmaLevel <- sigmaLevel %>% st_as_sf(coords=c("lon", "lat"))
    sigmaLayer <- sigmaLayer %>% st_as_sf(coords=c("lon", "lat"))
  }
  
  nc_close(ncin)
  
  return(list(n=nrow(mesh), mesh=mesh, hour=hour, 
              sigLev=sigmaLevel, sigLay=sigmaLayer))
}



day2 <- load_WeStCOMS("node", 
                      dir(hiRes_1h.dir, "0002.nc", full.names=T), 
                      make_sf=T, n_hours=24, n_sigLay=30, n_sigLev=31)
mesh.sf <- st_union(day2$mesh)
mesh.hull <- st_convex_hull(mesh.sf)

gridRes <- 0.005
lice_grid <- st_bbox(mesh.hull) %>%
  st_make_grid(cellsize=c(gridRes, gridRes), what="centers") %>%
  st_sf(id=1:length(.)) %>%
  mutate(inMesh=st_is_within_distance(., mesh.sf, 0.005, sparse=F)[,1]) %>%
  filter(inMesh)



# Hi-res
day2_5m <- nc_open(dir(hiRes_5min.dir, "0002.nc", full.names=T))

ncvar_get(day2_5m, "time")
# These still have 24 slots per file. So the easiest way to deal with this is
# probably (unfortunately) post-hoc processing to adjust the recorded times. 
# This is somewhat annoying since the days and hours are used in the outputs
# and in the filenames, and so it will require some work to get things right.
# I'll also need to adjust the recording intervals to reflect this difference.
times <- numeric()
for(i in seq_along(dir(hiRes_5min.dir))) {
  nc.i <- nc_open(dir(hiRes_5min.dir, full.names=T)[i])
  times <- c(times, ncvar_get(nc.i, "time")*24*60*60)
  nc_close(nc.i)
}
times <- as.POSIXct(times, origin="1858-11-17")
times - lag(times)


hiRes <- (0:255)*337.5
loRes <- (0:23)*3600


