# Mesh exploration
# Off Aqua
# Tim Szewczyk

# This script is just to get a handle on the high resolution mesh


library(tidyverse); library(ncdf4); library(sf); library(WeStCOMS); library(glue)

mesh.fp <- st_read("data/linnhe_mesh_footprint.gpkg")

hiRes.dir <- "D:/hydroOut/linnhe7/"
hiRes_1h.dir <- paste0(hiRes.dir, "linnhe7_tides_met_tsobc_riv/netcdf_2021/")
hiRes_5min.dir <- paste0(hiRes.dir, "linnhe7_tides_met_tsobc_riv_5min/netcdf_2021/")
wc.dir <- "D:/hydroOut/WeStCOMS2/"
wc_1h.dir <- paste0(wc.dir, "Archive/netcdf_2021/")




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
                 net_heat_flux=c(ncvar_get(ncin,"net_heat_flux")))#,
                 # precip=c(ncvar_get(ncin,"precip")),
                 # evap=c(ncvar_get(ncin,"evap")))
  } else {
    hour <- hour %>%
      add_column(uwind_speed=c(ncvar_get(ncin,"uwind_speed")),
                 vwind_speed=c(ncvar_get(ncin,"vwind_speed")),
                 # tauc=c(ncvar_get(ncin,"tauc")),
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
    # sigmaLevel <- sigmaLevel %>%
    #   add_column(omega=c(ncvar_get(ncin, "omega")))
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
                      dir(hiRes_1h.dir, "20211102", full.names=T), 
                      make_sf=T, n_hours=24, n_sigLay=30, n_sigLev=31)
day2.c <- load_WeStCOMS("elem", 
                        dir(hiRes_1h.dir, "20211102", full.names=T), 
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
day2_5m <- nc_open(dir(hiRes_5min.dir, "20211102", full.names=T))

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






load_WeStCOMS("node", 
              dir(hiRes_1h.dir, "20211102", full.names=T), 
              make_sf=T, n_hours=24, n_sigLay=30, n_sigLev=31)$sigLay %>%
  filter(sigLayer==1) %>% 
  st_set_crs(4326) %>%
  st_transform(27700) %>%
  group_by(hour) %>% 
  group_split() %>% 
  imap(~st_write(.x, glue("temp/linnhe_1h_nodeLay1_{str_pad(.y, 2, 'left', 0)}.gpkg")))

load_WeStCOMS("elem", 
              dir(hiRes_1h.dir, "20211102", full.names=T), 
              make_sf=T, n_hours=24, n_sigLay=30, n_sigLev=31)$sigLay %>%
  filter(sigLayer==1) %>%
  st_set_crs(4326) %>%
  st_transform(27700) %>%
  group_by(hour) %>% 
  group_split() %>% 
  imap(~st_write(.x, glue("temp/linnhe_1h_elemLay1_{str_pad(.y, 2, 'left', 0)}.gpkg")))

load_WeStCOMS("node", 
              dir(hiRes_1h.dir, "20211102", full.names=T), 
              make_sf=T, n_hours=24, n_sigLay=30, n_sigLev=31)$hour %>%
  st_set_crs(4326) %>%
  st_transform(27700) %>%
  group_by(hour) %>% 
  group_split() %>% 
  imap(~st_write(.x, glue("temp/linnhe_1h_nodeHour_{str_pad(.y, 2, 'left', 0)}.gpkg")))





load_WeStCOMS("node", 
              dir(wc_1h.dir, "20211102", full.names=T), 
              make_sf=T, n_hours=24, n_sigLay=10, n_sigLev=11)$sigLay %>%
  filter(sigLayer==1) %>% 
  st_set_crs(4326) %>%
  st_transform(27700) %>%
  filter(st_within(., mesh.fp, sparse=F)[,1]) %>%
  group_by(hour) %>% 
  group_split() %>% 
  imap(~st_write(.x, glue("temp/WeStCOMS2_1h_nodeLay1_{str_pad(.y, 2, 'left', 0)}.gpkg")))

load_WeStCOMS("elem", 
              dir(wc_1h.dir, "20211102", full.names=T), 
              make_sf=T, n_hours=24, n_sigLay=10, n_sigLev=11)$sigLay %>%
  filter(sigLayer==1) %>%
  st_set_crs(4326) %>%
  st_transform(27700) %>%
  filter(st_within(., mesh.fp, sparse=F)[,1]) %>%
  group_by(hour) %>% 
  group_split() %>% 
  imap(~st_write(.x, glue("temp/WeStCOMS2_1h_elemLay1_{str_pad(.y, 2, 'left', 0)}.gpkg")))

load_WeStCOMS("node", 
              dir(wc_1h.dir, "20211102", full.names=T), 
              make_sf=T, n_hours=24, n_sigLay=10, n_sigLev=11)$hour %>%
  st_set_crs(4326) %>%
  st_transform(27700) %>%
  filter(st_within(., mesh.fp, sparse=F)[,1]) %>%
  group_by(hour) %>% 
  group_split() %>% 
  imap(~st_write(.x, glue("temp/WeStCOMS2_1h_nodeHour_{str_pad(.y, 2, 'left', 0)}.gpkg")))





library(gganimate)
theme_set(theme_classic())


salinity.df <- map_dfr(dir("temp", "nodeLay1"), 
                       ~st_read(glue("temp/{.x}")) %>%
                         mutate(mesh=str_split_fixed(.x, "_1h_", 2)[,1])) %>%
  mutate(x=st_coordinates(.)[,1], 
         y=st_coordinates(.)[,2]) %>% 
  st_drop_geometry()
interp <- 1
anim <- salinity.df %>%
  ggplot() +
  geom_sf(data=mesh.fp, fill="grey20") +
  geom_point(aes(x, y, colour=salinity, size=area)) + 
  scale_colour_distiller("Surface\nsalinity", type="seq", palette="Blues") + 
  scale_radius(range=c(0.1, 2.5), guide="none") +
  transition_states(hour, wrap=F, transition_length=0, state_length=1) +
  facet_grid(.~mesh) +
  ggtitle(paste0("Hour: ", "{closest_state}"))
anim_save(glue("figs/mesh_surface_salinity.gif"), 
          anim, nframes=interp*24,
          fps=8, width=9, height=6, res=300, units="in")
anim <- salinity.df %>%
  ggplot() +
  geom_sf(data=mesh.fp, fill="grey20") +
  geom_point(aes(x, y, colour=temp, size=area)) + 
  scale_colour_distiller("Surface\ntemperature", type="div", palette="RdBu", direction=-1) + 
  scale_radius(range=c(0.1, 2.5), guide="none") +
  transition_states(hour, wrap=F, transition_length=0, state_length=1) +
  facet_grid(.~mesh) +
  ggtitle(paste0("Hour: ", "{closest_state}"))
anim_save(glue("figs/mesh_surface_temperature.gif"), 
          anim, nframes=interp*24,
          fps=8, width=9, height=6, res=300, units="in")
rm(salinity.df); rm(anim)
gc()

velocity.df <- map_dfr(dir("temp", "elemLay1"), 
                       ~st_read(glue("temp/{.x}")) %>%
                         mutate(mesh=str_split_fixed(.x, "_1h_", 2)[,1])) %>%
  mutate(x=st_coordinates(.)[,1], 
         y=st_coordinates(.)[,2]) %>% 
  st_drop_geometry()
interp <- 1
anim <- velocity.df %>%
  ggplot() +
  geom_sf(data=mesh.fp, fill="grey20") +
  geom_point(aes(x, y, colour=velocity2D), size=0.3) + 
  scale_colour_viridis_c("Surface\nwaterspeed") + 
  transition_states(hour, wrap=F, transition_length=0, state_length=1) +
  facet_grid(.~mesh) +
  ggtitle(paste0("Hour: ", "{closest_state}"))
anim_save(glue("figs/mesh_surface_speed.gif"), 
          anim, nframes=interp*24,
          fps=8, width=9, height=6, res=300, units="in")
rm(velocity.df); rm(anim)
gc()


shortwave.df <- map_dfr(dir("temp", "nodeHour"), 
                       ~st_read(glue("temp/{.x}")) %>%
                         mutate(mesh=str_split_fixed(.x, "_1h_", 2)[,1])) %>%
  mutate(x=st_coordinates(.)[,1], 
         y=st_coordinates(.)[,2]) %>% 
  st_drop_geometry()
interp <- 1
anim <- shortwave.df %>%
  ggplot() +
  geom_sf(data=mesh.fp, fill="grey20") +
  geom_point(aes(x, y, colour=short_wave, size=area)) + 
  scale_colour_distiller("Shortwave\nradiation", type="seq", palette="YlOrRd") + 
  scale_radius(range=c(0.1, 2.5), guide="none") +
  transition_states(hour, wrap=F, transition_length=0, state_length=1) +
  facet_grid(.~mesh) +
  ggtitle(paste0("Hour: ", "{closest_state}"))
anim_save(glue("figs/mesh_surface_shortwave.gif"), 
          anim, nframes=interp*24,
          fps=8, width=9, height=6, res=300, units="in")
anim <- shortwave.df %>%
  ggplot() +
  geom_sf(data=mesh.fp, fill="grey20") +
  geom_point(aes(x, y, colour=zeta, size=area)) + 
  scale_colour_distiller("Zeta", type="seq", palette="Purples") + 
  scale_radius(range=c(0.1, 2.5), guide="none") +
  transition_states(hour, wrap=F, transition_length=0, state_length=1) +
  facet_grid(.~mesh) +
  ggtitle(paste0("Hour: ", "{closest_state}"))
anim_save(glue("figs/mesh_surface_zeta.gif"), 
          anim, nframes=interp*24,
          fps=8, width=9, height=6, res=300, units="in")
rm(shortwave.df); rm(anim)
gc()
