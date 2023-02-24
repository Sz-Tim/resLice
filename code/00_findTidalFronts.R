

library(tidyverse); library(ncdf4); library(sf); library(glue); library(lubridate)
library(gganimate); theme_set(theme_bw())


identify_pulse <- F
make_gif <- T
dirs <- list(
  proj=getwd(),
  mesh="D:/hydroOut/",
  linnhe_1h="D:/hydroOut/linnhe7/linnhe7_tides_met_tsobc_riv/netcdf_2021",
  linnhe_5min="D:/hydroOut/linnhe7/linnhe7_tides_met_tsobc_riv_5min/orig",
  WeStCOMS2_1h="D:/hydroOut/WeStCOMS2/Archive/netcdf_2021",
  WeStCOMS2_5min="D:/hydroOut/WeStCOMS2/Archive_5min/orig"
)
mesh.fp <- st_read("data/linnhe_mesh_footprint.gpkg")

par.df <- expand_grid(mesh=c("WeStCOMS2", "linnhe"),
                      tRes=c("1h", "5min")) %>%
  # Delta(salinity) across timesteps & buffer on FVCOM mesh nodes to define fronts
  # NB: This was subjective through trial & error...
  mutate(thresh_deltaSal=case_when(mesh=="WeStCOMS2" & tRes=="1h" ~ -2.25,
                                   mesh=="WeStCOMS2" & tRes=="5min" ~ -0.7,
                                   mesh=="linnhe" & tRes=="1h" ~ -2.25,
                                   mesh=="linnhe" & tRes=="5min" ~ -1),
         nodeBuffer=if_else(mesh=="WeStCOMS2", 500, 500), 
         mesh_gpkg=glue("{ifelse(mesh=='WeStCOMS2', 'WeStCOMS2_','')}linnhe_mesh.gpkg")) %>%
  arrange(tRes, mesh)


for(j in 1:nrow(par.df)) {
  
  if(identify_pulse) {
    # Define mesh nodes in main body of Loch Linnhe
    pulseArea.sf <- st_read("data/temp/tidal_pulse_focal_area.gpkg")
    mesh.sf <- glue("{dirs$mesh}/{par.df$mesh[j]}_mesh.nc") %>%
      nc_open() %>%
      ncvar_get("nodexy_os") %>%
      as_tibble() %>%
      st_as_sf(coords=c("V1", "V2"), crs=27700) %>%
      mutate(i=row_number(),
             inPulseArea=st_is_within_distance(., pulseArea.sf, 50, sparse=F)[,1])
    
    # Find extreme salinity change across each time step
    salSurf.mx <- map(dir(dirs[[paste(par.df$mesh[j], par.df$tRes[j], sep="_")]], full.names=T), 
                      ~ncvar_get(nc_open(.x), "salinity")[,1,]) %>%
      do.call('cbind', .)
    salProf.mx <- map(dir(dirs[[paste(par.df$mesh[j], par.df$tRes[j], sep="_")]], full.names=T), 
                      ~ncvar_get(nc_open(.x), "salinity")[ifelse(par.df$mesh[j]=="WeStCOMS2", 71462, 12464),,]) %>%
      do.call('cbind', .)
    sigLay <- ncvar_get(nc_open(dir(dirs[[paste(par.df$mesh[j], par.df$tRes[j], sep="_")]], full.names=T)[1]),
                        "siglay")[1,]
    sigLev <- ncvar_get(nc_open(dir(dirs[[paste(par.df$mesh[j], par.df$tRes[j], sep="_")]], full.names=T)[1]),
                        "siglev")[1,]
    salProf.df <- expand_grid(time=1:ncol(salProf.mx), sigLayer=sigLay) %>% 
      mutate(salinity=c(salProf.mx),
             timeStart=as_datetime("2021-11-01 00:00:00") + (time-1)*ifelse(par.df$tRes[j]=="1h", 60, 5)*60,
             timeEnd=timeStart + ifelse(par.df$tRes[j]=="1h", 60, 5)*60,
             depth=60*sigLayer,
             sigLevTop=rep(sigLev[1:length(sigLay)], ncol(salProf.mx)),
             sigLevBottom=rep(sigLev[2:length(sigLev)], ncol(salProf.mx)),
             depthTop=60*sigLevTop,
             depthBottom=60*sigLevBottom) %>%
      filter(timeEnd <= "2021-11-08 00:00:00")
    write_csv(salProf.df, glue("out/salinity_profile_{par.df$mesh[j]}_{par.df$tRes[j]}.csv"))
    p <- ggplot(salProf.df, aes(fill=salinity)) +
      geom_rect(aes(xmin=timeStart, xmax=timeEnd, ymin=depthBottom, ymax=depthTop)) + 
      scale_fill_distiller(palette="Blues", direction=1, limits=c(17, 32)) +
      scale_y_continuous("Depth (m)") +
      scale_x_datetime("")
    ggsave(glue("figs/mesh/salinity_profile_{par.df$mesh[j]}_{par.df$tRes[j]}.png"), p,
           width=7, height=5, dpi=300)
    p <- salProf.df %>%
      filter(timeEnd <= "2021-11-03 00:00:00") %>%
      ggplot(aes(fill=salinity)) +
      geom_rect(aes(xmin=timeStart, xmax=timeEnd, ymin=depthBottom, ymax=depthTop)) + 
      scale_fill_distiller(palette="Blues", direction=1, limits=c(17, 32)) +
      scale_y_continuous("Depth (m)") +
      scale_x_datetime("")
    ggsave(glue("figs/mesh/salinity_profile_2days_{par.df$mesh[j]}_{par.df$tRes[j]}.png"), p,
           width=7, height=5, dpi=300)
    deltaSal <- salSurf.mx[,2:ncol(salSurf.mx)] - salSurf.mx[,1:(ncol(salSurf.mx)-1)]
    extremes.i <- map(1:ncol(deltaSal), ~which(deltaSal[,.x] <= par.df$thresh_deltaSal[j] &
                                                 deltaSal[,.x] < quantile(deltaSal[,.x], 0.01)))
    
    # Buffer extreme nodes to a single buffered polygon per time step
    extremes.sf <- map_dfr(1:length(extremes.i),
                           ~mesh.sf[extremes.i[[.x]],] %>%
                             mutate(trans=.x)) %>%
      filter(inPulseArea) %>%
      st_buffer(dist=par.df$nodeBuffer[j]) %>%
      group_by(trans) %>%
      summarise(geometry=st_union(geometry)) %>%
      st_cast("MULTIPOLYGON") %>%
      st_cast("POLYGON") %>%
      mutate(area=st_area(.)) %>%
      group_by(trans) %>%
      arrange(desc(area)) %>%
      slice_head(n=1) %>%
      select(-area) %>%
      mutate(timeStart=as_datetime("2021-11-01 00:00:00") + (trans-1)*ifelse(par.df$tRes[j]=="1h", 60, 5)*60) %>%
      st_write(glue("data/salinityPulse_{par.df$mesh[j]}_{par.df$tRes[j]}.gpkg"), append=FALSE)  
  }
  if(make_gif) {
    extremes.sf <- st_read(glue("data/salinityPulse_{par.df$mesh[j]}_{par.df$tRes[j]}.gpkg"))
    # invisible point at each time step to include steps with no fronts
    time.df <- tibble(trans=1:max(extremes.sf$trans)) %>%
      mutate(timeStart=as_datetime("2021-11-01 00:00:00") + (trans-1)*ifelse(par.df$tRes[j]=="1h", 60, 5)*60,
             x=201850,
             y=763450)
    
    # Animate
    anim <- ggplot() +
      geom_sf(data=mesh.fp, fill="grey85") +
      geom_point(data=time.df, aes(x, y), size=0.25, col="grey85") +
      geom_sf(data=extremes.sf, fill="red3", colour=NA) +
      transition_states(timeStart, transition_length=0, state_length=1) +
      ggtitle("{closest_state}")
    anim_save(glue("figs/mesh/tidal_pulse_{par.df$mesh[j]}_{par.df$tRes[j]}.gif"),
              anim, width=5, height=6, unit="in", res=300,
              fps=ifelse(par.df$tRes[j]=="1h", 3, 12),
              # fps=3,
              nframes=nrow(time.df))
  }
  
}

# 1h animations
extremes.sf <- map_dfr(unique(par.df$mesh), 
                       ~st_read(glue("data/salinityPulse_{.x}_1h.gpkg")) %>%
                         mutate(mesh=.x))
time.df <- tibble(trans=1:max(extremes.sf$trans)) %>%
  mutate(timeStart=as_datetime("2021-11-01 00:00:00") + (trans-1)*60*60,
         x=201850,
         y=763450)
anim <- ggplot() +
  geom_sf(data=mesh.fp, fill="grey85") +
  geom_point(data=time.df, aes(x, y), size=0.25, col="grey85") +
  geom_sf(data=extremes.sf, fill="blue3", colour=NA) +
  transition_states(timeStart, transition_length=0, state_length=1) +
  ggtitle("{closest_state}") +
  facet_grid(.~mesh) +
  theme(axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank())
anim_save(glue("figs/mesh/tidal_pulse_1h.gif"),
          anim, width=9, height=6, unit="in", res=300,
          fps=3, nframes=nrow(time.df))


# 5min animations
extremes.sf <- map_dfr(unique(par.df$mesh), 
                       ~st_read(glue("data/salinityPulse_{.x}_5min.gpkg")) %>%
                         mutate(mesh=.x))
time.df <- tibble(trans=1:max(extremes.sf$trans)) %>%
  mutate(timeStart=as_datetime("2021-11-01 00:00:00") + (trans-1)*5*60,
         x=201850,
         y=763450)
anim <- ggplot() +
  geom_sf(data=mesh.fp, fill="grey85") +
  geom_point(data=time.df, aes(x, y), size=0.25, col="grey85") +
  geom_sf(data=extremes.sf, fill="blue3", colour=NA) +
  transition_states(timeStart, transition_length=0, state_length=1) +
  ggtitle("{closest_state}") +
  facet_grid(.~mesh) +
  theme(axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank())
anim_save(glue("figs/mesh/tidal_pulse_5min.gif"),
          anim, width=9, height=6, unit="in", res=300,
          fps=36, nframes=nrow(time.df))

