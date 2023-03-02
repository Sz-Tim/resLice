


rm(list=ls())
gc()


library(tidyverse); library(glue); library(lubridate); library(sf); library(gganimate)
source("code/00_fn.R")
theme_set(theme_classic())

initDensity <- c("Scaled", "Uniform")[2]

dirs <- switch(get_os(),
               windows=list(proj=getwd(),
                            mesh="D:/hydroOut/",
                            out=glue("{getwd()}/out/siteRelease_init{initDensity}/")),
               linux=list(proj=getwd(),
                          mesh="/home/sa04ts/FVCOM_meshes",
                          out=glue("{getwd()}/out/siteRelease_init{initDensity}/")))

mesh.fp <- st_read("data/linnhe_mesh_footprint.gpkg")
mesh.sf <- list(linnhe7=st_read(glue("{dirs$mesh}/linnhe_mesh.gpkg")) %>%
                  mutate(mesh="linnhe7"),
                westcoms2=st_read(glue("{dirs$mesh}/WeStCOMS2_linnhe_mesh.gpkg")) %>%
                  mutate(mesh="WeStCOMS2")) %>%
  do.call('rbind', .) %>%
  select(i, area, depth, mesh, lochRegion, geom) %>%
  full_join(., read_csv("data/loch_regions.csv"), by="lochRegion")
sim_i <- read_csv(glue("{dirs$out}/sim_i.csv")) %>%
  mutate(liceSpeedF=factor(liceSpeed, levels=c(0.0001, 0.0005, 0.001), 
                           labels=c("Slow", "Medium", "Fast")))
site.i <- read_tsv("data/fishFarmSites.tsv", col_names=c("site", "x", "y"))
col_meshRes <- c("linnhe7, 1h"="#a6611a", "linnhe7, 5min"="#dfc27d",
                 "WeStCOMS2, 1h"="#018571", "WeStCOMS2, 5min"="#80cdc1")

loc.df <- readRDS(glue("out/00_processed/locations_site_init{initDensity}.rds")) %>%
  filter(status != 66) %>%
  mutate(liceSpeedF=factor(liceSpeedF, levels=levels(sim_i$liceSpeedF)),
         meshRes=factor(paste0(mesh, ", ", timeRes),
                        levels=c("WeStCOMS2, 1h", "linnhe7, 1h", "linnhe7, 5min", "WeStCOMS2, 5min")),
         hour=hour(timeCalculated),
         date=date(timeCalculated))
days <- unique(loc.df$date)
hours <- unique(loc.df$hour)
particles <- unique(filter(loc.df, date %in% days & hour %in% hours)$ID)
part.sample <- sample(particles, min(length(particles), 100e3))

loc.df <- loc.df %>%
  filter(date %in% days,
         hour %in% hours,
         ID %in% part.sample)
interp <- 3



# Etive focus
etive.bbox <- list(xmin=185276, xmax=212845, ymin=728357, ymax=746971)
loc.etive <- loc.df %>% 
  filter(x >= (etive.bbox$xmin - 2000),
         x <= (etive.bbox$xmax + 2000),
         y >= (etive.bbox$ymin - 2000),
         y <= (etive.bbox$ymax + 2000)) %>%
  mutate(meshRes=factor(meshRes, levels=names(col_meshRes)))
rm(loc.df)
interp <- 3
days <- unique(loc.etive$date)
hours <- unique(loc.etive$hour)
anim <- loc.etive %>%
  filter(status==2) %>%
  filter(liceSpeedF=="Medium") %>%
  # filter(date == "2021-11-06") %>%
  mutate(startLocation=if_else(startLocation %in% c("FS1067", "FS1101", "FS1112"),
                               "in Loch Etive", "outside Loch Etive")) %>%
  ggplot() +
  geom_sf(data=mesh.fp, fill="grey10", colour=NA) +
  geom_point(aes(x, y, colour=startLocation, alpha=density, group=ID), size=0.8) +
  geom_point(data=site.i, aes(x, y), colour="red", shape=13, size=4) +
  scale_colour_manual(values=c("#a6cee3", "#b2df8a")) +
  transition_states(timeCalculated, wrap=F, transition_length=1, state_length=0) +
  facet_wrap(~meshRes) +
  xlim(etive.bbox$xmin, etive.bbox$xmax) +
  ylim(etive.bbox$ymin, etive.bbox$ymax) +
  ggtitle("Copepodid particles: {closest_state}") +
  theme(axis.text=element_blank(),
        axis.title=element_blank())
anim_save(glue("figs/init{initDensity}_reps/tracks_etive.gif"),
          anim, nframes=interp*length(hours)*length(days),
          fps=32, width=10, height=7, res=500, units="in")


loc.etive %>%
  filter(status==2) %>%
  filter(liceSpeedF=="Medium") %>%
  filter(date == "2021-11-06") %>%
  filter(hour < 12) %>%
  filter(meshRes=="linnhe7, 1h") %>%
  mutate(startLocation=if_else(startLocation %in% c("FS1067", "FS1101", "FS1112"),
                               "in Loch Etive", "outside Loch Etive")) %>%
  ggplot() +
  geom_sf(data=mesh.fp, fill="grey10", colour=NA) +
  geom_point(aes(x, y, colour=startLocation, alpha=density, group=ID), size=1) +
  geom_point(data=site.i, aes(x, y), colour="#7570b3", shape=1, size=4) +
  scale_colour_manual(values=c("#d95f02", "#1b9e77")) +
  facet_wrap(~timeCalculated, nrow=3) +
  xlim(etive.bbox$xmin, etive.bbox$xmax) +
  ylim(etive.bbox$ymin, etive.bbox$ymax) +
  ggtitle("Copepodid particles: {closest_state}") +
  theme(axis.text=element_blank(),
        axis.title=element_blank())


interp <- 1
dens.sf <- loc.df %>% 
  filter(liceSpeedF=="Medium") %>%
  st_as_sf(coords=c("x", "y"), crs=27700) %>%
  st_join(mesh.sf %>% filter(mesh=="WeStCOMS2") %>% select(i)) %>%
  st_drop_geometry() %>%
  group_by(sim, meshRes, liceSpeedF, timeCalculated, i) %>%
  summarise(density=sum(density), N=n()) %>%
  ungroup %>%
  inner_join(mesh.sf %>% filter(mesh=="WeStCOMS2"), .) %>%
  mutate(density=log1p(density/(area*depth)))
anim <- dens.sf %>%
  ggplot() + 
  geom_sf(data=mesh.fp, fill=NA, colour="grey30") +
  geom_sf(aes(fill=density), colour=NA) + 
  scale_fill_viridis_c() +
  transition_states(timeCalculated, wrap=F, transition_length=0, state_length=1) +
  facet_grid(liceSpeedF~meshRes) +
  theme(axis.text=element_blank(),
        legend.position="bottom") +
  ggtitle( "Lice density: {closest_state}")
anim_save(glue("figs/initScaled/density_slower_site.gif"),
          anim, nframes=interp*length(hours)*length(days),
          fps=12, width=8, height=4, res=300, units="in")
rm(dens.sf); gc()

for(i in seq_along(sim_i)) {
  dens.sf <- loc.df %>% 
    filter(sim==i) %>%
    st_as_sf(coords=c("x", "y"), crs=27700) %>%
    st_join(mesh.sf %>% filter(mesh=="WeStCOMS2") %>% select(i)) %>%
    st_drop_geometry() %>%
    group_by(sim, meshRes, liceSpeedF, timeCalculated, i) %>%
    summarise(density=sum(density), N=n()) %>%
    ungroup %>%
    inner_join(mesh.sf %>% filter(mesh=="WeStCOMS2"), .) %>%
    mutate(density=log1p(density/(area*depth)))
  anim <- dens.sf %>%
    ggplot() + 
    geom_sf(data=mesh.fp, fill=NA, colour="grey30") +
    geom_sf(aes(fill=density), colour=NA) + 
    scale_fill_viridis_c() +
    transition_states(timeCalculated, wrap=F, transition_length=0, state_length=1) +
    facet_grid(liceSpeedF~meshRes) +
    theme(axis.text=element_blank(),
          legend.position="bottom") +
    ggtitle( "Lice density: {closest_state}")
  anim_save(glue("figs/initScaled/density_slower_site_sim{i}.gif"),
            anim, nframes=interp*length(hours)*length(days),
            fps=12, width=4, height=4, res=300, units="in")
  rm(dens.sf); gc()
}






anim <- loc.df %>%
  filter(age >= 6) %>%
  filter(depth < 20) %>%
  ggplot(aes(-depth, colour=meshRes)) +
  geom_density(adjust=1.5) + 
  scale_colour_brewer(type="div") +
  transition_states(timeCalculated, wrap=F, transition_length=1, state_length=0) +
  ggtitle( "{closest_state}") +
  facet_grid(.~liceSpeedF) +
  theme(legend.position="bottom") +
  coord_flip()
anim_save(glue("figs/initScaled/vertDistribution_site.gif"),
          anim, nframes=interp*length(hours)*length(days),
          fps=24, width=9, height=4, res=300, units="in")

# anim <- loc.df %>%
#   filter(age >= 6) %>%
#   filter(depth < 20) %>%
#   ggplot(aes(-depth, colour=liceSpeedF)) +
#   geom_density(adjust=1.5) + 
#   scale_colour_viridis_d(end=0.9) +
#   transition_states(timeCalculated, wrap=F, transition_length=1, state_length=0) +
#   ggtitle( "{closest_state}") +
#   facet_grid(.~meshRes) +
#   theme(legend.position="bottom") +
#   coord_flip()
# anim_save(glue("figs/initScaled/vertDistribution2_site.gif"),
#           anim, nframes=interp*length(hours)*length(days),
#           fps=24, width=9, height=4, res=300, units="in")

anim <- loc.df %>%
  # filter(timeRes=="1h") %>%
  filter(age >= 12) %>%
  ggplot(aes(xyTot, colour=meshRes)) +
  geom_density(adjust=1.3) +
  scale_colour_brewer(type="div") +
  scale_x_continuous("Cumulative xy movement (m)") +
  transition_states(timeCalculated, wrap=F, transition_length=1, state_length=0) +
  ggtitle( "{closest_state}") +
  facet_grid(.~liceSpeedF) +
  theme(legend.position="bottom")
anim_save(glue("figs/initScaled/xy_movement_site.gif"),
          anim, nframes=interp*length(hours)*length(days),
          fps=24, width=9, height=4, res=300, units="in")

anim <- loc.df %>%
  filter(age >= 12) %>%
  ggplot(aes(xyTot, colour=liceSpeedF)) +
  geom_density(adjust=1.3) +
  scale_colour_brewer(type="qual", palette=2) +
  scale_x_continuous("Cumulative xy movement (m)") +
  transition_states(timeCalculated, wrap=F, transition_length=1, state_length=0) +
  ggtitle( "{closest_state}") +
  facet_grid(.~meshRes) +
  theme(legend.position="bottom")
anim_save(glue("figs/initScaled/xy_movement2_site.gif"),
          anim, nframes=interp*length(hours)*length(days),
          fps=24, width=9, height=4, res=300, units="in")

anim <- loc.df %>%
  filter(timeRes=="1h") %>%
  filter(age >= 12) %>%
  ggplot(aes(xyTot, colour=meshRes)) +
  geom_density(adjust=1.3) +
  scale_colour_brewer(type="div") +
  scale_x_continuous("Cumulative xy movement (m)", trans="log1p") +
  transition_states(timeCalculated, wrap=F, transition_length=1, state_length=0) +
  ggtitle("{closest_state}") +
  facet_grid(.~liceSpeedF) +
  theme(legend.position="bottom")
anim_save(glue("figs/initScaled/xy_ln_movement_site.gif"),
          anim, nframes=interp*length(hours)*length(days),
          fps=24, width=9, height=4, res=300, units="in")

anim <- loc.df %>%
  filter(age >= 12) %>%
  ggplot(aes(xyTot, colour=liceSpeedF)) +
  geom_density() +
  scale_colour_brewer(type="qual", palette=2) +
  scale_x_continuous("Cumulative xy movement (m)", trans="log1p") +
  transition_states(timeCalculated, wrap=F, transition_length=1, state_length=0) +
  ggtitle("{closest_state}") +
  facet_grid(.~meshRes) +
  theme(legend.position="bottom")
anim_save(glue("figs/initScaled/xy_ln_movement2_site.gif"),
          anim, nframes=interp*length(hours)*length(days),
          fps=24, width=9, height=4, res=300, units="in")

# 
# anim <- loc.df %>%
#   ggplot() + 
#   geom_sf(data=mesh.fp, fill="grey10", colour=NA) +
#   geom_point(aes(x, y, colour=depth, alpha=density, group=ID), size=0.25) + 
#   scale_colour_viridis_c(direction=-1, option="D", limits=c(0,20), na.value="#000033") +
#   transition_states(timeCalculated, wrap=F, transition_length=1, state_length=0) +
#   facet_grid(meshRes~liceSpeedF) +
#   ggtitle( "{closest_state}")
# anim_save(glue("figs/initScaled/tracks_AllPart_new_site.gif"), 
#           anim, nframes=interp*length(hours)*length(days),
#           fps=32, width=6, height=9, res=300, units="in")
# 
# 
# 
# loc.df <- loc.df %>%
#   filter(liceSpeedF!="Passive") %>%
#   filter(date %in% days,
#          hour %in% hours,
#          ID %in% part.sample)
# interp <- 8
# anim <- loc.df %>%
#   ggplot() + 
#   geom_sf(data=mesh.fp, fill="grey10", colour=NA) +
#   geom_point(aes(x, y, colour=depth, alpha=density, group=ID), size=0.25) + 
#   scale_colour_viridis_c(direction=-1, option="D", limits=c(0,20), na.value="#000033") +
#   transition_states(timeCalculated, wrap=F, transition_length=1, state_length=0) +
#   facet_grid(.~meshRes) +
#   ggtitle( "{closest_state}")
# anim_save(glue("figs/initScaled/tracks_AllPart_nonPassive_site.gif"), 
#           anim, nframes=interp*length(hours)*length(days),
#           fps=32, width=9, height=4, res=300, units="in")







interp <- 8
velocity.df <- loc.df %>%
  filter(liceSpeedF=="Medium") %>%
  filter(date %in% days,
         hour %in% hours,
         ID %in% part.sample) %>%
  select(ID, timeCalculated, age, status, x, y, density, depth, 
         sim, mesh, meshRes, elem) %>%
  arrange(sim, timeCalculated, ID) %>%
  group_by(ID, sim) %>%
  mutate(x_m1=lag(x), 
         y_m1=lag(y),
         depth_m1=lag(depth)) %>%
  ungroup %>%
  mutate(dx=x-x_m1,
         dy=y-y_m1,
         dz=depth-depth_m1,
         dXY=sqrt(dx^2 + dy^2),
         bearing=atan2(dy, dx)) %>%
  filter(status != 66) %>%
  left_join(., sim_i %>% mutate(sim=as.numeric(i)) %>% select(sim, mesh, liceSpeedF)) %>%
  left_join(., mesh.sf %>% st_drop_geometry() %>% 
              select(i, mesh, bearing, depth) %>% 
              rename(elem=i, lochDir=bearing, meshDepth=depth),
            by=c("mesh", "elem")) %>%
  mutate(downloch=cos(lochDir-bearing),
         depth=case_when(depth < 2 ~ "0-2",
                         between(depth, 2, 10) ~ "2-10",
                         depth > 10 ~ "10-X"), 
         depth=factor(depth, levels=c("0-2", "2-10", "10-X")))

# anim <- velocity.df %>%
#   ggplot() + 
#   geom_sf(data=mesh.fp, fill="grey10", colour=NA) +
#   geom_point(aes(x, y, colour=downloch, alpha=density, group=ID), size=0.25) + 
#   scale_colour_gradient2("Direction", limits=c(-1,1), 
#                          breaks=c(-1,1), 
#                          labels=c("uploch", "downloch")) +
#   transition_states(timeCalculated, wrap=F, transition_length=1, state_length=0) +
#   facet_grid(depth~meshRes) +
#   ggtitle( "{closest_state}")
# anim_save(glue("figs/initScaled/tracks_Dir_nonPassive_site.gif"), 
#           anim, nframes=interp*length(hours)*length(days),
#           fps=32, width=9, height=10, res=300, units="in")
# 
# rm(velocity.df)
# gc()
# elem.Depth <- loc.df %>%
#   filter(status != 66, liceSpeedF=="Medium") %>% 
#   group_by(sim, meshRes, timeCalculated, elem) %>% 
#   summarise(mdDepth=median(depth)) %>%
#   ungroup %>%
#   left_join(., sim_i %>% mutate(sim=as.numeric(i)) %>% select(sim, mesh, liceSpeedF)) %>%
#   left_join(mesh.sf %>% 
#               select(i, mesh, depth) %>% 
#               rename(elem=i, meshDepth=depth),
#             ., by=c("mesh", "elem"))
# anim <- elem.Depth %>%
#   mutate(date=date(timeCalculated), hour=hour(timeCalculated)) %>%
#   filter(date %in% days,
#          hour %in% hours) %>%
#   ggplot() + 
#   geom_sf(data=mesh.fp, fill="grey10", colour=NA) +
#   geom_sf(aes(fill=mdDepth, group=elem), colour=NA) +
#   scale_fill_viridis_c("Median\nparticle\ndepth", direction=-1, option="D", 
#                        limits=c(0,20), na.value="#000033") +
#   transition_states(timeCalculated, wrap=F, transition_length=1, state_length=0) +
#   facet_grid(.~meshRes) +
#   ggtitle( "{closest_state}")
# anim_save(glue("figs/initScaled/mdDepth_map_nonPassive_site.gif"), 
#           anim, nframes=interp*length(hours)*length(days),
#           fps=32, width=9, height=4, res=300, units="in")
# anim <- elem.Depth %>%
#   mutate(date=date(timeCalculated), hour=hour(timeCalculated)) %>%
#   filter(date %in% days,
#          hour %in% hours) %>%
#   ggplot() + 
#   geom_sf(data=mesh.fp, fill="grey10", colour=NA) +
#   geom_sf(aes(fill=mdDepth/meshDepth, group=elem), colour=NA) +
#   scale_fill_viridis_c("Median\nparticle\ndepth (%)", direction=-1, option="D", 
#                        limits=c(0,1)) +
#   transition_states(timeCalculated, wrap=F, transition_length=1, state_length=0) +
#   facet_grid(.~meshRes) +
#   ggtitle( "{closest_state}")
# anim_save(glue("figs/initScaled/mdDepthPct_map_nonPassive_site.gif"), 
#           anim, nframes=interp*length(hours)*length(days),
#           fps=32, width=9, height=4, res=300, units="in")
# rm(elem.Depth)
# gc()
# 
loc.ls <- loc.df %>%
  filter(date %in% days,
         hour %in% hours,
         ID %in% part.sample) %>%
  group_by(sim) %>%
  group_split()
rm(loc.df)
gc()



library(doSNOW)
cl <- makeCluster(min(nrow(sim_i), 1))
registerDoSNOW(cl)
foreach(i=1:nrow(sim_i),
        .packages=c("tidyverse", "glue", "lubridate", "sf", "gganimate"),
        .export=c("sim_i", "mesh.fp", "loc.ls", "days", "hours"),
        .errorhandling="pass") %dopar% {

  interp <- 8  # interpolation frames
  theme_set(theme_classic())

  # anim <- loc.ls[[i]] %>%
  #   st_as_sf(coords=c("x", "y"), crs=27700) %>%
  #   ggplot() +
  #   geom_sf(data=mesh.fp, fill="grey10", colour=NA) +
  #   geom_sf(aes(colour=depth, alpha=density, group=ID), size=0.5) +
  #   scale_colour_viridis_c(direction=-1, option="D", limits=c(0,20), na.value="#000033") +
  #   transition_states(timeCalculated, wrap=F, transition_length=1, state_length=0) +
  #   ggtitle(paste0(glue("{sim_i$mesh[i]}, {sim_i$timeRes[i]}, {sim_i$liceSpeedF[i]}"),
  #                  ". {closest_state}"))
  
  bbox <- st_bbox(mesh.fp)
  
  anim <- loc.ls[[i]] %>%
    filter(between(x, bbox$xmin, bbox$xmax), between(y, bbox$ymin, bbox$ymax)) %>%
    ggplot() +
    geom_sf(data=mesh.fp, fill="grey10", colour=NA) +
    geom_point(aes(x, y, colour=depth, alpha=density, group=ID), size=0.5) +
    scale_colour_viridis_c(direction=-1, option="D", limits=c(0,20), na.value="#000033") +
    transition_states(timeCalculated, wrap=F, transition_length=1, state_length=0) +
    ggtitle(paste0(glue("{sim_i$mesh[i]}, {sim_i$timeRes[i]}, {sim_i$liceSpeedF[i]}"),
                   ". {closest_state}"))
  anim_save(glue("figs/initScaled/tracks_sim{str_pad(i,2,'left','0')}_10kPart_site.gif"),
            anim, nframes=interp*length(hours)*length(days),
            fps=32, width=7, height=6, res=300, units="in")
  paste("Finished", i)
}
stopCluster(cl)
# 
# 
# 
# 
# vel.ls <- velocity.df %>%
#   group_by(meshRes, liceSpeedF, depth) %>%
#   group_split()
# rm(velocity.df)
# rm(loc.df)
# gc()
# 
# library(doSNOW)
# cl <- makeCluster(5)
# registerDoSNOW(cl)
# foreach(i=seq_along(vel.ls),
#         .packages=c("tidyverse", "glue", "lubridate", "sf", "gganimate"),
#         .export=c("sim_i", "mesh.fp", "vel.ls", "days", "hours"),
#         .errorhandling="pass") %dopar% {
# 
#           interp <- 8  # interpolation frames
#           theme_set(theme_classic())
#           info <- tibble(meshRes=vel.ls[[i]]$meshRes[1],
#                          liceSpeedF=vel.ls[[i]]$liceSpeedF[1],
#                          depth=vel.ls[[i]]$depth[1]) %>%
#             mutate(meshRes=str_replace(meshRes, ", ", "_"))
# 
#           anim <- vel.ls[[i]] %>%
#             ggplot() +
#             geom_sf(data=mesh.fp, fill="grey10", colour=NA) +
#             geom_point(aes(x, y, colour=downloch, alpha=density, group=ID), size=0.25) +
#             scale_colour_gradient2("Direction", limits=c(-1,1),
#                                    breaks=c(-1,1),
#                                    labels=c("uploch", "downloch")) +
#             transition_states(timeCalculated, wrap=F, transition_length=1, state_length=0) +
#             ggtitle(paste0(glue("{info$meshRes}, {info$liceSpeedF}, {info$depth}m"),
#                            ". {closest_state}"))
#           anim_save(glue("figs/initScaled/tracks_dir_{info$meshRes}_{info$depth}m_site.gif"),
#                     anim, nframes=interp*length(hours)*length(days),
#                     fps=32, width=7, height=6, res=300, units="in")
#           paste("Finished", i)
#         }
# stopCluster(cl)

vel.ls <- velocity.df %>%
  group_by(liceSpeedF, depth) %>%
  group_split()
rm(velocity.df)
rm(loc.df)
gc()

library(doSNOW)
cl <- makeCluster(3)
registerDoSNOW(cl)
foreach(i=seq_along(vel.ls),
        .packages=c("tidyverse", "glue", "lubridate", "sf", "gganimate"),
        .export=c("sim_i", "mesh.fp", "vel.ls", "days", "hours"),
        .errorhandling="pass") %dopar% {
          
          interp <- 8  # interpolation frames
          theme_set(theme_classic())
          info <- tibble(liceSpeedF=vel.ls[[i]]$liceSpeedF[1],
                         depth=vel.ls[[i]]$depth[1])
          
          bbox <- st_bbox(mesh.fp)
          
          anim <- vel.ls[[i]] %>%
            filter(between(x, bbox$xmin, bbox$xmax), between(y, bbox$ymin, bbox$ymax)) %>%
            ggplot() +
            geom_sf(data=mesh.fp, fill="grey10", colour=NA) +
            geom_point(aes(x, y, colour=downloch, alpha=density, group=ID), size=0.25) +
            scale_colour_gradient2("Direction", limits=c(-1,1),
                                   breaks=c(-1,1),
                                   labels=c("uploch", "downloch")) +
            transition_states(timeCalculated, wrap=F, transition_length=1, state_length=0) +
            facet_grid(.~meshRes) +
            ggtitle(paste0(glue("{info$liceSpeedF}, {info$depth}m"),
                           ". {closest_state}"))
          anim_save(glue("figs/initScaled/tracks_dir_{info$depth}m_{info$liceSpeedF}_site.gif"),
                    anim, nframes=interp*length(hours)*length(days),
                    fps=32, width=9, height=4, res=300, units="in")
          paste("Finished", i)
        }
stopCluster(cl)


# 
# 
# vel.ls <- velocity.df %>%
#   group_by(meshRes) %>%
#   group_split()
# rm(velocity.df)
# rm(loc.df)
# gc()
# 
# library(doSNOW)
# cl <- makeCluster(3)
# registerDoSNOW(cl)
# foreach(i=seq_along(vel.ls), 
#         .packages=c("tidyverse", "glue", "lubridate", "sf", "gganimate"),
#         .export=c("sim_i", "mesh.fp", "vel.ls", "days", "hours"),
#         .errorhandling="pass") %dopar% {
#           
#           interp <- 8  # interpolation frames
#           theme_set(theme_classic())
#           info <- tibble(meshRes=vel.ls[[i]]$meshRes[1],
#                          liceSpeedF=vel.ls[[i]]$liceSpeedF[1],
#                          depth=vel.ls[[i]]$depth[1]) %>%
#             mutate(meshRes=str_replace(meshRes, ", ", "_"))
#           
#           anim <- vel.ls[[i]] %>%
#             ggplot() + 
#             geom_sf(data=mesh.fp, fill="grey10", colour=NA) +
#             geom_point(aes(x, y, colour=downloch, alpha=density, group=ID), size=0.25) + 
#             scale_colour_gradient2("Direction", limits=c(-1,1), 
#                                    breaks=c(-1,1), 
#                                    labels=c("uploch", "downloch")) +
#             transition_states(timeCalculated, wrap=F, transition_length=1, state_length=0) +
#             ggtitle(paste0(glue("{info$meshRes}, {info$liceSpeedF}"), 
#                            ". {closest_state}"))
#           anim_save(glue("figs/initScaled/tracks_dir_{info$meshRes}_site.gif"), 
#                     anim, nframes=interp*length(hours)*length(days),
#                     fps=32, width=7, height=6, res=300, units="in")
#           paste("Finished", i)
#         }
# stopCluster(cl)
