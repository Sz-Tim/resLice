


library(tidyverse); library(glue); library(lubridate); library(sf); library(gganimate)
source("code/00_fn.R")
theme_set(theme_classic())

dirs <- switch(get_os(),
               windows=list(proj=getwd(),
                            mesh="D:/hydroOut/",
                            out=glue("{getwd()}/out/")),
               linux=list(proj=getwd(),
                          mesh="/home/sa04ts/FVCOM_meshes",
                          out=glue("{getwd()}/out/")))

mesh.fp <- st_read("data/linnhe_mesh_footprint.gpkg")
mesh.sf <- list(linnhe7=st_read(glue("{dirs$mesh}/linnhe_mesh.gpkg")) %>%
                  mutate(mesh="linnhe7"),
                westcoms2=st_read(glue("{dirs$mesh}/WeStCOMS2_linnhe_mesh.gpkg")) %>%
                  mutate(mesh="WeStCOMS2")) %>%
  do.call('rbind', .) %>%
  select(i, area, depth, mesh, lochRegion, geom) %>%
  full_join(., read_csv("data/loch_regions.csv"), by="lochRegion")
sim_i <- read_csv(glue("{dirs$out}/sim_i.csv")) %>%
  mutate(liceSpeedF=factor(liceSpeed, labels=c("Passive", "Medium")))

loc.df <- readRDS("out/00_processed/locations.rds") %>%
  filter(status != 66) %>%
  mutate(liceSpeedF=factor(liceSpeedF, levels=levels(sim_i$liceSpeedF)),
         meshRes=factor(paste0(mesh, ", ", timeRes),
                        levels=c("linnhe7, 1h", "linnhe7, 5min",
                                 "WeStCOMS2, 1h", "WeStCOMS2, 5min")),
         hour=hour(timeCalculated),
         date=date(timeCalculated))
days <- unique(loc.df$date)
hours <- unique(loc.df$hour)
particles <- unique(filter(loc.df, date %in% days & hour %in% hours)$ID)
part.sample <- sample(particles, min(length(particles), 1e5))

loc.df <- loc.df %>%
  filter(date %in% days,
         hour %in% hours,
         ID %in% part.sample)
interp <- 8


anim <- loc.df %>%
  filter(liceSpeedF=="Medium") %>%
  filter(age > 6) %>%
  ggplot(aes(-depth, colour=meshRes)) +
  geom_density() + xlim(-20, 0) +
  scale_colour_brewer(type="qual", palette=2) +
  transition_states(timeCalculated, wrap=F, transition_length=1, state_length=0) +
  ggtitle( "{closest_state}") +
  coord_flip()
anim_save(glue("figs/vertDistribution.gif"),
          anim, nframes=interp*length(hours)*length(days),
          fps=32, width=5, height=4, res=300, units="in")

anim <- loc.df %>%
  filter(liceSpeedF=="Medium") %>%
  filter(timeRes=="1h") %>%
  filter(age > 6) %>%
  ggplot(aes(xyTot, colour=meshRes)) +
  geom_density() +
  scale_colour_brewer(type="qual", palette=2) +
  scale_x_continuous("Cumulative xy movement (m)") +
  transition_states(timeCalculated, wrap=F, transition_length=1, state_length=0) +
  ggtitle( "{closest_state}")
anim_save(glue("figs/xy_movement.gif"),
          anim, nframes=interp*length(hours)*length(days),
          fps=32, width=6, height=4, res=300, units="in")

anim <- loc.df %>%
  filter(liceSpeedF=="Medium") %>%
  filter(timeRes=="1h") %>%
  filter(age > 6) %>%
  ggplot(aes(xyTot, colour=meshRes)) +
  geom_density() +
  scale_colour_brewer(type="qual", palette=2) +
  scale_x_continuous("Cumulative xy movement (m)", trans="log1p") +
  transition_states(timeCalculated, wrap=F, transition_length=1, state_length=0) +
  ggtitle( "{closest_state}")
anim_save(glue("figs/xy_ln_movement.gif"),
          anim, nframes=interp*length(hours)*length(days),
          fps=32, width=5, height=4, res=300, units="in")

# 
# anim <- loc.df %>%
#   ggplot() + 
#   geom_sf(data=mesh.fp, fill="grey10", colour=NA) +
#   geom_point(aes(x, y, colour=depth, alpha=density, group=ID), size=0.25) + 
#   scale_colour_viridis_c(direction=-1, option="D", limits=c(0,20), na.value="#000033") +
#   transition_states(timeCalculated, wrap=F, transition_length=1, state_length=0) +
#   facet_grid(meshRes~liceSpeedF) +
#   ggtitle( "{closest_state}")
# anim_save(glue("figs/tracks_AllPart_new.gif"), 
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
# anim_save(glue("figs/tracks_AllPart_nonPassive.gif"), 
#           anim, nframes=interp*length(hours)*length(days),
#           fps=32, width=9, height=4, res=300, units="in")



interp <- 8
velocity.df <- loc.df %>%
  filter(liceSpeedF!="Passive") %>%
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
# anim_save(glue("figs/tracks_Dir_nonPassive.gif"), 
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
# anim_save(glue("figs/mdDepth_map_nonPassive.gif"), 
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
# anim_save(glue("figs/mdDepthPct_map_nonPassive.gif"), 
#           anim, nframes=interp*length(hours)*length(days),
#           fps=32, width=9, height=4, res=300, units="in")
# rm(elem.Depth)
# gc()
# 
# loc.ls <- loc.df %>%
#   filter(date %in% days,
#          hour %in% hours,
#          ID %in% part.sample) %>%
#   group_by(sim) %>%
#   group_split()
# rm(loc.df)
# gc()
# 
# 
# 
# library(doSNOW)
# cl <- makeCluster(min(nrow(sim_i), 10))
# registerDoSNOW(cl)
# foreach(i=1:nrow(sim_i), 
#         .packages=c("tidyverse", "glue", "lubridate", "sf", "gganimate"),
#         .export=c("sim_i", "mesh.fp", "loc.ls", "days", "hours"),
#         .errorhandling="pass") %dopar% {
#           
#   interp <- 8  # interpolation frames
#   theme_set(theme_classic())
#   
#   # anim <- loc.ls[[i]] %>%
#   #   st_as_sf(coords=c("x", "y"), crs=27700) %>%
#   #   ggplot() + 
#   #   geom_sf(data=mesh.fp, fill="grey10", colour=NA) +
#   #   geom_sf(aes(colour=depth, alpha=density, group=ID), size=0.5) + 
#   #   scale_colour_viridis_c(direction=-1, option="D", limits=c(0,20), na.value="#000033") +
#   #   transition_states(timeCalculated, wrap=F, transition_length=1, state_length=0) +
#   #   ggtitle(paste0(glue("{sim_i$mesh[i]}, {sim_i$timeRes[i]}, {sim_i$liceSpeedF[i]}"), 
#   #                  ". {closest_state}"))
#   anim <- loc.ls[[i]] %>%
#     ggplot() + 
#     geom_sf(data=mesh.fp, fill="grey10", colour=NA) +
#     geom_point(aes(x, y, colour=depth, alpha=density, group=ID), size=0.5) + 
#     scale_colour_viridis_c(direction=-1, option="D", limits=c(0,20), na.value="#000033") +
#     transition_states(timeCalculated, wrap=F, transition_length=1, state_length=0) +
#     ggtitle(paste0(glue("{sim_i$mesh[i]}, {sim_i$timeRes[i]}, {sim_i$liceSpeedF[i]}"), 
#                    ". {closest_state}"))
#   anim_save(glue("figs/tracks_sim{str_pad(i,2,'left','0')}_10kPart.gif"), 
#             anim, nframes=interp*length(hours)*length(days),
#             fps=32, width=7, height=6, res=300, units="in")
#   paste("Finished", i)
# }
# stopCluster(cl)
# 
# 
# 
# 
# vel.ls <- velocity.df %>%
#   group_by(meshRes, depth) %>%
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
#             ggtitle(paste0(glue("{info$meshRes}, {info$liceSpeedF}, {info$depth}m"),
#                            ". {closest_state}"))
#           anim_save(glue("figs/tracks_dir_{info$meshRes}_{info$depth}m.gif"),
#                     anim, nframes=interp*length(hours)*length(days),
#                     fps=32, width=7, height=6, res=300, units="in")
#           paste("Finished", i)
#         }
# stopCluster(cl)

vel.ls <- velocity.df %>%
  group_by(depth) %>%
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
          
          anim <- vel.ls[[i]] %>%
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
          anim_save(glue("figs/tracks_dir_{info$depth}m.gif"),
                    anim, nframes=interp*length(hours)*length(days),
                    fps=32, width=9, height=4, res=300, units="in")
          paste("Finished", i)
        }
stopCluster(cl)



# TODO: In QGIS, self-NNJoin each layer, calculated abs(salinity - join_salinity)
# then set threshold + buffer to define tidal pulse. Might work... Not sure how 
# to automate in QGIS, but there could be a way in sf

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
#           anim_save(glue("figs/tracks_dir_{info$meshRes}.gif"), 
#                     anim, nframes=interp*length(hours)*length(days),
#                     fps=32, width=7, height=6, res=300, units="in")
#           paste("Finished", i)
#         }
# stopCluster(cl)
