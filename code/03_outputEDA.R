# Output EDA
# Off Aqua
# Tim Szewczyk


# This script does some exploratory analysis on the output


# setup -------------------------------------------------------------------

library(tidyverse); library(glue); library(lubridate); library(sf)
source("code/00_fn.R")

dirs <- switch(get_os(),
               windows=list(proj=getwd(),
                            mesh="D:/hydroOut/",
                            out=glue("{getwd()}/out/")),
               linux=list(proj=getwd(),
                          mesh="/home/sa04ts/FVCOM_meshes",
                          out=glue("{getwd()}/out/")))




# load files --------------------------------------------------------------

mesh.sf <- list(linnhe7=st_read(glue("{dirs$mesh}/linnhe_mesh.gpkg")) %>%
                  mutate(mesh="linnhe7"),
                westcoms2=st_read(glue("{dirs$mesh}/WeStCOMS2_linnhe_mesh.gpkg")) %>%
                  mutate(mesh="WeStCOMS2")) %>%
  do.call('rbind', .) %>%
  select(i, area, depth, mesh, lochRegion, geom) %>%
  full_join(., read_csv("data/loch_regions.csv"), by="lochRegion")
mesh.fp <- st_read("data/linnhe_mesh_footprint.gpkg")
sim_i <- read_csv(glue("{dirs$out}/sim_i.csv")) %>%
  mutate(liceSpeedF=factor(liceSpeed, labels=c("Slow", "Medium", "Fast")),
         sim=as.numeric(i))
elemAct.sf <- st_read("out/00_processed/elementActivity.gpkg") %>%
  mutate(liceSpeedF=factor(liceSpeedF, levels=levels(sim_i$liceSpeedF)),
         meshRes=factor(paste0(mesh, ", ", timeRes),
                        levels=c("linnhe7, 1h", "linnhe7, 5min",
                                 "WeStCOMS2, 1h", "WeStCOMS2, 5min")))
loc.df <- readRDS("out/00_processed/locations.rds") %>%
  mutate(liceSpeedF=factor(liceSpeedF, levels=levels(sim_i$liceSpeedF)),
         meshRes=factor(paste0(mesh, ", ", timeRes),
                        levels=c("linnhe7, 1h", "linnhe7, 5min",
                                 "WeStCOMS2, 1h", "WeStCOMS2, 5min")),
         hour=hour(timeCalculated),
         date=date(timeCalculated))
# loc.sf <- st_as_sf(loc.df, coords=c("x", "y"), crs=27700)



# plot elements -----------------------------------------------------------

theme_set(theme_classic())
cmr <- readRDS("../../00_misc/cmr_cmaps.RDS")
col_hour <- scale_colour_gradientn(colours=cmr$seasons, limits=c(0,24),
                                   breaks=c(0,6,12,18,24))
fill_NSEW <- scale_fill_gradientn(colours=cmr$infinity, limits=c(-pi, pi),
                                  breaks=c(-pi, -pi/2, 0, pi/2, pi),
                                  labels=c("W", "S", "E", "N", "W"))
col_NSEW <- scale_colour_gradientn(colours=cmr$infinity, limits=c(-pi, pi),
                                   breaks=c(-pi, -pi/2, 0, pi/2, pi),
                                   labels=c("W", "S", "E", "N", "W"))
col_downloch <- scale_colour_gradient2("Direction", limits=c(-1,1), 
                                       breaks=c(-1,1), 
                                       labels=c("uploch", "downloch"))
fill_downloch <- scale_fill_gradient2("Direction", midpoint=0.5, 
                                      breaks=c(0, 1), labels=c("downloch", "uploch"),
                                      low=scales::muted("blue"), 
                                      mid="grey90",
                                      high=scales::muted("red"))



# mesh plots --------------------------------------------------------------

ggplot(mesh.sf, aes(fill=bearing)) + 
  geom_sf(colour=NA) + 
  fill_NSEW + 
  facet_grid(.~mesh)
ggsave("figs/mesh_bearing.png", width=6, height=4, dpi=300)
ggplot(mesh.sf, aes(fill=depth)) + 
  geom_sf(colour=NA) + 
  scale_fill_viridis_c() +
  facet_grid(.~mesh)
ggsave("figs/mesh_depth.png", width=6, height=4, dpi=300)



# vertical distribution ---------------------------------------------------

# tendency to be deeper in WeStCOMS, but need to filter by bathymetric depth
# since linnhe includes shallower areas
loc.df %>% filter(timeCalculated > "2021-11-01 ", status != 66) %>%
  ggplot(aes(depth, colour=meshRes)) +
  geom_density() +
  scale_colour_brewer(type="qual", palette=2) +
  xlim(0, 25) +
  facet_grid(liceSpeedF~.)
ggsave("figs/vertDist.png", width=5, height=6, dpi=300)

loc.df %>% filter(age>12 & status != 66) %>%
  ggplot(aes(depth, colour=liceSpeedF)) +
  geom_density() +
  scale_colour_viridis_d() +
  xlim(0, 25) +
  facet_grid(meshRes~.)

# By hour
loc.df %>% filter(age>12 & status != 66) %>%
  ggplot(aes(depth, colour=hour, group=hour)) +
  geom_density() +
  xlim(0, 25) +
  col_hour +
  facet_grid(meshRes~liceSpeedF)
ggsave("figs/vertDist_by_hour.png", width=9, height=4, dpi=300)

loc.df %>% filter(age>12 & status != 66) %>%
  ggplot(aes(depth, colour=meshRes)) +
  geom_density() +
  xlim(0, 25) +
  facet_grid(liceSpeedF~hour)
ggsave("figs/vertDist_by_hour2.png", width=15, height=6, dpi=300)

loc.df %>% filter(age>12 & status != 66) %>%
  ggplot(aes(depth, linetype=meshRes, colour=liceSpeedF, group=sim)) +
  geom_density() +
  xlim(0, 25) +
  scale_colour_viridis_d(end=0.9) +
  facet_wrap(~hour, ncol=6)
ggsave("figs/vertDist_by_hour3.png", width=12, height=9, dpi=300)





# tracks ------------------------------------------------------------------

part.samp <- sample(unique(loc.df$ID), 500)
loc.df %>% filter(ID %in% part.samp) %>%
  ggplot(aes(x, y, group=ID)) +
  geom_path(alpha=0.25) + 
  facet_grid(meshRes~liceSpeedF)






# mesh exits --------------------------------------------------------------

# lice more likely to exit loch in WeStCOMS, at faster sink/swim speeds
loc.df %>% 
  group_by(meshRes, liceSpeedF) %>% 
  summarise(propExit=sum(status==66)/n_distinct(ID)) %>% 
  ggplot(aes(liceSpeedF, propExit, colour=meshRes, group=meshRes)) + 
  geom_point() + geom_path()
ggsave("figs/exits.png", width=5, height=3, dpi=300)

# confirming exits occur at open boundaries only
loc.df %>%
  filter(status==66) %>%
  ggplot(aes(x, y, colour=age)) + 
  geom_point(alpha=0.2, size=0.5, shape=1) +
  scale_colour_viridis_c() +
  facet_grid(meshRes~liceSpeedF)

# Not sure...
loc.df %>%
  filter(status==66) %>%
  ggplot(aes(age, colour=liceSpeedF)) + 
  geom_density() +
  scale_colour_viridis_d() +
  facet_grid(meshRes~.)

# Not sure... lice exit more slowly at an older age at high res
loc.df %>%
  filter(status==66) %>%
  ggplot(aes(age, colour=meshRes)) + 
  geom_density() +
  facet_grid(liceSpeedF~.)




# total distance ----------------------------------------------------------

loc.df %>% 
  filter(age>12 & status != 66) %>%
  group_by(ID, sim, meshRes, liceSpeedF) %>%
  slice_tail(n=1) %>%
  ggplot(aes(xyTot/age, colour=meshRes, linetype=liceSpeedF, group=sim)) + 
  geom_density() + 
  scale_colour_brewer(type="qual") +
  scale_x_log10("Mean xy hourly displacement (m/h)")
ggsave("figs/xy_mean_log.png", width=6, height=4, dpi=300)
loc.df %>% 
  filter(status!=66) %>%
  group_by(ID, sim, meshRes, liceSpeedF) %>%
  slice_tail(n=1) %>%
  ggplot(aes(xyTot/age, colour=meshRes, linetype=liceSpeedF, group=sim)) + 
  geom_density() + 
  scale_colour_brewer(type="qual") +
  scale_x_continuous("Mean xy hourly displacement (m/h)")
ggsave("figs/xy_mean.png", width=6, height=4, dpi=300)

loc.df %>%
  group_by(ID, sim, meshRes, liceSpeedF) %>%
  arrange(age) %>%
  mutate(dZ=depth-lag(depth)) %>%
  filter(age>12 & status != 66) %>%
  ggplot(aes(abs(dZ), colour=meshRes, linetype=liceSpeedF, group=sim)) + 
  geom_density() + 
  scale_colour_brewer(type="qual") +
  scale_x_continuous("Hourly depth displacement (abs(m))", trans="log1p")
ggsave("figs/z_mean_log.png", width=6, height=4, dpi=300)
loc.df %>%
  group_by(ID, sim, meshRes, liceSpeedF) %>%
  arrange(age) %>%
  mutate(dZ=depth-lag(depth)) %>%
  filter(age>12 & status != 66) %>%
  ggplot(aes(dZ, colour=meshRes, linetype=liceSpeedF, group=sim)) + 
  geom_density() + 
  scale_colour_brewer(type="qual") +
  scale_x_continuous("Hourly depth displacement (m)")
ggsave("figs/z_mean.png", width=6, height=4, dpi=300)



# element activity --------------------------------------------------------

p <- ggplot(elemAct.sf, aes(fill=sink/total)) + geom_sf(colour=NA) +
  scale_fill_viridis_c(limits=c(0,1)) +
  facet_grid(meshRes~liceSpeed)
ggsave("figs/prSink_mesh_by_speed.png", p, width=8, height=10, dpi=300)

p <- ggplot(elemAct.sf, aes(fill=swim/total)) + geom_sf(colour=NA) +
  scale_fill_viridis_c(limits=c(0,1)) +
  facet_grid(meshRes~liceSpeed)
ggsave("figs/prSwim_mesh_by_speed.png", p, width=8, height=10, dpi=300)

p <- ggplot(elemAct.sf, aes(fill=float/total)) + geom_sf(colour=NA) +
  scale_fill_viridis_c(limits=c(0,1)) +
  facet_grid(meshRes~liceSpeed)
ggsave("figs/prFloat_mesh_by_speed.png", p, width=8, height=10, dpi=300)







# maps --------------------------------------------------------------------

for(i in 1:nrow(sim_i)) {
  p <- loc.df %>%
    filter(age>12 & status != 66) %>%
    filter(sim==i) %>%
    ggplot(aes(x,y, alpha=density, colour=depth)) + 
    geom_point(size=0.1, shape=1) + 
    scale_colour_viridis_c(direction=-1, option="D", limits=c(0,20), na.value="#000033") +
    ggtitle(glue("{sim_i$mesh[i]}, {sim_i$timeRes[i]}, {sim_i$liceSpeedF[i]}"))
  ggsave(glue("figs/part_locs_sim{str_pad(i,2,'left','0')}.png"), 
         p, width=6, height=6, dpi=900)
}







# velocities --------------------------------------------------------------

velocity.df <- loc.df %>%
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
  filter(age>12 & status != 66) %>%
  left_join(., sim_i %>% select(sim, mesh, liceSpeedF)) %>%
  left_join(., mesh.sf %>% st_drop_geometry() %>% 
              select(i, mesh, bearing, depth) %>% 
              rename(elem=i, lochDir=bearing, meshDepth=depth),
            by=c("mesh", "elem")) %>%
  mutate(downloch=cos(lochDir-bearing)) # -1: uploch, 1: downloch


velocity.df %>%
  filter(ID %% 100 == 0) %>%
  ggplot(aes(x, y, colour=bearing)) +
  geom_point(size=0.5) +
  col_NSEW +
  facet_grid(meshRes~liceSpeedF)

velocity.df %>%
  filter(ID %% 100 == 0) %>%
  ggplot(aes(x, y, colour=downloch)) +
  geom_point() +
  col_downloch +
  facet_wrap(~hour(timeCalculated), ncol=6)

velocity.df %>%
  mutate(depthCat=if_else(depth < 20, round(depth), 20)) %>%
  ggplot(aes(downloch, colour=depthCat, group=depthCat)) +
  geom_density() + 
  scale_colour_viridis_c(direction=-1) +
  facet_grid(~sim)


# compare across depths: absolute, or relative to meshDepth
velocity.df %>%
  mutate(depthCat=if_else(depth < 20, round(depth), 20)) %>%
  ggplot(aes(downloch, colour=depthCat, group=depthCat)) +
  geom_density() + 
  scale_colour_viridis_c(direction=-1) +
  facet_grid(~sim)
velocity.df %>%
  mutate(depthCat=if_else(depth<=5, "<= 5m", "> 5m")) %>%
  ggplot(aes(downloch, colour=depthCat, group=depthCat)) +
  geom_density() + 
  facet_grid(~sim)
velocity.df %>%
  mutate(depthCat=if_else(depth<=5, "<= 5m", "> 5m")) %>%
  ggplot(aes(downloch, colour=hour(timeCalculated), 
             group=hour(timeCalculated))) +
  geom_density() + 
  col_hour +
  facet_grid(depthCat~sim)


# p <- velocity.df %>%
#   filter(sim==4) %>%
#   filter(ID %% 20 == 0) %>%
#   mutate(depthCat=if_else(depth_m1<=5, "<= 5m", "> 5m")) %>%
#   ggplot() +
#   geom_sf(data=mesh.fp, fill=NA, colour="grey30") +
#   geom_path(aes(x, y, group=ID, colour=downloch), alpha=0.25, size=0.1) + 
#   col_downloch +
#   facet_wrap(~depthCat) +
#   ggtitle("High resolution, 1h, medium lice speed")
# ggsave("figs/direction_map_sim4.png", p, width=8, height=5, dpi=300)
# 
# p <- velocity.df %>%
#   filter(sim==3) %>%
#   filter(ID %% 20 == 0) %>%
#   mutate(depthCat=if_else(depth_m1<=5, "<= 5m", "> 5m")) %>%
#   ggplot() +
#   geom_sf(data=mesh.fp, fill=NA, colour="grey30") +
#   geom_path(aes(x, y, group=ID, colour=downloch), alpha=0.25, size=0.1) + 
#   col_downloch +
#   facet_wrap(~depthCat) +
#   ggtitle("High resolution, 1h, passive lice")
# ggsave("figs/direction_map_sim3.png", p, width=8, height=5, dpi=300)

p <- velocity.df %>%
  filter(ID %% 200 == 0) %>%
  ggplot() +
  geom_sf(data=mesh.fp, fill=NA, colour="grey30") +
  geom_path(aes(x, y, group=ID, colour=downloch), alpha=0.25, size=0.1) + 
  col_downloch +
  facet_grid(meshRes~liceSpeedF)
ggsave("figs/direction_map.png", p, width=8, height=10, dpi=300)



meshDir.sum <- left_join(mesh.sf,
                         velocity.df %>%
                           filter(liceSpeedF=="Medium") %>%
                           group_by(mesh, meshRes, elem) %>%
                           summarise(N=n(),
                                     totDens=sum(density),
                                     prUploch=mean(downloch<0)) %>%
                           ungroup %>% rename(i=elem),
                         by=c("mesh", "i")) %>%
  mutate(N=replace_na(N, 0),
         lnDens_m3=log(totDens/(area*depth))) %>%
  filter(!is.na(meshRes))
p <- ggplot(meshDir.sum) +
  geom_sf(data=mesh.fp, colour="grey30", fill=NA) +
  geom_sf(colour=NA, aes(fill=lnDens_m3)) + 
  scale_fill_viridis_c() +
  facet_grid(.~meshRes)
ggsave("figs/density_map.png", p, width=11, height=5, dpi=300)
p <- ggplot(meshDir.sum %>% filter(N>10)) +
  geom_sf(data=mesh.fp, colour="grey30", fill=NA) +
  geom_sf(colour=NA, aes(fill=prUploch)) + 
  fill_downloch +
  facet_grid(.~meshRes)
ggsave("figs/prUploch_map.png", p, width=11, height=5, dpi=300)
p <- ggplot(meshDir.sum %>% filter(N>10)) +
  geom_sf(data=mesh.fp, colour="grey30", fill=NA) +
  geom_sf(colour=NA, aes(fill=prUploch, alpha=lnDens_m3)) + 
  fill_downloch +
  facet_grid(.~meshRes)
ggsave("figs/prUploch_map_ALT.png", p, width=11, height=5, dpi=300)


meshDir.deep <- left_join(mesh.sf,
                          velocity.df %>%
                            filter(liceSpeedF=="Medium") %>%
                            mutate(deep=case_when(depth_m1 <= 2 ~ "0-2",
                                                  between(depth_m1, 2, 5) ~ "2-5",
                                                  between(depth_m1, 5, 10) ~ "5-10",
                                                  depth_m1 > 10 ~ ">10")) %>%
                            group_by(mesh, meshRes, elem, deep) %>%
                            summarise(N=n(),
                                      totDens=sum(density),
                                      prUploch=mean(downloch<0)) %>%
                            ungroup %>% rename(i=elem),
                          by=c("mesh", "i")) %>%
  mutate(N=replace_na(N, 0),
         lnDens_m3=case_when(deep=="0-2" ~ log(totDens/(area*pmin(depth, 2))),
                             deep=="2-5" ~ log(totDens/(area*pmin(depth-2, 3))),
                             deep=="5-10" ~ log(totDens/(area*pmin(depth-5, 5))),
                             deep==">10" ~ log(totDens/(area*pmin(depth-10, 1)))),
         deep=factor(deep, levels=c("0-2", "2-5", "5-10", ">10"))) %>%
  filter(!is.na(meshRes))

p <- ggplot(meshDir.deep %>% filter(!is.na(deep))) +
  geom_sf(data=mesh.fp, colour="grey30", fill=NA) +
  geom_sf(colour=NA, aes(fill=lnDens_m3)) + 
  scale_fill_viridis_c() +
  facet_grid(deep~meshRes)
ggsave("figs/density_map_depthStrat.png", p, width=11, height=14, dpi=300)
p <- ggplot(meshDir.deep %>% filter(!is.na(deep), N>10)) +
  geom_sf(data=mesh.fp, colour="grey30", fill=NA) +
  geom_sf(colour=NA, aes(fill=prUploch)) + 
  fill_downloch +
  facet_grid(deep~meshRes)
ggsave("figs/prUploch_map_depthStrat.png", p, width=11, height=14, dpi=300)
p <- ggplot(meshDir.deep %>% filter(!is.na(deep), N>10)) +
  geom_sf(data=mesh.fp, colour="grey30", fill=NA) +
  geom_sf(colour=NA, aes(fill=prUploch, alpha=lnDens_m3)) + 
  fill_downloch +
  facet_grid(deep~meshRes)
ggsave("figs/prUploch_map_depthStrat_ALT.png", p, width=11, height=14, dpi=300)
