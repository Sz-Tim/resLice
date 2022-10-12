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
                  mutate(mesh="WeStCOMS2", lochRegion=NA)) %>%
  do.call('rbind', .) %>%
  select(i, area, depth, mesh, lochRegion, geom)
mesh.fp <- st_read("data/linnhe_mesh_footprint.gpkg")
sim_i <- read_csv(glue("{dirs$out}/sim_i.csv")) %>%
  mutate(liceSpeedF=factor(liceSpeed, labels=c("Passive", "Medium")))
elemAct.sf <- st_read("out/00_processed/elementActivity.gpkg") %>%
  mutate(liceSpeedF=factor(liceSpeedF, levels=levels(sim_i$liceSpeedF)))
loc.df <- read_csv("out/00_processed/locations.csv") %>%
  mutate(liceSpeedF=factor(liceSpeedF, levels=levels(sim_i$liceSpeedF)))
loc.sf <- st_as_sf(loc.df, coords=c("x", "y"), crs=27700)



# plot elements -----------------------------------------------------------

theme_set(theme_classic())
cmr <- readRDS("../../00_misc/cmr_cmaps.RDS")
col_hour <- scale_colour_gradientn(colours=cmr$seasons, limits=c(0,24),
                                   breaks=c(0,6,12,18,24))




# vertical distribution ---------------------------------------------------

# tendency to be deeper in WeStCOMS, but need to filter by bathymetric depth
# since linnhe includes shallower areas
loc.df %>% filter(status==2) %>%
  ggplot(aes(depth, colour=mesh)) +
  geom_density() +
  xlim(0, 25) +
  facet_grid(liceSpeedF~.)

loc.df %>% filter(status==2) %>%
  ggplot(aes(depth, colour=liceSpeedF)) +
  geom_density() +
  scale_colour_viridis_d() +
  xlim(0, 25) +
  facet_grid(mesh~.)

# By hour
loc.df %>% filter(status==2) %>%
  ggplot(aes(depth, colour=hour, group=hour)) +
  geom_density() +
  xlim(0, 25) +
  col_hour +
  facet_grid(mesh~liceSpeedF)
ggsave("figs/vertDist_by_hour.png", width=9, height=4, dpi=300)

loc.df %>% filter(status==2) %>%
  ggplot(aes(depth, colour=mesh)) +
  geom_density() +
  xlim(0, 25) +
  facet_grid(liceSpeedF~hour)
ggsave("figs/vertDist_by_hour2.png", width=15, height=6, dpi=300)

loc.df %>% filter(status==2) %>%
  ggplot(aes(depth, linetype=mesh, colour=liceSpeedF, group=sim)) +
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
  facet_grid(mesh~liceSpeedF)






# mesh exits --------------------------------------------------------------

# lice more likely to exit loch in WeStCOMS, at faster sink/swim speeds
loc.df %>% 
  group_by(mesh, liceSpeed) %>% 
  summarise(propExit=sum(status==66)/n_distinct(ID)) %>% 
  ggplot(aes(liceSpeed, propExit, colour=mesh)) + 
  geom_point() + geom_line()
ggsave("figs/exits.png", width=5, height=3, dpi=300)

# confirming exits occur at open boundaries only
loc.df %>%
  filter(status==66) %>%
  ggplot(aes(x, y, colour=age)) + 
  geom_point(alpha=0.2, size=0.5, shape=1) +
  scale_colour_viridis_c() +
  facet_grid(mesh~liceSpeedF)

# Not sure...
loc.df %>%
  filter(status==66, startDate=="20211101") %>%
  ggplot(aes(age, colour=liceSpeedF)) + 
  geom_density() +
  scale_colour_viridis_d() +
  facet_grid(mesh~.)

# Not sure... lice exit more slowly at an older age at high res
loc.df %>%
  filter(status==66, startDate=="20211101") %>%
  ggplot(aes(age, colour=mesh)) + 
  geom_density() +
  facet_grid(liceSpeedF~.)




# total distance ----------------------------------------------------------

loc.df %>% 
  group_by(ID, sim, mesh, liceSpeedF) %>%
  arrange(age) %>%
  mutate(xy=sqrt((x-lag(x))^2 + (y-lag(y))^2)) %>%
  filter(!is.na(xy)) %>%
  mutate(xy_cumul=cumsum(xy)) %>%
  filter(status==2) %>%
  slice_tail(n=1) %>%
  ggplot(aes(xy_cumul/age, linetype=mesh, colour=liceSpeedF, group=sim)) + 
  geom_density() + 
  scale_colour_viridis_d() +
  scale_x_log10("Mean xy hourly displacement (m/h)")
ggsave("figs/xy_mean_log.png", width=6, height=4, dpi=300)
loc.df %>% 
  group_by(ID, sim, mesh, liceSpeedF) %>%
  arrange(age) %>%
  mutate(xy=sqrt((x-lag(x))^2 + (y-lag(y))^2)) %>%
  filter(!is.na(xy)) %>%
  mutate(xy_cumul=cumsum(xy)) %>%
  filter(status==2) %>%
  slice_tail(n=1) %>%
  ggplot(aes(xy_cumul/age, linetype=mesh, colour=liceSpeedF, group=sim)) + 
  geom_density() + 
  scale_colour_viridis_d() +
  scale_x_continuous("Mean xy hourly displacement (m/h)")
ggsave("figs/xy_mean.png", width=6, height=4, dpi=300)

loc.df %>%
  group_by(ID, sim, mesh, liceSpeedF) %>%
  arrange(age) %>%
  filter(status==2) %>%
  slice_tail(n=1) %>%
  ggplot(aes(zTot/age, linetype=mesh, colour=liceSpeedF, group=sim)) + 
  geom_density() + 
  scale_colour_viridis_d() +
  scale_x_log10("Mean z hourly displacement (m)")
ggsave("figs/z_mean_log.png", width=6, height=4, dpi=300)
loc.df %>%
  group_by(ID, sim, mesh, liceSpeedF) %>%
  arrange(age) %>%
  filter(status==2) %>%
  slice_tail(n=1) %>%
  ggplot(aes(zTot/age, linetype=mesh, colour=liceSpeedF, group=sim)) + 
  geom_density() + 
  scale_colour_viridis_d() +
  scale_x_continuous("Mean z hourly displacement (m)")
ggsave("figs/z_mean.png", width=6, height=4, dpi=300)



# element activity --------------------------------------------------------

p <- ggplot(elemAct.sf, aes(fill=sink/total)) + geom_sf(colour=NA) +
  scale_fill_viridis_c(limits=c(0,1)) +
  facet_grid(mesh~liceSpeed)
ggsave("figs/prSink_mesh_by_speed.png", p, width=10, height=5, dpi=300)

p <- ggplot(elemAct.sf, aes(fill=swim/total)) + geom_sf(colour=NA) +
  scale_fill_viridis_c(limits=c(0,1)) +
  facet_grid(mesh~liceSpeed)
ggsave("figs/prSwim_mesh_by_speed.png", p, width=10, height=5, dpi=300)

p <- ggplot(elemAct.sf, aes(fill=float/total)) + geom_sf(colour=NA) +
  scale_fill_viridis_c(limits=c(0,1)) +
  facet_grid(mesh~liceSpeed)
ggsave("figs/prFloat_mesh_by_speed.png", p, width=10, height=5, dpi=300)







# maps --------------------------------------------------------------------

for(i in 1:nrow(sim_i)) {
  p <- loc.df %>%
    filter(status==2) %>%
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
  select(ID, dateTime, age, status, x, y, sim, mesh, elem, density, depth) %>%
  arrange(sim, dateTime, ID) %>%
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
  filter(status==2)


velocity.df %>%
  filter(sim==2) %>%
  filter(ID %% 100 == 0) %>%
  ggplot(aes(x, y, colour=bearing)) +
  geom_point(size=0.5) +
  scale_colour_gradientn(colours=cmr$infinity, breaks=c(-pi, -pi/2, 0, pi/2, pi),
                         labels=c("W", "S", "E", "N", "W"), limits=c(-pi, pi))

# uploch (from mouth to head, but confirm) bearing in radians
lochDir.df <- st_read("d:/hydroOut/linnhe_skeleton.shp") %>%
  mutate(x=st_coordinates(.)[,1],
         y=st_coordinates(.)[,2]) %>%
  st_drop_geometry() %>%
  arrange(lochRegion, lochHead) %>%
  group_by(lochRegion) %>%
  summarise(bearing=atan2(first(y)-last(y), first(x)-last(x)))
mesh.sf <- mesh.sf %>%
  full_join(lochDir.df, by="lochRegion")
velocity.df <- velocity.df %>%
  filter(mesh=="linnhe7") %>%
  left_join(., mesh.sf %>% st_drop_geometry() %>% 
              select(i, mesh, bearing, depth) %>% 
              rename(elem=i, lochDir=bearing, meshDepth=depth),
            by=c("mesh", "elem")) %>%
  mutate(downloch=cos(lochDir-bearing)) # -1: uploch, 1: downloch
# cos(lochDir - partDir)
# positive = towards same direction (aligned at 1)
# negative = opposite direction (aligned at -1)

velocity.df %>%
  filter(sim==4) %>%
  filter(ID %% 100 == 0) %>%
  ggplot(aes(x, y, colour=downloch)) +
  geom_point() +
  scale_colour_gradient2("Direction", limits=c(-1,1), 
                         breaks=c(-1,1), labels=c("uploch", "downloch")) +
  facet_wrap(~hour(dateTime), ncol=6)

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
  ggplot(aes(downloch, colour=hour(dateTime), group=hour(dateTime))) +
  geom_density() + 
  scale_colour_gradientn(colours=cmr$seasons) +
  facet_grid(depthCat~sim)


p <- velocity.df %>%
  filter(sim==4) %>%
  filter(ID %% 20 == 0) %>%
  mutate(depthCat=if_else(depth<=5, "<= 5m", "> 5m")) %>%
  ggplot() +
  geom_sf(data=mesh.fp, fill=NA, colour="grey30") +
  geom_path(aes(x, y, group=ID, colour=downloch), alpha=0.25, size=0.1) + 
  scale_colour_gradient2("Direction", limits=c(-1,1), 
                         breaks=c(-1,1), labels=c("uploch", "downloch")) +
  facet_wrap(~depthCat) +
  ggtitle("High resolution, 1h, medium lice speed")
ggsave("figs/direction_map_sim4.png", p, width=8, height=5, dpi=300)

p <- velocity.df %>%
  filter(sim==3) %>%
  filter(ID %% 20 == 0) %>%
  mutate(depthCat=if_else(depth<=5, "<= 5m", "> 5m")) %>%
  ggplot() +
  geom_sf(data=mesh.fp, fill=NA, colour="grey30") +
  geom_path(aes(x, y, group=ID, colour=downloch), alpha=0.25, size=0.1) + 
  scale_colour_gradient2("Direction", limits=c(-1,1), 
                         breaks=c(-1,1), labels=c("uploch", "downloch")) +
  facet_wrap(~depthCat) +
  ggtitle("High resolution, 1h, passive lice")
ggsave("figs/direction_map_sim3.png", p, width=8, height=5, dpi=300)

p <- velocity.df %>%
  filter(ID %% 20 == 0) %>%
  ggplot() +
  geom_sf(data=mesh.fp, fill=NA, colour="grey30") +
  geom_path(aes(x, y, group=ID, colour=downloch), alpha=0.25, size=0.1) + 
  scale_colour_gradient2("Direction", limits=c(-1,1), 
                         breaks=c(-1,1), labels=c("uploch", "downloch")) +
  facet_wrap(~sim) +
  ggtitle("High resolution, 1h, passive (3) or medium (4) lice speed")
ggsave("figs/direction_map.png", p, width=8, height=5, dpi=300)
