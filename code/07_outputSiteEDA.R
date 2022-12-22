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
                            out=glue("{getwd()}/out/siteRelease/")),
               linux=list(proj=getwd(),
                          mesh="/home/sa04ts/FVCOM_meshes",
                          out=glue("{getwd()}/out/siteRelease/")))




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
  mutate(liceSpeedF=factor(liceSpeed, levels=c(0.0001, 0.0005, 0.001), 
                           labels=c("Slow", "Medium", "Fast")),
         sim=as.numeric(i))
elemAct.sf <- st_read("out/00_processed/elementActivity_site.gpkg") %>%
  mutate(liceSpeedF=factor(liceSpeedF, levels=levels(sim_i$liceSpeedF)),
         meshRes=factor(paste0(mesh, ", ", timeRes),
                        levels=c("WeStCOMS2, 1h", "linnhe7, 1h", "linnhe7, 5min", "WeStCOMS2, 5min")))
loc.df <- readRDS("out/00_processed/locations_site.rds") %>%
  mutate(liceSpeedF=factor(liceSpeedF, levels=levels(sim_i$liceSpeedF)),
         meshRes=factor(paste0(mesh, ", ", timeRes),
                        levels=c("WeStCOMS2, 1h", "linnhe7, 1h", "linnhe7, 5min", "WeStCOMS2, 5min")),
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
ggsave("figs/mesh_bearing_site.png", width=6, height=4, dpi=300)
ggplot(mesh.sf, aes(fill=depth)) + 
  geom_sf(colour=NA) + 
  scale_fill_viridis_c() +
  facet_grid(.~mesh)
ggsave("figs/mesh_depth_site.png", width=6, height=4, dpi=300)



# vertical distribution ---------------------------------------------------

# tendency to be deeper in WeStCOMS, but need to filter by bathymetric depth
# since linnhe includes shallower areas
loc.df %>% filter(age>=12, status != 66) %>%
  ggplot(aes(depth, colour=meshRes)) +
  geom_density() +
  scale_colour_brewer(type="div") +
  xlim(0, 25) +
  facet_grid(liceSpeedF~.)
ggsave("figs/vertDist_site.png", width=5, height=6, dpi=300)

loc.df %>% filter(age>=12 & status != 66) %>%
  ggplot(aes(depth, colour=liceSpeedF)) +
  geom_density() +
  scale_colour_viridis_d() +
  xlim(0, 25) +
  facet_grid(meshRes~.)

# By hour
loc.df %>% filter(age>=12 & status != 66) %>%
  ggplot(aes(depth, colour=hour, group=hour)) +
  geom_density() +
  xlim(0, 25) +
  col_hour +
  facet_grid(meshRes~liceSpeedF)
ggsave("figs/vertDist_by_hour_site.png", width=8, height=6, dpi=300)

loc.df %>% filter(age>=12 & status != 66) %>%
  ggplot(aes(depth, colour=meshRes)) +
  geom_density() +
  xlim(0, 25) +
  scale_colour_brewer(type="div") +
  facet_grid(liceSpeedF~hour)
ggsave("figs/vertDist_by_hour2_site.png", width=17, height=5, dpi=300)

loc.df %>% filter(age>=12 & status != 66) %>%
  ggplot(aes(depth, linetype=meshRes, colour=liceSpeedF, group=sim)) +
  geom_density() +
  xlim(0, 25) +
  scale_colour_viridis_d(end=0.9) +
  facet_wrap(~hour, ncol=6)
ggsave("figs/vertDist_by_hour3_site.png", width=12, height=9, dpi=300)





# tracks ------------------------------------------------------------------

# bbox <- st_bbox(mesh.fp)
# part.samp <- sample(unique(loc.df$ID), min(length(unique(loc.df$ID)), 1e3))
# loc.df %>% filter(ID %in% part.samp) %>%
#   filter(between(x, bbox$xmin, bbox$xmax), between(y, bbox$ymin, bbox$ymax)) %>%
#   ggplot(aes(x, y, group=ID)) +
#   geom_path(alpha=0.25) + 
#   facet_grid(meshRes~liceSpeedF)






# mesh exits --------------------------------------------------------------

# lice more likely to exit loch in WeStCOMS, at faster sink/swim speeds
loc.df %>% 
  group_by(meshRes, liceSpeedF, ID) %>% 
  summarise(exited=any(meshParticle==1)) %>% 
  group_by(meshRes, liceSpeedF) %>%
  summarise(propExit=mean(exited)) %>% 
  ggplot(aes(liceSpeedF, propExit, colour=meshRes, group=meshRes)) + 
  geom_point() + geom_path() +
  scale_colour_brewer(type="div") +
  ylim(0, NA) + labs(y="Proportion leaving Loch Linnhe")
ggsave("figs/exits_site.png", width=5, height=3, dpi=300)

# confirming exits occur at open boundaries only
loc.df %>% 
  filter(meshParticle==1) %>%
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
  scale_colour_brewer(type="div") + 
  geom_density() +
  facet_grid(liceSpeedF~.)




# total distance ----------------------------------------------------------

loc.df %>% 
  filter(age>=12 & status != 66) %>%
  group_by(ID, sim, meshRes, liceSpeedF) %>%
  slice_tail(n=1) %>%
  ggplot(aes(xyTot/age, colour=meshRes, linetype=liceSpeedF, group=sim)) + 
  geom_density() + 
  scale_colour_brewer(type="div") +
  scale_x_log10("Mean xy hourly displacement (m/h)")
ggsave("figs/xy_mean_log_site.png", width=6, height=4, dpi=300)
loc.df %>% 
  filter(status!=66) %>%
  group_by(ID, sim, meshRes, liceSpeedF) %>%
  slice_tail(n=1) %>%
  ggplot(aes(xyTot/age, colour=meshRes, linetype=liceSpeedF, group=sim)) + 
  geom_density() + 
  scale_colour_brewer(type="div") +
  scale_x_continuous("Mean xy hourly displacement (m/h)")
ggsave("figs/xy_mean_site.png", width=6, height=4, dpi=300)

loc.df %>%
  group_by(ID, sim, meshRes, liceSpeedF) %>%
  arrange(age) %>%
  mutate(dZ=depth-lag(depth)) %>%
  filter(age>=12 & status != 66) %>%
  ggplot(aes(abs(dZ), colour=meshRes, linetype=liceSpeedF, group=sim)) + 
  geom_density() + 
  scale_colour_brewer(type="div") +
  scale_x_continuous("Hourly depth displacement (abs(m))", trans="log1p")
ggsave("figs/z_mean_log_site.png", width=6, height=4, dpi=300)
loc.df %>%
  group_by(ID, sim, meshRes, liceSpeedF) %>%
  arrange(age) %>%
  mutate(dZ=depth-lag(depth)) %>%
  filter(age>=12 & status != 66) %>%
  ggplot(aes(dZ, colour=meshRes, linetype=liceSpeedF, group=sim)) + 
  geom_density() + 
  scale_colour_brewer(type="div") +
  scale_x_continuous("Hourly depth displacement (m)")
ggsave("figs/z_mean_site.png", width=6, height=4, dpi=300)



# element activity --------------------------------------------------------

p <- ggplot(elemAct.sf, aes(fill=sink/total)) + geom_sf(colour=NA) +
  scale_fill_viridis_c(limits=c(0,1)) +
  facet_grid(meshRes~liceSpeed)
ggsave("figs/prSink_mesh_by_speed_site.png", p, width=8, height=8, dpi=300)

p <- ggplot(elemAct.sf, aes(fill=swim/total)) + geom_sf(colour=NA) +
  scale_fill_viridis_c(limits=c(0,1)) +
  facet_grid(meshRes~liceSpeed)
ggsave("figs/prSwim_mesh_by_speed_site.png", p, width=8, height=8, dpi=300)

p <- ggplot(elemAct.sf, aes(fill=float/total)) + geom_sf(colour=NA) +
  scale_fill_viridis_c(limits=c(0,1)) +
  facet_grid(meshRes~liceSpeed)
ggsave("figs/prFloat_mesh_by_speed_site.png", p, width=8, height=8, dpi=300)







# maps --------------------------------------------------------------------

bbox <- st_bbox(mesh.fp)
for(i in 1:nrow(sim_i)) {
  p <- loc.df %>%
    filter(age>=12 & status != 66) %>%
    filter(sim==i) %>%
    filter(between(x, bbox$xmin, bbox$xmax), between(y, bbox$ymin, bbox$ymax)) %>%
    ggplot(aes(x,y, alpha=density, colour=depth)) + 
    geom_point(size=0.1, shape=1) + 
    scale_colour_viridis_c(direction=-1, option="D", limits=c(0,20), na.value="#000033") +
    ggtitle(glue("{sim_i$mesh[i]}, {sim_i$timeRes[i]}, {sim_i$liceSpeedF[i]}"))
  ggsave(glue("figs/part_locs_sim{str_pad(i,2,'left','0')}_site.png"), 
         p, width=6, height=6, dpi=900)
}






# element summaries -------------------------------------------------------

# part.samp <- sample(unique(loc.df$ID), 1000)
# loc.df2 <- loc.df %>% filter(ID %in% part.samp) %>%
#   filter(status != 66)
# 
# loc.df2 %>% 
#   filter(liceSpeedF=="Medium") %>%
#   filter(date=="2021-11-02", hour >= 20) %>%
#   st_as_sf(coords=c("x", "y"), crs=27700) %>%
#   st_join(mesh.sf %>% filter(mesh=="WeStCOMS2") %>% select(i)) %>%
#   st_drop_geometry() %>%
#   group_by(sim, meshRes, liceSpeedF, timeCalculated, i) %>%
#   summarise(density=sum(density), N=n()) %>%
#   ungroup %>%
#   inner_join(mesh.sf %>% filter(mesh=="WeStCOMS2"), .) %>%
#   mutate(density=log(density/(area*depth))) %>%
#   ggplot() + 
#   geom_sf(data=mesh.fp, fill=NA, colour="grey30") +
#   geom_sf(aes(fill=density), colour=NA) + scale_fill_viridis_c() +
#   facet_grid(meshRes~timeCalculated)








# velocities --------------------------------------------------------------

velocity.df <- loc.df %>%
  select(ID, timeCalculated, age, status, x, y, density, depth, meshParticle,
         sim, mesh, meshRes, elem) %>%
  mutate(meshParticle=if_else(meshParticle==0, mesh, "WeStCOMS2")) %>%
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
  filter(age>=12 & status != 66) %>%
  left_join(., sim_i %>% select(sim, mesh, liceSpeedF)) %>%
  left_join(., mesh.sf %>% st_drop_geometry() %>% 
              select(i, mesh, bearing) %>% 
              rename(meshParticle=mesh, elem=i, lochDir=bearing),
            by=c("meshParticle", "elem")) %>%
  mutate(downloch=cos(lochDir-bearing)) %>% # -1: uploch, 1: downloch
  st_as_sf(coords=c("x", "y"), crs=27700, remove=F) %>%
  st_join(mesh.sf %>% filter(mesh=="WeStCOMS2") %>% 
            select(i, depth) %>% rename(elem_WC=i, meshDepth=depth)) %>%
  st_drop_geometry() %>% 
  filter(complete.cases(.))


# velocity.df %>%
#   filter(ID %% 500 == 0) %>%
#   ggplot(aes(x, y, colour=bearing)) +
#   geom_point(size=0.5) +
#   col_NSEW +
#   facet_grid(meshRes~liceSpeedF)
# 
# velocity.df %>%
#   filter(ID %% 500 == 0) %>%
#   ggplot(aes(x, y, colour=downloch)) +
#   geom_point() +
#   col_downloch +
#   facet_wrap(~hour(timeCalculated), ncol=6)



# compare across depths: absolute, or relative to meshDepth
# velocity.df %>%
#   mutate(depthCat=if_else(depth < 20, round(depth), 20)) %>%
#   ggplot(aes(downloch, colour=depthCat, group=depthCat)) +
#   geom_density() + 
#   scale_colour_viridis_c(direction=-1) +
#   facet_grid(~sim)
# velocity.df %>%
#   mutate(depthCat=if_else(depth<=5, "<= 5m", "> 5m")) %>%
#   ggplot(aes(downloch, colour=depthCat, group=depthCat)) +
#   geom_density() + 
#   facet_grid(~sim)
# velocity.df %>%
#   mutate(depthCat=if_else(depth<=5, "<= 5m", "> 5m")) %>%
#   ggplot(aes(downloch, colour=hour(timeCalculated), 
#              group=hour(timeCalculated))) +
#   geom_density() + 
#   col_hour +
#   facet_grid(depthCat~sim)


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
# ggsave("figs/direction_map_sim4_site.png", p, width=8, height=5, dpi=300)
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
# ggsave("figs/direction_map_sim3_site.png", p, width=8, height=5, dpi=300)

p <- velocity.df %>%
  filter(ID %% 2000 == 0) %>%
  ggplot() +
  geom_sf(data=mesh.fp, fill=NA, colour="grey30") +
  geom_path(aes(x, y, group=ID, colour=downloch), alpha=0.25, size=0.1) + 
  col_downloch +
  facet_grid(meshRes~liceSpeedF)
ggsave("figs/direction_map_site.png", p, width=9, height=9, dpi=300)



meshDir.sum <- velocity.df %>%
  group_by(mesh, meshRes, liceSpeedF, elem_WC) %>%
  summarise(N=n(), 
            totDens=sum(density), 
            prUploch=mean(downloch<0)) %>%
  ungroup %>%
  inner_join(mesh.sf %>% 
               filter(mesh=="WeStCOMS2") %>% 
               select(i, depth, area) %>% 
               rename(elem_WC=i), .) %>%
  mutate(N=replace_na(N, 0),
         lnDens_m3=log(totDens/(area*depth)),
         lnN_m3=log(N/(area*depth)))

p <- ggplot(meshDir.sum %>% filter(liceSpeedF=="Medium")) +
  geom_sf(data=mesh.fp, colour="grey30", fill=NA) +
  geom_sf(colour=NA, aes(fill=lnDens_m3)) + 
  scale_fill_viridis_c() +
  facet_grid(.~meshRes)
ggsave("figs/density_map_site.png", p, width=13, height=5, dpi=300)
meshDiff.df <- meshDir.sum %>%
  mutate(meshRes=droplevels(meshRes)) %>%
  full_join(meshDir.sum %>%
              st_drop_geometry() %>%
              filter(meshRes=="WeStCOMS2, 1h" & liceSpeedF=="Medium") %>%
              rename(N_ref=N, totDens_ref=totDens, prUp_ref=prUploch, 
                     lnDens_ref=lnDens_m3, lnN_ref=lnN_m3) %>%
              select(elem_WC, ends_with("ref"))) %>%
  filter(meshRes != "WeStCOMS2, 1h") %>%
  filter(liceSpeedF=="Medium")%>%
  mutate(N_ref=if_else(is.na(N_ref), 0L, N_ref),
         lnN_ref=if_else(is.na(lnN_ref), 0, lnN_ref),
         totDens_ref=if_else(is.na(totDens_ref), 0, totDens_ref)) %>%
  mutate(N_diff=N-N_ref, 
         lnN_diff=lnN_m3-lnN_ref,
         totDens_diff=totDens-totDens_ref,
         prUp_diff=prUploch-prUp_ref,
         lnDens_diff=lnDens_m3-lnDens_ref)
p <- meshDiff.df %>%
  ggplot() +
  geom_sf(data=mesh.fp, colour="grey30", fill=NA) +
  geom_sf(colour=NA, aes(fill=lnDens_diff)) + 
  scale_fill_gradient2() +
  facet_grid(.~meshRes)
ggsave("figs/density_diff_map_site.png", p, width=11, height=5, dpi=300)
p <- meshDiff.df %>%
  ggplot() +
  geom_sf(data=mesh.fp, colour="grey30", fill=NA) +
  geom_sf(colour=NA, aes(fill=lnN_diff/lnN_ref)) + 
  scale_fill_gradient2() +
  facet_grid(.~meshRes)
ggsave("figs/densityN_diff_map_site.png", p, width=11, height=5, dpi=300)
p <- meshDiff.df %>%
  ggplot() +
  geom_sf(data=mesh.fp, colour="grey30", fill=NA) +
  geom_sf(colour=NA, aes(fill=prUp_diff)) + 
  scale_fill_gradient2() +
  facet_grid(.~meshRes)
ggsave("figs/prUploch_diff_map_site.png", p, width=11, height=5, dpi=300)
p <- ggplot(meshDir.sum %>% filter(liceSpeedF=="Medium") %>% filter(N>10)) +
  geom_sf(data=mesh.fp, colour="grey30", fill=NA) +
  geom_sf(colour=NA, aes(fill=prUploch)) + 
  fill_downloch +
  facet_grid(.~meshRes)
ggsave("figs/prUploch_map_site.png", p, width=13, height=5, dpi=300)
p <- ggplot(meshDir.sum %>% filter(liceSpeedF=="Medium") %>% filter(N>10)) +
  geom_sf(data=mesh.fp, colour="grey30", fill=NA) +
  geom_sf(colour=NA, aes(fill=prUploch, alpha=lnDens_m3)) + 
  fill_downloch +
  facet_grid(.~meshRes)
ggsave("figs/prUploch_map_ALT_site.png", p, width=13, height=5, dpi=300)


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
ggsave("figs/density_map_depthStrat_site.png", p, width=11, height=14, dpi=300)
p <- ggplot(meshDir.deep %>% filter(!is.na(deep), N>3)) +
  geom_sf(data=mesh.fp, colour="grey30", fill=NA) +
  geom_sf(colour=NA, aes(fill=prUploch)) + 
  fill_downloch +
  facet_grid(deep~meshRes)
ggsave("figs/prUploch_map_depthStrat_site.png", p, width=11, height=14, dpi=300)
p <- ggplot(meshDir.deep %>% filter(!is.na(deep), N>3)) +
  geom_sf(data=mesh.fp, colour="grey30", fill=NA) +
  geom_sf(colour=NA, aes(fill=prUploch, alpha=lnDens_m3)) + 
  fill_downloch +
  facet_grid(deep~meshRes)
ggsave("figs/prUploch_map_depthStrat_ALT_site.png", p, width=11, height=14, dpi=300)
