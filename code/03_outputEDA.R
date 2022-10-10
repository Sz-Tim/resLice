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
  select(i, area, depth, mesh, geom)
mesh.fp <- st_read("data/linnhe_mesh_footprint.gpkg")
sim_i <- read_csv(glue("{dirs$out}/sim_i.csv")) %>%
  mutate(liceSpeedF=factor(liceSpeed, labels=c("Passive", "Slow", "Medium", "Fast")))
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
  facet_grid(liceSpeed~.)

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
  ggplot(aes(xy_cumul/1e3, linetype=mesh, colour=liceSpeedF, group=sim)) + 
  geom_density() + 
  scale_colour_viridis_d() +
  scale_x_log10("Total xy hourly displacement (km)")
ggsave("figs/xy_total_log.png", width=6, height=4, dpi=300)
loc.df %>% 
  group_by(ID, sim, mesh, liceSpeedF) %>%
  arrange(age) %>%
  mutate(xy=sqrt((x-lag(x))^2 + (y-lag(y))^2)) %>%
  filter(!is.na(xy)) %>%
  mutate(xy_cumul=cumsum(xy)) %>%
  filter(status==2) %>%
  slice_tail(n=1) %>%
  ggplot(aes(xy_cumul/1e3, linetype=mesh, colour=liceSpeedF, group=sim)) + 
  geom_density() + 
  scale_colour_viridis_d() +
  scale_x_continuous("Total xy hourly displacement (km)")
ggsave("figs/xy_total.png", width=6, height=4, dpi=300)

loc.df %>%
  group_by(ID, sim, mesh, liceSpeedF) %>%
  arrange(age) %>%
  filter(status==2) %>%
  slice_tail(n=1) %>%
  ggplot(aes(zTot, linetype=mesh, colour=liceSpeedF, group=sim)) + 
  geom_density() + 
  scale_colour_viridis_d() +
  scale_x_log10("Total z hourly displacement (m)")
ggsave("figs/z_total_log.png", width=6, height=4, dpi=300)
loc.df %>%
  group_by(ID, sim, mesh, liceSpeedF) %>%
  arrange(age) %>%
  filter(status==2) %>%
  slice_tail(n=1) %>%
  ggplot(aes(zTot, linetype=mesh, colour=liceSpeedF, group=sim)) + 
  geom_density() + 
  scale_colour_viridis_d() +
  scale_x_continuous("Total z hourly displacement (m)")
ggsave("figs/z_total.png", width=6, height=4, dpi=300)



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

