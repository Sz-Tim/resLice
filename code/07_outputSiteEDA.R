# Output EDA
# Off Aqua
# Tim Szewczyk


# This script does some exploratory analysis on the output

rm(list=ls()); gc()
# setup -------------------------------------------------------------------

library(tidyverse); library(glue); library(lubridate); library(sf)
source("code/00_fn.R")

initDensity <- c("Scaled", "Uniform")[2]

dirs <- switch(get_os(),
               windows=list(proj=getwd(),
                            mesh="D:/hydroOut/",
                            out=glue("{getwd()}/out/siteRelease_init{initDensity}/")),
               linux=list(proj=getwd(),
                          mesh="/home/sa04ts/FVCOM_meshes",
                          out=glue("{getwd()}/out/siteRelease_init{initDensity}/")))

loch_order <- c("Loch Linnhe (main)", "Loch Linnhe (head)", "Loch Eil", "Loch Leven", 
                "Loch Creran", "Loch Etive (mouth)", "Loch Etive (head)",
                "Loch a' Choire", "Loch Feochan", "Beyond Linnhe")
corran <- st_buffer(st_as_sfc("POINT(201850 763450)", crs=27700), 17.5e3)
etive.bbox <- list(xmin=185276, xmax=212845, ymin=728357, ymax=746971)



# load files --------------------------------------------------------------

mesh.sf <- list(linnhe7=st_read(glue("{dirs$mesh}/linnhe_mesh.gpkg")) %>%
                  mutate(mesh="linnhe7", meshParticle="linnhe7"),
                westcoms2=st_read(glue("{dirs$mesh}/WeStCOMS2_linnhe_mesh.gpkg")) %>%
                  mutate(mesh="WeStCOMS2", meshParticle="WeStCOMS2"),
                westcoms2_full=st_read(glue("{dirs$mesh}/WeStCOMS2_mesh.gpkg")) %>%
                  mutate(mesh="WeStCOMS2", meshParticle="Beyond Linnhe",
                         lochRegion=10, mainLoch=0)) %>%
  do.call('rbind', .) %>%
  select(i, area, depth, mesh, meshParticle, lochRegion, mainLoch, geom) %>%
  full_join(., read_csv("data/loch_regions.csv"), by="lochRegion") %>%
  left_join(., st_read("data/lochRegions_linnhe7.gpkg") %>% st_drop_geometry()) %>%
  rename(area_elem=area, depth_elem=depth) %>%
  mutate(vol_elem=area_elem*depth_elem,
         vol30_elem=area_elem*pmin(depth_elem, 30))
linnhe_elems <- mesh.sf %>% st_drop_geometry() %>%
  filter(meshParticle!="Beyond Linnhe") %>% select(mesh, i) %>% mutate(inLinnhe=T)
mesh.lochRegions <- mesh.sf %>% 
  group_by(mesh, mainLoch, lochRegion, lochRegionF) %>%
  summarise(bearing=mean(bearing),
            vol_region=sum(vol_elem),
            vol30_region=sum(vol30_elem)) 
mesh.fp <- st_read("data/linnhe_mesh_footprint.gpkg")
sim_i <- read_csv(glue("{dirs$out}/sim_i.csv")) %>%
  mutate(liceSpeedF=factor(liceSpeed, levels=c(0.0001, 0.0005, 0.001), 
                           labels=c("Slow", "Medium", "Fast")),
         sim=as.numeric(i))
elemAct.sf <- st_read(glue("out/00_processed/elementActivity_site_init{initDensity}.gpkg")) %>%
  mutate(liceSpeedF=factor(liceSpeedF, levels=levels(sim_i$liceSpeedF)),
         meshRes=factor(paste0(mesh, ", ", timeRes),
                        levels=c("WeStCOMS2, 1h", "WeStCOMS2, 5min", "linnhe7, 1h", "linnhe7, 5min")))
loc.df <- readRDS(glue("out/00_processed/locations_site_init{initDensity}.rds")) %>%
  mutate(liceSpeedF=factor(liceSpeedF, levels=levels(sim_i$liceSpeedF)),
         meshRes=factor(paste0(mesh, ", ", timeRes),
                        levels=c("WeStCOMS2, 1h", "WeStCOMS2, 5min", "linnhe7, 1h", "linnhe7, 5min")),
         hour=hour(timeCalculated),
         date=date(timeCalculated)) %>%
  left_join(linnhe_elems %>% rename(elem=i)) %>%
  mutate(meshParticle=if_else(is.na(inLinnhe), "Beyond Linnhe", mesh)) %>%
  left_join(mesh.sf %>% 
              st_drop_geometry() %>% 
              rename(elem=i) %>%
              select(meshParticle, elem, lochRegionF, mainLoch, 
                     area_elem, depth_elem, vol_elem, vol30_elem)) %>%
  left_join(mesh.lochRegions %>% st_drop_geometry()) %>%
  mutate(lochRegionF=if_else(is.na(lochRegionF), "Beyond Linnhe", lochRegionF),
         lochRegionF=factor(lochRegionF, levels=loch_order),
         mainLoch=if_else(lochRegionF=="Beyond Linnhe", 2, mainLoch),
         mainLoch=factor(mainLoch, levels=c(2,0,1), 
                         labels=c("Beyond Linnhe", "Side lochs", "Main loch"))) 
beyondLinnhe_vol <- loc.df %>% 
  filter(lochRegionF=="Beyond Linnhe") %>%
  group_by(elem) %>%
  slice_head(n=1) %>%
  ungroup %>%
  summarise(vol_region=sum(vol_elem, na.rm=T),
            vol30_region=sum(vol30_elem, na.rm=T))
mainSide_vol <- mesh.lochRegions %>% st_drop_geometry() %>%
  filter(lochRegion != 10) %>%
  group_by(mesh, mainLoch) %>%
  summarise(vol_mainSide=sum(vol_region),
            vol30_mainSide=sum(vol30_region)) %>%
  ungroup %>%
  rename(meshParticle=mesh) %>%
  mutate(mainLoch=if_else(mainLoch==0, "Side lochs", "Main loch")) %>%
  bind_rows(beyondLinnhe_vol %>% 
              mutate(meshParticle="Beyond Linnhe", 
                     mainLoch="Beyond Linnhe") %>%
              rename(vol_mainSide=vol_region, vol30_mainSide=vol30_region))
loc.df <- loc.df %>%
  mutate(vol_region=if_else(lochRegionF=="Beyond Linnhe", beyondLinnhe_vol$vol_region, vol_region),
         vol30_region=if_else(lochRegionF=="Beyond Linnhe", beyondLinnhe_vol$vol30_region, vol30_region)) %>%
  left_join(mainSide_vol)
# loc.sf <- st_as_sf(loc.df, coords=c("x", "y"), crs=27700)
gc()


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
col_meshRes <- c("linnhe7, 1h"="#a6611a", "linnhe7, 5min"="#dfc27d",
                 "WeStCOMS2, 1h"="#018571", "WeStCOMS2, 5min"="#80cdc1")








# particle attributes -----------------------------------------------------
# part.samp <- sample(n_distinct(loc.df$ID), 8*14)
# 
# loc.df %>% 
#   filter(liceSpeedF=="Medium") %>% 
#   filter(ID %in% part.samp) %>% 
#   filter(density > 0 & density < 1) %>% 
#   ggplot(aes(age, brms::logit_scaled(density), colour=meshRes)) + 
#   geom_line() + 
#   scale_colour_manual(values=col_meshRes) + 
#   facet_wrap(~ID, ncol=14)
# loc.df %>% 
#   group_by(ID, meshRes, liceSpeedF) %>%
#   filter(age==max(age)) %>%
#   ggplot(aes(degreeDays, colour=meshRes)) + 
#   geom_density() + 
#   scale_colour_manual(values=col_meshRes) + 
#   facet_wrap(~liceSpeedF)
# loc.df %>% 
#   filter(startDate == 20211101) %>%
#   group_by(ID, meshRes, liceSpeedF) %>%
#   filter(age==max(age)) %>%
#   ggplot(aes(brms::logit_scaled(density), colour=meshRes)) + 
#   geom_density() + 
#   scale_colour_manual(values=col_meshRes) + 
#   facet_wrap(~liceSpeedF)
# library(brms)
# part.samp <- sample(unique(filter(loc.df, timeRes=="1h", startDate==20211101)$ID), 20)
# aging.df <- loc.df %>% 
#   # filter(startDate == 20211101) %>%
#   filter(ID %in% part.samp) %>%
#   filter(status != 666) %>%
#   filter(density > 0 & density < 1) %>%
#   mutate(dens_logit=logit_scaled(density),
#          age_sqrt=sqrt(age))
# ageDensity.brm <- brm(dens_logit ~ age_sqrt*meshRes*liceSpeedF + (1+age_sqrt|startLocation/ID),
#                       data=aging.df, cores=4)

library(ggridges)
loc.df %>% 
  filter(status==2) %>%
  filter(density > 0) %>%
  mutate(density=brms::logit_scaled(density),
         mainLoch=if_else(is.na(mainLoch), 2L, mainLoch),
         mainLoch=factor(mainLoch, levels=c(2,0,1), 
                         labels=c("Beyond Linnhe", "Side lochs", "Main loch"))) %>%
  ggplot(aes(density, y=mainLoch, colour=meshRes, fill=meshRes)) + 
  ggridges::geom_density_ridges(alpha=0.25, scale=0.9) +
  scale_fill_manual("Mesh", values=col_meshRes) +
  scale_colour_manual("Mesh", values=col_meshRes) +
  facet_grid(.~liceSpeedF) +
  labs(x="logit(Per-particle copepodid density)", y="")
ggsave(glue("figs/init{initDensity}/copepodDensity_mainLoch.png"), height=5, width=7)
loc.df %>% 
  filter(status==2) %>%
  filter(density > 0) %>%
  mutate(density=brms::logit_scaled(density),
         lochRegion=if_else(is.na(lochRegion), 8, lochRegion),
         lochRegion=factor(lochRegion, levels=1:8, 
                           labels=c(unique(mesh.lochRegions$lochRegionF), "Beyond Linnhe"))) %>%
  ggplot(aes(density, y=lochRegion, colour=meshRes, fill=meshRes)) + 
  ggridges::geom_density_ridges(alpha=0.25, scale=0.9) +
  scale_fill_manual("Mesh", values=col_meshRes) +
  scale_colour_manual("Mesh", values=col_meshRes) +
  facet_grid(.~liceSpeedF) +
  labs(x="logit(Per-particle copepodid density)", y="")
ggsave(glue("figs/init{initDensity}/copepodDensity_lochRegion.png"), height=7, width=7)
loc.df %>% 
  filter(status==2) %>%
  group_by(mesh, timeRes, meshRes, liceSpeedF, lochRegion) %>% 
  summarise(nPart=n(), 
            dens_tot=sum(density)) %>%
  mutate(lochRegion=if_else(is.na(lochRegion), 8, lochRegion),
         lochRegion=factor(lochRegion, levels=1:8, 
                           labels=c(unique(mesh.lochRegions$lochRegionF), "Beyond Linnhe")),
         lochRegion=factor(lochRegion, levels=levels(lochRegion)[c(6,1:3,4,5,7,8)])) %>%
  ggplot(aes(liceSpeedF, dens_tot, fill=meshRes)) +
  geom_bar(stat="identity", position="dodge", colour="grey30") + 
  scale_fill_manual("Mesh", values=col_meshRes) +
  facet_wrap(~lochRegion, nrow=2) +
  labs(x="Lice vertical speed", y="Integrated copepodid density") +
  theme(axis.text.x=element_text(angle=310, hjust=0, vjust=0.5),
        legend.position="bottom")
ggsave(glue("figs/init{initDensity}/copepodDensityTotal_lochRegion.png"), height=6, width=8)
loc.df %>% 
  filter(status==2) %>%
  filter(density > 0) %>%
  group_by(mesh, timeRes, meshRes, liceSpeedF, lochRegion) %>% 
  summarise(nPart=n()) %>%
  group_by(mesh, timeRes, meshRes, liceSpeedF) %>%
  mutate(propPart=nPart/sum(nPart),
         lochRegion=if_else(is.na(lochRegion), 8, lochRegion),
         lochRegion=factor(lochRegion, levels=1:8, 
                           labels=c(unique(mesh.lochRegions$lochRegionF), "Beyond Linnhe")),
         lochRegion=factor(lochRegion, levels=levels(lochRegion)[c(6,1:3,4,5,7,8)])) %>%
  ggplot(aes(liceSpeedF, propPart, fill=meshRes)) +
  geom_bar(stat="identity", position="dodge", colour="grey30") + 
  scale_fill_manual("Mesh", values=col_meshRes) +
  facet_wrap(~lochRegion, nrow=2) +
  labs(x="Lice vertical speed", y="Proportion of particle-hours") +
  theme(axis.text.x=element_text(angle=310, hjust=0, vjust=0.5),
        legend.position="bottom")
ggsave(glue("figs/init{initDensity}/copepodPropPartHr_lochRegion.png"), height=6, width=8)

loc.df %>% 
  filter(status==2) %>%
  group_by(mesh, timeRes, meshRes, liceSpeedF, mainLoch) %>% 
  summarise(nPart=n(), 
            dens_tot=sum(density)) %>%
  mutate(mainLoch=if_else(is.na(mainLoch), 2L, mainLoch),
         mainLoch=factor(mainLoch, levels=0:2, labels=c("Side lochs", "Main loch", "Beyond Linnhe"))) %>%
  ggplot(aes(liceSpeedF, dens_tot, fill=meshRes)) +
  geom_bar(stat="identity", position="dodge", colour="grey30") + 
  scale_fill_manual("Mesh", values=col_meshRes) +
  facet_grid(.~mainLoch) +
  labs(x="Lice vertical speed", y="Integrated copepodid density") +
  theme(axis.text.x=element_text(angle=310, hjust=0, vjust=0.5),
        legend.position="bottom")
ggsave(glue("figs/init{initDensity}/copepodDensityTotal_mainLoch.png"), height=4, width=6)
loc.df %>% 
  filter(status==2) %>%
  # filter(density > 0) %>%
  group_by(mesh, timeRes, meshRes, liceSpeedF, mainLoch) %>% 
  summarise(nPart=n()) %>%
  group_by(mesh, timeRes, meshRes, liceSpeedF) %>%
  mutate(propPart=nPart/sum(nPart),
         mainLoch=if_else(is.na(mainLoch), 2L, mainLoch),
         mainLoch=factor(mainLoch, levels=0:2, labels=c("Side lochs", "Main loch", "Beyond Linnhe"))) %>%
  ggplot(aes(liceSpeedF, propPart, fill=meshRes)) +
  geom_bar(stat="identity", position="dodge", colour="grey30") + 
  scale_fill_manual("Mesh", values=col_meshRes) +
  facet_grid(.~mainLoch) +
  labs(x="Lice vertical speed", y="Proportion of particle-hours") +
  theme(axis.text.x=element_text(angle=310, hjust=0, vjust=0.5),
        legend.position="bottom")
ggsave(glue("figs/init{initDensity}/copepodPropPartHr_mainLoch.png"), height=4, width=6)


loc.df %>% 
  filter(status==2) %>%
  group_by(mesh, timeRes, meshRes, liceSpeedF, mainLoch, timeCalculated) %>% 
  summarise(dens_tot=sum(density)) %>%
  mutate(mainLoch=if_else(is.na(mainLoch), 2L, mainLoch),
         mainLoch=factor(mainLoch, levels=0:2, labels=c("Side lochs", "Main loch", "Beyond Linnhe"))) %>%
  ggplot(aes(timeCalculated, dens_tot, colour=meshRes, linetype=liceSpeedF)) + 
  geom_line() + 
  scale_linetype_manual("Lice\nspeed", values=1:3) +
  scale_colour_manual("Mesh", values=col_meshRes) +
  guides(colour=guide_legend(ncol=2)) +
  facet_grid(.~mainLoch) +
  labs(x="Date", y="Copepodid density") +
  theme(legend.position="bottom",
        axis.text.x=element_text(angle=310, hjust=0, vjust=0.5))
ggsave(glue("figs/init{initDensity}/copepodTimeDensity_mainLoch.png"), height=4, width=8)

loc.df %>% 
  filter(status==2) %>%
  group_by(mesh, timeRes, meshRes, liceSpeedF, lochRegion, timeCalculated) %>% 
  summarise(dens_tot=sum(density)) %>%
  mutate(lochRegion=if_else(is.na(lochRegion), 8, lochRegion),
         lochRegion=factor(lochRegion, levels=1:8, 
                           labels=c(unique(mesh.lochRegions$lochRegionF), "Beyond Linnhe")),
         lochRegion=factor(lochRegion, levels=levels(lochRegion)[c(6,1:3,4,5,7,8)])) %>%
  ggplot(aes(timeCalculated, dens_tot, colour=meshRes, linetype=liceSpeedF)) + 
  geom_line() + 
  scale_linetype_manual("Lice\nspeed", values=1:3) +
  scale_colour_manual("Mesh", values=col_meshRes) +
  guides(colour=guide_legend(ncol=2)) +
  facet_wrap(~lochRegion, nrow=2) +
  labs(x="Date", y="Copepodid density") +
  theme(legend.position="bottom",
        axis.text.x=element_text(angle=310, hjust=0, vjust=0.5))
ggsave(glue("figs/init{initDensity}/copepodTimeDensity_lochRegion.png"), height=6, width=8)


loc.df %>% 
  filter(status==2) %>%
  group_by(mesh, timeRes, meshRes, liceSpeedF, mainLoch, timeCalculated) %>% 
  summarise(nPart=n()) %>%
  group_by(mesh, timeRes, meshRes, liceSpeedF) %>%
  mutate(propPart=nPart/sum(nPart),
         mainLoch=if_else(is.na(mainLoch), 2L, mainLoch),
         mainLoch=factor(mainLoch, levels=0:2, labels=c("Side lochs", "Main loch", "Beyond Linnhe"))) %>%
  ggplot(aes(timeCalculated, nPart, colour=meshRes, linetype=liceSpeedF)) + 
  geom_line() + 
  scale_linetype_manual("Lice\nspeed", values=1:3) +
  scale_colour_manual("Mesh", values=col_meshRes) +
  guides(colour=guide_legend(ncol=2)) +
  facet_grid(.~mainLoch) +
  labs(x="Date", y="Proportion of copepodid particles") +
  theme(legend.position="bottom",
        axis.text.x=element_text(angle=310, hjust=0, vjust=0.5))
ggsave(glue("figs/init{initDensity}/copepodTimePropPartHr_mainLoch.png"), height=4, width=8)

loc.df %>% 
  filter(status==2) %>%
  group_by(mesh, timeRes, meshRes, liceSpeedF, lochRegion, timeCalculated) %>% 
  summarise(nPart=n()) %>%
  group_by(mesh, timeRes, meshRes, liceSpeedF) %>%
  mutate(propPart=nPart/sum(nPart),lochRegion=if_else(is.na(lochRegion), 8, lochRegion),
         lochRegion=factor(lochRegion, levels=1:8, 
                           labels=c(unique(mesh.lochRegions$lochRegionF), "Beyond Linnhe")),
         lochRegion=factor(lochRegion, levels=levels(lochRegion)[c(6,1:3,4,5,7,8)])) %>%
  ggplot(aes(timeCalculated, nPart, colour=meshRes, linetype=liceSpeedF)) + 
  geom_line() + 
  scale_linetype_manual("Lice\nspeed", values=1:3) +
  scale_colour_manual("Mesh", values=col_meshRes) +
  guides(colour=guide_legend(ncol=2)) +
  facet_wrap(~lochRegion, nrow=2) +
  labs(x="Date", y="Copepodid density") +
  theme(legend.position="bottom",
        axis.text.x=element_text(angle=310, hjust=0, vjust=0.5))
ggsave(glue("figs/init{initDensity}/copepodTimePropPartHr_lochRegion.png"), height=6, width=8)


loc.df %>% 
  filter(status==2) %>%
  filter(density > 0) %>%
  filter(timeCalculated >= "2021-11-04 00:00:00") %>%
  mutate(density=brms::logit_scaled(density),
         mainLoch=if_else(is.na(mainLoch), 2L, mainLoch),
         mainLoch=factor(mainLoch, levels=c(2,0,1), 
                         labels=c("Beyond Linnhe", "Side lochs", "Main loch"))) %>%
  group_by(meshRes, liceSpeedF, mainLoch, timeCalculated) %>%
  summarise(dens_mn=mean(density),
            dens_lo=quantile(density, 0.25),
            dens_hi=quantile(density, 0.75)) %>%
  ggplot(aes(timeCalculated, dens_mn, ymin=dens_lo, ymax=dens_hi, colour=meshRes, fill=meshRes)) + 
  geom_ribbon(alpha=0.5, colour=NA) +
  geom_line() + 
  # geom_errorbar(position=position_dodge(width=0.4), width=0.2) +
  scale_fill_manual("Mesh", values=col_meshRes) +
  scale_colour_manual("Mesh", values=col_meshRes) +
  facet_grid(liceSpeedF~mainLoch) +
  theme(legend.position="bottom",
        axis.text.x=element_text(angle=310, hjust=0, vjust=0.5)) +
  labs(x="", y="logit(Per-particle copepodid density)")
ggsave(glue("figs/init{initDensity}/copepodTimePartDensity_mainLoch.png"), height=5, width=8)

loc.df %>% 
  filter(status==2) %>%
  filter(density > 0) %>%
  filter(timeCalculated >= "2021-11-04 00:00:00") %>%
  mutate(density=brms::logit_scaled(density),
         lochRegion=if_else(is.na(lochRegion), 8, lochRegion),
         lochRegion=factor(lochRegion, levels=1:8, 
                           labels=c(unique(mesh.lochRegions$lochRegionF), "Beyond Linnhe")),
         lochRegion=factor(lochRegion, levels=levels(lochRegion)[c(6,1:3,4,5,7,8)])) %>%
  group_by(meshRes, liceSpeedF, lochRegion, timeCalculated) %>%
  summarise(dens_mn=mean(density),
            dens_lo=quantile(density, 0.25),
            dens_hi=quantile(density, 0.75)) %>%
  ggplot(aes(timeCalculated, dens_mn, ymin=dens_lo, ymax=dens_hi, colour=meshRes, fill=meshRes)) + 
  geom_ribbon(alpha=0.5, colour=NA) +
  geom_line() + 
  # geom_errorbar(position=position_dodge(width=0.4), width=0.2) +
  scale_fill_manual("Mesh", values=col_meshRes) +
  scale_colour_manual("Mesh", values=col_meshRes) +
  facet_grid(liceSpeedF~lochRegion) +
  theme(legend.position="bottom",
        axis.text.x=element_text(angle=310, hjust=0, vjust=0.5)) +
  labs(x="", y="logit(Per-particle copepodid density)")
ggsave(glue("figs/init{initDensity}/copepodTimePart_lochRegion.png"), height=5, width=11)





# tidal pulses ------------------------------------------------------------

pulse.sf <- st_read(glue("out/00_processed/locations_sitePulse_init{initDensity}.gpkg")) %>%
  mutate(downloch_m=cos(-2.225858-bearing_m), # bearing of main body (loch region 7)
         downloch_p=cos(-2.225858-bearing_p),
         meshRes=paste(mesh, timeRes, sep=", "),
         meshRes=factor(meshRes, levels=names(col_meshRes)),
         liceSpeedF=factor(liceSpeedF, levels=levels(sim_i$liceSpeedF)))

pulse.dir <- pulse.sf %>%
  st_drop_geometry() %>%
  select(ID, sim, mesh, timeRes, meshRes, liceSpeedF, timeCalculated, downloch_m, downloch_p) %>%
  pivot_longer(starts_with("downloch"), names_to="Time", values_to="downloch") %>%
  mutate(Time=if_else(Time=="downloch_m", "before", "after"),
         Time=factor(Time, levels=c("before", "after")),
         particleDirection=if_else(downloch>0, "downloch", "uploch"))
pulse.dZ <- pulse.sf %>%
  st_drop_geometry() %>%
  select(ID, sim, mesh, timeRes, meshRes, liceSpeedF, timeCalculated, z, dmz, dpz) %>%
  pivot_longer(c("dmz", "dpz"), names_to="Time", values_to="dZ") %>%
  mutate(Time=if_else(Time=="dmz", "before", "after"),
         Time=factor(Time, levels=c("before", "after")))
pulse.dXY <- pulse.sf %>%
  st_drop_geometry() %>%
  select(ID, sim, mesh, timeRes, meshRes, liceSpeedF, timeCalculated, dmXY, dpXY) %>%
  pivot_longer(c("dmXY", "dpXY"), names_to="Time", values_to="dXY") %>%
  mutate(Time=if_else(Time=="dmXY", "before", "after"),
         Time=factor(Time, levels=c("before", "after")))
pulse.dDens <- pulse.sf %>%
  st_drop_geometry() %>%
  mutate(across(starts_with("density"), ~brms::logit_scaled(.x, lb=-0.0001, ub=1.0001))) %>%
  mutate(dmDens=density-density_m1, dpDens=density_p1-density) %>% 
  select(ID, sim, mesh, timeRes, meshRes, liceSpeedF, timeCalculated, density, dmDens, dpDens) %>%
  pivot_longer(c("dmDens", "dpDens"), names_to="Time", values_to="dDens") %>%
  mutate(Time=if_else(Time=="dmDens", "before", "after"),
         Time=factor(Time, levels=c("before", "after")))
pulse.dDens$dDens[pulse.dDens$dDens > 0] <- NA
pulse.sum <- pulse.dir %>%
  mutate(z0=pulse.dZ$z,
         dZ=pulse.dZ$dZ,
         dXY=pulse.dXY$dXY,
         dDens=pulse.dDens$dDens) %>%
  filter(!is.na(particleDirection))

pulse.Z <- pulse.sf %>%
  st_drop_geometry() %>%
  select(ID, sim, mesh, timeRes, meshRes, liceSpeedF, timeCalculated, z, z_m1, z_p1) %>%
  pivot_longer(starts_with("z"), names_to="Time", values_to="z") %>%
  mutate(Time=case_when(Time=="z_m1" ~ -1,
                        Time=="z" ~ 0,
                        Time=="z_p1" ~ 1),
         Time=factor(Time))
pulse.Dens <- pulse.sf %>%
  st_drop_geometry() %>%
  select(ID, sim, mesh, timeRes, meshRes, liceSpeedF, timeCalculated, density, density_m1, density_p1) %>%
  pivot_longer(starts_with("density"), names_to="Time", values_to="density") %>%
  mutate(Time=case_when(Time=="density_m1" ~ -1,
                        Time=="density" ~ 0,
                        Time=="density_p1" ~ 1),
         Time=factor(Time))


pulse.dir %>%
  filter(!is.na(downloch))%>%
  ggplot(aes(Time, fill=particleDirection)) +
  geom_bar(position="fill", colour="grey30") +
  facet_grid(liceSpeedF~meshRes) +
  labs(title="Direction of movement pre- and post-tidal pulse",
       x="Time re: tidal pulse", y="Proportion of particles")
ggsave(glue("figs/init{initDensity}/pulse_downlochBar.png"), width=8, height=6, dpi=300)

glm.dat <- pulse.dir %>%
  mutate(particleDirection=if_else(downloch>0, 1, 0),
         Time=factor(Time)) %>%
  filter(!is.na(particleDirection)) %>%
  group_by(mesh, timeRes, liceSpeedF, Time, timeCalculated) %>%
  summarise(down=sum(particleDirection),
            Npart=n()) %>%
  mutate(up=Npart-down)
library(brms)
out <- brm(up | trials(Npart) ~ mesh * timeRes * liceSpeedF * Time,
           data=glm.dat, family=binomial(), cores=4)
summary(out)
pred.df <- glm.dat %>%
  mutate(obs_pr=up/Npart) %>%
  group_by(mesh, timeRes, liceSpeedF, Time) %>%
  summarise(obs_mn=mean(obs_pr),
            obs_lo1=HDInterval::hdi(obs_pr)[1],
            obs_hi1=HDInterval::hdi(obs_pr)[2],
            obs_lo2=HDInterval::hdi(obs_pr, 0.8)[1],
            obs_hi2=HDInterval::hdi(obs_pr, 0.8)[2],
            obs_lo3=HDInterval::hdi(obs_pr, 0.5)[1],
            obs_hi3=HDInterval::hdi(obs_pr, 0.5)[2]) %>%
  ungroup %>%
  mutate(Npart=1e2)
preds <- posterior_predict(out, newdata=pred.df)
pred.df <- pred.df %>%
  mutate(pred_mn=colMeans(preds),
         pred_lo1=HDInterval::hdi(preds, 0.99)[1,],
         pred_hi1=HDInterval::hdi(preds, 0.99)[2,],
         pred_lo2=HDInterval::hdi(preds, 0.8)[1,],
         pred_hi2=HDInterval::hdi(preds, 0.8)[2,],
         pred_lo3=HDInterval::hdi(preds, 0.5)[1,],
         pred_hi3=HDInterval::hdi(preds, 0.5)[2,]) %>%
  mutate(across(starts_with("pred_"), ~.x/Npart),
         meshRes=paste(mesh, timeRes, sep=", "))
pred.df %>%
  mutate(Time=factor(Time, levels=rev(levels(Time)))) %>%
  ggplot(aes(pred_mn, liceSpeedF, colour=meshRes, shape=Time)) +
  geom_point(position=position_dodge(width=0.5), size=3) +
  geom_errorbarh(aes(xmin=pred_lo1, xmax=pred_hi1),
                 position=position_dodge(width=0.5), size=0.5, height=0.25) +
  scale_colour_manual("Mesh", values=col_meshRes) +
  scale_shape_manual("Time re:\ntidal pulse", values=c(19, 1)) +
  guides(shape=guide_legend(reverse=T),
         linetype=guide_legend(reverse=T)) +
  labs(#title="Direction of movement pre- and post-tidal pulse",
       x="Pr(uploch)", y="Lice swim/sink speed") +
  theme_bw() +
  theme(panel.grid.minor=element_blank())
ggsave(glue("figs/init{initDensity}/pulse_uplochPts.png"), width=6, height=4, dpi=300)
pred.df %>%
  ggplot(aes(Time, pred_mn, colour=meshRes, group=meshRes)) +
  geom_point(size=3) +
  geom_line() +
  geom_errorbar(aes(ymin=pred_lo3, ymax=pred_hi3), size=0.5, width=0.25) +
  scale_colour_manual("Mesh", values=col_meshRes) +
  labs(x="Time re: tidal pulse", y="Pr(uploch)") +
  facet_grid(.~liceSpeedF) +
  theme(panel.grid.minor=element_blank())
ggsave(glue("figs/init{initDensity}/pulse_uplochBump.png"), width=6, height=4, dpi=300)

pred.df %>%
  mutate(Time=factor(Time, levels=rev(levels(Time)))) %>%
  ggplot(aes(liceSpeedF, pred_mn, colour=meshRes, shape=Time)) +
  geom_point(size=2) +
  geom_line(aes(group=paste(Time, meshRes), linetype=Time)) +
  geom_errorbar(aes(ymin=pred_lo1, ymax=pred_hi1), size=0.25, width=0.1) +
  scale_colour_manual("Mesh", values=col_meshRes) +
  scale_shape_manual("Time re:\ntidal pulse", values=c(19, 1)) +
  scale_linetype_manual("Time re:\ntidal pulse", values=c(1, 2)) +
  guides(shape=guide_legend(reverse=T),
         linetype=guide_legend(reverse=T)) +
  labs(title="Mean direction of movement pre- and post-tidal pulse",
       x="Lice swim/sink speed", y="Pr(uploch)")
ggsave(glue("figs/init{initDensity}/pulse_uplochLines.png"), width=6, height=5, dpi=300)


# glm_zCat.dat <- pulse.sum %>%
#   mutate(particleDirection=if_else(downloch>0, 1, 0),
#          Time=factor(Time),
#          dZ=pmax(pmin(-dZ, 20), -20),
#          dZ_cat=factor(round(dZ, -1))) %>%
#   filter(!is.na(particleDirection)) %>%
#   group_by(mesh, timeRes, liceSpeedF, Time, dZ_cat, timeCalculated) %>%
#   summarise(down=sum(particleDirection),
#             Npart=n())
# library(brms)
# out.Z <- brm(down | trials(Npart) ~ mesh * timeRes * liceSpeedF * Time * dZ_cat,
#            data=glm_zCat.dat, family=binomial(), cores=4,
#            control=list(max_treedepth=20))
# summary(out.Z)
# preds <- posterior_predict(out)
# pred.df <- glm.dat %>%
#   ungroup %>%
#   mutate(pred_mn=colMeans(preds),
#          pred_lo1=HDInterval::hdi(preds, 0.999)[1,],
#          pred_hi1=HDInterval::hdi(preds, 0.999)[2,],
#          pred_lo2=HDInterval::hdi(preds, 0.8)[1,],
#          pred_hi2=HDInterval::hdi(preds, 0.8)[2,],
#          pred_lo3=HDInterval::hdi(preds, 0.5)[1,],
#          pred_hi3=HDInterval::hdi(preds, 0.5)[2,]) %>%
#   mutate(across(starts_with("pred_"), ~1-.x/Npart),
#          meshRes=paste(mesh, timeRes, sep=", "))


# more downward movement post-tidal pulse, esp. 5min
# out.dZ <- lm(dZ ~ meshRes * liceSpeedF * Time * particleDirection,
              # data=pulse.sum %>% mutate(dZ=-dZ))
pulse.sum %>%
  mutate(dZ=-dZ) %>%
  group_by(meshRes, liceSpeedF, Time, particleDirection) %>%
  summarise(dZ_mn=median(dZ, na.rm=T),
            dZ_lo1=HDInterval::hdi(dZ)[1],
            dZ_hi1=HDInterval::hdi(dZ)[2],
            dZ_lo2=HDInterval::hdi(dZ, 0.8)[1],
            dZ_hi2=HDInterval::hdi(dZ, 0.8)[2],
            dZ_lo3=HDInterval::hdi(dZ, 0.5)[1],
            dZ_hi3=HDInterval::hdi(dZ, 0.5)[2]) %>%
  ggplot(aes(particleDirection, dZ_mn, colour=meshRes, shape=Time)) +
  geom_point(position=position_dodge(width=0.5), size=2) +
  # geom_linerange(aes(ymin=dZ_lo1, ymax=dZ_hi1), position=position_dodge(width=0.5), size=0.25) +
  geom_linerange(aes(ymin=dZ_lo2, ymax=dZ_hi2), position=position_dodge(width=0.5), size=0.5) +
  geom_linerange(aes(ymin=dZ_lo3, ymax=dZ_hi3), position=position_dodge(width=0.5), size=1) +
  scale_colour_manual(values=col_meshRes) +
  scale_shape_manual(values=c(1, 19)) +
  facet_grid(.~liceSpeedF) +
  labs(x="Lice sink/swim speed", y="\u0394z (+ up, - down)") +
  ggtitle("Vertical movement pre- and post-tidal pulse") +
  theme(legend.position="bottom")
ggsave(glue("figs/init{initDensity}/pulse_dZ_hdi.png"), width=9, height=4, dpi=300)

pulse.sum %>%
  group_by(meshRes, liceSpeedF, Time, particleDirection) %>%
  summarise(dDens_mn=median(dDens, na.rm=T),
            dDens_lo1=HDInterval::hdi(dDens)[1],
            dDens_hi1=HDInterval::hdi(dDens)[2],
            dDens_lo2=HDInterval::hdi(dDens, 0.8)[1],
            dDens_hi2=HDInterval::hdi(dDens, 0.8)[2],
            dDens_lo3=HDInterval::hdi(dDens, 0.5)[1],
            dDens_hi3=HDInterval::hdi(dDens, 0.5)[2]) %>%
  ggplot(aes(particleDirection, dDens_mn, colour=meshRes, shape=Time)) +
  geom_point(position=position_dodge(width=0.5), size=2) +
  # geom_linerange(aes(ymin=dDens_lo1, ymax=dDens_hi1), position=position_dodge(width=0.5), size=0.25) +
  geom_linerange(aes(ymin=dDens_lo2, ymax=dDens_hi2), position=position_dodge(width=0.5), size=0.5) +
  geom_linerange(aes(ymin=dDens_lo3, ymax=dDens_hi3), position=position_dodge(width=0.5), size=1) +
  scale_colour_manual(values=col_meshRes) +
  scale_shape_manual(values=c(1, 19)) +
  facet_grid(.~liceSpeedF) +
  labs(x="Lice sink/swim speed", y="\u0394 density (mortality)") +
  ggtitle("Mortality pre- and post-tidal pulse") +
  theme(legend.position="bottom")
ggsave(glue("figs/init{initDensity}/pulse_dDens_hdi.png"), width=9, height=4, dpi=300)

pulse.sum %>%
  mutate(dZ=-dZ) %>%
  group_by(meshRes, liceSpeedF, Time, particleDirection, timeCalculated) %>%
  summarise(dZ_mn=median(dZ, na.rm=T)) %>%
  ggplot(aes(timeCalculated, dZ_mn, colour=particleDirection, shape=Time)) +
  geom_point() +
  # scale_colour_manual(values=col_meshRes) +
  scale_shape_manual(values=c(1, 19)) +
  facet_grid(meshRes~liceSpeedF)


# more xy movement post-tidal pulse, esp. 5min
pulse.sum %>%
  mutate(dXY=log(dXY)) %>%
  group_by(meshRes, liceSpeedF, Time, particleDirection) %>%
  summarise(dXY_mn=median(dXY, na.rm=T),
            dXY_lo1=HDInterval::hdi(dXY)[1],
            dXY_hi1=HDInterval::hdi(dXY)[2],
            dXY_lo2=HDInterval::hdi(dXY, 0.8)[1],
            dXY_hi2=HDInterval::hdi(dXY, 0.8)[2],
            dXY_lo3=HDInterval::hdi(dXY, 0.5)[1],
            dXY_hi3=HDInterval::hdi(dXY, 0.5)[2]) %>%
  ggplot(aes(particleDirection, dXY_mn, colour=meshRes, shape=Time)) +
  geom_point(position=position_dodge(width=0.5), size=2) +
  # geom_linerange(aes(ymin=dXY_lo1, ymax=dXY_hi1), position=position_dodge(width=0.5), size=0.25) +
  geom_linerange(aes(ymin=dXY_lo2, ymax=dXY_hi2), position=position_dodge(width=0.5), size=0.5) +
  geom_linerange(aes(ymin=dXY_lo3, ymax=dXY_hi3), position=position_dodge(width=0.5), size=1) +
  scale_colour_manual(values=col_meshRes) +
  scale_shape_manual(values=c(1, 19)) +
  facet_grid(.~liceSpeedF) +
  labs(x="Lice sink/swim speed", y="\u0394ln(XY)") +
  theme(legend.position="bottom")
ggsave(glue("figs/init{initDensity}/pulse_lndXY_hdi.png"), width=9, height=4, dpi=300)


pulse.sum %>%
  group_by(meshRes, liceSpeedF, Time, particleDirection) %>%
  sample_frac(0.1) %>%
  ungroup %>%
  mutate(dZ=-dZ,
         dXY=log(dXY)) %>%
  ggplot(aes(dZ, dXY, colour=Time)) +
  geom_point(alpha=0.25, size=0.5) +
  # scale_shape_manual(values=c(1, 19)) +
  facet_grid(particleDirection~meshRes*liceSpeedF)
pulse.sum %>%
  ungroup %>%
  mutate(dZ=-dZ,
         dXY=log(dXY)) %>%
  ggplot(aes(dZ, dXY)) +
  geom_density_2d_filled() +
  facet_grid(Time*particleDirection~meshRes*liceSpeedF)
pulse.sum %>%
  mutate(dZ=pmax(pmin(-dZ, 20), -20),
         dZ_cat=factor(round(dZ, -1)),
         dXY=log1p(dXY)) %>%
  ggplot(aes(dXY, dZ_cat, colour=Time, fill=Time)) +
  ggridges::geom_density_ridges(alpha=0.25, scale=1) +
  facet_grid(particleDirection~meshRes*liceSpeedF)
pulse.sum %>%
  mutate(dZ=-dZ,
         dZ_cat=case_when(dZ < -20 ~ "downward > 20m",
                          between(dZ, -20, 0) ~ "downward 0-20m",
                          dZ > 0 ~ "upward"),
         dZ_cat=factor(dZ_cat, levels=rev(sort(unique(dZ_cat)))),
         dXY=log1p(dXY)) %>%
  group_by(meshRes, liceSpeedF, Time, particleDirection, dZ_cat) %>%
  summarise(dXY_mn=median(dXY, na.rm=T),
            dXY_lo1=HDInterval::hdi(dXY)[1],
            dXY_hi1=HDInterval::hdi(dXY)[2],
            dXY_lo2=HDInterval::hdi(dXY, 0.8)[1],
            dXY_hi2=HDInterval::hdi(dXY, 0.8)[2],
            dXY_lo3=HDInterval::hdi(dXY, 0.5)[1],
            dXY_hi3=HDInterval::hdi(dXY, 0.5)[2]) %>%
  ggplot(aes(particleDirection, dXY_mn, colour=meshRes, shape=Time)) +
  geom_point(position=position_dodge(width=0.5), size=2) +
  geom_linerange(aes(ymin=dXY_lo1, ymax=dXY_hi1), position=position_dodge(width=0.5), size=0.25) +
  geom_linerange(aes(ymin=dXY_lo2, ymax=dXY_hi2), position=position_dodge(width=0.5), size=0.5) +
  geom_linerange(aes(ymin=dXY_lo3, ymax=dXY_hi3), position=position_dodge(width=0.5), size=1) +
  scale_colour_manual(values=col_meshRes) +
  scale_shape_manual(values=c(1, 19)) +
  facet_grid(dZ_cat~liceSpeedF) +
  labs(x="Particle direction", y="\u0394ln(XY)") +
  theme(legend.position="bottom")
ggsave(glue("figs/init{initDensity}/pulse_lndXY_by_dZcat_hdi.png"), width=8.5, height=7, dpi=300)

pulse.sum %>%
  mutate(dZ=-dZ,
         dZ_cat=case_when(dZ < -20 ~ "downward > 20m",
                          between(dZ, -20, 0) ~ "downward 0-20m",
                          dZ > 0 ~ "upward"),
         dZ_cat=factor(dZ_cat, levels=rev(sort(unique(dZ_cat)))),
         dXY=log1p(dXY)) %>%
  group_by(meshRes, liceSpeedF, Time, particleDirection, dZ_cat) %>%
  summarise(dXY_mn=median(dXY, na.rm=T),
            dXY_lo1=HDInterval::hdi(dXY)[1],
            dXY_hi1=HDInterval::hdi(dXY)[2],
            dXY_lo2=HDInterval::hdi(dXY, 0.8)[1],
            dXY_hi2=HDInterval::hdi(dXY, 0.8)[2],
            dXY_lo3=HDInterval::hdi(dXY, 0.5)[1],
            dXY_hi3=HDInterval::hdi(dXY, 0.5)[2]) %>%
  ggplot(aes(Time, dXY_mn, colour=meshRes, fill=meshRes,
             shape=particleDirection,
             group=paste(meshRes, particleDirection))) +
  geom_point(position=position_dodge(width=0.5), size=2) +
  geom_errorbar(aes(ymin=dXY_lo2, ymax=dXY_hi2), size=0.25, width=0.2,
                position=position_dodge(width=0.5), alpha=0.4) +
  geom_errorbar(aes(ymin=dXY_lo3, ymax=dXY_hi3), size=1, width=0,
                position=position_dodge(width=0.5), alpha=0.4) +
  geom_line(position=position_dodge(width=0.5), size=1) +
  scale_colour_manual(values=col_meshRes) +
  scale_fill_manual(values=col_meshRes) +
  scale_linetype_manual(values=c(4,1)) +
  scale_shape_manual(values=c(1, 19)) +
  facet_grid(dZ_cat~liceSpeedF*particleDirection) +
  guides(shape=guide_legend(reverse=T),
         linetype=guide_legend(reverse=T)) +
  labs(x="Hour re: tidal pulse", y="\u0394ln(XY)") +
  theme(panel.spacing.x=unit(c(-0.1, 0.3, -0.1, 0.3, -0.1), "cm"))
ggsave(glue("figs/init{initDensity}/pulse_lndXY_by_dZcat_bumpPlot.png"), width=8, height=6, dpi=300)

pulse.sum %>%
  mutate(mnZ=if_else(Time=="before", z0 - (dZ/2), z0 + (dZ/2)),
         mnZ_cat=case_when(mnZ > 20 ~ "20m + ",
                          between(mnZ, 10, 20) ~ "10-20m",
                          mnZ < 10 ~ "0-10m"),
         mnZ_cat=factor(mnZ_cat, levels=(sort(unique(mnZ_cat)))),
         dXY=log1p(dXY)) %>%
  group_by(meshRes, liceSpeedF, Time, particleDirection, mnZ_cat) %>%
  summarise(dXY_mn=median(dXY, na.rm=T),
            dXY_lo1=HDInterval::hdi(dXY)[1],
            dXY_hi1=HDInterval::hdi(dXY)[2],
            dXY_lo2=HDInterval::hdi(dXY, 0.8)[1],
            dXY_hi2=HDInterval::hdi(dXY, 0.8)[2],
            dXY_lo3=HDInterval::hdi(dXY, 0.5)[1],
            dXY_hi3=HDInterval::hdi(dXY, 0.5)[2]) %>%
  ggplot(aes(particleDirection, dXY_mn, colour=meshRes, shape=Time)) +
  geom_point(position=position_dodge(width=0.5), size=2) +
  geom_linerange(aes(ymin=dXY_lo1, ymax=dXY_hi1), position=position_dodge(width=0.5), size=0.25) +
  geom_linerange(aes(ymin=dXY_lo2, ymax=dXY_hi2), position=position_dodge(width=0.5), size=0.5) +
  geom_linerange(aes(ymin=dXY_lo3, ymax=dXY_hi3), position=position_dodge(width=0.5), size=1) +
  scale_colour_manual(values=col_meshRes) +
  scale_shape_manual(values=c(1, 19)) +
  facet_grid(mnZ_cat~liceSpeedF) +
  labs(x="Particle direction", y="\u0394ln(XY)") +
  theme(legend.position="bottom")
ggsave(glue("figs/init{initDensity}/pulse_lndXY_by_Zcat_hdi.png"), width=8.5, height=7, dpi=300)

pulse.sum %>%
  mutate(mnZ=if_else(Time=="before", z0 - (dZ/2), z0 + (dZ/2)),
         mnZ_cat=case_when(mnZ > 20 ~ "20m + ",
                           between(mnZ, 10, 20) ~ "10-20m",
                           mnZ < 10 ~ "0-10m"),
         mnZ_cat=factor(mnZ_cat, levels=(sort(unique(mnZ_cat)))),
         dXY=log1p(dXY)) %>%
  group_by(meshRes, liceSpeedF, Time, particleDirection, mnZ_cat) %>%
  summarise(dXY_mn=median(dXY, na.rm=T),
            dXY_lo1=HDInterval::hdi(dXY)[1],
            dXY_hi1=HDInterval::hdi(dXY)[2],
            dXY_lo2=HDInterval::hdi(dXY, 0.8)[1],
            dXY_hi2=HDInterval::hdi(dXY, 0.8)[2],
            dXY_lo3=HDInterval::hdi(dXY, 0.5)[1],
            dXY_hi3=HDInterval::hdi(dXY, 0.5)[2]) %>%
  ggplot(aes(Time, dXY_mn, colour=meshRes, fill=meshRes,
             shape=particleDirection,
             group=paste(meshRes, particleDirection))) +
  geom_point(position=position_dodge(width=0.5), size=2) +
  geom_errorbar(aes(ymin=dXY_lo2, ymax=dXY_hi2), size=0.25, width=0.2,
                position=position_dodge(width=0.5), alpha=0.4) +
  geom_errorbar(aes(ymin=dXY_lo3, ymax=dXY_hi3), size=1, width=0,
                position=position_dodge(width=0.5), alpha=0.4) +
  geom_line(position=position_dodge(width=0.5), size=1) +
  scale_colour_manual(values=col_meshRes) +
  scale_fill_manual(values=col_meshRes) +
  scale_shape_manual(values=c(1, 19)) +
  facet_grid(mnZ_cat~liceSpeedF*particleDirection) +
  guides(shape=guide_legend(reverse=T),
         linetype=guide_legend(reverse=T)) +
  labs(x="Hour re: tidal pulse", y="\u0394ln(XY)") +
  theme(panel.spacing.x=unit(c(-0.1, 0.3, -0.1, 0.3, -0.1), "cm"))
ggsave(glue("figs/init{initDensity}/pulse_lndXY_by_Zcat_bumpPlot.png"), width=8, height=6, dpi=300)

pulse.sum %>%
  mutate(mnZ=if_else(Time=="before", z0 - (dZ/2), z0 + (dZ/2)),
         mnZ_cat=case_when(mnZ > 20 ~ "20m + ",
                           between(mnZ, 10, 20) ~ "10-20m",
                           mnZ < 10 ~ "0-10m"),
         mnZ_cat=factor(mnZ_cat, levels=(sort(unique(mnZ_cat)))),
         dXY=log1p(dXY)) %>%
  group_by(meshRes, liceSpeedF, Time, particleDirection, mnZ_cat) %>%
  summarise(dXY_mn=median(dXY, na.rm=T),
            dXY_lo1=HDInterval::hdi(dXY)[1],
            dXY_hi1=HDInterval::hdi(dXY)[2],
            dXY_lo2=HDInterval::hdi(dXY, 0.8)[1],
            dXY_hi2=HDInterval::hdi(dXY, 0.8)[2],
            dXY_lo3=HDInterval::hdi(dXY, 0.5)[1],
            dXY_hi3=HDInterval::hdi(dXY, 0.5)[2]) %>%
  ggplot(aes(Time, dXY_mn, colour=meshRes, fill=meshRes,
             linetype=particleDirection, shape=particleDirection,
             group=paste(meshRes, particleDirection))) +
  geom_point() +
  geom_line() +
  scale_colour_manual(values=col_meshRes) +
  scale_fill_manual(values=col_meshRes) +
  scale_linetype_manual(values=c(2,1)) +
  scale_shape_manual(values=c(1, 19)) +
  facet_grid(mnZ_cat~liceSpeedF*particleDirection) +
  guides(shape=guide_legend(reverse=T),
         linetype=guide_legend(reverse=T)) +
  labs(x="Hour re: tidal pulse", y="\u0394ln(XY)") +
  theme(panel.spacing.x=unit(c(-0.1, 0.3, -0.1, 0.3, -0.1), "cm"))
ggsave(glue("figs/init{initDensity}/pulse_lndXY_by_Zcat_bumpPlot2.png"), width=8, height=6, dpi=300)


ggplot(pulse.Z, aes(z, y=Time, colour=meshRes)) +
  ggridges::geom_density_ridges(fill=NA, scale=0.99) + 
  scale_colour_manual(values=col_meshRes) +
  facet_grid(.~liceSpeedF) +
  labs(x="Depth (m)", y="Time re: tidal pulse") +
  theme(legend.position="bottom", 
        panel.grid.major.y=element_line(colour="grey80", size=0.5)) 

ggplot(pulse.Z, aes(z, y=mesh, fill=Time)) +
  ggridges::geom_density_ridges(scale=0.99, alpha=0.5, colour="grey30") + 
  scale_fill_viridis_d("Time (h) re:\ntidal pulse", option="C", end=0.9) +
  facet_grid(timeRes~liceSpeedF) +
  xlim(0, 40) +
  labs(x="Depth (m)", y="Mesh") +
  theme(legend.position="bottom", 
        panel.grid.major.y=element_line(colour="grey80", size=0.5)) 
ggsave(glue("figs/init{initDensity}/pulse_z_densities.png"), width=8, height=6, dpi=300)

ggplot(pulse.sum, aes(dZ, y=mesh, fill=Time)) +
  ggridges::geom_density_ridges(scale=0.99, alpha=0.5, colour="grey30") + 
  scale_fill_viridis_d("Time re:\ntidal pulse", option="C", end=0.9) +
  facet_grid(timeRes~liceSpeedF) +
  # xlim(0, 40) +
  labs(x="\u0394 z (m)", y="Mesh") +
  theme(legend.position="bottom", 
        panel.grid.major.y=element_line(colour="grey80", size=0.5)) 

ggplot(pulse.Dens, aes(density, y=mesh, fill=Time)) +
  ggridges::geom_density_ridges(scale=0.99, alpha=0.5, colour="grey30") + 
  scale_fill_viridis_d("Time (h) re:\ntidal pulse", option="C", end=0.9) +
  facet_grid(timeRes~liceSpeedF) +
  labs(x="Density", y="Mesh") +
  theme(legend.position="bottom", 
        panel.grid.major.y=element_line(colour="grey80", size=0.5)) 
ggsave(glue("figs/init{initDensity}/pulse_density_densities.png"), width=8, height=6, dpi=300)

pulse.sum %>%
  group_by(meshRes, liceSpeedF, ID) %>%
  filter(!any(is.na(dDens))) %>%
  ggplot(aes(dDens, y=mesh, fill=Time)) +
  ggridges::geom_density_ridges(scale=0.99, alpha=0.5, colour="grey30") + 
  scale_fill_viridis_d("Time re:\ntidal pulse", option="C", end=0.9) +
  facet_grid(timeRes~liceSpeedF) +
  # xlim(0, 40) +
  labs(x="\u0394 logit(density)", y="Mesh") +
  theme(legend.position="bottom", 
        panel.grid.major.y=element_line(colour="grey80", size=0.5)) 


pulse.sf %>% st_drop_geometry() %>% 
  mutate(z_mid=(z + z_m1)/2,
         mnZ_cat=case_when(z_mid > 20 ~ "20m + ",
                           between(z_mid, 10, 20) ~ "10-20m",
                           z_mid < 10 ~ "0-10m"),
         mnZ_cat=factor(mnZ_cat, levels=(sort(unique(mnZ_cat))))) %>% 
  ggplot(aes(bearing_m, y=..ndensity.., fill=mnZ_cat)) + 
  geom_vline(xintercept=-2.225858) +
  geom_histogram(bins=24, colour="grey30") +
  annotate("text", x=pi*3/4, y=2.8, label="pre", size=3) +
  facet_grid(liceSpeedF~meshRes) + 
  scale_fill_viridis_d("Mean depth", direction=-1) +
  coord_polar(start=pi/2, direction=-1) +
  ylim(0, 2.8) +
  theme(panel.grid.major.y=element_line(colour="grey80"),
        axis.title.x=element_blank(),
        axis.text=element_blank()) +
  labs(y="density", title="Hour before pulse")
ggsave(glue("figs/init{initDensity}/pulse_windrose_-1.png"), width=8, height=5, dpi=600)

pulse.sf %>% st_drop_geometry() %>% 
  mutate(z_mid=(z + z_p1)/2,
         mnZ_cat=case_when(z_mid > 20 ~ "20m + ",
                           between(z_mid, 10, 20) ~ "10-20m",
                           z_mid < 10 ~ "0-10m"),
         mnZ_cat=factor(mnZ_cat, levels=(sort(unique(mnZ_cat))))) %>% 
  ggplot(aes(bearing_p, y=..ndensity.., fill=mnZ_cat)) + 
  geom_vline(xintercept=-2.225858) +
  geom_histogram(bins=24, colour="grey30") + 
  annotate("text", x=pi*3/4, y=2.8, label="post", size=3) +
  facet_grid(liceSpeedF~meshRes) + 
  scale_fill_viridis_d("Mean depth", direction=-1) +
  coord_polar(start=pi/2, direction=-1) +
  ylim(0, 2.8) +
  theme(panel.grid.major.y=element_line(colour="grey80"),
        axis.title.x=element_blank(),
        axis.text=element_blank()) +
  labs(y="density", title="Hour after pulse")
ggsave(glue("figs/init{initDensity}/pulse_windrose_+1.png"), width=8, height=5, dpi=600)
library(magick) 
list.files(path=glue("figs/init{initDensity}"), pattern="pulse_windrose.*png", full.names=TRUE) %>% 
  image_read() %>% # reads each path file
  image_join() %>% # joins image
  image_animate(fps=1) %>% # animates, can opt for number of loops
  image_write(glue("figs/init{initDensity}/pulse_windrose.gif")) # write




# mesh plots --------------------------------------------------------------

ggplot(mesh.fp) + geom_sf() +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.line=element_blank())
ggsave("figs/mesh/mesh_fp.png", width=3, height=3.5, dpi=300)

# ggplot(mesh.sf, aes(fill=bearing)) + 
#   geom_sf(colour=NA) + 
#   fill_NSEW + 
#   facet_grid(.~mesh)
# ggsave(glue("figs/mesh/mesh_bearing.png", width=6, height=4, dpi=300)
# ggplot(mesh.sf, aes(fill=depth)) + 
#   geom_sf(colour=NA) + 
#   scale_fill_viridis_c() +
#   facet_grid(.~mesh)
# ggsave(glue("figs/mesh/mesh_depth.png", width=6, height=4, dpi=300)



# vertical distribution ---------------------------------------------------

# tendency to be deeper in WeStCOMS, but need to filter by bathymetric depth
# since linnhe includes shallower areas
p <- loc.df %>% filter(age>=12, status != 66) %>%
  ggplot(aes(depth, colour=meshRes)) +
  geom_density() +
  scale_colour_manual(values=col_meshRes) +
  xlim(0, 25) +
  facet_grid(liceSpeedF~.)
ggsave(glue("figs/init{initDensity}/vertDist_site.png"), p, width=5, height=6, dpi=300)

# p <- loc.df %>% filter(age>=12 & status != 66) %>%
#   ggplot(aes(depth, colour=liceSpeedF)) +
#   geom_density() +
#   scale_colour_viridis_d() +
#   xlim(0, 25) +
#   facet_grid(meshRes~.)

# By hour
p <- loc.df %>% filter(age>=12 & status != 66) %>%
  ggplot(aes(depth, colour=hour, group=hour)) +
  geom_density() +
  xlim(0, 25) +
  col_hour +
  facet_grid(meshRes~liceSpeedF)
ggsave(glue("figs/init{initDensity}/vertDist_by_hour_site.png"), p, width=8, height=6, dpi=300)

p <- loc.df %>% filter(age>=12 & status != 66) %>%
  ggplot(aes(depth, colour=meshRes)) +
  geom_density() +
  xlim(0, 25) +
  scale_colour_manual(values=col_meshRes) +
  facet_grid(liceSpeedF~hour)
ggsave(glue("figs/init{initDensity}/vertDist_by_hour2_site.png"), p, width=17, height=5, dpi=300)

p <- loc.df %>% filter(age>=12 & status != 66) %>%
  ggplot(aes(depth, linetype=meshRes, colour=liceSpeedF, group=sim)) +
  geom_density() +
  xlim(0, 25) +
  scale_colour_viridis_d(end=0.9) +
  facet_wrap(~hour, ncol=6)
ggsave(glue("figs/init{initDensity}/vertDist_by_hour3_site.png"), p, width=12, height=9, dpi=300)





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
# p <- loc.df %>% 
#   group_by(meshRes, liceSpeedF, ID) %>% 
#   summarise(exited=any(meshParticle==1)) %>% 
#   group_by(meshRes, liceSpeedF) %>%
#   summarise(propExit=mean(exited)) %>% 
#   ggplot(aes(liceSpeedF, propExit, colour=meshRes, group=meshRes)) + 
#   geom_point() + geom_path() +
#   scale_colour_manual(values=col_meshRes) +
#   ylim(0, NA) + labs(y="Proportion leaving Loch Linnhe")
# ggsave(glue("figs/init{initDensity}/exits_site.png"), p, width=5, height=3, dpi=300)

# confirming exits occur at open boundaries only
# p <- loc.df %>% 
#   filter(meshParticle==1) %>%
#   ggplot(aes(x, y, colour=age)) + 
#   geom_point(alpha=0.2, size=0.5, shape=1) +
#   scale_colour_viridis_c() +
#   facet_grid(meshRes~liceSpeedF)

# Not sure...
# p <- loc.df %>%
#   filter(status==66) %>%
#   ggplot(aes(age, colour=liceSpeedF)) + 
#   geom_density() +
#   scale_colour_viridis_d() +
#   facet_grid(meshRes~.)

# Not sure... lice exit more slowly at an older age at high res
# p <- loc.df %>%
#   filter(status==66) %>%
#   ggplot(aes(age, colour=meshRes)) +
#   scale_colour_manual(values=col_meshRes) + 
#   geom_density() +
#   facet_grid(liceSpeedF~.)




# total distance ----------------------------------------------------------

p <- loc.df %>%
  filter(age>=12 & status != 66) %>%
  group_by(ID, sim, meshRes, liceSpeedF) %>%
  slice_tail(n=1) %>%
  ggplot(aes(xyTot/age, colour=meshRes, linetype=liceSpeedF, group=sim)) +
  geom_density() +
  scale_colour_manual(values=col_meshRes) +
  scale_x_log10("Mean xy hourly displacement (m/h)")
ggsave(glue("figs/init{initDensity}/xy_mean_log_site.png"), p, width=6, height=4, dpi=300)
p <- loc.df %>%
  filter(age>=12 & status!=66) %>%
  group_by(ID, sim, meshRes, liceSpeedF) %>%
  slice_tail(n=1) %>%
  ggplot(aes(xyTot/age, colour=meshRes, linetype=liceSpeedF, group=sim)) +
  geom_density() +
  scale_colour_manual(values=col_meshRes) +
  scale_x_continuous("Mean xy hourly displacement (m/h)")
ggsave(glue("figs/init{initDensity}/xy_mean_site.png"), p, width=6, height=4, dpi=300)

p <- loc.df %>%
  group_by(ID, sim, meshRes, liceSpeedF) %>%
  arrange(age) %>%
  mutate(dZ=depth-lag(depth)) %>%
  filter(age>=12 & status != 66) %>%
  ggplot(aes(abs(dZ), colour=meshRes, linetype=liceSpeedF, group=sim)) +
  geom_density() +
  scale_colour_manual(values=col_meshRes) +
  scale_x_continuous("Hourly depth displacement (abs(m))", trans="log1p")
ggsave(glue("figs/init{initDensity}/z_mean_log_site.png"), p, width=6, height=4, dpi=300)
p <- loc.df %>%
  group_by(ID, sim, meshRes, liceSpeedF) %>%
  arrange(age) %>%
  mutate(dZ=depth-lag(depth)) %>%
  filter(age>=12 & status != 66) %>%
  ggplot(aes(abs(dZ), colour=meshRes, linetype=liceSpeedF, group=sim)) +
  geom_density() +
  scale_colour_manual(values=col_meshRes) +
  scale_x_continuous("Hourly depth displacement (m)", limits=c(0,20))
ggsave(glue("figs/init{initDensity}/z_mean_site.png"), p, width=6, height=4, dpi=300)



# element activity --------------------------------------------------------

# p <- ggplot(elemAct.sf, aes(fill=sink/total)) + geom_sf(colour=NA) +
#   scale_fill_viridis_c(limits=c(0,1)) +
#   facet_grid(meshRes~liceSpeed)
# ggsave(glue("figs/init{initDensity}/prSink_mesh_by_speed_site.png"), p, width=8, height=8, dpi=300)
# 
# p <- ggplot(elemAct.sf, aes(fill=swim/total)) + geom_sf(colour=NA) +
#   scale_fill_viridis_c(limits=c(0,1)) +
#   facet_grid(meshRes~liceSpeed)
# ggsave(glue("figs/init{initDensity}/prSwim_mesh_by_speed_site.png"), p, width=8, height=8, dpi=300)
# 
# p <- ggplot(elemAct.sf, aes(fill=float/total)) + geom_sf(colour=NA) +
#   scale_fill_viridis_c(limits=c(0,1)) +
#   facet_grid(meshRes~liceSpeed)
# ggsave(glue("figs/init{initDensity}/prFloat_mesh_by_speed_site.png"), p, width=8, height=8, dpi=300)







# maps --------------------------------------------------------------------

# bbox <- st_bbox(mesh.fp)
# for(i in 1:nrow(sim_i)) {
#   p <- loc.df %>%
#     filter(age>=12 & status != 66) %>%
#     filter(sim==i) %>%
#     filter(between(x, bbox$xmin, bbox$xmax), between(y, bbox$ymin, bbox$ymax)) %>%
#     ggplot(aes(x,y, alpha=density, colour=depth)) +
#     geom_point(size=0.1, shape=1) +
#     scale_colour_viridis_c(direction=-1, option="D", limits=c(0,20), na.value="#000033") +
#     ggtitle(glue("{sim_i$mesh[i]}, {sim_i$timeRes[i]}, {sim_i$liceSpeedF[i]}"))
#   ggsave(glue("figs/init{initDensity}/part_locs_sim{str_pad(i,2,'left','0')}_site.png"),
#          p, width=6, height=6, dpi=900)
# }

bbox <- st_bbox(mesh.fp)
part.samp <- sample(unique(loc.df$ID), 5000)
for(i in 1:nrow(sim_i)) {
  p <- loc.df %>%
    filter(status != 66) %>%
    filter(sim==i) %>%
    filter(between(x, bbox$xmin, bbox$xmax), between(y, bbox$ymin, bbox$ymax)) %>%
    arrange(sim, ID, age) %>%
    ggplot(aes(x, y, alpha=density, group=ID)) +
    geom_path() +
    scale_colour_viridis_c(direction=-1, option="D", limits=c(0,20), na.value="#000033") +
    ggtitle(glue("{sim_i$mesh[i]}, {sim_i$timeRes[i]}, {sim_i$liceSpeedF[i]}"))
  ggsave(glue("figs/init{initDensity}/part_paths_sim{str_pad(i,2,'left','0')}_site.png"),
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






# density maps ------------------------------------------------------------

# dens.sf <- loc.df %>% 
#   filter(status==2) %>%
#   filter(depth < 30) %>%
#   select(sim, meshRes, liceSpeedF, ID, x, y, density) %>%
#   st_as_sf(coords=c("x", "y"), crs=27700) %>%
#   st_join(mesh.sf %>% filter(mesh=="WeStCOMS2") %>% select(i)) %>%
#   st_drop_geometry() %>%
#   group_by(sim, meshRes, liceSpeedF, i) %>%
#   summarise(part_density=sum(density), N=n()) %>%
#   ungroup %>%
#   inner_join(mesh.sf %>% filter(mesh=="WeStCOMS2"), .) %>%
#   mutate(depth_30=pmin(depth, 30),
#          dens_m3=part_density/(area*depth_30),
#          ln_densm3=log(dens_m3),
#          lnDens_m3=log(part_density)/(area*depth_30))
# st_write(dens.sf, "temp/density_30m_initScaled.gpkg", append=F)
# 
# p <- ggplot(dens.sf) + 
#   geom_sf(data=mesh.fp) +
#   geom_sf(aes(fill=ln_densm3), colour=NA) + 
#   scale_fill_viridis_c(expression("ln(lice/m"^3~")"), option="magma") +
#   facet_grid(liceSpeedF~meshRes) +
#   theme(axis.text=element_blank())
# ggsave(glue("figs/init{initDensity}/density_map_site.png"), p, width=9, height=8, dpi=300)
# 
# p <- ggplot(dens.sf) + 
#   geom_sf(data=mesh.fp) +
#   geom_sf(aes(fill=dens_m3), colour=NA) + 
#   scale_fill_viridis_c(expression("lice/m"^3), option="magma") +
#   facet_grid(liceSpeedF~meshRes) +
#   theme(axis.text=element_blank())
# ggsave(glue("figs/init{initDensity}/density2_map_site.png"), p, width=9, height=8, dpi=300)
# 
# p <- ggplot(dens.sf) + 
#   geom_sf(data=mesh.fp) +
#   geom_sf(aes(fill=log(N/(area*depth_30))), colour=NA) + 
#   scale_fill_viridis_c(expression("ln(particles/m"^3~")"), option="magma") +
#   facet_grid(liceSpeedF~meshRes) +
#   theme(axis.text=element_blank())
# ggsave(glue("figs/init{initDensity}/density3_map_site.png"), p, width=9, height=8, dpi=300)
# 




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
# ggsave(glue("figs/init{initDensity}/direction_map_sim4_site.png"), p, width=8, height=5, dpi=300)
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
# ggsave(glue("figs/init{initDensity}/direction_map_sim3_site.png"), p, width=8, height=5, dpi=300)

p <- velocity.df %>%
  filter(ID %% 2000 == 0) %>%
  ggplot() +
  geom_sf(data=mesh.fp, fill=NA, colour="grey30") +
  geom_path(aes(x, y, group=ID, colour=downloch), alpha=0.25, size=0.1) + 
  col_downloch +
  facet_grid(meshRes~liceSpeedF)
ggsave(glue("figs/init{initDensity}/direction_map_site.png"), p, width=9, height=9, dpi=300)



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
ggsave(glue("figs/init{initDensity}/density_map_site.png"), p, width=13, height=5, dpi=300)
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
ggsave(glue("figs/init{initDensity}/density_diff_map_site.png"), p, width=11, height=5, dpi=300)
p <- meshDiff.df %>%
  ggplot() +
  geom_sf(data=mesh.fp, colour="grey30", fill=NA) +
  geom_sf(colour=NA, aes(fill=lnN_diff/lnN_ref)) + 
  scale_fill_gradient2() +
  facet_grid(.~meshRes)
ggsave(glue("figs/init{initDensity}/densityN_diff_map_site.png"), p, width=11, height=5, dpi=300)
p <- meshDiff.df %>%
  ggplot() +
  geom_sf(data=mesh.fp, colour="grey30", fill=NA) +
  geom_sf(colour=NA, aes(fill=prUp_diff)) + 
  scale_fill_gradient2() +
  facet_grid(.~meshRes)
ggsave(glue("figs/init{initDensity}/prUploch_diff_map_site.png"), p, width=11, height=5, dpi=300)
p <- ggplot(meshDir.sum %>% filter(liceSpeedF=="Medium") %>% filter(N>10)) +
  geom_sf(data=mesh.fp, colour="grey30", fill=NA) +
  geom_sf(colour=NA, aes(fill=prUploch)) + 
  fill_downloch +
  facet_grid(.~meshRes)
ggsave(glue("figs/init{initDensity}/prUploch_map_site.png"), p, width=13, height=5, dpi=300)
p <- ggplot(meshDir.sum %>% filter(liceSpeedF=="Medium") %>% filter(N>10)) +
  geom_sf(data=mesh.fp, colour="grey30", fill=NA) +
  geom_sf(colour=NA, aes(fill=prUploch, alpha=lnDens_m3)) + 
  fill_downloch +
  facet_grid(.~meshRes)
ggsave(glue("figs/init{initDensity}/prUploch_map_ALT_site.png"), p, width=13, height=5, dpi=300)


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
ggsave(glue("figs/init{initDensity}/density_map_depthStrat_site.png"), p, width=11, height=14, dpi=300)
p <- ggplot(meshDir.deep %>% filter(!is.na(deep), N>3)) +
  geom_sf(data=mesh.fp, colour="grey30", fill=NA) +
  geom_sf(colour=NA, aes(fill=prUploch)) + 
  fill_downloch +
  facet_grid(deep~meshRes)
ggsave(glue("figs/init{initDensity}/prUploch_map_depthStrat_site.png"), p, width=11, height=14, dpi=300)
p <- ggplot(meshDir.deep %>% filter(!is.na(deep), N>3)) +
  geom_sf(data=mesh.fp, colour="grey30", fill=NA) +
  geom_sf(colour=NA, aes(fill=prUploch, alpha=lnDens_m3)) + 
  fill_downloch +
  facet_grid(deep~meshRes)
ggsave(glue("figs/init{initDensity}/prUploch_map_depthStrat_ALT_site.png"), p, width=11, height=14, dpi=300)






