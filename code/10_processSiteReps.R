# Process output
# Off Aqua
# Tim Szewczyk


# This script processes and summarises the simulation output


# setup -------------------------------------------------------------------

library(tidyverse); library(glue); library(lubridate); library(sf)
source("code/00_fn.R")

initDensity <- c("Scaled", "Uniform")[2]

dirs <- switch(get_os(),
               windows=list(proj=getwd(),
                            out=glue("{getwd()}/out/siteConnectReps_init{initDensity}/")),
               linux=list(proj=getwd(),
                          out=glue("{getwd()}/out/siteConnectReps_init{initDensity}/")))




# load files --------------------------------------------------------------

sim_i <- read_csv(glue("{dirs$out}/sim_i.csv")) %>%
  mutate(liceSpeedF=factor(liceSpeed, levels=c(0.0001, 0.0005, 0.001), 
                           labels=c("Slow", "Medium", "Fast")))
sims <- dir(dirs$out, "sim_[0-9][0-9]$")
time.key <- read_csv("data/timeRes_key.csv") %>% 
  rename(fileHour=layer) %>%
  mutate(fileDate=as.character(fileDate)) %>%
  select(mesh, timeRes, fileDate, fileHour, timeCalculated)




# connectivity ------------------------------------------------------------

site.i <- read_tsv("data/fishFarmSites_all3yr.tsv", col_names=c("site", "x", "y"))
connect.ls <- map(sims, ~vector("list"))
for(inum in seq_along(sims)) {
  i <- sims[inum]
  sim_num <- as.numeric(str_sub(i, -2, -1))
  rep_seq <- dir(glue("{dirs$out}/{i}"), pattern="[0-9][0-9]")
  for(jnum in seq_along(rep_seq)) {
    j <- rep_seq[jnum]
    dir.ij <- glue("{dirs$out}/{i}/{j}/connectivity/")
    connect.f <- dir(glue("{dir.ij}"))
    connect.ls[[inum]][[jnum]] <- map_dfr(
      connect.f,
      ~read_delim(glue("{dir.ij}/{.x}"), 
                  delim=" ", col_names=site.i$site, show_col_types=F) %>%
        mutate(source=site.i$site) %>%
        pivot_longer(-"source", names_to="dest", values_to="density") %>%
        mutate(sim=i,
               rep=j,
               fdate=str_split_fixed(.x, "_", 3)[,2],
               steps=str_remove(str_split_fixed(.x, "_", 3)[,3], ".dat"))) %>%
      mutate(cumulHours=as.numeric(steps)/if_else(sim_i$timeRes[sim_num]=="1h", 24, 24),
             timeSim=as_datetime("2021-11-01 00:00:00") + cumulHours*60*60) %>%
      arrange(timeSim)
  }
  connect.ls[[inum]] <- do.call('rbind', connect.ls[[inum]])
}
do.call('rbind', connect.ls) %>%
  left_join(sim_i %>% mutate(sim=paste0("sim_", i)) %>%
              select(sim, mesh, timeRes, liceSpeedF)) %>%
  saveRDS(glue("out/00_processed/connectivity_siteReps_init{initDensity}.rds"))










# TO MOVE -----------------------------------------------------------------

connect.df <- readRDS(glue("out/00_processed/connectivity_siteReps_init{initDensity}.rds")) %>%
  mutate(meshRes=paste(mesh, timeRes, sep=", "),
         sim=as.numeric(str_sub(sim, -2, -1)))

theme_set(theme_bw())
col_meshRes <- c("linnhe7, 1h"="#a6611a", "linnhe7, 5min"="#dfc27d",
                 "WeStCOMS2, 1h"="#018571", "WeStCOMS2, 5min"="#80cdc1")
tot_part <- 16*(24*7-65)*5


connect.df %>%
  filter(timeSim > "2021-11-03 00:00:00") %>%
  group_by(sim, meshRes, liceSpeedF, timeSim, rep) %>%
  summarise(density=sum(density)) %>%
  ggplot(aes(timeSim, density, group=paste0(sim, rep), colour=meshRes)) + 
  geom_line(alpha=0.5) + 
  scale_colour_manual(values=col_meshRes) +
  ggtitle("Total infections") +
  facet_wrap(~liceSpeedF)

connect.df %>%
  filter(timeSim > "2021-11-03 00:00:00") %>%
  group_by(sim, meshRes, liceSpeedF, timeSim, source, dest) %>%
  summarise(density_mn=mean(density),
            density_lo=quantile(density, probs=0.1),
            density_hi=quantile(density, probs=0.9)) %>%
  ggplot(aes(timeSim, density_mn, ymin=density_lo, ymax=density_hi, colour=meshRes, fill=meshRes)) + 
  geom_line() + geom_ribbon(colour=NA, alpha=0.5) +
  facet_grid(source~dest) +
  scale_colour_manual(values=col_meshRes) +
  scale_fill_manual(values=col_meshRes) +
  ggtitle("Pairwise infections")

connect.df %>%
  filter(timeSim > "2021-11-03 00:00:00") %>%
  filter(source!=dest) %>%
  group_by(sim, meshRes, liceSpeedF, timeSim, source, rep) %>%
  summarise(density=sum(density)) %>%
  group_by(sim, meshRes, liceSpeedF, timeSim, source) %>%
  summarise(density_mn=mean(density),
            density_lo=quantile(density, probs=0.1),
            density_hi=quantile(density, probs=0.9)) %>%
  ggplot(aes(timeSim, density_mn, ymin=density_lo, ymax=density_hi, colour=meshRes, fill=meshRes)) + 
  geom_line() + geom_ribbon(colour=NA, alpha=0.5) +
  facet_grid(liceSpeedF~source) +
  scale_colour_manual(values=col_meshRes) +
  scale_fill_manual(values=col_meshRes) +
  ggtitle("Outflux to other sites")

connect.df %>%
  filter(timeSim > "2021-11-03 00:00:00") %>%
  group_by(sim, meshRes, liceSpeedF, timeSim, source, rep) %>%
  summarise(density=sum(density)) %>%
  group_by(sim, meshRes, liceSpeedF, timeSim, source) %>%
  summarise(density_mn=mean(density),
            density_lo=quantile(density, probs=0.1),
            density_hi=quantile(density, probs=0.9)) %>%
  ggplot(aes(timeSim, density_mn, ymin=density_lo, ymax=density_hi, colour=meshRes, fill=meshRes)) + 
  geom_line() + geom_ribbon(colour=NA, alpha=0.5) +
  facet_grid(liceSpeedF~source) +
  scale_colour_manual(values=col_meshRes) +
  scale_fill_manual(values=col_meshRes) +
  ggtitle("Outflux + self: Pressure caused")

connect.df %>%
  filter(timeSim > "2021-11-03 00:00:00") %>%
  filter(source!=dest) %>%
  group_by(sim, meshRes, liceSpeedF, timeSim, dest, rep) %>%
  summarise(density=sum(density)) %>%
  group_by(sim, meshRes, liceSpeedF, timeSim, dest) %>%
  summarise(density_mn=mean(density),
            density_lo=quantile(density, probs=0.1),
            density_hi=quantile(density, probs=0.9)) %>%
  ggplot(aes(timeSim, density_mn, ymin=density_lo, ymax=density_hi, colour=meshRes, fill=meshRes)) + 
  geom_line() + geom_ribbon(colour=NA, alpha=0.5) +
  facet_grid(liceSpeedF~dest) +
  scale_colour_manual(values=col_meshRes) +
  scale_fill_manual(values=col_meshRes) +
  ggtitle("Influx from other sites")

connect.df %>%
  filter(timeSim > "2021-11-03 00:00:00") %>%
  group_by(sim, meshRes, liceSpeedF, timeSim, dest, rep) %>%
  summarise(density=sum(density)) %>%
  group_by(sim, meshRes, liceSpeedF, timeSim, dest) %>%
  summarise(density_mn=mean(density),
            density_lo=quantile(density, probs=0.1),
            density_hi=quantile(density, probs=0.9)) %>%
  ggplot(aes(timeSim, density_mn, ymin=density_lo, ymax=density_hi, colour=meshRes, fill=meshRes)) + 
  geom_line() + geom_ribbon(colour=NA, alpha=0.5) +
  facet_grid(liceSpeedF~dest) +
  scale_colour_manual(values=col_meshRes) +
  scale_fill_manual(values=col_meshRes) +
  ggtitle("Influx + self: Pressure received")

connect.df %>%
  filter(timeSim > "2021-11-03 00:00:00") %>%
  filter(source==dest) %>%
  group_by(sim, meshRes, liceSpeedF, timeSim, dest, rep) %>%
  summarise(density=sum(density)) %>%
  group_by(sim, meshRes, liceSpeedF, timeSim, dest) %>%
  summarise(density_mn=mean(density),
            density_lo=quantile(density, probs=0.1),
            density_hi=quantile(density, probs=0.9)) %>%
  ggplot(aes(timeSim, density_mn, ymin=density_lo, ymax=density_hi, colour=meshRes, fill=meshRes)) + 
  geom_line() + geom_ribbon(colour=NA, alpha=0.5) +
  facet_grid(liceSpeedF~dest) +
  scale_colour_manual(values=col_meshRes) +
  scale_fill_manual(values=col_meshRes) +
  ggtitle("Self infection")

connect.df %>%
  filter(timeSim > "2021-11-03 00:00:00") %>%
  group_by(sim, meshRes, liceSpeedF, timeSim, rep) %>%
  summarise(density=sum(density)) %>%
  group_by(sim, meshRes, liceSpeedF, timeSim) %>%
  summarise(density_mn=mean(density),
            density_lo=quantile(density, probs=0.1),
            density_hi=quantile(density, probs=0.9)) %>%
  ggplot(aes(timeSim, density_mn, ymin=density_lo, ymax=density_hi, colour=meshRes, fill=meshRes)) + 
  geom_line() + geom_ribbon(colour=NA, alpha=0.5) + 
  facet_wrap(~liceSpeedF) +
  scale_colour_manual(values=col_meshRes) +
  scale_fill_manual(values=col_meshRes) +
  labs(x="Date", y="Total infection pressure")
connect.df %>%
  filter(timeSim > "2021-11-03 00:00:00") %>%
  group_by(sim, meshRes, liceSpeedF, timeSim, rep) %>%
  summarise(density=sum(density)) %>%
  group_by(sim, meshRes, liceSpeedF, timeSim) %>%
  summarise(density_mn=mean(density),
            density_lo=quantile(density, probs=0.1),
            density_hi=quantile(density, probs=0.9)) %>%
  ggplot(aes(timeSim, density_mn, ymin=density_lo, ymax=density_hi, colour=liceSpeedF, fill=liceSpeedF)) + 
  geom_line() + geom_ribbon(colour=NA, alpha=0.5) + 
  facet_wrap(~meshRes) +
  scale_colour_viridis_d("Lice speed", end=0.9) +
  scale_fill_viridis_d("Lice speed", end=0.9) +
  labs(x="Date", y="Total infection pressure")
  


connect.df %>%
  mutate(density=density/tot_part) %>%
  group_by(sim, meshRes, liceSpeedF, rep) %>%
  summarise(density=sum(density)) %>%
  group_by(sim, meshRes, liceSpeedF) %>%
  summarise(density_mn=mean(density),
            density_lo=quantile(density, probs=0.1),
            density_hi=quantile(density, probs=0.9)) %>%
  ggplot(aes(liceSpeedF, density_mn, ymin=density_lo, ymax=density_hi,
             colour=meshRes, fill=meshRes, group=meshRes)) + 
  geom_point() + geom_line() + geom_linerange() + geom_ribbon(alpha=0.5, colour=NA) +
  scale_colour_manual(values=col_meshRes) +
  scale_fill_manual(values=col_meshRes) +
  labs(x="Lice speed", y="Pr(Infection | particle)")
connect.df %>%
  mutate(density=density/tot_part) %>%
  group_by(sim, meshRes, liceSpeedF, rep) %>%
  summarise(density=sum(density)) %>%
  group_by(sim, meshRes, liceSpeedF) %>%
  summarise(density_mn=mean(density),
            density_lo=quantile(density, probs=0.1),
            density_hi=quantile(density, probs=0.9)) %>%
  ggplot(aes(meshRes, density_mn, ymin=density_lo, ymax=density_hi,
             colour=liceSpeedF, fill=liceSpeedF, group=liceSpeedF)) + 
  geom_point() + geom_line() + geom_linerange() + geom_ribbon(alpha=0.5, colour=NA) +
  scale_colour_viridis_d("Lice speed", end=0.9) +
  scale_fill_viridis_d("Lice speed", end=0.9) +
  labs(x="Mesh", y="Pr(Infection | particle)")


connect.df %>%
  group_by(sim, meshRes, liceSpeedF, source, rep) %>%
  summarise(density=sum(density)) %>%
  group_by(sim, meshRes, liceSpeedF, source) %>%
  summarise(density_mn=mean(density),
            density_lo=quantile(density, probs=0.1),
            density_hi=quantile(density, probs=0.9)) %>%
  ggplot(aes(source, density_mn, ymin=density_lo, ymax=density_hi, colour=meshRes,
             fill=meshRes, group=meshRes)) + 
  geom_point() + geom_line() + geom_linerange() + geom_ribbon(alpha=0.5, colour=NA) +
  facet_grid(liceSpeedF~.) +
  ggtitle("Outflux + self: Pressure caused")

connect.df %>%
  group_by(sim, meshRes, liceSpeedF, source, rep) %>%
  summarise(density=sum(density)) %>%
  group_by(sim, meshRes, liceSpeedF, source) %>%
  summarise(density_mn=mean(density),
            density_lo=quantile(density, probs=0.1),
            density_hi=quantile(density, probs=0.9)) %>%
  ungroup %>%
  mutate(sortOrder=factor(sim, levels=c(8,1:7,9:12))) %>%
  arrange(sortOrder, density_mn) %>%
  mutate(source=factor(source, levels=unique(source))) %>%
  ggplot(aes(source, density_mn, ymin=density_lo, ymax=density_hi, colour=meshRes,
             fill=meshRes, group=meshRes)) + 
  geom_point() + geom_line() + geom_linerange() + geom_ribbon(alpha=0.5, colour=NA) +
  facet_grid(liceSpeedF~.) +
  ggtitle("Outflux + self: Pressure caused")

connect.df %>%
  group_by(sim, meshRes, liceSpeedF, dest, rep) %>%
  summarise(density=sum(density)) %>%
  group_by(sim, meshRes, liceSpeedF, dest) %>%
  summarise(density_mn=mean(density),
            density_lo=quantile(density, probs=0.1),
            density_hi=quantile(density, probs=0.9)) %>%
  ungroup %>%
  mutate(sortOrder=factor(sim, levels=c(8,1:7,9:12))) %>%
  arrange(sortOrder, density_mn) %>%
  mutate(dest=factor(dest, levels=unique(dest))) %>%
  ggplot(aes(dest, density_mn, ymin=density_lo, ymax=density_hi, colour=meshRes,
             fill=meshRes, group=meshRes)) + 
  geom_point() + geom_line() + geom_linerange() + geom_ribbon(alpha=0.5, colour=NA) +
  facet_grid(liceSpeedF~.) +
  ggtitle("Influx + self: Pressure received")

connect.df %>%
  group_by(sim, meshRes, liceSpeedF, dest, rep) %>%
  summarise(density=sum(density)) %>%
  group_by(sim, meshRes, liceSpeedF, dest) %>%
  summarise(density_mn=mean(density),
            density_lo=quantile(density, probs=0.1),
            density_hi=quantile(density, probs=0.9)) %>%
  ungroup %>%
  mutate(sortOrder=factor(sim, levels=c(8,1:7,9:12))) %>%
  arrange(sortOrder, density_mn) %>%
  mutate(dest=factor(dest, levels=unique(dest))) %>%
  group_by(dest) %>%
  mutate(density_delta=density_mn-first(density_mn),
         density_pct_chg=density_delta/first(density_mn)) %>%
  ggplot(aes(dest, density_delta, colour=meshRes,
             fill=meshRes, group=sim)) + 
  geom_point() + geom_line() +
  facet_grid(liceSpeedF~.) +
  ggtitle("Influx + self: Pressure received") + 
  coord_flip()
connect.df %>%
  group_by(sim, meshRes, liceSpeedF, source, rep) %>%
  summarise(density=sum(density)) %>%
  group_by(sim, meshRes, liceSpeedF, source) %>%
  summarise(density_mn=mean(density),
            density_lo=quantile(density, probs=0.1),
            density_hi=quantile(density, probs=0.9)) %>%
  ungroup %>%
  mutate(sortOrder=factor(sim, levels=c(8,1:7,9:12))) %>%
  arrange(sortOrder, density_mn) %>%
  mutate(source=factor(source, levels=unique(source))) %>%
  group_by(source) %>%
  mutate(density_delta=density_mn-first(density_mn),
         density_pct_chg=density_delta/first(density_mn)) %>%
  ggplot(aes(source, density_delta, colour=meshRes,
             fill=meshRes, group=sim)) + 
  geom_point() + geom_line() +
  facet_grid(liceSpeedF~.) +
  ggtitle("Outflux + self: Pressure caused") + 
  coord_flip()

connect.df %>%
  group_by(sim, meshRes, liceSpeedF, dest, rep) %>%
  summarise(density=sum(density)) %>%
  group_by(sim, meshRes, liceSpeedF, dest) %>%
  summarise(density_mn=mean(density),
            density_lo=quantile(density, probs=0.1),
            density_hi=quantile(density, probs=0.9)) %>%
  ungroup %>%
  mutate(sortOrder=factor(sim, levels=c(8,1:7,9:12))) %>%
  arrange(sortOrder, density_mn) %>%
  mutate(dest=factor(dest, levels=unique(dest))) %>%
  group_by(dest) %>%
  mutate(density_delta=density_mn-first(density_mn),
         density_pct_chg=density_delta/first(density_mn)) %>%
  ggplot(aes(dest, density_delta, colour=liceSpeedF,
             fill=liceSpeedF, group=sim)) + 
  geom_point() + geom_line() +
  facet_grid(meshRes~.) +
  ggtitle("Influx + self: Pressure received") + 
  coord_flip()
connect.df %>%
  group_by(sim, meshRes, liceSpeedF, source, rep) %>%
  summarise(density=sum(density)) %>%
  group_by(sim, meshRes, liceSpeedF, source) %>%
  summarise(density_mn=mean(density),
            density_lo=quantile(density, probs=0.1),
            density_hi=quantile(density, probs=0.9)) %>%
  ungroup %>%
  mutate(sortOrder=factor(sim, levels=c(8,1:7,9:12))) %>%
  arrange(sortOrder, density_mn) %>%
  mutate(source=factor(source, levels=unique(source))) %>%
  group_by(source) %>%
  mutate(density_delta=density_mn-first(density_mn),
         density_pct_chg=density_delta/first(density_mn)) %>%
  ggplot(aes(source, density_delta, colour=liceSpeedF,
             fill=liceSpeedF, group=sim)) + 
  geom_point() + geom_line() +
  facet_grid(meshRes~.) +
  ggtitle("Outflux + self: Pressure caused") + 
  coord_flip()

connect.df %>%
  group_by(sim, meshRes, liceSpeedF, source, rep) %>%
  summarise(density=sum(density)) %>%
  group_by(sim, meshRes, liceSpeedF, source) %>%
  summarise(density_mn=mean(density),
            density_lo=quantile(density, probs=0.1),
            density_hi=quantile(density, probs=0.9)) %>%
  ungroup %>%
  mutate(sortOrder=factor(sim, levels=c(8,1:7,9:12))) %>%
  arrange(sortOrder, density_mn) %>%
  mutate(source=factor(source, levels=unique(source))) %>%
  group_by(source) %>%
  mutate(density_delta=density_mn-first(density_mn),
         density_pct_chg=density_delta/first(density_mn)) %>%
  filter(sim!=8) %>%
  ggplot(aes(density_delta, colour=meshRes,
             group=sim, shape=liceSpeedF, linetype=liceSpeedF)) + 
  geom_vline(xintercept=0) +
  geom_density() +
  ggtitle("Outflux + self: Pressure caused") 


mesh.fp <- st_read("data/linnhe_mesh_footprint.gpkg")
site.sf <- site.i %>% st_as_sf(coords=c("x", "y"), crs=27700)

influxPct.sf <- connect.df %>%
  group_by(sim, meshRes, liceSpeedF, dest, rep) %>%
  summarise(density=sum(density)) %>%
  group_by(sim, meshRes, liceSpeedF, dest) %>%
  summarise(density_mn=mean(density),
            density_lo=quantile(density, probs=0.1),
            density_hi=quantile(density, probs=0.9)) %>%
  ungroup %>%
  mutate(sortOrder=factor(sim, levels=c(8,1:7,9:12))) %>%
  arrange(sortOrder, density_mn) %>%
  mutate(dest=factor(dest, levels=unique(dest))) %>%
  group_by(dest) %>%
  mutate(density_delta=density_mn-first(density_mn),
         density_pct_chg=density_delta/first(density_mn)) %>%
  ungroup %>%
  left_join(site.sf, ., by=c("site"="dest"))

influxNoSelfPct.sf <- connect.df %>%
  filter(source!=dest) %>%
  group_by(sim, meshRes, liceSpeedF, dest, rep) %>%
  summarise(density=sum(density)) %>%
  group_by(sim, meshRes, liceSpeedF, dest) %>%
  summarise(density_mn=mean(density),
            density_lo=quantile(density, probs=0.1),
            density_hi=quantile(density, probs=0.9)) %>%
  ungroup %>%
  mutate(sortOrder=factor(sim, levels=c(8,1:7,9:12))) %>%
  arrange(sortOrder, density_mn) %>%
  mutate(dest=factor(dest, levels=unique(dest))) %>%
  group_by(dest) %>%
  mutate(density_delta=density_mn-first(density_mn),
         density_pct_chg=density_delta/first(density_mn)) %>%
  ungroup %>%
  left_join(site.sf, ., by=c("site"="dest"))

outfluxPct.sf <- connect.df %>%
  group_by(sim, meshRes, liceSpeedF, source, rep) %>%
  summarise(density=sum(density)) %>%
  group_by(sim, meshRes, liceSpeedF, source) %>%
  summarise(density_mn=mean(density),
            density_lo=quantile(density, probs=0.1),
            density_hi=quantile(density, probs=0.9)) %>%
  ungroup %>%
  mutate(sortOrder=factor(sim, levels=c(8,1:7,9:12))) %>%
  arrange(sortOrder, density_mn) %>%
  mutate(source=factor(source, levels=unique(source))) %>%
  group_by(source) %>%
  mutate(density_delta=density_mn-first(density_mn),
         density_pct_chg=density_delta/first(density_mn)) %>%
  ungroup %>%
  left_join(site.sf, ., by=c("site"="source"))

outfluxNoSelfPct.sf <- connect.df %>%
  filter(source!=dest) %>%
  group_by(sim, meshRes, liceSpeedF, source, rep) %>%
  summarise(density=sum(density)) %>%
  group_by(sim, meshRes, liceSpeedF, source) %>%
  summarise(density_mn=mean(density),
            density_lo=quantile(density, probs=0.1),
            density_hi=quantile(density, probs=0.9)) %>%
  ungroup %>%
  mutate(sortOrder=factor(sim, levels=c(8,1:7,9:12))) %>%
  arrange(sortOrder, density_mn) %>%
  mutate(source=factor(source, levels=unique(source))) %>%
  group_by(source) %>%
  mutate(density_delta=density_mn-first(density_mn),
         density_pct_chg=density_delta/first(density_mn)) %>%
  ungroup %>%
  left_join(site.sf, ., by=c("site"="source"))

ggplot(influxPct.sf) + geom_sf(data=mesh.fp) +
  geom_sf(aes(colour=density_pct_chg), size=4) +
  scale_colour_gradient2(low=scales::muted("blue"), high=scales::muted("red")) +
  facet_grid(liceSpeedF~meshRes) +
  ggtitle("Influx+self: \u0394 vs. WeStCOMS 1h (proportion)")
ggplot(influxNoSelfPct.sf) + geom_sf(data=mesh.fp) +
  geom_sf(aes(colour=density_pct_chg), size=4) +
  scale_colour_gradient2(low=scales::muted("blue"), high=scales::muted("red")) +
  facet_grid(liceSpeedF~meshRes) +
  ggtitle("Influx: \u0394 vs. WeStCOMS 1h (proportion)")
ggplot(outfluxPct.sf) + geom_sf(data=mesh.fp) +
  geom_sf(aes(colour=density_pct_chg), size=4) +
  scale_colour_gradient2(low=scales::muted("blue"), high=scales::muted("red")) +
  facet_grid(liceSpeedF~meshRes) +
  ggtitle("Outflux+self: \u0394 vs. WeStCOMS 1h (proportion)")
ggplot(outfluxNoSelfPct.sf) + geom_sf(data=mesh.fp) +
  geom_sf(aes(colour=density_pct_chg), size=4) +
  scale_colour_gradient2(low=scales::muted("blue"), high=scales::muted("red")) +
  facet_grid(liceSpeedF~meshRes) +
  ggtitle("Outflux: \u0394 vs. WeStCOMS 1h (proportion)")

ggplot(influxPct.sf) + geom_sf(data=mesh.fp) +
  geom_sf(aes(colour=density_delta), size=4) +
  scale_colour_gradient2(low=scales::muted("blue"), high=scales::muted("red")) +
  facet_grid(liceSpeedF~meshRes) +
  ggtitle("Influx+self: \u0394 vs. WeStCOMS 1h")
ggplot(influxNoSelfPct.sf) + geom_sf(data=mesh.fp) +
  geom_sf(aes(colour=density_delta), size=4) +
  scale_colour_gradient2(low=scales::muted("blue"), high=scales::muted("red")) +
  facet_grid(liceSpeedF~meshRes) +
  ggtitle("Influx: \u0394 vs. WeStCOMS 1h")
ggplot(outfluxPct.sf) + geom_sf(data=mesh.fp) +
  geom_sf(aes(colour=density_delta), size=4) +
  scale_colour_gradient2(low=scales::muted("blue"), high=scales::muted("red")) +
  facet_grid(liceSpeedF~meshRes) +
  ggtitle("Outflux+self: \u0394 vs. WeStCOMS 1h")
ggplot(outfluxNoSelfPct.sf) + geom_sf(data=mesh.fp) +
  geom_sf(aes(colour=density_delta), size=4) +
  scale_colour_gradient2(low=scales::muted("blue"), high=scales::muted("red")) +
  facet_grid(liceSpeedF~meshRes) +
  ggtitle("Outflux: \u0394 vs. WeStCOMS 1h")









# mesh fp -----------------------------------------------------------------

wc.fp <- st_read("../../WeStCOMS/data/WeStCOMS2_mesh.gpkg") %>% st_union()

ggplot(wc.fp) + 
  geom_sf(colour=NA, fill="#377eb8") + 
  geom_sf(data=mesh.fp, fill="#ff7f00", colour="#ff7f00", size=0.2) + 
  theme_classic() +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(), 
        axis.line=element_blank())
ggsave("figs/mesh_overview.png", height=18, width=8, dpi=300)
