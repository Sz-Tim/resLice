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
sims <- dir(dirs$out, "sim_[0-9][0-9]")

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



connect.df <- do.call('rbind', connect.ls) %>%
  left_join(sim_i %>% rename(sim=paste0("sim_", i)) %>%
              select(sim, mesh, timeRes, liceSpeedF))

connect.df <- readRDS(glue("out/00_processed/connectivity_siteReps_init{initDensity}.rds")) %>%
  mutate(meshRes=paste(mesh, timeRes, sep=", "))

theme_set(theme_bw())

connect.df %>%
  group_by(sim, meshRes, liceSpeedF, timeSim, rep) %>%
  summarise(density=sum(density)) %>%
  ggplot(aes(timeSim, density, group=paste0(liceSpeedF, rep), colour=liceSpeedF)) + 
  geom_line(alpha=0.5) + 
  facet_wrap(~meshRes) +
  ggtitle("Total infections") +
  guides(colour=guide_legend(override.aes=list(alpha=1)))
connect.df %>%
  group_by(sim, meshRes, liceSpeedF, timeSim, rep) %>%
  summarise(density=sum(density)) %>%
  ggplot(aes(timeSim, density, group=paste0(meshRes, rep), colour=meshRes)) + 
  geom_line(alpha=0.5) + 
  facet_wrap(~liceSpeedF) +
  ggtitle("Total infections") +
  guides(colour=guide_legend(override.aes=list(alpha=1)))

connect.df %>%
  group_by(sim, mesh, timeRes, liceSpeedF, timeSim, source, dest) %>%
  summarise(density_mn=mean(density),
            density_lo=quantile(density, probs=0.1),
            density_hi=quantile(density, probs=0.9)) %>%
  ggplot(aes(timeSim, density_mn, ymin=density_lo, ymax=density_hi, colour=liceSpeedF, fill=liceSpeedF)) + 
  geom_line() + geom_ribbon(colour=NA, alpha=0.5) +
  facet_grid(source~dest) +
  ggtitle("Pairwise infections")

connect.df %>%
  filter(source!=dest) %>%
  group_by(sim, mesh, timeRes, liceSpeedF, timeSim, source, rep) %>%
  summarise(density=sum(density)) %>%
  group_by(sim, mesh, timeRes, liceSpeedF, timeSim, source) %>%
  summarise(density_mn=mean(density),
            density_lo=quantile(density, probs=0.1),
            density_hi=quantile(density, probs=0.9)) %>%
  ggplot(aes(timeSim, density_mn, ymin=density_lo, ymax=density_hi, colour=liceSpeedF, fill=liceSpeedF)) + 
  geom_line() + geom_ribbon(colour=NA, alpha=0.5) +
  facet_wrap(~source) +
  ggtitle("Outflux to other sites")

connect.df %>%
  filter(source!=dest) %>%
  group_by(sim, mesh, timeRes, liceSpeedF, timeSim, dest, rep) %>%
  summarise(density=sum(density)) %>%
  group_by(sim, mesh, timeRes, liceSpeedF, timeSim, dest) %>%
  summarise(density_mn=mean(density),
            density_lo=quantile(density, probs=0.1),
            density_hi=quantile(density, probs=0.9)) %>%
  ggplot(aes(timeSim, density_mn, ymin=density_lo, ymax=density_hi, colour=liceSpeedF, fill=liceSpeedF)) + 
  geom_line() + geom_ribbon(colour=NA, alpha=0.5) +
  facet_wrap(~dest) +
  ggtitle("Influx from other sites")

connect.df %>%
  filter(source==dest) %>%
  group_by(sim, mesh, timeRes, liceSpeedF, timeSim, dest, rep) %>%
  summarise(density=sum(density)) %>%
  group_by(sim, mesh, timeRes, liceSpeedF, timeSim, dest) %>%
  summarise(density_mn=mean(density),
            density_lo=quantile(density, probs=0.1),
            density_hi=quantile(density, probs=0.9)) %>%
  ggplot(aes(timeSim, density_mn, ymin=density_lo, ymax=density_hi, colour=liceSpeedF, fill=liceSpeedF)) + 
  geom_line() + geom_ribbon(colour=NA, alpha=0.5) +
  facet_wrap(~dest) +
  ggtitle("Self-infection")

connect.df %>%
  group_by(sim, mesh, timeRes, liceSpeedF, timeSim, rep) %>%
  summarise(density=sum(density)) %>%
  group_by(sim, mesh, timeRes, liceSpeedF, timeSim) %>%
  summarise(density_mn=mean(density),
            density_lo=quantile(density, probs=0.1),
            density_hi=quantile(density, probs=0.9)) %>%
  ggplot(aes(timeSim, density_mn, ymin=density_lo, ymax=density_hi, colour=sim, fill=sim)) + 
  geom_line() + geom_ribbon(colour=NA, alpha=0.5)
