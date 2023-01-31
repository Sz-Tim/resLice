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
                            mesh="D:/hydroOut/",
                            out=glue("{getwd()}/out/siteRelease_init{initDensity}/")),
               linux=list(proj=getwd(),
                          mesh="/home/sa04ts/FVCOM_meshes",
                          out=glue("{getwd()}/out/siteRelease_init{initDensity}/")))




# load files --------------------------------------------------------------

mesh.sf <- list(linnhe7=st_read(glue("{dirs$mesh}/linnhe_mesh.gpkg")) %>%
                  mutate(mesh="linnhe7"),
                westcoms2=st_read(glue("{dirs$mesh}/WeStCOMS2_linnhe_mesh.gpkg")) %>%
                  mutate(mesh="WeStCOMS2")) %>%
  do.call('rbind', .) %>%
  select(i, area, depth, mesh, geom)
sim_i <- read_csv(glue("{dirs$out}/sim_i.csv")) %>%
  mutate(liceSpeedF=factor(liceSpeed, levels=c(0.0001, 0.0005, 0.001), 
                           labels=c("Slow", "Medium", "Fast")))
sim_seq <- 1:nrow(sim_i)
time.key <- read_csv("data/timeRes_key.csv") %>% 
  rename(fileHour=layer) %>%
  mutate(fileDate=as.character(fileDate)) %>%
  select(mesh, timeRes, fileDate, fileHour, timeCalculated)




# element activity --------------------------------------------------------

elemAct.df <- map_dfr(sim_seq,
                      ~glue("{sim_i$outDir[.x]}/elementActivity.dat") %>% 
                        read_delim(col_names=c("sink", "swim", "float", "total")) %>%
                        mutate(total=sink+swim+float,
                               i=row_number(),
                               sim=.x)) %>%
  left_join(sim_i %>% mutate(sim=as.numeric(i)) %>%
              select(mesh, timeRes, liceSpeed, liceSpeedF, sim))
elemAct.sf <- left_join(mesh.sf, elemAct.df, by=c("mesh", "i"))
write_sf(elemAct.sf, glue("out/00_processed/elementActivity_site_init{initDensity}.gpkg"))
rm(elemAct.df); rm(elemAct.sf); gc()



# particle locations ------------------------------------------------------

# 5 min files are too big. Need to process using temp storage, then combine.
for(i in sim_seq) {
  dir(glue("{sim_i$outDir[i]}locs"), "locations_", full.names=T) %>%
    map_dfr(~read_delim(.x, delim=" ", col_names=T, show_col_types=F) %>%
              select(-startLocation, -startDate, -depthLayer, -degreeDays) %>%
              rename(meshParticle=mesh) %>%
              mutate(date=str_sub(.x,-12,-5))) %>%
    mutate(sim=i) %>%
    left_join(sim_i %>% mutate(sim=as.numeric(i)) %>%
                select(mesh, timeRes, sim)) %>%
    left_join(time.key, by=c("mesh", "timeRes", 
                             "date"="fileDate", "hour"="fileHour")) %>%
    filter(minute(timeCalculated)==0) %>%
    select(-mesh, -timeRes, -hour, -date) %>%
    saveRDS(glue("out/00_processed/temp_site_init{initDensity}_sim_{i}.rds"))
  gc()
}
map_dfr(dir("out/00_processed", glue("temp_site_init{initDensity}_sim.*.rds"), full.names=T), readRDS) %>%
  left_join(sim_i %>% mutate(sim=as.numeric(i)) %>%
              select(sim, mesh, timeRes, liceSpeedF)) %>%
  saveRDS(glue("out/00_processed/locations_site_init{initDensity}.rds"))






# connectivity ------------------------------------------------------------

site.i <- read_tsv("data/fishFarmSites.tsv", col_names=c("site", "x", "y"))
connect.df <- vector("list", length(sim_seq))
for(i in sim_seq) {
  connect.df[[i]] <- dir(glue("{sim_i$outDir[i]}connectivity"), "connectivity_") %>%
    map_dfr(~read_delim(glue("{sim_i$outDir[i]}connectivity/{.x}"), 
                        delim=" ", col_names=site.i$site, show_col_types=F) %>%
              mutate(source=site.i$site) %>%
              pivot_longer(-"source", names_to="dest", values_to="density") %>%
              mutate(sim=i,
                     fdate=str_split_fixed(.x, "_", 3)[,2],
                     steps=str_remove(str_split_fixed(.x, "_", 3)[,3], ".dat"))) %>%
    mutate(#density=density/ifelse(sim_i$timeRes[i]=="1h", 1, 12),
           cumulHours=as.numeric(steps)/if_else(sim_i$timeRes[i]=="1h", 24, 24),
           timeSim=as_datetime("2021-11-01 00:00:00") + cumulHours*60*60) %>%
    arrange(timeSim)
}
do.call('rbind', connect.df) %>%
  left_join(sim_i %>% mutate(sim=as.numeric(i)) %>%
              select(sim, mesh, timeRes, liceSpeedF)) %>%
  saveRDS(glue("out/00_processed/connectivity_site_init{initDensity}.rds"))


