# Process output
# Off Aqua
# Tim Szewczyk


# This script processes and summarises the simulation output


# setup -------------------------------------------------------------------

library(tidyverse); library(glue); library(lubridate); library(sf)
source("code/00_fn.R")

dirs <- switch(get_os(),
               windows=list(proj=getwd(),
                            mesh="D:/hydroOut/",
                            out=glue("{getwd()}/out/gridRelease/")),
               linux=list(proj=getwd(),
                          mesh="/home/sa04ts/FVCOM_meshes",
                          out=glue("{getwd()}/out/gridRelease/")))




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
write_sf(elemAct.sf, "out/00_processed/elementActivity_grid.gpkg")



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
    saveRDS(glue("out/00_processed/temp_grid_sim_{i}.rds"))
  gc()
}
map_dfr(dir("out/00_processed", "temp_grid_sim.*.rds", full.names=T), readRDS) %>%
  left_join(sim_i %>% mutate(sim=as.numeric(i)) %>%
              select(sim, mesh, timeRes, liceSpeedF)) %>%
  saveRDS("out/00_processed/locations_grid.rds")







