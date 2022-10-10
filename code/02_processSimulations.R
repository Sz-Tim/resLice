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
sim_i <- read_csv(glue("{dirs$out}/sim_i.csv")) %>%
  mutate(liceSpeedF=factor(liceSpeed, labels=c("Passive", "Slow", "Medium", "Fast")))
sim_seq <- 1:nrow(sim_i)




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
write_sf(elemAct.sf, "out/00_processed/elementActivity.gpkg")



# particle locations ------------------------------------------------------


loc.df <- map_dfr(
  sim_seq, 
  ~dir(glue("{sim_i$outDir[.x]}locs"), "locations_", full.names=T) %>%
    map_dfr(~read_delim(.x, delim=" ", col_names=T, show_col_types=F) %>%
              select(-mesh) %>%
              mutate(date=ymd(str_sub(.x,-12,-5)),
                     dateTime=as_datetime(paste0(str_sub(.x,-12,-5), "-", hour),
                                          format="%Y%m%d-%H"))) %>%
    mutate(sim=.x)
  ) %>%
  left_join(sim_i %>% mutate(sim=as.numeric(i)) %>%
              select(mesh, timeRes, liceSpeed, liceSpeedF, sim))
write_csv(loc.df, "out/00_processed/locations.csv")







