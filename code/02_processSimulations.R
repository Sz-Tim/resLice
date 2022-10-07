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

sim_i <- read_csv(glue("{dirs$out}/sim_i.csv"))
sim_seq <- 1:nrow(sim_i)
elemAct.df <- map_dfr(sim_seq, ~glue("{sim_i$outDir[.x]}/elementActivity.dat") %>% 
                        read_delim(col_names=c("sink", "swim", "float", "total")) %>%
                        mutate(total=sink+swim+float,
                               i=row_number(),
                               sim=.x))


# elemActivity <- read_delim(glue("{out.dir}elementActivity.dat"), 
#                            col_names=c("sink", "swim", "float", "total")) %>%
#   mutate(total=sink+swim+float,
#          i=row_number())