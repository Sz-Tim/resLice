# Run simulations
# Off Aqua
# Tim Szewczyk


# This script prepares directories and properties files, then sequentially runs
# the specified simulations.




# setup -------------------------------------------------------------------

library(tidyverse); library(glue); library(lubridate); library(sf)
source("code/00_fn.R")



# define parameters -------------------------------------------------------

overwrite_jar <- T
cores <- 8

dirs <- switch(get_os(),
               windows=list(proj=getwd(),
                            mesh="D:/hydroOut/",
                            hydro.linnhe="D:/hydroOut/linnhe7/",
                            hydro.westcoms="D:/hydroOut/WeStCOMS2/Archive/",
                            jar="C:/Users/sa04ts/OneDrive - SAMS/Projects/OffAqua/particle_track/out/",
                            out=glue("{getwd()}/out/")),
               linux=list(proj=getwd(),
                          mesh="/home/sa04ts/FVCOM_meshes",
                          hydro.linnhe="/media/archiver/common/sa04ts-temp/linnhe7/",
                          hydro.westcoms="/media/archiver/common/sa01da-work/WeStCOMS2/Archive/",
                          jar=glue("{getwd()}/jar/"),
                          out=glue("{getwd()}/out/")))

sim.i <- expand_grid(mesh=c("WeStCOMS2", "linnhe7"),
                     timeRes=c("1h"),
                     liceSpeed=c(0, 0.0005)) %>%
  mutate(i=str_pad(row_number(), 2, "left", "0"),
         meshFile=if_else(mesh=="WeStCOMS2", 
                          glue("{dirs$mesh}WeStCOMS2_linnhe_mesh.nc"),
                          glue("{dirs$mesh}linnhe_mesh.nc")),
         hydroDir=if_else(mesh=="WeStCOMS2",
                          glue("{dirs$hydro.westcoms}"),
                          glue("{dirs$hydro.linnhe}linnhe7_tides_met_tsobc")),
         hydroDir=if_else(timeRes=="5min", glue("{hydroDir}_5min"), hydroDir),
         outDir=glue("{dirs$out}/sim_{i}/"),
         nDays=if_else(timeRes=="1h", 2, 85),
         dt=if_else(timeRes=="1h", 3600, 337.5),
         releaseInterval=if_else(timeRes=="1h", 100, 32),
         viabletime=if_else(timeRes=="1h", 12, 128),
         maxParticleAge=if_else(timeRes=="1h", 500, 5000))
write_csv(sim.i, glue("{dirs$out}/sim_i.csv"))  
sim_seq <- 1:nrow(sim.i)



# set properties ----------------------------------------------------------

walk(sim_seq, ~dir.create(sim.i$outDir[.x], showWarnings=F))

sep <- ifelse(get_os()=="windows", "\\", "/")
properties.ls <- map(
  sim_seq,
  ~setPartTrackProperties(parallelThreads=cores,
                          destinationDirectory=paste0(normalizePath(sim.i$outDir[.x]), sep),
                          datadir=paste0(normalizePath(sim.i$hydroDir[.x]), sep),
                          mesh1=paste0(normalizePath(sim.i$meshFile[.x]), sep),
                          location=str_to_lower(str_sub(sim.i$mesh[.x], 1, -2)),
                          minchVersion=str_sub(sim.i$mesh[.x], -1, -1),
                          sitefile="..\\..\\data\\linnhe_start_100m_corran_20km.tsv",
                          numberOfDays=sim.i$nDays[.x],
                          dt=sim.i$dt[.x],
                          releaseInterval=sim.i$releaseInterval[.x],
                          nparts=5,
                          viabletime=sim.i$viabletime[.x],
                          maxParticleAge=sim.i$maxParticleAge[.x],
                          vertSwimSpeedMean=sim.i$liceSpeed[.x],
                          vertSwimSpeedStd=sim.i$liceSpeed[.x]/5,
                          sinkingRateMean=sim.i$liceSpeed[.x],
                          sinkingRateStd=sim.i$liceSpeed[.x]/5,
                          recordMovement=F))
walk(sim_seq, 
     ~cat(properties.ls[[.x]] %>% 
            str_replace_all("\\\\", "\\\\\\\\") %>%
            str_replace_all("\\ ", "\\\\\\\\ "), 
          "\n", file=glue("{dirs$out}/sim_{sim.i$i[.x]}.properties")))




# make local copies -------------------------------------------------------

file.copy(glue("{dirs$jar}/particle_track.jar"), "jar/particle_track.jar", 
          overwrite=overwrite_jar)
file.copy(if_else(get_os()=="windows", "code/runSims_pc.sh", "code/runSims_smn.sh"), 
          glue("{dirs$out}/run_local.sh"), overwrite=T)
walk(sim_seq, ~file.copy("jar/lib", sim.i$outDir[.x], recursive=T, overwrite=T))



# run simulations ---------------------------------------------------------

for(i in sim_seq) {
  setwd(dirs$out)
  system2("bash", c("run_local.sh", 
                    glue("sim_{str_pad(i, 2, 'left', '0')}.properties"),
                    glue("sim_{str_pad(i, 2, 'left', '0')}/")))
  setwd(dirs$proj)
}

