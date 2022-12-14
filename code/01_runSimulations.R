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
cores <- 11
nDays <- 7

dirs <- switch(get_os(),
               windows=list(proj=getwd(),
                            mesh="D:/hydroOut",
                            hydro.linnhe="D:/hydroOut/linnhe7/linnhe7_tides_met_tsobc_riv",
                            hydro.westcoms="D:/hydroOut/WeStCOMS2/Archive",
                            jar="C:/Users/sa04ts/OneDrive - SAMS/Projects/OffAqua/particle_track/out",
                            out=glue("{getwd()}/out/gridRelease")),
               linux=list(proj=getwd(),
                          mesh="/home/sa04ts/FVCOM_meshes",
                          hydro.linnhe="/media/archiver/common/sa04ts-temp/linnhe7",
                          hydro.westcoms="/media/archiver/common/sa01da-work/WeStCOMS2/Archive",
                          jar=glue("{getwd()}/jar"),
                          out=glue("{getwd()}/out/gridRelease")))

sim.i <- expand_grid(mesh=c("linnhe7", "WeStCOMS2"),
                     timeRes=c("1h", "5min"),
                     liceSpeed=c(0.0001, 0.0005, 0.001)) %>%
  mutate(i=str_pad(row_number(), 2, "left", "0"),
         meshFile=if_else(mesh=="WeStCOMS2", 
                          glue("{dirs$mesh}/WeStCOMS2_linnhe_mesh.nc"),
                          glue("{dirs$mesh}/linnhe_mesh.nc")),
         hydroDir=if_else(mesh=="WeStCOMS2", dirs$hydro.westcoms, dirs$hydro.linnhe),
         hydroDir=glue("{hydroDir}{ifelse(timeRes=='5min','_5min','')}"),
         hydroDir2=glue("{dirs$hydro.westcoms}{ifelse(timeRes=='5min','_5min','')}"),
         outDir=glue("{dirs$out}/sim_{i}/"),
         nDays=if_else(timeRes=="1h", nDays, nDays*12),
         dt=if_else(timeRes=="1h", 3600, 300),
         stepsPerStep=if_else(timeRes=="1h", 24, 2),
         releaseInterval=100,
         viabletime=if_else(timeRes=="1h", 100, 1200),
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
                          datadir2=paste0(normalizePath(sim.i$hydroDir2[.x]), sep),
                          mesh2=normalizePath(paste0(dirs$mesh, "/WeStCOMS2_mesh.nc")),
                          sitefile="..\\..\\..\\data\\linnhe_start_100m_corran_20km.tsv",
                          numberOfDays=sim.i$nDays[.x],
                          dt=sim.i$dt[.x],
                          stepsPerStep=sim.i$stepsPerStep[.x],
                          releaseInterval=sim.i$releaseInterval[.x],
                          nparts=3,
                          viabletime=sim.i$viabletime[.x],
                          maxParticleAge=sim.i$maxParticleAge[.x],
                          vertSwimSpeedMean=sim.i$liceSpeed[.x],
                          vertSwimSpeedStd=sim.i$liceSpeed[.x]/5,
                          sinkingRateMean=sim.i$liceSpeed[.x],
                          sinkingRateStd=sim.i$liceSpeed[.x]/5,
                          variableDiffusion="false",
                          recordMovement="false"))
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

