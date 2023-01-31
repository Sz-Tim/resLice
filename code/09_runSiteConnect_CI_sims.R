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
cores <- ifelse(get_os()=="windows", 11, 50)
nDays <- 7
initDensity <- c("Scaled", "Uniform")[2]
nReps <- 20
nParts_hr_site <- 20

dirs <- switch(get_os(),
               windows=list(proj=getwd(),
                            mesh="D:/hydroOut",
                            hydro.linnhe="D:/hydroOut/linnhe7/linnhe7_tides_met_tsobc_riv",
                            hydro.westcoms="D:/hydroOut/WeStCOMS2/Archive",
                            jar="C:/Users/sa04ts/OneDrive - SAMS/Projects/OffAqua/particle_track/out",
                            out=glue("{getwd()}/out/siteConnectReps_init{initDensity}")),
               linux=list(proj=getwd(),
                          mesh="/home/sa04ts/FVCOM_meshes",
                          hydro.linnhe="/media/archiver/common/sa04ts-temp/linnhe7",
                          hydro.westcoms="/media/archiver/common/sa04ts-temp/WeStCOMS2",
                          jar=glue("/home/sa04ts/biotracker/"),
                          out=glue("{getwd()}/out/siteConnectReps_init{initDensity}")))

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
         siteDensity=glue("liceScale_daily_{timeRes}.tsv"),
         nDays=if_else(timeRes=="1h", nDays, nDays*12),
         dt=if_else(timeRes=="1h", 3600, 300),
         stepsPerStep=if_else(timeRes=="1h", 24, 2),
         releaseInterval=1,
         viableDegreeDays=40,
         maxDegreeDays=150)
if(initDensity=="Uniform") {
  sim.i$siteDensity <- ""
}
write_csv(sim.i, glue("{dirs$out}/sim_i.csv"))  
sim_seq <- 1:nrow(sim.i)
rep_seq <- str_pad(1:nReps, 2, "l", "0")


# set properties ----------------------------------------------------------

sep <- ifelse(get_os()=="windows", "\\", "/")

properties.ls <- vector("list", length(sim_seq)*length(rep_seq))
ij <- 1

for(i in sim_seq) {
  dir.create(sim.i$outDir[i], showWarnings=F)
  for(j in rep_seq) {
    dir.create(glue("{sim.i$outDir[i]}/{j}"), showWarnings=F)
    
    properties.ls[[ij]] <- setPartTrackProperties(
      parallelThreads=cores,
      destinationDirectory=paste0(normalizePath(paste0(sim.i$outDir[i], "/", j)), sep),
      datadir=paste0(normalizePath(sim.i$hydroDir[i]), sep),
      mesh1=paste0(normalizePath(sim.i$meshFile[i]), sep),
      location=str_to_lower(str_sub(sim.i$mesh[i], 1, -2)),
      minchVersion=str_sub(sim.i$mesh[i], -1, -1),
      datadir2=paste0(normalizePath(sim.i$hydroDir2[i]), sep),
      mesh2=normalizePath(paste0(dirs$mesh, "/WeStCOMS2_mesh.nc")),
      sitefile="..\\..\\..\\..\\data\\fishFarmSites.tsv",
      siteDensityPath=ifelse(sim.i$siteDensity[i]=="", "",
                             paste0("..\\..\\..\\..\\data\\", sim.i$siteDensity[i])), 
      numberOfDays=sim.i$nDays[i],
      dt=sim.i$dt[i],
      stepsPerStep=sim.i$stepsPerStep[i],
      releaseInterval=sim.i$releaseInterval[i],
      nparts=nParts_hr_site,
      viableDegreeDays=sim.i$viableDegreeDays[i],
      maxDegreeDays=sim.i$maxDegreeDays[i],
      vertSwimSpeedMean=sim.i$liceSpeed[i],
      vertSwimSpeedStd=sim.i$liceSpeed[i]/5,
      sinkingRateMean=sim.i$liceSpeed[i],
      sinkingRateStd=sim.i$liceSpeed[i]/5,
      variableDiffusion="false",
      recordMovement="false",
      recordElemActivity="false",
      recordLocations="false",
      recordConnectivity="true",
      connectivityInterval=ifelse(sim.i$timeRes[i]=="1h", 1, 12),
      verboseSetUp="true"
    )
    if(get_os()=="windows") {
      cat(properties.ls[[ij]] %>% 
            str_replace_all("\\\\", "\\\\\\\\") %>%
            str_replace_all("\\ ", "\\\\\\\\ "), 
          "\n", file=glue("{dirs$out}/sim_{sim.i$i[i]}_{j}.properties"))
    } else {
      cat(properties.ls[[ij]] %>% 
            str_replace_all("\\\\", "//") %>%
            str_replace_all("\\ ", "\\\\\\\\ "), 
          "\n", file=glue("{dirs$out}/sim_{sim.i$i[i]}_{j}.properties"))
    }
    
    
    ij <- ij + 1 
  }
}



# make local copies -------------------------------------------------------

file.copy(glue("{dirs$jar}/particle_track.jar"), "jar/particle_track.jar", 
          overwrite=overwrite_jar)
file.copy(if_else(get_os()=="windows", "code/runSimReps_pc.sh", "code/runSimReps_smn.sh"), 
          glue("{dirs$out}/run_local.sh"), overwrite=T)
for(i in sim_seq) {
  for(j in rep_seq) {
    file.copy("jar/lib", paste0(sim.i$outDir[i], j), recursive=T, overwrite=T)
  }
}
walk(sim_seq, ~file.copy("jar/lib", sim.i$outDir[.x], recursive=T, overwrite=T))



# run simulations ---------------------------------------------------------

for(i in sim_seq) {
  for(j in rep_seq) {
    setwd(dirs$out)
    system2("bash", c("run_local.sh", 
                      glue("sim_{str_pad(i, 2, 'left', '0')}_{j}.properties"),
                      glue("sim_{str_pad(i, 2, 'left', '0')}/{j}")))
    setwd(dirs$proj)
  }
}

