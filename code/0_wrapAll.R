


# Run and process site release
source("code/05_runSiteReleases.R")
rm(list=ls()); gc()
source("code/06_processSiteSimulations.R")
rm(list=ls()); gc()

# Run and process general release
source("code/01_runSimulations.R")
rm(list=ls()); gc()
source("code/02_processSimulations.R")
rm(list=ls()); gc()



# Make figures
source("code/0000_temp_connectivityEDA.R")
rm(list=ls()); gc()
source("code/03_outputEDA.R")
rm(list=ls()); gc()
source("code/07_outputSiteEDA.R")
rm(list=ls()); gc()


# Make animations
source("code/04_animTracks.R")
rm(list=ls()); gc()
source("code/08_animSiteTracks.R")
rm(list=ls()); gc()




