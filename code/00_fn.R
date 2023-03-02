# Functions
# Off Aqua
# Tim Szewczyk




# simulation set up -------------------------------------------------------

get_os <- function() {
  # https://www.r-bloggers.com/2015/06/identifying-the-os-from-r/
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}




extract_surface <- function(nc, dir, subdir, meshRes, westcoms_i=NULL) {
  nc.i <- nc_open(glue("{dir}/{subdir}/{nc}"))
  u <- ncvar_get(nc.i, "u")[,1,]
  v <- ncvar_get(nc.i, "v")[,1,]
  ww <- ncvar_get(nc.i, "ww")[,1,]
  temp <- ncvar_get(nc.i, "temp")[,1,]
  salinity <- ncvar_get(nc.i, "salinity")[,1,]
  nc_close(nc.i)
  
  currentSpeed <- sqrt(u^2 + v^2 + ww^2)
  elem.df <- expand_grid(fileHour=(1:ncol(u))-1, i=1:nrow(u)) %>% 
    mutate(fileDate=str_split_fixed(nc, "_", 3)[1,2],
           meshRes=meshRes,
           currentSpeed=c(currentSpeed),
           vertSpeed=c(ww))
  node.df <- expand_grid(fileHour=(1:ncol(temp))-1, i=1:nrow(temp)) %>% 
    mutate(fileDate=str_split_fixed(nc, "_", 3)[1,2],
           meshRes=meshRes,
           temp=c(temp),
           salinity=c(salinity))
  if(!is.null(westcoms_i)) {
    elem.df <- filter(elem.df, i %in% westcoms_i$elem)
    node.df <- filter(node.df, i %in% westcoms_i$node)
  }
  return(list(elem=elem.df, node=node.df))
}


setPartTrackProperties <- function(
  destinationDirectory="out/",
  coordRef="OSGB1936",
  datadir="/media/archiver/common/sa01da-work/WeStCOMS2/Archive/",
  datadirPrefix="netcdf_",
  datadirSuffix="",
  mesh1="/home/sa04ts/FVCOM_meshes/WeStCOMS2_linnhe_mesh.nc",
  mesh1Type="FVCOM",
  location="westcoms",
  minchVersion=2,
  datadir2="/media/archiver/common/sa01da-work/WeStCOMS2/Archive/",
  datadir2Prefix="netcdf_",
  datadir2Suffix="",
  mesh2="/home/sa04ts/FVCOM_meshes/WeStCOMS2_mesh.nc",
  mesh2Type="FVCOM",
  location2="westcoms",
  minchVersion2=2,
  sitefile="../../data/linnhe_start_1000m_grid.tsv",
  siteDensityPath="",
  verboseSetUp="false",
  start_ymd=20211101,
  numberOfDays=7,
  checkOpenBoundaries="true",
  openBoundaryThresh=500,
  duplicateLastDay="true",
  recordsPerFile1=25,
  dt=3600,
  verticalDynamics="true",
  maxDepth=10000,
  parallelThreads=4,
  releaseScenario=1,
  releaseInterval=3,
  nparts=1,
  setStartDepth="true",
  startDepth=1,
  stepsPerStep=30,
  variableDiffusion="false",
  D_h=0.1,
  D_hVert=0.001,
  salinityThreshMin=23,
  salinityThreshMax=31,
  mortalityRate=0.01,
  salinityMort="true",
  swimLightLevel="true",
  vertSwimSpeedMean=-0.0005,
  vertSwimSpeedStd=0.0001,
  sinkingRateMean=0.0005,
  sinkingRateStd=0.0001,
  viabletime=12,
  maxParticleAge=500,
  viableDegreeDays=-1,
  maxDegreeDays=-1,
  recordPsteps="false",
  splitPsteps="false",
  pstepsInterval=1,
  recordMovement="true",
  recordElemActivity="true",
  recordConnectivity="false",
  connectivityInterval=24,
  recordLocations="true",
  recordArrivals="false"
) {
  params <- c(
    destinationDirectory=destinationDirectory,
    coordRef=coordRef,
    datadir=datadir,
    datadirPrefix=datadirPrefix,
    datadirSuffix=datadirSuffix,
    mesh1=mesh1,
    mesh1Type=mesh1Type,
    location=location,
    minchVersion=minchVersion,
    datadir2=datadir2,
    datadir2Prefix=datadir2Prefix,
    datadir2Suffix=datadir2Suffix,
    mesh2=mesh2,
    mesh2Type=mesh2Type,
    location2=location2,
    minchVersion2=minchVersion2,
    sitefile=sitefile,
    siteDensityPath=siteDensityPath,
    verboseSetUp=verboseSetUp,
    start_ymd=start_ymd,
    numberOfDays=numberOfDays,
    checkOpenBoundaries=checkOpenBoundaries,
    openBoundaryThresh=openBoundaryThresh,
    duplicateLastDay=duplicateLastDay,
    recordsPerFile1=recordsPerFile1,
    dt=dt,
    verticalDynamics=verticalDynamics,
    maxDepth=maxDepth,
    parallelThreads=parallelThreads,
    releaseScenario=releaseScenario,
    releaseInterval=releaseInterval,
    nparts=nparts,
    setStartDepth=setStartDepth,
    startDepth=startDepth,
    stepsPerStep=stepsPerStep,
    variableDiffusion=variableDiffusion,
    D_h=D_h,
    D_hVert=D_hVert,
    salinityThreshMin=salinityThreshMin,
    salinityThreshMax=salinityThreshMax,
    mortalityRate=mortalityRate,
    salinityMort=salinityMort,
    swimLightLevel=swimLightLevel,
    vertSwimSpeedMean=vertSwimSpeedMean,
    vertSwimSpeedStd=vertSwimSpeedStd,
    sinkingRateMean=sinkingRateMean,
    sinkingRateStd=sinkingRateStd,
    viabletime=viabletime,
    maxParticleAge=maxParticleAge,
    viableDegreeDays=viableDegreeDays,
    maxDegreeDays=maxDegreeDays,
    recordPsteps=recordPsteps,
    splitPsteps=splitPsteps,
    pstepsInterval=pstepsInterval,
    recordMovement=recordMovement,
    recordElemActivity=recordElemActivity,
    recordConnectivity=recordConnectivity,
    connectivityInterval=connectivityInterval,
    recordLocations=recordLocations,
    recordArrivals=recordArrivals
  )
  return(paste(names(params), params, sep="=", collapse="\n"))
}


