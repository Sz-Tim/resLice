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



setPartTrackProperties <- function(
  destinationDirectory="out/",
  datadir="/media/archiver/common/sa01da-work/WeStCOMS2/Archive/",
  datadirPrefix="netcdf_",
  datadirSuffix="",
  mesh1="/home/sa04ts/FVCOM_meshes/WeStCOMS2_linnhe_mesh.nc",
  mesh1Type="FVCOM",
  coordRef="OSGB1936",
  location="westcoms",
  minchVersion=2,
  sitefile="../../data/linnhe_start_1000m_grid.tsv",
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
  mortalityRate=0,
  salinityMort="true",
  swimLightLevel="true",
  vertSwimSpeedMean=-0.0005,
  vertSwimSpeedStd=0.0001,
  sinkingRateMean=0.0005,
  sinkingRateStd=0.0001,
  viabletime=12,
  maxParticleAge=500,
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
    datadir=datadir,
    datadirPrefix=datadirPrefix,
    datadirSuffix=datadirSuffix,
    mesh1=mesh1,
    mesh1Type=mesh1Type,
    coordRef=coordRef,
    location=location,
    minchVersion=minchVersion,
    sitefile=sitefile,
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


