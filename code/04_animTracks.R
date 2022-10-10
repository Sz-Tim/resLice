


library(tidyverse); library(glue); library(lubridate); library(sf); library(gganimate)
source("code/00_fn.R")
theme_set(theme_classic())

dirs <- switch(get_os(),
               windows=list(proj=getwd(),
                            mesh="D:/hydroOut/",
                            out=glue("{getwd()}/out/")),
               linux=list(proj=getwd(),
                          mesh="/home/sa04ts/FVCOM_meshes",
                          out=glue("{getwd()}/out/")))

mesh.fp <- st_read("data/linnhe_mesh_footprint.gpkg")
sim_i <- read_csv(glue("{dirs$out}/sim_i.csv")) %>%
  mutate(liceSpeedF=factor(liceSpeed, labels=c("Passive", "Slow", "Medium", "Fast")))

loc.df <- read_csv("out/00_processed/locations.csv") %>%
  filter(status==2) %>%
  select(mesh, liceSpeedF, sim, date, hour, dateTime, ID, x, y, depth, density) %>%
  mutate(liceSpeedF=factor(liceSpeedF, levels=levels(sim_i$liceSpeedF)))
days <- unique(loc.df$date)[5:7]
hours <- unique(loc.df$hour)
particles <- unique(filter(loc.df, date %in% days & hour %in% hours)$ID)
part.sample <- sample(particles, min(length(particles), 1e6))

loc.df <- loc.df %>%
  filter(date %in% days,
         hour %in% hours,
         ID %in% part.sample)
interp <- 8
anim <- loc.df %>%
  ggplot() + 
  geom_sf(data=mesh.fp, fill="grey10", colour=NA) +
  geom_point(aes(x, y, colour=depth, alpha=density, group=ID), size=0.25) + 
  scale_colour_viridis_c(direction=-1, option="D", limits=c(0,20), na.value="#000033") +
  transition_states(dateTime, wrap=F, transition_length=1, state_length=0) +
  facet_grid(mesh~liceSpeedF) +
  ggtitle( "{closest_state}")
anim_save(glue("figs/tracks_AllPart.gif"), 
          anim, nframes=interp*length(hours)*length(days),
          fps=24, width=10, height=6, res=300, units="in")


loc.ls <- loc.df %>%
  filter(date %in% days,
         hour %in% hours,
         ID %in% part.sample) %>%
  group_by(sim) %>%
  group_split()
rm(loc.df)
gc()



library(doSNOW)
cl <- makeCluster(min(nrow(sim_i), 10))
registerDoSNOW(cl)
foreach(i=1:nrow(sim_i), 
        .packages=c("tidyverse", "glue", "lubridate", "sf", "gganimate"),
        .export=c("sim_i", "mesh.fp", "loc.ls", "days", "hours"),
        .errorhandling="pass") %dopar% {
          
  interp <- 4  # interpolation frames
  theme_set(theme_classic())
  
  # anim <- loc.ls[[i]] %>%
  #   st_as_sf(coords=c("x", "y"), crs=27700) %>%
  #   ggplot() + 
  #   geom_sf(data=mesh.fp, fill="grey10", colour=NA) +
  #   geom_sf(aes(colour=depth, alpha=density, group=ID), size=0.5) + 
  #   scale_colour_viridis_c(direction=-1, option="D", limits=c(0,20), na.value="#000033") +
  #   transition_states(dateTime, wrap=F, transition_length=1, state_length=0) +
  #   ggtitle(paste0(glue("{sim_i$mesh[i]}, {sim_i$timeRes[i]}, {sim_i$liceSpeedF[i]}"), 
  #                  ". {closest_state}"))
  anim <- loc.ls[[i]] %>%
    ggplot() + 
    geom_sf(data=mesh.fp, fill="grey10", colour=NA) +
    geom_point(aes(x, y, colour=depth, alpha=density, group=ID), size=0.5) + 
    scale_colour_viridis_c(direction=-1, option="D", limits=c(0,20), na.value="#000033") +
    transition_states(dateTime, wrap=F, transition_length=1, state_length=0) +
    ggtitle(paste0(glue("{sim_i$mesh[i]}, {sim_i$timeRes[i]}, {sim_i$liceSpeedF[i]}"), 
                   ". {closest_state}"))
  anim_save(glue("figs/tracks_sim{str_pad(i,2,'left','0')}_10kPart.gif"), 
            anim, nframes=interp*length(hours)*length(days),
            fps=20, width=7, height=6, res=300, units="in")
  paste("Finished", i)
}
stopCluster(cl)
