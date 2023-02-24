


rm(list=ls())
gc()



# setup -------------------------------------------------------------------

library(tidyverse); library(glue); library(lubridate); library(sf); library(ncdf4); library(gganimate)
source("code/00_fn.R")
theme_set(theme_classic())

re_extract <- T

dirs <- switch(get_os(),
               windows=list(proj=getwd(),
                            mesh="D:/hydroOut/",
                            hydro.linnhe="D:/hydroOut/linnhe7/linnhe7_tides_met_tsobc_riv",
                            hydro.westcoms="D:/hydroOut/WeStCOMS2/Archive"),
               linux=list(proj=getwd(),
                          mesh="/home/sa04ts/FVCOM_meshes",
                          out=glue("{getwd()}/out/")))

corran <- st_buffer(st_as_sfc("POINT(201850 763450)", crs=27700), 17.5e3)
mesh.fp <- st_read("data/linnhe_mesh_footprint.gpkg")
mesh.sf <- list(linnhe7=st_read(glue("{dirs$mesh}/linnhe_mesh.gpkg")) %>%
                  mutate(mesh="linnhe7"),
                westcoms2=st_read(glue("{dirs$mesh}/WeStCOMS2_linnhe_mesh.gpkg")) %>%
                  mutate(mesh="WeStCOMS2")) %>%
  do.call('rbind', .) %>%
  select(i, mesh, geom) %>%
  filter(st_within(., corran, sparse=F)[,1])
mesh.nodes <- list(linnhe7=st_read("temp/linnhe_1h_nodeLay1_01.gpkg") %>%
                     mutate(mesh="linnhe7"),
                   westcoms2=st_read("temp/WeStCOMS2_1h_nodeLay1_01.gpkg") %>%
                     mutate(mesh="WeStCOMS2")) %>%
  do.call('rbind', .) %>%
  select(i, mesh, geom) %>%
  filter(st_within(., corran, sparse=F)[,1])
wc_i.node <- read_csv("data/westcoms2-linnhe_nodes.csv")$i
wc_i.elem <- read_csv("data/westcoms2-linnhe_elems.csv")$i
time.key <- read_csv("data/timeRes_key.csv") %>% 
  rename(fileHour=layer) %>%
  mutate(fileDate=as.character(fileDate)) %>%
  select(mesh, timeRes, fileDate, fileHour, timeCalculated)



# extract & save surface conditions ---------------------------------------

if(re_extract) {
  # linnhe 1h
  nc.f <- dir(glue("{dirs$hydro.linnhe}/netcdf_2021"))
  linnhe_1h <- map(nc.f,
                   ~extract_surface(nc=.x, dir=dirs$hydro.linnhe, 
                                    subdir="netcdf_2021", meshRes="linnhe7, 1h"))
  # linnhe 5min
  nc.f <- dir(glue("{dirs$hydro.linnhe}_5min/netcdf_2021"))
  linnhe_5min_1 <- map(nc.f,
                       ~extract_surface(nc=.x, dir=paste0(dirs$hydro.linnhe, "_5min"), 
                                        subdir="netcdf_2021", meshRes="linnhe7, 5min"))
  nc.f <- dir(glue("{dirs$hydro.linnhe}_5min/netcdf_2022"))
  linnhe_5min_2 <- map(nc.f,
                       ~extract_surface(nc=.x, dir=paste0(dirs$hydro.linnhe, "_5min"), 
                                        subdir="netcdf_2022", meshRes="linnhe7, 5min"))
  # WeStCOMS2 1h
  nc.f <- dir(glue("{dirs$hydro.westcoms}/netcdf_2021"))
  westcoms_1h <- map(nc.f,
                     ~extract_surface(nc=.x, dir=dirs$hydro.westcoms, 
                                      subdir="netcdf_2021", meshRes="WeStCOMS2, 1h",
                                      westcoms_i=list(node=wc_i.node, elem=wc_i.elem)))
  
  # WeStCOMS2 5min
  nc.f <- dir(glue("{dirs$hydro.westcoms}_5min/netcdf_2021"))
  westcoms_5min_1 <- map(nc.f,
                         ~extract_surface(nc=.x, dir=paste0(dirs$hydro.westcoms, "_5min"), 
                                          subdir="netcdf_2021", meshRes="WeStCOMS2, 5min",
                                          westcoms_i=list(node=wc_i.node, elem=wc_i.elem)))
  nc.f <- dir(glue("{dirs$hydro.westcoms}_5min/netcdf_2022"))
  westcoms_5min_2 <- map(nc.f,
                         ~extract_surface(nc=.x, dir=paste0(dirs$hydro.westcoms, "_5min"), 
                                          subdir="netcdf_2022", meshRes="WeStCOMS2, 5min",
                                          westcoms_i=list(node=wc_i.node, elem=wc_i.elem)))
  
  elem.df <- bind_rows(map_dfr(linnhe_1h, ~.x$elem),
                       map_dfr(linnhe_5min_1, ~.x$elem),
                       map_dfr(linnhe_5min_2, ~.x$elem),
                       map_dfr(westcoms_1h, ~.x$elem),
                       map_dfr(westcoms_5min_1, ~.x$elem),
                       map_dfr(westcoms_5min_2, ~.x$elem)) %>%
    left_join(time.key %>% 
                mutate(meshRes=paste(mesh, timeRes, sep=", ")) %>%
                select(meshRes, fileDate, fileHour, timeCalculated)) %>%
    select(meshRes, timeCalculated, i, currentSpeed)
  saveRDS(elem.df, "temp/elem_df.rds")
  rm(elem.df); gc()
  node.df <- bind_rows(map_dfr(linnhe_1h, ~.x$node),
                       map_dfr(linnhe_5min_1, ~.x$node),
                       map_dfr(linnhe_5min_2, ~.x$node),
                       map_dfr(westcoms_1h, ~.x$node),
                       map_dfr(westcoms_5min_1, ~.x$node),
                       map_dfr(westcoms_5min_2, ~.x$node)) %>%
    left_join(time.key %>% 
                mutate(meshRes=paste(mesh, timeRes, sep=", ")) %>%
                select(meshRes, fileDate, fileHour, timeCalculated)) %>%
    select(meshRes, timeCalculated, i, temp, salinity)
  saveRDS(node.df, "temp/node_df.rds")
  rm(node.df)
  rm(list=c("linnhe_1h", "linnhe_5min_1", "linnhe_5min_2",
            "westcoms_1h", "westcoms_5min_1", "westcoms_5min_2"))
  gc() 
}



# animations --------------------------------------------------------------

hours_per_s <- 4
node.sf <- left_join(mesh.nodes, 
                     readRDS("temp/node_df.rds") %>%
                       filter(grepl("1h", meshRes)) %>%
                       filter(timeCalculated <= "2021-11-02 00:00:00") %>%
                       mutate(mesh=str_split_fixed(meshRes, ",", 2)[,1]))
anim <- node.sf %>%
  filter(!is.na(timeCalculated)) %>%
  ggplot(aes(colour=salinity)) + 
  geom_sf(size=0.5) + 
  scale_colour_distiller(palette="Blues", breaks=seq(10,30,by=5)) +
  transition_time(timeCalculated) +
  facet_grid(.~meshRes) +
  ggtitle( "{frame_time}")
anim_save(glue("figs/mesh_salinity_1h.gif"),
          anim, nframes=25,
          fps=hours_per_s, width=9, height=4, res=300, units="in")
anim <- node.sf %>%
  filter(!is.na(timeCalculated)) %>%
  ggplot(aes(colour=temp)) + 
  geom_sf(size=0.5) + 
  scale_colour_distiller(type="div", palette="RdBu") +
  transition_time(timeCalculated) +
  facet_grid(.~meshRes) +
  ggtitle( "{frame_time}")
anim_save(glue("figs/mesh_temperature_1h.gif"),
          anim, nframes=25,
          fps=hours_per_s, width=9, height=4, res=300, units="in")
rm(anim); rm(node.sf)

elem.sf <- left_join(mesh.sf, 
                     readRDS("temp/elem_df.rds") %>%
                       filter(grepl("1h", meshRes)) %>%
                       filter(timeCalculated <= "2021-11-02 00:00:00") %>%
                       mutate(mesh=str_split_fixed(meshRes, ",", 2)[,1]))
anim <- elem.sf %>%
  ggplot(aes(fill=currentSpeed)) + 
  geom_sf(colour=NA) + 
  scale_fill_viridis_c() +
  transition_time(timeCalculated) +
  facet_grid(.~meshRes) +
  ggtitle( "{frame_time}")
anim_save(glue("figs/mesh_currentSpeed_linnhe7_1h.gif"),
          anim, nframes=25,
          fps=hours_per_s, width=9, height=4, res=300, units="in")
rm(anim); rm(elem.sf)

node.sf <- left_join(mesh.nodes, 
                     readRDS("temp/node_df.rds") %>%
                       filter(grepl("5min", meshRes)) %>%
                       filter(timeCalculated <= "2021-11-02 00:00:00") %>%
                       mutate(mesh=str_split_fixed(meshRes, ",", 2)[,1]))
anim <- node.sf %>%
  filter(!is.na(timeCalculated)) %>%
  ggplot(aes(colour=salinity)) + 
  geom_sf(size=0.5) + 
  scale_colour_distiller(palette="Blues") +
  transition_time(timeCalculated) +
  facet_grid(.~meshRes) +
  ggtitle( "{frame_time}")
anim_save(glue("figs/mesh_salinity_5min.gif"),
          anim, nframes=289,
          fps=hours_per_s*12, width=9, height=4, res=300, units="in")
anim <- node.sf %>%
  filter(!is.na(timeCalculated)) %>%
  ggplot(aes(colour=temp)) + 
  geom_sf(size=0.5) + 
  scale_colour_distiller(type="div", palette="RdBu") +
  transition_time(timeCalculated) +
  facet_grid(.~meshRes) +
  ggtitle( "{frame_time}")
anim_save(glue("figs/mesh_temperature_5min.gif"),
          anim, nframes=289,
          fps=hours_per_s*12, width=9, height=4, res=300, units="in")
rm(anim); rm(node.sf)

elem.sf <- left_join(mesh.sf, 
                     readRDS("temp/elem_df.rds") %>%
                       filter(grepl("5min", meshRes)) %>%
                       filter(timeCalculated <= "2021-11-02 00:00:00") %>%
                       mutate(mesh=str_split_fixed(meshRes, ",", 2)[,1]))
anim <- elem.sf %>%
  filter(timeCalculated < "2021-11-02 00:00:00") %>%
  ggplot(aes(fill=currentSpeed)) + 
  geom_sf(colour=NA) + 
  scale_fill_viridis_c() +
  transition_time(timeCalculated) +
  facet_grid(.~meshRes) +
  ggtitle( "{frame_time}")
anim_save(glue("figs/mesh_currentSpeed_linnhe7_5min.gif"),
          anim, nframes=289,
          fps=24, width=9, height=4, res=300, units="in")






# salinity profiles -------------------------------------------------------

glue("out/salinity_profile_{par.df$mesh[j]}_{par.df$tRes[j]}.csv")

prof.df <- dir("out", "salinity_profile") %>% 
  map_dfr(~read_csv(paste0("out/", .x)) %>% 
            mutate(meshRes=str_sub(str_remove(.x, "salinity_profile_"),1,-5))) %>%
  mutate(meshRes=str_replace("_", ", "))

prof.df %>%
  ggplot(aes(fill=salinity)) +
  geom_rect(aes(xmin=timeStart, xmax=timeEnd, ymin=depthBottom, ymax=depthTop)) + 
  scale_fill_distiller(palette="Blues", direction=1, limits=c(17, 32)) +
  scale_y_continuous("Depth (m)") +
  scale_x_datetime("", date_breaks="1 day", date_labels="%d-%b") +
  facet_wrap(~meshRes) +
  theme(axis.text.x=element_text(angle=300, hjust=0, vjust=0.5),
        axis.title.x=element_blank())
ggsave(glue("figs/mesh/salinity_profiles.png"), 
       width=7, height=6, dpi=300)

prof.df %>%
  filter(timeEnd <= "2021-11-03 00:00:00") %>%
  ggplot(aes(fill=salinity)) +
  geom_rect(aes(xmin=timeStart, xmax=timeEnd, ymin=depthBottom, ymax=depthTop)) + 
  scale_fill_distiller(palette="Blues", direction=1, limits=c(17, 32)) +
  scale_y_continuous("Depth (m)", limits=c(-30, 0)) +
  scale_x_datetime("", date_breaks="1 day", date_labels="%d-%b") +
  facet_wrap(~meshRes) +
  theme(axis.text.x=element_text(angle=300, hjust=0, vjust=0.5),
        axis.title.x=element_blank())
ggsave(glue("figs/mesh/salinity_profiles_2days.png"), 
       width=7, height=6, dpi=300)
