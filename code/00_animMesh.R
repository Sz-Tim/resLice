


rm(list=ls())
gc()



# setup -------------------------------------------------------------------

library(tidyverse); library(glue); library(lubridate); library(sf); library(ncdf4); library(magick)
source("code/00_fn.R")
theme_set(theme_classic())

re_extract <- F
focus <- c("all", "corran", "etive")[3]

dirs <- switch(get_os(),
               windows=list(proj=getwd(),
                            mesh="D:/hydroOut/",
                            hydro.linnhe="D:/hydroOut/linnhe7/linnhe7_tides_met_tsobc_riv",
                            hydro.westcoms="D:/hydroOut/WeStCOMS2/Archive"),
               linux=list(proj=getwd(),
                          mesh="/home/sa04ts/FVCOM_meshes",
                          out=glue("{getwd()}/out/")))

corran <- st_buffer(st_as_sfc("POINT(201850 763450)", crs=27700), 17.5e3)
etive.bbox <- list(xmin=185276, xmax=212845, ymin=728357, ymax=746971)
mesh.fp <- st_read("data/linnhe_mesh_footprint.gpkg")
mesh.sf <- list(linnhe7=st_read(glue("{dirs$mesh}/linnhe_mesh.gpkg")) %>%
                  mutate(mesh="linnhe7"),
                westcoms2=st_read(glue("{dirs$mesh}/WeStCOMS2_linnhe_mesh.gpkg")) %>%
                  mutate(mesh="WeStCOMS2")) %>%
  do.call('rbind', .) %>%
  select(-lochRegion, -mainLoch, -depth, -area)
mesh.nodes <- list(linnhe7=st_read("temp/linnhe_1h_nodeLay1_01.gpkg") %>%
                     mutate(mesh="linnhe7"),
                   westcoms2=st_read("temp/WeStCOMS2_1h_nodeLay1_01.gpkg") %>%
                     mutate(mesh="WeStCOMS2")) %>%
  do.call('rbind', .) %>%
  select(i, mesh, geom)
if(focus=="corran") {
  mesh.sf <- mesh.sf %>%
    filter(st_within(., corran, sparse=F)[,1])
  mesh.nodes <- mesh.nodes %>%
    filter(st_within(., corran, sparse=F)[,1])
} else if(focus=="etive") {
  mesh.sf <- mesh.sf %>%
    st_crop(st_bbox(unlist(etive.bbox), crs=27700))
  mesh.nodes <- mesh.nodes %>%
    st_crop(st_bbox(unlist(etive.bbox), crs=27700))
}
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
    select(meshRes, timeCalculated, i, currentSpeed, vertSpeed)
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

hours_per_s <- 8



# 1h mesh animations ------------------------------------------------------

node.df <- readRDS("temp/node_df.rds") %>%
  filter(grepl("1h", meshRes)) %>%
  filter(!is.na(timeCalculated)) %>%
  rename(trinode=i)

elem.sf <- left_join(mesh.sf,
                     readRDS("temp/elem_df.rds") %>%
                       filter(grepl("1h", meshRes)) %>%
                       mutate(mesh=str_split_fixed(meshRes, ",", 2)[,1])) %>%
  filter(!is.na(timeCalculated)) %>%
  left_join(node.df, by=c("meshRes", "timeCalculated", "trinode_1"="trinode"), suffix=c("", "_1")) %>%
  left_join(node.df, by=c("meshRes", "timeCalculated", "trinode_2"="trinode"), suffix=c("", "_2")) %>%
  left_join(node.df, by=c("meshRes", "timeCalculated", "trinode_3"="trinode"), suffix=c("", "_3")) %>%
  rowwise() %>%
  mutate(temp_i=sum(temp, temp_2, temp_3, na.rm=T)/
           sum(!is.na(temp), !is.na(temp_2), !is.na(temp_3)),
         salinity_i=sum(salinity, salinity_2, salinity_3, na.rm=T)/
           sum(!is.na(salinity), !is.na(salinity_2), !is.na(salinity_3))) %>%
  select(-temp, -temp_2, -temp_3, -salinity, -salinity_2, -salinity_3, -starts_with("trinode"))

rm(node.df); gc()

elem.times <- sort(unique(elem.sf$timeCalculated))
sal.rng <- range(elem.sf$salinity_i)
temp.rng <- range(elem.sf$temp_i)
currentSpeed.rng <- range(elem.sf$currentSpeed)
vertSpeed.rng <- range(elem.sf$vertSpeed)
# for(i in 1:length(elem.times)) {
#   i.int <- as.integer(i)
#   date_i <- date(elem.times[i])
#   hour <- str_pad(hour(elem.times[i]), 2, "l", "0")
# 
#   p.i <- elem.sf %>%
#     filter(timeCalculated==elem.times[i.int]) %>%
#     ggplot(aes(fill=salinity_i)) +
#     geom_sf(colour=NA) +
#     scale_fill_distiller("Salinity (psu)", palette="Blues", breaks=seq(10,30,by=5), direction=1,
#                          limits=sal.rng) +
#     facet_grid(.~meshRes) +
#     ggtitle(paste0(date_i, " ", hour, ":00"))
#   ggsave(glue("figs/mesh/temp/mesh_salinity_1h_{focus}_{str_pad(i,4,'l','0')}.png"),
#          p.i, width=9, height=4, dpi=300)
# }
# for(i in 1:length(elem.times)) {
#   i.int <- as.integer(i)
#   date_i <- date(elem.times[i])
#   hour <- str_pad(hour(elem.times[i]), 2, "l", "0")
# 
#   p.i <- elem.sf %>%
#     filter(timeCalculated==elem.times[i.int]) %>%
#     ggplot(aes(fill=temp_i)) +
#     geom_sf(colour=NA) +
#     scale_fill_distiller("Temperature (C)", type="div", palette="RdBu", limits=temp.rng) +
#     facet_grid(.~meshRes) +
#     ggtitle(paste0(date_i, " ", hour, ":00"))
#   ggsave(glue("figs/mesh/temp/mesh_temp_1h_{focus}_{str_pad(i,4,'l','0')}.png"),
#          p.i, width=9, height=4, dpi=300)
# }
# for(i in 1:length(elem.times)) {
#   i.int <- as.integer(i)
#   date_i <- date(elem.times[i])
#   hour <- str_pad(hour(elem.times[i]), 2, "l", "0")
# 
#   p.i <- elem.sf %>%
#     filter(timeCalculated==elem.times[i.int]) %>%
#     ggplot(aes(fill=currentSpeed)) +
#     geom_sf(colour=NA) +
#     scale_fill_viridis_c("uv speed (m/s)", limits=currentSpeed.rng) +
#     facet_grid(.~meshRes) +
#     ggtitle(paste0(date_i, " ", hour, ":00"))
#   ggsave(glue("figs/mesh/temp/mesh_currentSpeed_1h_{focus}_{str_pad(i,4,'l','0')}.png"),
#          p.i, width=9, height=4, dpi=300)
# }
for(i in 111:length(elem.times)) {
  i.int <- as.integer(i)
  date_i <- date(elem.times[i])
  hour <- str_pad(hour(elem.times[i]), 2, "l", "0")

  p.i <- elem.sf %>%
    filter(timeCalculated==elem.times[i.int]) %>%
    ggplot(aes(fill=vertSpeed)) +
    geom_sf(colour=NA) +
    scale_fill_gradient2("ww speed (m/s)", limits=vertSpeed.rng, mid="grey95") +
    facet_grid(.~meshRes) +
    ggtitle(paste0(date_i, " ", hour, ":00"))
  ggsave(glue("figs/mesh/temp/mesh_vertSpeed_1h_{focus}_{str_pad(i,4,'l','0')}.png"),
         p.i, width=9, height=4, dpi=300)
}

# cat("Making salinity 1h\n")
# dir("figs/mesh/temp", glue("mesh_salinity_1h_{focus}_"), full.names=T) %>%
#   image_read() %>%
#   image_join() %>%
#   image_animate(delay=1/hours_per_s, optimize=T) %>%
#   image_write(glue("figs/mesh/mesh_salinity_1h_{focus}.gif"))
# 
# cat("Making temp 1h\n")
# dir("figs/mesh/temp", glue("mesh_temp_1h_{focus}_"), full.names=T) %>%
#   image_read() %>%
#   image_join() %>%
#   image_animate(delay=1/hours_per_s, optimize=) %>%
#   image_write(glue("figs/mesh/mesh_temperature_1h_{focus}.gif"))

cat("Making currentSpeed 1h\n")
dir("figs/mesh/temp", glue("mesh_currentSpeed_1h_{focus}_"), full.names=T) %>% 
  image_read() %>% 
  image_join() %>% 
  image_animate(delay=1/hours_per_s) %>% 
  image_write(glue("figs/mesh/mesh_currentSpeed_1h_{focus}.gif"))

cat("Making vertSpeed 1h\n")
dir("figs/mesh/temp", glue("mesh_vertSpeed_1h_{focus}_"), full.names=T) %>% 
  image_read() %>% 
  image_join() %>% 
  image_animate(delay=1/hours_per_s) %>% 
  image_write(glue("figs/mesh/mesh_vertSpeed_1h_{focus}.gif"))

rm(elem.sf)




# 5 min mesh animations ---------------------------------------------------

node.df <- readRDS("temp/node_df.rds") %>%
  filter(grepl("5min", meshRes)) %>%
  filter(!is.na(timeCalculated)) %>% 
  rename(trinode=i)

elem.sf <- left_join(mesh.sf,
                     readRDS("temp/elem_df.rds") %>%
                       filter(grepl("5min", meshRes)) %>%
                       mutate(mesh=str_split_fixed(meshRes, ",", 2)[,1])) %>%
  filter(!is.na(timeCalculated)) %>% 
  left_join(node.df, by=c("meshRes", "timeCalculated", "trinode_1"="trinode"), suffix=c("", "_1")) %>%
  left_join(node.df, by=c("meshRes", "timeCalculated", "trinode_2"="trinode"), suffix=c("", "_2")) %>%
  left_join(node.df, by=c("meshRes", "timeCalculated", "trinode_3"="trinode"), suffix=c("", "_3")) %>%
  rowwise() %>%
  mutate(temp_i=sum(temp, temp_2, temp_3, na.rm=T)/
           sum(!is.na(temp), !is.na(temp_2), !is.na(temp_3)),
         salinity_i=sum(salinity, salinity_2, salinity_3, na.rm=T)/
           sum(!is.na(salinity), !is.na(salinity_2), !is.na(salinity_3))) %>%
  select(-temp, -temp_2, -temp_3, -salinity, -salinity_2, -salinity_3, -starts_with("trinode"))

rm(node.df); gc()

elem.times <- sort(unique(elem.sf$timeCalculated))
sal.rng <- range(elem.sf$salinity_i)
temp.rng <- range(elem.sf$temp_i)
currentSpeed.rng <- range(elem.sf$currentSpeed)
vertSpeed.rng <- range(elem.sf$vertSpeed)
for(i in 1:length(elem.times)) {
  i.int <- as.integer(i)
  date_i <- date(elem.times[i])
  hour <- str_pad(hour(elem.times[i]), 2, "l", "0")
  minutes <- str_pad(minute(elem.time[i]), 2, "l", "0")
  
  p.i <- elem.sf %>%
    filter(timeCalculated==elem.times[i.int]) %>%
    ggplot(aes(fill=salinity_i)) +
    geom_sf(colour=NA) +
    scale_fill_distiller("Salinity (psu)", palette="Blues", breaks=seq(10,30,by=5), direction=1,
                         limits=sal.rng) +
    facet_grid(.~meshRes) +
    ggtitle(paste0(date_i, " ", hour, ":00:", minutes))
  ggsave(glue("figs/mesh/temp/mesh_salinity_5min_{focus}_{str_pad(i,4,'l','0')}.png"), 
         p.i, width=9, height=4, dpi=300)
}
for(i in 1:length(elem.times)) {
  i.int <- as.integer(i)
  date_i <- date(elem.times[i])
  hour <- str_pad(hour(elem.times[i]), 2, "l", "0")
  minutes <- str_pad(minute(elem.time[i]), 2, "l", "0")
  
  p.i <- elem.sf %>%
    filter(timeCalculated==elem.times[i.int]) %>%
    ggplot(aes(fill=temp_i)) +
    geom_sf(colour=NA) +
    scale_fill_distiller("Temperature (C)", type="div", palette="RdBu", limits=temp.rng) +
    facet_grid(.~meshRes) +
    ggtitle(paste0(date_i, " ", hour, ":00:", minutes))
  ggsave(glue("figs/mesh/temp/mesh_temp_5min_{focus}_{str_pad(i,4,'l','0')}.png"), 
         p.i, width=9, height=4, dpi=300)
}

for(i in 1:length(elem.times)) {
  i.int <- as.integer(i)
  date_i <- date(elem.times[i])
  hour <- str_pad(hour(elem.times[i]), 2, "l", "0")
  minutes <- str_pad(minute(elem.time[i]), 2, "l", "0")
  
  p.i <- elem.sf %>%
    filter(timeCalculated==elem.times[i.int]) %>%
    ggplot(aes(fill=currentSpeed)) +
    geom_sf(colour=NA) +
    scale_fill_viridis_c("uv speed (m/s)", limits=currentSpeed.rng) +
    facet_grid(.~meshRes) +
    ggtitle(paste0(date_i, " ", hour, ":00:", minutes))
  ggsave(glue("figs/mesh/temp/mesh_currentSpeed_5min_{focus}_{str_pad(i,4,'l','0')}.png"), 
         p.i, width=9, height=4, dpi=300)
}

for(i in 1:length(elem.times)) {
  i.int <- as.integer(i)
  date_i <- date(elem.times[i])
  hour <- str_pad(hour(elem.times[i]), 2, "l", "0")
  minutes <- str_pad(minute(elem.time[i]), 2, "l", "0")
  
  p.i <- elem.sf %>%
    filter(timeCalculated==elem.times[i.int]) %>%
    ggplot(aes(fill=vertSpeed)) +
    geom_sf(colour=NA) +
    scale_fill_gradient2("ww speed (m/s)", limits=vertSpeed.rng, mid="grey95") +
    facet_grid(.~meshRes) +
    ggtitle(paste0(date_i, " ", hour, ":00:", minutes))
  ggsave(glue("figs/mesh/temp/mesh_vertSpeed_5min_{focus}_{str_pad(i,4,'l','0')}.png"), 
         p.i, width=9, height=4, dpi=300)
}

cat("Making salinity 5min\n")
dir("figs/mesh/temp", glue("mesh_salinity_5min_{focus}_"), full.names=T) %>%
  image_read() %>%
  image_join() %>%
  image_animate(delay=1/(12*hours_per_s)) %>%
  image_write(glue("figs/mesh/mesh_salinity_5min_{focus}.gif"))

cat("Making temp 5min\n")
dir("figs/mesh/temp", glue("mesh_temp_5min_{focus}_"), full.names=T) %>%
  image_read() %>%
  image_join() %>%
  image_animate(delay=1/(12*hours_per_s)) %>%
  image_write(glue("figs/mesh/mesh_temp_5min_{focus}.gif"))

cat("Making currentSpeed 5min\n")
dir("figs/mesh/temp", glue("mesh_currentSpeed_5min_{focus}_"), full.names=T) %>% 
  image_read() %>% 
  image_join() %>% 
  image_animate(delay=1/(12*hours_per_s)) %>% 
  image_write(glue("figs/mesh/mesh_currentSpeed_5min_{focus}.gif"))

cat("Making vertSpeed 5min\n")
dir("figs/mesh/temp", glue("mesh_vertSpeed_5min_{focus}_"), full.names=T) %>% 
  image_read() %>% 
  image_join() %>% 
  image_animate(delay=1/(12*hours_per_s)) %>% 
  image_write(glue("figs/mesh/mesh_vertSpeed_5min_{focus}.gif"))

rm(elem.sf)


# 
# 
# 
# 
# # salinity profiles -------------------------------------------------------
# 
# glue("out/salinity_profile_{par.df$mesh[j]}_{par.df$tRes[j]}.csv")
# 
# site.i <- read_tsv("data/fishFarmSites.tsv", col_names=c("site", "x", "y"))
# for(i in 1:nrow(site.i)) {
#   prof.df <- dir("out", glue("salinity_profile_{site.i$site[i]}")) %>% 
#     map_dfr(~read_csv(paste0("out/", .x)) %>% 
#               mutate(meshRes=str_sub(.x, 25, -5))) %>%
#     mutate(meshRes=str_replace(meshRes, "_", ", "))
#   
#   prof.df %>%
#     ggplot(aes(fill=salinity)) +
#     geom_rect(aes(xmin=timeStart, xmax=timeEnd, ymin=depthBottom, ymax=depthTop)) + 
#     scale_fill_distiller(palette="Blues", direction=1, limits=c(10, 32)) +
#     scale_y_continuous("Depth (m)") +
#     scale_x_datetime("", date_breaks="1 day", date_labels="%d-%b") +
#     facet_wrap(~meshRes) +
#     ggtitle(site.i$site[i]) +
#     theme(axis.text.x=element_text(angle=300, hjust=0, vjust=0.5),
#           axis.title.x=element_blank())
#   ggsave(glue("figs/mesh/salinity_profiles_{site.i$site[i]}.png"), 
#          width=7, height=6, dpi=300)
#   
#   prof.df %>%
#     filter(timeStart >= "2021-11-05 00:00:00" & 
#              timeEnd <= "2021-11-07 00:00:00") %>%
#     ggplot(aes(fill=salinity)) +
#     geom_rect(aes(xmin=timeStart, xmax=timeEnd, ymin=depthBottom, ymax=depthTop)) + 
#     scale_fill_distiller(palette="Blues", direction=1, limits=c(10, 33)) +
#     scale_y_continuous("Depth (m)", limits=c(-30, 0)) +
#     scale_x_datetime("", date_breaks="1 day", date_labels="%d-%b") +
#     facet_wrap(~meshRes) +
#     ggtitle(site.i$site[i]) +
#     theme(axis.text.x=element_text(angle=300, hjust=0, vjust=0.5),
#           axis.title.x=element_blank())
#   ggsave(glue("figs/mesh/salinity_profiles_2days_{site.i$site[i]}.png"), 
#          width=7, height=6, dpi=300)
#   
#   prof.df %>%
#     filter(timeStart >= "2021-11-06 00:00:00" & 
#              timeEnd <= "2021-11-07 00:00:00") %>%
#     ggplot(aes(fill=salinity)) +
#     geom_rect(aes(xmin=timeStart, xmax=timeEnd, ymin=depthBottom, ymax=depthTop)) + 
#     scale_fill_distiller(palette="Blues", direction=1, limits=c(10, 33)) +
#     scale_y_continuous("Depth (m)", limits=c(-30, 0)) +
#     scale_x_datetime("", date_breaks="1 day", date_labels="%d-%b") +
#     facet_wrap(~meshRes) +
#     ggtitle(site.i$site[i]) +
#     theme(axis.text.x=element_text(angle=300, hjust=0, vjust=0.5),
#           axis.title.x=element_blank())
#   ggsave(glue("figs/mesh/salinity_profiles_1day_{site.i$site[i]}.png"), 
#          width=7, height=6, dpi=300)
# }
# 
