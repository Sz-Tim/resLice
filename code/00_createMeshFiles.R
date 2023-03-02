# Mesh creation
# Off Aqua
# Tim Szewczyk

# This script creates a mesh file for the hi resolution Loch Linnhe domain. This
# is probably not the cleanest or most efficient way to do this, but I need a
# mesh file that aligns with the expected names and variables loaded by the
# biotracker code, even if all of the information could in theory be read from a
# regular model output file instead.



# setup -------------------------------------------------------------------

library(tidyverse); library(ncdf4); library(sf); library(glue)

domain <- c("linnhe", "WeStCOMS2_linnhe", "WeStCOMS2_noLinnhe")[3]

westcoms.dir <- "D:/hydroOut/WeStCOMS2/Archive/netcdf_2021/"
hiRes.dir <- "D:hydroOut/linnhe7/"
hiRes_1h.dir <- paste0(hiRes.dir, "linnhe7_tides_met_tsobc_riv/")
hiRes_5min.dir <- paste0(hiRes.dir, "linnhe7_tides_met_tsobc_riv_5min/")

if(domain=="linnhe") {
  day2 <- nc_open(dir(hiRes_1h.dir, "0002.nc", full.names=T))
  open_elems <- read_csv("data/linnhe_mesh_openBoundaryElems.csv")
} else if(domain=="WeStCOMS2_linnhe") {
  day2 <- nc_open(dir(westcoms.dir, "20211102", full.names=T))
  open_elems <- read_csv("data/westcoms2-linnhe_mesh_openBoundaryElems.csv")
} else if(domain=="WeStCOMS2_noLinnhe") {
  day2 <- nc_open(dir(westcoms.dir, "20211102", full.names=T))
  open_elems <- read_csv("data/westcoms2-linnhe_elems.csv")
}



# extract dims ------------------------------------------------------------

# Variables to extract, following Tom's naming convention for clarity (see
# WeStCOMS2_mesh.nc) and to avoid rewriting a bunch of Java code
nc <- list(
  uvnode=cbind(ncvar_get(day2, "lonc"), ncvar_get(day2, "latc")),
  nodexy=cbind(ncvar_get(day2, "lon"), ncvar_get(day2, "lat")),
  trinodes=ncvar_get(day2, "nv"),
  nbe=ncvar_get(day2, "nbe"),
  depthNodexy=ncvar_get(day2, "h"),
  depthUvnode=ncvar_get(day2, "h_center"),
  siglay=ncvar_get(day2, "siglay")[1,], 
  siglev=ncvar_get(day2, "siglev")[1,]
)
nc_close(day2)
nc$boundaryNodesAll <- which(rowSums(nc$nbe==0) > 0)
nc$boundaryNodesOpen <- sort(unique(unlist(open_elems %>% select(starts_with("trinode")))))
nc$uvnode_os <- as_tibble(nc$uvnode) %>%
  st_as_sf(coords=c("V1", "V2"), crs=4326) %>%
  st_transform(27700) %>% st_coordinates
nc$nodexy_os <- as_tibble(nc$nodexy) %>%
  st_as_sf(coords=c("V1", "V2"), crs=4326) %>%
  st_transform(27700) %>% st_coordinates


# assign dims -------------------------------------------------------------

elems <- ncdim_def("elems", "element", 1:nrow(nc$uvnode))
two <- ncdim_def("two", "two", 1:2)
node <- ncdim_def("node", "nodes", 1:nrow(nc$nodexy))
nele <- elems
three <- ncdim_def("three", "three", 1:3)
siglayers <- ncdim_def("siglayers", "sigma layers", 1:length(nc$siglay))
siglevels <- ncdim_def("siglevels", "sigma levels", 1:length(nc$siglev))
bnode <- ncdim_def("bnode", "boundary", 1:length(nc$boundaryNodesAll))
obcnode <- ncdim_def("obcnode", "boundary", 1:length(nc$boundaryNodesOpen))



# create vars -------------------------------------------------------------

mesh.vars <- list(
  uvnode=ncvar_def("uvnode", "degrees", list(elems, two)),
  nodexy=ncvar_def("nodexy", "degrees", list(node, two)),
  uvnode_os=ncvar_def("uvnode_os", "m", list(elems, two)),
  nodexy_os=ncvar_def("nodexy_os", "m", list(node, two)),
  trinodes=ncvar_def("trinodes", "", list(nele, three), prec="integer"),
  nbe=ncvar_def("nbe", "", list(nele, three), prec="integer"),
  depthUvnode=ncvar_def("depthUvnode", "m", elems),
  depthNodexy=ncvar_def("depthNodexy", "m", node),
  boundaryNodesAll=ncvar_def("boundaryNodesAll", "", bnode, prec="integer"),
  boundaryNodesOpen=ncvar_def("boundaryNodesOpen", "", obcnode, prec="integer"),
  siglay=ncvar_def("siglay", "proportion", siglayers),
  siglev=ncvar_def("siglev", "proportion", siglevels)
)



# create mesh.nc ----------------------------------------------------------

mesh.nc <- nc_create(glue("data/{domain}_mesh.nc"), mesh.vars)
iwalk(names(mesh.vars), ~ncvar_put(mesh.nc, .x, nc[[.x]]))
nc_close(mesh.nc)



# create mesh.gpkg --------------------------------------------------------

mesh.sf <- tibble(i=1:nrow(nc$uvnode_os),
                  depth=nc$depthUvnode,
                  lonc=nc$uvnode_os[,1], 
                  latc=nc$uvnode_os[,2],
                  trinode_1=nc$trinodes[,1],
                  trinode_2=nc$trinodes[,2],
                  trinode_3=nc$trinodes[,3]) %>%
  rowwise() %>%
  mutate(coords=list(nc$nodexy_os[c(trinode_1, trinode_2, trinode_3, trinode_1),]),
         geom=list(st_polygon(list(coords)))) %>%
  ungroup %>%
  st_as_sf(crs=27700) %>%
  mutate(area=st_area(.)) %>%
  select(i, area, depth, trinode_1, trinode_2, trinode_3, geom)
write_sf(mesh.sf, "data/linnhe_mesh.gpkg")

mesh.footprint <- mesh.sf %>% st_union
write_sf(mesh.footprint, "data/linnhe_mesh_footprint.gpkg")



# confirm -----------------------------------------------------------------

mesh.nc <- nc_open(glue("data/{domain}_mesh.nc"))
mesh.nc
str(ncvar_get(mesh.nc, "uvnode"))
str(ncvar_get(mesh.nc, "uvnode_os"))
str(ncvar_get(mesh.nc, "nodexy_os"))
head(ncvar_get(mesh.nc, "trinodes"))
head(ncvar_get(mesh.nc, "nbe"))
head(ncvar_get(mesh.nc, "boundaryNodesOpen"))
ncvar_get(mesh.nc, "siglay")

nc_close(mesh.nc)





# loch bearings -----------------------------------------------------------

loch_pts <- st_read("data/temp/linnhe_skeleton.shp") %>% 
  arrange(lochRegion, lochHead)
loch_dirs <- loch_pts %>%
  mutate(x=st_coordinates(.)[,1],
         y=st_coordinates(.)[,2]) %>%
  st_drop_geometry() %>%
  group_by(lochRegion) %>%
  summarise(dx=last(x)-first(x),
            dy=last(y)-first(y)) %>%
  mutate(bearing=atan2(dy, dx))
loch_dirs %>% select(lochRegion, bearing) %>%
  write_csv("data/loch_regions.csv")

# loch_lines <- loch_pts %>%
#   group_by(lochRegion) %>%
#   summarise() %>%
#   st_cast("LINESTRING") %>% 
#   mutate(bearing=loch_dirs$bearing)
# ggplot(loch_lines, aes(colour=bearing)) + geom_sf(size=2) +
#   scale_colour_gradientn(colours=cmr$infinity, limits=c(-pi, pi),
#                          breaks=c(-pi, -pi/2, 0, pi/2, pi),
#                          labels=c("W", "S", "E", "N", "W"))

