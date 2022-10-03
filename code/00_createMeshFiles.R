# Mesh creation
# Off Aqua
# Tim Szewczyk

# This script creates a mesh file for the hi resolution Loch Linnhe domain. This
# is probably not the cleanest or most efficient way to do this, but I need a
# mesh file that aligns with the expected names and variables loaded by the
# biotracker code, even if all of the information could in theory be read from a
# regular model output file instead.



# setup -------------------------------------------------------------------

library(tidyverse); library(ncdf4); library(sf)

hiRes.dir <- "W:/common/sa04ts-temp/linnhe7/"
hiRes_1h.dir <- paste0(hiRes.dir, "linnhe7_tides_met_tsobc_mf/")

day2 <- nc_open(dir(hiRes_1h.dir, "0002.nc", full.names=T))
westcoms <- st_read("../../WeStCOMS/data/WeStCOMS2_mesh.gpkg")




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
nc$boundaryNodesAll <- which(rowSums(nc$nbe==0) > 0)
nc$boundaryNodesOpen <- nc$boundaryNodesAll # not used, but must be loaded...
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
obcnode <- ncdim_def("obcnode", "boundary", 1:length(nc$boundaryNodesAll))



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

mesh.nc <- nc_create("data/linnhe_mesh.nc", mesh.vars)
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
  mutate(coords=list(nc$nodexy_os[c(trinode_1, trinode_2, trinode_3, trinode_1),])) %>%
  ungroup %>%
  mutate(geometry=list(st_polygon(coords))) %>%
  st_as_sf(crs=27700)
mesh.sf_ <- mesh.sf %>%
  mutate(area=st_area()) %>%
  select(i, area, depth, trinode_1, trinode_2, trinode_3, geom)
write_sf(mesh.sf_, "data/linnhe_mesh.gpkg")



# confirm -----------------------------------------------------------------

mesh.nc <- nc_open("data/linnhe_mesh.nc")
mesh.nc
str(ncvar_get(mesh.nc, "uvnode"))
str(ncvar_get(mesh.nc, "uvnode_os"))
str(ncvar_get(mesh.nc, "nodexy_os"))
head(ncvar_get(mesh.nc, "trinodes"))
head(ncvar_get(mesh.nc, "nbe"))
ncvar_get(mesh.nc, "siglay")

nc_close(mesh.nc)
