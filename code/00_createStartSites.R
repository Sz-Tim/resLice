# Start sites
# Off Aqua
# Tim Szewczyk



# setup -------------------------------------------------------------------

library(tidyverse); library(ncdf4); library(sf)

hiRes.dir <- "W:/common/sa04ts-temp/linnhe7/"
hiRes_1h.dir <- paste0(hiRes.dir, "linnhe7_tides_met_tsobc_mf/")

mesh.nc <- nc_open("data/linnhe_mesh.nc")
mesh.sf <- as_tibble(ncvar_get(mesh.nc, "nodexy_os")) %>%
  st_as_sf(coords=c("V1", "V2"), crs=27700)

mesh.pts <- st_union(mesh.sf)

gridRes <- 500
lice_grid <- st_bbox(mesh.pts) %>%
  st_make_grid(cellsize=c(gridRes, gridRes), what="centers") %>%
  st_sf(id=1:length(.)) %>%
  mutate(inMesh=st_is_within_distance(., mesh.pts, 1000, sparse=F)[,1]) %>%
  filter(inMesh)

st_coordinates(lice_grid) %>% as_tibble() %>%
  mutate(i=row_number(), depth=0) %>%
  select(i, X, Y, depth) %>%
  write_tsv(paste0("data/linnhe_start_", gridRes, "m_grid.tsv"),
            col_names=F)
