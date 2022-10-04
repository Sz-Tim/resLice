# Start sites
# Off Aqua
# Tim Szewczyk



# setup -------------------------------------------------------------------

library(tidyverse); library(sf); library(glue)
theme_set(theme_classic())

mesh.fp <- st_read("data/linnhe_mesh_footprint.gpkg")


gridRes <- 1000
lice_grid <- st_bbox(mesh.fp) %>%
  st_make_grid(cellsize=c(gridRes, gridRes), what="centers") %>%
  st_sf(id=1:length(.)) %>%
  mutate(inMesh=st_within(., mesh.fp, sparse=F)[,1]) %>%
  filter(inMesh)

st_coordinates(lice_grid) %>% as_tibble() %>%
  mutate(i=row_number(), depth=0) %>%
  select(i, X, Y, depth) %>%
  write_tsv(paste0("data/linnhe_start_", gridRes, "m_grid.tsv"),
            col_names=F)

ggplot(lice_grid) + 
  geom_sf(data=mesh.fp, fill="cadetblue") + 
  geom_sf(size=0.5) +
  ggtitle(glue("Start sites: {gridRes}m grid ({nrow(lice_grid)} pts)"))
ggsave(glue("figs/startSites_{gridRes}m_grid.png"), width=5, height=6)
