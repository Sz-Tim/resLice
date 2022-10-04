# Mesh exploration
# Off Aqua
# Tim Szewczyk

# This script is to debug / test output options



# setup -------------------------------------------------------------------

library(tidyverse); library(ncdf4); library(sf); library(glue); library(lubridate)
theme_set(theme_bw())

# directories
out.dir <- "../sealice_runs/00_test/out_debug/linnhe7_1hr/"

# files
# mesh.sf <- st_read("data/linnhe_mesh.gpkg")
mesh.fp <- st_read("data/linnhe_mesh_footprint.gpkg")
# start.i <- read_delim(glue("{out.dir}startSitesUsed.dat"))
# elemActivity <- read_delim(glue("{out.dir}elementActivity.dat"), 
#                            col_names=c("sink", "swim", "float", "total")) %>%
#   mutate(total=sink+swim+float,
#          i=row_number())
mvmt <- read_delim(glue("{out.dir}movementFile.dat"))
mvmt.df <- bind_rows(mvmt %>% mutate(res="1h"),
                     mvmt2 %>% mutate(res="5min"))


# start sites -------------------------------------------------------------

# ggplot(start.i, aes(x, y, colour=salinity)) + 
#   geom_point() + scale_colour_viridis_c()





# element activity --------------------------------------------------------

# elemAct.sf <- mesh.sf %>%
#   left_join(elemActivity, by="i")
# ggplot(elemAct.sf) + 
#   geom_sf(colour=NA, aes(fill=log1p(total))) +
#   scale_fill_viridis_c()
# ggplot(elemAct.sf) + 
#   geom_sf(colour=NA, aes(fill=sink/total)) +
#   scale_fill_viridis_c()
# ggplot(elemAct.sf) + 
#   geom_sf(colour=NA, aes(fill=swim/total)) +
#   scale_fill_viridis_c()
# ggplot(elemAct.sf) + 
#   geom_sf(colour=NA, aes(fill=float/total)) +
#   scale_fill_viridis_c()



# movement ----------------------------------------------------------------

# part.sample <- sample(unique(mvmt$ID), 1)
# part.sample <- sample(filter(mvmt, startDate=="20211101")$ID, 50)
# 
# mvmt %>% 
#   filter(ID %in% part.sample) %>%
#   st_as_sf(coords=c("x", "y"), crs=27700) %>%
#   mutate(uploch=dX>0 & dY>0,
#          upward=dZ<0) %>%
#   ggplot() + 
#   geom_sf(data=mesh.fp, fill="grey70") +
#   geom_sf(size=0.5, aes(colour=density, group=ID)) + 
#   scale_colour_viridis_c()
# 
# mvmt %>% 
#   filter(ID %in% part.sample) %>%
#   st_as_sf(coords=c("x", "y"), crs=27700) %>%
#   mutate(uploch=dX>0 & dY>0,
#          upward=dZ<0) %>%
#   ggplot() + 
#   geom_sf(data=mesh.fp, fill="grey70") +
#   geom_sf(size=0.5, aes(colour=uploch, alpha=uploch)) +
#   scale_alpha_manual(values=c(0.2, 1))
# 
# mvmt %>% 
#   filter(ID %in% part.sample) %>%
#   ggplot(aes(x, y, colour=-z, group=ID)) + 
#   geom_path() + scale_colour_viridis_c()
# 
# mvmt %>% 
#   filter(ID %in% part.sample) %>%
#   ggplot(aes(x, y, colour=temp-tempSurface, group=ID)) + 
#   geom_path() + scale_colour_gradient2(low=scales::muted("blue"), high=scales::muted("red"))
# 
# mvmt %>% 
#   filter(ID %in% part.sample) %>%
#   ggplot(aes(x, y, colour=dX>0 & dY>0, group=ID)) + 
#   geom_path() +
#   scale_colour_manual("Uploch", values=c("grey", "red3"))
# 
# mvmt %>% 
#   filter(ID %in% part.sample) %>%
#   ggplot(aes(age, -z, colour=dX>0 & dY>0, group=ID)) + 
#   geom_line() +
#   scale_colour_manual("Uploch", values=c("grey", "red3")) +
#   facet_wrap(~ID)
# 
# mvmt %>% 
#   filter(ID %in% part.sample) %>%
#   ggplot(aes(age, -z, colour=dZ>0, group=ID)) + 
#   geom_line() +
#   facet_wrap(~ID)
# mvmt %>% 
#   filter(ID %in% part.sample) %>%
#   ggplot(aes(age, -z, colour=factor(sink), group=ID)) + 
#   geom_line() +
#   facet_wrap(~ID)
# 
# mvmt %>%
#   mutate(uploch=dX>0 & dY>0) %>%
#   ggplot(aes(dZ>0, fill=uploch)) + geom_bar()
# 
# mvmt %>%
#   mutate(uploch=dX>0 & dY>0) %>%
#   ggplot(aes(uploch, fill=factor(sink))) + geom_bar(position="fill")
# 
# mvmt %>%
#   mutate(uploch=dX>0 & dY>0) %>%
#   group_by(uploch) %>%
#   summarise(prSink=mean(sink))
# 
# mvmt %>%
#   mutate(uploch=dX>0 & dY>0) %>%
#   ggplot(aes(-z, colour=uploch)) + 
#   xlim(-20, 0) +
#   geom_density()






# psteps ------------------------------------------------------------------

# psteps.f <- dir(out.dir, "pstepsMature")[90:91]
# psteps.df <- map_dfr(psteps.f, 
#                      ~paste0(out.dir, .x) %>%
#                        read_delim(delim=" ", col_names=FALSE) %>%
#                        rename(i=X1) %>%
#                        rowwise() %>%
#                        summarise(i=i,
#                                  Tot=sum(c_across(starts_with("X")))) %>%
#                        mutate(date=str_split_fixed(.x, "_", 3)[,2],
#                               time=str_sub(str_split_fixed(.x, "_", 3)[,3], 1, -5)))
# psteps.sf <- mesh.sf %>%
#   inner_join(psteps.df, by="i") %>%
#   mutate(Lice_per_m3=Tot/(area*depth),
#          ln_Lice_per_m3=log(Lice_per_m3))
# ggplot(psteps.sf, aes(fill=ln_Lice_per_m3)) + 
#   geom_sf(data=mesh.fp, fill="grey") +
#   geom_sf(colour=NA) + 
#   scale_fill_viridis_c(limits=c(-15, -7)) +
#   facet_wrap(~time)
# ggplot(psteps.sf, aes(fill=ln_Lice_per_m3)) + 
#   geom_sf(data=mesh.fp, fill="grey") +
#   geom_sf(colour=NA) + 
#   scale_fill_viridis_c() +
#   facet_wrap(~time)
# ggplot(psteps.sf, aes(fill=Lice_per_m3)) + 
#   geom_sf(data=mesh.fp, fill="grey") +
#   geom_sf(colour=NA) + 
#   scale_fill_viridis_c() +
#   facet_wrap(~time)







# locations ---------------------------------------------------------------

loc.f <- dir(out.dir, "locations_")
loc.df <- map_dfr(loc.f, 
                  ~paste0(out.dir, .x) %>%
                    read_delim(delim=" ", col_names=T, show_col_types=F) %>%
                    mutate(date=ymd(str_sub(.x, 11, -5)),
                           dateTime=as_datetime(paste0(str_sub(.x,11,-5), "-", hour),
                                                format="%Y%m%d-%H")))

loc.df <- bind_rows(loc.df %>% mutate(res="1h"),
                    loc.df2 %>% mutate(res="5min"))
ggplot(loc.df, aes(x, y, group=ID, colour=zTot)) + 
  geom_path() +
  scale_colour_viridis_c()


# 
# loc.df %>%
#   ggplot(aes(depth, colour=hour)) + 
#   geom_density() + 
#   facet_wrap(~date) +
#   scale_colour_viridis_c()
# 
# elem.sum <- full_join(
#   mesh.sf %>% select(-starts_with("trinode")),
#   full_join(
#     expand_grid(i=unique(mesh.sf$i),
#                 date=unique(loc.df$date),
#                 hour=unique(loc.df$hour)),
#     loc.df %>% group_by(date, hour, elem) %>%
#       summarise(N=n(), 
#                 DensTot=sum(density),
#                 mnDepth=mean(depth)) %>%
#       rename(i=elem), 
#     by=c("i", "date", "hour")),
#   by="i")
# 
# elem.sum %>%
#   filter(date=="2021-11-05", hour<5) %>%
#   ggplot() +
#   geom_sf(colour=NA, aes(fill=mnDepth)) +
#   scale_fill_viridis_c() +
#   facet_wrap(~hour)


library(gganimate)
anim <- loc.df %>%
  # filter(date=="2021-11-05") %>%
  # filter(hour < 4) %>%
  filter(ID %% 10 == 0) %>%
  st_as_sf(coords=c("x", "y"), crs=27700) %>%
  ggplot() + 
  geom_sf(data=mesh.fp, fill="grey10", colour=NA) +
  geom_sf(aes(colour=depth, alpha=density), size=0.5) + 
  scale_colour_viridis_c(direction=-1, end=0.9, option="D") +
  transition_states(dateTime, wrap=F, transition_length=100, state_length=1) +
  ggtitle(paste0("Hi-res, 1 hour. {closest_state}"))
anim_save("locs_test.gif", anim, nframes=8*24*5, fps=10, 
          width=7, height=6, res=300, units="in")
