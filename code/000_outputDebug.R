# Mesh exploration
# Off Aqua
# Tim Szewczyk

# This script is to debug / test output options



# setup -------------------------------------------------------------------

library(tidyverse); library(ncdf4); library(sf); library(glue)
theme_set(theme_bw())

# directories
out.dir <- "../sealice_runs/00_test/out_debug/"

# files
mesh.nc <- 
start.i <- read_delim(glue("{out.dir}startSitesUsed.dat"))
elemActivity <- read_delim(glue("{out.dir}elementActivity.dat"), 
                           col_names=c("sink", "swim", "float", "total")) %>%
  mutate(total=sink+swim+float,
         i=row_number())
mvmt <- read_delim(glue("{out.dir}movementFile.dat"))



# start sites -------------------------------------------------------------

ggplot(start.i, aes(x, y, colour=salinity)) + 
  geom_point() + scale_colour_viridis_c()




# movement ----------------------------------------------------------------

part.sample <- sample(unique(mvmt$ID), 36)
part.sample <- sample(filter(mvmt, startDate=="20211101")$ID, 200)

mvmt %>% 
  filter(ID %in% part.sample) %>%
  ggplot(aes(x, y, colour=density, group=ID)) + 
  geom_path() + scale_colour_viridis_c()

mvmt %>% 
  filter(ID %in% part.sample) %>%
  ggplot(aes(x, y, colour=temp-tempSurface, group=ID)) + 
  geom_path() + scale_colour_gradient2(low=scales::muted("blue"), high=scales::muted("red"))

mvmt %>% 
  filter(ID %in% part.sample) %>%
  ggplot(aes(x, y, colour=dX>0 & dY>0, group=ID)) + 
  geom_path() +
  scale_colour_manual("Uploch", values=c("grey", "red3"))

mvmt %>% 
  filter(ID %in% part.sample) %>%
  ggplot(aes(age, -z, colour=dX>0 & dY>0, group=ID)) + 
  geom_line() +
  scale_colour_manual("Uploch", values=c("grey", "red3")) +
  facet_wrap(~ID)

mvmt %>% 
  filter(ID %in% part.sample) %>%
  ggplot(aes(age, -z, colour=dZ>0, group=ID)) + 
  geom_line() +
  facet_wrap(~ID)
mvmt %>% 
  filter(ID %in% part.sample) %>%
  ggplot(aes(age, -z, colour=factor(sink), group=ID)) + 
  geom_line() +
  facet_wrap(~ID)

mvmt %>%
  mutate(uploch=dX>0 & dY>0) %>%
  ggplot(aes(dZ>0, fill=uploch)) + geom_bar()

mvmt %>%
  mutate(uploch=dX>0 & dY>0) %>%
  ggplot(aes(uploch, fill=factor(sink))) + geom_bar(position="fill")

mvmt %>%
  mutate(uploch=dX>0 & dY>0) %>%
  group_by(uploch) %>%
  summarise(prSink=mean(sink))

mvmt %>%
  mutate(uploch=dX>0 & dY>0) %>%
  ggplot(aes(-z, colour=uploch)) + 
  xlim(-20, 0) +
  geom_density()
