# Initial lice densities
# Off Aqua
# Tim Szewczyk

# This script uses data gathered from http://aquaculture.scotland.gov.uk/ to 
# estimate the daily release of lice from each site

# 



# setup -------------------------------------------------------------------

library(tidyverse); library(sf); library(lubridate); library(zoo)

mesh.fp <- st_read("../resLice/data/linnhe_mesh_footprint.gpkg")

sites.i <- read_csv("../resLice/data/ms_site_details.csv") %>%
  janitor::clean_names(case="small_camel") %>%
  st_as_sf(coords=c("easting", "northing"), crs=27700, remove=F) %>%
  filter(st_within(., mesh.fp, sparse=F))



# initial datasets --------------------------------------------------------

lice.i <- read_csv("../resLice/data/Sealice - SEALICE.csv") %>%
  janitor::clean_names(case="small_camel") %>%
  filter(siteNo %in% sites.i$marineScotlandSiteId &
           weekBeginning > "2021-10-01" &
           weekBeginning < "2021-12-31") %>%
  right_join(sites.i, ., by=c("siteName", "marineScotlandSiteId"="siteNo")) %>%
  group_by(marineScotlandSiteId) %>%
  arrange(marineScotlandSiteId, weekBeginning) %>%
  mutate(weeklyAverageAf=as.numeric(na.fill(zoo(weeklyAverageAf), "extend"))) %>%
  ungroup

fish.i <- read_csv("../resLice/data/biomass_monthly_reports.csv") %>%
  janitor::clean_names(case="small_camel") %>%
  filter(year %in% paste0("01-", c("Oct", "Nov", "Dec"), "-2021")) %>%
  st_as_sf(coords=c("easting", "northing"), crs=27700, remove=F) %>%
  filter(st_within(., mesh.fp, sparse=F)) %>%
  mutate(marineScotlandSiteId=sites.i$marineScotlandSiteId[st_nearest_feature(., sites.i)])

combo.df <- full_join(lice.i %>% st_drop_geometry() %>%
                        mutate(year=year(weekBeginning), 
                               month=month(weekBeginning),
                               day=day(weekBeginning)) %>%
                        select(marineScotlandSiteId, siteName, 
                               year, month, day, 
                               weeklyAverageAf), 
                      fish.i %>% st_drop_geometry() %>%
                        mutate(date=dmy(year),
                               year=year(date),
                               month=month(date)) %>%
                        select(marineScotlandSiteId, sepaSite, siteName,
                               year, month,
                               actualBiomassOnSiteTonnes)) %>%
  mutate(week=week(ymd(paste(year, month, day, sep="-")))) %>%
  group_by(marineScotlandSiteId) %>%
  filter(any(actualBiomassOnSiteTonnes > 0)) %>%
  ungroup



# dates -------------------------------------------------------------------

dates.df <- expand_grid(mssID=unique(combo.df$marineScotlandSiteId),
                        date=seq(ymd("2021-10-01"), ymd("2021-12-31"), by=1)) %>%
  mutate(year=year(date), month=month(date), day=day(date))

# If no values for specific date, use weekly average, then monthly average
# If no value for specific site, use global average
out.df <- full_join(dates.df, combo.df, 
                    by=c("mssID"="marineScotlandSiteId", "year", "month", "day")) %>%
  arrange(mssID, siteName) %>%
  group_by(mssID) %>%
  mutate(siteName=first(siteName),
         sepaSite=first(sepaSite)) %>%
  arrange(mssID, year, month, day) %>%
  group_by(mssID) %>%
  mutate(weeklyAverageAf=as.numeric(na.fill(zoo(weeklyAverageAf), c(NA,"extend",NA)))) %>% 
  group_by(year, month, day) %>%
  mutate(weeklyAverageAf=if_else(is.na(weeklyAverageAf), mean(weeklyAverageAf, na.rm=T), weeklyAverageAf)) %>%
  group_by(year, week) %>%
  mutate(weeklyAverageAf=if_else(is.na(weeklyAverageAf), mean(weeklyAverageAf, na.rm=T), weeklyAverageAf)) %>%
  group_by(year, month) %>%
  mutate(weeklyAverageAf=if_else(is.na(weeklyAverageAf), mean(weeklyAverageAf, na.rm=T), weeklyAverageAf)) %>%
  group_by(mssID, year, month) %>%
  mutate(actualBiomassOnSiteTonnes=mean(actualBiomassOnSiteTonnes, na.rm=T)) %>%
  ungroup %>%
  mutate(weeklyAverageAf=if_else(weeklyAverageAf==0, 0.0001, weeklyAverageAf),
         liceScale=actualBiomassOnSiteTonnes*28.2*240 * weeklyAverageAf) 



# merge and save ----------------------------------------------------------

out.df %>%
  filter(date > "2021-11-09", date < "2021-12-01") %>%
  group_by(mssID) %>%
  filter(any(actualBiomassOnSiteTonnes > 0)) %>%
  ungroup %>%
  mutate(liceScale=round(liceScale/24, 3)) %>%
  mutate(date.c=str_remove_all(date, "-")) %>%
  select(mssID, date.c, liceScale) %>%
  pivot_wider(names_from=date.c, values_from=liceScale) %>%
  write_tsv("../resLice/data/liceScale_daily_1h.tsv", col_names=F)

out.df %>%
  filter(date > "2021-10-31", date < "2021-11-10") %>%
  group_by(mssID) %>%
  filter(any(actualBiomassOnSiteTonnes > 0)) %>%
  ungroup %>%
  mutate(liceScale=round(liceScale/24, 3)) %>%
  mutate(date.c=str_remove_all(date, "-")) %>%
  select(mssID, date.c, liceScale) %>%
  pivot_wider(names_from=date.c, values_from=liceScale) %>%
  write_tsv("../resLice/data/liceScale_daily_1h.tsv", col_names=F)

map_dfr(1:12, 
        ~out.df %>%
          filter(date > "2021-10-31", date < "2021-11-10") %>%
          group_by(mssID) %>%
          filter(any(actualBiomassOnSiteTonnes > 0)) %>%
          ungroup %>%
          mutate(liceScale=round(liceScale/(24*12), 3)) %>%
          mutate(date.c=paste0(str_remove_all(date, "-"), ".", .x)) %>%
          select(mssID, date.c, liceScale)) %>%
  group_by(mssID) %>%
  mutate(date=ymd("2021-10-31")+row_number()) %>%
  ungroup %>%
  mutate(date.c=str_remove_all(date, "-")) %>%
  select(-date) %>%
  pivot_wider(names_from=date.c, values_from=liceScale) %>%
  write_tsv("../resLice/data/liceScale_daily_5min.tsv", col_names=F)

sites.i %>% st_drop_geometry() %>%
  filter(marineScotlandSiteId %in% unique(read_tsv("../resLice/data/liceScale_daily_1h.tsv", col_names=F)[[1]])) %>%
  select(marineScotlandSiteId, easting, northing) %>%
  write_tsv("../resLice/data/fishFarmSites.tsv", col_names=F)



