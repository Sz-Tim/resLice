# Hydrofile names
# Off Aqua
# Tim Szewczyk



# setup -------------------------------------------------------------------

library(tidyverse); library(glue); library(lubridate); library(ncdf4)

dirs <- list(
  linnhe_1h="D:/hydroOut/linnhe7/linnhe7_tides_met_tsobc/netcdf_2021",
  linnhe_5min="D:/hydroOut/linnhe7/linnhe7_tides_met_tsobc_5min/netcdf_2021",
  WeStCOMS2_1h="D:/hydroOut/WeStCOMS2/Archive/netcdf_2021"
)

startDate <- ymd("2021-11-01")




# rename ------------------------------------------------------------------

f.orig <- dir(dirs$linnhe_5min)
dateSeq <- str_remove_all(startDate + (seq_along(f.orig)-1), "-")

# RUN TO RENAME FILES
# f.renamed <- str_replace(f.orig, "_", glue("_{dateSeq}_"))
# file.rename(glue("{dir.5min}/{f.orig}"), glue("{dir.5min}/{f.renamed}"))




# key ---------------------------------------------------------------------

getTimes <- function(nc, unit="s") {
  library(ncdf4)
  nc_temp <- nc_open(nc)
  times <- ncvar_get(nc_temp, "time")
  nc_close(nc_temp)
  if(unit=="s") {
    times <- times*24*60*60
  }
  return(times)
}


time.df <- list(
  map_dfr(dir(dirs$linnhe_5min), 
          ~tibble(file=.x,
                  fileDate=str_split_fixed(.x, "_", 3)[,2], 
                  time_raw=getTimes(glue("{dirs$linnhe_5min}/{.x}"))) %>%
            mutate(layer=row_number()-1)) %>%
    mutate(mesh="linnhe7", 
           timeRes="5min"),
  map_dfr(dir(dirs$linnhe_1h), 
          ~tibble(file=.x,
                  fileDate=str_split_fixed(.x, "_", 3)[,2], 
                  time_raw=getTimes(glue("{dirs$linnhe_1h}/{.x}"))) %>%
            mutate(layer=row_number()-1)) %>%
    mutate(mesh="linnhe7", 
           timeRes="1h"),
  map_dfr(dir(dirs$WeStCOMS2_1h), 
          ~tibble(file=.x,
                  fileDate=str_split_fixed(.x, "_", 3)[,2], 
                  time_raw=getTimes(glue("{dirs$WeStCOMS2_1h}/{.x}"))) %>%
            mutate(layer=row_number()-1)) %>%
    mutate(mesh="WeStCOMS2", 
           timeRes="1h")
) %>%
  do.call('rbind', .) %>%
  mutate(time_raw_dt=as.POSIXct(time_raw, origin="1858-11-17"), 
         date=date(time_raw_dt)) %>%
  group_by(mesh, timeRes) %>%
  arrange(time_raw_dt) %>%
  mutate(timeStep=difftime(time_raw_dt, lag(time_raw_dt), units="mins"),
         i=row_number()-1,
         timeCalculated=as_datetime("2021-11-01 00:00:00") + round(mean(timeStep, na.rm=T))*i)

write_csv(time.df, "data/timeRes_key.csv")

