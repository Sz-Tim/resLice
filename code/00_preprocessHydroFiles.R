# Hydrofile names
# Off Aqua
# Tim Szewczyk



# setup -------------------------------------------------------------------

library(tidyverse); library(glue); library(lubridate); library(ncdf4)

dirs <- list(
  linnhe_1h="D:/hydroOut/linnhe7/linnhe7_tides_met_tsobc_riv/netcdf_2021",
  linnhe_5min="D:/hydroOut/linnhe7/linnhe7_tides_met_tsobc_riv_5min/orig",
  linnhe_5min_new="D:/hydroOut/linnhe7/linnhe7_tides_met_tsobc_riv_5min/netcdf_2021",
  WeStCOMS2_1h="D:/hydroOut/WeStCOMS2/Archive/netcdf_2021"
)

startDate <- ymd("2021-11-01")




# rename ------------------------------------------------------------------

f.orig <- dir(dirs$linnhe_5min)
dateSeq <- str_remove_all(startDate + (seq_along(f.orig)-1), "-")

# RUN TO RENAME MAX'S FILES
# f.renamed <- str_replace(f.orig, "_", glue("_{dateSeq}_"))
# file.rename(glue("{dir.5min}/{f.orig}"), glue("{dir.5min}/{f.renamed}"))

# RUN TO RENAME DIMA'S RE-RUN FILES
# f.renamed <- str_replace(f.orig, "_202111[0-1][0-9]_", glue("_{dateSeq}_"))
# file.copy(glue("{dirs$linnhe_5min}/{f.orig}"), glue("{dirs$linnhe_5min_new}/{f.renamed}"))




# filter 5min to 1h -------------------------------------------------------
extract.i <- tibble(file_5min=dir(dirs$linnhe_5min_new)) %>%
  mutate(cumul_hours=row_number()*2,
         day=(cumul_hours-1) %/% 24,
         date=str_remove_all(startDate + day, "-"),
         layer_start_1h=(cumul_hours-1) %% 24,
         layer_end_1h=layer_start_1h+1,
         file_1h=glue("linnhe7_{date}_0001.nc"))

# In Dima's files, layer 0 and layer 12 are rounded hours
# I need to create a new file for each day
init.nc <- nc_open(glue("{dirs$linnhe_5min_new}/{extract.i$file_5min[1]}"))
# ncdf attributes and mesh info
nc <- list(
  uvnode=cbind(ncvar_get(init.nc, "lonc"), ncvar_get(init.nc, "latc")),
  nodexy=cbind(ncvar_get(init.nc, "lon"), ncvar_get(init.nc, "lat")),
  siglay=ncvar_get(init.nc, "siglay")[1,], 
  siglev=ncvar_get(init.nc, "siglev")[1,]
)
elems <- ncdim_def("elems", "element", 1:nrow(nc$uvnode))
two <- ncdim_def("two", "two", 1:2)
node <- ncdim_def("node", "nodes", 1:nrow(nc$nodexy))
three <- ncdim_def("three", "three", 1:3)
siglayers <- ncdim_def("siglayers", "sigma layers", 1:length(nc$siglay))
siglevels <- ncdim_def("siglevels", "sigma levels", 1:length(nc$siglev))

dims.i <- tribble(~name, ~len, ~obj,
                  "elems", elems$len, elems,
                  "node", node$len, node,
                  "siglayers", siglayers$len, siglayers,
                  "siglevels", siglevels$len, siglevels,
                  "Time", 24, ncdim_def("Time", "Time", 1:24), 
                  "three", 3, three,
                  "nine", 9, ncdim_def("nine", "nine", 1:9))

# hydro variables
var.needed <- c("lon", "lat", "lonc", "latc", "art1",
                "h_center", "h", "Times", "zeta", "siglay", "siglev",
                "u", "v", "ua", "va", "ww", "omega", "temp", "salinity", 
                "km", "short_wave", "uwind_speed", "vwind_speed",
                "net_heat_flux")
var.toUpdate <- c("Times", "zeta", "u", "v", "ww", "temp", "salinity", "km", "short_wave")
var.init.ls <- map(var.needed, ~ncvar_get(init.nc, .x)) %>% setNames(var.needed)
var.att.ls <- map(var.needed, ~ncatt_get(init.nc, .x)) %>% setNames(var.needed)
nc_close(init.nc)
for(i in 1:n_distinct(extract.i$file_1h)) {
  new.f <- unique(extract.i$file_1h)[i]
  extract.i_i <- extract.i %>% filter(file_1h==new.f)
  vars.hour_i <- var.init.ls
  ncvar.i <- vector("list", length(vars.hour_i)) %>% setNames(names(vars.hour_i))
  
  for(j in 1:nrow(extract.i_i)) {
    nc.j <- nc_open(glue("{dirs$linnhe_5min_new}/{extract.i_i$file_5min[j]}"))
    i_j <- extract.i_i$layer_start_1h[j]:extract.i_i$layer_end_1h[j]
    for(v in var.needed) {
      dims.v <- dim(vars.hour_i[[v]])
      ndims <- length(dims.v)
      time.dim <- which(dims.v==24)
      if(v %in% var.toUpdate) {
        if(ndims==1) {
          vars.hour_i[[v]][i_j] <- ncvar_get(nc.j, v)[c(1,13)]
        }
        if(ndims==2) {
          if(time.dim==1) {
            vars.hour_i[[v]][i_j,] <- ncvar_get(nc.j, v)[c(1,13),]
          } else {
            vars.hour_i[[v]][,i_j] <- ncvar_get(nc.j, v)[,c(1,13)]
          }
        }
        if(ndims==3) {
          if(time.dim==1) {
            vars.hour_i[[v]][i_j,,] <- ncvar_get(nc.j, v)[c(1,13),,]
          } else if(time.dim==2) {
            vars.hour_i[[v]][,i_j,] <- ncvar_get(nc.j, v)[,c(1,13),]
          } else {
            vars.hour_i[[v]][,,i_j] <- ncvar_get(nc.j, v)[,,c(1,13)]
          }
        }
      }
    }
    nc_close(nc.j)
  }
  vars.hour_i[["Times"]] <- str_replace(vars.hour_i[["Times"]], "T", " ") %>%
    as_datetime() %>% as.numeric() %>% array()
  for(v in var.needed) {
    dims.v <- dim(vars.hour_i[[v]])
    ncvar.i[[v]] <- ncvar_def(name=v, 
                              units=ifelse(is.null(var.att.ls[[v]]$units), 
                                           "none", 
                                           var.att.ls[[v]]$units), 
                              dim=dims.i$obj[match(dims.v, dims.i$len)])
  }
  
  hour.nc <- nc_create(glue("{dirs$linnhe_1h}/{new.f}"), ncvar.i)
  walk(seq_along(ncvar.i), ~ncvar_put(hour.nc, var.needed[.x], vars.hour_i[[.x]]))
  nc_close(hour.nc)
}








# key ---------------------------------------------------------------------

getTimes <- function(nc, unit="s") {
  library(ncdf4)
  nc_temp <- nc_open(nc)
  times <- c(ncvar_get(nc_temp, "Times"))
  nc_close(nc_temp)
  return(times)
}


time.df <- list(
  map_dfr(dir(dirs$linnhe_5min_new), 
          ~tibble(file=.x,
                  fileDate=str_split_fixed(.x, "_", 3)[,2], 
                  time_raw=getTimes(glue("{dirs$linnhe_5min_new}/{.x}"))) %>%
            mutate(layer=row_number()-1,
                   time_dt=str_replace(time_raw, "T", " ") %>%
                     as_datetime())) %>%
    mutate(mesh="linnhe7", 
           timeRes="5min"),
  map_dfr(dir(dirs$linnhe_1h), 
          ~tibble(file=.x,
                  fileDate=str_split_fixed(.x, "_", 3)[,2], 
                  time_raw=getTimes(glue("{dirs$linnhe_1h}/{.x}"))) %>%
            mutate(layer=row_number()-1,
                   time_dt=as_datetime(time_raw))) %>%
    mutate(mesh="linnhe7", 
           timeRes="1h"),
  map_dfr(dir(dirs$WeStCOMS2_1h), 
          ~tibble(file=.x,
                  fileDate=str_split_fixed(.x, "_", 3)[,2], 
                  time_raw=getTimes(glue("{dirs$WeStCOMS2_1h}/{.x}"))) %>%
            mutate(layer=row_number()-1,
                   time_dt=str_replace(time_raw, "T", " ") %>%
                     as_datetime())) %>%
    mutate(mesh="WeStCOMS2", 
           timeRes="1h")
) %>%
  do.call('rbind', .) %>%
  mutate(date=date(time_dt)) %>%
  group_by(mesh, timeRes) %>%
  arrange(time_dt) %>%
  mutate(timeStep=difftime(time_dt, lag(time_dt), units="mins"),
         i=row_number()-1,
         timeCalculated=round_date(time_dt, "5 minute"))

write_csv(time.df, "data/timeRes_key.csv")

