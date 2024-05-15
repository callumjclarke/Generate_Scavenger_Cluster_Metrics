library('move2')
library('lubridate')
library('magrittr')
library('dplyr')
library('sf')
library('units')

# Shortened 'not in':
`%!in%` <- Negate(`%in%`)

# Efficient interval function
int_overlaps_numeric <- function (int1, int2) {
  stopifnot(c(is.interval(int1), is.interval(int2)))
  
  x <- dplyr::intersect(int1, int2)@.Data
  x[is.na(x)] <- 0
  as.duration(x)
}

# Function to calculate geometric medians:
calcGMedianSF <- function(data) {
  
  if (st_geometry_type(data[1,]) == "POINT") {
    med <- data %>% 
      st_coordinates()
    med <- Gmedian::Weiszfeld(st_coordinates(data))$median %>% as.data.frame() %>%
      rename(x = V1, y = V2) %>%
      st_as_sf(coords = c("x", "y"), crs = st_crs(data)) %>%
      st_geometry()
  }
  
  
  if (st_geometry_type(data[1,]) == "MULTIPOINT") {
    med <- data %>% 
      st_coordinates() %>%
      as.data.frame() %>%
      group_by(L1) %>%
      group_map(
        ~ Gmedian::Weiszfeld(.)$median 
      ) %>%
      do.call(rbind, .) %>%
      as.data.frame() %>%
      st_as_sf(coords = colnames(.), crs = st_crs(data)) %>%
      st_geometry()
  }
  
  return(med)
  
}


# Showcase injecting app setting (parameter `year`)
rFunction = function(data, 
                     sleepduration,
                     sleepdist,
                     stationary_dist
                     ) {
  

  
  # ///////////////////////////////////////////////////
  # 1. Check that xy.clust is in data (write later) + other input checks  ----
  # ///////////////////////////////////////////////////
  
  logger.info("[1] Running input checks")
  logger.trace(paste0(
    "   Data contains ",
    nrow(data), 
    " entries"
  ))
  
  
  if (any(c("sunrise_timestamp", "sunset_timestamp") %!in% colnames(data))) {
    logger.fatal("Input data does not contain sunrise-sunset data. Please use 'Add Local and Solar Time' MoveApp before this stage in the workflow")
    stop()
  }
  
  if (!mt_is_move2(data)) {
    logger.fatal("Input data is not a move2 object. Please ensure input is in move2 format")
    stop()
  }
  
  
  
  # ///////////////////////////////////////////////////
  # 2. Group by cluster, calculate geometric median ----
  # ///////////////////////////////////////////////////
  
  logger.info("[2] Calculating geometric median")
  
  newdat <- data %>%
    mutate(ID = mt_track_id(.), 
           timestamp = mt_time(.),
           timediff_hrs = mt_time_lags(.) %>% 
             units::set_units("minutes") %>% 
             units::set_units("hours"),
           
           # Sunrise calculations:
           sunint = interval(sunrise_timestamp, sunset_timestamp),
           # Find hours of overlap with daytime for each location
           dayhrs = int_overlaps_numeric(
             int1 = interval(timestamp, lead(timestamp)),
             int2 = sunint
           ) %>%
             units::set_units("seconds") %>%
             units::set_units("hours")
    ) 
  
  clustertable <- newdat %>%
    filter(!is.na(xy.clust)) %>%
    group_by(xy.clust) %>%
    summarise(
      geometry = st_union(geometry),
      all_cluster_points = geometry,
      
      ntags = length(unique(ID)),
      tags = list(unique(ID)),
      
      # Time calculations:
      firstdatetime = min(timestamp),
      lastdatetime = max(timestamp),
      days = length(unique(lubridate::date(timestamp))),
      totaldays = (ceiling_date(lastdatetime, unit = "days") - floor_date(firstdatetime, unit = "days")) %>% as.integer(),
      daysempty = as.numeric(totaldays - days),
      daysemptyprop = daysempty / totaldays,
      
      # Day-night calculations:
      # Calculate dayprop by calculating the proportion
      # of each event occuring in daytime
      Total = n(),
      totalhrs = sum(timediff_hrs, na.rm = T),
      totalhrs_daytime = sum(dayhrs, na.rm = T),
      dayprop = totalhrs_daytime / totalhrs,
      
    ) %>%
    
    # Set cluster location to geometric median:
    st_set_geometry(
      calcGMedianSF(.) %>%
        st_geometry()
    )
  
  
  
  # ///////////////////////////////////////////////////
  # 3. Identify daytime-sleeping locations ----
  # ///////////////////////////////////////////////////
  
  logger.info("[3] Generating sleeping locations")
  
  stat <- newdat %>%
    mutate(stationary = ifelse(
      as.numeric(dist_m) < stationary_dist, 1, 0
    ),
    statruns = data.table::rleid(stationary),
    statruns = ifelse(stationary == 0, NA, statruns),
    statruns = ifelse(
      lag(stationary) == 0 & lead(stationary) == 0,
      NA, statruns
    )
    ) %>%
    filter(!is.na(statruns)) %>%
    group_by(statruns) %>%
    summarise(
      geometry = st_union(geometry),
      geometry = st_centroid(geometry),
      
      statstart = min(timestamp, na.rm = T),
      statend = max(timestamp, na.rm = T),
      statdur = difftime(statend, statstart, units = "hours") 
    ) %>%
    filter(date(statstart) == date(statend),
           statdur > 4 # change this to a user setting
    )
  
  clustcent <- clustertable %>%
    select(c("xy.clust", "geometry")) %>%
    st_buffer(250) %>%
    st_intersects(stat) %>%
    as.data.frame() %>%
    left_join(
      clustertable %>% 
        mutate(row.id = row_number()) %>%
        select(c("xy.clust", "row.id")),
      by = "row.id"
    ) %>%
    as.data.frame() %>%
    group_by(xy.clust) %>%
    summarise(
      sleepspots = n()
    ) 
  clustertable %<>% 
    left_join(clustcent, by = "xy.clust") %>%
    mutate(
      sleepspots = ifelse(is.na(sleepspots), 0, sleepspots)
    )
  
  
  
  # ///////////////////////////////////////////////////
  # 4. Generate visit duration data -----
  # ///////////////////////////////////////////////////
  
  logger.info("[4] Generating visit-duration data")
  
  clusttimes <- newdat %>%
    as.data.frame() %>%
    group_by(ID) %>%
    mutate(
      clustrun = data.table::rleid(xy.clust)
    ) %>%
    filter(!is.na(xy.clust)) %>%
    group_by(ID, clustrun) %>%
    summarise(
      xy.clust = first(xy.clust),
      visitdur = sum(timediff_hrs, na.rm = T),
      .groups = "drop"
    ) %>%
    group_by(xy.clust) %>%
    summarise(
      mean_visit_duration = mean(visitdur, na.rm = T)
    )
  clustertable %<>% 
    left_join(clusttimes, by = "xy.clust")
  
  
  # ///////////////////////////////////////////////////
  # 5. Generate revisit data -----
  # ///////////////////////////////////////////////////
  
  logger.info("[5] Generating revisit-frequency data")
  
  newrevs <- newdat %>%
    group_by(ID) %>%
    mutate(firstloc = case_when(
      # Identify 'first location' in each cluster visit
      # i.e. in cluster but previous location isn't 
      !is.na(xy.clust) &
        (lag(xy.clust) != xy.clust | 
           is.na(lag(xy.clust))) ~ 1,
      TRUE ~ 0
    )) %>%
    
    # Now calculate revisit times 
    # Remove non-cluster locations:
    filter(!is.na(xy.clust)) %>%
    group_by(ID, xy.clust) %>%
    mutate(
      newlags = difftime(timestamp, lag(timestamp), units = "hours")
    ) %>%
    filter(firstloc == 1) %>%
    as.data.frame() %>%
    group_by(xy.clust) %>%
    summarise(
      nrevisits = n() - 1,
      mean_revisit_time = mean(newlags, na.rm = T),
    )
  clustertable %<>% 
    left_join(newrevs, by = "xy.clust")
  
  
  # ///////////////////////////////////////////////////
  # 6. Returning output -----
  # ///////////////////////////////////////////////////
  
  logger.info("[6] Returning output")
  logger.trace(paste0("   Clusters identified: ", nrow(clustertable)))
  
  # Return clustertable:
  attr(data, "cluster_tbl") <- clustertable
  return(data)
  
  
}
