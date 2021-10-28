#### NEFSC Lobster Zone Re-stratification
#### 10/21/2021


####  Packages  ####
library(sf)
library(tidyverse)
library(targets)
library(gmRi)

# paths
box_paths <- research_access_paths()
res_path <- box_paths$res

####  Data  ####

# Load clean trawl data, add length weight data
nefsc_clean <- gmRi::gmri_survdat_prep(survdat_source = "most recent")
nefsc_clean <- add_lw_info(nefsc_clean, cutoff = TRUE)


# Lobstah strata
# lobstrata <- read_sf(paste0(res_path, "Shapefiles/Statistical_Areas/Statistical_Areas_2010.shp"))
lobstrata <- read_sf(paste0(res_path, "Shapefiles/Statistical_Areas/Statistical_Areas_2010_withNames.shp"))
lobstrata <- rename_all(lobstrata, tolower)



####  Re-stratify trawl data  ####



# From box/Res_data/NMFS_trawl...
# modified for this repurposing exercise

# Assign statistical zone from new sf for stat zones
assign_stat_zones <- function(survdat, zone_sf, strata_col_in = "id", strata_col_out = "stat_zone", keep_NA = FALSE){
  
  # # Testing:
  # survdat <- nefsc_clean
  # zone_sf <- lobstrata
  # strata_col_in <- "short_name"
  # strata_col_out <- "lobster_strata"
  # keep_NA = FALSE
  
  # transfer to shorthands
  x <- as.data.frame(survdat)
  out_name_sym <- sym(strata_col_out)
  
  # use only station data for overlay/intersection
  stations <- distinct(x, cruise6, stratum, station, decdeg_beglat, decdeg_beglon)
  
  # Convert stations to sf
  stations_sf <- st_as_sf(stations, coords = c("decdeg_beglon", "decdeg_beglat"), crs = 4326, remove = FALSE)
  
  # Project to Lambert Conformal Conic
  lcc <- st_crs("+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-72 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0 ") 
  stations_sf <- st_transform(stations_sf, crs = lcc)
  
  # Prepare statistical zones in same CRS
  stratum <- st_transform(lobstrata, crs = lcc)
  
  # rename stratum column to match desired label
  names(stratum)[which(names(stratum) == strata_col_in)] <- strata_col_out
  
  # Identify points within each polygon/strata
  stations_sf <- st_join(stations_sf, stratum, join = st_within, left = TRUE)
  
  # good to chuck
  # intersect_lst <- st_intersects(stations_sf, stratum) # returns list with unknown number
  # 
  # # checking names? $name key
  # name_key <- stratum %>% filter(rowid %in% as.character(intersect_lst)) %>% 
  #   st_drop_geometry()
  # 
  # 
  # # which points are in each polygon
  # stations_sf <-  stations_sf %>% 
  #   mutate({{out_name_sym}} := as.character(intersect_lst))
  # 
  # # Convert back to lat/lon epsg: 4326 
  # (don't need to because we kept original coords as column)
  # stations_wgs_sf <- st_transform(stations_sf, crs = 4326)
  
  # Don't need to convert back b/c we kept coordinates
  stations_wgs <- st_drop_geometry(stations_sf)
  
  # Keep NA's or not?
  if(keep_NA == F){ stations_wgs <- filter(stations_wgs, is.na({{out_name_sym}}) == FALSE)}
  
  # Join station assignments back into full data
  out <- right_join(stations_wgs, x, by = c('cruise6', 'stratum', 'station', "decdeg_beglat", "decdeg_beglon")) %>% 
    mutate({{out_name_sym}} := as.character({{out_name_sym}}))

  # return the table
  return(out)
  
  }


# Assign those zones!
nefsc_lobsta_zones <- assign_stat_zones(survdat = nefsc_clean, 
                                        zone_sf = select(lobstrata, id, short_name, geometry), 
                                        strata_col_in = "id", 
                                        strata_col_out = "lobster_strata",
                                        keep_NA = FALSE)

# Tidy up?
# nefsc_lobsta_zones <- nefsc_lobsta_zones %>% filter(lobster_strata != "integer(0)")

# Map test
new_england <- rnaturalearth::ne_states("united states of america", returnclass = "sf") 
canada <- rnaturalearth::ne_states("canada", returnclass = "sf") 



####  Filter  Areas to GoM  ####
gom_dat <- nefsc_lobsta_zones %>% 
  filter(survey_area == "GoM",
         lobster_strata %not in% c(561, 522, 551, 526, 466, 467),
         is.na(lobster_strata) == FALSE)

# make sf to check
gom_dat_sf <- gom_dat %>% 
  distinct(decdeg_beglon, decdeg_beglat, lobster_strata) %>% 
  st_as_sf(coords = c("decdeg_beglon", "decdeg_beglat"), crs = 4326) 

# map check
ggplot() +
  geom_sf(data = gom_dat_sf, aes(color = lobster_strata)) +
  geom_sf(data = new_england, size = 0.3) +
  geom_sf(data = canada, size = 0.3) +
  coord_sf(xlim = c(-71, -65.8), ylim = c(41, 44.5)) +
  theme(legend.position = "bottom")



#### Get Stratum Area in km2  ####

# Take lobster strata and get the area of each zone

# Lambert Conformal Conic
lcc <- st_crs("+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-72 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0 ")
lobstrata_lcc <- st_transform(lobstrata, lcc)

# Get areas
# rename the id column to match nefsc_lobsta
# drop geometry
zone_areas <- lobstrata_lcc %>% 
  rename(lobster_strata = id) %>% 
  mutate(lobster_strata = as.character(lobster_strata),
         area_m2 = st_area(lobstrata_lcc),
         area_km2 = units::set_units(area_m2, "km^2")) %>% 
  select(lobster_strata, area_km2) %>% 
  st_drop_geometry()


# merge in the areas
gom_weights <- left_join(gom_dat, zone_areas, by = "lobster_strata")




####  Get Stratified Abundances  ####

stratify_lobster_strata_catch <- function(survdat_weights, area_label, area_col){
  
  # https://noaa-edab.github.io/survdat/articles/calc_strat_mean.html
  
 # # Testing:
 #  area_label <- "lobster_strata"
 #  area_col <- "area_km2"
 label_sym <- sym(area_label)
 areas_sym <- sym(area_col)
  
  ####  1. Set Constants: 
  
  # Area covered by an albatross standard tow in km2
  alb_tow_km2 <- 0.0384
  
  # catchability coefficient - ideally should change for species guilds
  q <- 1
  
  
  
  ####  2. Stratum Area & Effort Ratios 
  # Get Annual Stratum Effort, and Area Ratios
  # The number of tows in each stratum by year
  # area of a stratum relative to total area of all stratum sampled that year
  
  # Get Total area of all strata sampled in each year
  total_stratum_areas <- dplyr::group_by(survdat_weights, est_year)
  total_stratum_areas <- dplyr::distinct(total_stratum_areas, {{label_sym}}, .keep_all = T)
  total_stratum_areas <- dplyr::summarise(total_stratum_areas,
                                          tot_s_area =  sum({{areas_sym}}, na.rm = T),
                                          .groups = "drop")
  
  
  # Calculate strata area relative to total area i.e. stratio or stratum weights
  survdat_weights <- dplyr::left_join(survdat_weights, total_stratum_areas, by = "est_year")
  survdat_weights <- dplyr::mutate(survdat_weights, 
                                   st_ratio = {{areas_sym}} / tot_s_area,
                                   st_ratio = as.numeric(st_ratio))
  
  
  # We have total areas, now we want effort within each
  # Number of unique tows per stratum, within each season
  yr_strat_effort <- dplyr::group_by(survdat_weights, est_year, season, {{label_sym}})
  yr_strat_effort <- dplyr::summarise(yr_strat_effort, strat_ntows = dplyr::n_distinct(id), .groups = "drop")
  
  
  
  # Add those yearly effort counts back for later
  # (area stratified abundance)
  survdat_weights <- dplyr::left_join(survdat_weights, yr_strat_effort, by = c("est_year", "season", area_label))
  
  
  
  
  ####  4. Derived Stratum Area Estimates
  
  # a. Catch / tow, for that year & season
  survdat_weights <-  dplyr::mutate(survdat_weights,
                                    
                                    # Abundance
                                    abund_tow_s   = numlen_adj / strat_ntows,
                                    
                                    # All size biomass
                                    # Biomass is repeated across length classes at each station by species
                                    # the number of length classes is tallied where the conversion factor is done
                                    biom_per_lclass = (biomass_kg / n_len_class),
                                    biom_tow_s = biom_per_lclass / strat_ntows,
                                    
                                    # Length specific biomass
                                    lwbio_tow_s = sum_weight_kg / strat_ntows)
  
  
  # b. Stratified Mean Catch Rates
  survdat_weights <-  dplyr::mutate(survdat_weights,
                                    
                                    # Stratified mean abundance CPUE, weighted by the stratum areas
                                    strat_mean_abund_s = abund_tow_s * st_ratio,
                                    
                                    # Stratified mean BIOMASS CPUE
                                    strat_mean_biom_s = biom_tow_s * st_ratio,
                                    
                                    # Stratified mean LW Biomass
                                    strat_mean_lwbio_s = lwbio_tow_s * st_ratio)
  
  
  # c. Stratified Total Abundance/Biomass
  # convert from catch rate by area swept to total catch for entire stratum
  survdat_weights <-  dplyr::mutate(survdat_weights,
                                    
                                    # Total Abundance
                                    strat_total_abund_s = round((strat_mean_abund_s * tot_s_area / alb_tow_km2) / q),
    
                                    # Total BIOMASS from the biomass of all lengths
                                    strat_total_biom_s = (strat_mean_biom_s * tot_s_area / alb_tow_km2) / q,
                                    
                                    # Two options for to estimate lw biomass
                                    # Result is the same 4/20/2021
                                    # Option 1: Individual LW Biomass * expanded abundance at length
                                    strat_total_lwbio_s = (ind_weight_kg * strat_total_abund_s) / q
                                    
                                    # # Option 2: Size specific lw biomass / tow, expanded to total area
                                    # strat_total_lwbio_s  = (strat_mean_lwbio_s * tot_s_area / alb_tow_km2) / q
  )
  
  
  
  
 
  
  
  
}


# Run stratification
gom_strat <- stratify_lobster_strata_catch(survdat_weights = gom_weights, 
                                           area_label = "lobster_strata", 
                                           area_col = "area_km2")


# Tidy up?
# there are now two different "stratum columns" floating around
gom_strat <- gom_strat %>% 
  select(-c(strat_num, stratum, nafodiv, full_name, st_ratio, strat_ntows, biom_per_lclass)) 

gom_strat <- gom_strat %>% 
  mutate(strat_total_abund_s = as.numeric(strat_total_abund_s), 
         strat_total_biom_s = as.numeric(strat_total_biom_s), 
         strat_total_lwbio_s = as.numeric(strat_total_lwbio_s))



#### Pull Predators  ####


# Pull out species we care about
lobster_predators <- c(
  "atlantic halibut",
  "atlantic wolffish",
  "barndoor skate",
  "black sea bass",
  "atlantic cod",
  "fourspot flounder",
  "haddock",
  "little skate",
  "longhorn sculpin",
  "ocean pout",
  "red hake",
  "sea raven",
  "silver hake",
  "smooth skate",
  "spiny dogfish",
  "spotted hake",
  "thorny skate",
  "white hake",
  "winter flounder"
)


# Filter
gom_predators <- gom_strat %>% 
  filter(comname %in% lobster_predators) %>% 
  mutate(lobster_pred = "lobster predator")

# check
gom_predators %>% distinct(comname) %>% arrange(comname)



####  Assess Predator Indices  ####
# get stratified mean abundances
predator_summs <- gom_predators %>% 
  bind_rows(gom_predators %>% mutate(season = "Both")) %>% 
  group_by(est_year, season, lobster_pred) %>% 
  summarise(
    # station details
    `station total` = n_distinct(id),
    `total fish` = sum(numlen_adj),
    # Survey totals
    `total biomass` = sum(sum_weight_kg),
    `average abundance` = mean(numlen_adj),
    #  Stratified Means
    # Should these averages be weighted by the number of fish caught?
    `stratified mean abundance` = mean(strat_mean_abund_s),
    `stratified mean biomass (kg)` = mean(strat_mean_lwbio_s),
    `strat-abund weighted mean abundance` = weighted.mean(strat_mean_abund_s, strat_total_abund_s),
    `strat-abund weighted mean biomass (kg)` = weighted.mean(strat_mean_lwbio_s, strat_total_abund_s),
    # Weighted Average Sizes
    `stratified mean length (cm)` = weighted.mean(length_cm, strat_total_abund_s),
    `stratified mean weight (kg)` = weighted.mean(ind_weight_kg, strat_total_abund_s),
    # Stratified Totals
    `stratified total abundance` = sum(strat_total_abund_s),
    `stratified total biomass (kg)` = sum(strat_total_lwbio_s),
            .groups = "drop") %>% 
  mutate(season = factor(season, levels = c("Spring", "Fall", "Both")))


# standardize? No, Matt will have to anyways



####  Data Exploration  ####

# Stop re-typing code
plot_pred_index <- function(pred_dat, col_select = "stratified total biomass (kg)"){
  col_sym <- sym(col_select)
  ggplot(pred_dat, aes(x = est_year, y = {{col_sym}}, color = season)) +
    geom_line(alpha = 0.2) +
    geom_point(alpha = 0.6) +
    geom_smooth(formula = y~x, method = "loess") + 
    scale_color_gmri() +
    facet_wrap(~season, ncol = 1, scales = "free") + 
    scale_y_continuous(labels = scales::comma_format()) +
    labs(x = "", y = str_replace_all(col_select, "_", " "),
         subtitle = str_c("Lobster Predator Complex: ", col_select)) +
    theme(legend.position = "right") +
    guides(color = guide_legend(title = "Season", title.position = "top", title.hjust = 0.5))
  
}


# Plots
plot_pred_index(predator_summs, col_select = "stratified total biomass (kg)")
plot_pred_index(predator_summs, col_select = "strat_mean_len_cm")

# Fish silhouettes?
