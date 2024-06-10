# IWD_analysis.R
#------------------------------------------------------------------------------#
#           Inverse-distance weighted smoothing - air pollution data           #
#                                                                              #
# This script generates weekly smoothed values of air pollution data from      #
# Mexico City using a IDW analysis.                                            #
#                                                                              #
# Source:                                                                      #
# - Air pollution data on PM2.5:                                               #
#     PM25_data-wrangling.R                                                    #
# - Data on monitoring stations:                                               #
#     Gobierno de la Ciudad de México                                          #   
#     https://datos.cdmx.gob.mx/dataset/red-automatica-de-monitoreo-atmosferico#
# - Vector data files for geospatial analysis:                                 #
#     Instituto Nacional de Estadística y Geografía (INEGI)                    #
#     https://en.www.inegi.org.mx/app/mapas/                                   #
#                                                                              #
# Authors:                                                                     #
#   Valeria Gracia Olvera, <vgracia@stanford.edu>.                             #
#                                                                              #
# Created: 06/01/24                                                            #
# Last update: 06/10/24                                                        #
#------------------------------------------------------------------------------#
rm(list = ls())  # clean environment

# Load libraries ----------------------------------------------------------
library(sf)
library(raster)
library(ggplot2)
library(tmap)
library(dplyr)
library(cartogram)
library(readxl)
library(spatstat)
library(sparr)
library(gstat)

# Set up ------------------------------------------------------------------
v_years <- 2014:2019

# RMSE function
RMSE_resid <- function(x){
  return(sqrt(mean(x^2)))
}

# Load monitoring stations data
df_monitoring_stations_raw <- read_xlsx("data_raw/monitoring_stations_MCMA_byMUN.xlsx")

# Municipios data
df_mun <- mxmaps::df_mxmunicipio %>% 
  filter(state_code == "09") %>% 
  mutate(MUN_NAME = toupper(municipio_name)) %>% 
  dplyr::select(municipio_code, MUN_NAME) %>% 
  rename(CVE_MUN = municipio_code) %>% 
  mutate(MUN = paste(CVE_MUN, "-", MUN_NAME))

# Analysis ----------------------------------------------------------------
# Shapefiles
poligonos_cdmx <- st_read(dsn = "data_raw/shapefiles/09mun.shp", layer="09mun")
poligonos_cdmx_transf <- st_transform(poligonos_cdmx, 
                                      crs = CRS("+proj=longlat +datum=WGS84"))

poligonos_cdmx_transf_names <- left_join(poligonos_cdmx_transf,df_mun, by = "CVE_MUN")

poligonos_edomx <- st_read(dsn = "data_raw/shapefiles/15mun.shp", layer="15mun")
poligonos_edomx_transf <- st_transform(poligonos_edomx, 
                                       crs = CRS("+proj=longlat +datum=WGS84"))

poligonos_morel <- read_sf(dsn = "data_raw/shapefiles/17mun.shp", layer="17mun")
poligonos_morel_transf <- st_transform(poligonos_morel, 
                                       crs = CRS("+proj=longlat +datum=WGS84"))

tm_shape(poligonos_cdmx_transf_names) + tm_borders() + 
  tm_fill(col = "MUN", alpha = 0.2, title = "Municipality") +
  tm_text("CVE_MUN", size=0.6, just = "center", col = "black") +
  tm_layout(inner.margins = c(0.12,0.4,0.12,0.12),
            legend.width = 0.65) +
  tm_compass()

tmap_save(tmap_last(),"figs/fig_map_CDMX_by_mun.jpg", width = 7, height = 6) 

# Air pollution data 
# PM2.5
load("data/df_conc_PM25.rda")

# Join data with monitoring stations
df_conc_PM25 <- left_join(df_conc_PM25, df_monitoring_stations_raw, by = "station")

# Filtering monitoring stations
v_stations <- unique(df_conc_PM25$station)

df_monitoring_stations <- df_monitoring_stations_raw %>% 
  filter(station %in% v_stations)

# Map monitoring stations
sf_monitoring_stations <- st_as_sf(df_monitoring_stations,
                                   coords = c("x","y"),
                                   crs = CRS("+proj=longlat +datum=WGS84")
)

tm_shape(poligonos_cdmx_transf) + tm_borders() + 
  tm_fill(col = "CVE_MUN", alpha = 0.2, legend.show = F) +
  tm_shape(poligonos_edomx_transf) + tm_borders() +
  tm_fill(col = "white", alpha = 0.2) +
  tm_shape(poligonos_morel_transf) + tm_borders() +
  tm_fill(col = "white", alpha = 0.2) +
  tm_shape(sf_monitoring_stations) +
  tm_dots(size = 0.2) +
  tm_text("station", size=0.6, ymod = -0.5, just = "center") +
  tm_layout(inner.margins = rep(0.12,4)) + 
  # tm_add_legend(type = "fill", 
  #               labels = c("Mexico City","State of Mexico","Morelia"),
  #               col = c("red", "forestgreen", "steelblue"),
  #               alpha = 0.5) +
  tm_compass()

tmap_save(tmap_last(),"figs/fig_map_stations.jpg", width = 7, height = 7) 

# For every week perform a spatial smoothing using IDW  
v_pollutants <- unique(df_conc_PM25$pollutant)
v_cve_mun <- unique(poligonos_cdmx_transf$CVE_MUN)

df_idw_values <- data.frame(NULL)
for(pollutant_i in v_pollutants){ # pollutant_i = "PM25"
  
  print(pollutant_i)
  
  df_temp <- df_conc_PM25 %>% 
    filter(pollutant == pollutant_i)
  
  max_value <- max(df_temp$mean_week)
  min_value <- min(df_temp$mean_week)
  
  v_weeks <- as.character(unique(df_temp$week))
  
  for(week_i in v_weeks){ # week_i = "2019-12-16"
    
    print(week_i)
    
    # Filter
    df_temp_data <- df_temp %>% 
      filter(week == week_i)
    
    if(dim(df_temp_data)[1] < 4){
      next
    }
    
    # Convert to spatial object
    coordinates(df_temp_data) <- ~x+y
    
    # IDW
    study_area <- poligonos_cdmx_transf
    extend_study_area <- extent(study_area)
    
    grid_study_area <- expand.grid(x = seq(from = extend_study_area@xmin,
                                           to = extend_study_area@xmax,
                                           length.out = 500),
                                   y = seq(from = extend_study_area@ymin,
                                           to = extend_study_area@ymax,
                                           length.out = 500))
    
    coordinates(grid_study_area) <- ~ x + y
    gridded(grid_study_area) <- TRUE
    
    # plot(grid_study_area,
    #      main = paste("Air quality monitoring Station\n and Interpolation Grid"),
    #      col = "grey",
    #      cex.main = 0.9)
    # 
    # plot(study_area, add = TRUE, border = "red")
    # plot(df_temp_data, add = TRUE)
    
    # Choose the best model
    v_rmse <- c()
    for(idp_i in 1:6){ # idp_i = 1
      res_idw_temp <- gstat(formula = mean_week ~ 1, # intercept only model
                            data = df_temp_data, 
                            nmax = length(df_temp_data), 
                            set = list(idp = idp_i))
      
      # Cross validation
      res_cv <- gstat.cv(res_idw_temp,
                         nfold = nrow(df_temp_data), # Set the number of folds to the number of rows
                         set = list(idp = idp_i),
                         verbose = F,
      )
      
      v_rmse <- c(v_rmse,RMSE_resid(res_cv@data$residual))
    }
    
    # IDW
    idp_op <- which(v_rmse == min(v_rmse))
    res_idw <- gstat(formula = mean_week ~ 1, # intercept only model
                     data = df_temp_data, 
                     nmax = length(df_temp_data), 
                     set = list(idp = idp_op))
    
    predict_res_idw <- predict(object = res_idw,
                               newdata = grid_study_area)
    
    # Masking
    raster_breaks <- seq(floor(min_value), 
                         ceiling(max_value),
                         length.out = 10)
    
    for(cve_mun_i in v_cve_mun){ # cve_mun_i = "011"
      map_raster <- mask(raster(predict_res_idw), 
                         study_area %>% filter(CVE_MUN == cve_mun_i))
      
      df_idw_values <- rbind(df_idw_values,
                             data.frame(pollutant = pollutant_i,
                                        week = week_i,
                                        CVE_mun = cve_mun_i,
                                        value = mean(map_raster@data@values, na.rm = T),
                                        idp = idp_op))
    }
    
    # # Map
    # map_raster <- mask(raster(predict_res_idw), poligonos_cdmx_transf)
    # 
    # tm_plot <- tm_shape(map_raster) + tm_raster(colorNA = NULL,
    #                                             title = paste0(pollutant_i,"\n",
    #                                                            week_i),
    #                                             palette = "YlOrRd") +
    #   tm_shape(poligonos_cdmx_transf) + tm_borders() + tm_text("CVE_MUN") +
    #   tm_shape(df_temp_data) + tm_dots(size = 0.1, shape = 4, col = "black")
    # 
    # tmap_save(tm = tm_plot,
    #           filename = paste0("figs/IDW/fig_map_",pollutant_i,"_",week_i,".jpg"),
    #           width = 7, height = 7)
  }
  
}

# Save data
save(df_idw_values, file = "output/df_idw_values_week.rda")
