# Limpiar entorno
rm(list = ls())
gc() 

library(tidyverse)
library(readr)
library(stringr)
library(sf)
library(ggplot2)
library(ggrepel)

df <- read_csv(
  file = "precipitacion_acumulado_2025-09-01_21-13.csv")

df_2 <- read_csv(
  file = "station_coordinates.csv")

departamentos_sf <- st_read("SHP/departamentos.shp")

departamentos_sf <- departamentos_sf %>%
  filter(PROVINCIA == "CORDOBA")


df_2 <- df_2[df_2$original_id %in% df$station_code, ]


# Normalizar a texto para que comparen igual
norm_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)                  # quita espacios al inicio/final
  x <- sub("\\.0+$", "", x)       # quita ".0"
  x <- sub("^0+", "", x)          # quita ceros a la izquierda
  x
}

df <- df %>% mutate(station_code_norm = norm_id(station_code))
df_2 <- df_2 %>% mutate(original_id_norm = norm_id(original_id))

df_Acumulados <- df %>%
  left_join(df_2, by = c("station_code_norm" = "original_id_norm"))

# 3) Convertir tu data frame df (acumulados) a puntos sf
#    Ajustá los nombres de columnas si se llaman distinto a "latitude"/"longitude".

df_Acumulados <- df_Acumulados %>%  # por las dudas
  st_as_sf(coords = c("longitude","latitude"), crs = 4326)

ggplot() +
  geom_sf(data = departamentos_sf, fill = "grey95", color = "grey50", linewidth = 0.2) +
  geom_sf(data = df_Acumulados, size = 2) +
  geom_text_repel(
    data = df_Acumulados %>% st_drop_geometry(), # para etiquetar con ggrepel hay que pasar sin geometría
    aes(x = st_coordinates(df_Acumulados)[,1],
        y = st_coordinates(df_Acumulados)[,2],
        label = station_code),
    size = 3, max.overlaps = 50, min.segment.length = 0
  ) +
  coord_sf() +
  theme_minimal() +
  labs(title = "Ubicacion de estaciones",
       subtitle = "Etiquetas: Codigo de estacion",
       x = NULL, y = NULL)
