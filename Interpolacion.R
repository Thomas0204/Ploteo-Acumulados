# Limpiar entorno
rm(list = ls())
gc() 
library(sf)
library(dplyr)
library(ggplot2)
library(automap)
library(gstat)
library(terra)
library(stars)
library(ggrepel)
library(readr)
library(sp)  # Añadir sp para funciones coordinates

# Cargar datos
df <- read_csv("precipitacion_acumulado_2025-09-01_21-13.csv")
df_2 <- read_csv("station_coordinates.csv")
df_2 <- df_2[df_2$original_id %in% df$station_code, ]
departamentos_sf <- st_read("SHP/departamentos.shp")
departamentos_sf <- departamentos_sf %>%
  filter(PROVINCIA == "CORDOBA")


# Definir CRS
utm_crs <- 32720
departamentos_utm <- st_transform(departamentos_sf, utm_crs)

# Normalizar IDs
norm_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("\\.0+$", "", x)
  x <- sub("^0+", "", x)
  x
}

df <- df %>% mutate(station_code_norm = norm_id(station_code))
df_2 <- df_2 %>% mutate(original_id_norm = norm_id(original_id))

# Unir datos y crear objeto sf
df_Acumulados <- df %>%
  left_join(df_2, by = c("station_code_norm" = "original_id_norm")) %>%
  filter(!is.na(longitude) & !is.na(latitude) & !is.na(accumulated_rain_mm)) %>%
  st_as_sf(coords = c("longitude","latitude"), crs = 4326)

# Extraer valores únicos
first_reading_val <- unique(df_Acumulados$first_reading)
last_reading_val  <- unique(df_Acumulados$last_reading)

# Como todos son iguales, nos quedamos con el primero
first_reading_val <- first_reading_val[1]
last_reading_val  <- last_reading_val[1]

# Si querés darles formato bonito
first_reading_str <- format(first_reading_val, "%Y-%m-%d %H:%M")
last_reading_str  <- format(last_reading_val, "%Y-%m-%d %H:%M")

# TRANSFORMAR A UTM PARA KRIGING
#--------------------------------
df_Acumulados_utm <- st_transform(df_Acumulados, utm_crs)

# CREAR GRILLA PARA INTERPOLACIÓN
#---------------------------------
# Obtener los límites de Córdoba
bbox_cordoba <- st_bbox(departamentos_utm)

# Crear grilla regular (ajusta la resolución según necesites)
# 1000m = 1km de resolución, puedes cambiar a 500 o 2000
grid_size <- 2000

# Crear secuencias de coordenadas
x_seq <- seq(bbox_cordoba[1], bbox_cordoba[3], by = grid_size)
y_seq <- seq(bbox_cordoba[2], bbox_cordoba[4], by = grid_size)

# Crear grilla

grid <- expand.grid(x = x_seq, y = y_seq)
coordinates(grid) <- ~x + y
proj4string(grid) <- CRS(paste0("+init=epsg:", utm_crs))

# Convertir a sf y filtrar puntos dentro de Córdoba

grid_sf <- st_as_sf(grid)
grid_points <- st_intersection(grid_sf, st_union(departamentos_utm))

# PREPARAR DATOS PARA GSTAT
#---------------------------
# Convertir df_Acumulados a formato sp (requerido por gstat)

df_sp <- as(df_Acumulados_utm, "Spatial")

# ANÁLISIS EXPLORATORIO DE DATOS
#--------------------------------
# Verificar distribución de datos

summary(df_Acumulados$accumulated_rain_mm)
hist(df_Acumulados$accumulated_rain_mm, 
     main = "Distribución de Precipitación Acumulada",
     xlab = "Precipitación (mm)")

# VARIOGRAMA
#------------
# Calcular variograma experimental

v_exp <- variogram(accumulated_rain_mm ~ 1, df_sp)

# Ajustar modelo de variograma automáticamente
v_fit <- autofitVariogram(accumulated_rain_mm ~ 1, df_sp)

# Visualizar variograma
plot(v_fit)

# KRIGING ORDINARIO
#-------------------
# Convertir grid_points a sp
grid_sp <- as(grid_points, "Spatial")

# Realizar kriging
kriging_result <- krige(accumulated_rain_mm ~ 1, 
                        df_sp, 
                        grid_sp, 
                        model = v_fit$var_model)

# Convertir resultado a sf
kriging_sf <- st_as_sf(kriging_result)

# VISUALIZACIÓN

#---------------
# Limpiar datos interpolados (remover NAs)
kriging_sf_clean <- kriging_sf %>%
  filter(!is.na(var1.pred))

library(terra)
library(stars)
library(scales)

# --- (si no lo tenías aún) rasterizar el kriging a la malla y enmascarar por Córdoba
r <- rast(
  extent = ext(as_Spatial(st_union(departamentos_utm))),
  resolution = grid_size,                    # p.ej. 5000 m
  crs = paste0("EPSG:", utm_crs)
)
krig_spatvect <- vect(kriging_sf_clean)      # sf -> SpatVector
r_var <- rasterize(krig_spatvect, r, field = "var1.pred", fun = mean)
r_var_mask <- mask(r_var, vect(departamentos_utm))
r_stars <- st_as_stars(r_var_mask)

mapa_kriging_log <- ggplot() +
  geom_stars(data = r_stars) +
  scale_fill_gradientn(
    name   = "Precipitación\n(mm)",
    colours = c("white", "lightblue", "blue", "darkblue", "purple", "red"),
    values  = scales::rescale(c(0, 10, 50, 100, 150, 180, 240, 300)),
    breaks  = c(0, 10, 25, 50, 75, 100, 150, 180, 240, 300),
    limits  = c(0, 300),
    labels  = scales::label_number(accuracy = 1),
    oob     = scales::squish
  ) +
  geom_sf(data = departamentos_utm, fill = NA, color = "black", linewidth = 0.4) +
  coord_sf(crs = st_crs(utm_crs)) +
  theme_void() +
  labs(
    title    = "Precipitación Acumulada",
    subtitle = paste0("Periodo de observación: ",
                      first_reading_str, " a ", last_reading_str)) +
  theme(
    plot.title    = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "right",
    legend.title    = element_text(size = 12, face = "bold"),
    legend.text     = element_text(size = 11)
  ) +
  guides(
    fill = guide_colorbar(
      barwidth  = 1.5,   # ancho de la barra
      barheight = 20,    # largo de la barra
      ticks     = TRUE,
      frame.colour = "black"
    )
  )

print(mapa_kriging_log)

ggsave(
  "mapa_kriging.png",
  plot = mapa_kriging_log,
  width = 16, height = 9, units = "in", dpi = 600,
  bg = "white"    # <- evita gris en el PNG
)

# MAPA CON CONTORNOS (ALTERNATIVO)
#----------------------------------
# Crear mapa alternativo con contornos suaves
mapa_alternativo <- ggplot() +
  geom_sf(data = departamentos_utm, fill = "white", color = "black", linewidth = 0.5) +
  geom_sf(data = kriging_sf_clean, 
          aes(fill = var1.pred), 
          color = NA) +
  scale_fill_gradient2(name = "Precipitación\n(mm)", 
                       low = "blue", 
                       mid = "yellow", 
                       high = "red",
                       midpoint = median(kriging_sf_clean$var1.pred, na.rm = TRUE),
                       na.value = "transparent") +
  geom_sf(data = df_Acumulados_utm, 
          color = "black", 
          size = 2, 
          alpha = 0.8,
          shape = 21,
          fill = "white",
          stroke = 1) +
  labs(title = "Mapa de Precipitación con Estaciones",
       subtitle = "Interpolación Kriging - Córdoba") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        legend.position = "bottom")

print(mapa_alternativo)

# MAPA CON CONTORNOS
#-------------------
# Crear isolíneas
iso_lines <- st_contour(kriging_sf, contour_lines = TRUE, breaks = 10)

mapa_contornos <- mapa_base +
  geom_sf(data = kriging_sf, aes(fill = var1.pred), color = NA) +
  geom_sf(data = iso_lines, color = "white", alpha = 0.7, size = 0.3) +
  scale_fill_viridis_c(name = "Precipitación\n(mm)", 
                       option = "plasma") +
  geom_sf(data = df_Acumulados_utm, 
          color = "red", 
          size = 1.5, 
          alpha = 0.8)

print(mapa_contornos)