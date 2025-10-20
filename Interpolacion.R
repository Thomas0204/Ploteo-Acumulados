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
library(jpeg)
library(grid)# Para agregar logo

# Cargar datos
df <- read_csv("datos/precipitacion_acumulado_2025-10-05_11-38.csv")
df_2 <- read_csv("datos/station_coordinates.csv")
df_2 <- df_2[df_2$original_id %in% df$station_code, ]
departamentos_sf <- st_read("SHP/departamentos.shp")
departamentos_sf <- departamentos_sf %>%
  filter(PROVINCIA == "CORDOBA")

# Coordenadas de ciudades
ciudades <- data.frame(
  nombre = c("Córdoba", "Villa María", "Río Cuarto", "Mina Clavero", 
             "V.Mackenna", "Laboulaye", "Cruz del Eje", "Villa de María", 
             "S.Francisco"),
  longitude = c(-64.1914, -63.2438, -64.3497, -65.0044, -64.3907, 
                -63.3885, -64.8076, -63.7234, -62.0839),
  latitude = c(-31.4193, -32.4465, -33.1230, -31.7296, -33.9183, 
               -34.1276, -30.7209, -29.9043, -31.4256)
)

# Convertir ciudades a sf
ciudades_sf <- st_as_sf(ciudades, coords = c("longitude", "latitude"), crs = 4326)

# Definir CRS
utm_crs <- 32720
departamentos_utm <- st_transform(departamentos_sf, utm_crs)
ciudades_utm <- st_transform(ciudades_sf, utm_crs)

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
#grid <- expand.grid(x = x_seq, y = y_seq)
#coordinates(grid) <- ~x + y
#proj4string(grid) <- CRS(paste0("+init=epsg:", utm_crs))

# Convertir a sf y filtrar puntos dentro de Córdoba
#grid_sf <- st_as_sf(grid)
#grid_points <- st_intersection(grid_sf, st_union(departamentos_utm))

# --- GRILLA REGULAR DE 2 km SOLO DENTRO DE CÓRDOBA (UTM)
grid_points <- st_make_grid(departamentos_utm, cellsize = grid_size, what = "centers") %>%
  st_as_sf() %>%
  st_intersection(st_union(departamentos_utm))
st_crs(grid_points) <- st_crs(departamentos_utm)

# Para gstat:
grid_sp <- as(grid_points, "Spatial")
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

# Rasterizar el kriging a la malla y enmascarar por Córdoba
r <- rast(
  extent = ext(as_Spatial(st_union(departamentos_utm))),
  resolution = grid_size,                    # p.ej. 2000 m
  crs = paste0("EPSG:", utm_crs)
  
)
krig_spatvect <- vect(kriging_sf_clean)      # sf -> SpatVector

r_var <- rasterize(krig_spatvect, r, field = "var1.pred", fun = mean)

r_var_mask <- mask(r_var, vect(departamentos_utm))

# Dejar valores <0.1 como NA para pintarlos en blanco
r_var_mask[r_var_mask < 0.1] <- NA

r_stars <- st_as_stars(r_var_mask)

# Extraer coordenadas de las ciudades para usar con geom_text_repel
ciudades_coords <- ciudades_utm %>%
  mutate(
    x = st_coordinates(.)[,1],
    y = st_coordinates(.)[,2]
  ) %>%
  st_drop_geometry()

mapa_kriging_log <- ggplot() +
  geom_stars(data = r_stars) +
  scale_fill_stepsn(
    name    = "\n(mm)",
    colours = c(
      "#c7e3ff", "#417eba", "#3e4887", "#253494",
      "#fdd0a2", "#fcae91", "#fb6a4a"
    ),
    breaks  = c(0.1,5,15,30,50,80,120), # arranca en 0.1
    limits  = c(0,120),
    labels  = scales::label_number(accuracy = 1),
    na.value = "white",
    oob = scales::squish
  ) +
  geom_sf(data = departamentos_utm, fill = NA, color = "black", linewidth = 0.4) +
  geom_sf(data = ciudades_utm, color = "black", size = 2, shape = 19) +
  geom_text_repel(
    data = ciudades_coords,
    aes(x = x, y = y, label = nombre),
    size = 3, color = "black",
    bg.color = "white", bg.r = 0.1,
    force = 2, max.overlaps = Inf, seed = 123
  ) +
  coord_sf(crs = st_crs(utm_crs)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title    = element_text(hjust = 0.5, size = 13, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    plot.caption  = element_text(hjust = 0.5, size = 8, colour = "grey40"), # caption pequeño
    legend.position = "right",
    legend.title    = element_text(size = 12, face = "bold"),
    legend.text     = element_text(size = 11),
    plot.margin     = margin(20, 20, 20, 20)
  ) +
  labs(
    title    = "Precipitación Acumulada",
    subtitle = paste0(first_reading_str, " a ", last_reading_str),
    caption  = "Fuente: Estaciones meteorológicas pertenecientes a la provincia de Córdoba 
    Elaborado a partir del método de interpolación de Kriging"   # 👈 caption
  ) +
  guides(
    fill = guide_legend(            # leyenda discreta con “cajitas” iguales
      keyheight = unit(30, "pt"),
      keywidth  = unit(16, "pt")
    )
  )

print(mapa_kriging_log)

mapa_con_watermark <- mapa_kriging_log +
  annotate(
    "text",
    x = mean(st_bbox(departamentos_utm)[c("xmin","xmax")]),  # centro en X
    y = mean(st_bbox(departamentos_utm)[c("ymin","ymax")]),  # centro en Y
    label = "OHMC",              # 🔹 tu marca
    color = "grey70",            # color
    alpha = 0.3,                 # transparencia
    angle = 45,                  # inclinación
    size = 30,                   # tamaño (ajustá a gusto)
    fontface = "bold"
  )

print(mapa_con_watermark)


ggsave(
  "Interpolacion_10_05.png",
  plot = mapa_con_watermark,
  width = 16, height = 9, units = "in", dpi = 600, bg = "white"
)
