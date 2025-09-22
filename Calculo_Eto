######################### CALCULO EVAPOTRANSPIRACION
##############################################################################
# Autor: Marco Vinicio Mogro
# Fecha: 2025-09-22
# Descripción: Calculo de la evapotranspiracion utilizando el metodo de Thornthwaite
# para la cuenca del Paute (Ecuador) utilizando datos de temperatura promedio.
##############################################################################

# Configurar codificación para evitar problemas con caracteres especiales
options(encoding = "UTF-8")
Sys.setlocale("LC_ALL", "C")

# Cargar librerías
library(data.table)
library(dplyr)
library(Evapotranspiration)
library(SPEI)
library(terra)
library(ggplot2)
library(gridExtra)
library(tibble)
library(tidyr)
library(lubridate)
library(scales)

cat("Librerias cargadas exitosamente\n")

# Datos de entrada
Temperatura <- fread("C:/9.CLASES/Database_Complete/Paute/Ensamble_Paute(2000-2015)_Temp_Promedio.csv")

# Asegurarse que la columna de fecha esté en formato correcto
Temperatura$Fecha <- as.Date(Temperatura$Fecha, format="%Y-%m-%d")

# Datos zona de estudio
Shape <- terra::vect("C:/9.CLASES/INTERPOLACION/Paute/Subcuencas_CRP_50k_UTM_SAM56.shp")

# Pasar a sistema de coordenadas 4326
Shape <- terra::project(Shape, "EPSG:4326")
plot(Shape)

# Extraer el centroide del shapefile, para conocer la latitud que se
# incluirá en el cálculo de evapotranspiración
centroid <- terra::centroids(Shape, TRUE)

# Extraer latitud del centroide (diferencia entre ymin y max de las coordenadas)
coords <- terra::crds(centroid)
latitud <- mean(coords[, "y"])  # Promedio de todas las latitudes de los centroides
print(paste("Latitud promedio es:", latitud))

# Opcional: Cargar observaciones
# Obs <- terra::vect("C:/9.CLASES/INTERPOLACION/Paute/Observaciones.shp")

# Visualizar datos
plot(Shape)
print(Shape$Subcuenca)

# Preparar datos para el cálculo de evapotranspiración
# Crear objeto climático para el método de Thornthwaite
clim <- data.frame(
  Year = as.numeric(format(Temperatura$Fecha, "%Y")),
  Month = as.numeric(format(Temperatura$Fecha, "%m")),
  Tavg = Temperatura$TempAvg
)

# Calcular evapotranspiración usando el método de Thornthwaite
ET <- thornthwaite(
  clim$Tavg,
  lat = latitud,  # Ajustar según la latitud de la cuenca
  na.rm = TRUE,
  verbose = TRUE
)

# Agregar resultados de evapotranspiración al data.table original
Temperatura$ET <- ET
summary(Temperatura$ET)

# Visualizar resultados básicos
plot(Temperatura$Fecha, Temperatura$ET, 
     type = "l", 
     xlab = "Fecha", 
     ylab = "Evapotranspiracion (mm/Month)", 
     main = "Evapotranspiracion - Metodo Thornthwaite")

######################## VISUALIZACION DE RESULTADOS ########################

# Agregar variables adicionales para análisis
Temperatura$Anio <- year(Temperatura$Fecha)
Temperatura$Month <- month(Temperatura$Fecha)
Temperatura$MonthNombre <- month.name[Temperatura$Month]

# Verificar que las variables se crearon correctamente
cat("Rango de Years:", min(Temperatura$Anio, na.rm=TRUE), "-", max(Temperatura$Anio, na.rm=TRUE), "\n")

# =============================================================================
# VISUALIZACIONES AVANZADAS
# =============================================================================

# 1. Gráfico principal de serie temporal con tendencia
p1 <- ggplot(Temperatura, aes(x = Fecha, y = ET)) +
  geom_line(color = "steelblue", linewidth = 0.8, alpha = 0.8) +
  geom_smooth(method = "loess", se = TRUE, color = "red", linewidth = 1.2, alpha = 0.3) +
  labs(
    title = "Serie Temporal de Evapotranspiracion - Metodo Thornthwaite",
    subtitle = paste("Cuenca del Paute (1981-2015) - Latitud:", round(latitud, 3), "grados"),
    x = "Date",
    y = "Evapotranspiracion (mm/Month)",
    caption = "Fuente: Datos climaticos Paute | Metodo: Thornthwaite"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  ) +
  scale_x_date(date_labels = "%Y", date_breaks = "2 years")

# 2. Gráfico de promedios anuales
et_anual <- Temperatura %>%
  group_by(Anio) %>%
  summarise(
    ET_promedio = mean(ET, na.rm = TRUE),
    ET_max = max(ET, na.rm = TRUE),
    ET_min = min(ET, na.rm = TRUE),
    .groups = 'drop'
  )

p3 <- ggplot(et_anual, aes(x = Anio, y = ET_promedio)) +
  geom_line(linewidth = 1.2, color = "darkgreen") +
  geom_point(size = 3, color = "darkgreen") +
  geom_ribbon(aes(ymin = ET_min, ymax = ET_max), alpha = 0.2, fill = "darkgreen") +
  labs(
    title = "Evapotranspiracion Promedio Anual",
    subtitle = "Linea: Promedio | Banda: Rango min-max",
    x = "Year",
    y = "Evapotranspiracion (mm/Month)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5)
  ) +
  scale_x_continuous(breaks = seq(1981, 2015, 2))

# 3. Gráfico de correlación ET vs Temperatura
p4 <- ggplot(Temperatura, aes(x = TempAvg, y = ET)) +
  geom_point(alpha = 0.6, color = "coral") +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1) +
  labs(
    title = "Relacion ET vs Temperatura",
    x = "Temperatura Promedio (C)",
    y = "Evapotranspiracion (mm/Month)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
  ) +
  annotate("text", x = Inf, y = -Inf,
           label = paste("R2 =", round(summary(lm(ET ~ TempAvg, data = Temperatura))$r.squared, 3)),
           hjust = 1.1, vjust = -0.5, size = 4, color = "black")

# 4. Heatmap mensual por Year
et_long <- Temperatura %>%
  group_by(Anio, Month) %>%
  summarise(ET = mean(ET, na.rm = TRUE), .groups = 'drop') %>%
  complete(Anio, Month, fill = list(ET = NA))

p5 <- ggplot(et_long, aes(x = factor(Month), y = factor(Anio), fill = ET)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = mean(Temperatura$ET, na.rm = TRUE),
                       name = "ET\n(mm/Month)") +
  labs(
    title = "Heatmap de Evapotranspiracion Mensual",
    x = "Month",
    y = "Year"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 0),
    legend.position = "right"
  ) +
  scale_x_discrete(labels = month.abb)

# Mostrar gráficos
cat("=== GRAFICO PRINCIPAL ===\n")
suppressWarnings(print(p1))

cat("=== ANALISIS COMPLEMENTARIOS ===\n")
# Combinar gráficos complementarios
suppressWarnings({
  grid.arrange(p3, p4, ncol = 1, nrow = 2)
})
# Mostrar heatmap mensual 
print(p5)

# Estadísticas descriptivas
cat("\n=== ESTADISTICAS DESCRIPTIVAS ===\n")
cat("Evapotranspiracion (mm/Month):\n")
print(summary(Temperatura$ET))

cat("\nEvapotranspiracion promedio por Year:\n")
print(et_anual)

# Guardar resultados
# Descomenta las siguientes líneas si quieres guardar los archivos
# write.csv(Temperatura, "C:/9.CLASES/Database_Complete/Paute/Resultados_ET_Paute_Completo.csv", row.naMonth = FALSE)
# ggsave("ET_Serie_Temporal.png", p1, width = 12, height = 8, dpi = 300)
# ggsave("ET_Analisis_Complementario.png", arrangeGrob(p3, p4, p5, ncol = 1), width = 12, height = 16, dpi = 300)


