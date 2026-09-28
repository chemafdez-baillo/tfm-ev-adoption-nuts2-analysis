# ============================================================
# TFM - MUMADE (UCLM)
# Introducción del VE en regiones europeas (NUTS-2)
# Script 04b: Clustering revisado
#             DBSCAN + K-means (sin Grecia)
# ============================================================
# Motivación del script:
#   El Script 04 reveló que las 10 regiones griegas (EL) forman
#   outliers multivariantes por valores extremos en IN_2 y pol_1yr,
#   resultado de imputación al 100% a nivel país. Se excluyen del
#   clustering y se documentan como limitación metodológica.
#   Se complementa DBSCAN con K-means como análisis de robustez.
# ============================================================

# ── 0. PAQUETES ──────────────────────────────────────────────

library(dplyr)
library(ggplot2)
library(patchwork)
library(dbscan)
library(factoextra)
library(cluster)
library(tidyr)
library(scales)


# ── 1. CARGA Y FILTRADO ──────────────────────────────────────
df_full <- readRDS("df_dbscan.rds")

vars_cluster <- c("SC_1", "tert_share", "IN_1", "IN_2",
                  "pol_1yr_gdp", "own_tax_gdp")

# Excluir Grecia: imputación 100% a nivel país → outlier estadístico
df <- df_full |> filter(country != "EL")

cat("Regiones totales:      ", nrow(df_full), "\n")
cat("Regiones tras filtro:  ", nrow(df), "\n")
cat("Regiones griegas excl.:", nrow(df_full) - nrow(df), "\n\n")

X <- as.matrix(df[, vars_cluster])


# ══════════════════════════════════════════════════════════════
#   BLOQUE A: DBSCAN
# ══════════════════════════════════════════════════════════════

# ── A1. Gráfico k-distancias ──────────────────────────────────
minPts_val <- 4
kd         <- kNNdist(X, k = minPts_val - 1)
kd_sorted  <- sort(as.vector(kd), decreasing = TRUE)

df_codo <- data.frame(indice    = seq_along(kd_sorted),
                      distancia = kd_sorted)

p_codo <- ggplot(df_codo, aes(x = indice, y = distancia)) +
  geom_line(color = "#2166AC", linewidth = 0.8) +
  geom_hline(yintercept = 1.5, linetype = "dashed",
             color = "#D6604D", linewidth = 0.6) +
  geom_hline(yintercept = 2.0, linetype = "dashed",
             color = "#F4A582", linewidth = 0.6) +
  annotate("text", x = 15, y = 1.55,
           label = "ε = 1.5", color = "#D6604D", size = 3.5) +
  annotate("text", x = 15, y = 2.05,
           label = "ε = 2.0", color = "#F4A582", size = 3.5) +
  labs(title    = "Gráfico k-distancias — muestra sin Grecia",
       subtitle = paste0("k = ", minPts_val - 1,
                         " | minPts = ", minPts_val,
                         " | n = ", nrow(df), " regiones"),
       x = "Puntos ordenados por distancia (desc.)",
       y = paste0("Distancia al ", minPts_val - 1,
                  "º vecino más cercano")) +
  theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(color = "gray40"))

ggsave("dbscan_codo_sinEL.png", p_codo,
       width = 9, height = 5, dpi = 200)
cat("✓ Gráfico k-distancias guardado: dbscan_codo_sinEL.png\n")

# ── A2. Sensibilidad ε ────────────────────────────────────────
cat("\n--- Sensibilidad ε (minPts = ", minPts_val,
    ", sin Grecia) ---\n")
cat(sprintf("%-6s %-10s %-10s %-12s\n",
            "eps", "clusters", "ruido", "%ruido"))

eps_vals <- c(1.0, 1.5, 2.0, 2.5, 3.0)
for (eps in eps_vals) {
  db  <- dbscan(X, eps = eps, minPts = minPts_val)
  n_c <- max(db$cluster)
  n_r <- sum(db$cluster == 0)
  pct <- round(n_r / nrow(X) * 100, 1)
  cat(sprintf("%-6.1f %-10d %-10d %-12.1f\n", eps, n_c, n_r, pct))
}

# ── A3. DBSCAN definitivo ─────────────────────────────────────
# Ajusta eps_final según la tabla de sensibilidad
eps_final <- 1.5

db_final <- dbscan(X, eps = eps_final, minPts = minPts_val)

cat("\n--- DBSCAN definitivo (eps =", eps_final,
    "| minPts =", minPts_val, ", sin Grecia) ---\n")
cat("Distribución por clúster (0 = ruido):\n")
print(table(db_final$cluster))

df <- df |> mutate(cluster_db = factor(db_final$cluster))

# ── A4. Silueta DBSCAN ────────────────────────────────────────
df_asig <- df |> filter(cluster_db != 0)
X_asig  <- as.matrix(df_asig[, vars_cluster])
cl_asig <- as.integer(as.character(df_asig$cluster_db))

if (n_distinct(cl_asig) >= 2) {
  sil_db <- silhouette(cl_asig, dist(X_asig))
  cat("\nÍndice de silueta DBSCAN:",
      round(mean(sil_db[, 3]), 3), "\n")
} else {
  cat("\n⚠ DBSCAN produce 1 clúster con estos parámetros.\n")
  cat("Revisa eps_final o considera usar solo K-means.\n")
}


# ══════════════════════════════════════════════════════════════
#   BLOQUE B: K-MEANS
# ══════════════════════════════════════════════════════════════

# ── B1. Número óptimo de clústeres: método del codo + silueta ─
cat("\n", strrep("═", 55), "\n")
cat("  BLOQUE B: K-MEANS\n")
cat(strrep("═", 55), "\n\n")

set.seed(123)   # reproducibilidad

# Método del codo (within-cluster sum of squares)
p_codo_km <- fviz_nbclust(X, kmeans, method = "wss",
                          k.max = 8, nstart = 25) +
  geom_vline(xintercept = 3, linetype = "dashed",
             color = "#D6604D") +
  geom_vline(xintercept = 4, linetype = "dashed",
             color = "#F4A582") +
  labs(title    = "K-means — Método del codo (WSS)",
       subtitle = "Muestra sin Grecia | n = 226 regiones") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

# Método de la silueta media
p_sil_km <- fviz_nbclust(X, kmeans, method = "silhouette",
                         k.max = 8, nstart = 25) +
  labs(title    = "K-means — Silueta media por k",
       subtitle = "Muestra sin Grecia | n = 226 regiones") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

(p_codo_km | p_sil_km) +
  plot_annotation(
    title    = "Selección del número de clústeres — K-means",
    subtitle = "TFM MUMADE · Adopción del VE en regiones NUTS-2 europeas (2023)",
    theme    = theme(plot.title    = element_text(size = 13, face = "bold"),
                     plot.subtitle = element_text(size = 9, color = "gray40"))
  )
ggsave("kmeans_seleccion_k.png", width = 12, height = 5, dpi = 200)
cat("✓ Gráfico selección k guardado: kmeans_seleccion_k.png\n")

# ── B2. K-means definitivo ────────────────────────────────────
# Ajusta k_final según los gráficos anteriores
k_final <- 4   # <── modifica si los gráficos sugieren otro valor

set.seed(123)
km_final <- kmeans(X, centers = k_final,
                   nstart = 50, iter.max = 100)

cat("\n--- K-means definitivo (k =", k_final, ") ---\n")
cat("Distribución por clúster:\n")
print(table(km_final$cluster))

# Silueta K-means
sil_km <- silhouette(km_final$cluster, dist(X))
cat("\nÍndice de silueta K-means:",
    round(mean(sil_km[, 3]), 3), "\n")

df <- df |> mutate(cluster_km = factor(km_final$cluster))

# ── B3. Comparación DBSCAN vs K-means ────────────────────────
cat("\n--- Comparación de asignaciones DBSCAN vs K-means ---\n")
print(table(DBSCAN = df$cluster_db, Kmeans = df$cluster_km))


# ══════════════════════════════════════════════════════════════
#   BLOQUE C: CARACTERIZACIÓN (sobre K-means como principal)
# ══════════════════════════════════════════════════════════════
cat("\n", strrep("═", 55), "\n")
cat("  BLOQUE C: CARACTERIZACIÓN DE CLÚSTERES\n")
cat(strrep("═", 55), "\n\n")

perfil_km <- df |>
  group_by(cluster_km) |>
  summarise(
    n          = n(),
    ev_share   = round(mean(ev_share, na.rm = TRUE), 4),
    SC_1       = round(mean(SC_1, na.rm = TRUE), 3),
    tert_share = round(mean(tert_share, na.rm = TRUE), 3),
    IN_1       = round(mean(IN_1, na.rm = TRUE), 3),
    IN_2       = round(mean(IN_2, na.rm = TRUE), 3),
    pol_1yr    = round(mean(pol_1yr_gdp, na.rm = TRUE), 4),
    own_tax    = round(mean(own_tax_gdp, na.rm = TRUE), 4),
    .groups    = "drop"
  ) |>
  arrange(desc(ev_share))

cat("Perfil medio por clúster K-means (ordenado por ev_share):\n")
print(perfil_km)

# Diagnóstico pol_1yr_gdp (encargo Miguel Ángel)
cat("\n--- pol_1yr_gdp por clúster K-means ---\n")
cat("(Examinar si alto pol_1yr coincide con bajo ev_share → endogeneidad)\n\n")
pol_km <- df |>
  group_by(cluster_km) |>
  summarise(
    n            = n(),
    pol_media    = round(mean(pol_1yr_gdp, na.rm = TRUE), 4),
    pol_mediana  = round(median(pol_1yr_gdp, na.rm = TRUE), 4),
    ev_share_med = round(mean(ev_share, na.rm = TRUE), 4),
    .groups      = "drop"
  ) |>
  arrange(desc(ev_share_med))
print(pol_km)

# Países dominantes por clúster
cat("\n--- Países más frecuentes por clúster K-means ---\n")
paises_km <- df |>
  count(cluster_km, country) |>
  group_by(cluster_km) |>
  slice_max(n, n = 6) |>
  arrange(cluster_km, desc(n))
print(paises_km, n = 60)


# ══════════════════════════════════════════════════════════════
#   BLOQUE D: VISUALIZACIONES
# ══════════════════════════════════════════════════════════════

colores_k <- c("1" = "#2166AC", "2" = "#D6604D", "3" = "#1A9850",
               "4" = "#8B4513", "5" = "#7B2D8B")

# ── D1. Proyección PCA 2D ─────────────────────────────────────
pca_vis  <- prcomp(X, scale. = FALSE)
df_pca   <- as.data.frame(pca_vis$x[, 1:2]) |>
  mutate(cluster_km = df$cluster_km,
         cluster_db = df$cluster_db,
         ev_share   = df$ev_share)

var_exp <- summary(pca_vis)$importance[2, 1:2] * 100

p_proj_km <- ggplot(df_pca,
                    aes(x = PC1, y = PC2, color = cluster_km)) +
  geom_point(alpha = 0.7, size = 2) +
  stat_ellipse(aes(group = cluster_km), type = "norm",
               linetype = "dashed", linewidth = 0.5) +
  scale_color_manual(values = colores_k, name = "Clúster") +
  labs(title    = "K-means — Proyección en espacio PCA",
       subtitle = paste0("k = ", k_final,
                         " | n = ", nrow(df), " regiones (sin Grecia)"),
       x = paste0("CP1 (", round(var_exp[1], 1), "%)"),
       y = paste0("CP2 (", round(var_exp[2], 1), "%)")) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(color = "gray40"))

ggsave("kmeans_proyeccion_pca.png", p_proj_km,
       width = 9, height = 6, dpi = 200)
cat("\n✓ Proyección PCA K-means guardada: kmeans_proyeccion_pca.png\n")

# ── D2. Perfil de clústeres ───────────────────────────────────
perfil_long <- perfil_km |>
  select(cluster_km, SC_1, tert_share, IN_1, IN_2,
         pol_1yr, own_tax) |>
  pivot_longer(-cluster_km,
               names_to  = "variable",
               values_to = "valor")

p_perfil <- ggplot(perfil_long,
                   aes(x = variable, y = valor,
                       fill = factor(cluster_km))) +
  geom_col(position = "dodge", alpha = 0.85) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "gray40") +
  scale_fill_manual(values = colores_k, name = "Clúster") +
  labs(title    = "Perfil medio de clústeres K-means",
       subtitle = "Variables estandarizadas | muestra sin Grecia",
       x = NULL, y = "Valor medio estandarizado",
       caption  = "TFM MUMADE · UCLM 2023") +
  theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(color = "gray40"),
        axis.text.x   = element_text(angle = 30, hjust = 1))

ggsave("kmeans_perfil_clusters.png", p_perfil,
       width = 10, height = 6, dpi = 200)
cat("✓ Perfil clústeres guardado: kmeans_perfil_clusters.png\n")

# ── D3. Boxplot ev_share por clúster ─────────────────────────
p_box <- ggplot(df, aes(x = cluster_km, y = ev_share,
                        fill = cluster_km)) +
  geom_boxplot(alpha = 0.75, outlier.size = 1.5) +
  scale_fill_manual(values = colores_k) +
  labs(title    = "Distribución de ev_share por clúster K-means",
       subtitle = "Variable dependiente principal | muestra sin Grecia",
       x = "Clúster", y = "ev_share (proporción)",
       caption  = "TFM MUMADE · UCLM 2023") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none",
        plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(color = "gray40"))

ggsave("kmeans_boxplot_evshare.png", p_box,
       width = 8, height = 5, dpi = 200)
cat("✓ Boxplot ev_share guardado: kmeans_boxplot_evshare.png\n")

# ── D4. Comparación silueta DBSCAN vs K-means ────────────────
p_sil_final <- fviz_silhouette(sil_km, print.summary = FALSE) +
  labs(title    = paste0("Silueta K-means (k = ", k_final, ")"),
       subtitle = paste0("Silueta media = ",
                         round(mean(sil_km[, 3]), 3))) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggsave("kmeans_silueta.png", p_sil_final,
       width = 9, height = 5, dpi = 200)
cat("✓ Silueta K-means guardada: kmeans_silueta.png\n")


# ── EXPORTAR ─────────────────────────────────────────────────
df_export <- df |>
  select(region, country, region_name, quality_flag,
         ev_share, cluster_db, cluster_km,
         all_of(vars_cluster))

write.csv(df_export, "resultados_clusters_v2.csv",
          row.names = FALSE, fileEncoding = "UTF-8")
saveRDS(df_export, "df_resultados_finales.rds")

cat("\n✓ Script 04b completado.\n")
cat("Archivos generados:\n")
cat("  dbscan_codo_sinEL.png\n")
cat("  kmeans_seleccion_k.png\n")
cat("  kmeans_proyeccion_pca.png\n")
cat("  kmeans_perfil_clusters.png\n")
cat("  kmeans_boxplot_evshare.png\n")
cat("  kmeans_silueta.png\n")
cat("  resultados_clusters_v2.csv\n")
cat("  df_resultados_finales.rds\n")


df_full <- readRDS("df_dbscan.rds")
df_full |>
  filter(country == "EL") |>
  select(region, region_name, country) |>
  print(n = 20)

write.csv(
  df_export |>
    select(region, region_name, country, quality_flag) |>
    arrange(country, region),
  "anexo_A.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)