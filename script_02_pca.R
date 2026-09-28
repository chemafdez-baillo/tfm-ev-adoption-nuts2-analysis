# ============================================================
# TFM - MUMADE (UCLM)
# Introducción del VE en regiones europeas (NUTS-2)
# Script 02: Análisis de Componentes Principales (PCA)
#            Bloques socioeconómico e infraestructura
# ============================================================

# ── 0. PAQUETES ──────────────────────────────────────────────

library(readxl)
library(dplyr)
library(psych)         # KMO, Bartlett, fa.parallel
library(factoextra)    # fviz_eig, fviz_pca_var, fviz_pca_biplot
library(ggplot2)
library(patchwork)
library(corrplot)


# ── 1. IMPORTACIÓN Y PREPARACIÓN ─────────────────────────────
df_raw <- read_excel("final_ev_dataset_2023_with_policy.xlsx", sheet = "final_clean_dataset")

df <- df_raw |>
  transmute(
    region       = region,
    country      = country,
    region_name  = region_name,
    quality_flag = quality_flag,
    
    # Dependientes (para correlacionar con componentes al final)
    ev_share     = ev_share_vehicle_stock,
    ev_p1000h    = ev_per_1000_people,
    ev_p1000v    = ev_per_1000_vehicles,
    
    # Bloque socioeconómico
    ln_gdp_pc    = log(gdp_per_capita_eur),
    ln_density   = log(population_density_km2),
    tert_share   = tertiarisation_share,
    emp_pc       = employment_per_capita,
    readiness    = ev_adoption_readiness_index,
    
    # Bloque infraestructura (SIN ev_per_charging_point — duplicado inverso)
    cp_km2       = charging_points_per_km2,
    cp_p1000h    = charging_points_per_1000_people,
    cp_p1000v    = charging_points_per_1000_vehicles,
    cp_ev        = charging_points_per_ev
    # ev_per_charging_point EXCLUIDA: r = -0.98 con cp_ev → redundancia matemática
  )


# ── 2. BLOQUE SOCIOECONÓMICO ──────────────────────────────────
cat("\n", strrep("═", 60), "\n")
cat("  BLOQUE SOCIOECONÓMICO\n")
cat(strrep("═", 60), "\n\n")

socio_vars <- c("ln_gdp_pc", "ln_density", "tert_share", "emp_pc", "readiness")
bloque_socio <- df |> select(all_of(socio_vars)) |> na.omit()

# ── 2a. Retest KMO y Bartlett con variables log-transformadas ─
cat("--- Tests de adecuación (variables finales) ---\n")
cor_socio <- cor(bloque_socio, method = "spearman")
kmo_socio <- KMO(cor_socio)
print(kmo_socio)
bart_socio <- cortest.bartlett(cor_socio, n = nrow(bloque_socio))
cat("Bartlett chi2 =", round(bart_socio$chisq, 2),
    " | p =", format(bart_socio$p.value, scientific = TRUE), "\n\n")

# ── 2b. Número óptimo de componentes: regla Kaiser + scree ───
cat("--- Criterio de selección de componentes ---\n")
pca_socio_all <- prcomp(bloque_socio, scale. = TRUE)
cat("Eigenvalues:\n")
eigenvalues_s <- pca_socio_all$sdev^2
print(round(eigenvalues_s, 3))
cat("\nVarianza explicada acumulada:\n")
print(round(cumsum(eigenvalues_s / sum(eigenvalues_s)) * 100, 1))

# Análisis paralelo (criterio más riguroso que Kaiser)
cat("\n--- Análisis paralelo (Horn) ---\n")
parallel_socio <- fa.parallel(bloque_socio, fm = "minres", fa = "pc",
                              n.iter = 200, show.legend = FALSE,
                              main = "Análisis paralelo — Bloque socioeconómico")
# ── 2c. PCA final con componentes retenidos ───────────────────
# Decisión: retener los componentes con eigenvalue > 1 (Kaiser) confirmados por análisis paralelo
n_comp_socio <- sum(eigenvalues_s > 1)
cat("\nComponentes retenidos (Kaiser, eigenvalue > 1):", n_comp_socio, "\n\n")

pca_socio <- principal(bloque_socio, nfactors = n_comp_socio,
                       rotate = "varimax", scores = TRUE)

cat("--- Cargas factoriales (varimax) ---\n")
print(pca_socio$loadings, cutoff = 0.30)

cat("\n--- Varianza explicada por componente ---\n")
print(pca_socio$Vaccounted)

cat("\n--- Comunalidades ---\n")
print(round(pca_socio$communality, 3))

# ── 2d. Diagnóstico de tert_share ────────────────────────────
cat("\n--- Diagnóstico de tert_share ---\n")
cat("Comunalidad de tert_share:", round(pca_socio$communality["tert_share"], 3), "\n")
cat("Carga máxima de tert_share:", round(max(abs(pca_socio$loadings["tert_share",])), 3), "\n")
cat("\nDecisión: si comunalidad < 0.40 o carga máxima < 0.50,\n")
cat("tert_share se tratará como variable independiente en la regresión.\n\n")

# ── 2e. Puntuaciones de las componentes socioeconómicas ───────
scores_socio <- as.data.frame(pca_socio$scores)
colnames(scores_socio) <- paste0("SC_", seq_len(n_comp_socio))

# Añadir identificadores
scores_socio <- bind_cols(
  df |> select(region, country, region_name, quality_flag,
               ev_share, ev_p1000h, ev_p1000v) |> na.omit(),
  scores_socio
)

cat("Primeras filas de las puntuaciones socioeconómicas:\n")
print(head(scores_socio))


# ── 3. BLOQUE INFRAESTRUCTURA ─────────────────────────────────
cat("\n", strrep("═", 60), "\n")
cat("  BLOQUE INFRAESTRUCTURA\n")
cat(strrep("═", 60), "\n\n")

infra_vars <- c("cp_km2", "cp_p1000h", "cp_p1000v", "cp_ev")
# Nota: ev_per_charging_point ya fue excluida en la importación

bloque_infra <- df |> select(all_of(infra_vars)) |> na.omit()

# ── 3a. Retest KMO y Bartlett sin ev_cp ──────────────────────
cat("--- Tests de adecuación (4 variables, sin EV/CP) ---\n")
cor_infra <- cor(bloque_infra, method = "spearman")
kmo_infra <- KMO(cor_infra)
print(kmo_infra)
bart_infra <- cortest.bartlett(cor_infra, n = nrow(bloque_infra))
cat("Bartlett chi2 =", round(bart_infra$chisq, 2),
    " | p =", format(bart_infra$p.value, scientific = TRUE), "\n\n")

# ── 3b. Número óptimo de componentes ─────────────────────────
cat("--- Criterio de selección de componentes ---\n")
pca_infra_all <- prcomp(bloque_infra, scale. = TRUE)
eigenvalues_i <- pca_infra_all$sdev^2
cat("Eigenvalues:\n")
print(round(eigenvalues_i, 3))
cat("\nVarianza explicada acumulada:\n")
print(round(cumsum(eigenvalues_i / sum(eigenvalues_i)) * 100, 1))

cat("\n--- Criterio de selección: Kaiser + varianza acumulada ---\n")
cat("(Análisis paralelo descartado: caso Heywood por alta colinealidad)\n")
cat("Eigenvalues > 1:", sum(eigenvalues_i > 1), "\n")
cat("Varianza explicada por CP1:",
    round(eigenvalues_i[1] / sum(eigenvalues_i) * 100, 1), "%\n")
cat("Si CP1 > 70% y eigenvalue[2] < 1 → retener 1 componente.\n")

# ── 3c. PCA final infraestructura ────────────────────────────
n_comp_infra <- sum(eigenvalues_i > 1)
cat("\nComponentes retenidos (Kaiser):", n_comp_infra, "\n\n")

pca_infra <- principal(bloque_infra, nfactors = n_comp_infra,
                       rotate = "varimax", scores = TRUE)

cat("--- Cargas factoriales (varimax) ---\n")
print(pca_infra$loadings, cutoff = 0.30)

cat("\n--- Varianza explicada por componente ---\n")
print(pca_infra$Vaccounted)

cat("\n--- Comunalidades ---\n")
print(round(pca_infra$communality, 3))

# ── 3d. Puntuaciones infraestructura ─────────────────────────
scores_infra <- as.data.frame(pca_infra$scores)
colnames(scores_infra) <- paste0("IN_", seq_len(n_comp_infra))

scores_infra <- bind_cols(
  df |> select(region, country, region_name, quality_flag,
               ev_share, ev_p1000h, ev_p1000v) |> na.omit(),
  scores_infra
)

cat("Primeras filas de las puntuaciones de infraestructura:\n")
print(head(scores_infra))


# ── 4. VARIABLES DE POLÍTICA PÚBLICA (directas, sin PCA) ─────
cat("\n", strrep("═", 60), "\n")
cat("  POLÍTICA PÚBLICA — Variables directas\n")
cat(strrep("═", 60), "\n\n")

# Según el diagnóstico del Script 01, KMO = 0.46 → no PCA
# Se usan directamente:
#   pol_1yr_gdp  = ev_policy_total_first_year_over_gdppc  (esfuerzo total)
#   own_tax_gdp  = annual_ownership_tax_saving_over_gdppc (componente autónoma)

df_policy <- df_raw |>
  transmute(
    region      = region,
    country     = country,
    region_name = region_name,
    pol_1yr_gdp = ev_policy_total_first_year_over_gdppc,
    own_tax_gdp = annual_ownership_tax_saving_over_gdppc
  )

cat("Estadísticos de las variables de política:\n")
print(summary(df_policy[, c("pol_1yr_gdp", "own_tax_gdp")]))
cat("\nNúmero de valores distintos (confirma variación between-country):\n")
cat("pol_1yr_gdp:", n_distinct(df_policy$pol_1yr_gdp), "valores únicos\n")
cat("own_tax_gdp:", n_distinct(df_policy$own_tax_gdp), "valores únicos\n")


# ── 5. DATASET CONSOLIDADO PARA REGRESIÓN ────────────────────
cat("\n", strrep("═", 60), "\n")
cat("  DATASET CONSOLIDADO\n")
cat(strrep("═", 60), "\n\n")

# Unir componentes socioeconómicas, de infraestructura,
# variables de política y variables dependientes
df_regresion <- df |>
  select(region, country, region_name, quality_flag,
         ev_share, ev_p1000h, ev_p1000v,
         tert_share) |>       # tert_share: pendiente de decisión
  left_join(
    scores_socio |> select(region, starts_with("SC_")),
    by = "region"
  ) |>
  left_join(
    scores_infra |> select(region, starts_with("IN_")),
    by = "region"
  ) |>
  left_join(
    df_policy |> select(region, pol_1yr_gdp, own_tax_gdp),
    by = "region"
  ) |>
  na.omit()

cat("Dimensiones del dataset para regresión:", dim(df_regresion), "\n")
cat("Variables disponibles:\n")
print(names(df_regresion))

# Guardar para el Script de regresiones
saveRDS(df_regresion, "df_regresion.rds")
cat("\n✓ Dataset guardado como 'df_regresion.rds'\n")


# ── 6. VISUALIZACIONES ───────────────────────────────────────

# ── 6a. Scree plots comparativos ─────────────────────────────
p_scree_s <- fviz_eig(pca_socio_all,
                      addlabels = TRUE, ylim = c(0, 65),
                      barfill = "#5B8DB8", barcolor = "#5B8DB8",
                      linecolor = "#C0392B") +
  geom_hline(yintercept = 20, linetype = "dashed",
             color = "gray40", linewidth = 0.5) +
  labs(title = "Bloque socioeconómico",
       x = "Componente", y = "% varianza explicada") +
  theme_minimal(base_size = 11)

p_scree_i <- fviz_eig(pca_infra_all,
                      addlabels = TRUE, ylim = c(0, 90),
                      barfill = "#E8974A", barcolor = "#E8974A",
                      linecolor = "#C0392B") +
  geom_hline(yintercept = 25, linetype = "dashed",
             color = "gray40", linewidth = 0.5) +
  labs(title = "Bloque infraestructura",
       x = "Componente", y = "% varianza explicada") +
  theme_minimal(base_size = 11)

(p_scree_s | p_scree_i) +
  plot_annotation(
    title    = "Scree plots — Varianza explicada por componente",
    subtitle = "TFM MUMADE · Adopción del VE en regiones NUTS-2 europeas (2023)",
    theme    = theme(plot.title    = element_text(size = 13, face = "bold"),
                     plot.subtitle = element_text(size = 9, color = "gray40"))
  )
ggsave("pca_scree_plots.png", width = 12, height = 5, dpi = 200)

# ── 6b. Círculo de correlaciones (biplots de cargas) ─────────
p_circ_s <- fviz_pca_var(pca_socio_all,
                         axes = c(1, 2),
                         col.var = "contrib",
                         gradient.cols = c("#2166AC", "#F7F7F7", "#D6604D"),
                         repel = TRUE,
                         title = "Bloque socioeconómico — CP1 vs CP2") +
  theme_minimal(base_size = 11)

p_circ_i <- fviz_pca_var(pca_infra_all,
                         axes = c(1, 2),
                         col.var = "contrib",
                         gradient.cols = c("#1A9641", "#F7F7F7", "#D73027"),
                         repel = TRUE,
                         title = "Bloque infraestructura — CP1 vs CP2") +
  theme_minimal(base_size = 11)

(p_circ_s | p_circ_i) +
  plot_annotation(
    title    = "Círculos de correlación de variables — PCA",
    subtitle = "TFM MUMADE · Adopción del VE en regiones NUTS-2 europeas (2023)",
    theme    = theme(plot.title    = element_text(size = 13, face = "bold"),
                     plot.subtitle = element_text(size = 9, color = "gray40"))
  )
ggsave("pca_circulos_correlacion.png", width = 12, height = 5, dpi = 200)

# ── 6c. Correlación de componentes con variables dependientes ─
cat("\n--- Correlación de componentes con variables dependientes ---\n")

# Socioeconómico
comp_socio_names <- paste0("SC_", seq_len(n_comp_socio))
cor_dep_socio <- cor(
  df_regresion[, comp_socio_names],
  df_regresion[, c("ev_share", "ev_p1000h", "ev_p1000v")],
  method = "spearman",
  use = "complete.obs"
)
cat("\nCorrelación componentes socioeconómicas ~ dependientes (Spearman):\n")
print(round(cor_dep_socio, 3))

# Infraestructura
comp_infra_names <- paste0("IN_", seq_len(n_comp_infra))
cor_dep_infra <- cor(
  df_regresion[, comp_infra_names],
  df_regresion[, c("ev_share", "ev_p1000h", "ev_p1000v")],
  method = "spearman",
  use = "complete.obs"
)
cat("\nCorrelación componentes infraestructura ~ dependientes (Spearman):\n")
print(round(cor_dep_infra, 3))

cat("\n✓ Script 02 completado.\n")
cat("Figuras guardadas: pca_scree_plots.png | pca_circulos_correlacion.png\n")
cat("Dataset para regresión: df_regresion.rds\n")
