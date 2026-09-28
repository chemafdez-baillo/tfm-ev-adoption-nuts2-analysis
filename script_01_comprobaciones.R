
# ============================================================
# TFM - MUMADE (UCLM)
# Introducción del VE en regiones europeas (NUTS-2)
# Script 01: Importación, selección de variables y
#            matrices de correlación intra-bloque
# ============================================================

# ── 0. LIMPIAR ENTORNO Y PAQUETES ──────────────────────────────────────────────
rm(list = ls())
library(readxl)
library(dplyr)
library(corrplot)
library(ggcorrplot)
library(psych)               # KMO y Bartlett
library(PerformanceAnalytics) # chart.Correlation


# ── 1. IMPORTACIÓN ───────────────────────────────────────────
df_raw <- read_excel("final_ev_dataset_2023_with_policy.xlsx", sheet = "final_clean_dataset")

cat("Dimensiones:", nrow(df_raw), "filas x", ncol(df_raw), "columnas\n")
cat("Países:", n_distinct(df_raw$country), "\n")
cat("Regiones:", n_distinct(df_raw$region), "\n")


# ── 2. SELECCIÓN Y TRANSFORMACIONES ──────────────────────────

df <- df_raw |>
  transmute(
    # Identificadores
    region      = region,
    country     = country,
    region_name = region_name,
    nuts2       = nuts2,
    quality_flag = quality_flag,
    
    # ── DEPENDIENTES ──────────────────────────────────────────
    # Principal: cuota de VE sobre parque total de vehículos
    ev_share    = ev_share_vehicle_stock,
    # Alternativas (para análisis de robustez)
    ev_p1000h   = ev_per_1000_people,
    ev_p1000v   = ev_per_1000_vehicles,
    # Nota: bev_share_ev_stock = 1 en todas las obs. → descartada
    
    # ── BLOQUE SOCIOECONÓMICO ─────────────────────────────────
    # Log de GDP per cápita (fuerte asimetría positiva)
    ln_gdp_pc   = log(gdp_per_capita_eur),
    # Log de densidad de población (rango 0.3–5035, muy asimétrica)
    ln_density  = log(population_density_km2),
    # Grado de terciarización de la economía (% GVA servicios)
    tert_share  = tertiarisation_share,
    # Tasa de empleo
    emp_pc      = employment_per_capita,
    # Índice de preparación para la adopción del VE (ya tipificado)
    readiness   = ev_adoption_readiness_index,
    
    # ── BLOQUE INFRAESTRUCTURA ────────────────────────────────
    # Densidad geográfica de puntos de carga
    cp_km2      = charging_points_per_km2,
    # Acceso poblacional a puntos de carga
    cp_p1000h   = charging_points_per_1000_people,
    # Acceso relativo al parque de vehículos
    cp_p1000v   = charging_points_per_1000_vehicles,
    # Disponibilidad de carga por EV existente
    cp_ev       = charging_points_per_ev,
    # Presión sobre la red (inverso conceptual de cp_ev)
    ev_cp       = ev_per_charging_point,
    
    # ── BLOQUE POLÍTICA PÚBLICA ───────────────────────────────
    # Versiones en euros (solo 26 valores distintos, nivel país)
    sub_eur     = bev_purchase_subsidy_eur,
    reg_tax_eur = registration_tax_saving_eur,
    own_tax_eur = annual_ownership_tax_saving_eur,
    pol_1yr_eur = ev_policy_total_first_year_eur,
    # Versiones normalizadas por GDP per cápita (variación regional real)
    sub_gdp     = bev_purchase_subsidy_over_gdppc,
    reg_tax_gdp = registration_tax_saving_over_gdppc,
    own_tax_gdp = annual_ownership_tax_saving_over_gdppc,
    pol_1yr_gdp = ev_policy_total_first_year_over_gdppc
  )

cat("\nObservaciones con datos completos:", sum(complete.cases(df)), "\n")
glimpse(df)


# ── 3. ESTADÍSTICOS DESCRIPTIVOS ─────────────────────────────

vars_analiticas <- df |>
  select(ev_share, ev_p1000h, ev_p1000v,
         ln_gdp_pc, ln_density, tert_share, emp_pc, readiness,
         cp_km2, cp_p1000h, cp_p1000v, cp_ev, ev_cp,
         sub_gdp, reg_tax_gdp, own_tax_gdp, pol_1yr_gdp)

desc <- describe(vars_analiticas)   # psych::describe
print(round(desc[, c("n","mean","sd","median","min","max","skew","kurtosis")], 3))


# ── 4. MATRICES DE CORRELACIÓN INTRA-BLOQUE ──────────────────
# Usamos Spearman: más robusto ante outliers y no exige normalidad.
# Si prefieres Pearson, cambia method = "pearson".

metodo_cor <- "spearman"

# Subconjuntos por bloque (solo las variables analíticas, sin NAs)
bloque_socio  <- df |> select(ln_gdp_pc, ln_density, tert_share, emp_pc, readiness) |> na.omit()
bloque_infra  <- df |> select(cp_km2, cp_p1000h, cp_p1000v, cp_ev, ev_cp)           |> na.omit()
bloque_policy <- df |> select(sub_gdp, reg_tax_gdp, own_tax_gdp, pol_1yr_gdp)       |> na.omit()

cor_socio  <- cor(bloque_socio,  method = metodo_cor)
cor_infra  <- cor(bloque_infra,  method = metodo_cor)
cor_policy <- cor(bloque_policy, method = metodo_cor)

# ── 4a. Imprimir matrices ────────────────────────────────────
cat("\n=== Correlaciones Bloque Socioeconómico ===\n")
print(round(cor_socio, 3))

cat("\n=== Correlaciones Bloque Infraestructura ===\n")
print(round(cor_infra, 3))

cat("\n=== Correlaciones Bloque Política Pública ===\n")
print(round(cor_policy, 3))


# ── 5. TESTS FORMALES: KMO Y BARTLETT ────────────────────────
# KMO > 0.6 y Bartlett p < 0.05 justifican el PCA en cada bloque.

cat("\n=== Test KMO y Bartlett — Bloque Socioeconómico ===\n")
KMO(cor_socio)
cortest.bartlett(cor_socio, n = nrow(bloque_socio))

cat("\n=== Test KMO y Bartlett — Bloque Infraestructura ===\n")
KMO(cor_infra)
cortest.bartlett(cor_infra, n = nrow(bloque_infra))

cat("\n=== Test KMO y Bartlett — Bloque Política Pública ===\n")
KMO(cor_policy)
cortest.bartlett(cor_policy, n = nrow(bloque_policy))


# ── 6. VISUALIZACIÓN ─────────────────────────────────────────

# Paleta común y etiquetas legibles
etiq_socio  <- c("ln(GDP pc)", "ln(Densidad)", "Terciarización", "Empleo pc", "Readiness")
etiq_infra  <- c("CP/km²", "CP/1000 hab.", "CP/1000 veh.", "CP/EV", "EV/CP")
etiq_policy <- c("Subsidio/GDPpc", "Impto.Matr./GDPpc", "Impto.Prop./GDPpc", "Total 1er año/GDPpc")

colnames(cor_socio)  <- rownames(cor_socio)  <- etiq_socio
colnames(cor_infra)  <- rownames(cor_infra)  <- etiq_infra
colnames(cor_policy) <- rownames(cor_policy) <- etiq_policy

# ── corrplot (versión clásica) ───────────────────────────────
par(mfrow = c(1, 3))

corrplot(cor_socio,
         method  = "color", type = "upper", order = "AOE",
         addCoef.col = "black", number.cex = 0.75,
         tl.col  = "black", tl.cex = 0.8, cl.cex = 0.7,
         title   = "Socioeconómico", mar = c(0, 0, 2, 0))

corrplot(cor_infra,
         method  = "color", type = "upper", order = "AOE",
         addCoef.col = "black", number.cex = 0.75,
         tl.col  = "black", tl.cex = 0.8, cl.cex = 0.7,
         title   = "Infraestructura", mar = c(0, 0, 2, 0))

corrplot(cor_policy,
         method  = "color", type = "upper", order = "AOE",
         addCoef.col = "black", number.cex = 0.8,
         tl.col  = "black", tl.cex = 0.8, cl.cex = 0.7,
         title   = "Política pública", mar = c(0, 0, 2, 0))

par(mfrow = c(1, 1))

# ── ggcorrplot (versión ggplot2, más fácil de exportar) ──────
library(ggplot2)
library(patchwork)   # instalar si no está: install.packages("patchwork")

p1 <- ggcorrplot(cor_socio,
                 method = "square", type = "upper",
                 lab = TRUE, lab_size = 3,
                 colors = c("#d73027", "white", "#1a9850"),
                 title = "Bloque socioeconómico",
                 ggtheme = theme_minimal(base_size = 10))

p2 <- ggcorrplot(cor_infra,
                 method = "square", type = "upper",
                 lab = TRUE, lab_size = 3,
                 colors = c("#d73027", "white", "#1a9850"),
                 title = "Bloque infraestructura",
                 ggtheme = theme_minimal(base_size = 10))

p3 <- ggcorrplot(cor_policy,
                 method = "square", type = "upper",
                 lab = TRUE, lab_size = 3,
                 colors = c("#d73027", "white", "#1a9850"),
                 title = "Bloque política pública",
                 ggtheme = theme_minimal(base_size = 10))

# Composición en una sola figura
(p1 | p2 | p3) +
  plot_annotation(
    title    = "Matrices de correlación de Spearman intra-bloque",
    subtitle = "TFM MUMADE · Adopción del VE en regiones NUTS-2 europeas (2023)",
    theme    = theme(plot.title = element_text(size = 13, face = "bold"),
                     plot.subtitle = element_text(size = 9, color = "gray40"))
  )

# Guardar figura
ggsave("correlaciones_intra_bloque.png", width = 14, height = 5, dpi = 200)
cat("\nFigura guardada como 'correlaciones_intra_bloque.png'\n")


# ── 7. REGLA DE DECISIÓN PARA EL PCA ─────────────────────────
# Si algún par de variables supera |r| >= 0.70, el PCA está justificado
# en ese bloque. Imprimimos un resumen automático.

decisión_pca <- function(mat_cor, nombre_bloque) {
  mat <- mat_cor
  diag(mat) <- NA
  max_r <- max(abs(mat), na.rm = TRUE)
  pares_altos <- sum(abs(mat) >= 0.70, na.rm = TRUE) / 2  # triangulo superior
  cat(sprintf("\n[%s] Max |r| = %.3f | Pares con |r| >= 0.70: %d → PCA %s\n",
              nombre_bloque, max_r, pares_altos,
              ifelse(max_r >= 0.70, "RECOMENDADO ✓", "no prioritario")))
}

decisión_pca(cor_socio,  "Socioeconómico")
decisión_pca(cor_infra,  "Infraestructura")
decisión_pca(cor_policy, "Política Pública")

