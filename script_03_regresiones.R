# ============================================================
# TFM - MUMADE (UCLM)
# Introducción del VE en regiones europeas (NUTS-2)
# Script 03: Modelos de regresión
#            Componentes PCA → penetración del VE
# ============================================================

# ── 0. PAQUETES ──────────────────────────────────────────────

library(dplyr)
library(ggplot2)
library(patchwork)
library(lmtest)     # bptest (Breusch-Pagan), coeftest
library(sandwich)   # vcovHC (errores robustos)
library(car)        # vif, ncvTest, durbinWatsonTest
library(broom)      # tidy, glance → tablas limpias de resultados
library(kableExtra) # formateo de tablas para exportar


# ── 1. CARGA DEL DATASET ─────────────────────────────────────
df <- readRDS("df_regresion.rds")

cat("Observaciones disponibles:", nrow(df), "\n")
cat("Variables:\n")
print(names(df))

# Crear versión log de la dependiente principal
# ev_share es una proporción [0,1] con fuerte asimetría derecha
df <- df |>
  mutate(
    ln_ev_share  = log(ev_share),
    ln_ev_p1000h = log(ev_p1000h + 0.001),  # +0.001 por si hay ceros
    ln_ev_p1000v = log(ev_p1000v + 0.001)
  )

# Estadísticos descriptivos de las dependientes
cat("\n--- Asimetría de las variables dependientes ---\n")
cat("ev_share     — skewness:", round(psych::skew(df$ev_share), 3), "\n")
cat("ln_ev_share  — skewness:", round(psych::skew(df$ln_ev_share), 3), "\n")
cat("(Si |skewness| baja con log, usar versión transformada)\n\n")


# ── 2. ESPECIFICACIÓN DE LOS MODELOS ─────────────────────────
# Se estiman tres modelos por bloque + modelo conjunto:
#
#   M1: solo bloque socioeconómico
#   M2: solo bloque infraestructura
#   M3: solo política pública
#   M4: modelo conjunto (todos los bloques)
#
# Variable dependiente principal: ln_ev_share
# Variables de robustez: ln_ev_p1000h, ln_ev_p1000v

dep_principal  <- "ln_ev_share"
dep_robustez   <- c("ln_ev_p1000h", "ln_ev_p1000v")

# Detectar nombres de componentes automáticamente
comp_socio <- names(df)[grepl("^SC_", names(df))]
comp_infra <- names(df)[grepl("^IN_", names(df))]


# ── 3. FUNCIÓN AUXILIAR: ESTIMAR + DIAGNOSTICAR ──────────────
estimar_modelo <- function(formula_str, data, nombre) {
  
  f   <- as.formula(formula_str)
  mod <- lm(f, data = data)
  
  cat("\n", strrep("─", 55), "\n")
  cat(" ", nombre, "\n")
  cat(strrep("─", 55), "\n")
  cat("Fórmula:", formula_str, "\n\n")
  
  # Resumen estándar
  print(summary(mod))
  
  # ── Diagnósticos ─────────────────────────────────────────
  cat("\n--- Diagnósticos ---\n")
  
  # 1. VIF (multicolinealidad) — solo si hay > 1 regresor
  tryCatch({
    v <- vif(mod)
    cat("VIF:\n")
    print(round(v, 2))
    if (any(v > 5))
      cat("⚠ Algún VIF > 5: posible multicolinealidad moderada\n")
    if (any(v > 10))
      cat("✗ Algún VIF > 10: multicolinealidad severa\n")
    else
      cat("✓ VIF < 5 en todos los regresores\n")
  }, error = function(e) cat("VIF no calculable (modelo simple)\n"))
  
  # 2. Breusch-Pagan (heterocedasticidad)
  bp <- bptest(mod)
  cat("\nBreusch-Pagan: chi2 =", round(bp$statistic, 3),
      "| p =", round(bp$p.value, 4), "\n")
  if (bp$p.value < 0.05) {
    cat("⚠ Heterocedasticidad detectada → se usarán errores robustos HC3\n")
    het <- TRUE
  } else {
    cat("✓ No se rechaza homocedasticidad\n")
    het <- FALSE
  }
  
  # 3. Normalidad de residuos (Shapiro-Wilk, muestra ≤ 5000)
  sw <- shapiro.test(residuals(mod))
  cat("\nShapiro-Wilk: W =", round(sw$statistic, 4),
      "| p =", round(sw$p.value, 4), "\n")
  if (sw$p.value < 0.05)
    cat("⚠ Residuos no normales (puede aceptarse con n > 100 por TLC)\n")
  else
    cat("✓ Residuos compatibles con normalidad\n")
  
  # 4. Coeficientes con errores robustos si hay heterocedasticidad
  if (het) {
    cat("\n--- Coeficientes con errores robustos HC3 ---\n")
    print(coeftest(mod, vcov = vcovHC(mod, type = "HC3")))
  }
  
  return(mod)
}


# ── 4. MODELOS POR BLOQUE ────────────────────────────────────

# ── M1: Bloque socioeconómico ─────────────────────────────────
# SC1 + tert_share (excluida del PCA por comunalidad < 0.40)
f_M1 <- paste(dep_principal, "~",
              paste(c(comp_socio, "tert_share"), collapse = " + "))
M1 <- estimar_modelo(f_M1, df, "M1 — Bloque socioeconómico")

# ── M2: Bloque infraestructura ────────────────────────────────
f_M2 <- paste(dep_principal, "~",
              paste(comp_infra, collapse = " + "))
M2 <- estimar_modelo(f_M2, df, "M2 — Bloque infraestructura")

# ── M3: Bloque política pública ───────────────────────────────
f_M3 <- paste(dep_principal, "~ pol_1yr_gdp + own_tax_gdp")
M3 <- estimar_modelo(f_M3, df, "M3 — Política pública")

# ── M4: Modelo conjunto ───────────────────────────────────────
f_M4 <- paste(dep_principal, "~",
              paste(c(comp_socio, "tert_share",
                      comp_infra,
                      "pol_1yr_gdp", "own_tax_gdp"),
                    collapse = " + "))
M4 <- estimar_modelo(f_M4, df, "M4 — Modelo conjunto")


# ── 5. TABLA COMPARATIVA DE MODELOS ──────────────────────────
cat("\n", strrep("═", 60), "\n")
cat("  TABLA COMPARATIVA DE MODELOS\n")
cat(strrep("═", 60), "\n\n")

extraer_stats <- function(mod, nombre) {
  s <- summary(mod)
  data.frame(
    Modelo   = nombre,
    R2       = round(s$r.squared, 3),
    R2_adj   = round(s$adj.r.squared, 3),
    F_stat   = round(s$fstatistic[1], 2),
    p_F      = round(pf(s$fstatistic[1], s$fstatistic[2],
                        s$fstatistic[3], lower.tail = FALSE), 4),
    AIC      = round(AIC(mod), 1),
    BIC      = round(BIC(mod), 1),
    n        = nobs(mod)
  )
}

tabla_comp <- bind_rows(
  extraer_stats(M1, "M1 Socioeconómico"),
  extraer_stats(M2, "M2 Infraestructura"),
  extraer_stats(M3, "M3 Política"),
  extraer_stats(M4, "M4 Conjunto")
)
print(tabla_comp)


# ── 6. ANÁLISIS DE ROBUSTEZ ───────────────────────────────────
# Reestimar M4 con las otras dos dependientes
cat("\n", strrep("═", 60), "\n")
cat("  ANÁLISIS DE ROBUSTEZ — M4 con dependientes alternativas\n")
cat(strrep("═", 60), "\n\n")

vars_M4 <- paste(c(comp_socio, "tert_share",
                   comp_infra,
                   "pol_1yr_gdp", "own_tax_gdp"),
                 collapse = " + ")

M4b <- estimar_modelo(
  paste("ln_ev_p1000h ~", vars_M4), df,
  "M4b — Robustez: ln(EV por 1000 hab.)"
)

M4c <- estimar_modelo(
  paste("ln_ev_p1000v ~", vars_M4), df,
  "M4c — Robustez: ln(EV por 1000 veh.)"
)

# Comparativa de robustez: ¿cambia el signo o significatividad de los coefs?
cat("\n--- Signos y significatividad comparados (M4 vs M4b vs M4c) ---\n")
tidy_M4  <- tidy(M4)  |> select(term, estimate, p.value) |>
  rename(est_M4 = estimate, p_M4 = p.value)
tidy_M4b <- tidy(M4b) |> select(term, estimate, p.value) |>
  rename(est_M4b = estimate, p_M4b = p.value)
tidy_M4c <- tidy(M4c) |> select(term, estimate, p.value) |>
  rename(est_M4c = estimate, p_M4c = p.value)

robustez <- tidy_M4 |>
  left_join(tidy_M4b, by = "term") |>
  left_join(tidy_M4c, by = "term") |>
  mutate(
    sig_M4  = ifelse(p_M4  < 0.05, "✓", ""),
    sig_M4b = ifelse(p_M4b < 0.05, "✓", ""),
    sig_M4c = ifelse(p_M4c < 0.05, "✓", ""),
    signo_consistente = sign(est_M4) == sign(est_M4b) &
      sign(est_M4) == sign(est_M4c)
  )
print(robustez |> select(term, est_M4, sig_M4,
                         est_M4b, sig_M4b,
                         est_M4c, sig_M4c,
                         signo_consistente))


# ── 7. VARIABLES SIGNIFICATIVAS PARA DBSCAN ──────────────────
# Criterio: significativas (p < 0.05) en M4 Y signo consistente
#           en los tres modelos de robustez
cat("\n", strrep("═", 60), "\n")
cat("  VARIABLES SIGNIFICATIVAS PARA DBSCAN\n")
cat(strrep("═", 60), "\n\n")

vars_sig <- robustez |>
  filter(p_M4 < 0.05, signo_consistente, term != "(Intercept)") |>
  pull(term)

cat("Variables que entrarán en el clustering DBSCAN:\n")
print(vars_sig)

# Dataset para DBSCAN (estandarizado)
df_dbscan <- df |>
  select(region, country, region_name, quality_flag,
         ev_share, all_of(vars_sig)) |>
  na.omit()

# Estandarizar las variables explicativas (DBSCAN es sensible a escala)
df_dbscan_scaled <- df_dbscan |>
  mutate(across(all_of(vars_sig), scale))

saveRDS(df_dbscan_scaled, "df_dbscan.rds")
cat("\n✓ Dataset para DBSCAN guardado como 'df_dbscan.rds'\n")
cat("Dimensiones:", dim(df_dbscan_scaled), "\n")


# ── 8. GRÁFICOS DE DIAGNÓSTICO ────────────────────────────────
# Gráficos clásicos de residuos para M4
cat("\n--- Generando gráficos de diagnóstico de M4 ---\n")
png("diagnostico_M4.png", width = 1400, height = 1200, res = 150)
par(mfrow = c(2, 2))
plot(M4, which = 1:4,
     main = c("Residuos vs Ajustados — M4",
              "Q-Q Normal — M4",
              "Scale-Location — M4",
              "Cook's Distance — M4"))
dev.off()
cat("✓ Figura guardada: diagnostico_M4.png\n")

# Gráfico de coeficientes del modelo conjunto (con IC al 95%)
coef_plot <- tidy(M4, conf.int = TRUE) |>
  filter(term != "(Intercept)") |>
  mutate(
    significativo = p.value < 0.05,
    term = factor(term, levels = rev(term))
  )

p_coef <- ggplot(coef_plot,
                 aes(x = estimate, y = term,
                     color = significativo)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "gray50", linewidth = 0.5) +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.25, linewidth = 0.7) +
  geom_point(size = 3) +
  scale_color_manual(values = c("FALSE" = "gray60",
                                "TRUE"  = "#2166AC"),
                     labels = c("No significativo (p ≥ 0.05)",
                                "Significativo (p < 0.05)"),
                     name   = NULL) +
  labs(title    = "Coeficientes del modelo conjunto (M4)",
       subtitle = paste0("Variable dependiente: ln(EV share) | n = ",
                         nobs(M4), " regiones NUTS-2"),
       x = "Coeficiente (IC 95%)",
       y = NULL,
       caption = "TFM MUMADE · UCLM 2023") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(color = "gray40"))

ggsave("coeficientes_M4.png", p_coef,
       width = 9, height = 6, dpi = 200)
cat("✓ Figura guardada: coeficientes_M4.png\n")

cat("\n✓ Script 03 completado.\n")
cat("Próximo paso: Script 04 — Clustering DBSCAN\n")
cat("Input necesario: df_dbscan.rds\n")
