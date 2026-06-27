################################################################################
# R SCRIPT FOR: Article 1 - CAM molecules (VCAM/ICAM/selectins) in HIV
# Analyses:
# 1) Data cleaning & group definition (AIDS vs NO-AIDS)
# 2) PCA + K-means clustering (PCA vars: CAMs + age + CBC) + Confusion matrix
# 3) Logistic regression (AIDS) with univariate screening (p <= 0.20) + forced sex
# 4) Model diagnostics and metrics (no ROC figure)
# 5) Save figures and tables in "artigo1_results/"
################################################################################

#### 0) SETUP ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(FactoMineR)
  library(factoextra)
  library(broom)
  library(patchwork)
  library(writexl)
  library(ResourceSelection) # hoslem.test
  library(car)              # vif
  library(pscl)             # pR2
  library(pROC)             # AUC (no plot)
  library(scales)           # pretty percent labels (optional)
})

cat("=== Article 1 (CAM) analysis started ===\n")

#### 1) OUTPUT FOLDER ----
results_path_artigo1 <- "artigo1_results/"
if (!dir.exists(results_path_artigo1)) dir.create(results_path_artigo1, recursive = TRUE)

#### 2) INPUT DATA ----
stopifnot(exists("banco_CAM"))
cat("Rows in banco_CAM:", nrow(banco_CAM), "\n")

#### 3) DATA CLEANING & GROUPS ----
# age column = i
# sex column = s
# cd4 column = cd4
# viral load log column = log
# CAM molecules: vcam, icam, cd62e, cd62l, cd62p

vars_cam <- c("vcam", "icam", "cd62e", "cd62l", "cd62p")

vars_cbc <- c(
  "wbc", "linf", "neu", "eos", "bas", "mon",
  "rbc", "hgb", "hct", "vcm", "hcm", "chcm", "rdw", "plt"
)

dados_cam <- banco_CAM %>%
  transmute(
    codigo = codigo,
    idade  = suppressWarnings(as.numeric(i)),
    sexo   = as.factor(s),
    cd4    = suppressWarnings(as.numeric(cd4)),
    Log    = suppressWarnings(as.numeric(log)),
    vcam   = suppressWarnings(as.numeric(vcam)),
    icam   = suppressWarnings(as.numeric(icam)),
    cd62e  = suppressWarnings(as.numeric(cd62e)),
    cd62l  = suppressWarnings(as.numeric(cd62l)),
    cd62p  = suppressWarnings(as.numeric(cd62p)),
    across(all_of(vars_cbc), ~ suppressWarnings(as.numeric(.x)))
  ) %>%
  mutate(
    aids_status = factor(if_else(cd4 <= 350, "AIDS", "NO-AIDS"),
                         levels = c("NO-AIDS", "AIDS"))
  ) %>%
  filter(!is.na(aids_status), !is.na(idade), !is.na(sexo))

cat("After basic filtering, N =", nrow(dados_cam), "\n")
cat("Outcome counts:\n")
print(table(dados_cam$aids_status))

################################################################################
#### 4) PCA + CLUSTERING (CAMs + age + CBC) + CONFUSION MATRIX ----
################################################################################
cat("\n--- 4) PCA + clustering (CAMs + age + CBC) ---\n")

# PCA variables: ONLY CAMs + age + CBC (all numeric)
vars_pca <- c(vars_cam, "idade", vars_cbc)

dados_pca <- dados_cam %>%
  select(codigo, aids_status, sexo, all_of(vars_pca)) %>%
  drop_na()

cat("PCA dataset N =", nrow(dados_pca), "\n")

X_pca <- dados_pca %>% select(all_of(vars_pca))

# Run PCA with scaling
pca_res <- PCA(X_pca, graph = FALSE, scale.unit = TRUE)

# ---------- THEME: larger fonts + axis lines in black + no grid ----------
base_theme_clean <- theme_minimal(base_size = 16) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    plot.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 14)
  )

# p1: Individuals colored by AIDS group
p1 <- fviz_pca_ind(
  pca_res,
  geom.ind = "point",
  col.ind = dados_pca$aids_status,
  addEllipses = TRUE,
  legend.title = "Group"
) + labs(
  title = "PCA of CAM markers + age + CBC",
  subtitle = "Individuals colored by AIDS status (CD4 ≤ 350)"
) + base_theme_clean

# p2: Biplot
p2 <- fviz_pca_biplot(
  pca_res,
  repel = TRUE
) + labs(
  title = "PCA biplot",
  subtitle = "Variables and individuals (scaled)"
) + base_theme_clean

# Clustering input (scaled matrix)
X_scaled <- scale(X_pca)

# p3: Elbow
p3 <- fviz_nbclust(X_scaled, kmeans, method = "wss") +
  labs(
    title = "Elbow method (K-means)",
    subtitle = "Within-cluster sum of squares"
  ) + base_theme_clean

# p4: K-means clusters (k=2 default)
set.seed(123)
k_chosen <- 2
km_res <- kmeans(X_scaled, centers = k_chosen, nstart = 25)

p4 <- fviz_cluster(
  km_res,
  data = X_scaled,
  geom = "point",
  ellipse.type = "convex"
) + labs(
  title = paste0("K-means clustering (k = ", k_chosen, ")"),
  subtitle = "Unsupervised clustering on scaled PCA variables"
) + base_theme_clean

################################################################################
#### 4B) PERMANOVA, PERMDISP AND CLUSTER STABILITY CHECKS ----
################################################################################

suppressPackageStartupMessages({
  library(vegan)
  library(cluster)
})

cat("\n--- 4B) PERMANOVA, PERMDISP and silhouette analyses ---\n")

# Distance matrix using the same scaled data used for PCA/K-means
dist_pca <- dist(X_scaled, method = "euclidean")

# PERMANOVA: tests whether AIDS and NO-AIDS differ in multivariate space
set.seed(123)
permanova_aids <- vegan::adonis2(
  dist_pca ~ aids_status,
  data = dados_pca,
  permutations = 999
)

cat("\nPERMANOVA by AIDS status:\n")
print(permanova_aids)

# PERMDISP: tests whether group dispersions differ
permdisp_aids <- vegan::betadisper(
  dist_pca,
  group = dados_pca$aids_status
)

cat("\nPERMDISP by AIDS status:\n")
print(anova(permdisp_aids))

# PERMANOVA for K-means clusters
set.seed(123)
permanova_cluster <- vegan::adonis2(
  dist_pca ~ cluster,
  data = cluster_df,
  permutations = 999
)

cat("\nPERMANOVA by K-means cluster:\n")
print(permanova_cluster)

# PERMDISP for K-means clusters
permdisp_cluster <- vegan::betadisper(
  dist_pca,
  group = cluster_df$cluster
)

cat("\nPERMDISP by K-means cluster:\n")
print(anova(permdisp_cluster))

# Silhouette width for K-means clusters
sil <- cluster::silhouette(km_res$cluster, dist_pca)

cat("\nSilhouette summary:\n")
print(summary(sil))

mean_sil <- mean(sil[, "sil_width"], na.rm = TRUE)

cat("\nMean silhouette width:", round(mean_sil, 3), "\n")

# Export results
sink(file.path(results_path_artigo1, "permanova_permdisp_silhouette_results.txt"))
cat("PERMANOVA by AIDS status:\n")
print(permanova_aids)

cat("\nPERMDISP by AIDS status:\n")
print(anova(permdisp_aids))

cat("\nPERMANOVA by K-means cluster:\n")
print(permanova_cluster)

cat("\nPERMDISP by K-means cluster:\n")
print(anova(permdisp_cluster))

cat("\nSilhouette summary:\n")
print(summary(sil))

cat("\nMean silhouette width:", round(mean_sil, 3), "\n")
sink()

# p5: Confusion matrix (cluster vs AIDS status) with % within cluster (row %)
cluster_df <- tibble(
  aids_status = dados_pca$aids_status,
  cluster = factor(km_res$cluster, levels = sort(unique(km_res$cluster)))
)

tab_counts <- table(cluster_df$cluster, cluster_df$aids_status)
tab_perc <- prop.table(tab_counts, margin = 1) * 100  # within cluster

conf_df <- as.data.frame(tab_counts) %>%
  rename(cluster = Var1, aids_status = Var2, n = Freq) %>%
  mutate(
    perc = as.vector(tab_perc),
    label = sprintf("%d\n(%.1f%%)", n, perc)
  )

p5 <- ggplot(conf_df, aes(x = aids_status, y = cluster, fill = perc)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = label), size = 5, fontface = "bold") +
  scale_fill_continuous(name = "% within cluster") +
  labs(
    title = "Cluster composition by AIDS status",
    subtitle = "Row percentages within each cluster",
    x = "Clinical group",
    y = "K-means cluster"
  ) + base_theme_clean

# Patchwork panel with labels A) B) C) D) E)
# NOTE: tag_levels="A" gives A, B, C...; we style it as "A)" using tag_prefix/tag_suffix
pca_panel <- ((p1 + p2) / (p3 + p4) / p5) +
  plot_annotation(tag_levels = "A", tag_suffix = ")") &
  theme(plot.tag = element_text(size = 18, face = "bold"))

print(pca_panel)

ggsave(
  filename = file.path(results_path_artigo1, "figure_pca_panel_p1_to_p5.png"),
  plot = pca_panel,
  width = 16, height = 16, dpi = 300
)

cat("PCA panel saved.\n\n")

################################################################################
#### 5) LOGISTIC REGRESSION (AIDS) WITH UNIVARIATE SCREENING ----
################################################################################
cat("\n--- 5) Logistic regression with univariate screening (p <= 0.20) ---\n")

# IMPORTANT: Fix the "OR ~ 1 (1-1)" issue:
# - Model CAM effects per 1 SD increase (z-score) => interpretable ORs.
# Also: ensure we keep Log even if it has NA (drop_na later handles it).

dados_model <- dados_cam %>%
  mutate(
    # z-scores for CAMs (only among non-missing values)
    across(all_of(vars_cam), ~ as.numeric(scale(.x)), .names = "z_{.col}"),
    aids_status = factor(aids_status, levels = c("NO-AIDS", "AIDS")),
    sexo = as.factor(sexo)
  )

candidate_predictors_logistic <- c(
  "idade", "sexo", "Log",
  "z_vcam", "z_icam", "z_cd62e", "z_cd62l", "z_cd62p"
)

modelo_log_data <- dados_model %>%
  select(aids_status, all_of(candidate_predictors_logistic)) %>%
  drop_na()

cat("Logistic dataset N =", nrow(modelo_log_data), "\n")
cat("Outcome counts:\n")
print(table(modelo_log_data$aids_status))

# ---- Helper: univariate model returning OR + CI + p(Wald) and p(LRT) ----
get_univ_glm_stats <- function(df, outcome, predictor) {
  f <- as.formula(paste0(outcome, " ~ ", predictor))
  m <- glm(f, data = df, family = binomial)
  
  tt <- broom::tidy(m, exponentiate = TRUE, conf.int = TRUE)
  row_pred <- tt %>% filter(term == predictor)
  
  lrt <- anova(m, test = "Chisq")
  p_lrt <- lrt$`Pr(>Chi)`[2]
  
  tibble(
    predictor = predictor,
    OR = as.numeric(row_pred$estimate),
    CI_low = as.numeric(row_pred$conf.low),
    CI_high = as.numeric(row_pred$conf.high),
    p_wald = as.numeric(row_pred$p.value),
    p_lrt = as.numeric(p_lrt)
  )
}

# Run univariate screening for all candidates EXCEPT sexo (forced)
univ_candidates <- setdiff(candidate_predictors_logistic, "sexo")

univ_screen <- purrr::map_dfr(
  univ_candidates,
  ~ get_univ_glm_stats(modelo_log_data, "aids_status", .x)
) %>%
  arrange(p_lrt)

cat("\nUnivariate screening (sorted by LRT p):\n")
print(univ_screen)

# Select by p_lrt <= 0.20
selected_predictors <- univ_screen %>%
  filter(p_lrt <= 0.20) %>%
  pull(predictor)

# Force sex
selected_predictors <- unique(c("sexo", selected_predictors))

# Safety fallback
if (length(selected_predictors) == 1 && selected_predictors[1] == "sexo") {
  selected_predictors <- c("sexo", "idade", "Log")
  cat("No predictors met p<=0.20; using fallback: sexo + idade + Log\n")
}

cat("Selected predictors (p<=0.20) + forced sexo:\n")
print(selected_predictors)

final_formula <- as.formula(paste("aids_status ~", paste(selected_predictors, collapse = " + ")))
modelo_logistico <- glm(final_formula, data = modelo_log_data, family = binomial)

cat("\n--- Final logistic model summary ---\n")
print(summary(modelo_logistico))

# Final adjusted OR table (drop intercept for reporting)
tabela_or <- broom::tidy(modelo_logistico, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  transmute(
    Variavel = term,
    OR = round(estimate, 3),
    CI_low = round(conf.low, 3),
    CI_high = round(conf.high, 3),
    p_value = signif(p.value, 3)
  )

cat("\n--- Adjusted ORs (final model) ---\n")
print(tabela_or)

################################################################################
#### 6) MODEL DIAGNOSTICS & METRICS ----
################################################################################
cat("\n--- 6) Model diagnostics & metrics ---\n")

anova_lrt <- anova(modelo_logistico, test = "Chisq")
cat("\nSequential LRT (added terms):\n")
print(anova_lrt)

r2 <- pscl::pR2(modelo_logistico)
cat("\nPseudo-R2:\n")
print(r2)

y01 <- as.numeric(modelo_log_data$aids_status) - 1
hl <- ResourceSelection::hoslem.test(y01, fitted(modelo_logistico), g = 10)
cat("\nHosmer-Lemeshow GOF:\n")
print(hl)

vifs <- car::vif(modelo_logistico)
cat("\nVIF:\n")
print(vifs)

roc_obj <- pROC::roc(
  response = modelo_log_data$aids_status,
  predictor = fitted(modelo_logistico),
  levels = c("NO-AIDS", "AIDS"),
  direction = "<"
)
auc_val <- as.numeric(pROC::auc(roc_obj))
cat("\nAUC =", round(auc_val, 3), "\n")

diag_table <- tibble(
  metric = c("N", "AIC", "Null deviance", "Residual deviance",
             "McFadden_R2", "Nagelkerke_r2CU", "AUC",
             "HL_chi_square", "HL_df", "HL_p_value"),
  value = c(
    nrow(modelo_log_data),
    AIC(modelo_logistico),
    modelo_logistico$null.deviance,
    modelo_logistico$deviance,
    unname(r2["McFadden"]),
    unname(r2["r2CU"]),
    auc_val,
    unname(hl$statistic),
    unname(hl$parameter),
    unname(hl$p.value)
  )
)

cat("\nDiagnostics table:\n")
print(diag_table)

################################################################################
#### 7) EXPORT TABLES (Excel) ----
################################################################################
cat("\n--- 7) Exporting tables ---\n")

univ_export <- univ_screen %>%
  mutate(
    OR = round(OR, 3),
    CI_low = round(CI_low, 3),
    CI_high = round(CI_high, 3),
    p_wald = signif(p_wald, 3),
    p_lrt = signif(p_lrt, 3)
  ) %>%
  transmute(
    Variavel = predictor,
    OR = OR,
    CI_low = CI_low,
    CI_high = CI_high,
    p_wald = p_wald,
    p_lrt = p_lrt
  )

vif_export <- tibble(term = names(vifs), VIF = as.numeric(vifs))

# IMPORTANT: Sheet names must be <= 31 chars and cannot include []:*?/\
write_xlsx(
  list(
    "Univariate_OR_CI" = univ_export,
    "Final_model_AOR"  = tabela_or,
    "Diagnostics"      = diag_table,
    "VIF"              = vif_export
  ),
  path = file.path(results_path_artigo1, "tables_logistic_CAM_article1.xlsx")
)

cat("Excel saved to:", file.path(results_path_artigo1, "tables_logistic_CAM_article1.xlsx"), "\n")
cat("\n=== Done! Check folder:", results_path_artigo1, "===\n")
