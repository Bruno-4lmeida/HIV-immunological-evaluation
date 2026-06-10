################################################################################
#
# R SCRIPT FOR:
# IMMUNE SIGNATURES IN HIV INFECTION: sIL-6 AND VIRAL LOAD PREDICT ADVANCED IMMUNODEFICIENCY, 
# WHILE ANTIRETROVIRAL THERAPY RAPIDLY REDUCES sIL-2R LEVELS
#
# Author: Bruno Almeida Silva
# Date: 2026-01-21
#
################################################################################

#### 0. SETUP ####

# ---- Packages ----
# install.packages(c("readxl","tidyverse","ggpubr","FactoMineR","factoextra","corrplot","broom","writexl","patchwork","rstatix","scales"))

library(readxl)
library(tidyverse)
library(ggpubr)
library(FactoMineR)
library(factoextra)
library(corrplot)
library(broom)
library(writexl)
library(patchwork)
library(rstatix)
library(scales)

# ---- Global options ----
set.seed(123)

# ---- Paths (EDIT THESE TWO LINES) ----
data_path_main      <- "C:/Users/bruno/OneDrive/Documentos/BRUNO/DOUTORADO/EXPERIMENTOS_DOC/hiv_data.xlsx"
data_path_treatment <- "C:/Users/bruno/OneDrive/Documentos/BRUNO/DOUTORADO/EXPERIMENTOS_DOC/analise_tratados.xlsx"

# ---- Output folder ----
results_path <- "artigo2_results"
if (!dir.exists(results_path)) dir.create(results_path, recursive = TRUE)

cat("=== Script started. Outputs will be saved in:", results_path, "===\n\n")

base_theme <- theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 14)
  )

################################################################################
#### 1. DATA CLEANING AND PREPARATION ####
################################################################################

cat("--- 1. Data Cleaning and Preparation ---\n")

hiv_data <- read_excel(data_path_main)

dados_limpos <- hiv_data %>%
  filter(!is.na(cd4)) %>%
  mutate(
    cv   = suppressWarnings(as.numeric(as.character(cv))),
    Log  = suppressWarnings(as.numeric(as.character(Log))),
    sexo = as.factor(sexo)
  ) %>%
  filter(idade < 100)

cat("Data cleaned. Kept", nrow(dados_limpos), "rows.\n\n")


################################################################################
#### 2. FEATURE ENGINEERING ####
################################################################################

cat("--- 2. Feature Engineering ---\n")

dados_analise <- dados_limpos %>%
  mutate(
    cd4_group = factor(
      case_when(
        cd4 <= 200 ~ "<= 200",
        cd4 > 200 & cd4 < 350 ~ ">200 & <350",
        cd4 >= 350 ~ ">= 350"
      ),
      levels = c("<= 200", ">200 & <350", ">= 350")
    ),
    IL2_IL2R_ratio = IL2 / il2r,
    IL6_IL6R_ratio = IL6 / il6r
  )

cat("Created cd4_group and cytokine ratios.\n\n")


################################################################################
#### 3. COMPARATIVE ANALYSIS: BIOMARKERS BY CD4 GROUP ####
################################################################################

cat("--- 3. Comparative Plots by CD4 Group ---\n")

my_comparisons <- list(
  c("<= 200", ">200 & <350"),
  c("<= 200", ">= 350"),
  c(">200 & <350", ">= 350")
)


# p1: il2r
p1 <- ggplot(dados_analise, aes(x = cd4_group, y = il2r, fill = cd4_group)) +
  geom_boxplot() +
  geom_jitter(width = 0.12, alpha = 0.5) +
  stat_compare_means(method = "kruskal.test", label.y = 4500) +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    p.adjust.method = "BH",
    label = "p.signif"
  ) +
  scale_y_continuous(breaks = seq(0, 5000, by = 500), limits = c(0, 5000)) +
  labs(
    title = "sIL-2R across CD4 strata",
    x = "CD4+ group (cells/mm³)",
    y = "sIL-2R (pg/mL)"
  ) +
  base_theme

# p2: il6r
wilcox_y_positions <- c(1200, 1500, 1700)
p2 <- ggplot(dados_analise, aes(x = cd4_group, y = il6r, fill = cd4_group)) +
  geom_boxplot() +
  geom_jitter(width = 0.12, alpha = 0.5) +
  stat_compare_means(method = "kruskal.test", label.y = 1750) +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    p.adjust.method = "BH",
    label = "p.signif",
    label.y = wilcox_y_positions
  ) +
  scale_y_continuous(breaks = seq(0, 2000, by = 300), limits = c(0, 2000)) +
  labs(
    title = "sIL-6R across CD4 strata",
    x = "CD4+ group (cells/mm³)",
    y = "sIL-6R (pg/mL)"
  ) +
  base_theme

# p3: cd40l_pg_ml
p3 <- ggplot(dados_analise, aes(x = cd4_group, y = cd40l_pg_ml, fill = cd4_group)) +
  geom_boxplot() +
  geom_jitter(width = 0.12, alpha = 0.5) +
  stat_compare_means(method = "kruskal.test", label.y.npc = 0.9) +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    p.adjust.method = "BH",
    label = "p.signif"
  ) +
  scale_y_continuous(breaks = seq(0, 40000, by = 10000), limits = c(0, 40000)) +
  labs(
    title = "sCD40L across CD4 strata",
    x = "CD4+ group (cells/mm³)",
    y = "sCD40L (pg/mL)"
  ) +
  base_theme

# p4: CD4/CD8
p4 <- ggplot(dados_analise, aes(x = cd4_group, y = `CD4/CD8`, fill = cd4_group)) +
  geom_boxplot() +
  geom_jitter(width = 0.12, alpha = 0.5) +
  stat_compare_means(method = "kruskal.test", label.y = 3) +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    p.adjust.method = "BH",
    label = "p.signif"
  ) +
  labs(
    title = "CD4/CD8 ratio across CD4 strata",
    x = "CD4+ group (cells/mm³)",
    y = "CD4/CD8 ratio"
  ) +
  base_theme

# p5: IL2
p5 <- ggplot(dados_analise, aes(x = cd4_group, y = IL2, fill = cd4_group)) +
  geom_boxplot() +
  geom_jitter(width = 0.12, alpha = 0.5) +
  stat_compare_means(method = "kruskal.test", label.y = 9.5) +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    p.adjust.method = "BH",
    label = "p.signif"
  ) +
  labs(
    title = "sIL-2 across CD4 strata",
    x = "CD4+ group (cells/mm³)",
    y = "sIL-2 (pg/mL)"
  ) +
  base_theme

# p6: IL6
p6 <- ggplot(dados_analise, aes(x = cd4_group, y = IL6, fill = cd4_group)) +
  geom_boxplot() +
  geom_jitter(width = 0.12, alpha = 0.5) +
  stat_compare_means(method = "kruskal.test", label.y = 200) +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    p.adjust.method = "BH",
    label = "p.signif"
  ) +
  labs(
    title = "sIL-6 across CD4 strata",
    x = "CD4+ group (cells/mm³)",
    y = "sIL-6 (pg/mL)"
  ) +
  base_theme

# Combined grid for CD4 group plots
plot_cd4_groups <- patchwork::wrap_plots(
  list(p1, p2, p5, p6, p3, p4),
  ncol = 2
) +
  patchwork::plot_annotation(
    tag_levels = list(c("A)", "B)", "C)", "D)", "E)", "F)"))
  )

print(plot_cd4_groups)

ggsave(
  filename = file.path(results_path, "figure_cd4_groups_p1_to_p6.png"),
  plot = plot_cd4_groups,
  width = 18, height = 16, dpi = 300
)

cat("CD4-group plots saved.\n\n")


################################################################################
#### 4. PCA ####
################################################################################

cat("--- 4. PCA ---\n")

# Build PCA dataframe keeping cd4_group aligned (no rowname hacks)
pca_df <- dados_analise %>%
  select(cd4_group, il2r, il6r, cd4, cd8, `CD4/CD8`, IL2, IL6) %>%
  drop_na()

pca_data <- pca_df %>% select(-cd4_group)

pca_result <- PCA(pca_data, graph = FALSE, scale.unit = TRUE)

plot_pca_individuals <- fviz_pca_ind(
  pca_result,
  geom.ind = "point",
  col.ind = pca_df$cd4_group,
  addEllipses = TRUE,
  legend.title = "CD4 Group"
) + 
  labs(title = "PCA of HIV patients (immune markers)") +
  base_theme

plot_pca_biplot <- fviz_pca_biplot(
  pca_result,
  repel = TRUE
) + labs(title = "PCA biplot (variables + individuals)") + theme_minimal(base_size = 14)

print(plot_pca_individuals)
print(plot_pca_biplot)

ggsave(file.path(results_path, "figure_pca_individuals.png"), plot_pca_individuals, width = 9, height = 7, dpi = 300)
ggsave(file.path(results_path, "figure_pca_biplot.png"), plot_pca_biplot, width = 9, height = 7, dpi = 300)

cat("PCA complete.\n\n")


################################################################################
#### 5. CLUSTERING ####
################################################################################

cat("--- 5. Clustering ---\n")

plot_elbow <- fviz_nbclust(pca_data, kmeans, method = "wss") + labs(subtitle = "Elbow method")
print(plot_elbow)
ggsave(file.path(results_path, "figure_elbow_method.png"), plot_elbow, width = 8, height = 6, dpi = 300)

# K-means (choose k=3 as in your original script; you can change after elbow inspection)
k <- 3
kmeans_result <- kmeans(pca_data, centers = k, nstart = 25)

plot_kmeans <- fviz_cluster(
  kmeans_result,
  data = pca_data,
  geom = "point",
  ellipse.type = "convex",
  ggtheme = theme_bw()
) + labs(title = paste0("K-means clustering (k=", k, ")"))+
  base_theme

print(plot_kmeans)
ggsave(file.path(results_path, "figure_kmeans_clusters.png"), plot_kmeans, width = 9, height = 7, dpi = 300)

# Hierarchical
hclust_result <- hclust(dist(pca_data), method = "ward.D2")

plot_hclust <- fviz_dend(
  hclust_result,
  k = k,
  cex = 0.6,
  rect = TRUE,
  rect_fill = TRUE
) + labs(title = paste0("Hierarchical clustering (Ward.D2), k=", k))

print(plot_hclust)
ggsave(file.path(results_path, "figure_hclust_dendrogram.png"), plot_hclust, width = 10, height = 6, dpi = 300)

# Unified clustering panel
plot_cluster_panel <- ggarrange(
  plot_pca_individuals, plot_pca_biplot, plot_elbow, plot_kmeans,
  ncol = 2, nrow = 2,
  labels = c("A", "B", "C", "D")
)

print(plot_cluster_panel)
ggsave(file.path(results_path, "figure_cluster_panel.png"), plot_cluster_panel, width = 14, height = 12, dpi = 300)

cat("Clustering complete.\n\n")

################################################################################
#### 5B. CLUSTERING (Option B: cluster on PCA scores / PCs) ####
################################################################################

cat("--- 5B. Clustering on PCs (PCA scores) ---\n")

# Choose how many PCs to use for clustering.
# Common choices: 2–10 PCs, or enough PCs to explain ~80–90% variance.
m <- 5  # <-- adjust if needed

# FactoMineR stores individual coordinates (scores) here:
pc_scores <- as.data.frame(pca_result$ind$coord[, 1:m, drop = FALSE])

# Optional: keep group labels aligned (not used in clustering, only for coloring/inspection)
pc_scores$cd4_group <- pca_df$cd4_group

# --- Elbow method on PC space ---
plot_elbow_pc <- fviz_nbclust(pc_scores %>% select(-cd4_group), kmeans, method = "wss") +
  labs(subtitle = paste0("Elbow method (PC1–PC", m, ")"))

print(plot_elbow_pc)
ggsave(file.path(results_path, "figure_elbow_method_PCs.png"),
       plot_elbow_pc, width = 8, height = 6, dpi = 300)

# --- K-means on PC space ---
set.seed(42)  # for reproducibility
k <- 3  # <-- change after elbow inspection if needed
kmeans_pc <- kmeans(pc_scores %>% select(-cd4_group), centers = k, nstart = 25)

plot_kmeans_pc <- fviz_cluster(
  kmeans_pc,
  data = pc_scores %>% select(-cd4_group),
  geom = "point",
  ellipse.type = "convex",
  ggtheme = theme_bw()
) + labs(title = paste0("K-means clustering on PCs (k=", k, ", PC1–PC", m, ")"))

print(plot_kmeans_pc)
ggsave(file.path(results_path, "figure_kmeans_clusters_PCs.png"),
       plot_kmeans_pc, width = 9, height = 7, dpi = 300)

# --- Hierarchical clustering on PC space ---
hclust_pc <- hclust(dist(pc_scores %>% select(-cd4_group)), method = "ward.D2")

plot_hclust_pc <- fviz_dend(
  hclust_pc,
  k = k,
  cex = 0.6,
  rect = TRUE,
  rect_fill = TRUE
) + labs(title = paste0("Hierarchical clustering on PCs (Ward.D2), k=", k, ", PC1–PC", m, ")"))

print(plot_hclust_pc)
ggsave(file.path(results_path, "figure_hclust_dendrogram_PCs.png"),
       plot_hclust_pc, width = 10, height = 6, dpi = 300)

# --- Unified clustering panel (PCA + clustering-on-PCs) ---
plot_cluster_panel_pc <- ggarrange(
  plot_pca_individuals, plot_pca_biplot, plot_elbow_pc, plot_kmeans_pc,
  ncol = 2, nrow = 2,
  labels = c("A", "B", "C", "D")
)

print(plot_cluster_panel_pc)
ggsave(file.path(results_path, "figure_cluster_panel_TESTE.png"),
       plot_cluster_panel_pc, width = 14, height = 12, dpi = 300)

cat("Clustering on PCs complete.\n\n")



################################################################################
#### 6. MULTIVARIATE LINEAR MODEL (CD4) WITH UNIVARIATE SCREENING (p <= 0.20) ####
################################################################################

cat("--- 6. Linear regression (CD4) with univariate screening ---\n")

# Candidate predictors (edit here if you want to add/remove)
candidate_predictors_linear <- c(
  "idade", "sexo", "Log",
  "il2r", "il6r",
  "IL2_IL2R_ratio", "IL6_IL6R_ratio",
  "TNF", "IFN"
)

modelo_linear_data <- dados_analise %>%
  select(cd4, all_of(candidate_predictors_linear)) %>%
  drop_na()

cat("Linear model dataset N =", nrow(modelo_linear_data), "\n")

# --- helper: get univariate p-value using ANOVA for each predictor ---
get_univ_p_lm <- function(df, outcome, predictor) {
  f <- as.formula(paste0(outcome, " ~ ", predictor))
  m <- lm(f, data = df)
  a <- anova(m)
  # predictor is always row 1 in this simple model
  p <- a$`Pr(>F)`[1]
  tibble(predictor = predictor, p_value = p)
}

univ_screen_linear <- map_dfr(candidate_predictors_linear, ~get_univ_p_lm(modelo_linear_data, "cd4", .x)) %>%
  arrange(p_value)

print(univ_screen_linear)

selected_predictors_linear <- univ_screen_linear %>%
  filter(p_value <= 0.20) %>%
  pull(predictor)

cat("Selected predictors (p<=0.20):", paste(selected_predictors_linear, collapse = ", "), "\n")

# If nothing passes, fall back to a minimal model
if (length(selected_predictors_linear) == 0) {
  selected_predictors_linear <- c("idade", "sexo", "Log")
  cat("No predictors met p<=0.20. Using fallback:", paste(selected_predictors_linear, collapse = ", "), "\n")
}

final_formula_linear <- as.formula(paste("cd4 ~", paste(selected_predictors_linear, collapse = " + ")))
modelo_linear <- lm(final_formula_linear, data = modelo_linear_data)

cat("\n--- Linear Model Summary ---\n")
print(summary(modelo_linear))

tabela_modelo_linear <- tidy(modelo_linear, conf.int = TRUE)
print(tabela_modelo_linear)

write_xlsx(
  list("Univariate screening (linear)" = univ_screen_linear,
       "Final linear model (CD4)" = tabela_modelo_linear),
  path = file.path(results_path, "tables_linear_model_cd4.xlsx")
)

cat("Linear model tables saved.\n\n")


################################################################################
#### 7. LOGISTIC MODEL (AIDS: CD4 < 350) WITH UNIVARIATE SCREENING + SEX FORCED ####
################################################################################

cat("--- 7. Logistic regression (AIDS) with univariate screening (sexo forced) ---\n")

candidate_predictors_logistic <- c(
  "idade", "sexo", "Log",
  "IL2", "IL6", "il2r", "il6r", "CD40L")

# --- 7.1 Data prep ---
modelo_log_data <- dados_analise %>%
  mutate(aids_status = factor(if_else(cd4 < 350, "Yes", "No"),
                              levels = c("No", "Yes"))) %>%
  select(aids_status, all_of(candidate_predictors_logistic)) %>%
  drop_na()

cat("Logistic model dataset N =", nrow(modelo_log_data), "\n")
cat("Outcome counts:\n")
print(table(modelo_log_data$aids_status))

# --- 7.2 Univariate logistic regressions ---
run_univ_logistic <- function(df, outcome, predictor) {
  f <- as.formula(paste0(outcome, " ~ ", predictor))
  m <- glm(f, data = df, family = binomial)
  
  broom::tidy(m, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      predictor = predictor,
      OR = round(estimate, 3),
      CI_low = round(conf.low, 3),
      CI_high = round(conf.high, 3),
      p_value = signif(p.value, 3)
    ) %>%
    select(predictor, term, OR, CI_low, CI_high, p_value)
}

univ_table_logistic <- purrr::map_dfr(
  setdiff(candidate_predictors_logistic, "sexo"),  # sexo excluded from screening
  ~run_univ_logistic(modelo_log_data, "aids_status", .x)
) %>%
  arrange(as.numeric(p_value))

cat("\n--- Univariate logistic results (crude OR) ---\n")
print(univ_table_logistic)

# --- 7.3 Variable selection (p <= 0.20), excluding sexo ---
selected_predictors_logistic <- univ_table_logistic %>%
  mutate(p_num = suppressWarnings(as.numeric(p_value))) %>%
  group_by(predictor) %>%
  summarise(min_p = min(p_num, na.rm = TRUE), .groups = "drop") %>%
  filter(min_p <= 0.20) %>%
  pull(predictor)

# Force sexo into the model
selected_predictors_logistic <- unique(c("sexo", selected_predictors_logistic))

cat("\nSelected predictors (p<=0.20) + forced sexo:\n",
    paste(selected_predictors_logistic, collapse = ", "), "\n")

# Fallback safety (if only sexo survives)
if (length(selected_predictors_logistic) == 1) {
  selected_predictors_logistic <- c("sexo", "idade", "Log")
  cat("Only sexo selected. Using fallback:", paste(selected_predictors_logistic, collapse = ", "), "\n")
}

# --- 7.4 Fit multivariable logistic model ---
final_formula_logistic <- as.formula(
  paste("aids_status ~", paste(selected_predictors_logistic, collapse = " + "))
)

modelo_logistico <- glm(
  final_formula_logistic,
  data = modelo_log_data,
  family = binomial
)

cat("\n--- Summary of the Multivariable Logistic Model ---\n")
print(summary(modelo_logistico))

# Adjusted OR table
tabela_or_ajustada <- broom::tidy(
  modelo_logistico,
  exponentiate = TRUE,
  conf.int = TRUE
) %>%
  filter(term != "(Intercept)") %>%
  transmute(
    term,
    AOR = round(estimate, 3),
    CI_low = round(conf.low, 3),
    CI_high = round(conf.high, 3),
    p_value = signif(p.value, 3)
  )

cat("\n--- Adjusted Odds Ratios (Final Model) ---\n")
print(tabela_or_ajustada)

# --- 7.5 Save tables ---
write_xlsx(
  list(
    "Univariate logistic (crude OR)" = univ_table_logistic,
    "Final logistic model (AOR, sexo forced)" = tabela_or_ajustada
  ),
  path = file.path(results_path, "tables_logistic_model_aids.xlsx")
)

cat("Logistic model tables saved.\n\n")

################################################################################
#### TABLE 1A - CHARACTERIZATION OF THE CROSS-SECTIONAL POPULATION ####
################################################################################

cat("--- Table 1A: Cross-sectional population characteristics ---\n")


# Função para formatar mediana (Q1–Q3)
fmt_median_iqr <- function(x, digits = 1) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) return(NA_character_)
  q <- quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  paste0(
    format(round(q[2], digits), nsmall = digits),
    " (",
    format(round(q[1], digits), nsmall = digits),
    "–",
    format(round(q[3], digits), nsmall = digits),
    ")"
  )
}

# Função para sexo F/M
fmt_sex <- function(x) {
  x <- as.factor(x)
  f <- sum(x %in% c("F", "Female", "female", "f"), na.rm = TRUE)
  m <- sum(x %in% c("M", "Male", "male", "m"), na.rm = TRUE)
  paste0(f, "/", m)
}

# Base para a tabela
dados_tabela1 <- dados_analise %>%
  mutate(
    idade = suppressWarnings(as.numeric(idade)),
    cv = suppressWarnings(as.numeric(cv)),
    cd4 = suppressWarnings(as.numeric(cd4)),
    cd8 = suppressWarnings(as.numeric(cd8)),
    cd4_cd8 = suppressWarnings(as.numeric(`CD4/CD8`)),
    il2 = suppressWarnings(as.numeric(IL2)),
    il6 = suppressWarnings(as.numeric(IL6)),
    sil2r = suppressWarnings(as.numeric(il2r)),
    sil6r = suppressWarnings(as.numeric(il6r)),
    scd40l = suppressWarnings(as.numeric(cd40l_pg_ml))
  )

# Resumo por grupo
tabela1A_long <- dados_tabela1 %>%
  group_by(cd4_group) %>%
  summarise(
    N = n(),
    `Sex (F/M)` = fmt_sex(sexo),
    `Age (years)` = fmt_median_iqr(idade, 1),
    `Viral load (copies/mL)` = fmt_median_iqr(cv, 0),
    `CD4+ T cells (cells/mm³)` = fmt_median_iqr(cd4, 0),
    `CD8+ T cells (cells/mm³)` = fmt_median_iqr(cd8, 0),
    `CD4/CD8 ratio` = fmt_median_iqr(cd4_cd8, 2),
    `sIL-2 (pg/mL)` = fmt_median_iqr(il2, 2),
    `sIL-6 (pg/mL)` = fmt_median_iqr(il6, 2),
    `sIL-2R (pg/mL)` = fmt_median_iqr(sil2r, 0),
    `sIL-6R (pg/mL)` = fmt_median_iqr(sil6r, 0),
    `sCD40L (pg/mL)` = fmt_median_iqr(scd40l, 0),
    .groups = "drop"
  )

print(tabela1A_long)

################################################################################
#### P-VALUES ####
################################################################################

# Função para Kruskal-Wallis
kw_p <- function(var) {
  f <- as.formula(paste0(var, " ~ cd4_group"))
  p <- kruskal.test(f, data = dados_tabela1)$p.value
  signif(p, 3)
}

# Qui-quadrado para sexo
sexo_tab <- table(dados_tabela1$sexo, dados_tabela1$cd4_group)
sexo_p <- tryCatch(chisq.test(sexo_tab)$p.value, error = function(e) NA_real_)

# Tabela de p-values
tabela1A_p <- tibble(
  Variable = c(
    "N",
    "Sex (F/M)",
    "Age (years)",
    "Viral load (copies/mL)",
    "CD4+ T cells (cells/mm³)",
    "CD8+ T cells (cells/mm³)",
    "CD4/CD8 ratio",
    "sIL-2 (pg/mL)",
    "sIL-6 (pg/mL)",
    "sIL-2R (pg/mL)",
    "sIL-6R (pg/mL)",
    "sCD40L (pg/mL)"
  ),
  `P-value` = c(
    NA,
    signif(sexo_p, 3),
    kw_p("idade"),
    kw_p("cv"),
    kw_p("cd4"),
    kw_p("cd8"),
    kw_p("cd4_cd8"),
    kw_p("il2"),
    kw_p("il6"),
    kw_p("sil2r"),
    kw_p("sil6r"),
    kw_p("scd40l")
  )
)

################################################################################
#### TRANSFORM TO ARTICLE FORMAT ####
################################################################################

# Converter tudo para character antes do pivot_longer
tabela1A <- tabela1A_long %>%
  mutate(across(-cd4_group, as.character)) %>%
  pivot_longer(cols = -cd4_group, names_to = "Variable", values_to = "Value") %>%
  pivot_wider(names_from = cd4_group, values_from = Value)

# Ordenar variáveis
ordem_variaveis <- c(
  "N",
  "Sex (F/M)",
  "Age (years)",
  "Viral load (copies/mL)",
  "CD4+ T cells (cells/mm³)",
  "CD8+ T cells (cells/mm³)",
  "CD4/CD8 ratio",
  "sIL-2 (pg/mL)",
  "sIL-6 (pg/mL)",
  "sIL-2R (pg/mL)",
  "sIL-6R (pg/mL)",
  "sCD40L (pg/mL)"
)

tabela1A <- tabela1A %>%
  mutate(Variable = factor(Variable, levels = ordem_variaveis)) %>%
  arrange(Variable)

# Juntar p-values
tabela1A_final <- tabela1A %>%
  left_join(tabela1A_p, by = "Variable")

print(tabela1A_final)

################################################################################
#### EXPORT ####
################################################################################

write_xlsx(
  list("Table1A_cross_sectional" = tabela1A_final),
  path = file.path(results_path, "table1A_cross_sectional_characteristics.xlsx")
)

cat("Table 1A saved to:", file.path(results_path, "table1A_cross_sectional_characteristics.xlsx"), "\n")


################################################################################
#### 8. TREATMENT ANALYSIS (BT vs AT at 2m and 4m) + ADD IL2 and IL6 ####
################################################################################

cat("--- 8. Treatment analysis (paired BT vs AT at 2m and 4m) ---\n")

# ---- Helper: robust numeric parsing for mixed PT/EN formats (cell-by-cell) ----
parse_num_robust <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x == ""] <- NA
  
  out <- sapply(x, function(s) {
    if (is.na(s)) return(NA_real_)
    
    s <- gsub("\\s+", "", s)
    
    has_comma <- grepl(",", s)
    has_dot   <- grepl("\\.", s)
    
    # Case: both comma and dot present -> decide decimal by last separator
    if (has_comma && has_dot) {
      last_comma <- max(gregexpr(",", s)[[1]])
      last_dot   <- max(gregexpr("\\.", s)[[1]])
      
      if (last_dot > last_comma) {
        # decimal = dot, comma = thousands (e.g., 4,832.74)
        s2 <- gsub(",", "", s)
        return(suppressWarnings(as.numeric(s2)))
      } else {
        # decimal = comma, dot = thousands (e.g., 4.832,74)
        s2 <- gsub("\\.", "", s)
        s2 <- gsub(",", ".", s2)
        return(suppressWarnings(as.numeric(s2)))
      }
    }
    
    # Case: only comma -> decimal comma
    if (has_comma && !has_dot) {
      s2 <- gsub(",", ".", s)
      return(suppressWarnings(as.numeric(s2)))
    }
    
    # Case: only dot OR no separators -> standard numeric
    suppressWarnings(as.numeric(s))
  })
  
  as.numeric(out)
}

# Read ALL as text first to avoid Excel auto-conversions that differ by locale
hiv_tratamento_wide_raw <- readxl::read_excel(data_path_treatment, col_types = "text")

# Create Patient_ID and parse numeric columns robustly
hiv_tratamento_wide <- hiv_tratamento_wide_raw %>%
  dplyr::mutate(Patient_ID = dplyr::row_number(), .before = 1) %>%
  dplyr::mutate(dplyr::across(-Patient_ID, parse_num_robust))

# ---- Identify timepoint columns (2m and 4m) ----
cols_2m <- stringr::str_subset(names(hiv_tratamento_wide), "(bt|at)2m$")
cols_4m <- stringr::str_subset(names(hiv_tratamento_wide), "(bt|at)4m$")

# ---- Build long format (robust to biomarker names with underscores) ----
make_long_time <- function(df_wide, cols_time, suffix) {
  df_wide %>%
    dplyr::select(Patient_ID, dplyr::all_of(cols_time)) %>%
    tidyr::pivot_longer(
      cols = -Patient_ID,
      names_to = "var",
      values_to = "Dosage",
      values_drop_na = FALSE
    ) %>%
    dplyr::mutate(
      Biomarker = stringr::str_replace(var, paste0("_(bt|at)", suffix, "$"), ""),
      Timecode  = stringr::str_extract(var, paste0("(bt|at)", suffix, "$")),
      Timepoint = factor(
        dplyr::case_when(
          stringr::str_starts(Timecode, "bt") ~ "Before (BT)",
          stringr::str_starts(Timecode, "at") ~ "After (AT)",
          TRUE ~ NA_character_
        ),
        levels = c("Before (BT)", "After (AT)")
      )
    ) %>%
    dplyr::select(-var, -Timecode) %>%
    dplyr::filter(!is.na(Timepoint))
}

dados_2m <- make_long_time(hiv_tratamento_wide, cols_2m, "2m")
dados_4m <- make_long_time(hiv_tratamento_wide, cols_4m, "4m")

cat("Available treated biomarkers (2m): ", paste(sort(unique(dados_2m$Biomarker)), collapse = ", "), "\n")
cat("Available treated biomarkers (4m): ", paste(sort(unique(dados_4m$Biomarker)), collapse = ", "), "\n")

# ---- OPTIONAL sanity check: verify BT != AT for cd40l if present ----
if (all(c("cd40l_bt2m", "cd40l_at2m") %in% names(hiv_tratamento_wide))) {
  chk <- hiv_tratamento_wide %>%
    dplyr::transmute(
      bt = cd40l_bt2m,
      at = cd40l_at2m,
      ok_pair = !is.na(bt) & !is.na(at),
      equal = ok_pair & (bt == at)
    ) %>%
    dplyr::summarise(
      n_pairs = sum(ok_pair),
      n_equal = sum(equal),
      prop_equal = ifelse(n_pairs > 0, n_equal / n_pairs, NA_real_)
    )
  cat("\nSanity check cd40l 2m (pairs/equal):\n")
  print(chk)
}

# ---- Paired test (complete pairs only; handles all-zero diffs) ----
run_single_paired_test_no_impute <- function(data, marker_name) {
  data_filtered <- data %>% dplyr::filter(Biomarker == marker_name)
  
  if (nrow(data_filtered) == 0) {
    return(tibble(
      group1 = "Before (BT)", group2 = "After (AT)",
      n = NA_integer_, statistic = NA_real_, p = NA_real_, p.signif = "NA"
    ))
  }
  
  wide <- data_filtered %>%
    tidyr::pivot_wider(id_cols = Patient_ID, names_from = Timepoint, values_from = Dosage)
  
  if (!all(c("Before (BT)", "After (AT)") %in% names(wide))) {
    return(tibble(
      group1 = "Before (BT)", group2 = "After (AT)",
      n = NA_integer_, statistic = NA_real_, p = NA_real_, p.signif = "NA"
    ))
  }
  
  wide_complete <- wide %>% tidyr::drop_na(`Before (BT)`, `After (AT)`)
  
  if (nrow(wide_complete) < 2) {
    return(tibble(
      group1 = "Before (BT)", group2 = "After (AT)",
      n = nrow(wide_complete), statistic = NA_real_, p = NA_real_, p.signif = "NA"
    ))
  }
  
  # If all paired differences are exactly zero, Wilcoxon returns NA; handle explicitly
  diffs <- wide_complete$`After (AT)` - wide_complete$`Before (BT)`
  if (all(diffs == 0)) {
    return(tibble(
      group1 = "Before (BT)", group2 = "After (AT)",
      n = nrow(wide_complete), statistic = 0, p = 1, p.signif = "ns"
    ))
  }
  
  long_complete <- wide_complete %>%
    tidyr::pivot_longer(
      cols = c(`Before (BT)`, `After (AT)`),
      names_to = "Timepoint",
      values_to = "Dosage"
    )
  
  long_complete %>%
    rstatix::wilcox_test(Dosage ~ Timepoint, paired = TRUE, detailed = TRUE, ref.group = "Before (BT)") %>%
    rstatix::add_significance()
}

# ---- General plotting function: choose transform for y-axis ----
create_facet_plot_general <- function(data_2m, data_4m, marker,
                                      y_trans = c("identity", "log10"),
                                      y_label = NULL,
                                      title = NULL) {
  y_trans <- match.arg(y_trans)
  
  data_plot <- dplyr::bind_rows(
    data_2m %>% dplyr::filter(Biomarker == marker) %>% dplyr::mutate(Time = "2 Months"),
    data_4m %>% dplyr::filter(Biomarker == marker) %>% dplyr::mutate(Time = "4 Months")
  )
  
  if (nrow(data_plot) == 0) {
    warning(paste0("Marker '", marker, "' not found in treated dataset. Returning NULL."))
    return(NULL)
  }
  
  data_complete <- data_plot %>% tidyr::drop_na(Dosage)
  
  test_results <- dplyr::bind_rows(
    run_single_paired_test_no_impute(data_2m, marker) %>% dplyr::mutate(Time = "2 Months"),
    run_single_paired_test_no_impute(data_4m, marker) %>% dplyr::mutate(Time = "4 Months")
  )
  
  # y.position for p-value labels
  if (y_trans == "log10") {
    y_max <- log10(max(data_plot$Dosage, na.rm = TRUE))
    test_results <- test_results %>% dplyr::mutate(y.position = y_max + 0.5)
  } else {
    y_max <- max(data_plot$Dosage, na.rm = TRUE)
    test_results <- test_results %>% dplyr::mutate(y.position = y_max * 1.10)
  }
  
  if (is.null(y_label)) y_label <- marker
  if (is.null(title)) title <- paste0(marker, " across treatment time")
  
  ggplot2::ggplot(data_plot, ggplot2::aes(x = Timepoint, y = Dosage, fill = Timepoint)) +
    ggplot2::geom_boxplot(na.rm = TRUE) +
    ggplot2::geom_line(
      data = data_complete,
      ggplot2::aes(group = Patient_ID),
      color = "gray",
      alpha = 0.5,
      position = ggplot2::position_dodge(width = 0.2)
    ) +
    ggplot2::geom_point(
      data = data_complete,
      ggplot2::aes(group = Patient_ID),
      position = ggplot2::position_dodge(width = 0.2)
    ) +
    {
      if (y_trans == "log10") {
        ggplot2::scale_y_continuous(
          trans = "log10",
          breaks = scales::trans_breaks("log10", function(x) 10^x),
          labels = scales::trans_format("log10", scales::math_format(10^.x))
        )
      } else {
        ggplot2::scale_y_continuous()
      }
    } +
    ggpubr::stat_pvalue_manual(
      test_results,
      label = "p.signif",
      y.position = "y.position",
      tip.length = 0.01,
      size = 5
    ) +
    ggplot2::facet_wrap(~Time) +
    base_theme +
    ggplot2::labs(title = title, x = "Treatment status", y = y_label)
}


# ---- Labels for treated markers (extend as needed) ----
biomarker_labels <- c(
  "il2r"  = "sIL-2R (pg/mL)",
  "il6r"  = "sIL-6R (pg/mL)",
  "cd40l" = "sCD40L (pg/mL)",
  "il2"   = "sIL-2 (pg/mL)",
  "il6"   = "sIL-6 (pg/mL)",
  "cd4"   = "CD4 (cells/mm³)",
  "cv"    = "Viral load (copies/mL)"
)

# ---- Core treated plots (existing) ----
p10 <- create_facet_plot_general(
  dados_2m, dados_4m, marker = "il2r",
  y_trans = "log10",
  y_label = biomarker_labels["il2r"],
  title = "sIL-2R before vs after treatment"
)

p11 <- create_facet_plot_general(
  dados_2m, dados_4m, marker = "il6r",
  y_trans = "log10",
  y_label = biomarker_labels["il6r"],
  title = "sIL-6R before vs after treatment"
)

p12 <- create_facet_plot_general(
  dados_2m, dados_4m, marker = "cd40l",
  y_trans = "log10",
  y_label = biomarker_labels["cd40l"],
  title = "sCD40L before vs after treatment"
)

# ---- NEW: IL2 and IL6 (you added these columns) ----
p15 <- create_facet_plot_general(
  dados_2m, dados_4m, marker = "il2",
  y_trans = "identity",
  y_label = biomarker_labels["il2"],
  title = "sIL-2 before vs after treatment"
)

p16 <- create_facet_plot_general(
  dados_2m, dados_4m, marker = "il6",
  y_trans = "identity",
  y_label = biomarker_labels["il6"],
  title = "sIL-6 before vs after treatment"
)

# ---- Optional: CD4 and CV (only if present in treated dataset) ----
available_markers <- sort(unique(c(unique(dados_2m$Biomarker), unique(dados_4m$Biomarker))))

p13 <- if ("cd4" %in% available_markers) {
  create_facet_plot_general(
    dados_2m, dados_4m, marker = "cd4",
    y_trans = "identity",
    y_label = biomarker_labels["cd4"],
    title = "CD4 before vs after treatment"
  )
} else {
  cat("NOTE: 'cd4' not found in treated dataset. Skipping p13.\n")
  NULL
}

p14 <- if ("cv" %in% available_markers) {
  create_facet_plot_general(
    dados_2m, dados_4m, marker = "cv",
    y_trans = "log10",
    y_label = biomarker_labels["cv"],
    title = "Viral load before vs after treatment"
  )
} else {
  cat("NOTE: 'cv' not found in treated dataset. Skipping p14.\n")
  NULL
}

# ---- Final panel (robust): patchwork ----
plot_list_treat <- list(p10, p11, p15, p16,  p12, p13, p14)

# Ensure no NULL (safety)
plot_list_treat <- plot_list_treat[!sapply(plot_list_treat, is.null)]

# Build panel with patchwork
final_plots_treatment <- patchwork::wrap_plots(
  plot_list_treat,
  ncol = 2
) +
  patchwork::plot_annotation(
    tag_levels = "A"
  )

print(final_plots_treatment)

# Save reliably with ggsave (patchwork object is ggplot-compatible)
ggsave(
  filename = file.path(results_path, "figure_treatment_panel.png"),
  plot = final_plots_treatment,
  width = 16, height = 18, dpi = 300
)

cat("Treatment panel saved:", file.path(results_path, "figure_treatment_panel.png"), "\n\n")
