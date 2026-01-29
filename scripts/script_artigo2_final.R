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

# Helper: consistent theme
base_theme <- theme_minimal(base_size = 16) + theme(legend.position = "none")

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
  stat_compare_means(method = "kruskal.test", label.y.npc = 0.9) +
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
  stat_compare_means(method = "kruskal.test", label.y.npc = 0.95) +
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
  stat_compare_means(method = "kruskal.test", label.y = 125) +
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
) + labs(title = "PCA of HIV patients (immune markers)") + theme_minimal(base_size = 14)

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
) + labs(title = paste0("K-means clustering (k=", k, ")"))

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
#### 8. TREATMENT ANALYSIS (BT vs AT at 2m and 4m) + ADD IL2 and IL6 ####
################################################################################

cat("--- 8. Treatment analysis (paired BT vs AT at 2m and 4m) ---\n")

# ---- Helper: robust numeric parsing for mixed PT/EN formats ----
parse_num_smart <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x == ""] <- NA
  
  has_comma <- any(grepl(",", x), na.rm = TRUE)
  if (has_comma) {
    # likely PT-BR style: 1.234,56
    x2 <- gsub("\\.", "", x)   # remove thousands dots
    x2 <- gsub(",", ".", x2)   # decimal comma -> dot
    suppressWarnings(as.numeric(x2))
  } else {
    # likely EN style: 1234.56
    suppressWarnings(as.numeric(x))
  }
}

hiv_tratamento_wide_raw <- read_excel(data_path_treatment)

# Create Patient_ID (if your file already has an ID column, you can edit this)
hiv_tratamento_wide <- hiv_tratamento_wide_raw %>%
  rowid_to_column(var = "Patient_ID") %>%
  mutate(across(-Patient_ID, parse_num_smart))

# ---- Identify timepoint columns (2m and 4m) ----
cols_2m <- str_subset(names(hiv_tratamento_wide), "(bt|at)2m$")
cols_4m <- str_subset(names(hiv_tratamento_wide), "(bt|at)4m$")

# ---- Build long format (robust to biomarker names with underscores) ----
make_long_time <- function(df_wide, cols_time, suffix) {
  df_wide %>%
    select(Patient_ID, all_of(cols_time)) %>%
    pivot_longer(
      cols = -Patient_ID,
      names_to = "var",
      values_to = "Dosage",
      values_drop_na = FALSE
    ) %>%
    mutate(
      Biomarker = str_replace(var, paste0("_(bt|at)", suffix, "$"), ""),
      Timecode  = str_extract(var, paste0("(bt|at)", suffix, "$")),
      Timepoint = factor(
        case_when(
          str_starts(Timecode, "bt") ~ "Before (BT)",
          str_starts(Timecode, "at") ~ "After (AT)",
          TRUE ~ NA_character_
        ),
        levels = c("Before (BT)", "After (AT)")
      )
    ) %>%
    select(-var, -Timecode) %>%
    filter(!is.na(Timepoint))
}

dados_2m <- make_long_time(hiv_tratamento_wide, cols_2m, "2m")
dados_4m <- make_long_time(hiv_tratamento_wide, cols_4m, "4m")

cat("Available treated biomarkers (2m): ", paste(sort(unique(dados_2m$Biomarker)), collapse = ", "), "\n")
cat("Available treated biomarkers (4m): ", paste(sort(unique(dados_4m$Biomarker)), collapse = ", "), "\n")

# ---- Paired test (complete pairs only) ----
run_single_paired_test_no_impute <- function(data, marker_name) {
  data_filtered <- data %>% filter(Biomarker == marker_name)
  
  if (nrow(data_filtered) == 0) {
    return(tibble(
      group1 = "Before (BT)", group2 = "After (AT)",
      n = NA_integer_, statistic = NA_real_, p = NA_real_, p.signif = "NA"
    ))
  }
  
  wide <- data_filtered %>%
    pivot_wider(id_cols = Patient_ID, names_from = Timepoint, values_from = Dosage)
  
  if (!all(c("Before (BT)", "After (AT)") %in% names(wide))) {
    return(tibble(
      group1 = "Before (BT)", group2 = "After (AT)",
      n = NA_integer_, statistic = NA_real_, p = NA_real_, p.signif = "NA"
    ))
  }
  
  wide_complete <- wide %>% drop_na(`Before (BT)`, `After (AT)`)
  
  if (nrow(wide_complete) < 2) {
    return(tibble(
      group1 = "Before (BT)", group2 = "After (AT)",
      n = nrow(wide_complete), statistic = NA_real_, p = NA_real_, p.signif = "NA"
    ))
  }
  
  long_complete <- wide_complete %>%
    pivot_longer(
      cols = c(`Before (BT)`, `After (AT)`),
      names_to = "Timepoint",
      values_to = "Dosage"
    )
  
  long_complete %>%
    wilcox_test(Dosage ~ Timepoint, paired = TRUE, detailed = TRUE, ref.group = "Before (BT)") %>%
    add_significance()
}

# ---- General plotting function: choose transform for y-axis ----
create_facet_plot_general <- function(data_2m, data_4m, marker,
                                      y_trans = c("identity", "log10"),
                                      y_label = NULL,
                                      title = NULL) {
  y_trans <- match.arg(y_trans)
  
  data_plot <- bind_rows(
    data_2m %>% filter(Biomarker == marker) %>% mutate(Time = "2 Months"),
    data_4m %>% filter(Biomarker == marker) %>% mutate(Time = "4 Months")
  )
  
  # If marker not present in both timepoints, still try plot with what exists
  if (nrow(data_plot) == 0) {
    warning(paste0("Marker '", marker, "' not found in treated dataset. Returning NULL."))
    return(NULL)
  }
  
  data_complete <- data_plot %>% drop_na(Dosage)
  
  test_results <- bind_rows(
    run_single_paired_test_no_impute(data_2m, marker) %>% mutate(Time = "2 Months"),
    run_single_paired_test_no_impute(data_4m, marker) %>% mutate(Time = "4 Months")
  )
  
  # y.position for pvalue labels
  if (y_trans == "log10") {
    y_max <- log10(max(data_plot$Dosage, na.rm = TRUE))
    test_results <- test_results %>% mutate(y.position = y_max + 0.5)
  } else {
    y_max <- max(data_plot$Dosage, na.rm = TRUE)
    test_results <- test_results %>% mutate(y.position = y_max * 1.10)
  }
  
  if (is.null(y_label)) y_label <- marker
  if (is.null(title)) title <- paste0(marker, " across treatment time")
  
  p <- ggplot(data_plot, aes(x = Timepoint, y = Dosage, fill = Timepoint)) +
    geom_boxplot(na.rm = TRUE) +
    geom_line(
      data = data_complete,
      aes(group = Patient_ID),
      color = "gray",
      alpha = 0.5,
      position = position_dodge(width = 0.2)
    ) +
    geom_point(
      data = data_complete,
      aes(group = Patient_ID),
      position = position_dodge(width = 0.2)
    ) +
    {
      if (y_trans == "log10") {
        scale_y_continuous(
          trans = "log10",
          breaks = trans_breaks("log10", function(x) 10^x),
          labels = trans_format("log10", math_format(10^.x))
        )
      } else {
        scale_y_continuous()
      }
    } +
    stat_pvalue_manual(
      test_results,
      label = "p.signif",
      y.position = "y.position",
      tip.length = 0.01,
      size = 5
    ) +
    facet_wrap(~Time) +
    theme_minimal(base_size = 18) +
    labs(title = title, x = "Treatment status", y = y_label) +
    theme(legend.position = "none", plot.title = element_text(face = "bold"))
  
  p
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
