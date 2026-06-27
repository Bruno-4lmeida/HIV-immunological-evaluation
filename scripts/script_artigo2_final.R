#### R SCRIPT FOR: ####
# IL-6 AND VIRAL LOAD ARE ASSOCIATED WITH CD4-DEFINED ADVANCED IMMUNODEFICIENCY,
# WHILE ART RAPIDLY REDUCES SIL-2R LEVELS IN HIV INFECTION
#
# Author: Bruno Almeida Silva

#### 0. SETUP ####

#### Packages
required_packages <- c(
  "readxl",
  "tidyverse",
  "ggpubr",
  "FactoMineR",
  "factoextra",
  "corrplot",
  "broom",
  "writexl",
  "patchwork",
  "rstatix",
  "scales",
  "car",
  "ResourceSelection",
  "pscl",
  "pROC",
  "caret",
  "cluster"
)

installed_packages <- rownames(installed.packages())

packages_to_install <- required_packages[!(required_packages %in% installed_packages)]

if (length(packages_to_install) > 0) {
  install.packages(packages_to_install)
}

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
library(car)
library(ResourceSelection)
library(pscl)
library(pROC)
library(caret)
library(cluster)

# ---- Global options ----
set.seed(123)

# ---- Paths: EDIT THESE TWO LINES IF NEEDED ----
data_path_main      <- "C:/Users/bruno/OneDrive/Documentos/BRUNO/DOUTORADO/EXPERIMENTOS_DOC/hiv_data.xlsx"
data_path_treatment <- "C:/Users/bruno/OneDrive/Documentos/BRUNO/DOUTORADO/EXPERIMENTOS_DOC/analise_tratados.xlsx"

# ---- Output folder ----
results_path <- "artigo2_results"

if (!dir.exists(results_path)) {
  dir.create(results_path, recursive = TRUE)
}

cat("=== Script started. Outputs will be saved in:", results_path, "===\n\n")


#### 0.1. THEMES AND HELPER FUNCTIONS ####

base_theme <- theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 16, hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 14)
  )

base_theme_with_legend <- theme_minimal(base_size = 16) +
  theme(
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 16, hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 14)
  )

# ---- Function to format p-values ----
format_p <- function(p) {
  case_when(
    is.na(p) ~ NA_character_,
    p < 0.0001 ~ "<0.0001",
    TRUE ~ formatC(p, format = "f", digits = 4)
  )
}

# ---- Function to generate significance labels ----
p_to_signif <- function(p) {
  case_when(
    is.na(p) ~ "NA",
    p < 0.0001 ~ "****",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

# ---- Function to format median and IQR ----
fmt_median_iqr <- function(x, digits = 1) {
  x <- suppressWarnings(as.numeric(x))
  
  if (all(is.na(x))) {
    return(NA_character_)
  }
  
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

# ---- Function to format sex as F/M ----
fmt_sex <- function(x) {
  x <- as.character(x)
  
  f <- sum(x %in% c("F", "Female", "female", "f", "Feminino", "feminino"), na.rm = TRUE)
  m <- sum(x %in% c("M", "Male", "male", "m", "Masculino", "masculino"), na.rm = TRUE)
  
  paste0(f, "/", m)
}

# ---- Robust numeric parser for mixed decimal formats ----
parse_num_robust <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x == ""] <- NA
  
  out <- sapply(x, function(s) {
    
    if (is.na(s)) {
      return(NA_real_)
    }
    
    s <- gsub("\\s+", "", s)
    
    has_comma <- grepl(",", s)
    has_dot   <- grepl("\\.", s)
    
    # Both comma and dot are present: decide decimal separator by last separator
    if (has_comma && has_dot) {
      
      last_comma <- max(gregexpr(",", s)[[1]])
      last_dot   <- max(gregexpr("\\.", s)[[1]])
      
      if (last_dot > last_comma) {
        # Decimal separator is dot, comma is thousands separator
        s2 <- gsub(",", "", s)
        return(suppressWarnings(as.numeric(s2)))
      } else {
        # Decimal separator is comma, dot is thousands separator
        s2 <- gsub("\\.", "", s)
        s2 <- gsub(",", ".", s2)
        return(suppressWarnings(as.numeric(s2)))
      }
    }
    
    # Only comma is present: decimal comma
    if (has_comma && !has_dot) {
      s2 <- gsub(",", ".", s)
      return(suppressWarnings(as.numeric(s2)))
    }
    
    # Only dot or no separator
    suppressWarnings(as.numeric(s))
  })
  
  as.numeric(out)
}


#### 1. DATA CLEANING AND PREPARATION ####

cat("--- 1. Data cleaning and preparation ---\n")

hiv_data <- read_excel(data_path_main)

dados_limpos <- hiv_data %>%
  filter(!is.na(cd4)) %>%
  mutate(
    cv    = suppressWarnings(as.numeric(as.character(cv))),
    Log   = suppressWarnings(as.numeric(as.character(Log))),
    idade = suppressWarnings(as.numeric(as.character(idade))),
    cd4   = suppressWarnings(as.numeric(as.character(cd4))),
    cd8   = suppressWarnings(as.numeric(as.character(cd8))),
    sexo  = as.factor(sexo)
  ) %>%
  filter(idade < 100)

cat("Data cleaned. Kept", nrow(dados_limpos), "rows.\n\n")


#### 2. FEATURE ENGINEERING ####

cat("--- 2. Feature engineering ---\n")

dados_analise <- dados_limpos %>%
  mutate(
    cd4_group = factor(
      case_when(
        cd4 <= 200 ~ "<= 200",
        cd4 > 200 & cd4 < 350 ~ ">200 & <350",
        cd4 >= 350 ~ ">= 350",
        TRUE ~ NA_character_
      ),
      levels = c("<= 200", ">200 & <350", ">= 350")
    ),
    IL2_IL2R_ratio = IL2 / il2r,
    IL6_IL6R_ratio = IL6 / il6r
  )

# ---- Define the preferred sCD40L variable ----
# The manuscript figures and Table 1 use cd40l_pg_ml.
# If cd40l_pg_ml exists, it is preferred as the final concentration in pg/mL.
# If not, the script falls back to CD40L.

if ("cd40l_pg_ml" %in% names(dados_analise)) {
  cd40l_var <- "cd40l_pg_ml"
} else if ("CD40L" %in% names(dados_analise)) {
  cd40l_var <- "CD40L"
} else if ("cd40l" %in% names(dados_analise)) {
  cd40l_var <- "cd40l"
} else {
  stop("No sCD40L variable found. Check column names in dados_analise.")
}

cat("sCD40L variable used:", cd40l_var, "\n")

cat("CD4 strata counts:\n")
print(table(dados_analise$cd4_group, useNA = "ifany"))
cat("\n")


#### 3. TABLE 1: BASELINE CHARACTERISTICS ####

cat("--- 3. Table 1: baseline characteristics ---\n")

dados_tabela1 <- dados_analise %>%
  mutate(
    idade    = suppressWarnings(as.numeric(idade)),
    cv       = suppressWarnings(as.numeric(cv)),
    Log      = suppressWarnings(as.numeric(Log)),
    cd4      = suppressWarnings(as.numeric(cd4)),
    cd8      = suppressWarnings(as.numeric(cd8)),
    cd4_cd8  = suppressWarnings(as.numeric(`CD4/CD8`)),
    il2      = suppressWarnings(as.numeric(IL2)),
    il6      = suppressWarnings(as.numeric(IL6)),
    sil2r    = suppressWarnings(as.numeric(il2r)),
    sil6r    = suppressWarnings(as.numeric(il6r)),
    scd40l   = suppressWarnings(as.numeric(.data[[cd40l_var]]))
  )

# ---- Summary by CD4 group ----
tabela1_long <- dados_tabela1 %>%
  group_by(cd4_group) %>%
  summarise(
    N = n(),
    `Sex (F/M)` = fmt_sex(sexo),
    `Age (years)` = fmt_median_iqr(idade, 1),
    `Viral load (copies/mL)` = fmt_median_iqr(cv, 0),
    `Viral load (log10 copies/mL)` = fmt_median_iqr(Log, 2),
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

print(tabela1_long)

# ---- P-values for Table 1 ----
kw_p <- function(var) {
  f <- as.formula(paste0(var, " ~ cd4_group"))
  p <- kruskal.test(f, data = dados_tabela1)$p.value
  p
}

sexo_tab <- table(dados_tabela1$sexo, dados_tabela1$cd4_group)

sexo_p <- tryCatch(
  chisq.test(sexo_tab)$p.value,
  error = function(e) NA_real_
)

tabela1_p <- tibble(
  Variable = c(
    "N",
    "Sex (F/M)",
    "Age (years)",
    "Viral load (copies/mL)",
    "Viral load (log10 copies/mL)",
    "CD4+ T cells (cells/mm³)",
    "CD8+ T cells (cells/mm³)",
    "CD4/CD8 ratio",
    "sIL-2 (pg/mL)",
    "sIL-6 (pg/mL)",
    "sIL-2R (pg/mL)",
    "sIL-6R (pg/mL)",
    "sCD40L (pg/mL)"
  ),
  p_raw = c(
    NA_real_,
    sexo_p,
    kw_p("idade"),
    kw_p("cv"),
    kw_p("Log"),
    kw_p("cd4"),
    kw_p("cd8"),
    kw_p("cd4_cd8"),
    kw_p("il2"),
    kw_p("il6"),
    kw_p("sil2r"),
    kw_p("sil6r"),
    kw_p("scd40l")
  )
) %>%
  mutate(
    `P-value` = format_p(p_raw)
  ) %>%
  select(Variable, `P-value`, p_raw)

# ---- Transform to article format ----
ordem_variaveis_tabela1 <- c(
  "N",
  "Sex (F/M)",
  "Age (years)",
  "Viral load (copies/mL)",
  "Viral load (log10 copies/mL)",
  "CD4+ T cells (cells/mm³)",
  "CD8+ T cells (cells/mm³)",
  "CD4/CD8 ratio",
  "sIL-2 (pg/mL)",
  "sIL-6 (pg/mL)",
  "sIL-2R (pg/mL)",
  "sIL-6R (pg/mL)",
  "sCD40L (pg/mL)"
)

tabela1_final <- tabela1_long %>%
  mutate(across(-cd4_group, as.character)) %>%
  pivot_longer(
    cols = -cd4_group,
    names_to = "Variable",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = cd4_group,
    values_from = Value
  ) %>%
  mutate(Variable = factor(Variable, levels = ordem_variaveis_tabela1)) %>%
  arrange(Variable) %>%
  left_join(tabela1_p %>% select(Variable, `P-value`), by = "Variable")

print(tabela1_final)

write_xlsx(
  list(
    "Table1_baseline" = tabela1_final,
    "Table1_pvalues_raw" = tabela1_p
  ),
  path = file.path(results_path, "table1_baseline_characteristics.xlsx")
)

cat("Saved Table 1:", file.path(results_path, "table1_baseline_characteristics.xlsx"), "\n\n")


#### 4. GLOBAL KRUSKAL-WALLIS TESTS WITH BH CORRECTION ACROSS BIOMARKERS ####

cat("--- 4. Global Kruskal-Wallis tests with BH correction across biomarkers ---\n")

biomarker_vars <- c(
  "IL2",
  "IL6",
  "il2r",
  "il6r",
  cd40l_var
)

biomarker_labels_global <- c(
  "IL2" = "sIL-2",
  "IL6" = "sIL-6",
  "il2r" = "sIL-2R",
  "il6r" = "sIL-6R"
)

biomarker_labels_global[cd40l_var] <- "sCD40L"

kw_global_biomarkers <- purrr::map_dfr(biomarker_vars, function(v) {
  
  formula_kw <- as.formula(paste0(v, " ~ cd4_group"))
  test_kw <- kruskal.test(formula_kw, data = dados_analise)
  
  tibble(
    variable = v,
    biomarker = biomarker_labels_global[v],
    n_complete = sum(!is.na(dados_analise[[v]]) & !is.na(dados_analise$cd4_group)),
    kruskal_chisq = unname(test_kw$statistic),
    df = unname(test_kw$parameter),
    p_raw = test_kw$p.value
  )
}) %>%
  mutate(
    p_BH_global = p.adjust(p_raw, method = "BH"),
    p_raw_formatted = format_p(p_raw),
    p_BH_global_formatted = format_p(p_BH_global),
    signif_raw = p_to_signif(p_raw),
    signif_BH = p_to_signif(p_BH_global)
  ) %>%
  arrange(p_BH_global)

print(kw_global_biomarkers)

write_xlsx(
  list("KW_global_BH_biomarkers" = kw_global_biomarkers),
  path = file.path(results_path, "table_kw_global_BH_biomarkers.xlsx")
)

cat("Saved:", file.path(results_path, "table_kw_global_BH_biomarkers.xlsx"), "\n\n")


#### 5. PAIRWISE POST-HOC WILCOXON TESTS WITH BH CORRECTION ####

cat("--- 5. Pairwise Wilcoxon post-hoc tests with BH correction ---\n")

# CD4/CD8 is included in pairwise analyses but not in the
# global BH correction across soluble biomarkers.
plot_variables <- c(biomarker_vars, "CD4/CD8")

biomarker_labels_global["CD4/CD8"] <- "CD4/CD8 ratio"

posthoc_biomarkers <- purrr::map_dfr(plot_variables, function(v) {
  
  df_tmp <- dados_analise %>%
    select(cd4_group, value = all_of(v)) %>%
    drop_na()
  
  pw <- pairwise.wilcox.test(
    x = df_tmp$value,
    g = df_tmp$cd4_group,
    p.adjust.method = "BH",
    exact = FALSE
  )
  
  as.data.frame(as.table(pw$p.value)) %>%
    filter(!is.na(Freq)) %>%
    transmute(
      variable = v,
      biomarker = unname(biomarker_labels_global[v]),
      group1 = as.character(Var1),
      group2 = as.character(Var2),
      p_adj_BH = Freq,
      p_adj_BH_formatted = format_p(p_adj_BH),
      significance = p_to_signif(p_adj_BH)
    )
})

print(posthoc_biomarkers)

write_xlsx(
  list("Posthoc_Wilcoxon_BH" = posthoc_biomarkers),
  path = file.path(results_path, "table_posthoc_wilcoxon_BH.xlsx")
)

cat("Saved:", file.path(
  results_path,
  "table_posthoc_wilcoxon_BH.xlsx"
), "\n\n")


#### 6. BOXPLOTS: BIOMARKERS BY CD4 STRATA ####

cat("--- 6. Comparative boxplots by CD4 strata ---\n")

my_comparisons <- list(
  c("<= 200", ">200 & <350"),
  c("<= 200", ">= 350"),
  c(">200 & <350", ">= 350")
)

comparison_key <- function(group1, group2) {
  paste(sort(c(group1, group2)), collapse = "|||")
}

comparison_lookup <- tibble(
  comparison_key = purrr::map_chr(
    my_comparisons,
    ~comparison_key(.x[1], .x[2])
  ),
  comparison_order = seq_along(my_comparisons)
)

make_p_label <- function(p, adjusted = FALSE) {
  
  test_name <- if (adjusted) {
    "Kruskal-Wallis, BH-adjusted "
  } else {
    "Kruskal-Wallis, "
  }
  
  if (p < 0.001) {
    paste0(test_name, "p < 0.001")
  } else {
    paste0(test_name, "p = ", sprintf("%.3f", p))
  }
}

create_cd4_boxplot <- function(data, yvar, title, ylab) {
  
  y_values <- data[[yvar]]
  y_max <- max(y_values, na.rm = TRUE)
  
  if (!is.finite(y_max) || y_max <= 0) {
    y_max <- 1
  }
  
  y_positions <- y_max * c(1.08, 1.20, 1.32)
  
  # Use the globally BH-adjusted KW p value for soluble biomarkers.
  kw_row <- kw_global_biomarkers %>%
    filter(variable == yvar)
  
  if (nrow(kw_row) == 1) {
    
    kw_label <- make_p_label(
      kw_row$p_BH_global[[1]],
      adjusted = TRUE
    )
    
  } else {
    
    # CD4/CD8 was not included in the BH family of soluble biomarkers.
    df_kw <- data %>%
      select(cd4_group, value = all_of(yvar)) %>%
      drop_na()
    
    kw_test <- kruskal.test(value ~ cd4_group, data = df_kw)
    
    kw_label <- make_p_label(
      kw_test$p.value,
      adjusted = FALSE
    )
  }
  
  # Use the previously calculated BH-adjusted pairwise results.
  pairwise_plot <- posthoc_biomarkers %>%
    filter(variable == yvar) %>%
    mutate(
      comparison_key = purrr::map2_chr(
        group1,
        group2,
        comparison_key
      )
    ) %>%
    left_join(comparison_lookup, by = "comparison_key") %>%
    arrange(comparison_order) %>%
    mutate(
      y.position = y_positions[comparison_order]
    )
  
  ggplot(
    data,
    aes(
      x = cd4_group,
      y = .data[[yvar]],
      fill = cd4_group
    )
  ) +
    geom_boxplot(
      outlier.shape = NA,
      alpha = 0.8,
      na.rm = TRUE
    ) +
    geom_jitter(
      width = 0.12,
      alpha = 0.55,
      size = 1.8,
      na.rm = TRUE
    ) +
    annotate(
      "text",
      x = 2,
      y = y_max * 1.48,
      label = kw_label,
      size = 5
    ) +
    ggpubr::stat_pvalue_manual(
      pairwise_plot,
      label = "significance",
      xmin = "group1",
      xmax = "group2",
      y.position = "y.position",
      tip.length = 0.01,
      bracket.size = 0.5,
      size = 5,
      hide.ns = FALSE,
      inherit.aes = FALSE
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.25))
    ) +
    labs(
      title = title,
      x = "CD4+ T-cell stratum (cells/mm³)",
      y = ylab
    ) +
    base_theme
}

p_il2r <- create_cd4_boxplot(
  dados_analise, "il2r",
  "sIL-2R across CD4 strata", "sIL-2R (pg/mL)"
)

p_il6r <- create_cd4_boxplot(
  dados_analise, "il6r",
  "sIL-6R across CD4 strata", "sIL-6R (pg/mL)"
)

p_il2 <- create_cd4_boxplot(
  dados_analise, "IL2",
  "sIL-2 across CD4 strata", "sIL-2 (pg/mL)"
)

p_il6 <- create_cd4_boxplot(
  dados_analise, "IL6",
  "sIL-6 across CD4 strata", "sIL-6 (pg/mL)"
)

p_cd40l <- create_cd4_boxplot(
  dados_analise, cd40l_var,
  "sCD40L across CD4 strata", "sCD40L (pg/mL)"
)

p_cd4_cd8 <- create_cd4_boxplot(
  dados_analise, "CD4/CD8",
  "CD4/CD8 ratio across CD4 strata", "CD4/CD8 ratio"
)

plot_cd4_groups <- patchwork::wrap_plots(
  list(p_il2r, p_il6r, p_il2, p_il6, p_cd40l, p_cd4_cd8),
  ncol = 2
) +
  patchwork::plot_annotation(
    tag_levels = list(c("A)", "B)", "C)", "D)", "E)", "F)"))
  )

print(plot_cd4_groups)

ggsave(
  file.path(results_path, "figure_cd4_groups_p1_to_p6.png"),
  plot = plot_cd4_groups,
  width = 18,
  height = 16,
  dpi = 300
)

#### 7. CREATE TRANSFORMED DATASET FOR LOGISTIC MODELS ####

cat("--- 7. Creating transformed dataset for logistic models ---\n")

dados_modelo2 <- dados_analise %>%
  mutate(
    outcome_cd4_350 = factor(
      if_else(cd4 < 350, "Yes", "No"),
      levels = c("No", "Yes")
    ),
    
    # Log10(x+1) transformation
    log_IL2   = log10(IL2 + 1),
    log_IL6   = log10(IL6 + 1),
    log_il2r  = log10(il2r + 1),
    log_il6r  = log10(il6r + 1),
    log_cd40l = log10(.data[[cd40l_var]] + 1),
    
    # Z-score after log10 transformation
    z_log_IL2   = as.numeric(scale(log_IL2)),
    z_log_IL6   = as.numeric(scale(log_IL6)),
    z_log_il2r  = as.numeric(scale(log_il2r)),
    z_log_il6r  = as.numeric(scale(log_il6r)),
    z_log_cd40l = as.numeric(scale(log_cd40l)),
    
    # Z-score for age and viral load
    z_Log   = as.numeric(scale(Log)),
    z_idade = as.numeric(scale(idade))
  )

cat("Dataset dados_modelo2 created.\n")
cat("Outcome counts for CD4-defined advanced immunodeficiency:\n")
print(table(dados_modelo2$outcome_cd4_350, useNA = "ifany"))
cat("\n")

cat("=== END OF BLOCK 1 ===\n")
#### 8. PRINCIPAL COMPONENT ANALYSIS (PCA) ####

cat("--- 8. Principal component analysis (PCA) ---\n")

# PCA variables:
# We include soluble immune markers, classical cellular immune markers,
# and plasma viral load in log10 scale.
# The data are scaled in the PCA function.

pca_df <- dados_analise %>%
  select(
    cd4_group,
    il2r,
    il6r,
    cd4,
    cd8,
    `CD4/CD8`,
    Log,
    IL2,
    IL6,
    all_of(cd40l_var)
  ) %>%
  drop_na()

# Rename variables to improve PCA outputs and figure labels
pca_df <- pca_df %>%
  rename(
    sIL2R = il2r,
    sIL6R = il6r,
    CD4 = cd4,
    CD8 = cd8,
    CD4_CD8_ratio = `CD4/CD8`,
    Viral_load_log10 = Log,
    sIL2 = IL2,
    sIL6 = IL6,
    sCD40L = all_of(cd40l_var)
  )

pca_data <- pca_df %>%
  select(-cd4_group)

cat("PCA dataset N =", nrow(pca_data), "\n")
cat("Variables included in PCA:\n")
print(names(pca_data))

# Run PCA
pca_result <- FactoMineR::PCA(
  pca_data,
  graph = FALSE,
  scale.unit = TRUE
)

# PCA eigenvalues and variance explained
pca_eigenvalues <- factoextra::get_eigenvalue(pca_result) %>%
  as.data.frame()

# Save the correct percentages for axis labels
dim1_percent <- round(pca_eigenvalues$variance.percent[1], 1)
dim2_percent <- round(pca_eigenvalues$variance.percent[2], 1)

cat("PCA eigenvalues and explained variance:\n")
print(pca_eigenvalues)

cat("Dim1 variance percent:", dim1_percent, "%\n")
cat("Dim2 variance percent:", dim2_percent, "%\n")
cat("Dim1 + Dim2 cumulative variance percent:",
    round(pca_eigenvalues$cumulative.variance.percent[2], 1), "%\n")

write_xlsx(
  list("PCA_eigenvalues" = pca_eigenvalues),
  path = file.path(results_path, "table_pca_eigenvalues.xlsx")
)

cat("Saved:", file.path(results_path, "table_pca_eigenvalues.xlsx"), "\n")

# Variable contributions to PCA dimensions
pca_var <- factoextra::get_pca_var(pca_result)

pca_variable_contributions <- as.data.frame(pca_var$contrib) %>%
  rownames_to_column("variable")

pca_variable_coordinates <- as.data.frame(pca_var$coord) %>%
  rownames_to_column("variable")

pca_variable_cos2 <- as.data.frame(pca_var$cos2) %>%
  rownames_to_column("variable")

write_xlsx(
  list(
    "PCA_variable_contributions" = pca_variable_contributions,
    "PCA_variable_coordinates" = pca_variable_coordinates,
    "PCA_variable_cos2" = pca_variable_cos2
  ),
  path = file.path(results_path, "tables_pca_variables.xlsx")
)

cat("Saved:", file.path(results_path, "tables_pca_variables.xlsx"), "\n")

cat("Variable contributions to Dim1:\n")
print(
  pca_variable_contributions %>%
    select(variable, Dim.1, Dim.2) %>%
    arrange(desc(Dim.1))
)

cat("Variable contributions to Dim2:\n")
print(
  pca_variable_contributions %>%
    select(variable, Dim.1, Dim.2) %>%
    arrange(desc(Dim.2))
)

# PCA individual plot
plot_pca_individuals <- factoextra::fviz_pca_ind(
  pca_result,
  geom.ind = "point",
  col.ind = pca_df$cd4_group,
  palette = NULL,
  addEllipses = TRUE,
  ellipse.level = 0.95,
  legend.title = "CD4 stratum"
) +
  labs(
    title = "PCA of HIV patients based on immune markers",
    x = paste0("Dim1 (", dim1_percent, "%)"),
    y = paste0("Dim2 (", dim2_percent, "%)")
  ) +
  base_theme_with_legend +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

print(plot_pca_individuals)

ggsave(
  filename = file.path(results_path, "figure_pca_individuals.png"),
  plot = plot_pca_individuals,
  width = 9,
  height = 7,
  dpi = 300
)

cat("Saved:", file.path(results_path, "figure_pca_individuals.png"), "\n")

# PCA biplot
plot_pca_biplot <- factoextra::fviz_pca_biplot(
  pca_result,
  repel = TRUE,
  col.var = "black",
  col.ind = "gray40",
  geom.ind = "point"
) +
  labs(
    title = "PCA biplot of immune markers",
    x = paste0("Dim1 (", dim1_percent, "%)"),
    y = paste0("Dim2 (", dim2_percent, "%)")
  ) +
  theme_minimal(base_size = 14)

print(plot_pca_biplot)

ggsave(
  filename = file.path(results_path, "figure_pca_biplot.png"),
  plot = plot_pca_biplot,
  width = 9,
  height = 7,
  dpi = 300
)

cat("Saved:", file.path(results_path, "figure_pca_biplot.png"), "\n\n")

#### 9. K-MEANS CLUSTERING ON ORIGINAL SCALED PCA VARIABLES ####

cat("--- 9. K-means clustering on original scaled PCA variables ---\n")

# The same numeric variables used in PCA are scaled before clustering.
pca_data_scaled <- scale(pca_data)

# Elbow method
plot_elbow <- factoextra::fviz_nbclust(
  pca_data_scaled,
  kmeans,
  method = "wss"
) +
  labs(
    title = "Elbow method for K-means clustering",
    subtitle = "Original scaled immune marker matrix"
  ) +
  theme_minimal(base_size = 14)

print(plot_elbow)

ggsave(
  filename = file.path(results_path, "figure_elbow_method.png"),
  plot = plot_elbow,
  width = 8,
  height = 6,
  dpi = 300
)

cat("Saved:", file.path(results_path, "figure_elbow_method.png"), "\n")

# K-means with k = 3
# We keep k = 3 because this was the solution suggested by the elbow method
# and corresponds to the number of CD4 strata, but this remains exploratory.
set.seed(123)

k_original <- 3

kmeans_result <- kmeans(
  pca_data_scaled,
  centers = k_original,
  nstart = 25
)

# K-means cluster plot
plot_kmeans <- factoextra::fviz_cluster(
  kmeans_result,
  data = pca_data_scaled,
  geom = "point",
  ellipse.type = "convex",
  ggtheme = theme_bw()
) +
  labs(
    title = paste0("K-means clustering, k = ", k_original),
    subtitle = "Clustering based on original scaled immune markers"
  ) +
  base_theme_with_legend

print(plot_kmeans)

ggsave(
  filename = file.path(results_path, "figure_kmeans_clusters.png"),
  plot = plot_kmeans,
  width = 9,
  height = 7,
  dpi = 300
)

cat("Saved:", file.path(results_path, "figure_kmeans_clusters.png"), "\n")

# Cluster composition by CD4 group
cluster_composition_original <- tibble(
  cd4_group = pca_df$cd4_group,
  cluster = factor(kmeans_result$cluster)
) %>%
  count(cluster, cd4_group) %>%
  group_by(cluster) %>%
  mutate(percent_within_cluster = 100 * n / sum(n)) %>%
  ungroup()

print(cluster_composition_original)

write_xlsx(
  list("Cluster_composition_original" = cluster_composition_original),
  path = file.path(results_path, "table_cluster_composition_original.xlsx")
)

cat("Saved:", file.path(results_path, "table_cluster_composition_original.xlsx"), "\n\n")


#### 10. HIERARCHICAL CLUSTERING ####

cat("--- 10. Hierarchical clustering ---\n")

hclust_result <- hclust(
  dist(pca_data_scaled),
  method = "ward.D2"
)

plot_hclust <- factoextra::fviz_dend(
  hclust_result,
  k = k_original,
  cex = 0.6,
  rect = TRUE,
  rect_fill = TRUE
) +
  labs(
    title = paste0("Hierarchical clustering, Ward.D2, k = ", k_original)
  ) +
  theme_minimal(base_size = 14)

print(plot_hclust)

ggsave(
  filename = file.path(results_path, "figure_hclust_dendrogram.png"),
  plot = plot_hclust,
  width = 10,
  height = 6,
  dpi = 300
)

cat("Saved:", file.path(results_path, "figure_hclust_dendrogram.png"), "\n\n")


#### 11. CLUSTERING ON PCA SCORES ####

cat("--- 11. K-means clustering on PCA scores ---\n")

# Choice of PCs:
# Here we use the first 5 principal components.
# This retains the dominant multivariate structure while reducing noise.
m_pcs <- 5

pc_scores <- as.data.frame(pca_result$ind$coord[, 1:m_pcs, drop = FALSE])

pc_scores$cd4_group <- pca_df$cd4_group

cat("Number of PCs used for clustering:", m_pcs, "\n")

# Elbow method on PC space
plot_elbow_pc <- factoextra::fviz_nbclust(
  pc_scores %>% select(-cd4_group),
  kmeans,
  method = "wss"
) +
  labs(
    title = "Elbow method for K-means clustering",
    subtitle = paste0("PC1–PC", m_pcs, " scores")
  ) +
  theme_minimal(base_size = 14)

print(plot_elbow_pc)

ggsave(
  filename = file.path(results_path, "figure_elbow_method_PCs.png"),
  plot = plot_elbow_pc,
  width = 8,
  height = 6,
  dpi = 300
)

cat("Saved:", file.path(results_path, "figure_elbow_method_PCs.png"), "\n")

# K-means on PCA scores
set.seed(123)

k_pc <- 3

kmeans_pc <- kmeans(
  pc_scores %>% select(-cd4_group),
  centers = k_pc,
  nstart = 25
)

plot_kmeans_pc <- factoextra::fviz_cluster(
  kmeans_pc,
  data = pc_scores %>% select(-cd4_group),
  geom = "point",
  ellipse.type = "convex",
  ggtheme = theme_bw()
) +
  labs(
    title = paste0("K-means clustering on PCA scores, k = ", k_pc),
    subtitle = paste0("PC1–PC", m_pcs)
  ) +
  base_theme_with_legend

print(plot_kmeans_pc)

ggsave(
  filename = file.path(results_path, "figure_kmeans_clusters_PCs.png"),
  plot = plot_kmeans_pc,
  width = 9,
  height = 7,
  dpi = 300
)

cat("Saved:", file.path(results_path, "figure_kmeans_clusters_PCs.png"), "\n")

# Cluster composition in PC space
cluster_composition_pc <- tibble(
  cd4_group = pc_scores$cd4_group,
  cluster = factor(kmeans_pc$cluster)
) %>%
  count(cluster, cd4_group) %>%
  group_by(cluster) %>%
  mutate(percent_within_cluster = 100 * n / sum(n)) %>%
  ungroup()

print(cluster_composition_pc)

write_xlsx(
  list("Cluster_composition_PC" = cluster_composition_pc),
  path = file.path(results_path, "table_cluster_composition_PC.xlsx")
)

cat("Saved:", file.path(results_path, "table_cluster_composition_PC.xlsx"), "\n")

# Hierarchical clustering on PCA scores
hclust_pc <- hclust(
  dist(pc_scores %>% select(-cd4_group)),
  method = "ward.D2"
)

plot_hclust_pc <- factoextra::fviz_dend(
  hclust_pc,
  k = k_pc,
  cex = 0.6,
  rect = TRUE,
  rect_fill = TRUE
) +
  labs(
    title = paste0("Hierarchical clustering on PCA scores, Ward.D2, k = ", k_pc)
  ) +
  theme_minimal(base_size = 14)

print(plot_hclust_pc)

ggsave(
  filename = file.path(results_path, "figure_hclust_dendrogram_PCs.png"),
  plot = plot_hclust_pc,
  width = 10,
  height = 6,
  dpi = 300
)

cat("Saved:", file.path(results_path, "figure_hclust_dendrogram_PCs.png"), "\n\n")


#### 12. SILHOUETTE VALIDATION OF CLUSTERS ####

cat("--- 12. Silhouette validation of clusters ---\n")

# Silhouette for original scaled variables
sil_original <- cluster::silhouette(
  kmeans_result$cluster,
  dist(pca_data_scaled)
)

sil_original_summary <- tibble(
  analysis = "K-means on original scaled immune markers",
  k = k_original,
  mean_silhouette = mean(sil_original[, "sil_width"]),
  median_silhouette = median(sil_original[, "sil_width"]),
  min_silhouette = min(sil_original[, "sil_width"]),
  max_silhouette = max(sil_original[, "sil_width"])
)

sil_original_by_cluster <- as.data.frame(sil_original) %>%
  group_by(cluster) %>%
  summarise(
    n = n(),
    mean_silhouette = mean(sil_width),
    median_silhouette = median(sil_width),
    min_silhouette = min(sil_width),
    max_silhouette = max(sil_width),
    .groups = "drop"
  )

# Silhouette for PCA-score clustering
sil_pc <- cluster::silhouette(
  kmeans_pc$cluster,
  dist(pc_scores %>% select(-cd4_group))
)

sil_pc_summary <- tibble(
  analysis = paste0("K-means on PCA scores PC1–PC", m_pcs),
  k = k_pc,
  mean_silhouette = mean(sil_pc[, "sil_width"]),
  median_silhouette = median(sil_pc[, "sil_width"]),
  min_silhouette = min(sil_pc[, "sil_width"]),
  max_silhouette = max(sil_pc[, "sil_width"])
)

sil_pc_by_cluster <- as.data.frame(sil_pc) %>%
  group_by(cluster) %>%
  summarise(
    n = n(),
    mean_silhouette = mean(sil_width),
    median_silhouette = median(sil_width),
    min_silhouette = min(sil_width),
    max_silhouette = max(sil_width),
    .groups = "drop"
  )

silhouette_summary_all <- bind_rows(
  sil_original_summary,
  sil_pc_summary
)

print(silhouette_summary_all)
print(sil_original_by_cluster)
print(sil_pc_by_cluster)

write_xlsx(
  list(
    "Silhouette_summary_all" = silhouette_summary_all,
    "Silhouette_original_by_cluster" = sil_original_by_cluster,
    "Silhouette_PC_by_cluster" = sil_pc_by_cluster,
    "Cluster_composition_original" = cluster_composition_original,
    "Cluster_composition_PC" = cluster_composition_pc
  ),
  path = file.path(results_path, "cluster_validation_silhouette.xlsx")
)

cat("Saved:", file.path(results_path, "cluster_validation_silhouette.xlsx"), "\n")

# Silhouette figure for original scaled variables
plot_silhouette_original <- factoextra::fviz_silhouette(sil_original) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Silhouette plot",
    subtitle = "K-means on original scaled immune markers"
  )

print(plot_silhouette_original)

ggsave(
  filename = file.path(results_path, "figure_silhouette_clusters_original.png"),
  plot = plot_silhouette_original,
  width = 8,
  height = 6,
  dpi = 300
)

cat("Saved:", file.path(results_path, "figure_silhouette_clusters_original.png"), "\n")

# Silhouette figure for PCA scores
plot_silhouette_pc <- factoextra::fviz_silhouette(sil_pc) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Silhouette plot",
    subtitle = paste0("K-means on PCA scores PC1–PC", m_pcs)
  )

print(plot_silhouette_pc)

ggsave(
  filename = file.path(results_path, "figure_silhouette_clusters_PC.png"),
  plot = plot_silhouette_pc,
  width = 8,
  height = 6,
  dpi = 300
)

cat("Saved:", file.path(results_path, "figure_silhouette_clusters_PC.png"), "\n\n")


#### 13. COMBINED PCA AND CLUSTERING PANELS ####

cat("--- 13. Combined PCA and clustering panels ---\n")

# Panel using original scaled variables
plot_cluster_panel_original <- ggarrange(
  plot_pca_individuals,
  plot_pca_biplot,
  plot_elbow,
  plot_kmeans,
  ncol = 2,
  nrow = 2,
  labels = c("A", "B", "C", "D")
)

print(plot_cluster_panel_original)

ggsave(
  filename = file.path(results_path, "figure_cluster_panel.png"),
  plot = plot_cluster_panel_original,
  width = 14,
  height = 12,
  dpi = 300
)

cat("Saved:", file.path(results_path, "figure_cluster_panel.png"), "\n")

# Panel using PCA scores
plot_cluster_panel_pc <- ggarrange(
  plot_pca_individuals,
  plot_pca_biplot,
  plot_elbow_pc,
  plot_kmeans_pc,
  ncol = 2,
  nrow = 2,
  labels = c("A", "B", "C", "D")
)

print(plot_cluster_panel_pc)

ggsave(
  filename = file.path(results_path, "figure_cluster_panel_PCs.png"),
  plot = plot_cluster_panel_pc,
  width = 14,
  height = 12,
  dpi = 300
)

cat("Saved:", file.path(results_path, "figure_cluster_panel_PCs.png"), "\n\n")

cat("=== END OF BLOCK 2 ===\n")
#### 14. MULTIVARIABLE LOGISTIC REGRESSION WITH LOG10-Z VARIABLES ####

cat("--- 14. Multivariable logistic regression with log10-z variables ---\n")

# Outcome:
# CD4-defined advanced immunodeficiency = CD4+ T-cell count <350 cells/mm³
#
# Predictors:
# Continuous soluble biomarkers are log10(x+1)-transformed and standardized.
# ORs therefore represent the change in odds per 1 standard deviation increase
# in each log-transformed predictor.

modelo_log_data2 <- dados_modelo2 %>%
  select(
    outcome_cd4_350,
    sexo,
    z_idade,
    z_Log,
    z_log_IL2,
    z_log_IL6,
    z_log_il2r,
    z_log_il6r,
    z_log_cd40l
  ) %>%
  drop_na()

cat("Logistic model dataset N =", nrow(modelo_log_data2), "\n")
cat("Outcome counts:\n")
print(table(modelo_log_data2$outcome_cd4_350))
cat("\n")


#### 14.1. UNIVARIATE LOGISTIC REGRESSION SCREENING ####

cat("--- 14.1. Univariate logistic regression screening ---\n")

candidate_predictors_logistic2 <- c(
  "z_idade",
  "sexo",
  "z_Log",
  "z_log_IL2",
  "z_log_IL6",
  "z_log_il2r",
  "z_log_il6r",
  "z_log_cd40l"
)

# ---- Function for univariate logistic regression
run_univ_logistic2 <- function(df, outcome, predictor) {
  
  f <- as.formula(paste0(outcome, " ~ ", predictor))
  
  m <- glm(
    f,
    data = df,
    family = binomial
  )
  
  broom::tidy(
    m,
    exponentiate = TRUE,
    conf.int = TRUE
  ) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      predictor = predictor,
      OR = estimate,
      CI_low = conf.low,
      CI_high = conf.high,
      p_value = p.value,
      p_value_formatted = format_p(p_value)
    ) %>%
    select(
      predictor,
      term,
      OR,
      CI_low,
      CI_high,
      p_value,
      p_value_formatted
    )
}

# Sex is not used for screening because it will be forced into the model
# Run univariable analyses for all predictors, including sex
univ_logistic2 <- purrr::map_dfr(
  candidate_predictors_logistic2,
  ~run_univ_logistic2(modelo_log_data2, "outcome_cd4_350", .x)
) %>%
  arrange(p_value)

print(univ_logistic2)

# ---- Variable selection using p <= 0.20
selected_predictors_logistic2 <- univ_logistic2 %>%
  group_by(predictor) %>%
  summarise(
    min_p = min(p_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(min_p <= 0.20) %>%
  pull(predictor)

# Force sex into the multivariable model as an adjustment covariate
selected_predictors_logistic2 <- unique(
  c("sexo", selected_predictors_logistic2)
)

cat("Selected predictors for final model:\n")
print(selected_predictors_logistic2)
cat("\n")

# Safety fallback
if (length(selected_predictors_logistic2) == 1) {
  selected_predictors_logistic2 <- c("sexo", "z_idade", "z_Log")
  cat("Only sex was selected. Using fallback model with sex, age, and viral load.\n")
}

final_formula_logistic2 <- as.formula(
  paste(
    "outcome_cd4_350 ~",
    paste(selected_predictors_logistic2, collapse = " + ")
  )
)

cat("Final logistic regression formula:\n")
print(final_formula_logistic2)
cat("\n")


#### 14.2. FINAL MULTIVARIABLE LOGISTIC REGRESSION MODEL ####

cat("--- 14.2. Final multivariable logistic regression model ---\n")

modelo_logistico2 <- glm(
  final_formula_logistic2,
  data = modelo_log_data2,
  family = binomial
)

cat("Model summary:\n")
print(summary(modelo_logistico2))
cat("\n")

# ---- Adjusted OR table
tabela_or_ajustada2 <- broom::tidy(
  modelo_logistico2,
  exponentiate = TRUE,
  conf.int = TRUE
) %>%
  filter(term != "(Intercept)") %>%
  transmute(
    term,
    AOR = estimate,
    CI_low = conf.low,
    CI_high = conf.high,
    p_value = p.value,
    p_value_formatted = format_p(p_value),
    AOR_CI95 = paste0(
      round(AOR, 2),
      " (",
      round(CI_low, 2),
      "–",
      round(CI_high, 2),
      ")"
    )
  )

print(tabela_or_ajustada2)

write_xlsx(
  list(
    "Univariate_logistic_zlog" = univ_logistic2,
    "Final_model_zlog" = tabela_or_ajustada2
  ),
  path = file.path(results_path, "tables_logistic_model_zlog.xlsx")
)

cat("Saved:", file.path(results_path, "tables_logistic_model_zlog.xlsx"), "\n\n")


#### 15. GENERAL MODEL STATISTICS ####

cat("--- 15. General statistics for final logistic model ---\n")

# ---- Null model
modelo_nulo2 <- glm(
  outcome_cd4_350 ~ 1,
  data = modelo_log_data2,
  family = binomial
)

# ---- Likelihood-ratio test
lr_test2 <- anova(
  modelo_nulo2,
  modelo_logistico2,
  test = "Chisq"
)

# ---- Pseudo-R2
pseudo_r2_2 <- pscl::pR2(modelo_logistico2)

# ---- Predicted probabilities
modelo_log_data2$pred_prob <- predict(
  modelo_logistico2,
  type = "response"
)

# ---- Hosmer-Lemeshow test
y_numeric2 <- ifelse(modelo_log_data2$outcome_cd4_350 == "Yes", 1, 0)

hl_test2 <- ResourceSelection::hoslem.test(
  x = y_numeric2,
  y = modelo_log_data2$pred_prob,
  g = 10
)

# ---- VIF
vif_values2 <- car::vif(modelo_logistico2)

vif_table2 <- tibble(
  variable = names(vif_values2),
  VIF = as.numeric(vif_values2)
)

# ---- Model statistics table
model_stats2 <- tibble(
  metric = c(
    "N",
    "AIC",
    "BIC",
    "LR Chi-square",
    "LR df",
    "LR p-value",
    "McFadden R2",
    "Nagelkerke/Cragg-Uhler R2",
    "Hosmer-Lemeshow Chi-square",
    "Hosmer-Lemeshow df",
    "Hosmer-Lemeshow p-value"
  ),
  value = c(
    nrow(modelo_log_data2),
    AIC(modelo_logistico2),
    BIC(modelo_logistico2),
    lr_test2$Deviance[2],
    lr_test2$Df[2],
    lr_test2$`Pr(>Chi)`[2],
    unname(pseudo_r2_2["McFadden"]),
    unname(pseudo_r2_2["r2CU"]),
    unname(hl_test2$statistic),
    unname(hl_test2$parameter),
    unname(hl_test2$p.value)
  )
) %>%
  mutate(
    value_formatted = case_when(
      metric == "LR p-value" ~ format_p(value),
      metric == "Hosmer-Lemeshow p-value" ~ format_p(value),
      TRUE ~ as.character(round(value, 4))
    )
  )

cat("Model statistics:\n")
print(model_stats2)
cat("\n")

cat("VIF values:\n")
print(vif_table2)
cat("\n")

write_xlsx(
  list(
    "Model_statistics" = model_stats2,
    "VIF" = vif_table2
  ),
  path = file.path(results_path, "model_statistics_logistic_zlog.xlsx")
)

cat("Saved:", file.path(results_path, "model_statistics_logistic_zlog.xlsx"), "\n\n")


#### 16. ROC CURVE AND AUC ####

cat("--- 16. ROC curve and AUC ---\n")

roc_model2 <- pROC::roc(
  response = modelo_log_data2$outcome_cd4_350,
  predictor = modelo_log_data2$pred_prob,
  levels = c("No", "Yes"),
  direction = "<",
  ci = TRUE
)

auc_value2 <- as.numeric(pROC::auc(roc_model2))
auc_ci2 <- as.numeric(pROC::ci.auc(roc_model2))

cat("AUC:", auc_value2, "\n")
cat("AUC 95% CI:", auc_ci2[1], "-", auc_ci2[3], "\n\n")

# ---- Best cutoff by Youden index
best_coords2 <- pROC::coords(
  roc_model2,
  x = "best",
  best.method = "youden",
  ret = c(
    "threshold",
    "sensitivity",
    "specificity",
    "ppv",
    "npv",
    "accuracy"
  ),
  transpose = FALSE
)

print(best_coords2)

youden_cutoff2 <- as.numeric(best_coords2$threshold)

# ---- Brier score
brier_score2 <- mean(
  (y_numeric2 - modelo_log_data2$pred_prob)^2
)

roc_summary2 <- tibble(
  AUC = auc_value2,
  AUC_CI_low = auc_ci2[1],
  AUC_CI_mid = auc_ci2[2],
  AUC_CI_high = auc_ci2[3],
  Youden_cutoff = youden_cutoff2,
  Brier_score = brier_score2
)

print(roc_summary2)

# ---- ROC figure
roc_plot2 <- ggroc(
  roc_model2,
  linewidth = 1.2
) +
  theme_minimal(base_size = 16) +
  labs(
    title = "ROC curve for the final multivariable logistic model",
    subtitle = paste0(
      "AUC = ",
      round(auc_value2, 3),
      " (95% CI: ",
      round(auc_ci2[1], 3),
      "–",
      round(auc_ci2[3], 3),
      ")"
    ),
    x = "1 - Specificity",
    y = "Sensitivity"
  )

print(roc_plot2)

ggsave(
  filename = file.path(results_path, "figure_ROC_logistic_model_zlog.png"),
  plot = roc_plot2,
  width = 8,
  height = 6,
  dpi = 300
)

cat("Saved:", file.path(results_path, "figure_ROC_logistic_model_zlog.png"), "\n\n")


#### 17. CONFUSION MATRICES AND CLASSIFICATION METRICS ####

cat("--- 17. Confusion matrices and classification metrics ---\n")

get_confusion_metrics <- function(data, cutoff) {
  
  pred_class <- factor(
    ifelse(data$pred_prob >= cutoff, "Yes", "No"),
    levels = c("No", "Yes")
  )
  
  cm <- caret::confusionMatrix(
    data = pred_class,
    reference = data$outcome_cd4_350,
    positive = "Yes"
  )
  
  tab <- as.data.frame(cm$table)
  
  metrics <- tibble(
    cutoff = cutoff,
    TN = cm$table["No", "No"],
    FP = cm$table["Yes", "No"],
    FN = cm$table["No", "Yes"],
    TP = cm$table["Yes", "Yes"],
    accuracy = unname(cm$overall["Accuracy"]),
    sensitivity = unname(cm$byClass["Sensitivity"]),
    specificity = unname(cm$byClass["Specificity"]),
    PPV = unname(cm$byClass["Pos Pred Value"]),
    NPV = unname(cm$byClass["Neg Pred Value"]),
    balanced_accuracy = unname(cm$byClass["Balanced Accuracy"])
  )
  
  list(
    confusion_table = tab,
    metrics = metrics,
    cm_object = cm
  )
}

# ---- Metrics using conventional cutoff = 0.5 ----
cm_05 <- get_confusion_metrics(
  data = modelo_log_data2,
  cutoff = 0.5
)

cat("Metrics using cutoff = 0.5:\n")
print(cm_05$metrics)
cat("\n")

# ---- Metrics using Youden cutoff ----
cm_youden <- get_confusion_metrics(
  data = modelo_log_data2,
  cutoff = youden_cutoff2
)

cat("Metrics using Youden cutoff:\n")
print(cm_youden$metrics)
cat("\n")

# ---- Export ROC and confusion outputs ----
write_xlsx(
  list(
    "ROC_summary" = roc_summary2,
    "Best_cutoff_Youden" = as.data.frame(best_coords2),
    "Metrics_cutoff_0.5" = cm_05$metrics,
    "Confusion_cutoff_0.5" = cm_05$confusion_table,
    "Metrics_cutoff_Youden" = cm_youden$metrics,
    "Confusion_cutoff_Youden" = cm_youden$confusion_table
  ),
  path = file.path(results_path, "roc_confusion_logistic_zlog.xlsx")
)

cat("Saved:", file.path(results_path, "roc_confusion_logistic_zlog.xlsx"), "\n\n")


#### 18. INTERNAL BOOTSTRAP VALIDATION FOR AUC ####

cat("--- 18. Internal bootstrap validation for AUC ---\n")

set.seed(123)

B <- 1000

formula_boot <- formula(modelo_logistico2)

boot_auc <- replicate(B, {
  
  idx <- sample(
    seq_len(nrow(modelo_log_data2)),
    replace = TRUE
  )
  
  boot_data <- modelo_log_data2[idx, ]
  
  # Skip samples with only one outcome class
  if (length(unique(boot_data$outcome_cd4_350)) < 2) {
    return(NA_real_)
  }
  
  fit_boot <- tryCatch(
    glm(
      formula_boot,
      data = boot_data,
      family = binomial
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit_boot)) {
    return(NA_real_)
  }
  
  pred_boot <- tryCatch(
    predict(
      fit_boot,
      newdata = modelo_log_data2,
      type = "response"
    ),
    error = function(e) rep(NA_real_, nrow(modelo_log_data2))
  )
  
  if (all(is.na(pred_boot))) {
    return(NA_real_)
  }
  
  roc_boot <- tryCatch(
    pROC::roc(
      response = modelo_log_data2$outcome_cd4_350,
      predictor = pred_boot,
      levels = c("No", "Yes"),
      direction = "<",
      quiet = TRUE
    ),
    error = function(e) NULL
  )
  
  if (is.null(roc_boot)) {
    return(NA_real_)
  }
  
  as.numeric(pROC::auc(roc_boot))
})

boot_auc <- boot_auc[!is.na(boot_auc)]

bootstrap_auc_summary <- tibble(
  n_bootstrap_valid = length(boot_auc),
  AUC_apparent = auc_value2,
  AUC_bootstrap_mean = mean(boot_auc),
  AUC_bootstrap_median = median(boot_auc),
  AUC_bootstrap_CI_low = as.numeric(quantile(boot_auc, 0.025)),
  AUC_bootstrap_CI_high = as.numeric(quantile(boot_auc, 0.975)),
  optimism = auc_value2 - mean(boot_auc)
)

cat("Bootstrap AUC summary:\n")
print(bootstrap_auc_summary, width = Inf)
cat("\n")

write_xlsx(
  list(
    "Bootstrap_AUC_summary" = bootstrap_auc_summary,
    "Bootstrap_AUC_values" = tibble(AUC = boot_auc)
  ),
  path = file.path(results_path, "bootstrap_auc_logistic_zlog.xlsx")
)

cat("Saved:", file.path(results_path, "bootstrap_auc_logistic_zlog.xlsx"), "\n\n")


#### 19. OPTIONAL: COMPACT FINAL MODEL TABLE FOR MANUSCRIPT ####

cat("--- 19. Compact final model table for manuscript ---\n")

# This table combines the main adjusted ORs and model performance metrics.

final_model_table_for_article <- tabela_or_ajustada2 %>%
  mutate(
    variable_label = case_when(
      term == "sexoM" ~ "Sex (male)",
      term == "z_Log" ~ "Viral load, log10 copies/mL, z-score",
      term == "z_idade" ~ "Age, z-score",
      term == "z_log_IL6" ~ "sIL-6, log10(x+1), z-score",
      term == "z_log_il6r" ~ "sIL-6R, log10(x+1), z-score",
      term == "z_log_IL2" ~ "sIL-2, log10(x+1), z-score",
      term == "z_log_il2r" ~ "sIL-2R, log10(x+1), z-score",
      term == "z_log_cd40l" ~ "sCD40L, log10(x+1), z-score",
      TRUE ~ term
    )
  ) %>%
  select(
    variable_label,
    term,
    AOR,
    CI_low,
    CI_high,
    AOR_CI95,
    p_value,
    p_value_formatted
  )

model_performance_table_for_article <- tibble(
  metric = c(
    "AUC",
    "AUC 95% CI lower",
    "AUC 95% CI upper",
    "Bootstrap mean AUC",
    "Bootstrap AUC 95% CI lower",
    "Bootstrap AUC 95% CI upper",
    "Brier score",
    "Youden cutoff",
    "Accuracy at Youden cutoff",
    "Sensitivity at Youden cutoff",
    "Specificity at Youden cutoff",
    "PPV at Youden cutoff",
    "NPV at Youden cutoff",
    "McFadden R2",
    "Nagelkerke/Cragg-Uhler R2",
    "Hosmer-Lemeshow p-value"
  ),
  value = c(
    auc_value2,
    auc_ci2[1],
    auc_ci2[3],
    bootstrap_auc_summary$AUC_bootstrap_mean,
    bootstrap_auc_summary$AUC_bootstrap_CI_low,
    bootstrap_auc_summary$AUC_bootstrap_CI_high,
    brier_score2,
    youden_cutoff2,
    cm_youden$metrics$accuracy,
    cm_youden$metrics$sensitivity,
    cm_youden$metrics$specificity,
    cm_youden$metrics$PPV,
    cm_youden$metrics$NPV,
    unname(pseudo_r2_2["McFadden"]),
    unname(pseudo_r2_2["r2CU"]),
    unname(hl_test2$p.value)
  )
) %>%
  mutate(
    value_formatted = round(value, 4)
  )

print(final_model_table_for_article)
print(model_performance_table_for_article)

write_xlsx(
  list(
    "Final_model_ORs_for_article" = final_model_table_for_article,
    "Model_performance_for_article" = model_performance_table_for_article
  ),
  path = file.path(results_path, "final_logistic_model_tables_for_article.xlsx")
)

cat("Saved:", file.path(results_path, "final_logistic_model_tables_for_article.xlsx"), "\n\n")

cat("=== END OF BLOCK 3 ===\n")
#### 20. TREATMENT ANALYSIS: PAIRED BT VS AT AT 2 AND 4 MONTHS ####

cat("--- 20. Treatment analysis: paired BT vs AT at 2 and 4 months ---\n")

# This section evaluates paired changes before and after ART.
# BT = before treatment
# AT = after treatment
# Separate analyses are performed for 2-month and 4-month paired subgroups.


#### 20.1. READ AND CLEAN TREATMENT DATASET ####

cat("--- 20.1. Reading and cleaning treatment dataset ---\n")

# Read all columns as text first to avoid problems with mixed decimal formats
hiv_tratamento_wide_raw <- readxl::read_excel(
  data_path_treatment,
  col_types = "text"
)

# Create Patient_ID and parse numeric columns robustly
hiv_tratamento_wide <- hiv_tratamento_wide_raw %>%
  mutate(
    Patient_ID = row_number(),
    .before = 1
  ) %>%
  mutate(
    across(
      -Patient_ID,
      parse_num_robust
    )
  )

cat("Treatment dataset dimensions:\n")
print(dim(hiv_tratamento_wide))

cat("Treatment dataset column names:\n")
print(names(hiv_tratamento_wide))
cat("\n")


#### 20.2. IDENTIFY 2-MONTH AND 4-MONTH COLUMNS ####

cat("--- 20.2. Identifying 2-month and 4-month columns ---\n")

cols_2m <- stringr::str_subset(
  names(hiv_tratamento_wide),
  "(bt|at)2m$"
)

cols_4m <- stringr::str_subset(
  names(hiv_tratamento_wide),
  "(bt|at)4m$"
)

cat("2-month columns:\n")
print(cols_2m)

cat("4-month columns:\n")
print(cols_4m)
cat("\n")


#### 20.3. CONVERT WIDE TREATMENT DATASET TO LONG FORMAT ####

cat("--- 20.3. Converting treatment dataset to long format ---\n")

make_long_time <- function(df_wide, cols_time, suffix) {
  
  df_wide %>%
    select(
      Patient_ID,
      all_of(cols_time)
    ) %>%
    pivot_longer(
      cols = -Patient_ID,
      names_to = "var",
      values_to = "Dosage",
      values_drop_na = FALSE
    ) %>%
    mutate(
      Biomarker = stringr::str_replace(
        var,
        paste0("_(bt|at)", suffix, "$"),
        ""
      ),
      Timecode = stringr::str_extract(
        var,
        paste0("(bt|at)", suffix, "$")
      ),
      Timepoint = factor(
        case_when(
          stringr::str_starts(Timecode, "bt") ~ "Before (BT)",
          stringr::str_starts(Timecode, "at") ~ "After (AT)",
          TRUE ~ NA_character_
        ),
        levels = c("Before (BT)", "After (AT)")
      )
    ) %>%
    select(
      -var,
      -Timecode
    ) %>%
    filter(
      !is.na(Timepoint)
    )
}

dados_2m <- make_long_time(
  df_wide = hiv_tratamento_wide,
  cols_time = cols_2m,
  suffix = "2m"
)

dados_4m <- make_long_time(
  df_wide = hiv_tratamento_wide,
  cols_time = cols_4m,
  suffix = "4m"
)

cat("Available biomarkers at 2 months:\n")
print(sort(unique(dados_2m$Biomarker)))

cat("Available biomarkers at 4 months:\n")
print(sort(unique(dados_4m$Biomarker)))
cat("\n")


#### 20.4. SANITY CHECKS FOR TREATMENT DATA ####

cat("--- 20.4. Sanity checks for treatment data ---\n")

# This check verifies if BT and AT values are accidentally identical for sCD40L.
# This is useful because spreadsheet copy/paste errors may occur in paired datasets.

if (all(c("cd40l_bt2m", "cd40l_at2m") %in% names(hiv_tratamento_wide))) {
  
  chk_cd40l_2m <- hiv_tratamento_wide %>%
    transmute(
      bt = cd40l_bt2m,
      at = cd40l_at2m,
      ok_pair = !is.na(bt) & !is.na(at),
      equal = ok_pair & (bt == at)
    ) %>%
    summarise(
      n_pairs = sum(ok_pair),
      n_equal = sum(equal),
      prop_equal = ifelse(n_pairs > 0, n_equal / n_pairs, NA_real_)
    )
  
  cat("Sanity check for sCD40L at 2 months:\n")
  print(chk_cd40l_2m)
}

if (all(c("cd40l_bt4m", "cd40l_at4m") %in% names(hiv_tratamento_wide))) {
  
  chk_cd40l_4m <- hiv_tratamento_wide %>%
    transmute(
      bt = cd40l_bt4m,
      at = cd40l_at4m,
      ok_pair = !is.na(bt) & !is.na(at),
      equal = ok_pair & (bt == at)
    ) %>%
    summarise(
      n_pairs = sum(ok_pair),
      n_equal = sum(equal),
      prop_equal = ifelse(n_pairs > 0, n_equal / n_pairs, NA_real_)
    )
  
  cat("Sanity check for sCD40L at 4 months:\n")
  print(chk_cd40l_4m)
}

cat("\n")


#### 20.5. PAIRED WILCOXON TEST FUNCTION ####

cat("--- 20.5. Defining paired Wilcoxon test function ---\n")

run_single_paired_test_no_impute <- function(data, marker_name) {
  
  data_filtered <- data %>%
    filter(
      Biomarker == marker_name
    )
  
  if (nrow(data_filtered) == 0) {
    return(
      tibble(
        Biomarker = marker_name,
        group1 = "Before (BT)",
        group2 = "After (AT)",
        n_pairs = NA_integer_,
        statistic = NA_real_,
        p = NA_real_,
        p_formatted = NA_character_,
        p_signif = "NA",
        median_before = NA_real_,
        median_after = NA_real_,
        median_difference = NA_real_
      )
    )
  }
  
  wide <- data_filtered %>%
    pivot_wider(
      id_cols = Patient_ID,
      names_from = Timepoint,
      values_from = Dosage
    )
  
  if (!all(c("Before (BT)", "After (AT)") %in% names(wide))) {
    return(
      tibble(
        Biomarker = marker_name,
        group1 = "Before (BT)",
        group2 = "After (AT)",
        n_pairs = NA_integer_,
        statistic = NA_real_,
        p = NA_real_,
        p_formatted = NA_character_,
        p_signif = "NA",
        median_before = NA_real_,
        median_after = NA_real_,
        median_difference = NA_real_
      )
    )
  }
  
  wide_complete <- wide %>%
    drop_na(
      `Before (BT)`,
      `After (AT)`
    )
  
  n_pairs <- nrow(wide_complete)
  
  if (n_pairs < 2) {
    return(
      tibble(
        Biomarker = marker_name,
        group1 = "Before (BT)",
        group2 = "After (AT)",
        n_pairs = n_pairs,
        statistic = NA_real_,
        p = NA_real_,
        p_formatted = NA_character_,
        p_signif = "NA",
        median_before = ifelse(n_pairs > 0, median(wide_complete$`Before (BT)`, na.rm = TRUE), NA_real_),
        median_after = ifelse(n_pairs > 0, median(wide_complete$`After (AT)`, na.rm = TRUE), NA_real_),
        median_difference = ifelse(
          n_pairs > 0,
          median(wide_complete$`After (AT)` - wide_complete$`Before (BT)`, na.rm = TRUE),
          NA_real_
        )
      )
    )
  }
  
  diffs <- wide_complete$`After (AT)` - wide_complete$`Before (BT)`
  
  median_before <- median(wide_complete$`Before (BT)`, na.rm = TRUE)
  median_after <- median(wide_complete$`After (AT)`, na.rm = TRUE)
  median_difference <- median(diffs, na.rm = TRUE)
  
  # If all differences are exactly zero, Wilcoxon may return NA.
  # In this case, p is set to 1.
  if (all(diffs == 0)) {
    return(
      tibble(
        Biomarker = marker_name,
        group1 = "Before (BT)",
        group2 = "After (AT)",
        n_pairs = n_pairs,
        statistic = 0,
        p = 1,
        p_formatted = "1.0000",
        p_signif = "ns",
        median_before = median_before,
        median_after = median_after,
        median_difference = median_difference
      )
    )
  }
  
  wt <- wilcox.test(
    wide_complete$`Before (BT)`,
    wide_complete$`After (AT)`,
    paired = TRUE,
    exact = FALSE
  )
  
  tibble(
    Biomarker = marker_name,
    group1 = "Before (BT)",
    group2 = "After (AT)",
    n_pairs = n_pairs,
    statistic = unname(wt$statistic),
    p = wt$p.value,
    p_formatted = format_p(wt$p.value),
    p_signif = p_to_signif(wt$p.value),
    median_before = median_before,
    median_after = median_after,
    median_difference = median_difference
  )
}


#### 20.6. RUN PAIRED TESTS FOR ALL AVAILABLE TREATMENT MARKERS ####

cat("--- 20.6. Running paired tests for treatment markers ---\n")

available_markers <- sort(
  unique(
    c(
      unique(dados_2m$Biomarker),
      unique(dados_4m$Biomarker)
    )
  )
)

cat("All available treatment markers:\n")
print(available_markers)

# Markers to prioritize in manuscript figure
markers_treatment_priority <- c(
  "il2r",
  "il6r",
  "il2",
  "il6",
  "cd40l",
  "cd4",
  "cv"
)

markers_treatment <- intersect(
  markers_treatment_priority,
  available_markers
)

# Add any extra available markers at the end
markers_treatment <- unique(
  c(
    markers_treatment,
    setdiff(available_markers, markers_treatment)
  )
)

cat("Markers included in paired treatment analysis:\n")
print(markers_treatment)

treatment_tests_2m <- purrr::map_dfr(
  markers_treatment,
  ~run_single_paired_test_no_impute(dados_2m, .x)
) %>%
  mutate(
    Time = "2 Months",
    .before = 1
  )

treatment_tests_4m <- purrr::map_dfr(
  markers_treatment,
  ~run_single_paired_test_no_impute(dados_4m, .x)
) %>%
  mutate(
    Time = "4 Months",
    .before = 1
  )

treatment_tests_all <- bind_rows(
  treatment_tests_2m,
  treatment_tests_4m
) %>%
  mutate(
    Biomarker_label = case_when(
      Biomarker == "il2r" ~ "sIL-2R",
      Biomarker == "il6r" ~ "sIL-6R",
      Biomarker == "il2" ~ "sIL-2",
      Biomarker == "il6" ~ "sIL-6",
      Biomarker == "cd40l" ~ "sCD40L",
      Biomarker == "cd4" ~ "CD4+ T cells",
      Biomarker == "cv" ~ "Viral load",
      TRUE ~ Biomarker
    )
  ) %>%
  select(
    Time,
    Biomarker,
    Biomarker_label,
    n_pairs,
    median_before,
    median_after,
    median_difference,
    statistic,
    p,
    p_formatted,
    p_signif
  )

print(treatment_tests_all)

write_xlsx(
  list(
    "Treatment_paired_tests_all" = treatment_tests_all,
    "Treatment_paired_tests_2m" = treatment_tests_2m,
    "Treatment_paired_tests_4m" = treatment_tests_4m
  ),
  path = file.path(results_path, "table_treatment_paired_wilcoxon.xlsx")
)

cat("Saved:", file.path(results_path, "table_treatment_paired_wilcoxon.xlsx"), "\n\n")


#### 20.7. TREATMENT PLOT FUNCTION ####

cat("--- 20.7. Defining treatment plot function ---\n")

biomarker_labels_treatment <- c(
  "il2r"  = "sIL-2R (pg/mL)",
  "il6r"  = "sIL-6R (pg/mL)",
  "cd40l" = "sCD40L (pg/mL)",
  "il2"   = "sIL-2 (pg/mL)",
  "il6"   = "sIL-6 (pg/mL)",
  "cd4"   = "CD4+ T cells (cells/mm³)",
  "cv"    = "Viral load (copies/mL)"
)

create_facet_plot_treatment <- function(data_2m, data_4m, marker,
                                        transform_y = c("none", "log10_plus1"),
                                        y_label = NULL,
                                        title = NULL) {
  
  transform_y <- match.arg(transform_y)
  
  data_plot <- bind_rows(
    data_2m %>%
      filter(Biomarker == marker) %>%
      mutate(Time = "2 Months"),
    data_4m %>%
      filter(Biomarker == marker) %>%
      mutate(Time = "4 Months")
  )
  
  if (nrow(data_plot) == 0) {
    warning(paste0("Marker '", marker, "' not found in treatment dataset. Returning NULL."))
    return(NULL)
  }
  
  # Create plotting variable
  data_plot <- data_plot %>%
    mutate(
      PlotValue = case_when(
        transform_y == "log10_plus1" ~ log10(Dosage + 1),
        TRUE ~ Dosage
      )
    )
  
  data_complete <- data_plot %>%
    drop_na(PlotValue)
  
  test_results <- bind_rows(
    run_single_paired_test_no_impute(data_2m, marker) %>%
      mutate(Time = "2 Months"),
    run_single_paired_test_no_impute(data_4m, marker) %>%
      mutate(Time = "4 Months")
  ) %>%
    mutate(
      group1 = "Before (BT)",
      group2 = "After (AT)",
      p.signif = p_signif
    )
  
  # Position of p-value labels based on transformed plotting scale
  y_max_by_time <- data_plot %>%
    group_by(Time) %>%
    summarise(
      y_max = max(PlotValue, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      y_max = ifelse(is.finite(y_max), y_max, 1),
      y.position = y_max + 0.12 * abs(y_max)
    )
  
  # If values are very compressed or y_max is close to zero
  y_max_by_time <- y_max_by_time %>%
    mutate(
      y.position = ifelse(y.position <= 0, y_max + 0.5, y.position)
    )
  
  test_results <- test_results %>%
    left_join(y_max_by_time %>% select(Time, y.position), by = "Time")
  
  if (is.null(y_label)) {
    y_label <- marker
  }
  
  if (transform_y == "log10_plus1") {
    y_label <- paste0("log10(", y_label, " + 1)")
  }
  
  if (is.null(title)) {
    title <- paste0(marker, " before vs after ART")
  }
  
  p <- ggplot(
    data_plot,
    aes(
      x = Timepoint,
      y = PlotValue,
      fill = Timepoint
    )
  ) +
    geom_boxplot(
      outlier.shape = NA,
      alpha = 0.8,
      na.rm = TRUE
    ) +
    geom_line(
      data = data_complete,
      aes(group = Patient_ID),
      color = "gray50",
      alpha = 0.5,
      position = position_dodge(width = 0.2)
    ) +
    geom_point(
      data = data_complete,
      aes(group = Patient_ID),
      position = position_dodge(width = 0.2),
      size = 2,
      alpha = 0.8
    ) +
    ggpubr::stat_pvalue_manual(
      test_results,
      label = "p.signif",
      y.position = "y.position",
      tip.length = 0.01,
      size = 5
    ) +
    facet_wrap(~Time, scales = "free_y") +
    scale_y_continuous(
      expand = expansion(mult = c(0.08, 0.22))
    ) +
    labs(
      title = title,
      x = "Treatment status",
      y = y_label
    ) +
    base_theme
  
  return(p)
}

#### 20.8. CREATE LONGITUDINAL TREATMENT FIGURE ####

cat("--- 20.8. Creating longitudinal treatment figure ---\n")

# For highly skewed markers, we plot log10(x+1) directly as the y variable.
# This avoids visual compression caused by scale_y_log10.

marker_y_transform <- c(
  "il2r" = "log10_plus1",
  "il6r" = "log10_plus1",
  "cd40l" = "log10_plus1",
  "il2" = "none",
  "il6" = "none",
  "cd4" = "none",
  "cv" = "log10_plus1"
)

p_treat_il2r <- if ("il2r" %in% available_markers) {
  create_facet_plot_treatment(
    dados_2m,
    dados_4m,
    marker = "il2r",
    transform_y = marker_y_transform["il2r"],
    y_label = biomarker_labels_treatment["il2r"],
    title = "sIL-2R before vs after ART"
  )
} else {
  NULL
}

p_treat_il6r <- if ("il6r" %in% available_markers) {
  create_facet_plot_treatment(
    dados_2m,
    dados_4m,
    marker = "il6r",
    transform_y = marker_y_transform["il6r"],
    y_label = biomarker_labels_treatment["il6r"],
    title = "sIL-6R before vs after ART"
  )
} else {
  NULL
}

p_treat_il2 <- if ("il2" %in% available_markers) {
  create_facet_plot_treatment(
    dados_2m,
    dados_4m,
    marker = "il2",
    transform_y = marker_y_transform["il2"],
    y_label = biomarker_labels_treatment["il2"],
    title = "sIL-2 before vs after ART"
  )
} else {
  NULL
}

p_treat_il6 <- if ("il6" %in% available_markers) {
  create_facet_plot_treatment(
    dados_2m,
    dados_4m,
    marker = "il6",
    transform_y = marker_y_transform["il6"],
    y_label = biomarker_labels_treatment["il6"],
    title = "sIL-6 before vs after ART"
  )
} else {
  NULL
}

p_treat_cd40l <- if ("cd40l" %in% available_markers) {
  create_facet_plot_treatment(
    dados_2m,
    dados_4m,
    marker = "cd40l",
    transform_y = marker_y_transform["cd40l"],
    y_label = biomarker_labels_treatment["cd40l"],
    title = "sCD40L before vs after ART"
  )
} else {
  NULL
}

p_treat_cd4 <- if ("cd4" %in% available_markers) {
  create_facet_plot_treatment(
    dados_2m,
    dados_4m,
    marker = "cd4",
    transform_y = marker_y_transform["cd4"],
    y_label = biomarker_labels_treatment["cd4"],
    title = "CD4+ T cells before vs after ART"
  )
} else {
  NULL
}

p_treat_cv <- if ("cv" %in% available_markers) {
  create_facet_plot_treatment(
    dados_2m,
    dados_4m,
    marker = "cv",
    transform_y = marker_y_transform["cv"],
    y_label = biomarker_labels_treatment["cv"],
    title = "Viral load before vs after ART"
  )
} else {
  NULL
}

plot_list_treatment <- list(
  p_treat_il2r,
  p_treat_il6r,
  p_treat_il2,
  p_treat_il6,
  p_treat_cd40l,
  p_treat_cd4,
  p_treat_cv
)

plot_list_treatment <- plot_list_treatment[
  !sapply(plot_list_treatment, is.null)
]

final_plots_treatment <- patchwork::wrap_plots(
  plot_list_treatment,
  ncol = 2
) +
  patchwork::plot_annotation(
    tag_levels = "A"
  )

print(final_plots_treatment)

ggsave(
  filename = file.path(results_path, "figure_treatment_panel.png"),
  plot = final_plots_treatment,
  width = 16,
  height = 18,
  dpi = 300
)

cat("Saved:", file.path(results_path, "figure_treatment_panel.png"), "\n\n")
#### 21. EXPORT CONSOLIDATED RESULTS ####

cat("--- 21. Exporting consolidated results ---\n")

# This file brings together the main statistical outputs into a single Excel file.
# It is useful for manuscript writing and checking all results in one place.

# Some objects may not exist if previous blocks were not run.
# The following checks avoid script interruption.

sheets_to_export <- list()

if (exists("tabela1_final")) {
  sheets_to_export[["Table1_baseline"]] <- tabela1_final
}

if (exists("kw_global_biomarkers")) {
  sheets_to_export[["KW_global_BH_biomarkers"]] <- kw_global_biomarkers
}

if (exists("posthoc_biomarkers")) {
  sheets_to_export[["Posthoc_Wilcoxon_BH"]] <- posthoc_biomarkers
}

if (exists("pca_eigenvalues")) {
  sheets_to_export[["PCA_eigenvalues"]] <- as.data.frame(pca_eigenvalues)
}

if (exists("pca_variable_contributions")) {
  sheets_to_export[["PCA_variable_contrib"]] <- pca_variable_contributions
}

if (exists("silhouette_summary_all")) {
  sheets_to_export[["Silhouette_summary"]] <- silhouette_summary_all
}

if (exists("cluster_composition_original")) {
  sheets_to_export[["Cluster_composition_original"]] <- cluster_composition_original
}

if (exists("cluster_composition_pc")) {
  sheets_to_export[["Cluster_composition_PC"]] <- cluster_composition_pc
}

if (exists("univ_logistic2")) {
  sheets_to_export[["Univariate_logistic_zlog"]] <- univ_logistic2
}

if (exists("tabela_or_ajustada2")) {
  sheets_to_export[["Final_logistic_model_zlog"]] <- tabela_or_ajustada2
}

if (exists("model_stats2")) {
  sheets_to_export[["Model_statistics"]] <- model_stats2
}

if (exists("vif_table2")) {
  sheets_to_export[["VIF"]] <- vif_table2
}

if (exists("roc_summary2")) {
  sheets_to_export[["ROC_summary"]] <- roc_summary2
}

if (exists("cm_05")) {
  sheets_to_export[["Metrics_cutoff_0.5"]] <- cm_05$metrics
  sheets_to_export[["Confusion_cutoff_0.5"]] <- cm_05$confusion_table
}

if (exists("cm_youden")) {
  sheets_to_export[["Metrics_cutoff_Youden"]] <- cm_youden$metrics
  sheets_to_export[["Confusion_cutoff_Youden"]] <- cm_youden$confusion_table
}

if (exists("bootstrap_auc_summary")) {
  sheets_to_export[["Bootstrap_AUC_summary"]] <- bootstrap_auc_summary
}

if (exists("treatment_tests_all")) {
  sheets_to_export[["Treatment_paired_tests"]] <- treatment_tests_all
}

if (length(sheets_to_export) > 0) {
  
  write_xlsx(
    sheets_to_export,
    path = file.path(results_path, "ALL_MAIN_STATISTICAL_OUTPUTS.xlsx")
  )
  
  cat("Saved:", file.path(results_path, "ALL_MAIN_STATISTICAL_OUTPUTS.xlsx"), "\n")
  
} else {
  
  cat("No objects available to export in consolidated file.\n")
}

cat("\n")


#### 22. SESSION INFO ####

cat("--- 22. Session info ---\n")

session_information <- capture.output(sessionInfo())

writeLines(
  session_information,
  con = file.path(results_path, "sessionInfo.txt")
)

cat("Saved:", file.path(results_path, "sessionInfo.txt"), "\n\n")


#### 23. FINAL MESSAGE ####

cat("=== SCRIPT FINISHED SUCCESSFULLY ===\n")
cat("All results were saved in:", results_path, "\n")
cat("Main outputs include:\n")
cat("- table1_baseline_characteristics.xlsx\n")
cat("- table_kw_global_BH_biomarkers.xlsx\n")
cat("- table_posthoc_wilcoxon_BH.xlsx\n")
cat("- tables_pca_variables.xlsx\n")
cat("- cluster_validation_silhouette.xlsx\n")
cat("- tables_logistic_model_zlog.xlsx\n")
cat("- model_statistics_logistic_zlog.xlsx\n")
cat("- roc_confusion_logistic_zlog.xlsx\n")
cat("- bootstrap_auc_logistic_zlog.xlsx\n")
cat("- table_treatment_paired_wilcoxon.xlsx\n")
cat("- ALL_MAIN_STATISTICAL_OUTPUTS.xlsx\n")
cat("- figure_cd4_groups_p1_to_p6.png\n")
cat("- figure_cluster_panel.png\n")
cat("- figure_cluster_panel_PCs.png\n")
cat("- figure_ROC_logistic_model_zlog.png\n")
cat("- figure_treatment_panel.png\n")
cat("=== END OF BLOCK 4 ===\n")
