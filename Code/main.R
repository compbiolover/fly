# Name: main.R
# Author: Andrew Willems <awillems@vols.utk.edu>
# Purpose: Main housing of all analysis for lethality_counts_072223.xlsx and Caffeine Lethality Test - Starvation.xlsx

# ----Loading packages----
pacman::p_load(
  ggpp, # Package for nice tables in our plots
  gplots, # Package to make heatmap.2 heatmaps
  HeatmapR, # Another package to make heatmaps
  patchwork, # Package to help with arranging ggplots
  readxl, # Package for reading Excel files
  survival, # Package for survival analysis
  survminer, # Package for Kaplan–Meier plots
  tidyverse # Package for data manipulation
)

# ----Loading raw data----
no_caf_df <- read_xlsx(path = "Data/lethality_counts_072223.xlsx", sheet = "Zero Caffeine")
mid_caf_df <- read_xlsx(path = "Data/lethality_counts_072223.xlsx", sheet = "7.5 mM Caffeine")
high_caf_df <- read_xlsx(path = "Data/lethality_counts_072223.xlsx", sheet = "15 mM Caffeine")


# ----Loading processed data (made from combining processed No Caf, Mid Caf, and High Caf dataframes)----
all_data <- read.csv("Data/all_data_reps_consolidated.csv", row.names = 1)
all_conds_75 <- read.csv("Data/all_conds_75.csv", row.names = 1)
all_conds_15 <- read.csv("Data/all_conds_15.csv", row.names = 1)

# ----Loading processed No Caf data----
no_caf_df <- read.csv("Data/low_caff_df.csv", sep = ",")
no_caf_df$genotype <- gsub("\\s+", "", no_caf_df$genotype)


# ----Survival analysis for Zero Caf Data----
no_caf_mod <- survfit(Surv(time, event) ~ genotype, data = no_caf_df)

# Log-rank test can't generate a p-value (as expected for 0 caffeine data)
# because their is no variation in survival
survdiff(Surv(time, event) ~ genotype, no_caf_df)

# ----Plotting KM curves for Zero Caf Data----
p <- ggsurvplot(
  fit = no_caf_mod, legend.title = "Genotype",
  legend.labs = c(
    "pdf/Df-1",
    "pdf/Df-2",
    "pdf/Df-3",
    "pdf/Df-4",
    "pdf01-1",
    "pdf01-2",
    "pdf01-3",
    "pdf01-4",
    "w1118-1",
    "w1118-2"
  ),
  font.main = c(18, "bold"),
  font.x = c(16, "bold"),
  font.y = c(16, "bold"),
  font.tickslab = c(14, "plain"),
  legend = "bottom", pval = TRUE,
  pval.method = TRUE, ggtheme = theme(
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    legend.key = element_blank()
  )
)

p <- p + ggtitle("0 mM Caffeine")
p


# ----Saving 0 Caf Plot----
ggsave(
  plot = p$plot, filename = "zero_caf_km.png",
  device = "png", path = "Outputs/Plots/",
  width = 8, height = 8, units = "in", dpi = 600
)

ggsave(
  plot = p$plot, filename = "zero_caf_km.svg",
  device = "svg", path = "Outputs/Plots/",
  width = 8, height = 8, units = "in", dpi = 600
)


# ----Loading processed Mid Caf data----
mid_caf_df <- read.csv("Data/mid_caff_df.csv", sep = ",")
mid_caf_df$genotype <- gsub("\\s+", "", mid_caf_df$genotype)

# ----Survival analysis for Mid Caf Data----
mid_caf_mod <- survfit(Surv(time = time, event = event) ~ genotype, data = mid_caf_df)

# ----Mid Caf Genotype Log-rank test. P-value: 4e-07----
mid_genotype_lrt <- survdiff(Surv(time, event) ~ genotype, mid_caf_df)
mid_genotype_pairwise <- pairwise_survdiff(Surv(time, event) ~ genotype, mid_caf_df)
mid_genotype_pairwise$p.value <- round(mid_genotype_pairwise$p.value, digits = 3)

# ----Plotting KM curves for 7.5 Caf Data----
p <- ggsurvplot(
  fit = mid_caf_mod,
  legend.title = "Genotype",
  legend.labs = c(
    "pdf/Df-1",
    "pdf/Df-2",
    "pdf/Df-3",
    "pdf/Df-4",
    "pdf01-1",
    "pdf01-2",
    "pdf01-3",
    "pdf01-4",
    "w1118-1",
    "w1118-2",
    "w1118-3"
  ),
  xlab = "Time (hours)",
  conf.int = FALSE,
  conf.int.alpha = 0.1,
  font.main = c(18, "bold"),
  font.x = c(16, "bold"),
  font.y = c(16, "bold"),
  font.tickslab = c(14, "plain"),
  legend = "bottom",
  pval = TRUE,
  pval.method = TRUE,
  ggtheme = theme(
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    legend.key = element_blank()
  )
)

p <- p + ggtitle("7.5 mM Caffeine")
p


# ----Saving 7.5 Caf Plot----
ggsave(
  plot = p$plot, filename = "mid_caf_km.png",
  device = "png", path = "Outputs/Plots/",
  width = 8, height = 8, units = "in", dpi = 600
)

ggsave(
  plot = p$plot, filename = "mid_caf_km.svg",
  device = "svg", path = "Outputs/Plots/",
  width = 8, height = 8, units = "in", dpi = 600
)

my_colors <- colorRampPalette(colors = c("red", "yellow", "orange"))(100)

# ----Doing pairwise heat map of Mid Caf results and saving them to file----
png(filename = "Outputs/Plots/mid_caf_heatmap_p_value.png", width = 8, height = 8, units = "in", res = 1080)
mid_heatmap <- heatmap.2(mid_genotype_pairwise$p.value,
  main = "7.5 mM\nPairwise Log-rank",
  Rowv = NA,
  Colv = NA, trace = "none",
  dendrogram = "none", cellnote = mid_genotype_pairwise$p.value,
  notecol = "black", key = TRUE, col = my_colors
)

dev.off()


# ----Loading processed High Caf data----
high_caf_df <- read.csv("Data/high_caff_df.csv", sep = ",")
high_caf_df$genotype <- gsub("\\s+", "", high_caf_df$genotype)

# ----Survival analysis for High Caf Data----
high_caf_mod <- survfit(Surv(time = time, event = event) ~ genotype, data = high_caf_df)

# ----High Caf Genotype Log-rank test. P-value: 7e-05 ----
high_genotype_lrt <- survdiff(Surv(time, event) ~ genotype, high_caf_df)
high_genotype_pairwise <- pairwise_survdiff(Surv(time, event) ~ genotype, high_caf_df)
high_genotype_pairwise$p.value <- round(high_genotype_pairwise$p.value, digits = 3)

# ----Plotting KM curves for 15 Caf Data----
p <- ggsurvplot(
  fit = high_caf_mod,
  legend.title = "Genotype",
  legend.labs = c(
    "pdf/Df-1",
    "pdf/Df-2",
    "pdf/Df-3",
    "pdf/Df-4",
    "pdf01-1",
    "pdf01-2",
    "pdf01-3",
    "pdf01-4",
    "w1118-1",
    "w1118-2"
  ),
  xlab = "Time (hours)",
  conf.int = FALSE,
  conf.int.alpha = 0.1,
  font.main = c(18, "bold"),
  font.x = c(16, "bold"),
  font.y = c(16, "bold"),
  font.tickslab = c(14, "plain"),
  legend = "bottom",
  pval = TRUE,
  pval.method = TRUE,
  ggtheme = theme(
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    legend.key = element_blank()
  )
)

p <- p + ggtitle("15 mM Caffeine")
p


# ----Saving 15 Caf KM Plot----
ggsave(
  plot = p$plot, filename = "high_caf_km.png",
  device = "png", path = "Outputs/Plots/",
  width = 8, height = 8, units = "in", dpi = 600
)

ggsave(
  plot = p$plot, filename = "high_caf_km.svg",
  device = "svg", path = "Outputs/Plots/",
  width = 8, height = 8, units = "in", dpi = 600
)

my_colors <- colorRampPalette(colors = c("red", "yellow", "orange"))(100)

# ----Doing pairwise heat map of High Caf results and saving them to file----
png(filename = "Outputs/Plots/high_caf_heatmap_p_value.png", width = 8, height = 8, units = "in", res = 1080)
high_heatmap <- heatmap.2(high_genotype_pairwise$p.value,
  main = "15 mM\nPairwise Log-rank",
  Rowv = NA,
  Colv = NA, trace = "none",
  dendrogram = "none", cellnote = high_genotype_pairwise$p.value,
  notecol = "black", key = TRUE, col = my_colors
)

dev.off()



# ----Now binding all caf conditions together into a single dataframe to do survival analysis across conditions----
all_data <- bind_rows(no_caf_df, mid_caf_df, high_caf_df)
all_data$fly_id <- 1:nrow(all_data)


# ----Survival analysis for Caf levels----
caf_mod <- survfit(Surv(time = time, event = event) ~ caff_amount, data = all_data)

# ----All Caf Log-rank test. P-value: <2e-16 ----
caf_lrt <- survdiff(Surv(time, event) ~ caff_amount, all_data)
caf_pairwise <- pairwise_survdiff(Surv(time, event) ~ caff_amount, all_data)

# ----Plotting KM curves for Caf Data----
p <- ggsurvplot(
  fit = caf_mod,
  legend.title = "Caffeine",
  legend.labs = c(0, 7.5, 15),
  xlab = "Time (hours)",
  conf.int = FALSE,
  conf.int.alpha = 0.1,
  font.main = c(18, "bold"),
  font.x = c(16, "bold"),
  font.y = c(16, "bold"),
  font.tickslab = c(14, "plain"),
  legend = "bottom",
  pval = TRUE,
  pval.method = TRUE,
  ggtheme = theme(
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    legend.key = element_blank()
  )
)

p <- p + ggtitle("All Caffeine")
p


# ----Saving All Caf Plot----
ggsave(
  plot = p$plot, filename = "all_caf_km.png",
  device = "png", path = "Outputs/Plots/",
  width = 8, height = 8, units = "in", dpi = 600
)

ggsave(
  plot = p$plot, filename = "all_caf_km.svg",
  device = "svg", path = "Outputs/Plots/",
  width = 8, height = 8, units = "in", dpi = 600
)

my_colors <- colorRampPalette(colors = c("red", "yellow", "orange"))(100)

# ----Doing pairwise heat map of all Caf results and saving them to file----
png(filename = "Outputs/Plots/all_caf_heatmap_p_value.png", width = 8, height = 8, units = "in", res = 1080)
caff_heatmap <- heatmap.2(caf_pairwise$p.value,
  main = "All Caffeine\nPairwise Log-rank",
  Rowv = NA,
  Colv = NA, trace = "none",
  dendrogram = "none", cellnote = caf_pairwise$p.value,
  notecol = "black", key = TRUE, col = my_colors
)

dev.off()



# ----Now pairing down to each genotype for comparisons across Caf levels----
perform_genotype_analysis <- function(genotype_condition, data, round_p_value = FALSE, plot_heatmaps = FALSE, display_ci = TRUE, display_pvalue = FALSE, display_pvalue_method = FALSE) {
  return_list <- list()

  # Filter the dataframe based on the genotype condition
  df_gene <- filter(data, genotype == genotype_condition)

  # Check to make sure there are at least 2 Caf levels
  if (length(unique(df_gene$caff_amount)) > 1) {
    # Perform survival analysis
    genotype_mod <- survfit(Surv(time, event) ~ caff_amount, data = df_gene)
    genotype_lrt <- survdiff(Surv(time, event) ~ caff_amount, df_gene)
    genotype_pairwise <- pairwise_survdiff(Surv(time, event) ~ caff_amount, df_gene)
    if (round_p_value) {
      genotype_pairwise$p.value <- round(genotype_pairwise$p.value, digits = 3)
    }

    # Plot KM curves
    p <- ggsurvplot(
      fit = genotype_mod,
      data = df_gene,
      legend.title = "Caffeine (mM)",
      legend.labs = c(0, 7.5, 15),
      xlab = "Time (hours)",
      conf.int = display_ci,
      conf.int.alpha = 0.2,
      font.main = c(18, "bold"),
      font.x = c(16, "bold"),
      font.y = c(16, "bold"),
      font.tickslab = c(14, "plain"),
      legend = "bottom",
      pval = display_pvalue,
      pval.method = display_pvalue_method,
      ggtheme = theme(
        plot.title = element_text(hjust = 0.5),
        panel.background = element_blank(),
        legend.key = element_blank()
      )
    )

    # Set title
    p <- p + ggtitle(paste(genotype_condition))

    return_list[["KM"]] <- p


    # Save plots
    filenames <- c("png", "svg")
    for (filename in filenames) {
      if (genotype_condition == "pdf/Df-1" | genotype_condition == "pdf/Df-2" | genotype_condition == "pdf/Df-3" | genotype_condition == "pdf/Df-4" | genotype_condition == "pdf/Df") {
        km_filename <- gsub("pdf/", "", genotype_condition)
      } else {
        km_filename <- genotype_condition
      }
      ggsave(
        plot = p$plot, filename = paste0(km_filename, "_caf_km.", filename),
        device = filename, path = "Outputs/Plots",
        width = 8, height = 8, units = "in", dpi = 600
      )
    }

    if (plot_heatmaps) {
      # Create heatmap
      my_colors <- colorRampPalette(colors = c("red", "yellow", "orange"))(100)
      heatmap_filename <- file.path("Outputs/Plots", paste0(genotype_condition, "_heatmap_p_value.png"))
      heatmap_filename_svg <- file.path("Outputs/Plots", paste0(genotype_condition, "_heatmap_p_value.svg"))
      
      if (genotype_condition == "pdf/Df-1" | genotype_condition == "pdf/Df-2" | genotype_condition == "pdf/Df-3" | genotype_condition == "pdf/Df-4" | genotype_condition == "pdf/Df") {
        heatmap_filename <- gsub("pdf/", "", heatmap_filename)
      } else {
        heatmap_filename <- heatmap_filename
      }

      # Saving heatmap of pairwise p-values in PNG format
      png(filename = heatmap_filename, width = 8, height = 8, units = "in", res = 1080)
      finished_heatmap <- heatmap.2(genotype_pairwise$p.value,
        main = paste(genotype_condition, "\nPairwise Log-rank"),
        Rowv = NA, Colv = NA, trace = "none", dendrogram = "none",
        cellnote = genotype_pairwise$p.value, notecol = "black",
        key = TRUE, col = my_colors
      )
      dev.off()

      # Saving heatmap of pairwise p-values in SVG format
      svg(filename = heatmap_filename_svg, width = 8, height = 8)
      finished_heatmap <- heatmap.2(genotype_pairwise$p.value,
        main = paste(genotype_condition, "\nPairwise Log-rank"),
        Rowv = NA, Colv = NA, trace = "none", dendrogram = "none",
        cellnote = genotype_pairwise$p.value, notecol = "black",
        key = TRUE, col = my_colors
      )
      dev.off()

      return_list[["Heatmap"]] <- finished_heatmap
    }

    print(paste0(genotype_condition, " analysis completed sucessfully"))

    return(p)
  } else {
    print(paste0(genotype_condition, " genotype is being skipped because it has only 1 caffeine level"))
  }
}


# Looping through all conditions
conds <- unique(all_data$genotype)
all_kms <- list()
for (cond in conds) {
  current_km <- perform_genotype_analysis(genotype_condition = cond, all_data)
  all_kms[[paste0(cond)]] <- current_km
}


# ----Now making combo plot with all genotypes across caf levels----
combo_plot <- arrange_ggsurvplots(all_kms, ncol = 5, nrow = 2, print = FALSE, )
ggsave("all_genotypes_across_caf.png", combo_plot, device = "png", path = "Outputs/Plots/", width = 32, height = 18, units = "in", dpi = 600)


# ----Now treating the 1,2,3,4 etc as replicates and consolidating the data ----
all_data <- all_data %>%
  mutate(genotype = str_remove_all(genotype, "-[1-4]"))


# ----KM plot with all genotypes across specific caffeine conditions with the genotypes consolidated ----
conds <- unique(all_data$genotype)
all_kms <- list()
for (cond in conds) {
  current_km <- perform_genotype_analysis(genotype_condition = cond, all_data, display_ci = TRUE, display_pvalue = FALSE, display_pvalue_method = FALSE)
  all_kms[[paste0(cond)]] <- current_km
}


# ----Now making combo plot with all consolidated genotypes across caf levels----
((all_kms$w1118$plot) / all_kms$pdf01$plot | (all_kms$`pdf/Df`$plot) / plot_spacer()) + plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 18, face = "bold"))

# Saving as SVG
ggsave(filename = "all_consolidated_genotypes_across_caf_global_p_value.svg", 
       plot = last_plot(), 
       device = "svg", 
       path = "Outputs/Plots/", 
       width = 8,
       height = 8, 
       units = "in",
       dpi = 600)



# Saving as PNG
ggsave(filename = "all_consolidated_genotypes_across_caf_global_p_value.png", 
       plot = last_plot(), 
       device = "png", 
       path = "Outputs/Plots/", 
       width = 8,
       height = 8, 
       units = "in",
       dpi = 600)



# ----Now trying to add pairwise log-rank test results to plot----
all_p_values <- list()
all_methods <- list()
conds <- unique(all_data$genotype)

for (cond in conds) {
  current_data <- filter(all_data, genotype == cond)
  res <- pairwise_survdiff(Surv(time, event) ~ caff_amount,
                           data = current_data)
  
  all_p_values[[cond]] <-  signif(res$p.value, digits = 3)
  all_methods[[cond]] <- res$p.adjust.method
}


# ----Making the heatmaps and saving them as SVGs----
svg(filename = "Outputs/Plots/w1118_heatmap.svg", width = 8, height = 8)
par(cex.main=4, cex.lab = 3.5)
heatmap.2(all_p_values$w1118, dendrogram = "none", Rowv = FALSE, trace = "none", main = "w1118", 
          cellnote = all_p_values$w1118, notecol = "black", notecex = 4, cexCol = 4,
          cexRow = 4, key = FALSE, xlab = "Caffeine (mM)", ylab = "Caffeine (mM)")
dev.off()

svg(filename = "Outputs/Plots/pdf01_heatmap.svg", width = 8, height = 8)
par(cex.main=4)
heatmap.2(all_p_values$pdf01, dendrogram = "none", Rowv = FALSE, trace = "none", main = "pdf01", 
          cellnote = all_p_values$pdf01, notecol = "black", notecex = 4, cexCol = 4,
          cexRow = 4, key = FALSE, xlab = "Caffeine (mM)", ylab = "Caffeine (mM)")
dev.off()

svg(filename = "Outputs/Plots/pdf_df_heatmap.svg", width = 8, height = 8)
par(cex.main=4)
heatmap.2(all_p_values$`pdf/Df`, dendrogram = "none", Rowv = FALSE, trace = "none", main = "pdf/Df", 
          cellnote = all_p_values$`pdf/Df`, notecol = "black", notecex = 4, cexCol = 4,
          cexRow = 4, key = FALSE, xlab = "Caffeine (mM)", ylab = "Caffeine (mM)")
dev.off()



# ----Making the heatmaps and saving them as PNGs----
png(filename = "Outputs/Plots/w1118_heatmap.png", width = 8, height = 8, units = "in", res = 1080)
par(cex.main=4)
heatmap.2(all_p_values$w1118, dendrogram = "none", Rowv = FALSE, trace = "none", main = "w1118", 
          cellnote = all_p_values$w1118, notecol = "black", notecex = 4, cexCol = 4,
          cexRow = 4, key = FALSE, xlab = "Caffeine (mM)", ylab = "Caffeine (mM)")
dev.off()

png(filename = "Outputs/Plots/pdf01_heatmap.png", width = 8, height = 8, units = "in", res = 1080)
par(cex.main=4)
heatmap.2(all_p_values$pdf01, dendrogram = "none", Rowv = FALSE, trace = "none", main = "pdf01", 
          cellnote = all_p_values$pdf01, notecol = "black", notecex = 4, cexCol = 4,
          cexRow = 4, key = FALSE, xlab = "Caffeine (mM)", ylab = "Caffeine (mM)")
dev.off()

png(filename = "Outputs/Plots/pdf_df_heatmap.png", width = 8, height = 8, units = "in", res = 1080)
par(cex.main=4)
heatmap.2(all_p_values$`pdf/Df`, dendrogram = "none", Rowv = FALSE, trace = "none", main = "pdf/Df", 
          cellnote = all_p_values$`pdf/Df`, notecol = "black", notecex = 4, cexCol = 4,
          cexRow = 4, key = FALSE, xlab = "Caffeine (mM)", ylab = "Caffeine (mM)")
dev.off()


# ----Making additional pairwise KM plots under 7.5 and 15 mM Caf levels----
all_conds_75 <- all_data %>%
  filter(caff_amount == 7.5)

all_conds_15 <- all_data %>%
  filter(caff_amount == 15)

# ----7.5 mM Caf all genotypes Log-rank test. P-value: 1e-08----
caf_75_genotype_mod <- survfit(Surv(time, event) ~ genotype, data = all_conds_75)
caf_75_genotype_lrt <- survdiff(Surv(time, event) ~ genotype, all_conds_75)
caf_75_genotype_pairwise <- pairwise_survdiff(Surv(time, event) ~ genotype, all_conds_75)

# ----Plotting KM curves for 7.5 Caf all genotype Data----
p75 <- ggsurvplot(
  fit = caf_75_genotype_mod,
  legend.title = "Genotype",
  legend.labs = c("pdf/Df", "pdf01", "w1118"),
  xlab = "Time (hours)",
  conf.int = TRUE,
  conf.int.alpha = 0.2,
  font.main = c(18, "bold"),
  font.x = c(16, "bold"),
  font.y = c(16, "bold"),
  font.tickslab = c(14, "plain"),
  legend = "bottom",
  pval = FALSE,
  pval.method = FALSE,
  ggtheme = theme(
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    legend.key = element_blank()
  )
)

p75 <- p75 + ggtitle("7.5 mM")
p75

# ----Saving 7.5 mM Caf All Genotype Plot----
ggsave(
  plot = p75$plot, filename = "7.5_caf_all_genotype_km.png",
  device = "png", path = "Outputs/Plots/",
  width = 8, height = 8, units = "in", dpi = 600
)

ggsave(
  plot = p75$plot, filename = "7.5_caf_all_genotype_km.svg",
  device = "svg", path = "Outputs/Plots/",
  width = 8, height = 8, units = "in", dpi = 600
)



# ----15 mM Caf all genotypes Log-rank test. P-value: 2e-07----
caf_15_genotype_mod <- survfit(Surv(time, event) ~ genotype, data = all_conds_15)
caf_15_genotype_lrt <- survdiff(Surv(time, event) ~ genotype, all_conds_15)
caf_15_genotype_pairwise <- pairwise_survdiff(Surv(time, event) ~ genotype, all_conds_15)

# ----Plotting KM curves for 15 Caf all genotype Data----
p15 <- ggsurvplot(
  fit = caf_15_genotype_mod,
  legend.title = "Genotype",
  legend.labs = c("pdf/Df", "pdf01", "w1118"),
  xlab = "Time (hours)",
  conf.int = TRUE,
  conf.int.alpha = 0.2,
  font.main = c(18, "bold"),
  font.x = c(16, "bold"),
  font.y = c(16, "bold"),
  font.tickslab = c(14, "plain"),
  legend = "bottom",
  pval = FALSE,
  pval.method = FALSE,
  ggtheme = theme(
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    legend.key = element_blank()
  )
)

p15 <- p15 + ggtitle("15 mM")
p15

# ----Saving 15 mM Caf All Genotype Plot----
ggsave(
  plot = p15$plot, filename = "15_caf_all_genotype_km.png",
  device = "png", path = "Outputs/Plots/",
  width = 8, height = 8, units = "in", dpi = 600
)

ggsave(
  plot = p15$plot, filename = "15_caf_all_genotype_km.svg",
  device = "svg", path = "Outputs/Plots/",
  width = 8, height = 8, units = "in", dpi = 600
)

(p75$plot + p15$plot)/(guide_area()) + plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 18, face = "bold")) 


(p75$plot + p15$plot)/(guide_area()) + plot_layout(ncol = 1, nrow = 2, guides = "collect")+
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 18, face = "bold"))

# Saving as SVG
ggsave(filename = "all_consolidated_genotypes_just_7.5_and_15_caf_separate.svg", 
       plot = last_plot(), 
       device = "svg", 
       path = "Outputs/Plots/", 
       width = 8,
       height = 8, 
       units = "in",
       dpi = 600)


# ----Making and saving pairwise heatmap for 7.5mM Caf level----
all_p_values <- list()
all_methods <- list()
cafs <- unique(all_data$caff_amount)
cafs <- cafs[cafs != 0.0]

for (caf in cafs) {
  current_data <- filter(all_data, caff_amount == caf)
  res <- pairwise_survdiff(Surv(time, event) ~ genotype,
                           data = current_data)
  
  all_p_values[[paste0(as.character(caf))]] <-  signif(res$p.value, digits = 3)
  all_methods[[paste0(as.character(caf))]] <- res$p.adjust.method
}

my_colors <- colorRampPalette(colors = c("red", "yellow", "orange"))(100)
heatmap_filename <- file.path("Outputs/Plots", paste0("7.5_caf_level_all_genotypes_pairwise_heatmap_p_value.svg"))


png(filename = heatmap_filename, width = 8, height = 8, units = "in", res = 1080)
par(cex.main=4)
finished_heatmap <- heatmap.2(all_p_values$`7.5`,
                              main = paste("7.5 mM"),
                              Rowv = NA, Colv = NA, trace = "none", dendrogram = "none",
                              cellnote = all_p_values$`7.5`, notecol = "black", notecex = 4, 
                              cexRow = 4, cexCol = 4,
                              key = FALSE, col = my_colors
)
dev.off()


svg(filename = heatmap_filename, width = 8, height = 8)
par(cex.main=4)
finished_heatmap <- heatmap.2(all_p_values$`7.5`,
                              main = paste("7.5 mM"),
                              Rowv = NA, Colv = NA, trace = "none", dendrogram = "none",
                              cellnote = all_p_values$`7.5`, notecol = "black", notecex = 4, 
                              cexRow = 4, cexCol = 4,
                              key = FALSE, col = my_colors
)
dev.off()

# ----Making and saving pairwise heatmap for 15mM Caf level----
my_colors <- colorRampPalette(colors = c("red", "yellow", "orange"))(100)
heatmap_filename <- file.path("Outputs/Plots", paste0("15_caf_level_all_genotypes_pairwise_heatmap_p_value.png"))


png(filename = heatmap_filename, width = 8, height = 8, units = "in", res = 1080)
par(cex.main=4)
finished_heatmap <- heatmap.2(all_p_values$`15`,
                              main = paste("15 mM"),
                              Rowv = NA, Colv = NA, trace = "none", dendrogram = "none",
                              cellnote = all_p_values$`15`, notecol = "black", notecex = 4, 
                              cexRow = 4, cexCol = 4,
                              key = FALSE, col = my_colors
)
dev.off()


svg(filename = heatmap_filename, width = 8, height = 8)
par(cex.main=4)
finished_heatmap <- heatmap.2(all_p_values$`15`,
                              main = paste("15 mM"),
                              Rowv = NA, Colv = NA, trace = "none", dendrogram = "none",
                              cellnote = all_p_values$`15`, notecol = "black", notecex = 4, 
                              cexRow = 4, cexCol = 4,
                              key = FALSE, col = my_colors
)
dev.off()



# ---- Now separating each genotype within 7.5 mM level to get pairwise KM plots----
# Get unique genotypes
all_pairwise_km_plots <- list()
counter <- 1
plot_tags <- c("E", "F", "G", "H", "I", "J")
genotypes <- unique(all_conds_75$genotype)

# Create all possible combinations of 2 genotypes
genotype_combinations <- combn(genotypes, 2, simplify = FALSE)

# Loop through each genotype combination
for (genotype_pair in genotype_combinations) {
  genotype_1 <- genotype_pair[1]
  genotype_2 <- genotype_pair[2]
  
  # Filter data for each genotype
  data_genotype_1 <- filter(all_conds_75, genotype == genotype_1)
  data_genotype_2 <- filter(all_conds_75, genotype == genotype_2)
  
  # Combine just the two conditions of interest together
  all_genotype_data <- bind_rows(data_genotype_1, data_genotype_2)
  
  # Perform survival analysis
  survfit_1 <- survfit(Surv(time, event) ~ genotype, data = all_genotype_data)
  
  # Create KM plots
  p <- ggsurvplot(
    fit = survfit_1,
    legend.title = "Genotype",
    xlab = "Time (hours)",
    conf.int = TRUE,
    conf.int.alpha = 0.2,
    font.main = c(18, "bold"),
    font.x = c(16, "bold"),
    font.y = c(16, "bold"),
    font.tickslab = c(14, "plain"),
    legend = "bottom",
    pval = TRUE,
    pval.method = TRUE,
    ggtheme = theme(
      plot.title = element_text(hjust = 0.5),
      panel.background = element_blank(),
      legend.key = element_blank()
    )
  )
  
  p <- p + ggtitle(paste("7.5 mM: ", genotype_1, "vs", genotype_2))
  
  # Save plots
  filename_png <- paste0("7.5_caf_genotype_comparison_", genotype_2, "_vs_", genotype_1, ".png")
  filename_svg <- paste0("7.5_caf_genotype_comparison_", genotype_2, "_vs_", genotype_1, ".svg")
  
  ggsave(
    plot = p$plot, filename = filename_png,
    device = "png", path = "Outputs/Plots/",
    width = 8, height = 8, units = "in", dpi = 600
  )
  
  ggsave(
    plot = p$plot, filename = filename_svg,
    device = "svg", path = "Outputs/Plots/",
    width = 8, height = 8, units = "in", dpi = 600
  )
  
  all_pairwise_km_plots[[counter]] <- p + labs(tag = paste0(plot_tags[counter]))
  counter <- counter + 1
}

# ---- Now separating each genotype within 15 mM level to get pairwise KM plots----
# Get unique genotypes
genotypes <- unique(all_conds_15$genotype)

# Create all possible combinations of 2 genotypes
genotype_combinations <- combn(genotypes, 2, simplify = FALSE)

# Loop through each genotype combination
for (genotype_pair in genotype_combinations) {
  genotype_1 <- genotype_pair[1]
  genotype_2 <- genotype_pair[2]
  
  # Filter data for each genotype
  data_genotype_1 <- filter(all_conds_15, genotype == genotype_1)
  data_genotype_2 <- filter(all_conds_15, genotype == genotype_2)
  
  # Combine just the two conditions of interest together
  all_genotype_data <- bind_rows(data_genotype_1, data_genotype_2)
  
  # Perform survival analysis
  survfit_1 <- survfit(Surv(time, event) ~ genotype, data = all_genotype_data)
  
  # Create KM plots
  p <- ggsurvplot(
    fit = survfit_1,
    legend.title = "Genotype",
    xlab = "Time (hours)",
    conf.int = TRUE,
    conf.int.alpha = 0.2,
    font.main = c(18, "bold"),
    font.x = c(16, "bold"),
    font.y = c(16, "bold"),
    font.tickslab = c(14, "plain"),
    legend = "bottom",
    pval = TRUE,
    pval.method = TRUE,
    ggtheme = theme(
      plot.title = element_text(hjust = 0.5),
      panel.background = element_blank(),
      legend.key = element_blank()
    )
  )
  
  p <- p + ggtitle(paste("15 mM: ", genotype_1, "vs", genotype_2))
  
  # Save plots
  filename_png <- paste0("15_caf_genotype_comparison_", genotype_1, "_vs_", genotype_2, ".png")
  filename_svg <- paste0("15_caf_genotype_comparison_", genotype_1, "_vs_", genotype_2, ".svg")
  
  ggsave(
    plot = p$plot, filename = filename_png,
    device = "png", path = "Outputs/Plots/",
    width = 8, height = 8, units = "in", dpi = 600
  )
  
  ggsave(
    plot = p$plot, filename = filename_svg,
    device = "svg", path = "Outputs/Plots/",
    width = 8, height = 8, units = "in", dpi = 600
  )
  
  all_pairwise_km_plots[[counter]] <-  p + 
    labs(tag = paste0(plot_tags[counter])) 
  
  counter <- counter + 1
}



# Making combo plot
combo_plot <- arrange_ggsurvplots(all_pairwise_km_plots, print = FALSE, ncol = 2, nrow = 3)
ggsave("pairwise_genotype_comps_by_caf_level_corrected.png", path = "Outputs/Plots/", 
       device = "png", width = 16, height = 16, units = "in", dpi = 600, plot = combo_plot)

ggsave("pairwise_genotype_comps_by_caf_level_corrected.svg", path = "Outputs/Plots/", 
       device = "svg", width = 16, height = 16, units = "in", dpi = 600, plot = combo_plot)

# ---- Starvation Analysis ----
# ----Loading processed starvation data----
starve_df_0 <- read.csv("Data/starvation_data_0_mM.csv")
starve_df_75 <- read.csv("Data/starvation_data_75_mM.csv")

# ----7.5 mM Genotype within same condition----
starved_75_df <- starve_df_75 %>%
  filter(condition == "starved")

colnames(starved_75_df)[4] <- "genotype"

caf_75_genotype_mod <- survfit(Surv(time, event) ~ genotype, data = starved_75_df)
caf_75_genotype_lrt <- survdiff(Surv(time, event) ~ genotype, data = starved_75_df)
caf_75_genotype_pairwise <- pairwise_survdiff(Surv(time, event) ~ genotype, starved_75_df)

# ----Plotting KM curves for 7.5 Caf all genotype Data----
p75 <- ggsurvplot(
  fit = caf_75_genotype_mod,
  legend.title = "Genotype",
  xlab = "Time (hours)",
  conf.int = TRUE,
  conf.int.alpha = 0.2,
  font.main = c(18, "bold"),
  font.x = c(16, "bold"),
  font.y = c(16, "bold"),
  font.tickslab = c(14, "plain"),
  legend = "bottom",
  pval = FALSE,
  pval.method = FALSE,
  ggtheme = theme(
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    legend.key = element_blank()
  )
)

p75 <- p75 + ggtitle("7.5 mM\nStarved")
p75

# ----Saving 7.5 mM Caf Starved All Genotype Plot----
ggsave(
  plot = p75$plot, filename = "7.5_caf_all_genotype_starved_km.png",
  device = "png", path = "Outputs/Plots/",
  width = 8, height = 8, units = "in", dpi = 600
)

ggsave(
  plot = p75$plot, filename = "7.5_caf_all_genotype_starved_km.svg",
  device = "svg", path = "Outputs/Plots/",
  width = 8, height = 8, units = "in", dpi = 600
)



# ----7.5 mM w1118 conditions compared----
w1118_conds <- 

colnames(starved_75_df)[4] <- "genotype"

caf_75_genotype_mod <- survfit(Surv(time, event) ~ genotype, data = starved_75_df)
caf_75_genotype_lrt <- survdiff(Surv(time, event) ~ genotype, data = starved_75_df)
caf_75_genotype_pairwise <- pairwise_survdiff(Surv(time, event) ~ genotype, starved_75_df)

# ----Plotting KM curves for 7.5 Caf all genotype Data----
p75 <- ggsurvplot(
  fit = caf_75_genotype_mod,
  legend.title = "Genotype",
  xlab = "Time (hours)",
  conf.int = TRUE,
  conf.int.alpha = 0.2,
  font.main = c(18, "bold"),
  font.x = c(16, "bold"),
  font.y = c(16, "bold"),
  font.tickslab = c(14, "plain"),
  legend = "bottom",
  pval = FALSE,
  pval.method = FALSE,
  ggtheme = theme(
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    legend.key = element_blank()
  )
)

p75 <- p75 + ggtitle("7.5 mM\nStarved")
p75

# ----Saving 7.5 mM Caf Starved All Genotype Plot----
ggsave(
  plot = p75$plot, filename = "7.5_caf_all_genotype_starved_km.png",
  device = "png", path = "Outputs/Plots/",
  width = 8, height = 8, units = "in", dpi = 600
)

ggsave(
  plot = p75$plot, filename = "7.5_caf_all_genotype_starved_km.svg",
  device = "svg", path = "Outputs/Plots/",
  width = 8, height = 8, units = "in", dpi = 600
)


