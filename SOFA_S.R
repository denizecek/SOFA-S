# SOFA-S Run Script 12.03.25


seed = 123
scale = TRUE
normalise = TRUE
bin_size = 150

run_sofa <- function(coord_table_path,
                     omic_data_path,
                     metadata_path,
                     survival_table_path,
                     hdf5_output = "sofa_trained_model.hdf5",
                     normalise = TRUE,
                     scale = TRUE,
                     seed = 123,
                     bin_size = 150) {
  
  # Load required packages
  library(MOFA2)
  library(tidyverse)
  library(patchwork)
  library(survival)
  library(survminer)
  library(broom)
  library(ComplexHeatmap)
  library(RColorBrewer)
  
  # Set seed
  set.seed(seed)
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Update full path for model
  hdf5_output_path <- file.path(output_dir, hdf5_output)
  
  # Read input data
  coords <- read.csv(coord_table_path)
  omic_data <- read.csv(omic_data_path)
  metadata <- read.csv(metadata_path)
  survival_metadata <- read.csv(survival_table_path)
  
  
  # Define marker columns as all columns in omics data
  marker_cols <- colnames(omic_data)
  metadata$Xnew <- ceiling(coords$X / bin_size)
  metadata$Ynew <- ceiling(coords$Y / bin_size)
  full_data <- cbind(metadata, omic_data)
  
  # Columns to keep from metadata (exclude X, Y, Xnew, Ynew, omics)
  metadata_cols <- setdiff(colnames(metadata), c("Xnew", "Ynew","ROI"))
  
  # Summarise omics
  gbm_data <- full_data %>%
    group_by(Xnew, Ynew, ROI) %>%
    summarise(
      across(all_of(metadata_cols), ~ first(.)),
      across(all_of(marker_cols), sum, na.rm = TRUE), .groups = "drop")
  
  # Plot and save scatter
  coord_plot <- ggplot(gbm_data, aes(x = Xnew, y = Ynew)) +
    geom_point(color = "blue", alpha = 0.6) +
    theme_minimal() +
    labs(
      title = "Scatter Plot of Binned X and Y Coordinates",
      x = "X Coordinate (Binned)",
      y = "Y Coordinate (Binned)"
    )
  
  ggsave(filename = file.path(output_dir, "scatter_plot.png"), plot = coord_plot, width = 6, height = 5)
  
  # Omics preprocessing
  omic_data_selected <- gbm_data[, marker_cols]
  if (normalise && scale) {
    omic_data_logit_apply <- log(omic_data_selected + 1)
    omic_data_norm_scaled <- scale(omic_data_logit_apply, center = TRUE, scale = TRUE)
    normalised_omic_data <- t(omic_data_norm_scaled)
  } else {
    normalised_omic_data <- t(omic_data_selected)
  }
  matrices <- list(omic = as.matrix(normalised_omic_data))
  
  coords_mat <- as.matrix(gbm_data[, c("Xnew", "Ynew")])
  coords_mat <- t(coords_mat)
  
  # MOFA model
  #MOFA_clin <- create_mofa(matrices,groups = clinical$ROI)
  MOFA_clin <- create_mofa(matrices)
  colnames(coords_mat) <- MOFA_clin@samples_metadata$sample
  MOFA_clin <- set_covariates(MOFA_clin, covariates = coords_mat)
  
  data_opts <- get_default_data_options(MOFA_clin)
  model_opts <- get_default_model_options(MOFA_clin)
  train_opts <- get_default_training_options(MOFA_clin)
  mefisto_opts <- get_default_mefisto_options(MOFA_clin)
  
  MOFA_clin <- prepare_mofa(
    object = MOFA_clin,
    data_options = data_opts,
    model_options = model_opts,
    mefisto_options = mefisto_opts,
    training_options = train_opts
  ) 
  
  gbm.trained <- run_mofa(MOFA_clin, hdf5_output_path, use_basilisk = TRUE)
  
  # Metadata
  
  gbm.trained@samples_metadata <- gbm_data[, c("sample", "Xnew", "Ynew", metadata_cols)]
  
  
  cov= as.data.frame(gbm.trained@covariates$group1)
  colnames(cov)
  gbm.trained@samples_metadata$sample = colnames(cov)
  
  # Variance explained plot
  var_exp_plot <- plot_variance_explained(gbm.trained, max_r2 = 15, y = "view", x = "factor")
  
  # Save the plot
  ggsave(
    filename = file.path(output_dir, "variance_explained.png"),
    plot = var_exp_plot,
    width = 7, height = 5
  )
  
  # Get number of factors
  n_factors <- gbm.trained@dimensions$K
  
  # Create subfolders for weights
  weights_dir <- file.path(output_dir, "weights")
  dir.create(weights_dir, showWarnings = FALSE)
  
  # Initialize lists to store ggplot objects
  weights_plots <- list()
  top_weights_plots <- list()
  
  for (i in 1:n_factors) {
    # Plot all weights for factor i
    w_plot <- plot_weights(gbm.trained, view = 1, factors = i)
    weights_plots[[paste0("Factor", i)]] <- w_plot
    
    ggsave(
      filename = file.path(weights_dir, paste0("weights_factor_", i, ".png")),
      plot = w_plot,
      width = 6, height = 4
    )
    
    # Plot top 50 weights
    top_plot <- plot_top_weights(gbm.trained, view = 1, factors = i, nfeatures = 50)
    top_weights_plots[[paste0("Factor", i)]] <- top_plot
    
    ggsave(
      filename = file.path(weights_dir, paste0("top_weights_factor_", i, ".png")),
      plot = top_plot,
      width = 6, height = 4
    )
  }
  
  weights_longdf <- get_weights(gbm.trained, as.data.frame = TRUE, scale = TRUE)
  weights_longdf <- as.data.frame(weights_longdf)
  
  significants <- as.character(weights_longdf$feature)
  weights_omic <- gbm.trained@expectations$W$omic
  weights <- rbind(weights_omic)
  
  fac <- weights[rownames(weights) %in% unique(significants), ]
  
  # Prepare heatmap PNG output path
  heatmap_path <- file.path(output_dir, "weights_heatmap.png")
  
  # Open PNG device
  png(filename = heatmap_path, width = 1200, height = 1200, res = 150)
  
  # Draw heatmap
  ht <- ComplexHeatmap::Heatmap(
    fac,
    name = "weights",
    gap = grid::unit(2, "mm"),
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    show_row_dend = TRUE,
    row_names_gp = grid::gpar(fontsize = 8)
  )
  
  ComplexHeatmap::draw(ht)
  dev.off()
  
  
  # Generate smoothness plot
  smoothness_plot <- plot_smoothness(gbm.trained)
  
  # Define path
  smoothness_path <- file.path(output_dir, "factor_smoothness.png")
  
  # Save it as PNG
  ggsave(
    filename = smoothness_path,
    plot = smoothness_plot,
    width = 7, height = 5, dpi = 300
  )
  
  # Factors vs Covariates plot
  factors_vs_cov_plot <- plot_factors_vs_cov(gbm.trained)
  
  # Save it
  factors_vs_cov_path <- file.path(output_dir, "factors_vs_cov.png")
  ggsave(
    filename = factors_vs_cov_path,
    plot = factors_vs_cov_plot,
    width = 7, height = 5, dpi = 300
  )
 
  # Create subfolder for data_vs_cov plots
  data_vs_cov_dir <- file.path(output_dir, "data_vs_cov")
  dir.create(data_vs_cov_dir, showWarnings = FALSE)
  
  # Number of factors
  n_factors <- gbm.trained@dimensions$K
  
  # Store plots in a list
  data_vs_cov_plots <- list()
  
  for (i in 1:n_factors) {
    plot_i <- plot_data_vs_cov(
      gbm.trained,
      factor = i,
      color_by = "group",
      dot_size = 2
    )
    
    # Save the plot
    plot_path <- file.path(data_vs_cov_dir, paste0("data_vs_cov_factor", i, ".png"))
    ggsave(
      filename = plot_path,
      plot = plot_i,
      width = 7, height = 5, dpi = 300
    )
    
    # Store in return list
    data_vs_cov_plots[[paste0("Factor", i)]] <- plot_i
  }
  
  # Factor scores
  factor_scores <- get_factors(gbm.trained, factors = "all") %>% as.data.frame()
  # factor_scores$PatientID <- 1:nrow(factor_scores)
  # factor_scores$Time <- 250:(249 + nrow(factor_scores))
  # factor_scores$Event <- 1
  # Ensure factor_scores has 'sample' column for joining
  
  write.csv(factor_scores, file.path(output_dir, "factor_scores.csv"), row.names = FALSE)
  factor_scores$sample <- paste0("sample_", 1:nrow(factor_scores))
  # Cox models
  cox_results <- lapply(1:5, function(i) {
    factor_name <- paste0("group1.Factor", i)
    coxph(Surv(Time, Event) ~ factor_scores[[factor_name]], data = factor_scores)
  })
  cox_summary <- broom::tidy(cox_results[[1]], conf.int = TRUE)
  write.csv(cox_summary, file.path(output_dir, "cox_summary.csv"), row.names = FALSE)
  
  km_plots <- list()
  
  for (i in 1:n_factors) {
    factor_name <- paste0("group1.Factor", i)
    group_col <- paste0(factor_name, "_group")
    
    # Binarize the factor scores
    factor_scores[[group_col]] <- ifelse(
      factor_scores[[factor_name]] > median(factor_scores[[factor_name]], na.rm = TRUE),
      "High", "Low"
    )
    
    # Create KM fit
    surv_object <- Surv(factor_scores[["Time"]], factor_scores[["Event"]])
    fit <- survfit(surv_object ~ factor_scores[[group_col]])
    
    # Generate plot
    km_plot <- ggsurvplot(
      fit,
      data = factor_scores,
      pval = TRUE,
      title = paste("Survival Analysis for", factor_name),
      xlab = "Time (days)",
      ylab = "Survival Probability"
    )
    
    # Save to list
    km_plots[[paste0("Factor", i)]] <- km_plot$plot
  }
  for (i in 1:n_factors) {
    ggsave(
      filename = file.path(output_dir, paste0("km_factor_", i, ".png")),
      plot = km_plots[[paste0("Factor", i)]],
      width = 6, height = 5, dpi = 300
    )
  }
  
  # Create combined data: Time, Event, and factor columns
  combined_data <- factor_scores[, c("Time", "Event", paste0("group1.Factor", 1:n_factors))]
  
  # Build Cox formula
  cox_formula <- as.formula(
    paste0("Surv(Time, Event) ~ ", paste0("group1.Factor", 1:n_factors, collapse = " + "))
  )
  
  # Fit Cox model
  cox_model <- coxph(cox_formula, data = combined_data)
  
  # Tidy summary with confidence intervals
  cox_summary <- broom::tidy(cox_model, conf.int = TRUE)
  
  # Save forest plot (grid-based)
  forest_path <- file.path(output_dir, "cox_forest.png")
  png(forest_path, width = 800, height = 600, res = 120)
  ggforest(cox_model, data = combined_data)
  dev.off()
  
  # Create ggplot HR plot
  hr_plot <- ggplot(cox_summary, aes(x = term, y = estimate, ymin = conf.low, ymax = conf.high)) +
    geom_pointrange() +
    coord_flip() +
    labs(
      title = "Hazard Ratios for MOFA Factors",
      x = "Factors",
      y = "Hazard Ratio (95% CI)"
    )
  
  # Save HR plot
  hr_path <- file.path(output_dir, "cox_hr_plot.png")
  ggsave(hr_path, plot = hr_plot, width = 7, height = 5, dpi = 300)
  
  
  # Return output
  list(
    # Core models and data
    mofa_model = gbm.trained,
    factor_scores = factor_scores,
    cox_model = cox_model,
    cox_results = cox_results,
    cox_summary = cox_summary,
    
    # ggplot objects
    scatter_plot = coord_plot,
    variance_explained_plot = var_exp_plot,
    smoothness_plot = smoothness_plot,
    factors_vs_cov_plot = factors_vs_cov_plot,
    hr_plot = hr_plot,
    weights_plots = weights_plots,
    top_weights_plots = top_weights_plots,
    data_vs_cov_plots = data_vs_cov_plots,
    km_plots = km_plots,
    
    # Output files and folders
    output_files = list(
      # Tables
      hdf5 = hdf5_output_path,
      factor_scores_csv = file.path(output_dir, "factor_scores.csv"),
      cox_summary_csv = file.path(output_dir, "cox_summary.csv"),
      
      # Core plots
      scatter = file.path(output_dir, "scatter_plot.png"),
      variance_explained = file.path(output_dir, "variance_explained.png"),
      smoothness = file.path(output_dir, "factor_smoothness.png"),
      factors_vs_cov = file.path(output_dir, "factors_vs_cov.png"),
      hr_plot = hr_path,
      cox_forest = forest_path,
      
      # Folders containing grouped plots
      weights_dir = weights_dir,
      data_vs_cov_dir = data_vs_cov_dir,
      
      # Individual KM plots
      km_plots = list.files(output_dir, pattern = "km_factor_.*\\.png$", full.names = TRUE),
      
      # Heatmap
      heatmap = heatmap_path
    )
  )
}