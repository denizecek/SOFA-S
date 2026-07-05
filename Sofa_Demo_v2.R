# ============================================================
# SOFA-S demo pipeline
# ------------------------------------------------------------
# Purpose:
#   Run a spatial omics factor analysis workflow using MOFA2/MEFISTO,
#   starting from:
#     1. cell/bin coordinates
#     2. marker expression matrix
#     3. ROI-level clinical metadata including survival information
#
# Key assumptions:
#   - coord_table, omic_data and metadata have the same row order.
#   - Each row represents one cell before binning.
#   - metadata contains one row per cell and includes an ROI/core identifier.
#   - Survival variables are named Time and Event, where Event is usually:
#       0 = alive
#       1 = death
#
# ============================================================

run_sofa <- function(coord_table_path,
                     omic_data_path,
                     metadata_path,
                     output_dir,
                     hdf5_output = "sofa_trained_model.hdf5",
                     normalise = TRUE,
                     scale = TRUE,
                     seed = 123,
                     bin_size = 500) {

  # ------------------------------------------------------------
  # 1. Load required packages
  # ------------------------------------------------------------
  suppressPackageStartupMessages({
    library(MOFA2)
    library(tidyverse)
    library(survival)
    library(survminer)
    library(broom)
    library(ComplexHeatmap)
  })

  set.seed(seed)

  # Use cairo graphics device where available.
  # This helps avoid X11-related PNG errors on HPC systems.
  options(bitmapType = "cairo")

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  hdf5_output_path <- file.path(output_dir, hdf5_output)

  # ------------------------------------------------------------
  # 2. Read input data
  # ------------------------------------------------------------
  coords <- read.csv(coord_table_path, check.names = FALSE)
  omic_data <- read.csv(omic_data_path, check.names = FALSE)
  metadata <- read.csv(metadata_path, check.names = FALSE)

  # Basic checks to ensure the three input tables are aligned.
  stopifnot(nrow(coords) == nrow(omic_data))
  stopifnot(nrow(metadata) == nrow(omic_data))

  # Check required coordinate columns.
  if (!all(c("X", "Y") %in% colnames(coords))) {
    stop("coord_table must contain columns named X and Y")
  }

  # Check bin size.
  if (!is.numeric(bin_size) || length(bin_size) != 1 || bin_size <= 0) {
    stop("bin_size must be a single positive numeric value")
  }

  # ------------------------------------------------------------
  # 3. Standardise metadata columns
  # ------------------------------------------------------------
  # The workflow expects an ROI column. If the metadata uses 'core',
  # copy it to ROI. If neither exists, create a single dummy ROI.
  if (!"ROI" %in% colnames(metadata)) {
    if ("core" %in% colnames(metadata)) {
      metadata$ROI <- metadata$core
    } else {
      warning("No ROI/core column found. All rows will be assigned to ROI1.")
      metadata$ROI <- "ROI1"
    }
  }

  # The workflow expects a sample column. Here sample is initially cell-level;
  # after binning we will replace it with a unique bin-level sample ID.
  if (!"sample" %in% colnames(metadata)) {
    metadata$sample <- paste0("sample_", seq_len(nrow(metadata)))
  }

  # ------------------------------------------------------------
  # 4. Spatial binning
  # ------------------------------------------------------------
  # Convert continuous X/Y coordinates into discrete spatial bins.
  # Smaller bin_size gives higher spatial resolution but increases runtime.
  metadata$Xnew <- ceiling(coords$X / bin_size)
  metadata$Ynew <- ceiling(coords$Y / bin_size)

  # Remove accidental overlap between marker columns and metadata columns.
  # This prevents metadata variables from being treated as omics features.
  omic_data <- omic_data[, !colnames(omic_data) %in% colnames(metadata), drop = FALSE]
  marker_cols <- colnames(omic_data)

  if (length(marker_cols) == 0) {
    stop("No marker columns remain after removing metadata overlaps.")
  }

  # Combine metadata and marker expression for bin-level aggregation.
  full_data <- bind_cols(metadata, omic_data)

  # Metadata columns to carry forward after aggregation.
  # Xnew/Ynew/ROI are grouping variables and should not be summarised.
  metadata_cols <- setdiff(colnames(metadata), c("Xnew", "Ynew", "ROI"))

  # Aggregate cell-level marker expression into spatial bins.
  # Currently using sum(), which captures total marker signal per bin.
  sofa_data <- full_data %>%
    group_by(Xnew, Ynew, ROI) %>%
    summarise(
      n_cells = n(),
      across(all_of(metadata_cols), ~ first(.x)),
      across(all_of(marker_cols), ~ sum(.x, na.rm = TRUE)),
      .groups = "drop"
    )

  # Unique sample name for each spatial bin.
  sofa_data$bin_sample <- paste0(
    sofa_data$ROI,
    "_X", sofa_data$Xnew,
    "_Y", sofa_data$Ynew
  )

  # ------------------------------------------------------------
  # 5. Quality-control plots before MOFA
  # ------------------------------------------------------------
  # Plot spatial bins by ROI. Useful for checking bin resolution.
  coord_plot <- ggplot(sofa_data, aes(x = Xnew, y = Ynew)) +
    geom_point(size = 1, alpha = 0.6) +
    facet_wrap(~ ROI) +
    coord_equal() +
    theme_bw() +
    theme(strip.text = element_text(size = 5)) +
    labs(
      title = "Spatial bins by ROI",
      x = "X bin",
      y = "Y bin"
    )

  ggsave(
    file.path(output_dir, "scatter_plot_by_ROI.png"),
    coord_plot,
    width = 14,
    height = 10,
    dpi = 300
  )

  # Plot number of cells per spatial bin.
  cell_count_plot <- ggplot(sofa_data, aes(x = Xnew, y = Ynew, fill = n_cells)) +
    geom_tile(color = "grey80") +
    facet_wrap(~ ROI) +
    coord_equal() +
    theme_minimal() +
    labs(
      title = "Number of cells per spatial bin",
      x = "Binned X coordinate",
      y = "Binned Y coordinate",
      fill = "n cells"
    )

  ggsave(
    file.path(output_dir, "cell_count_by_ROI.png"),
    cell_count_plot,
    width = 14,
    height = 10,
    dpi = 300
  )

  # Plot total number of cells per ROI.
  roi_counts <- full_data %>%
    count(ROI)

  roi_barplot <- ggplot(
    roi_counts,
    aes(x = reorder(ROI, n), y = n)
  ) +
    geom_col() +
    coord_flip() +
    theme_bw() +
    labs(
      title = "Cell count per ROI",
      x = "ROI",
      y = "Number of cells"
    )

  ggsave(
    file.path(output_dir, "cells_per_roi.png"),
    roi_barplot,
    width = 8,
    height = 10,
    dpi = 300
  )

  # ------------------------------------------------------------
  # 6. Prepare omics matrix for MOFA
  # ------------------------------------------------------------
  omic_data_selected <- sofa_data[, marker_cols, drop = FALSE]

  # Remove features with zero variance or no observations.
  # MOFA cannot learn meaningful factors from constant features.
  keep_features <- apply(omic_data_selected, 2, function(x) {
    sd(x, na.rm = TRUE) > 0 && sum(!is.na(x)) > 0
  })

  omic_data_selected <- omic_data_selected[, keep_features, drop = FALSE]
  marker_cols <- colnames(omic_data_selected)

  if (length(marker_cols) == 0) {
    stop("All marker features were removed during zero-variance filtering.")
  }

  # Normalisation:
  #   log1p transformation reduces skew from high marker counts.
  #   scale() centres and scales each marker across spatial bins.
  if (normalise && scale) {
    omic_data_log <- log1p(omic_data_selected)
    omic_data_norm_scaled <- scale(omic_data_log, center = TRUE, scale = TRUE)
    normalised_omic_data <- t(omic_data_norm_scaled)
  } else {
    normalised_omic_data <- t(as.matrix(omic_data_selected))
  }

  # MOFA expects features as rows and samples as columns.
  colnames(normalised_omic_data) <- sofa_data$bin_sample
  matrices <- list(omic = as.matrix(normalised_omic_data))

  # MEFISTO covariates: spatial coordinates for each bin.
  coords_mat <- as.matrix(sofa_data[, c("Xnew", "Ynew")])
  coords_mat <- t(coords_mat)
  colnames(coords_mat) <- sofa_data$bin_sample

  # ------------------------------------------------------------
  # 7. Build and train MOFA/MEFISTO model
  # ------------------------------------------------------------
  MOFA_clin <- create_mofa(matrices)
  MOFA_clin <- set_covariates(MOFA_clin, covariates = coords_mat)

  MOFA_clin <- prepare_mofa(
    object = MOFA_clin,
    data_options = get_default_data_options(MOFA_clin),
    model_options = get_default_model_options(MOFA_clin),
    mefisto_options = get_default_mefisto_options(MOFA_clin),
    training_options = get_default_training_options(MOFA_clin)
  )

  sofa_trained <- run_mofa(
    MOFA_clin,
    outfile = hdf5_output_path,
    use_basilisk = TRUE
  )

  # Attach bin-level metadata to the trained MOFA object.
  sofa_data$sample <- sofa_data$bin_sample
  sofa_trained@samples_metadata <- sofa_data

  # ------------------------------------------------------------
  # 8. MOFA model diagnostics and weights
  # ------------------------------------------------------------
  var_exp_plot <- plot_variance_explained(
    sofa_trained,
    max_r2 = 15,
    y = "view",
    x = "factor"
  )

  ggsave(
    file.path(output_dir, "variance_explained.png"),
    var_exp_plot,
    width = 7,
    height = 5
  )

  n_factors <- sofa_trained@dimensions$K

  weights_dir <- file.path(output_dir, "weights")
  dir.create(weights_dir, showWarnings = FALSE)

  # Save factor loading plots for each factor.
  for (i in seq_len(n_factors)) {

    w_plot <- plot_weights(
      sofa_trained,
      view = 1,
      factors = i
    )

    ggsave(
      file.path(weights_dir, paste0("weights_factor_", i, ".png")),
      w_plot,
      width = 6,
      height = 4
    )

    top_plot <- plot_top_weights(
      sofa_trained,
      view = 1,
      factors = i,
      nfeatures = min(50, length(marker_cols))
    )

    ggsave(
      file.path(weights_dir, paste0("top_weights_factor_", i, ".png")),
      top_plot,
      width = 6,
      height = 4
    )
  }

  # Heatmap of marker weights across all factors.
  # PDF avoids X11 graphics errors on HPC.
  weights_omic <- sofa_trained@expectations$W$omic

  pdf(file.path(output_dir, "weights_heatmap.pdf"), width = 8, height = 8)

  ht <- ComplexHeatmap::Heatmap(
    weights_omic,
    name = "weights",
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    show_row_dend = TRUE,
    row_names_gp = grid::gpar(fontsize = 8)
  )

  ComplexHeatmap::draw(ht)
  dev.off()

  smoothness_plot <- plot_smoothness(sofa_trained)

  ggsave(
    file.path(output_dir, "factor_smoothness.png"),
    smoothness_plot,
    width = 7,
    height = 5,
    dpi = 300
  )

  factors_vs_cov_plot <- plot_factors_vs_cov(sofa_trained)

  ggsave(
    file.path(output_dir, "factors_vs_cov.png"),
    factors_vs_cov_plot,
    width = 7,
    height = 5,
    dpi = 300
  )

  # ------------------------------------------------------------
  # 9. Extract factor scores
  # ------------------------------------------------------------
  factor_scores_long <- get_factors(
    sofa_trained,
    factors = "all",
    as.data.frame = TRUE
  ) %>%
    as.data.frame()

  # Add spatial-bin metadata back to each factor score.
  factor_scores_long <- factor_scores_long %>%
    left_join(
      sofa_trained@samples_metadata,
      by = "sample"
    )

  write.csv(
    factor_scores_long,
    file.path(output_dir, "factor_scores_long.csv"),
    row.names = FALSE
  )

  factor_names <- unique(factor_scores_long$factor)

  # Convert factor scores to wide bin-level table.
  # One row = one spatial bin.
  factor_scores_wide_bin <- factor_scores_long %>%
    select(sample, factor, value) %>%
    pivot_wider(
      names_from = factor,
      values_from = value
    ) %>%
    left_join(
      sofa_trained@samples_metadata,
      by = "sample"
    )

  write.csv(
    factor_scores_wide_bin,
    file.path(output_dir, "factor_scores_wide_bin_level.csv"),
    row.names = FALSE
  )

  # ------------------------------------------------------------
  # 10. Spatial factor maps
  # ------------------------------------------------------------
  spatial_factor_dir <- file.path(output_dir, "spatial_factor_maps_by_ROI")
  dir.create(spatial_factor_dir, showWarnings = FALSE)

  spatial_factor_plots <- list()

  # This combined plot can be very large if many ROIs are included.
  # It is useful for small subsets, but for large cohorts the per-factor
  # or per-ROI plots below are more readable.
  p_all <- ggplot(
    factor_scores_long,
    aes(x = Xnew, y = Ynew, fill = value)
  ) +
    geom_tile(color = "grey80") +
    facet_grid(factor ~ ROI) +
    coord_equal() +
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      na.value = "grey90"
    ) +
    theme_minimal() +
    labs(
      title = "Spatial distribution of SOFA factors by ROI",
      x = "Binned X coordinate",
      y = "Binned Y coordinate",
      fill = "Factor value"
    )

  ggsave(
    filename = file.path(spatial_factor_dir, "spatial_maps_all_factors_by_ROI.png"),
    plot = p_all,
    width = 14,
    height = max(8, length(factor_names) * 2.5),
    dpi = 300
  )

  # Save one faceted spatial map per factor.
  for (f in factor_names) {

    p <- factor_scores_long %>%
      filter(factor == f) %>%
      ggplot(aes(x = Xnew, y = Ynew, fill = value)) +
      geom_tile(color = "grey80") +
      facet_wrap(~ ROI) +
      coord_equal() +
      scale_fill_gradient2(
        low = "blue",
        mid = "white",
        high = "red",
        midpoint = 0,
        na.value = "grey90"
      ) +
      theme_minimal() +
      labs(
        title = paste("Spatial distribution of", f),
        x = "Binned X coordinate",
        y = "Binned Y coordinate",
        fill = "Factor value"
      )

    spatial_factor_plots[[f]] <- p

    safe_f <- gsub("[^A-Za-z0-9_]", "_", f)

    ggsave(
      filename = file.path(spatial_factor_dir, paste0("spatial_map_", safe_f, "_by_ROI.png")),
      plot = p,
      width = 12,
      height = 8,
      dpi = 300
    )

    # Save one separate map per factor and ROI.
    # This is useful for manual inspection but can create many files.
    for (current_roi in unique(factor_scores_long$ROI)) {

      p_roi <- factor_scores_long %>%
        filter(factor == f, ROI == current_roi) %>%
        ggplot(aes(x = Xnew, y = Ynew, fill = value)) +
        geom_tile(color = "grey80") +
        coord_equal() +
        scale_fill_gradient2(
          low = "blue",
          mid = "white",
          high = "red",
          midpoint = 0,
          na.value = "grey90"
        ) +
        theme_bw() +
        labs(
          title = paste(f, "-", current_roi),
          x = "Binned X coordinate",
          y = "Binned Y coordinate",
          fill = "Factor value"
        )

      safe_roi <- gsub("[^A-Za-z0-9_]", "_", current_roi)

      ggsave(
        file.path(
          spatial_factor_dir,
          paste0(safe_f, "_", safe_roi, ".png")
        ),
        p_roi,
        width = 5,
        height = 5,
        dpi = 300
      )
    }
  }

  # ------------------------------------------------------------
  # 11. ROI-level factor summaries for valid survival analysis
  # ------------------------------------------------------------
  #   MOFA was trained on spatial bins, but survival analysis must not treat
  #   bins from the same ROI as independent samples.
  #
  # Therefore, we summarise factor scores per ROI before survival testing.
  # Mean captures overall factor activity.
  # SD captures within-ROI spatial heterogeneity.
  roi_factor_scores_long <- factor_scores_long %>%
    group_by(ROI, Event, Time, factor) %>%
    summarise(
      mean_factor = mean(value, na.rm = TRUE),
      sd_factor = sd(value, na.rm = TRUE),
      max_factor = max(value, na.rm = TRUE),
      min_factor = min(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Event_status = ifelse(Event == 1, "Death", "Alive")
    )

  write.csv(
    roi_factor_scores_long,
    file.path(output_dir, "roi_factor_scores_long.csv"),
    row.names = FALSE
  )

  roi_factor_scores_wide <- roi_factor_scores_long %>%
    select(ROI, Event, Time, factor, mean_factor) %>%
    pivot_wider(
      names_from = factor,
      values_from = mean_factor
    )

  write.csv(
    roi_factor_scores_wide,
    file.path(output_dir, "roi_factor_scores_wide_mean.csv"),
    row.names = FALSE
  )

  # ------------------------------------------------------------
  # 12. Alive/Death comparison using ROI-level factor means
  # ------------------------------------------------------------
  if (all(c("Event") %in% colnames(roi_factor_scores_long))) {

    factor_event_roi_plot <- ggplot(
      roi_factor_scores_long,
      aes(x = Event_status, y = mean_factor)
    ) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
      facet_wrap(~ factor, scales = "free_y") +
      theme_bw() +
      labs(
        title = "Mean SOFA factor scores by survival status",
        x = "Survival status",
        y = "Mean factor score per ROI"
      )

    ggsave(
      file.path(output_dir, "mean_factor_scores_by_event_ROI_level.png"),
      factor_event_roi_plot,
      width = 10,
      height = 6,
      dpi = 300
    )

    # Optional: within-ROI spatial heterogeneity by survival status.
    factor_sd_roi_plot <- ggplot(
      roi_factor_scores_long,
      aes(x = Event_status, y = sd_factor)
    ) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
      facet_wrap(~ factor, scales = "free_y") +
      theme_bw() +
      labs(
        title = "Spatial heterogeneity of SOFA factors by survival status",
        x = "Survival status",
        y = "Within-ROI SD of factor score"
      )

    ggsave(
      file.path(output_dir, "factor_spatial_heterogeneity_by_event_ROI_level.png"),
      factor_sd_roi_plot,
      width = 10,
      height = 6,
      dpi = 300
    )
  }

  # ------------------------------------------------------------
  # 13. Survival analysis at ROI/patient level
  # ------------------------------------------------------------
  cox_results <- NULL
  cox_summary <- NULL
  km_plots <- list()

  if (all(c("Time", "Event") %in% colnames(roi_factor_scores_wide))) {

    roi_survival_data <- roi_factor_scores_wide %>%
      filter(!is.na(Time), !is.na(Event))

    # Single-factor Cox models.
    cox_results <- lapply(factor_names, function(factor_name) {
      coxph(
        as.formula(paste0("Surv(Time, Event) ~ `", factor_name, "`")),
        data = roi_survival_data
      )
    })

    cox_summary <- map2_dfr(
      cox_results,
      factor_names,
      ~ broom::tidy(.x, conf.int = TRUE, exponentiate = TRUE) %>%
        mutate(factor = .y)
    )

    write.csv(
      cox_summary,
      file.path(output_dir, "cox_summary_ROI_level.csv"),
      row.names = FALSE
    )

    # Kaplan-Meier plots based on median split of ROI-level factor score.
    # Median split is exploratory and should be interpreted with caution.
    for (factor_name in factor_names) {

      tmp <- roi_survival_data %>%
        mutate(
          Factor_group = ifelse(
            .data[[factor_name]] > median(.data[[factor_name]], na.rm = TRUE),
            "High",
            "Low"
          ),
          Factor_group = factor(Factor_group, levels = c("Low", "High"))
        )

      fit <- survfit(
        Surv(Time, Event) ~ Factor_group,
        data = tmp
      )

      km_plot <- ggsurvplot(
        fit,
        data = tmp,
        pval = TRUE,
        title = paste("Survival analysis for", factor_name),
        xlab = "Time",
        ylab = "Survival probability",
        legend.title = factor_name,
        legend.labs = c("Low", "High")
      )

      km_plots[[factor_name]] <- km_plot$plot

      safe_factor_name <- gsub("[^A-Za-z0-9_]", "_", factor_name)

      ggsave(
        file.path(output_dir, paste0("km_ROI_level_", safe_factor_name, ".png")),
        km_plot$plot,
        width = 6,
        height = 5,
        dpi = 300
      )
    }
  } else {
    warning("Time and/or Event columns not found. Survival analysis skipped.")
  }

  # ------------------------------------------------------------
  # 14. Return outputs
  # ------------------------------------------------------------
  list(
    mofa_model = sofa_trained,
    binned_data = sofa_data,
    factor_scores_long = factor_scores_long,
    factor_scores_wide_bin = factor_scores_wide_bin,
    roi_factor_scores_long = roi_factor_scores_long,
    roi_factor_scores_wide = roi_factor_scores_wide,
    factor_names = factor_names,
    spatial_factor_plots = spatial_factor_plots,
    cox_results = cox_results,
    cox_summary = cox_summary,
    km_plots = km_plots,
    output_files = list(
      hdf5 = hdf5_output_path,
      factor_scores_long_csv = file.path(output_dir, "factor_scores_long.csv"),
      factor_scores_wide_bin_csv = file.path(output_dir, "factor_scores_wide_bin_level.csv"),
      roi_factor_scores_long_csv = file.path(output_dir, "roi_factor_scores_long.csv"),
      roi_factor_scores_wide_csv = file.path(output_dir, "roi_factor_scores_wide_mean.csv"),
      cox_summary_ROI_level_csv = file.path(output_dir, "cox_summary_ROI_level.csv"),
      scatter_by_ROI = file.path(output_dir, "scatter_plot_by_ROI.png"),
      cell_count_by_ROI = file.path(output_dir, "cell_count_by_ROI.png"),
      cells_per_roi = file.path(output_dir, "cells_per_roi.png"),
      mean_factor_scores_by_event_ROI_level = file.path(output_dir, "mean_factor_scores_by_event_ROI_level.png"),
      spatial_heterogeneity_by_event_ROI_level = file.path(output_dir, "factor_spatial_heterogeneity_by_event_ROI_level.png"),
      variance_explained = file.path(output_dir, "variance_explained.png"),
      smoothness = file.path(output_dir, "factor_smoothness.png"),
      factors_vs_cov = file.path(output_dir, "factors_vs_cov.png"),
      weights_dir = weights_dir,
      spatial_factor_maps_by_ROI = spatial_factor_dir,
      heatmap = file.path(output_dir, "weights_heatmap.pdf")
    )
  )
}

# ------------------------------------------------------------
# Example run
# ------------------------------------------------------------
# Update paths before running on a different dataset or bin size.

results <- run_sofa(
  coord_table_path = "/mnt/iusers01/dw01/p09731dk/scratch/sofa/coord_table_basel_100cores_newname.csv",
  omic_data_path = "/mnt/iusers01/dw01/p09731dk/scratch/sofa/omic_data_basel_100cores_newname.csv",
  metadata_path = "/mnt/iusers01/dw01/p09731dk/scratch/sofa/metadata_table_basel_100cores_newname.csv",
  output_dir = "/mnt/iusers01/dw01/p09731dk/scratch/sofa/100cores_bin500_v2/",
  hdf5_output = "sofa_basel_100cores_bin500_v2.hdf5",
  normalise = TRUE,
  scale = TRUE,
  seed = 123,
  bin_size = 500
)
