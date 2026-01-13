#' Test Community-Level Mean-Variability Association
#'
#' @description This function tests community-level mean-variability association
#' from a sccomp result. It extracts and analyzes the relationship between mean
#' abundance and variability across all cell groups.
#'
#' @param .data A tibble. The result of sccomp_estimate.
#' @param contrasts A vector of character strings. For example if your formula is
#'   `~ 0 + treatment` and the factor treatment has values `yes` and `no`, your
#'   contrast could be "contrasts = c(treatmentyes - treatmentno)".
#' @param prob_threshold A real between 0 and 1. Posterior probability threshold
#'   for determining significance (default 0.95).
#' @param pass_fit A boolean. Whether to pass the Stan fit as attribute in the output.
#'   Because the Stan fit can be very large, setting this to FALSE can be used to
#'   lower the memory imprint to save the output.
#'
#' @return A tibble (`tbl`) with the following columns:
#' \itemize{
#'   \item parameter - The parameter from the design matrix (e.g., "(Intercept)", "DiseaseStateCRC")
#'   \item component - Component number (1 or 2) for bimodal models; not present for single models
#'   \item factor - The covariate factor in the formula, if applicable
#'   \item v_intercept_lower - Lower (2.5%) quantile of community-level variability intercept
#'   \item v_intercept_effect - Mean of posterior distribution for community intercept
#'   \item v_intercept_upper - Upper (97.5%) quantile of community-level variability intercept
#'   \item v_intercept_prob_gt_0 - Posterior probability that intercept > 0
#'   \item v_intercept_prob_lt_0 - Posterior probability that intercept < 0
#'   \item v_intercept_prob_direction - max(prob_gt_0, prob_lt_0)
#'   \item v_intercept_direction - "increase", "decrease", or "none"
#'   \item v_intercept_significant - TRUE if prob_direction > prob_threshold
#'   \item v_slope_lower - Lower (2.5%) quantile of community-level variability slope
#'   \item v_slope_effect - Mean of posterior distribution for community slope
#'   \item v_slope_upper - Upper (97.5%) quantile of community-level variability slope
#'   \item v_slope_prob_gt_0 - Posterior probability that slope > 0
#'   \item v_slope_prob_lt_0 - Posterior probability that slope < 0
#'   \item v_slope_prob_direction - max(prob_gt_0, prob_lt_0)
#'   \item v_slope_direction - "increase", "decrease", or "none"
#'   \item v_slope_significant - TRUE if prob_direction > prob_threshold
#'   \item rhat_intercept, ess_bulk_intercept, ess_tail_intercept - Convergence diagnostics for intercept
#'   \item rhat_slope, ess_bulk_slope, ess_tail_slope - Convergence diagnostics for slope
#' }
#'
#' @export
#'
#' @examples
#' \donttest{
#'   if (instantiate::stan_cmdstan_exists()) {
#'     data("counts_obj")
#'
#'     estimates = sccomp_estimate(
#'       counts_obj,
#'       ~ 0 + type, ~1, "sample", "cell_group", "count",
#'       cores = 1
#'     ) |>
#'     sccomp_test_community()
#'   }
#' }
#'
sccomp_test_community <- function(.data,
                                  contrasts = NULL,
                                  prob_threshold = 0.95,
                                  pass_fit = TRUE) {

  # Get metadata
  model_input = .data |> attr("model_input")
  noise_model = .data |> attr("noise_model")

  # Check if bimodal (multi_beta_binomial = bimodal model)
  is_bimodal = !is.null(noise_model) && noise_model == "multi_beta_binomial"

  # Get community-level slope draws
  community_slope_draws = get_community_slope_draws(.data, contrasts, is_bimodal)

  # If no parameters match, return empty
  if (nrow(community_slope_draws) == 0) {
    warning("No community-level parameters found matching the contrasts.")
    return(tibble())
  }

  # Calculate statistics
  result = community_slope_draws |>
    community_draws_to_statistics(prob_threshold, is_bimodal)

  # Add factor labels if available
  if ("factor" %in% colnames(model_input$factor_parameter_dictionary)) {
    factor_dict = model_input$factor_parameter_dictionary |>
      select(factor, design_matrix_col)

    result = result |>
      left_join(factor_dict, by = c("parameter" = "design_matrix_col")) |>
      select(parameter, factor, everything())
  }

  # Add attributes
  if (pass_fit) {
    result = result |> add_attr(.data |> attr("fit"), "fit")
  }

  result |>
    add_attr(.data |> attr("model_input"), "model_input") |>
    add_attr(.data |> attr(".sample"), ".sample") |>
    add_attr(.data |> attr(".cell_group"), ".cell_group") |>
    add_attr(.data |> attr(".count"), ".count") |>
    add_attr(.data |> attr("formula_composition"), "formula_composition") |>
    add_attr(.data |> attr("formula_variability"), "formula_variability") |>
    add_attr(.data |> attr("inference_method"), "inference_method") |>
    add_attr(.data |> attr("noise_model"), "noise_model") |>
    add_attr(prob_threshold, "prob_threshold") |>
    add_class("sccomp_community_tbl")
}

#' Extract Community-Level Slope Draws
#'
#' @description Helper function to extract prec_coeff draws from Stan fit and
#' reshape them for community-level analysis.
#'
#' @param .data A sccomp result object
#' @param contrasts Optional contrast specifications
#' @param is_bimodal Logical indicating if the model is bimodal
#'
#' @return A tibble with draws for intercept and slope parameters
#'
#' @noRd
get_community_slope_draws <- function(.data, contrasts, is_bimodal) {

  variability_factors = .data |> attr("model_input") %$% XA |> colnames()
  n_params = length(variability_factors)

  if (is_bimodal) {
    fit = .data |> attr("fit")
    draws_raw = as_draws_df(fit$draws("prec_coeff"))

    # Convert to long format
    draws_long = draws_raw |>
      pivot_longer(
        cols = starts_with("prec_coeff"),
        names_to = "param_name",
        values_to = "value"
      ) |>
      mutate(
        # Extract indices from prec_coeff[C, a]
        C = as.integer(str_extract(param_name, "(?<=\\[)\\d+")),
        a = as.integer(str_extract(param_name, "(?<=,)\\d+(?=\\])")),

        # Map C to component and type
        # C=1: i1 (component=1, intercept)
        # C=2: s1 (component=1, slope)
        # C=3: s2 (component=2, slope)
        # C=4: i2 (component=2, intercept)
        component = case_when(
          C %in% c(1, 2) ~ 1,
          C %in% c(3, 4) ~ 2
        ),
        param_type = case_when(
          C %in% c(1, 4) ~ "intercept",
          C %in% c(2, 3) ~ "slope"
        ),

        # a is the parameter index
        param_idx = a
      ) |>
      select(-param_name, -C, -a)

    # Pivot to have intercept and slope as columns
    draws = draws_long |>
      pivot_wider(
        names_from = param_type,
        values_from = value
      )

    # Add parameter names
    draws = draws |>
      left_join(
        variability_factors |> enframe(name = "param_idx", value = "parameter"),
        by = "param_idx"
      ) |>
      select(parameter, component, param_idx, .chain, .iteration, .draw, intercept, slope)

    # Get convergence diagnostics
    fit_summary = fit$summary("prec_coeff")

    convergence_df = fit_summary |>
      as_tibble() |>
      mutate(
        C = as.integer(str_extract(variable, "(?<=\\[)\\d+")),
        a = as.integer(str_extract(variable, "(?<=,)\\d+(?=\\])")),
        component = case_when(
          C %in% c(1, 2) ~ 1,
          C %in% c(3, 4) ~ 2
        ),
        param_type = case_when(
          C %in% c(1, 4) ~ "intercept",
          C %in% c(2, 3) ~ "slope"
        ),
        param_idx = a
      ) |>
      select(param_idx, component, param_type, any_of(c("rhat", "ess_bulk", "ess_tail"))) |>
      pivot_wider(
        names_from = param_type,
        values_from = any_of(c("rhat", "ess_bulk", "ess_tail")),
        names_sep = "_"
      )

    draws = draws |>
      left_join(convergence_df, by = c("param_idx", "component")) |>
      select(-param_idx)

  } else {
    # Single model case - keep unchanged
    fit = .data |> attr("fit")
    draws_raw = as_draws_df(fit$draws("prec_coeff"))

    draws_long = draws_raw |>
      pivot_longer(
        cols = starts_with("prec_coeff"),
        names_to = "param_name",
        values_to = "value"
      ) |>
      mutate(
        C = as.integer(str_extract(param_name, "(?<=\\[)\\d+")),
        a = as.integer(str_extract(param_name, "(?<=,)\\d+(?=\\])")),
        param_idx = a,
        param_type = if_else(C == 1, "intercept", "slope")
      ) |>
      select(-param_name, -C)

    draws = draws_long |>
      pivot_wider(
        names_from = param_type,
        values_from = value
      )

    draws = draws |>
      left_join(
        variability_factors |> enframe(name = "param_idx", value = "parameter"),
        by = "param_idx"
      ) |>
      select(parameter, param_idx, .chain, .iteration, .draw, intercept, slope)

    fit_summary = fit$summary("prec_coeff")

    convergence_df = fit_summary |>
      as_tibble() |>
      mutate(
        C = as.integer(str_extract(variable, "(?<=\\[)\\d+")),
        a = as.integer(str_extract(variable, "(?<=,)\\d+(?=\\])")),
        param_idx = a,
        param_type = if_else(C == 1, "intercept", "slope")
      ) |>
      select(param_idx, param_type, any_of(c("rhat", "ess_bulk", "ess_tail"))) |>
      pivot_wider(
        names_from = param_type,
        values_from = any_of(c("rhat", "ess_bulk", "ess_tail")),
        names_sep = "_"
      )

    draws = draws |>
      left_join(convergence_df, by = "param_idx") |>
      select(-param_idx)
  }

  # Filter by contrasts
  if (!is.null(contrasts)) {
    contrast_params = contrasts_to_parameter_list(contrasts)
    draws = draws |>
      filter(parameter %in% contrast_params)
  }

  return(draws)
}

#' Calculate Statistics for Community-Level Parameters
#'
#' @description Helper function to calculate posterior statistics from draws
#'
#' @param draws Tibble of posterior draws
#' @param prob_threshold Probability threshold for significance
#' @param is_bimodal Logical indicating if the model is bimodal
#'
#' @return A tibble with summary statistics
#'
#' @noRd
community_draws_to_statistics <- function(draws, prob_threshold = 0.95, is_bimodal) {

  if (is_bimodal) {
    # For bimodal, we want one row per (parameter × component) combination

    # Calculate statistics for intercept
    intercept_stats = draws |>
      group_by(parameter, component) |>
      summarise(
        v_intercept_lower = quantile(intercept, 0.025),
        v_intercept_effect = mean(intercept),
        v_intercept_upper = quantile(intercept, 0.975),
        v_intercept_prob_gt_0 = mean(intercept > 0),
        v_intercept_prob_lt_0 = mean(intercept < 0),
        .groups = "drop"
      ) |>
      mutate(
        v_intercept_prob_direction = pmax(v_intercept_prob_gt_0, v_intercept_prob_lt_0),
        v_intercept_direction = case_when(
          v_intercept_prob_gt_0 > v_intercept_prob_lt_0 & v_intercept_prob_direction > prob_threshold ~ "increase",
          v_intercept_prob_lt_0 > v_intercept_prob_gt_0 & v_intercept_prob_direction > prob_threshold ~ "decrease",
          TRUE ~ "none"
        ),
        v_intercept_significant = v_intercept_prob_direction > prob_threshold
      )

    # Calculate statistics for slope
    slope_stats = draws |>
      group_by(parameter, component) |>
      summarise(
        v_slope_lower = quantile(slope, 0.025),
        v_slope_effect = mean(slope),
        v_slope_upper = quantile(slope, 0.975),
        v_slope_prob_gt_0 = mean(slope > 0),
        v_slope_prob_lt_0 = mean(slope < 0),
        .groups = "drop"
      ) |>
      mutate(
        v_slope_prob_direction = pmax(v_slope_prob_gt_0, v_slope_prob_lt_0),
        v_slope_direction = case_when(
          v_slope_prob_gt_0 > v_slope_prob_lt_0 & v_slope_prob_direction > prob_threshold ~ "increase",
          v_slope_prob_lt_0 > v_slope_prob_gt_0 & v_slope_prob_direction > prob_threshold ~ "decrease",
          TRUE ~ "none"
        ),
        v_slope_significant = v_slope_prob_direction > prob_threshold
      )

    # Get convergence diagnostics - aggregate by parameter and component
    convergence_intercept = draws |>
      group_by(parameter, component) |>
      summarise(
        across(matches("(rhat|ess_bulk|ess_tail)_intercept$"), first),
        .groups = "drop"
      )

    convergence_slope = draws |>
      group_by(parameter, component) |>
      summarise(
        across(matches("(rhat|ess_bulk|ess_tail)_slope$"), first),
        .groups = "drop"
      )

    # Combine all together
    result = intercept_stats |>
      left_join(slope_stats, by = c("parameter", "component")) |>
      left_join(convergence_intercept, by = c("parameter", "component")) |>
      left_join(convergence_slope, by = c("parameter", "component")) |>
      select(
        parameter,
        component,
        v_intercept_lower:v_intercept_significant,
        v_slope_lower:v_slope_significant,
        everything()
      )

  } else {
    # Single model case - one row per parameter

    # Calculate statistics for intercept
    intercept_stats = draws |>
      group_by(parameter) |>
      summarise(
        v_intercept_lower = quantile(intercept, 0.025),
        v_intercept_effect = mean(intercept),
        v_intercept_upper = quantile(intercept, 0.975),
        v_intercept_prob_gt_0 = mean(intercept > 0),
        v_intercept_prob_lt_0 = mean(intercept < 0),
        .groups = "drop"
      ) |>
      mutate(
        v_intercept_prob_direction = pmax(v_intercept_prob_gt_0, v_intercept_prob_lt_0),
        v_intercept_direction = case_when(
          v_intercept_prob_gt_0 > v_intercept_prob_lt_0 & v_intercept_prob_direction > prob_threshold ~ "increase",
          v_intercept_prob_lt_0 > v_intercept_prob_gt_0 & v_intercept_prob_direction > prob_threshold ~ "decrease",
          TRUE ~ "none"
        ),
        v_intercept_significant = v_intercept_prob_direction > prob_threshold
      )

    # Calculate statistics for slope
    slope_stats = draws |>
      group_by(parameter) |>
      summarise(
        v_slope_lower = quantile(slope, 0.025),
        v_slope_effect = mean(slope),
        v_slope_upper = quantile(slope, 0.975),
        v_slope_prob_gt_0 = mean(slope > 0),
        v_slope_prob_lt_0 = mean(slope < 0),
        .groups = "drop"
      ) |>
      mutate(
        v_slope_prob_direction = pmax(v_slope_prob_gt_0, v_slope_prob_lt_0),
        v_slope_direction = case_when(
          v_slope_prob_gt_0 > v_slope_prob_lt_0 & v_slope_prob_direction > prob_threshold ~ "increase",
          v_slope_prob_lt_0 > v_slope_prob_gt_0 & v_slope_prob_direction > prob_threshold ~ "decrease",
          TRUE ~ "none"
        ),
        v_slope_significant = v_slope_prob_direction > prob_threshold
      )

    # Get convergence diagnostics
    convergence_intercept = draws |>
      group_by(parameter) |>
      summarise(
        across(matches("(rhat|ess_bulk|ess_tail)_intercept$"), first),
        .groups = "drop"
      )

    convergence_slope = draws |>
      group_by(parameter) |>
      summarise(
        across(matches("(rhat|ess_bulk|ess_tail)_slope$"), first),
        .groups = "drop"
      )

    # Combine all together
    result = intercept_stats |>
      left_join(slope_stats, by = "parameter") |>
      left_join(convergence_intercept, by = "parameter") |>
      left_join(convergence_slope, by = "parameter") |>
      select(
        parameter,
        v_intercept_lower:v_intercept_significant,
        v_slope_lower:v_slope_significant,
        everything()
      )
  }

  return(result)
}
