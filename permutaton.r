# ================================================================
# Monte Carlo Simulation for Confidential Global Treatment Permutation
# ================================================================
# Objectives:
#   1. Verify exact preservation of statistical information under a
#      single global bijective treatment-label permutation.
#   2. Evaluate posterior disclosure risk for the hidden permutation.
#
# Supported settings:
#   - G = 2, 3, or 4 treatment groups
#   - Continuous or binary outcomes
#   - Equal or unequal randomization
#   - One or multiple interim analyses
#   - External-information scenarios E0, E1, E2, and E3
#
# Required packages:
#   install.packages(c("MASS", "dplyr", "tidyr", "purrr", "tibble"))
# ================================================================

suppressPackageStartupMessages({
  library(MASS)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
})

# ----------------------------------------------------------------
# 1. General helper functions
# ----------------------------------------------------------------

expit <- function(x) {
  1 / (1 + exp(-x))
}

log_sum_exp <- function(x) {
  z <- max(x)
  z + log(sum(exp(x - z)))
}

normalize_log_weights <- function(log_w) {
  exp(log_w - log_sum_exp(log_w))
}

ar1_covariance <- function(p, rho_x) {
  outer(seq_len(p), seq_len(p), function(j, k) rho_x^abs(j - k))
}

validate_probabilities <- function(rho, G) {
  if (length(rho) != G) stop("rho must have length G.")
  if (any(rho <= 0)) stop("All randomization probabilities must be positive.")
  if (abs(sum(rho) - 1) > 1e-10) stop("Randomization probabilities must sum to 1.")
  invisible(TRUE)
}

make_interim_sizes <- function(n, J) {
  unique(pmax(1L, floor(seq_len(J) * n / J)))
}

# Generate all permutations of 1, ..., G without an extra package.
all_permutations <- function(x) {
  if (length(x) == 1L) return(matrix(x, nrow = 1L))
  do.call(
    rbind,
    lapply(seq_along(x), function(i) {
      cbind(x[i], all_permutations(x[-i]))
    })
  )
}

permutation_id <- function(perm) {
  paste(perm, collapse = "-")
}

inverse_permutation <- function(perm) {
  inv <- integer(length(perm))
  inv[perm] <- seq_along(perm)
  inv
}

apply_permutation <- function(labels, perm) {
  perm[labels]
}

# ----------------------------------------------------------------
# 2. Data-generating mechanisms
# ----------------------------------------------------------------

generate_covariates <- function(n, p = 5, rho_x = 0.3) {
  Sigma_x <- ar1_covariance(p, rho_x)
  X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma_x)
  colnames(X) <- paste0("X", seq_len(p))
  X
}

make_treatment_effects <- function(G, delta) {
  # mu = (0, delta, 2*delta, ..., (G-1)*delta)
  (seq_len(G) - 1) * delta
}

generate_trial <- function(
    n,
    G,
    rho = rep(1 / G, G),
    p = 5,
    rho_x = 0.3,
    outcome_type = c("continuous", "binary"),
    delta = 0.5,
    sigma = 1,
    baseline_event = 0.30,
    beta = rep(0.5 / sqrt(p), p),
    seed = NULL) {

  outcome_type <- match.arg(outcome_type)
  validate_probabilities(rho, G)
  if (!is.null(seed)) set.seed(seed)

  A <- sample.int(G, size = n, replace = TRUE, prob = rho)
  X <- generate_covariates(n, p, rho_x)

  if (length(beta) != p) stop("beta must have length p.")
  eta_x <- drop(X %*% beta)
  trt_effect <- make_treatment_effects(G, delta)

  if (outcome_type == "continuous") {
    Y <- trt_effect[A] + eta_x + rnorm(n, mean = 0, sd = sigma)
    event_prob <- rep(NA_real_, n)
  } else {
    alpha <- qlogis(baseline_event)
    linear_predictor <- alpha + trt_effect[A] + eta_x
    event_prob <- expit(linear_predictor)
    Y <- rbinom(n, size = 1, prob = event_prob)
  }

  dat <- data.frame(
    id = seq_len(n),
    A = factor(A, levels = seq_len(G)),
    Y = Y,
    X,
    check.names = FALSE
  )

  attr(dat, "event_probability") <- event_prob
  attr(dat, "treatment_effects") <- trt_effect
  attr(dat, "rho") <- rho
  dat
}

release_trial <- function(dat, perm) {
  A_num <- as.integer(as.character(dat$A))
  R_num <- apply_permutation(A_num, perm)

  dat %>%
    mutate(R = factor(R_num, levels = levels(A))) %>%
    select(id, Y, R, starts_with("X"), A)
}

# ----------------------------------------------------------------
# 3. Original and released statistical analyses
# ----------------------------------------------------------------

fit_trial_model <- function(dat, label_var, outcome_type) {
  x_names <- grep("^X", names(dat), value = TRUE)
  rhs <- paste(c(label_var, x_names), collapse = " + ")
  form <- as.formula(paste("Y ~", rhs))

  if (outcome_type == "continuous") {
    lm(form, data = dat)
  } else {
    glm(form, data = dat, family = binomial())
  }
}

omnibus_treatment_test <- function(dat, label_var, outcome_type) {
  x_names <- grep("^X", names(dat), value = TRUE)
  full_rhs <- paste(c(label_var, x_names), collapse = " + ")
  reduced_rhs <- paste(x_names, collapse = " + ")

  full_form <- as.formula(paste("Y ~", full_rhs))
  reduced_form <- as.formula(paste("Y ~", reduced_rhs))

  if (outcome_type == "continuous") {
    fit_full <- lm(full_form, data = dat)
    fit_reduced <- lm(reduced_form, data = dat)
    tab <- anova(fit_reduced, fit_full)

    list(
      statistic = unname(tab$F[2]),
      p_value = unname(tab$`Pr(>F)`[2]),
      logLik = as.numeric(logLik(fit_full)),
      deviance = deviance(fit_full),
      fitted = fitted(fit_full),
      residuals = residuals(fit_full),
      model = fit_full
    )
  } else {
    fit_full <- glm(full_form, data = dat, family = binomial())
    fit_reduced <- glm(reduced_form, data = dat, family = binomial())
    tab <- anova(fit_reduced, fit_full, test = "LRT")

    list(
      statistic = unname(tab$Deviance[2]),
      p_value = unname(tab$`Pr(>Chi)`[2]),
      logLik = as.numeric(logLik(fit_full)),
      deviance = deviance(fit_full),
      fitted = fitted(fit_full),
      residuals = residuals(fit_full, type = "response"),
      model = fit_full
    )
  }
}

# Recover released coefficients to the original treatment-label ordering.
# The model matrix columns differ only by a permutation/reparameterization,
# so fitted values and omnibus tests are the most direct invariants.
recover_treatment_means_continuous <- function(model, dat, label_var, G) {
  x_names <- grep("^X", names(dat), value = TRUE)
  newdat <- as.data.frame(matrix(0, nrow = G, ncol = length(x_names)))
  names(newdat) <- x_names
  newdat[[label_var]] <- factor(seq_len(G), levels = seq_len(G))
  newdat <- newdat[, c(label_var, x_names), drop = FALSE]
  as.numeric(predict(model, newdata = newdat))
}

compare_original_released <- function(dat_released, perm, outcome_type) {
  G <- nlevels(dat_released$A)

  original <- omnibus_treatment_test(dat_released, "A", outcome_type)
  released <- omnibus_treatment_test(dat_released, "R", outcome_type)

  fitted_diff <- max(abs(original$fitted - released$fitted))
  loglik_diff <- released$logLik - original$logLik
  stat_diff <- released$statistic - original$statistic
  p_diff <- released$p_value - original$p_value

  mean_diff <- NA_real_
  if (outcome_type == "continuous") {
    mu_A <- recover_treatment_means_continuous(original$model, dat_released, "A", G)
    mu_R <- recover_treatment_means_continuous(released$model, dat_released, "R", G)
    mu_R_recovered <- mu_R[perm]  # true g maps to released perm[g]
    mean_diff <- max(abs(mu_R_recovered - mu_A))
  }

  tibble(
    p_original = original$p_value,
    p_released = released$p_value,
    statistic_original = original$statistic,
    statistic_released = released$statistic,
    logLik_original = original$logLik,
    logLik_released = released$logLik,
    delta_p = p_diff,
    delta_statistic = stat_diff,
    delta_logLik = loglik_diff,
    max_abs_fitted_difference = fitted_diff,
    max_abs_recovered_mean_difference = mean_diff,
    reject_original = original$p_value < 0.05,
    reject_released = released$p_value < 0.05
  )
}

# ----------------------------------------------------------------
# 4. Posterior inference for the confidential permutation
# ----------------------------------------------------------------

# External information can include:
#   expected_outcome[g]: expected mean/event probability for true treatment g
#   outcome_sd: assumed SD for released-group sample means (continuous)
#   expected_ae[g]: expected adverse-event probability for true treatment g
#   known_rho[g]: known randomization probability for true treatment g
#
# Candidate permutation pi uses R = pi(A). Thus released group h corresponds
# to true group inverse(pi)[h].

allocation_loglik <- function(counts_by_released, candidate_perm, known_rho) {
  inv <- inverse_permutation(candidate_perm)
  probs_released <- known_rho[inv]
  dmultinom(counts_by_released, prob = probs_released, log = TRUE)
}

outcome_loglik_continuous <- function(dat, candidate_perm, expected_outcome,
                                      assumed_sigma = 1) {
  G <- length(candidate_perm)
  inv <- inverse_permutation(candidate_perm)
  ll <- 0

  for (h in seq_len(G)) {
    y_h <- dat$Y[as.integer(dat$R) == h]
    g <- inv[h]
    if (length(y_h) > 0) {
      ll <- ll + sum(dnorm(y_h, mean = expected_outcome[g],
                           sd = assumed_sigma, log = TRUE))
    }
  }
  ll
}

outcome_loglik_binary <- function(dat, candidate_perm, expected_outcome) {
  G <- length(candidate_perm)
  inv <- inverse_permutation(candidate_perm)
  ll <- 0

  for (h in seq_len(G)) {
    y_h <- dat$Y[as.integer(dat$R) == h]
    g <- inv[h]
    p_g <- min(max(expected_outcome[g], 1e-8), 1 - 1e-8)
    if (length(y_h) > 0) {
      ll <- ll + sum(dbinom(y_h, size = 1, prob = p_g, log = TRUE))
    }
  }
  ll
}

safety_loglik_binary <- function(dat, candidate_perm, expected_ae) {
  if (!"AE" %in% names(dat)) return(0)

  G <- length(candidate_perm)
  inv <- inverse_permutation(candidate_perm)
  ll <- 0

  for (h in seq_len(G)) {
    ae_h <- dat$AE[as.integer(dat$R) == h]
    g <- inv[h]
    p_g <- min(max(expected_ae[g], 1e-8), 1 - 1e-8)
    if (length(ae_h) > 0) {
      ll <- ll + sum(dbinom(ae_h, size = 1, prob = p_g, log = TRUE))
    }
  }
  ll
}

permutation_posterior <- function(
    dat,
    permutations,
    scenario = c("E0", "E1", "E2", "E3"),
    outcome_type,
    known_rho = NULL,
    expected_outcome = NULL,
    assumed_sigma = 1,
    expected_ae = NULL,
    prior = NULL) {

  scenario <- match.arg(scenario)
  K <- nrow(permutations)
  G <- ncol(permutations)

  if (is.null(prior)) prior <- rep(1 / K, K)
  if (length(prior) != K || any(prior <= 0)) stop("Invalid permutation prior.")

  counts <- tabulate(as.integer(dat$R), nbins = G)
  log_post <- log(prior)

  for (k in seq_len(K)) {
    pi_k <- permutations[k, ]

    if (scenario %in% c("E1", "E3")) {
      if (is.null(known_rho)) stop("known_rho is required for E1/E3.")
      log_post[k] <- log_post[k] + allocation_loglik(counts, pi_k, known_rho)
    }

    if (scenario %in% c("E2", "E3")) {
      if (is.null(expected_outcome)) stop("expected_outcome is required for E2/E3.")

      if (outcome_type == "continuous") {
        log_post[k] <- log_post[k] + outcome_loglik_continuous(
          dat, pi_k, expected_outcome, assumed_sigma
        )
      } else {
        log_post[k] <- log_post[k] + outcome_loglik_binary(
          dat, pi_k, expected_outcome
        )
      }
    }

    if (scenario == "E3" && !is.null(expected_ae)) {
      log_post[k] <- log_post[k] + safety_loglik_binary(dat, pi_k, expected_ae)
    }
  }

  posterior <- normalize_log_weights(log_post)

  tibble(
    permutation_index = seq_len(K),
    permutation = apply(permutations, 1, permutation_id),
    posterior = posterior
  )
}

posterior_mapping_matrix <- function(posterior, permutations) {
  G <- ncol(permutations)
  C <- matrix(0, nrow = G, ncol = G,
              dimnames = list(
                released_group = paste0("R", seq_len(G)),
                true_group = paste0("A", seq_len(G))
              ))

  for (k in seq_len(nrow(permutations))) {
    inv <- inverse_permutation(permutations[k, ])
    for (h in seq_len(G)) {
      C[h, inv[h]] <- C[h, inv[h]] + posterior[k]
    }
  }
  C
}

privacy_measures <- function(posterior_tbl, permutations, true_perm,
                             interim_look, n_interim) {
  p <- posterior_tbl$posterior
  G <- ncol(permutations)
  K <- nrow(permutations)

  entropy <- -sum(ifelse(p > 0, p * log(p), 0))
  d_max <- max(p)
  max_set <- which(abs(p - max(p)) < 1e-12)
  map_index <- sample(max_set, size = 1L)
  true_index <- which(apply(permutations, 1, function(x) all(as.integer(x) == as.integer(true_perm))))
  C <- posterior_mapping_matrix(p, permutations)
  group_disclosure <- apply(C, 1, max)

  tibble(
    interim_look = interim_look,
    n_interim = n_interim,
    max_posterior = d_max,
    classification_privacy = (1 - d_max) / (1 - 1 / K),
    entropy = entropy,
    entropy_privacy = if (K > 1) entropy / log(K) else 1,
    effective_permutations = exp(entropy),
    mutual_information = log(K) - entropy,
    normalized_information_disclosure = if (K > 1) 1 - entropy / log(K) else 0,
    map_recovery = as.integer(map_index == true_index),
    true_permutation_posterior = p[true_index],
    max_group_disclosure = max(group_disclosure),
    mean_group_disclosure = mean(group_disclosure)
  )
}

# ----------------------------------------------------------------
# 5. One Monte Carlo replication
# ----------------------------------------------------------------

simulate_one_replication <- function(
    replication,
    n,
    G,
    rho,
    p,
    rho_x,
    outcome_type,
    delta,
    sigma,
    baseline_event,
    J,
    external_scenario,
    efficacy_information_strength = 1,
    safety_event_prob = NULL,
    seed = 20260723) {

  set.seed(seed + replication)

  dat <- generate_trial(
    n = n,
    G = G,
    rho = rho,
    p = p,
    rho_x = rho_x,
    outcome_type = outcome_type,
    delta = delta,
    sigma = sigma,
    baseline_event = baseline_event
  )

  permutations <- all_permutations(seq_len(G))
  true_perm_index <- sample.int(nrow(permutations), 1)
  true_perm <- permutations[true_perm_index, ]
  released <- release_trial(dat, true_perm)

  # Optional safety outcome for E3.
  if (!is.null(safety_event_prob)) {
    A_num <- as.integer(as.character(released$A))
    released$AE <- rbinom(n, 1, safety_event_prob[A_num])
  }

  preservation <- compare_original_released(released, true_perm, outcome_type) %>%
    mutate(
      replication = replication,
      true_permutation = permutation_id(true_perm),
      .before = 1
    )

  # External efficacy information.
  true_effects <- make_treatment_effects(G, delta)
  if (outcome_type == "continuous") {
    expected_outcome <- efficacy_information_strength * true_effects
  } else {
    # Marginal approximation at X = 0, as described in the manuscript.
    expected_outcome <- expit(qlogis(baseline_event) +
                                efficacy_information_strength * true_effects)
  }

  interim_sizes <- make_interim_sizes(n, J)

  privacy <- map_dfr(seq_along(interim_sizes), function(j) {
    n_j <- interim_sizes[j]
    dat_j <- released[seq_len(n_j), , drop = FALSE]

    post <- permutation_posterior(
      dat = dat_j,
      permutations = permutations,
      scenario = external_scenario,
      outcome_type = outcome_type,
      known_rho = rho,
      expected_outcome = expected_outcome,
      assumed_sigma = sigma,
      expected_ae = safety_event_prob
    )

    privacy_measures(
      posterior_tbl = post,
      permutations = permutations,
      true_perm = true_perm,
      interim_look = j,
      n_interim = n_j
    ) %>%
      mutate(
        replication = replication,
        true_permutation = permutation_id(true_perm),
        .before = 1
      )
  })

  list(preservation = preservation, privacy = privacy)
}

# ----------------------------------------------------------------
# 6. Monte Carlo driver for one scenario
# ----------------------------------------------------------------

run_simulation_scenario <- function(
    M = 1000,
    n = 500,
    G = 3,
    rho = rep(1 / G, G),
    p = 5,
    rho_x = 0.3,
    outcome_type = c("continuous", "binary"),
    delta = 0.5,
    sigma = 1,
    baseline_event = 0.30,
    J = 3,
    external_scenario = c("E0", "E1", "E2", "E3"),
    efficacy_information_strength = 1,
    safety_event_prob = NULL,
    seed = 20260723,
    progress = TRUE) {

  outcome_type <- match.arg(outcome_type)
  external_scenario <- match.arg(external_scenario)
  validate_probabilities(rho, G)

  if (progress) {
    message("Running scenario: G=", G,
            ", n=", n,
            ", outcome=", outcome_type,
            ", delta=", delta,
            ", J=", J,
            ", external scenario=", external_scenario,
            ", M=", M)
  }

  simulations <- vector("list", M)

  for (m in seq_len(M)) {
    simulations[[m]] <- simulate_one_replication(
      replication = m,
      n = n,
      G = G,
      rho = rho,
      p = p,
      rho_x = rho_x,
      outcome_type = outcome_type,
      delta = delta,
      sigma = sigma,
      baseline_event = baseline_event,
      J = J,
      external_scenario = external_scenario,
      efficacy_information_strength = efficacy_information_strength,
      safety_event_prob = safety_event_prob,
      seed = seed
    )

    if (progress && (m %% max(1, floor(M / 10)) == 0 || m == M)) {
      message("  completed ", m, " / ", M)
    }
  }

  preservation_raw <- bind_rows(map(simulations, "preservation"))
  privacy_raw <- bind_rows(map(simulations, "privacy"))

  scenario_columns <- tibble(
    M = M,
    n = n,
    G = G,
    randomization = paste(round(rho, 4), collapse = ":"),
    p = p,
    rho_x = rho_x,
    outcome_type = outcome_type,
    delta = delta,
    sigma = sigma,
    baseline_event = baseline_event,
    J = J,
    external_scenario = external_scenario,
    efficacy_information_strength = efficacy_information_strength
  )

  preservation_summary <- preservation_raw %>%
    summarise(
      rejection_original = mean(reject_original),
      rejection_released = mean(reject_released),
      rejection_difference = rejection_released - rejection_original,
      discordant_rejection_rate = mean(reject_original != reject_released),
      max_abs_delta_p = max(abs(delta_p), na.rm = TRUE),
      max_abs_delta_statistic = max(abs(delta_statistic), na.rm = TRUE),
      max_abs_delta_logLik = max(abs(delta_logLik), na.rm = TRUE),
      max_abs_fitted_difference = max(max_abs_fitted_difference, na.rm = TRUE),
      max_abs_recovered_mean_difference = if (
        all(is.na(max_abs_recovered_mean_difference))
      ) NA_real_ else max(max_abs_recovered_mean_difference, na.rm = TRUE)
    ) %>%
    bind_cols(scenario_columns, .)

  privacy_summary <- privacy_raw %>%
    group_by(interim_look, n_interim) %>%
    summarise(
      mean_max_posterior = mean(max_posterior),
      mean_classification_privacy = mean(classification_privacy),
      mean_entropy = mean(entropy),
      mean_entropy_privacy = mean(entropy_privacy),
      mean_effective_permutations = mean(effective_permutations),
      mean_mutual_information = mean(mutual_information),
      mean_information_disclosure = mean(normalized_information_disclosure),
      map_recovery_probability = mean(map_recovery),
      mean_true_permutation_posterior = mean(true_permutation_posterior),
      mean_max_group_disclosure = mean(max_group_disclosure),
      mean_group_disclosure = mean(mean_group_disclosure),
      .groups = "drop"
    ) %>%
    mutate(
      M = M,
      n = n,
      G = G,
      randomization = paste(round(rho, 4), collapse = ":"),
      p = p,
      rho_x = rho_x,
      outcome_type = outcome_type,
      delta = delta,
      sigma = sigma,
      baseline_event = baseline_event,
      J = J,
      external_scenario = external_scenario,
      efficacy_information_strength = efficacy_information_strength,
      .before = 1
    )

  list(
    preservation_raw = preservation_raw,
    privacy_raw = privacy_raw,
    preservation_summary = preservation_summary,
    privacy_summary = privacy_summary
  )
}

# ----------------------------------------------------------------
# 7. Scenario-grid driver
# ----------------------------------------------------------------

run_scenario_grid <- function(scenario_grid, M = 1000,
                              seed = 20260723, progress = TRUE) {
  required <- c(
    "n", "G", "rho", "p", "rho_x", "outcome_type", "delta",
    "sigma", "baseline_event", "J", "external_scenario"
  )
  missing_cols <- setdiff(required, names(scenario_grid))
  if (length(missing_cols) > 0) {
    stop("Scenario grid is missing: ", paste(missing_cols, collapse = ", "))
  }

  results <- vector("list", nrow(scenario_grid))

  for (s in seq_len(nrow(scenario_grid))) {
    row <- scenario_grid[s, ]

    rho_s <- row$rho[[1]]
    info_strength <- if ("efficacy_information_strength" %in% names(row)) {
      row$efficacy_information_strength[[1]]
    } else 1
    safety_prob <- if ("safety_event_prob" %in% names(row)) {
      row$safety_event_prob[[1]]
    } else NULL

    results[[s]] <- run_simulation_scenario(
      M = M,
      n = row$n[[1]],
      G = row$G[[1]],
      rho = rho_s,
      p = row$p[[1]],
      rho_x = row$rho_x[[1]],
      outcome_type = row$outcome_type[[1]],
      delta = row$delta[[1]],
      sigma = row$sigma[[1]],
      baseline_event = row$baseline_event[[1]],
      J = row$J[[1]],
      external_scenario = row$external_scenario[[1]],
      efficacy_information_strength = info_strength,
      safety_event_prob = safety_prob,
      seed = seed + 100000 * s,
      progress = progress
    )
  }

  list(
    preservation_summary = bind_rows(map(results, "preservation_summary")),
    privacy_summary = bind_rows(map(results, "privacy_summary")),
    preservation_raw = bind_rows(map(results, "preservation_raw"), .id = "scenario_id"),
    privacy_raw = bind_rows(map(results, "privacy_raw"), .id = "scenario_id")
  )
}

# ----------------------------------------------------------------
# 8. Example analyses
# ----------------------------------------------------------------

# Example 1: quick continuous-outcome check
# Under delta = 0, the rejection rate estimates type-I error.
example_continuous_null <- run_simulation_scenario(
  M = 200,             # use M = 5000 for the manuscript
  n = 500,
  G = 3,
  rho = c(1/3, 1/3, 1/3),
  p = 5,
  rho_x = 0.3,
  outcome_type = "continuous",
  delta = 0,
  sigma = 1,
  J = 3,
  external_scenario = "E0",
  seed = 20260723
)

print(example_continuous_null$preservation_summary)
print(example_continuous_null$privacy_summary)

# Example 2: continuous-outcome power and efficacy-based disclosure
example_continuous_alt <- run_simulation_scenario(
  M = 200,
  n = 500,
  G = 3,
  rho = c(1/3, 1/3, 1/3),
  p = 5,
  rho_x = 0.3,
  outcome_type = "continuous",
  delta = 0.5,
  sigma = 1,
  J = 5,
  external_scenario = "E2",
  efficacy_information_strength = 1,
  seed = 20260723
)

print(example_continuous_alt$preservation_summary)
print(example_continuous_alt$privacy_summary)

# Example 3: binary outcome with unequal allocation
example_binary <- run_simulation_scenario(
  M = 200,
  n = 500,
  G = 3,
  rho = c(0.4, 0.4, 0.2),       # 2:2:1 allocation
  p = 5,
  rho_x = 0.3,
  outcome_type = "binary",
  delta = log(1.50),
  baseline_event = 0.30,
  J = 3,
  external_scenario = "E1",
  seed = 20260723
)

print(example_binary$preservation_summary)
print(example_binary$privacy_summary)

# ----------------------------------------------------------------
# 9. Manuscript-sized scenario grid
# ----------------------------------------------------------------

# This is a compact example. Expand it as needed.
scenario_grid <- tribble(
  ~n, ~G, ~rho,                  ~p, ~rho_x, ~outcome_type, ~delta,       ~sigma, ~baseline_event, ~J, ~external_scenario,
  100, 2, list(c(0.5, 0.5)),       5, 0.3,    "continuous",  0,            1,      0.30,            1,  "E0",
  500, 2, list(c(0.5, 0.5)),       5, 0.3,    "continuous",  0.5,          1,      0.30,            3,  "E2",
  500, 2, list(c(2/3, 1/3)),       5, 0.3,    "continuous",  0.5,          1,      0.30,            3,  "E1",
  500, 3, list(c(1/3, 1/3, 1/3)), 5, 0.3,    "continuous",  0.5,          1,      0.30,            5,  "E0",
  500, 3, list(c(0.4, 0.4, 0.2)), 5, 0.3,    "binary",      log(1.50),    1,      0.30,            3,  "E3",
 1000, 4, list(rep(0.25, 4)),       10, 0.6,   "binary",      log(1.25),    1,      0.30,            5,  "E2"
)

# Add optional E3 safety information as a list-column.
scenario_grid <- scenario_grid %>%
  mutate(
    efficacy_information_strength = 1,
    safety_event_prob = map(G, function(g) seq(0.10, 0.10 + 0.05 * (g - 1), length.out = g))
  )

# Full manuscript run:
# grid_results <- run_scenario_grid(
#   scenario_grid = scenario_grid,
#   M = 5000,
#   seed = 20260723,
#   progress = TRUE
# )
#
# write.csv(grid_results$preservation_summary,
#           "preservation_summary.csv", row.names = FALSE)
# write.csv(grid_results$privacy_summary,
#           "privacy_summary.csv", row.names = FALSE)
# write.csv(grid_results$preservation_raw,
#           "preservation_raw.csv", row.names = FALSE)
# write.csv(grid_results$privacy_raw,
#           "privacy_raw.csv", row.names = FALSE)

# ================================================================
# End of file
# ================================================================

# ================================================================
# 10. Additional existing regrouping methods
# ================================================================
# The methods below extend the confidential global-permutation
# simulation to include:
#   1. Pairwise global label swapping.
#   2. Exact inverse decoding of deterministic permutations.
#   3. Participant-level stochastic Markov regrouping.
#   4. Repeated applications of a Markov transition matrix.
#   5. Transition-matrix correction of marginal treatment means.
#   6. Power/type-I-error comparisons across release mechanisms.
#
# A deterministic global permutation is only a relabeling and therefore
# preserves the fitted model space and the omnibus treatment test exactly.
# In contrast, participant-level stochastic regrouping mixes true treatment
# groups. Its naive released-label analysis generally loses power. When the
# transition matrix and randomization probabilities are known, marginal
# treatment means can be recovered by matrix inversion, provided the induced
# mixing matrix is nonsingular and sufficiently well conditioned.

# ----------------------------------------------------------------
# 10.1 Pairwise global swapping and exact decoding
# ----------------------------------------------------------------

validate_permutation_vector <- function(perm, G = length(perm)) {
  perm <- as.integer(perm)
  if (length(perm) != G || !identical(sort(perm), seq_len(G))) {
    stop("perm must be a one-to-one permutation of 1, ..., G.")
  }
  invisible(TRUE)
}

pairwise_swap_permutation <- function(G, first, second) {
  if (length(first) != 1L || length(second) != 1L ||
      first == second ||
      !all(c(first, second) %in% seq_len(G))) {
    stop("first and second must be distinct indices in 1, ..., G.")
  }

  perm <- seq_len(G)
  perm[c(first, second)] <- perm[c(second, first)]
  perm
}

compose_permutation_vectors <- function(first_perm, second_perm) {
  # Apply first_perm first and second_perm second:
  # g -> first_perm[g] -> second_perm[first_perm[g]].
  G <- length(first_perm)
  validate_permutation_vector(first_perm, G)
  validate_permutation_vector(second_perm, G)
  second_perm[first_perm]
}

permutation_from_pairwise_swaps <- function(G, swaps) {
  # swaps is a two-column matrix/data frame. Each row gives one
  # transposition to be applied sequentially.
  if (is.null(swaps) || length(swaps) == 0L) return(seq_len(G))

  swaps <- as.matrix(swaps)
  if (ncol(swaps) != 2L) {
    stop("swaps must have exactly two columns.")
  }

  current_perm <- seq_len(G)

  for (j in seq_len(nrow(swaps))) {
    transposition <- pairwise_swap_permutation(
      G = G,
      first = swaps[j, 1L],
      second = swaps[j, 2L]
    )
    current_perm <- compose_permutation_vectors(
      current_perm,
      transposition
    )
  }

  current_perm
}

release_trial_pairwise <- function(dat, swaps) {
  G <- nlevels(dat$A)
  perm <- permutation_from_pairwise_swaps(G, swaps)
  released <- release_trial(dat, perm)
  attr(released, "pairwise_permutation") <- perm
  attr(released, "pairwise_swaps") <- swaps
  released
}

decode_released_labels <- function(released_labels, perm) {
  validate_permutation_vector(perm)
  inv <- inverse_permutation(perm)
  inv[as.integer(released_labels)]
}

verify_exact_decoding <- function(dat_released, perm) {
  decoded <- decode_released_labels(dat_released$R, perm)
  true_labels <- as.integer(as.character(dat_released$A))

  tibble(
    exact_recovery = all(decoded == true_labels),
    number_incorrect = sum(decoded != true_labels),
    incorrect_fraction = mean(decoded != true_labels)
  )
}

# ----------------------------------------------------------------
# 10.2 Stochastic Markov-kernel regrouping
# ----------------------------------------------------------------

validate_transition_kernel <- function(K, tolerance = 1e-10) {
  K <- as.matrix(K)

  if (nrow(K) != ncol(K)) {
    stop("K must be a square transition matrix.")
  }
  if (any(!is.finite(K)) || any(K < 0)) {
    stop("K must contain finite, nonnegative entries.")
  }
  if (any(abs(rowSums(K) - 1) > tolerance)) {
    stop("Every row of K must sum to one.")
  }

  invisible(TRUE)
}

matrix_power <- function(K, exponent) {
  validate_transition_kernel(K)
  if (length(exponent) != 1L || exponent < 0 || exponent != floor(exponent)) {
    stop("exponent must be a nonnegative integer.")
  }

  G <- nrow(K)
  result <- diag(G)
  base <- K
  exponent <- as.integer(exponent)

  while (exponent > 0L) {
    if (exponent %% 2L == 1L) result <- result %*% base
    base <- base %*% base
    exponent <- exponent %/% 2L
  }

  result
}

make_symmetric_markov_kernel <- function(G, stay_probability = 0.80) {
  if (G < 2L) stop("G must be at least two.")
  if (stay_probability < 0 || stay_probability > 1) {
    stop("stay_probability must lie between zero and one.")
  }

  off_diagonal <- (1 - stay_probability) / (G - 1)
  K <- matrix(off_diagonal, nrow = G, ncol = G)
  diag(K) <- stay_probability
  K
}

release_trial_markov <- function(dat, K, n_transitions = 1L, seed = NULL) {
  validate_transition_kernel(K)
  G <- nlevels(dat$A)
  if (nrow(K) != G) stop("K must have one row and column per treatment group.")
  if (!is.null(seed)) set.seed(seed)

  K_effective <- matrix_power(K, n_transitions)
  A_num <- as.integer(as.character(dat$A))

  R_num <- vapply(
    A_num,
    function(g) sample.int(G, size = 1L, prob = K_effective[g, ]),
    integer(1L)
  )

  released <- dat %>%
    mutate(R = factor(R_num, levels = seq_len(G))) %>%
    select(id, Y, R, starts_with("X"), A, everything())

  attr(released, "transition_kernel") <- K
  attr(released, "effective_transition_kernel") <- K_effective
  attr(released, "n_transitions") <- n_transitions
  released
}

markov_convergence_diagnostics <- function(K, powers = c(1, 2, 5, 10, 25, 50)) {
  validate_transition_kernel(K)
  powers <- sort(unique(as.integer(powers[powers >= 0])))

  eigen_values <- eigen(t(K), only.values = TRUE)$values
  stationary_index <- which.min(abs(eigen_values - 1))
  stationary_vector <- Re(eigen(t(K))$vectors[, stationary_index])
  stationary_vector <- stationary_vector / sum(stationary_vector)
  stationary_matrix <- matrix(
    rep(stationary_vector, each = nrow(K)),
    nrow = nrow(K),
    byrow = FALSE
  )

  map_dfr(powers, function(m) {
    Km <- matrix_power(K, m)
    tibble(
      number_of_transitions = m,
      determinant = det(Km),
      condition_number = kappa(Km),
      max_distance_from_identity = max(abs(Km - diag(nrow(K)))),
      max_distance_from_stationary = max(abs(Km - stationary_matrix))
    )
  })
}

# ----------------------------------------------------------------
# 10.3 Transition-matrix correction for marginal means
# ----------------------------------------------------------------

induced_released_mixing_matrix <- function(K, rho, n_transitions = 1L) {
  validate_transition_kernel(K)
  G <- nrow(K)
  validate_probabilities(rho, G)

  K_effective <- matrix_power(K, n_transitions)
  released_prob <- drop(rho %*% K_effective)

  if (any(released_prob <= 0)) {
    stop("Every released group must have positive probability.")
  }

  # B[h, g] = Pr(A = g | R = h).
  B <- matrix(0, nrow = G, ncol = G)
  for (h in seq_len(G)) {
    for (g in seq_len(G)) {
      B[h, g] <- rho[g] * K_effective[g, h] / released_prob[h]
    }
  }

  dimnames(B) <- list(
    released_group = paste0("R", seq_len(G)),
    true_group = paste0("A", seq_len(G))
  )
  B
}

released_group_moment_estimates <- function(dat_released, G, outcome_type) {
  R_num <- as.integer(dat_released$R)

  estimates <- variances <- sample_sizes <- numeric(G)

  for (h in seq_len(G)) {
    y_h <- dat_released$Y[R_num == h]
    sample_sizes[h] <- length(y_h)

    if (length(y_h) == 0L) {
      estimates[h] <- NA_real_
      variances[h] <- NA_real_
    } else {
      estimates[h] <- mean(y_h)

      if (outcome_type == "continuous") {
        variances[h] <- if (length(y_h) > 1L) {
          var(y_h) / length(y_h)
        } else {
          NA_real_
        }
      } else {
        p_hat <- estimates[h]
        variances[h] <- max(p_hat * (1 - p_hat), 1e-8) / length(y_h)
      }
    }
  }

  list(
    estimate = estimates,
    covariance = diag(variances, nrow = G),
    sample_size = sample_sizes
  )
}

omnibus_wald_from_means <- function(mean_estimate, covariance) {
  G <- length(mean_estimate)
  if (G < 2L) stop("At least two groups are required.")

  C <- cbind(-1, diag(G - 1L))
  contrast <- drop(C %*% mean_estimate)
  contrast_covariance <- C %*% covariance %*% t(C)

  statistic <- tryCatch(
    drop(t(contrast) %*% solve(contrast_covariance, contrast)),
    error = function(e) NA_real_
  )

  tibble(
    statistic = statistic,
    df = G - 1L,
    p_value = if (is.finite(statistic)) {
      pchisq(statistic, df = G - 1L, lower.tail = FALSE)
    } else {
      NA_real_
    }
  )
}

correct_markov_released_means <- function(
    dat_released,
    K,
    rho,
    outcome_type,
    n_transitions = 1L,
    condition_limit = 1e8) {

  G <- nlevels(dat_released$A)
  B <- induced_released_mixing_matrix(K, rho, n_transitions)
  moments <- released_group_moment_estimates(dat_released, G, outcome_type)

  if (any(!is.finite(moments$estimate)) ||
      any(!is.finite(moments$covariance))) {
    return(list(
      status = "Not estimable: an empty or singleton released group",
      mixing_matrix = B,
      condition_number = kappa(B),
      corrected_mean = rep(NA_real_, G),
      corrected_covariance = matrix(NA_real_, G, G),
      omnibus = tibble(statistic = NA_real_, df = G - 1L, p_value = NA_real_)
    ))
  }

  condition_number <- kappa(B)
  if (!is.finite(condition_number) || condition_number > condition_limit) {
    return(list(
      status = "Not estimable: mixing matrix is singular or ill conditioned",
      mixing_matrix = B,
      condition_number = condition_number,
      corrected_mean = rep(NA_real_, G),
      corrected_covariance = matrix(NA_real_, G, G),
      omnibus = tibble(statistic = NA_real_, df = G - 1L, p_value = NA_real_)
    ))
  }

  B_inverse <- solve(B)
  corrected_mean <- drop(B_inverse %*% moments$estimate)
  corrected_covariance <- B_inverse %*%
    moments$covariance %*%
    t(B_inverse)

  list(
    status = "Estimated",
    mixing_matrix = B,
    condition_number = condition_number,
    released_mean = moments$estimate,
    corrected_mean = corrected_mean,
    corrected_covariance = corrected_covariance,
    omnibus = omnibus_wald_from_means(
      corrected_mean,
      corrected_covariance
    )
  )
}

# ----------------------------------------------------------------
# 10.4 One-replication comparison of all regrouping methods
# ----------------------------------------------------------------

compare_regrouping_methods_one_replication <- function(
    replication,
    n,
    G,
    rho,
    p,
    rho_x,
    outcome_type,
    delta,
    sigma,
    baseline_event,
    K = make_symmetric_markov_kernel(G, 0.80),
    n_transitions = 1L,
    alpha = 0.05,
    seed = 20260723) {

  set.seed(seed + replication)

  dat <- generate_trial(
    n = n,
    G = G,
    rho = rho,
    p = p,
    rho_x = rho_x,
    outcome_type = outcome_type,
    delta = delta,
    sigma = sigma,
    baseline_event = baseline_event
  )

  original_test <- omnibus_treatment_test(dat, "A", outcome_type)

  # Existing method: one confidential global permutation.
  permutations <- all_permutations(seq_len(G))
  global_perm <- permutations[sample.int(nrow(permutations), 1L), ]
  global_release <- release_trial(dat, global_perm)
  global_test <- omnibus_treatment_test(global_release, "R", outcome_type)
  global_decode <- verify_exact_decoding(global_release, global_perm)

  # Existing method: sequential pairwise swaps producing a global permutation.
  number_of_swaps <- max(1L, G - 1L)
  swaps <- t(replicate(
    number_of_swaps,
    sample(seq_len(G), size = 2L, replace = FALSE)
  ))
  pairwise_release <- release_trial_pairwise(dat, swaps)
  pairwise_perm <- attr(pairwise_release, "pairwise_permutation")
  pairwise_test <- omnibus_treatment_test(pairwise_release, "R", outcome_type)
  pairwise_decode <- verify_exact_decoding(pairwise_release, pairwise_perm)

  # New stochastic participant-level Markov regrouping.
  markov_release <- release_trial_markov(
    dat = dat,
    K = K,
    n_transitions = n_transitions
  )
  markov_naive_test <- omnibus_treatment_test(
    markov_release,
    "R",
    outcome_type
  )

  markov_corrected <- correct_markov_released_means(
    dat_released = markov_release,
    K = K,
    rho = rho,
    outcome_type = outcome_type,
    n_transitions = n_transitions
  )

  tibble(
    replication = replication,
    method = c(
      "Original labels",
      "Global permutation",
      "Pairwise global swapping",
      "Markov naive released labels",
      "Markov matrix corrected"
    ),
    statistic = c(
      original_test$statistic,
      global_test$statistic,
      pairwise_test$statistic,
      markov_naive_test$statistic,
      markov_corrected$omnibus$statistic
    ),
    p_value = c(
      original_test$p_value,
      global_test$p_value,
      pairwise_test$p_value,
      markov_naive_test$p_value,
      markov_corrected$omnibus$p_value
    ),
    reject = p_value < alpha,
    exact_decoding = c(
      TRUE,
      global_decode$exact_recovery,
      pairwise_decode$exact_recovery,
      FALSE,
      FALSE
    ),
    transition_condition_number = c(
      NA_real_, NA_real_, NA_real_,
      kappa(induced_released_mixing_matrix(K, rho, n_transitions)),
      markov_corrected$condition_number
    ),
    correction_status = c(
      NA_character_, NA_character_, NA_character_, NA_character_,
      markov_corrected$status
    )
  )
}

# ----------------------------------------------------------------
# 10.5 Monte Carlo driver for method and power comparisons
# ----------------------------------------------------------------

run_regrouping_method_comparison <- function(
    M = 1000,
    n = 500,
    G = 3,
    rho = rep(1 / G, G),
    p = 5,
    rho_x = 0.3,
    outcome_type = c("continuous", "binary"),
    delta = 0.5,
    sigma = 1,
    baseline_event = 0.30,
    K = make_symmetric_markov_kernel(G, 0.80),
    n_transitions = 1L,
    alpha = 0.05,
    seed = 20260723,
    progress = TRUE) {

  outcome_type <- match.arg(outcome_type)
  validate_probabilities(rho, G)
  validate_transition_kernel(K)

  if (!all(dim(K) == c(G, G))) {
    stop("K must be a G by G transition matrix.")
  }

  result_list <- vector("list", M)

  for (m in seq_len(M)) {
    result_list[[m]] <- compare_regrouping_methods_one_replication(
      replication = m,
      n = n,
      G = G,
      rho = rho,
      p = p,
      rho_x = rho_x,
      outcome_type = outcome_type,
      delta = delta,
      sigma = sigma,
      baseline_event = baseline_event,
      K = K,
      n_transitions = n_transitions,
      alpha = alpha,
      seed = seed
    )

    if (progress && (m %% max(1L, floor(M / 10L)) == 0L || m == M)) {
      message("  regrouping-method comparison completed ", m, " / ", M)
    }
  }

  raw <- bind_rows(result_list)

  summary <- raw %>%
    group_by(method) %>%
    summarise(
      estimable_replications = sum(is.finite(p_value)),
      rejection_probability = mean(reject, na.rm = TRUE),
      mean_statistic = mean(statistic, na.rm = TRUE),
      mean_p_value = mean(p_value, na.rm = TRUE),
      exact_decoding_probability = mean(exact_decoding),
      mean_transition_condition_number = if (
        all(is.na(transition_condition_number))
      ) NA_real_ else mean(transition_condition_number, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      estimand = ifelse(delta == 0, "Type-I error", "Power"),
      M = M,
      n = n,
      G = G,
      outcome_type = outcome_type,
      delta = delta,
      n_transitions = n_transitions,
      alpha = alpha,
      .before = 1
    )

  invariance <- raw %>%
    select(replication, method, statistic, p_value) %>%
    pivot_wider(
      names_from = method,
      values_from = c(statistic, p_value)
    ) %>%
    summarise(
      max_global_statistic_difference = max(abs(
        .data[["statistic_Global permutation"]] -
          .data[["statistic_Original labels"]]
      ), na.rm = TRUE),
      max_global_p_value_difference = max(abs(
        .data[["p_value_Global permutation"]] -
          .data[["p_value_Original labels"]]
      ), na.rm = TRUE),
      max_pairwise_statistic_difference = max(abs(
        .data[["statistic_Pairwise global swapping"]] -
          .data[["statistic_Original labels"]]
      ), na.rm = TRUE),
      max_pairwise_p_value_difference = max(abs(
        .data[["p_value_Pairwise global swapping"]] -
          .data[["p_value_Original labels"]]
      ), na.rm = TRUE)
    )

  list(
    raw = raw,
    summary = summary,
    deterministic_invariance = invariance,
    convergence = markov_convergence_diagnostics(K),
    effective_transition_kernel = matrix_power(K, n_transitions),
    induced_mixing_matrix = induced_released_mixing_matrix(
      K,
      rho,
      n_transitions
    )
  )
}

# ----------------------------------------------------------------
# 10.6 Examples
# ----------------------------------------------------------------

# Pairwise swapping example for G = 3:
# swaps_example <- matrix(c(1, 2,
#                           1, 3),
#                         ncol = 2, byrow = TRUE)
# pairwise_perm_example <- permutation_from_pairwise_swaps(
#   G = 3,
#   swaps = swaps_example
# )
# print(pairwise_perm_example)

# Markov-kernel example:
# K_example <- make_symmetric_markov_kernel(
#   G = 3,
#   stay_probability = 0.80
# )
# print(markov_convergence_diagnostics(K_example))

# Type-I error comparison:
# method_null <- run_regrouping_method_comparison(
#   M = 1000,
#   n = 500,
#   G = 3,
#   rho = c(1/3, 1/3, 1/3),
#   outcome_type = "continuous",
#   delta = 0,
#   K = K_example,
#   n_transitions = 1,
#   seed = 20260723
# )
# print(method_null$summary)
# print(method_null$deterministic_invariance)

# Power comparison:
# method_alternative <- run_regrouping_method_comparison(
#   M = 1000,
#   n = 500,
#   G = 3,
#   rho = c(1/3, 1/3, 1/3),
#   outcome_type = "continuous",
#   delta = 0.5,
#   K = K_example,
#   n_transitions = 1,
#   seed = 20260723
# )
# print(method_alternative$summary)

# Save the added-method outputs:
# write.csv(method_alternative$raw,
#           "regrouping_method_comparison_raw.csv", row.names = FALSE)
# write.csv(method_alternative$summary,
#           "regrouping_method_comparison_summary.csv", row.names = FALSE)
# write.csv(method_alternative$deterministic_invariance,
#           "deterministic_permutation_invariance.csv", row.names = FALSE)
# write.csv(method_alternative$convergence,
#           "markov_convergence_diagnostics.csv", row.names = FALSE)

# ================================================================
# End of additional regrouping methods
# ================================================================


# ================================================================
# 11. Execute regrouping-method comparisons and save outputs
# ================================================================
# This section is active. Running the full script creates the requested
# regrouping_method_comparison_summary.csv file.
#
# Increase method_comparison_M to 5000 for manuscript-sized runs.

method_output_dir <- "regrouping_method_outputs"
dir.create(method_output_dir, showWarnings = FALSE, recursive = TRUE)

method_comparison_M <- 200L
method_comparison_n <- 500L
method_comparison_G <- 3L
method_comparison_rho <- rep(1 / method_comparison_G, method_comparison_G)
method_comparison_p <- 5L
method_comparison_rho_x <- 0.30
method_comparison_outcome <- "continuous"
method_comparison_sigma <- 1
method_comparison_baseline_event <- 0.30
method_comparison_alpha <- 0.05
method_comparison_transitions <- 1L
method_comparison_stay_probability <- 0.80
method_comparison_seed <- 20260723L

method_comparison_K <- make_symmetric_markov_kernel(
  G = method_comparison_G,
  stay_probability = method_comparison_stay_probability
)

cat("\nRunning regrouping-method comparison under the null...\n")
method_null <- run_regrouping_method_comparison(
  M = method_comparison_M,
  n = method_comparison_n,
  G = method_comparison_G,
  rho = method_comparison_rho,
  p = method_comparison_p,
  rho_x = method_comparison_rho_x,
  outcome_type = method_comparison_outcome,
  delta = 0,
  sigma = method_comparison_sigma,
  baseline_event = method_comparison_baseline_event,
  K = method_comparison_K,
  n_transitions = method_comparison_transitions,
  alpha = method_comparison_alpha,
  seed = method_comparison_seed,
  progress = TRUE
)

cat("\nRunning regrouping-method comparison under the alternative...\n")
method_alternative <- run_regrouping_method_comparison(
  M = method_comparison_M,
  n = method_comparison_n,
  G = method_comparison_G,
  rho = method_comparison_rho,
  p = method_comparison_p,
  rho_x = method_comparison_rho_x,
  outcome_type = method_comparison_outcome,
  delta = 0.50,
  sigma = method_comparison_sigma,
  baseline_event = method_comparison_baseline_event,
  K = method_comparison_K,
  n_transitions = method_comparison_transitions,
  alpha = method_comparison_alpha,
  seed = method_comparison_seed + 100000L,
  progress = TRUE
)

# Long-format summary containing both type-I error and power.
regrouping_method_comparison_summary <- bind_rows(
  method_null$summary %>% mutate(simulation_hypothesis = "Null", .before = 1),
  method_alternative$summary %>% mutate(simulation_hypothesis = "Alternative", .before = 1)
) %>%
  arrange(method, simulation_hypothesis)

# Wide-format table with type-I error and power side by side.
regrouping_method_type1_power_summary <- regrouping_method_comparison_summary %>%
  select(
    method,
    simulation_hypothesis,
    rejection_probability,
    estimable_replications,
    mean_statistic,
    mean_p_value,
    exact_decoding_probability,
    mean_transition_condition_number
  ) %>%
  pivot_wider(
    names_from = simulation_hypothesis,
    values_from = c(
      rejection_probability,
      estimable_replications,
      mean_statistic,
      mean_p_value,
      exact_decoding_probability,
      mean_transition_condition_number
    ),
    names_glue = "{.value}_{simulation_hypothesis}"
  ) %>%
  rename(
    type_I_error = rejection_probability_Null,
    power = rejection_probability_Alternative
  )

# Raw replication-level results.
regrouping_method_comparison_raw <- bind_rows(
  method_null$raw %>% mutate(simulation_hypothesis = "Null", .before = 1),
  method_alternative$raw %>% mutate(simulation_hypothesis = "Alternative", .before = 1)
)

# Deterministic permutation invariance checks.
deterministic_permutation_invariance <- bind_rows(
  method_null$deterministic_invariance %>%
    mutate(simulation_hypothesis = "Null", .before = 1),
  method_alternative$deterministic_invariance %>%
    mutate(simulation_hypothesis = "Alternative", .before = 1)
)

# Transition-matrix diagnostics and matrices.
markov_convergence_diagnostics_output <- method_alternative$convergence

effective_transition_kernel_output <- as.data.frame(
  method_alternative$effective_transition_kernel,
  check.names = FALSE
)
effective_transition_kernel_output <- cbind(
  true_group = paste0("A", seq_len(nrow(effective_transition_kernel_output))),
  effective_transition_kernel_output
)

induced_mixing_matrix_output <- as.data.frame(
  method_alternative$induced_mixing_matrix,
  check.names = FALSE
)
induced_mixing_matrix_output <- cbind(
  released_group = rownames(method_alternative$induced_mixing_matrix),
  induced_mixing_matrix_output
)
rownames(induced_mixing_matrix_output) <- NULL

transition_kernel_output <- as.data.frame(method_comparison_K, check.names = FALSE)
colnames(transition_kernel_output) <- paste0("R", seq_len(method_comparison_G))
transition_kernel_output <- cbind(
  true_group = paste0("A", seq_len(method_comparison_G)),
  transition_kernel_output
)

# Save every requested output.
write.csv(
  regrouping_method_comparison_summary,
  file = file.path(
    method_output_dir,
    "regrouping_method_comparison_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  regrouping_method_type1_power_summary,
  file = file.path(
    method_output_dir,
    "regrouping_method_type1_power_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  regrouping_method_comparison_raw,
  file = file.path(
    method_output_dir,
    "regrouping_method_comparison_raw.csv"
  ),
  row.names = FALSE
)

write.csv(
  deterministic_permutation_invariance,
  file = file.path(
    method_output_dir,
    "deterministic_permutation_invariance.csv"
  ),
  row.names = FALSE
)

write.csv(
  markov_convergence_diagnostics_output,
  file = file.path(
    method_output_dir,
    "markov_convergence_diagnostics.csv"
  ),
  row.names = FALSE
)

write.csv(
  transition_kernel_output,
  file = file.path(
    method_output_dir,
    "markov_transition_kernel.csv"
  ),
  row.names = FALSE
)

write.csv(
  effective_transition_kernel_output,
  file = file.path(
    method_output_dir,
    "effective_transition_kernel.csv"
  ),
  row.names = FALSE
)

write.csv(
  induced_mixing_matrix_output,
  file = file.path(
    method_output_dir,
    "induced_released_mixing_matrix.csv"
  ),
  row.names = FALSE
)

cat("\n============================================================\n")
cat("Regrouping-method comparison completed.\n")
cat("Requested summary file:\n")
cat(normalizePath(file.path(
  method_output_dir,
  "regrouping_method_comparison_summary.csv"
), mustWork = FALSE), "\n")
cat("============================================================\n")

print(regrouping_method_comparison_summary)
print(regrouping_method_type1_power_summary)

# ================================================================
# End of complete revised simulation program
# ================================================================
