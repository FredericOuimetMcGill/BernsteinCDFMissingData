
suppressPackageStartupMessages({
  library(parallel)
  library(ggplot2)
})

options(stringsAsFactors = FALSE, scipen = 999)

# ============================================================
# Paths
# ============================================================

PATH <- "C://Users//oufr6443//Downloads//test3"

if (!dir.exists(PATH)) {
  stop(sprintf("The path does not exist: %s", PATH))
}

setwd(PATH)

TABLE_DIR <- file.path(getwd(), "simulation_tables")
FIG_DIR   <- file.path(getwd(), "simulation_figures")

dir_create_safe <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

dir_create_safe(TABLE_DIR)
dir_create_safe(FIG_DIR)

# ============================================================
# Hyperparameters (edit here)
# ============================================================

SEED <- 12345

# Parallelization
N_CORES_AVAILABLE <- parallel::detectCores(logical = TRUE)
if (is.na(N_CORES_AVAILABLE) || N_CORES_AVAILABLE < 2L) {
  N_CORES_AVAILABLE <- 2L
}
N_CORES <- max(1L, N_CORES_AVAILABLE - 1L)

# Data-generating mechanism
ALPHA_Y <- 0.9
BETA_Y  <- 0.9
RHO     <- 0.6

# Benchmark logistic MAR mechanism
PI0_BENCHMARK <- 0.6
PI1_BENCHMARK <- 0.9
BETA0_BENCHMARK <- qlogis(PI0_BENCHMARK)
BETA1_BENCHMARK <- qlogis(PI1_BENCHMARK) - qlogis(PI0_BENCHMARK)

# Monte Carlo designs
SAMPLE_SIZE_GRID    <- c(25, 50, 100, 200, 400, 800, 1600, 3200, 6400)
N_FIXED_MISSING     <- 400
MISSING_RATE_GRID   <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40)
MC_REPS             <- 1000

# Numerical integration
GRID_POINTS <- 2^7

# Bernstein LSCV grid
M_MIN    <- 1L
M_CAP    <- 300L
M_FACTOR <- 5

# KDE LSCV grid
H_C_MIN       <- 0.05
H_C_MAX       <- 1.00
H_GRID_LENGTH <- 20L

# Numerical integration hyperparameters for ISE / BISE
BISE_HALF_WIDTH_DENOM <- 1

RISK_NUMINT_METHOD <- "adaptIntegrate"
# Allowed values:
#   "integrate"       : stats::integrate
#   "adaptIntegrate"  : cubature::adaptIntegrate
#   "adaptiveSimpson" : custom adaptive Simpson rule

RISK_NUMINT_REL_TOL       <- 1e-2
RISK_NUMINT_ABS_TOL       <- 0
RISK_NUMINT_SUBDIVISIONS  <- 200L   # used by stats::integrate only
RISK_NUMINT_STOP_ON_ERROR <- TRUE   # used by stats::integrate only
RISK_ADAPTINT_MAX_EVAL    <- 5000L  # used by cubature::adaptIntegrate only
RISK_ASIMPSON_MAX_DEPTH   <- 12L    # used by adaptiveSimpson only

# Numerical safeguards
PI_HAT_FLOOR <- 1e-8
PROB_CLIP    <- 1e-12

# Output formatting
TABLE_STAT_MULTIPLIER <- 1e8
TABLE_DIGITS          <- 0L

# LaTeX table interline spacing (\arraystretch; 1 = default)
TABLE_USE_ARRAYSTRETCH <- TRUE
TABLE_ARRAYSTRETCH     <- 1.15

# Boxplot aesthetics
BOXPLOT_WIDTH         <- 9
BOXPLOT_HEIGHT        <- 5.75
BOXPLOT_BASE_SIZE     <- 14
BOXPLOT_Y_PAD_FRACTION <- 0.04
ESTIMATOR_COLORS <- c(
  "Unsmoothed" = "#1f77b4",
  "I-IPW KDE"  = "#2ca02c",
  "Bernstein"  = "#d62728"
)

# Other toggles
VERBOSE <- TRUE

if (!is.numeric(TABLE_STAT_MULTIPLIER) ||
    length(TABLE_STAT_MULTIPLIER) != 1L ||
    !is.finite(TABLE_STAT_MULTIPLIER) ||
    TABLE_STAT_MULTIPLIER <= 0) {
  stop("TABLE_STAT_MULTIPLIER must be a single positive finite number.")
}

TABLE_DIGITS <- as.integer(TABLE_DIGITS)
if (!is.finite(TABLE_DIGITS) || length(TABLE_DIGITS) != 1L || TABLE_DIGITS < 0L) {
  stop("TABLE_DIGITS must be a single nonnegative integer.")
}

if (!is.logical(TABLE_USE_ARRAYSTRETCH) ||
    length(TABLE_USE_ARRAYSTRETCH) != 1L ||
    is.na(TABLE_USE_ARRAYSTRETCH)) {
  stop("TABLE_USE_ARRAYSTRETCH must be TRUE or FALSE.")
}

if (!is.numeric(TABLE_ARRAYSTRETCH) ||
    length(TABLE_ARRAYSTRETCH) != 1L ||
    !is.finite(TABLE_ARRAYSTRETCH) ||
    TABLE_ARRAYSTRETCH <= 0) {
  stop("TABLE_ARRAYSTRETCH must be a single positive finite number.")
}

# ============================================================
# Fixed labels / naming conventions
# ============================================================

ESTIMATOR_KEYS <- c("unsmoothed", "kde", "bernstein")
ESTIMATOR_LABELS <- c(
  unsmoothed = "Unsmoothed",
  kde        = "I-IPW KDE",
  bernstein  = "Bernstein"
)
ESTIMATOR_LEVELS <- unname(ESTIMATOR_LABELS[ESTIMATOR_KEYS])

REGIME_LEVELS  <- c("pseudo", "feasible")
MEASURE_LEVELS <- c("ISE", "BISE")

# ============================================================
# Timers
# ============================================================

timer_report <- function(label, start_time, end_time = Sys.time()) {
  elapsed_sec <- as.numeric(difftime(end_time, start_time, units = "secs"))
  elapsed_min <- elapsed_sec / 60
  message(sprintf("[TIMER] %s: %.2f seconds (%.2f minutes)", label, elapsed_sec, elapsed_min))
  invisible(list(seconds = elapsed_sec, minutes = elapsed_min))
}

# ============================================================
# True CDF and sampling
# ============================================================

F_true <- function(y) {
  pbeta(y, shape1 = ALPHA_Y, shape2 = BETA_Y)
}

rY <- function(n) {
  rbeta(n, shape1 = ALPHA_Y, shape2 = BETA_Y)
}

# ============================================================
# Latent-Gaussian construction for X
# ============================================================

generate_auxiliary_from_Y <- function(Y, rho = RHO, prob_clip = PROB_CLIP) {
  U <- pbeta(Y, shape1 = ALPHA_Y, shape2 = BETA_Y)
  U <- pmin(pmax(U, prob_clip), 1 - prob_clip)
  T <- qnorm(U)
  
  Z <- rnorm(length(Y))
  X_star <- rho * T + sqrt(1 - rho^2) * Z
  X <- as.integer(X_star > 0)
  
  list(U = U, T = T, Z = Z, X_star = X_star, X = X)
}

# ============================================================
# Logistic propensity utilities
# ============================================================

propensity_fun <- function(x, beta0, beta1) {
  plogis(beta0 + beta1 * x)
}

expected_observation_prob <- function(beta0, beta1) {
  0.5 * (plogis(beta0) + plogis(beta0 + beta1))
}

solve_beta0_for_missing_rate <- function(target_missing_rate,
                                         beta1 = BETA1_BENCHMARK,
                                         interval = c(-30, 30)) {
  target_observation_rate <- 1 - target_missing_rate
  
  objective <- function(beta0) {
    expected_observation_prob(beta0, beta1) - target_observation_rate
  }
  
  uniroot(objective, interval = interval, tol = 1e-12)$root
}

# ============================================================
# Dataset simulation
# ============================================================

simulate_dataset <- function(n, beta0, beta1, rho = RHO) {
  Y <- rY(n)
  
  aux <- generate_auxiliary_from_Y(Y = Y, rho = rho, prob_clip = PROB_CLIP)
  X  <- aux$X
  
  pi    <- propensity_fun(X, beta0 = beta0, beta1 = beta1)
  delta <- rbinom(n, size = 1, prob = pi)
  
  list(
    Y      = Y,
    X      = X,
    pi     = pi,
    delta  = delta,
    T      = aux$T,
    X_star = aux$X_star
  )
}

# ============================================================
# Propensity estimation from discrete X
# ============================================================

estimate_pi_hat <- function(X, delta, pi_floor = PI_HAT_FLOOR) {
  X <- as.integer(X)
  out <- numeric(length(X))
  
  for (val in sort(unique(X))) {
    idx <- which(X == val)
    den <- length(idx)
    pih <- if (den > 0L) sum(delta[idx]) / den else pi_floor
    pih <- max(pih, pi_floor)
    out[idx] <- pih
  }
  
  out
}

# ============================================================
# IPW empirical CDF
# ============================================================

ipw_ecdf_at <- function(Y, w, t_points, n_total) {
  obs <- which(w > 0)
  
  if (length(obs) == 0L) {
    return(rep(0, length(t_points)))
  }
  
  Yobs <- Y[obs]
  wobs <- w[obs]
  
  ord <- order(Yobs)
  Yso <- Yobs[ord]
  wso <- wobs[ord]
  
  cs <- cumsum(wso)
  
  idx <- findInterval(t_points, Yso, left.open = FALSE, rightmost.closed = TRUE)
  Fv  <- numeric(length(idx))
  
  pos <- idx > 0L
  Fv[pos] <- cs[idx[pos]]
  
  (1 / n_total) * Fv
}

# ============================================================
# Bernstein smoothing and Bernstein LSCV
# ============================================================

bernstein_smooth <- function(F_at_knots, y_points, m) {
  ks <- 0:m
  vapply(
    y_points,
    function(y) sum(F_at_knots * dbinom(ks, size = m, prob = y)),
    numeric(1)
  )
}

m_candidates_of_n <- function(n) {
  m_max <- min(floor(M_FACTOR * n^(2 / 3)), M_CAP, n)
  if (m_max < M_MIN) {
    m_max <- M_MIN
  }
  seq.int(M_MIN, m_max, by = 1L)
}

bernstein_term1_exact <- function(F_knots, m) {
  ks <- 0:m
  Ikl <- outer(
    ks, ks,
    function(k, l) choose(m, k) * choose(m, l) * beta(k + l + 1, 2 * m - k - l + 1)
  )
  sum(outer(F_knots, F_knots) * Ikl)
}

bernstein_tail_integrals <- function(y, m) {
  ks <- 0:m
  (1 - pbeta(y, ks + 1, m - ks + 1)) / (m + 1)
}

cv_score_bernstein <- function(m, Y, w, n_total) {
  if (n_total <= 1L) {
    return(Inf)
  }
  
  t_k     <- seq(0, 1, length.out = m + 1L)
  F_knots <- ipw_ecdf_at(Y, w, t_k, n_total = n_total)
  
  term1 <- bernstein_term1_exact(F_knots, m)
  
  obs <- which(w > 0)
  if (length(obs) == 0L) {
    return(term1)
  }
  
  sumW_leq_tk <- n_total * F_knots
  
  acc <- 0
  for (i in obs) {
    ind_i     <- as.numeric(t_k >= Y[i])
    F_minus_k <- (sumW_leq_tk - w[i] * ind_i) / (n_total - 1)
    tail_i    <- bernstein_tail_integrals(Y[i], m)
    acc <- acc + w[i] * sum(F_minus_k * tail_i)
  }
  
  term2 <- (2 / n_total) * acc
  
  term1 - term2
}

# ============================================================
# Integrated IPW Gaussian KDE CDF and KDE LSCV
# ============================================================

ipw_gaussian_kde_cdf <- function(Y, w, y_points, n_total, h) {
  obs <- which(w > 0 & is.finite(Y))
  
  if (length(obs) == 0L) {
    return(rep(0, length(y_points)))
  }
  
  Yo <- Y[obs]
  wo <- w[obs]
  
  vapply(
    y_points,
    function(y) {
      sum(wo * pnorm((y - Yo) / h)) / n_total
    },
    numeric(1)
  )
}

h_candidates_of_n <- function(n) {
  c_seq <- seq(H_C_MIN, H_C_MAX, length.out = H_GRID_LENGTH)
  c_seq * n^(-1 / 5)
}

cv_score_kde <- function(h, Y, w, n_total, y_grid) {
  Glen <- length(y_grid)
  
  if (!is.finite(h) || h <= 0 || n_total <= 1L) {
    return(Inf)
  }
  
  obs <- which(w > 0 & is.finite(Y))
  if (length(obs) <= 1L) {
    return(Inf)
  }
  
  Yo <- Y[obs]
  wo <- w[obs]
  
  S_vec <- numeric(Glen)
  for (i in seq_along(Yo)) {
    S_vec <- S_vec + wo[i] * pnorm((y_grid - Yo[i]) / h)
  }
  
  F_vec <- S_vec / n_total
  term1 <- mean(F_vec^2)
  
  term2_acc <- 0
  for (i in seq_along(Yo)) {
    add_i <- wo[i] * pnorm((y_grid - Yo[i]) / h)
    F_mi  <- (S_vec - add_i) / (n_total - 1L)
    
    sel  <- y_grid >= Yo[i]
    tail <- sum(F_mi[sel]) / Glen
    
    term2_acc <- term2_acc + wo[i] * tail
  }
  
  term2 <- (2 / n_total) * term2_acc
  
  term1 - term2
}

# ============================================================
# Risk criteria
# ============================================================

boundary_mask_of_n <- function(n, y_grid) {
  delta <- n^(-2 / 3) / BISE_HALF_WIDTH_DENOM
  as.numeric((y_grid <= delta) | (y_grid >= (1 - delta)))
}

normalize_numint_method <- function(method = RISK_NUMINT_METHOD) {
  key <- tolower(gsub("[^[:alnum:]]", "", method))
  
  if (key == "integrate") {
    return("integrate")
  }
  if (key %in% c("adaptintegrate", "adapt", "cubature")) {
    return("adaptIntegrate")
  }
  if (key %in% c("adaptivesimpson", "simpson")) {
    return("adaptiveSimpson")
  }
  
  stop("Unknown RISK_NUMINT_METHOD. Use 'integrate', 'adaptIntegrate', or 'adaptiveSimpson'.")
}

make_risk_integrand <- function(est_on_grid, F_true_on_grid) {
  if (length(est_on_grid) != length(F_true_on_grid)) {
    stop("est_on_grid and F_true_on_grid must have the same length.")
  }
  if (length(est_on_grid) < 2L) {
    stop("At least two grid points are required for numerical integration.")
  }
  if (any(!is.finite(est_on_grid)) || any(!is.finite(F_true_on_grid))) {
    stop("est_on_grid and F_true_on_grid must be finite.")
  }
  
  y_grid_local <- seq(0, 1, length.out = length(est_on_grid))
  
  est_fun <- approxfun(
    x = y_grid_local,
    y = est_on_grid,
    method = "linear",
    rule = 2,
    ties = "ordered"
  )
  
  true_fun <- approxfun(
    x = y_grid_local,
    y = F_true_on_grid,
    method = "linear",
    rule = 2,
    ties = "ordered"
  )
  
  function(y) {
    d <- est_fun(y) - true_fun(y)
    d * d
  }
}

adaptive_simpson_1d <- function(f, lower, upper,
                                rel_tol = RISK_NUMINT_REL_TOL,
                                abs_tol = RISK_NUMINT_ABS_TOL,
                                max_depth = RISK_ASIMPSON_MAX_DEPTH) {
  if (!is.finite(lower) || !is.finite(upper) || upper < lower) {
    stop("Invalid integration bounds.")
  }
  if (upper == lower) {
    return(0)
  }
  
  max_depth <- max(0L, as.integer(max_depth))
  
  simpson_panel <- function(a, b, fa, fm, fb) {
    (b - a) * (fa + 4 * fm + fb) / 6
  }
  
  recurse <- function(a, b, fa, fm, fb, S, tol_here, depth_left) {
    m  <- (a + b) / 2
    lm <- (a + m) / 2
    rm <- (m + b) / 2
    
    flm <- f(lm)
    frm <- f(rm)
    
    S_left  <- simpson_panel(a, m, fa, flm, fm)
    S_right <- simpson_panel(m, b, fm, frm, fb)
    S2 <- S_left + S_right
    
    err_est <- abs(S2 - S)
    
    if (depth_left <= 0L || err_est <= 15 * tol_here) {
      return(S2 + (S2 - S) / 15)
    }
    
    recurse(a, m, fa, flm, fm, S_left,  tol_here / 2, depth_left - 1L) +
      recurse(m, b, fm, frm, fb, S_right, tol_here / 2, depth_left - 1L)
  }
  
  a  <- lower
  b  <- upper
  m  <- (a + b) / 2
  fa <- f(a)
  fm <- f(m)
  fb <- f(b)
  
  S0   <- simpson_panel(a, b, fa, fm, fb)
  tol0 <- max(abs_tol, rel_tol * max(1, abs(S0)), .Machine$double.eps)
  
  recurse(a, b, fa, fm, fb, S0, tol0, max_depth)
}

numerical_integral_1d <- function(f, lower, upper) {
  method <- normalize_numint_method(RISK_NUMINT_METHOD)
  
  if (!is.finite(lower) || !is.finite(upper) || upper < lower) {
    stop("Invalid integration bounds.")
  }
  if (upper == lower) {
    return(0)
  }
  
  if (method == "integrate") {
    out <- stats::integrate(
      f = f,
      lower = lower,
      upper = upper,
      subdivisions = as.integer(RISK_NUMINT_SUBDIVISIONS),
      rel.tol = RISK_NUMINT_REL_TOL,
      abs.tol = RISK_NUMINT_ABS_TOL,
      stop.on.error = RISK_NUMINT_STOP_ON_ERROR
    )
    return(as.numeric(out$value))
  }
  
  if (method == "adaptIntegrate") {
    if (!requireNamespace("cubature", quietly = TRUE)) {
      stop("RISK_NUMINT_METHOD = 'adaptIntegrate' requires the 'cubature' package.")
    }
    
    f_cubature <- function(x) {
      if (is.matrix(x)) {
        return(vapply(seq_len(ncol(x)), function(j) f(x[1L, j]), numeric(1)))
      }
      f(as.numeric(x)[1L])
    }
    
    out <- cubature::adaptIntegrate(
      f = f_cubature,
      lowerLimit = lower,
      upperLimit = upper,
      tol = RISK_NUMINT_REL_TOL,
      absError = RISK_NUMINT_ABS_TOL,
      maxEval = as.integer(RISK_ADAPTINT_MAX_EVAL)
    )
    return(as.numeric(out$integral))
  }
  
  if (method == "adaptiveSimpson") {
    return(
      adaptive_simpson_1d(
        f = f,
        lower = lower,
        upper = upper,
        rel_tol = RISK_NUMINT_REL_TOL,
        abs_tol = RISK_NUMINT_ABS_TOL,
        max_depth = RISK_ASIMPSON_MAX_DEPTH
      )
    )
  }
  
  stop("Unknown numerical integration method.")
}

compute_ISE <- function(est_on_grid, F_true_on_grid) {
  err2_fun <- make_risk_integrand(est_on_grid, F_true_on_grid)
  numerical_integral_1d(err2_fun, lower = 0, upper = 1)
}

compute_BISE <- function(est_on_grid, F_true_on_grid, boundary_mask, n) {
  delta <- n^(-2 / 3) / BISE_HALF_WIDTH_DENOM
  
  if (!is.finite(delta) || delta <= 0) {
    stop("The boundary width delta must be positive and finite.")
  }
  
  err2_fun <- make_risk_integrand(est_on_grid, F_true_on_grid)
  
  boundary_integral <- if (delta >= 0.5) {
    numerical_integral_1d(err2_fun, lower = 0, upper = 1)
  } else {
    numerical_integral_1d(err2_fun, lower = 0, upper = delta) +
      numerical_integral_1d(err2_fun, lower = 1 - delta, upper = 1)
  }
  
  boundary_integral / (2 * delta)
}

# ============================================================
# One Monte Carlo replication
# ============================================================

one_replication <- function(n, beta0, beta1, y_grid, F_true_grid,
                            m_grid, h_grid, boundary_mask) {
  dat   <- simulate_dataset(n = n, beta0 = beta0, beta1 = beta1, rho = RHO)
  Y     <- dat$Y
  X     <- dat$X
  pi    <- dat$pi
  delta <- dat$delta
  
  # Pseudo weights
  w_pseudo <- delta / pi
  
  # Feasible weights
  pi_hat <- estimate_pi_hat(X, delta, pi_floor = PI_HAT_FLOOR)
  w_feasible <- delta / pi_hat
  
  # Pseudo: unsmoothed
  F_uns_pseudo <- ipw_ecdf_at(Y, w_pseudo, y_grid, n_total = n)
  
  # Pseudo: Bernstein
  cv_m_pseudo <- vapply(
    m_grid,
    cv_score_bernstein,
    numeric(1),
    Y = Y,
    w = w_pseudo,
    n_total = n
  )
  m_hat_pseudo <- m_grid[which.min(cv_m_pseudo)]
  
  t_k_pseudo <- seq(0, 1, length.out = m_hat_pseudo + 1L)
  F_knots_pseudo <- ipw_ecdf_at(Y, w_pseudo, t_k_pseudo, n_total = n)
  F_bern_pseudo <- bernstein_smooth(F_knots_pseudo, y_grid, m_hat_pseudo)
  
  # Pseudo: KDE
  cv_h_pseudo <- vapply(
    h_grid,
    cv_score_kde,
    numeric(1),
    Y = Y,
    w = w_pseudo,
    n_total = n,
    y_grid = y_grid
  )
  h_hat_pseudo <- h_grid[which.min(cv_h_pseudo)]
  F_kde_pseudo <- ipw_gaussian_kde_cdf(Y, w_pseudo, y_grid, n_total = n, h = h_hat_pseudo)
  
  # Feasible: unsmoothed
  F_uns_feasible <- ipw_ecdf_at(Y, w_feasible, y_grid, n_total = n)
  
  # Feasible: Bernstein
  cv_m_feasible <- vapply(
    m_grid,
    cv_score_bernstein,
    numeric(1),
    Y = Y,
    w = w_feasible,
    n_total = n
  )
  m_hat_feasible <- m_grid[which.min(cv_m_feasible)]
  
  t_k_feasible <- seq(0, 1, length.out = m_hat_feasible + 1L)
  F_knots_feasible <- ipw_ecdf_at(Y, w_feasible, t_k_feasible, n_total = n)
  F_bern_feasible <- bernstein_smooth(F_knots_feasible, y_grid, m_hat_feasible)
  
  # Feasible: KDE
  cv_h_feasible <- vapply(
    h_grid,
    cv_score_kde,
    numeric(1),
    Y = Y,
    w = w_feasible,
    n_total = n,
    y_grid = y_grid
  )
  h_hat_feasible <- h_grid[which.min(cv_h_feasible)]
  F_kde_feasible <- ipw_gaussian_kde_cdf(Y, w_feasible, y_grid, n_total = n, h = h_hat_feasible)
  
  c(
    pseudo_unsmoothed_ISE = compute_ISE(F_uns_pseudo, F_true_grid),
    pseudo_kde_ISE        = compute_ISE(F_kde_pseudo, F_true_grid),
    pseudo_bernstein_ISE  = compute_ISE(F_bern_pseudo, F_true_grid),
    
    pseudo_unsmoothed_BISE = compute_BISE(F_uns_pseudo, F_true_grid, boundary_mask, n),
    pseudo_kde_BISE        = compute_BISE(F_kde_pseudo, F_true_grid, boundary_mask, n),
    pseudo_bernstein_BISE  = compute_BISE(F_bern_pseudo, F_true_grid, boundary_mask, n),
    
    feasible_unsmoothed_ISE = compute_ISE(F_uns_feasible, F_true_grid),
    feasible_kde_ISE        = compute_ISE(F_kde_feasible, F_true_grid),
    feasible_bernstein_ISE  = compute_ISE(F_bern_feasible, F_true_grid),
    
    feasible_unsmoothed_BISE = compute_BISE(F_uns_feasible, F_true_grid, boundary_mask, n),
    feasible_kde_BISE        = compute_BISE(F_kde_feasible, F_true_grid, boundary_mask, n),
    feasible_bernstein_BISE  = compute_BISE(F_bern_feasible, F_true_grid, boundary_mask, n)
  )
}

# ============================================================
# Convert a results matrix to long format
# ============================================================

results_matrix_to_long <- function(res_mat, experiment_name, setting_value, setting_label, n_used) {
  out_list <- list()
  idx <- 1L
  
  for (regime in REGIME_LEVELS) {
    for (measure in MEASURE_LEVELS) {
      for (est_key in ESTIMATOR_KEYS) {
        col_name <- paste(regime, est_key, measure, sep = "_")
        out_list[[idx]] <- data.frame(
          experiment    = experiment_name,
          setting_value = rep(setting_value, nrow(res_mat)),
          setting_label = rep(setting_label, nrow(res_mat)),
          n             = rep(n_used, nrow(res_mat)),
          regime        = rep(regime, nrow(res_mat)),
          measure       = rep(measure, nrow(res_mat)),
          estimator     = rep(ESTIMATOR_LABELS[[est_key]], nrow(res_mat)),
          value         = res_mat[, col_name],
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
      }
    }
  }
  
  do.call(rbind, out_list)
}

# ============================================================
# Summaries for tables
# ============================================================

scale_summary_statistics <- function(summary_df, id_col_name,
                                     scale_factor = TABLE_STAT_MULTIPLIER) {
  stat_cols <- setdiff(names(summary_df), id_col_name)
  
  if (length(stat_cols) == 0L || abs(scale_factor - 1) <= .Machine$double.eps^0.5) {
    return(summary_df)
  }
  
  summary_df[stat_cols] <- lapply(summary_df[stat_cols], function(x) scale_factor * x)
  summary_df
}

make_summary_table <- function(long_df, experiment_name, regime, measure,
                               scale_factor = TABLE_STAT_MULTIPLIER) {
  sub_df <- long_df[
    long_df$experiment == experiment_name &
      long_df$regime == regime &
      long_df$measure == measure,
    ,
    drop = FALSE
  ]
  
  settings_df <- unique(sub_df[, c("setting_value", "setting_label")])
  settings_df <- settings_df[order(settings_df$setting_value), , drop = FALSE]
  
  id_col_name <- if (experiment_name == "sample_size") "n" else "missing_rate_percent"
  
  out <- data.frame(setting_value = settings_df$setting_value)
  names(out)[1] <- id_col_name
  
  for (est_key in ESTIMATOR_KEYS) {
    est_label <- ESTIMATOR_LABELS[[est_key]]
    
    mean_vals <- numeric(nrow(settings_df))
    sd_vals   <- numeric(nrow(settings_df))
    
    for (i in seq_len(nrow(settings_df))) {
      sv <- settings_df$setting_value[i]
      vals <- sub_df$value[sub_df$setting_value == sv & sub_df$estimator == est_label]
      
      mean_vals[i] <- mean(vals)
      sd_vals[i]   <- sd(vals)
    }
    
    out[[paste0("mean_", est_key)]] <- mean_vals
    out[[paste0("sd_", est_key)]]   <- sd_vals
  }
  
  ordered_cols <- c(
    id_col_name,
    paste0("mean_", ESTIMATOR_KEYS),
    paste0("sd_", ESTIMATOR_KEYS)
  )
  
  out <- out[, ordered_cols, drop = FALSE]
  scale_summary_statistics(out, id_col_name = id_col_name, scale_factor = scale_factor)
}

# ============================================================
# Table writers
# ============================================================

fmt_num <- function(x, digits = TABLE_DIGITS) {
  formatC(x, format = "f", digits = digits)
}

latex_scale_suffix <- function(scale_factor = TABLE_STAT_MULTIPLIER, tol = 1e-12) {
  if (abs(scale_factor - 1) <= tol) {
    return("")
  }
  
  exponent <- log10(scale_factor)
  if (is.finite(exponent) && abs(exponent - round(exponent)) <= tol) {
    return(paste0(" ($\\times 10^{", as.integer(round(exponent)), "}$)"))
  }
  
  factor_label <- formatC(scale_factor, format = "fg", digits = 10)
  paste0(" ($\\times ", factor_label, "$)")
}

stat_block_header <- function(stat_name, measure,
                              scale_factor = TABLE_STAT_MULTIPLIER) {
  measure_tag <- if (measure == "ISE") "ISE" else "BISE"
  paste0(stat_name, " ", measure_tag, latex_scale_suffix(scale_factor = scale_factor))
}

make_table_caption <- function(experiment_name, regime, measure) {
  measure_phrase <- if (measure == "ISE") "ISE statistics" else "Boundary ISE statistics"
  regime_phrase <- if (regime == "pseudo") "pseudo estimators" else "feasible estimators"
  setting_phrase <- if (experiment_name == "sample_size") {
    "as a function of the sample size"
  } else {
    "as a function of the missing rate"
  }
  
  paste0(measure_phrase, " for ", regime_phrase, " ", setting_phrase, ".")
}

write_summary_tex <- function(summary_df, experiment_name, regime, measure, file_tex,
                              digits = TABLE_DIGITS, tol = 1e-12,
                              scale_factor = TABLE_STAT_MULTIPLIER) {
  first_col_display <- if (experiment_name == "sample_size") "$n$" else "Missing rate (\\%)"
  caption_text <- make_table_caption(experiment_name = experiment_name,
                                     regime = regime,
                                     measure = measure)
  label_text <- paste0("tab:", regime, "_", measure, "_", experiment_name)
  
  mean_header <- stat_block_header("Mean", measure, scale_factor = scale_factor)
  sd_header   <- stat_block_header("Standard deviation", measure, scale_factor = scale_factor)
  
  con <- file(file_tex, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  
  writeLines("\\begin{table}[H]", con)
  writeLines("\\centering", con)
  writeLines("\\small", con)
  writeLines("\\begingroup", con)
  
  if (isTRUE(TABLE_USE_ARRAYSTRETCH)) {
    arraystretch_text <- formatC(TABLE_ARRAYSTRETCH, format = "fg", digits = 10)
    writeLines(paste0("\\renewcommand{\\arraystretch}{", arraystretch_text, "}"), con)
  }
  
  writeLines("\\begin{tabular}{|r|ccc|ccc|}", con)
  writeLines("\\hline", con)
  
  hdr1 <- c(
    first_col_display,
    paste0("\\multicolumn{3}{|c|}{", mean_header, "}"),
    paste0("\\multicolumn{3}{|c|}{", sd_header, "}")
  )
  writeLines(paste0(paste(hdr1, collapse = " & "), " \\\\"), con)
  writeLines("\\hline", con)
  
  hdr2 <- c(" ", ESTIMATOR_LEVELS, ESTIMATOR_LEVELS)
  writeLines(paste0(paste(hdr2, collapse = " & "), " \\\\"), con)
  writeLines("\\hline", con)
  
  mean_cols <- paste0("mean_", ESTIMATOR_KEYS)
  sd_cols   <- paste0("sd_", ESTIMATOR_KEYS)
  
  for (i in seq_len(nrow(summary_df))) {
    row_chars <- character(7L)
    
    first_value <- summary_df[i, 1]
    if (experiment_name == "sample_size") {
      row_chars[1] <- as.character(as.integer(first_value))
    } else {
      row_chars[1] <- as.character(as.integer(round(first_value)))
    }
    
    mean_vals <- as.numeric(summary_df[i, mean_cols])
    sd_vals   <- as.numeric(summary_df[i, sd_cols])
    
    mean_chars <- vapply(mean_vals, fmt_num, "", digits = digits)
    sd_chars   <- vapply(sd_vals, fmt_num, "", digits = digits)
    
    mean_min <- min(mean_vals)
    sd_min   <- min(sd_vals)
    
    mean_chars[abs(mean_vals - mean_min) <= tol] <-
      paste0("\\textbf{", mean_chars[abs(mean_vals - mean_min) <= tol], "}")
    sd_chars[abs(sd_vals - sd_min) <= tol] <-
      paste0("\\textbf{", sd_chars[abs(sd_vals - sd_min) <= tol], "}")
    
    row_chars[2:4] <- mean_chars
    row_chars[5:7] <- sd_chars
    
    writeLines(paste0(paste(row_chars, collapse = " & "), " \\\\"), con)
  }
  
  writeLines("\\hline", con)
  writeLines("\\end{tabular}", con)
  writeLines("\\endgroup", con)
  writeLines(paste0("\\caption{", caption_text, "}"), con)
  writeLines(paste0("\\label{", label_text, "}"), con)
  writeLines("\\end{table}", con)
}

save_summary_outputs <- function(summary_df, experiment_name, regime, measure) {
  file_stub <- paste(regime, measure, experiment_name, sep = "_")
  
  csv_file <- file.path(TABLE_DIR, paste0(file_stub, ".csv"))
  tex_file <- file.path(TABLE_DIR, paste0(file_stub, ".tex"))
  
  write.csv(summary_df, csv_file, row.names = FALSE)
  write_summary_tex(
    summary_df = summary_df,
    experiment_name = experiment_name,
    regime = regime,
    measure = measure,
    file_tex = tex_file,
    digits = TABLE_DIGITS,
    scale_factor = TABLE_STAT_MULTIPLIER
  )
  
  invisible(list(csv = csv_file, tex = tex_file))
}

# ============================================================
# Boxplots
# ============================================================

make_boxplot_caption <- function(experiment_name, regime, measure) {
  experiment_phrase <- if (experiment_name == "sample_size") {
    "sample-size experiment"
  } else {
    "missing-rate experiment"
  }
  
  measure_phrase <- if (measure == "ISE") "ISE" else "boundary ISE"
  regime_phrase <- if (regime == "pseudo") "pseudo estimators" else "feasible estimators"
  
  paste0("Boxplots of the ", measure_phrase, " for the ", regime_phrase,
         " in the ", experiment_phrase, ".")
}

compute_boxplot_ylim <- function(values, groups, pad_fraction = BOXPLOT_Y_PAD_FRACTION) {
  split_vals <- split(values, groups, drop = TRUE)
  
  whisker_bounds <- lapply(split_vals, function(v) {
    stats <- boxplot.stats(v)$stats
    c(lower = stats[1], upper = stats[5])
  })
  
  whisker_mat <- do.call(rbind, whisker_bounds)
  lower_whisker <- min(whisker_mat[, "lower"], na.rm = TRUE)
  upper_whisker <- max(whisker_mat[, "upper"], na.rm = TRUE)
  
  span <- upper_whisker - lower_whisker
  pad <- if (is.finite(span) && span > 0) {
    pad_fraction * span
  } else {
    max(pad_fraction * max(abs(c(lower_whisker, upper_whisker)), na.rm = TRUE), 1e-8)
  }
  
  c(max(0, lower_whisker - pad), upper_whisker + pad)
}

save_boxplot <- function(long_df, experiment_name, regime, measure) {
  sub_df <- long_df[
    long_df$experiment == experiment_name &
      long_df$regime == regime &
      long_df$measure == measure,
    ,
    drop = FALSE
  ]
  
  ord <- unique(sub_df[, c("setting_value", "setting_label")])
  ord <- ord[order(ord$setting_value), , drop = FALSE]
  
  sub_df$setting_factor <- factor(sub_df$setting_label, levels = ord$setting_label)
  sub_df$estimator <- factor(sub_df$estimator, levels = ESTIMATOR_LEVELS)
  
  x_lab <- if (experiment_name == "sample_size") "Sample size n" else "Target missing rate"
  y_lab <- if (measure == "ISE") "ISE" else "Boundary ISE"
  legend_title <- if (regime == "pseudo") "Pseudo estimator" else "Feasible estimator"
  
  whisker_groups <- interaction(sub_df$setting_factor, sub_df$estimator, drop = TRUE)
  y_limits <- compute_boxplot_ylim(sub_df$value, whisker_groups)
  
  p <- ggplot(sub_df, aes(x = setting_factor, y = value, fill = estimator)) +
    geom_boxplot(
      position = position_dodge2(width = 0.80, preserve = "single"),
      width = 0.70,
      outlier.shape = NA
    ) +
    scale_fill_manual(values = ESTIMATOR_COLORS, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    coord_cartesian(ylim = y_limits) +
    labs(
      x = x_lab,
      y = y_lab,
      fill = legend_title
    ) +
    theme_minimal(base_size = BOXPLOT_BASE_SIZE) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 0, hjust = 0.5)
    )
  
  file_stub <- paste(regime, measure, experiment_name, "boxplot", sep = "_")
  file_pdf  <- file.path(FIG_DIR, paste0(file_stub, ".pdf"))
  
  ggsave(
    filename = file_pdf,
    plot = p,
    width = BOXPLOT_WIDTH,
    height = BOXPLOT_HEIGHT
  )
  
  invisible(file_pdf)
}

# ============================================================
# Parallel setup
# ============================================================

WORKER_EXPORTS <- c(
  "ALPHA_Y", "BETA_Y", "RHO", "PROB_CLIP", "PI_HAT_FLOOR",
  "M_MIN", "M_CAP", "M_FACTOR",
  "H_C_MIN", "H_C_MAX", "H_GRID_LENGTH",
  "BISE_HALF_WIDTH_DENOM",
  "REGIME_LEVELS", "MEASURE_LEVELS", "ESTIMATOR_KEYS", "ESTIMATOR_LABELS",
  "F_true", "rY",
  "generate_auxiliary_from_Y",
  "propensity_fun",
  "simulate_dataset",
  "estimate_pi_hat",
  "ipw_ecdf_at",
  "bernstein_smooth",
  "m_candidates_of_n",
  "bernstein_term1_exact",
  "bernstein_tail_integrals",
  "cv_score_bernstein",
  "ipw_gaussian_kde_cdf",
  "h_candidates_of_n",
  "cv_score_kde",
  "boundary_mask_of_n",
  "compute_ISE",
  "compute_BISE",
  "one_replication"
)

WORKER_EXPORTS <- unique(c(
  WORKER_EXPORTS,
  "RISK_NUMINT_METHOD",
  "RISK_NUMINT_REL_TOL",
  "RISK_NUMINT_ABS_TOL",
  "RISK_NUMINT_SUBDIVISIONS",
  "RISK_NUMINT_STOP_ON_ERROR",
  "RISK_ADAPTINT_MAX_EVAL",
  "RISK_ASIMPSON_MAX_DEPTH",
  "normalize_numint_method",
  "make_risk_integrand",
  "adaptive_simpson_1d",
  "numerical_integral_1d"
))

setup_parallel_cluster <- function(num_cores = N_CORES) {
  cl <- parallel::makeCluster(num_cores)
  parallel::clusterExport(cl, varlist = WORKER_EXPORTS, envir = .GlobalEnv)
  cl
}

# ============================================================
# Generic experiment runner
# ============================================================

run_experiment <- function(cl, experiment_name, settings, n_fun, beta0_fun, beta1_fun,
                           seed_offset = 0L) {
  y_grid      <- seq(0, 1, length.out = GRID_POINTS)
  F_true_grid <- F_true(y_grid)
  
  long_list <- vector("list", length(settings))
  
  for (i in seq_along(settings)) {
    setting <- settings[[i]]
    
    n_here     <- n_fun(setting)
    beta0_here <- beta0_fun(setting)
    beta1_here <- beta1_fun(setting)
    
    m_grid_here        <- m_candidates_of_n(n_here)
    h_grid_here        <- h_candidates_of_n(n_here)
    boundary_mask_here <- boundary_mask_of_n(n_here, y_grid)
    
    if (experiment_name == "sample_size") {
      setting_value <- as.integer(n_here)
      setting_label <- as.character(n_here)
      progress_msg <- sprintf("Running %s experiment for n = %d", experiment_name, n_here)
    } else {
      setting_value <- as.integer(round(100 * setting))
      setting_label <- paste0(setting_value, "%")
      progress_msg <- sprintf(
        "Running %s experiment for target missing rate = %s",
        experiment_name,
        setting_label
      )
    }
    
    if (VERBOSE) {
      message(progress_msg)
    }
    
    seeds <- as.integer(SEED + seed_offset + 100000L * i + seq_len(MC_REPS))
    
    res_list <- parallel::parLapplyLB(
      cl,
      X = seeds,
      fun = function(seed, n_here, beta0_here, beta1_here,
                     y_grid, F_true_grid, m_grid_here, h_grid_here, boundary_mask_here) {
        set.seed(seed)
        one_replication(
          n = n_here,
          beta0 = beta0_here,
          beta1 = beta1_here,
          y_grid = y_grid,
          F_true_grid = F_true_grid,
          m_grid = m_grid_here,
          h_grid = h_grid_here,
          boundary_mask = boundary_mask_here
        )
      },
      n_here = n_here,
      beta0_here = beta0_here,
      beta1_here = beta1_here,
      y_grid = y_grid,
      F_true_grid = F_true_grid,
      m_grid_here = m_grid_here,
      h_grid_here = h_grid_here,
      boundary_mask_here = boundary_mask_here
    )
    
    res_mat <- do.call(rbind, res_list)
    
    long_list[[i]] <- results_matrix_to_long(
      res_mat = res_mat,
      experiment_name = experiment_name,
      setting_value = setting_value,
      setting_label = setting_label,
      n_used = n_here
    )
  }
  
  do.call(rbind, long_list)
}

# ============================================================
# Generate all tables and boxplots from one long results table
# ============================================================

generate_outputs <- function(long_df, experiment_name) {
  for (regime in REGIME_LEVELS) {
    for (measure in MEASURE_LEVELS) {
      summary_df <- make_summary_table(
        long_df = long_df,
        experiment_name = experiment_name,
        regime = regime,
        measure = measure,
        scale_factor = TABLE_STAT_MULTIPLIER
      )
      
      save_summary_outputs(
        summary_df = summary_df,
        experiment_name = experiment_name,
        regime = regime,
        measure = measure
      )
      
      save_boxplot(
        long_df = long_df,
        experiment_name = experiment_name,
        regime = regime,
        measure = measure
      )
    }
  }
}

# ============================================================
# LSCV diagnostics subsection
# ============================================================

# --- Hyperparameters for this subsection (edit here only) ---
RUN_LSCV_DIAGNOSTICS <- TRUE
LSCV_DIAG_SEED       <- 777
LSCV_FIG_WIDTH       <- 11
LSCV_FIG_HEIGHT      <- 8.5
LSCV_LINE_WIDTH      <- 2
LSCV_POINT_CEX       <- 0.8
LSCV_PAR_MAR         <- c(4.2, 4.2, 3.0, 1.2)
LSCV_MAIN_CEX        <- 1.0
LSCV_AXIS_CEX        <- 0.95
LSCV_LAB_CEX         <- 1.0

LSCV_FIG_DIR <- file.path(FIG_DIR, "LSCV_figures")
dir_create_safe(LSCV_FIG_DIR)

save_one_lscv_diagnostic <- function(config_name,
                                     file_stub,
                                     n_here,
                                     beta0_here,
                                     beta1_here,
                                     seed_here = LSCV_DIAG_SEED) {
  set.seed(seed_here)
  
  y_grid <- seq(0, 1, length.out = GRID_POINTS)
  
  dat   <- simulate_dataset(n = n_here, beta0 = beta0_here, beta1 = beta1_here, rho = RHO)
  Y     <- dat$Y
  X     <- dat$X
  pi    <- dat$pi
  delta <- dat$delta
  
  # Weights
  w_pseudo   <- delta / pi
  pi_hat     <- estimate_pi_hat(X, delta, pi_floor = PI_HAT_FLOOR)
  w_feasible <- delta / pi_hat
  
  # Candidate grids
  m_grid <- m_candidates_of_n(n_here)
  h_grid <- h_candidates_of_n(n_here)
  
  # LSCV values
  cv_m_pseudo <- vapply(
    m_grid,
    cv_score_bernstein,
    numeric(1),
    Y = Y,
    w = w_pseudo,
    n_total = n_here
  )
  m_hat_pseudo <- m_grid[which.min(cv_m_pseudo)]
  
  cv_m_feasible <- vapply(
    m_grid,
    cv_score_bernstein,
    numeric(1),
    Y = Y,
    w = w_feasible,
    n_total = n_here
  )
  m_hat_feasible <- m_grid[which.min(cv_m_feasible)]
  
  cv_h_pseudo <- vapply(
    h_grid,
    cv_score_kde,
    numeric(1),
    Y = Y,
    w = w_pseudo,
    n_total = n_here,
    y_grid = y_grid
  )
  h_hat_pseudo <- h_grid[which.min(cv_h_pseudo)]
  
  cv_h_feasible <- vapply(
    h_grid,
    cv_score_kde,
    numeric(1),
    Y = Y,
    w = w_feasible,
    n_total = n_here,
    y_grid = y_grid
  )
  h_hat_feasible <- h_grid[which.min(cv_h_feasible)]
  
  pdf(
    file = file.path(LSCV_FIG_DIR, paste0(file_stub, ".pdf")),
    width = LSCV_FIG_WIDTH,
    height = LSCV_FIG_HEIGHT
  )
  
  op <- par(no.readonly = TRUE)
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)
  
  par(mfrow = c(2, 2), mar = LSCV_PAR_MAR, cex.main = LSCV_MAIN_CEX,
      cex.axis = LSCV_AXIS_CEX, cex.lab = LSCV_LAB_CEX)
  
  plot(
    h_grid, cv_h_pseudo,
    type = "b", lwd = LSCV_LINE_WIDTH, pch = 16, cex = LSCV_POINT_CEX,
    xlab = "h (bandwidth)", ylab = "LSCV(h)",
    main = paste0(config_name, "\nPseudo I-IPW KDE")
  )
  abline(v = h_hat_pseudo, lty = 2)
  legend("topright", legend = sprintf("h* = %.6f", h_hat_pseudo), bty = "n")
  
  plot(
    h_grid, cv_h_feasible,
    type = "b", lwd = LSCV_LINE_WIDTH, pch = 16, cex = LSCV_POINT_CEX,
    xlab = "h (bandwidth)", ylab = "LSCV(h)",
    main = paste0(config_name, "\nFeasible I-IPW KDE")
  )
  abline(v = h_hat_feasible, lty = 2)
  legend("topright", legend = sprintf("h* = %.6f", h_hat_feasible), bty = "n")
  
  plot(
    m_grid, cv_m_pseudo,
    type = "b", lwd = LSCV_LINE_WIDTH, pch = 16, cex = LSCV_POINT_CEX,
    xlab = "m (Bernstein degree)", ylab = "LSCV(m)",
    main = paste0(config_name, "\nPseudo Bernstein")
  )
  abline(v = m_hat_pseudo, lty = 2)
  legend("topright", legend = sprintf("m* = %d", m_hat_pseudo), bty = "n")
  
  plot(
    m_grid, cv_m_feasible,
    type = "b", lwd = LSCV_LINE_WIDTH, pch = 16, cex = LSCV_POINT_CEX,
    xlab = "m (Bernstein degree)", ylab = "LSCV(m)",
    main = paste0(config_name, "\nFeasible Bernstein")
  )
  abline(v = m_hat_feasible, lty = 2)
  legend("topright", legend = sprintf("m* = %d", m_hat_feasible), bty = "n")
  
  invisible(
    list(
      m_grid = m_grid,
      h_grid = h_grid,
      cv_m_pseudo = cv_m_pseudo,
      cv_m_feasible = cv_m_feasible,
      cv_h_pseudo = cv_h_pseudo,
      cv_h_feasible = cv_h_feasible,
      m_hat_pseudo = m_hat_pseudo,
      m_hat_feasible = m_hat_feasible,
      h_hat_pseudo = h_hat_pseudo,
      h_hat_feasible = h_hat_feasible
    )
  )
}

run_lscv_diagnostics_section <- function() {
  if (!RUN_LSCV_DIAGNOSTICS) {
    return(invisible(NULL))
  }
  
  t_lscv <- Sys.time()
  message("\n=== LSCV diagnostics subsection ===")
  
  # Sample-size experiment configurations
  for (n_here in SAMPLE_SIZE_GRID) {
    config_name <- paste0("Sample-size experiment: n = ", n_here)
    file_stub   <- paste0("LSCV_sample_size_n_", n_here)
    
    save_one_lscv_diagnostic(
      config_name = config_name,
      file_stub   = file_stub,
      n_here      = n_here,
      beta0_here  = BETA0_BENCHMARK,
      beta1_here  = BETA1_BENCHMARK,
      seed_here   = LSCV_DIAG_SEED + n_here
    )
  }
  
  # Missing-rate experiment configurations
  for (miss in MISSING_RATE_GRID) {
    miss_pct   <- as.integer(round(100 * miss))
    beta0_here <- solve_beta0_for_missing_rate(
      target_missing_rate = miss,
      beta1 = BETA1_BENCHMARK
    )
    
    config_name <- paste0("Missing-rate experiment: ", miss_pct, "% missing")
    file_stub   <- paste0("LSCV_missing_rate_", miss_pct, "pct")
    
    save_one_lscv_diagnostic(
      config_name = config_name,
      file_stub   = file_stub,
      n_here      = N_FIXED_MISSING,
      beta0_here  = beta0_here,
      beta1_here  = BETA1_BENCHMARK,
      seed_here   = LSCV_DIAG_SEED + 1000 + miss_pct
    )
  }
  
  timer_report("LSCV diagnostics subsection", t_lscv)
  invisible(TRUE)
}

run_lscv_diagnostics_section()

# ============================================================
# Main
# ============================================================

main <- function() {
  global_start <- Sys.time()
  message("\n=== GLOBAL TIMER STARTED ===")
  message(sprintf("Detected cores: %d | Using cores: %d", N_CORES_AVAILABLE, N_CORES))
  message(sprintf(
    "Benchmark beta0 = %.10f | Benchmark beta1 = %.10f",
    BETA0_BENCHMARK, BETA1_BENCHMARK
  ))
  
  beta0_missing_grid <- vapply(
    MISSING_RATE_GRID,
    solve_beta0_for_missing_rate,
    numeric(1),
    beta1 = BETA1_BENCHMARK
  )
  
  missing_calibration <- data.frame(
    missing_rate = 100 * MISSING_RATE_GRID,
    beta0 = beta0_missing_grid,
    expected_observation_rate = vapply(
      beta0_missing_grid,
      expected_observation_prob,
      numeric(1),
      beta1 = BETA1_BENCHMARK
    ),
    expected_missing_rate = 1 - vapply(
      beta0_missing_grid,
      expected_observation_prob,
      numeric(1),
      beta1 = BETA1_BENCHMARK
    )
  )
  
  message("Calibrated intercepts for the missing-rate experiment:")
  print(missing_calibration, row.names = FALSE)
  
  cl <- NULL
  on.exit({
    if (!is.null(cl)) {
      try(parallel::stopCluster(cl), silent = TRUE)
    }
  }, add = TRUE)
  
  t_cluster <- Sys.time()
  cl <- setup_parallel_cluster(num_cores = N_CORES)
  timer_report("Parallel cluster initialization", t_cluster)
  
  # Sample-size experiment
  t_sample <- Sys.time()
  
  long_sample_size <- run_experiment(
    cl = cl,
    experiment_name = "sample_size",
    settings = as.list(SAMPLE_SIZE_GRID),
    n_fun = function(setting) as.integer(setting),
    beta0_fun = function(setting) BETA0_BENCHMARK,
    beta1_fun = function(setting) BETA1_BENCHMARK,
    seed_offset = 0L
  )
  
  generate_outputs(long_df = long_sample_size, experiment_name = "sample_size")
  timer_report("Sample-size experiment", t_sample)
  
  # Missing-rate experiment
  t_missing <- Sys.time()
  
  long_missing_rate <- run_experiment(
    cl = cl,
    experiment_name = "missing_rate",
    settings = as.list(MISSING_RATE_GRID),
    n_fun = function(setting) N_FIXED_MISSING,
    beta0_fun = function(setting) solve_beta0_for_missing_rate(setting, beta1 = BETA1_BENCHMARK),
    beta1_fun = function(setting) BETA1_BENCHMARK,
    seed_offset = 10000000L
  )
  
  generate_outputs(long_df = long_missing_rate, experiment_name = "missing_rate")
  timer_report("Missing-rate experiment", t_missing)
  
  t_stop <- Sys.time()
  parallel::stopCluster(cl)
  cl <- NULL
  timer_report("Parallel cluster shutdown", t_stop)
  
  timer_report("TOTAL runtime", global_start)
  message("=== GLOBAL TIMER ENDED ===\n")
  
  invisible(
    list(
      sample_size = long_sample_size,
      missing_rate = long_missing_rate
    )
  )
}

main()

