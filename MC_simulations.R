## ============================================================
## Monte Carlo comparison: Bernstein-smoothed vs. unsmoothed vs. IPW–KDE CDF
##  - Study 1: Pseudo-estimators (known propensity scores)
##  - Study 2: Feasible estimators (π estimated from discrete X)
##
## Outputs: two tables with, for each n,
##   mean ISE and var ISE for (Unsmoothed, Bernstein, KDE)
##   (KDE uses a pluggable kernel CDF; default = Gaussian)
## ============================================================

## -----------------------
## Hyperparameters (edit me)
## -----------------------

## Parallel settings
library(parallel)
n_cores_avail <- parallel::detectCores(logical = TRUE)
n_cores       <- max(1, n_cores_avail - 1)   # leave 1 core free

set.seed(2025)

n_grid        <- c(25, 50, 100, 200, 400, 800, 1600, 3200, 6400) # sample sizes
MC            <- 100                  # Monte Carlo replications per n
y_dist        <- list(name = "beta", shape1 = 2, shape2 = 5)  # Y ~ Beta(a,b) on [0,1]
x_levels      <- c(0, 1)               # discrete auxiliary variable takes values in {0,1}
px            <- 0.5                   # P(X=1); P(X=0)=1-px
pi_vals       <- c(`0` = 0.6, `1` = 0.9)  # true propensity π(x): x=0 -> 0.6, x=1 -> 0.9
pi_min        <- 0.05                  # lower bound for feasible π-hat to avoid division by ~0
grid_points   <- 512                   # integration grid for ISE, CV ∈ [0,1]
m_min         <- 1                     # minimum Bernstein degree in search
m_cap         <- 300                   # absolute cap on m
m_factor      <- 3                     # m_max(n) = min(floor(m_factor * n^(2/3)), m_cap, n)
verbose       <- TRUE                  # print progress

## KDE kernel choice (CDF form). Options: "gaussian","epanechnikov","triangular","uniform"
kde_kernel    <- "gaussian"

# Paths
path <- "C://Users//fred1//Dropbox//Gharbi_Jedidi_Khardani_Ouimet_2025_missing_data"
setwd(path)

TABLE_DIR  <- file.path(getwd(), "simulation_tables")
FIG_DIR    <- file.path(getwd(), "simulation_figures")

dir.create(TABLE_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR,   showWarnings = FALSE, recursive = TRUE)

## -----------------------
## Utilities
## -----------------------

## True CDF of Y on [0,1]
F_true <- function(y) {
  if (y_dist$name == "beta") {
    return(pbeta(y, shape1 = y_dist$shape1, shape2 = y_dist$shape2))
  } else {
    stop("Only Beta implemented; provide your F_true.")
  }
}

## Sample Y
rY <- function(n) {
  if (y_dist$name == "beta") {
    return(rbeta(n, shape1 = y_dist$shape1, shape2 = y_dist$shape2))
  } else {
    stop("Only Beta implemented; provide your rY.")
  }
}

## Draw X ∈ {0,1}
rX <- function(n) rbinom(n, size = 1, prob = px)

## Propensity π(x) (known for pseudo study)
pi_fun <- function(x) {
  vapply(x, function(xx) pi_vals[as.character(xx)], numeric(1))
}

## Generate one dataset
simulate_dataset <- function(n) {
  X  <- rX(n)
  Y  <- rY(n)
  pi <- pi_fun(X)
  delta <- rbinom(n, size = 1, prob = pi)
  list(Y = Y, X = X, pi = pi, delta = delta)
}

## Estimate π-hat_i for feasible study (discrete X)
estimate_pi_hat <- function(X, delta, pi_min = 0.05) {
  out <- numeric(length(X))
  for (val in unique(X)) {
    idx <- which(X == val)
    num <- sum(delta[idx] == 1)
    den <- length(idx)
    pih <- if (den > 0) num / den else 0
    pih <- max(pih, pi_min)  # floor to avoid division by ~0
    out[idx] <- pih
  }
  out
}

## Weighted IPW ECDF evaluated at a vector of points 't_points'
## Fhat(t) = (1/n) * sum_i w_i * 1{Y_i ≤ t}, where w_i already includes delta_i/π_i[hat]
ipw_ecdf_at <- function(Y, w, t_points, n_total) {
  ## Only observed contribute because w_i = 0 when delta_i = 0
  obs <- which(w > 0)
  if (length(obs) == 0L) return(rep(0, length(t_points)))
  Yobs <- Y[obs]
  wobs <- w[obs]
  ord  <- order(Yobs)
  Yso  <- Yobs[ord]
  wso  <- wobs[ord]
  cs   <- cumsum(wso)
  idx  <- findInterval(t_points, Yso, left.open = FALSE, rightmost.closed = TRUE)
  Fv   <- numeric(length(idx))
  pos  <- idx > 0L
  Fv[pos] <- cs[idx[pos]]
  (1 / n_total) * Fv
}

## Bernstein smoothing: B_m(F)(y) = sum_{k=0}^m F(k/m) * dbinom(k; m, y)
bernstein_smooth <- function(F_at_knots, y_points, m) {
  ks <- 0:m
  sapply(y_points, function(y) sum(F_at_knots * dbinom(ks, size = m, prob = y)))
}

## Build m-candidate set for a given n
m_candidates_of_n <- function(n) {
  m_max <- min(floor(m_factor * n^(2/3)), m_cap, n)
  if (m_max < m_min) m_max <- m_min
  seq.int(m_min, m_max, by = 1L)
}

## Numerical integral over [0,1] via uniform grid mean
integral_mean <- function(values) mean(values)

## Exact integral helpers for Bernstein LSCV (Lebesgue ISE)
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

## LSCV for pseudo (known π) under Lebesgue ISE
cv_score_pseudo <- function(m, Y, w, n_total, y_grid) {
  t_k     <- seq(0, 1, length.out = m + 1L)
  F_knots <- ipw_ecdf_at(Y, w, t_k, n_total = n_total)
  
  ## Term 1: ∫ B_m(F_n)^2 dy (exact)
  term1 <- bernstein_term1_exact(F_knots, m)
  
  ## Term 2: (2/n) * sum_i w_i * ∫_{Y_i}^1 B_m(F_{-i})(y) dy
  sumW_leq_tk <- n_total * F_knots
  obs <- which(w > 0)
  if (length(obs) == 0L) return(term1)
  
  acc <- 0
  for (i in obs) {
    ind_i     <- as.numeric(t_k >= Y[i])                 # 1{Y_i ≤ k/m}
    F_minus_k <- (sumW_leq_tk - w[i] * ind_i) / (n_total - 1)
    tail_i    <- bernstein_tail_integrals(Y[i], m)       # ∫_{Y_i}^1 b_{m,k}(y) dy
    acc <- acc + w[i] * sum(F_minus_k * tail_i)
  }
  term2 <- (2 / n_total) * acc
  
  term1 - term2
}

## Feasible case: identical structure, just different weights w = δ/π̂
cv_score_feasible <- function(m, Y, w, n_total, y_grid) {
  cv_score_pseudo(m, Y, w, n_total, y_grid)
}

## Compute ISE on grid for a given estimator (vector on y_grid)
compute_ISE <- function(est_on_grid, F_true_on_grid) {
  integral_mean((est_on_grid - F_true_on_grid)^2)
}

## -------- NEW: Kernel-CDF factory + generic IPW–KDE CDF on [0,1] --------

## Kernel CDFs (standardized): G(u) = ∫_{-∞}^u K(t) dt
get_kernel_cdf <- function(name = c("gaussian","epanechnikov","triangular","uniform")) {
  name <- match.arg(name)
  switch(name,
         gaussian = function(u) pnorm(u),
         epanechnikov = function(u) {
           out <- numeric(length(u))
           out[u <= -1] <- 0
           out[u >=  1] <- 1
           mid <- (u > -1) & (u < 1)
           um  <- u[mid]
           out[mid] <- 0.5 + 0.75*(um - um^3/3)  # = 0.75*(um - um^3/3 + 2/3)
           out
         },
         triangular = function(u) {
           out <- numeric(length(u))
           out[u <= -1] <- 0
           out[u >=  1] <- 1
           mid <- (u > -1) & (u < 1)
           um  <- u[mid]
           out[mid] <- 0.5 + um - 0.5*um*abs(um)
           out
         },
         uniform = function(u) {
           pmin(pmax((u + 1)/2, 0), 1)
         }
  )
}

## IPW–KDE CDF on [0,1] with boundary renormalization (generic kernel via its CDF)
## F_h(y) = [ (1/n) Σ w_i {G((y - Y_i)/h) - G(-Y_i/h)} ] / [ (1/n) Σ w_i {G((1 - Y_i)/h) - G(-Y_i/h)} ]
ipw_kde_cdf <- function(Y, w, y_points, n_total, h, kernel = kde_kernel) {
  obs <- which(w > 0 & is.finite(Y))
  if (!length(obs)) return(rep(0, length(y_points)))
  Yo <- Y[obs]; wo <- w[obs]
  G  <- get_kernel_cdf(kernel)
  base <- G((-Yo) / h)
  top1 <- G((1 - Yo) / h)
  denom <- sum(wo * (top1 - base)) / n_total
  if (denom <= 0) denom <- 1
  
  vapply(y_points, function(y) {
    zy   <- (y - Yo) / h
    num  <- sum(wo * (G(zy) - base)) / n_total
    val  <- num / denom
    if (val < 0) 0 else if (val > 1) 1 else val
  }, numeric(1))
}

## Bandwidth candidate grid: h = c * n^{-1/5}, c in [c_min, c_max]
h_candidates_of_n <- function(n, c_min = 0.2, c_max = 1.5, num = 25L) {
  c_seq <- seq(c_min, c_max, length.out = num)
  c_seq * n^(-1/5)
}

## LSCV for KDE-CDF (generic kernel; numerical on y_grid):
## LSCV(h) = ∫ F_h(y)^2 dy - (2/n) Σ_i w_i ∫_{Y_i}^1 F_{-i,h}(y) dy
cv_score_kde_core <- function(Y, w, n_total, y_grid, h, kernel = kde_kernel) {
  Glen <- length(y_grid)
  obs <- which(w > 0 & is.finite(Y))
  if (!length(obs)) return(0)
  
  Yo <- Y[obs]; wo <- w[obs]
  G  <- get_kernel_cdf(kernel)
  base <- G((-Yo)/h)
  top1 <- G((1 - Yo)/h)
  D    <- sum(wo * (top1 - base)) / n_total
  
  # Build N_vec(y_j) on the grid
  N_vec <- numeric(Glen)
  for (i in seq_along(Yo)) {
    add_j <- G((y_grid - Yo[i]) / h) - base[i]
    N_vec <- N_vec + wo[i] * add_j
  }
  N_vec <- N_vec / n_total
  
  F_vec <- N_vec / ifelse(D > 0, D, 1)
  F_vec <- pmin(pmax(F_vec, 0), 1)
  term1 <- mean(F_vec^2)
  
  # Term 2: leave-one-out tail integrals numerically
  term2_acc <- 0
  w_over_n  <- wo / n_total
  for (i in seq_along(Yo)) {
    b_i   <- w_over_n[i] * (top1[i] - base[i])
    den_i <- D - b_i
    add_j <- w_over_n[i] * (G((y_grid - Yo[i]) / h) - base[i])
    N_mi  <- N_vec - add_j
    F_mi  <- N_mi / ifelse(den_i > 0, den_i, 1)
    F_mi  <- pmin(pmax(F_mi, 0), 1)
    sel   <- y_grid >= Yo[i]
    tail  <- sum(F_mi[sel]) / Glen
    term2_acc <- term2_acc + wo[i] * tail
  }
  term2 <- (2 / n_total) * term2_acc
  
  term1 - term2
}
cv_score_kde_pseudo   <- function(h, Y, w, n_total, y_grid) cv_score_kde_core(Y, w, n_total, y_grid, h, kernel = kde_kernel)
cv_score_kde_feasible <- function(h, Y, w, n_total, y_grid) cv_score_kde_core(Y, w, n_total, y_grid, h, kernel = kde_kernel)

## -----------------------
## One replication (returns ISEs only)
## -----------------------

one_replication <- function(n, y_grid, F_true_grid) {
  dat   <- simulate_dataset(n)
  Y     <- dat$Y
  X     <- dat$X
  pi    <- dat$pi
  delta <- dat$delta
  
  ## Pseudo weights (known π)
  w_pseudo <- delta / pi
  
  ## Feasible weights (π-hat from discrete X)
  pi_hat   <- estimate_pi_hat(X, delta, pi_min = pi_min)
  w_feas   <- delta / pi_hat
  
  ## --- PSEUDO: UnsMoothed ---
  F_uns_pseudo_grid <- ipw_ecdf_at(Y, w_pseudo, y_grid, n_total = n)
  ISE_uns_pseudo    <- compute_ISE(F_uns_pseudo_grid, F_true_grid)
  
  ## --- PSEUDO: Bernstein w/ LSCV m ---
  m_grid  <- m_candidates_of_n(n)
  cv_vals <- vapply(m_grid, cv_score_pseudo, numeric(1),
                    Y = Y, w = w_pseudo, n_total = n, y_grid = y_grid)
  m_hat <- m_grid[which.min(cv_vals)]
  t_k   <- seq(0, 1, length.out = m_hat + 1L)
  F_knots <- ipw_ecdf_at(Y, w_pseudo, t_k, n_total = n)
  F_smooth_pseudo_grid <- bernstein_smooth(F_knots, y_grid, m_hat)
  ISE_sm_pseudo  <- compute_ISE(F_smooth_pseudo_grid, F_true_grid)
  
  ## --- PSEUDO: KDE-CDF (Gaussian, by default) w/ LSCV h ---
  h_grid <- h_candidates_of_n(n)
  cv_h   <- vapply(h_grid, cv_score_kde_pseudo, numeric(1),
                   Y = Y, w = w_pseudo, n_total = n, y_grid = y_grid)
  h_hat  <- h_grid[which.min(cv_h)]
  F_kde_pseudo_grid <- ipw_kde_cdf(Y, w_pseudo, y_grid, n_total = n, h = h_hat, kernel = kde_kernel)
  ISE_kde_pseudo    <- compute_ISE(F_kde_pseudo_grid, F_true_grid)
  
  ## --- FEASIBLE: UnsMoothed ---
  F_uns_feas_grid <- ipw_ecdf_at(Y, w_feas, y_grid, n_total = n)
  ISE_uns_feas    <- compute_ISE(F_uns_feas_grid, F_true_grid)
  
  ## --- FEASIBLE: Bernstein w/ LSCV m ---
  m_grid2 <- m_candidates_of_n(n)
  cv_vals2 <- vapply(m_grid2, cv_score_feasible, numeric(1),
                     Y = Y, w = w_feas, n_total = n, y_grid = y_grid)
  m_hat2 <- m_grid2[which.min(cv_vals2)]
  t_k2   <- seq(0, 1, length.out = m_hat2 + 1L)
  F_knots2 <- ipw_ecdf_at(Y, w_feas, t_k2, n_total = n)
  F_smooth_feas_grid <- bernstein_smooth(F_knots2, y_grid, m_hat2)
  ISE_sm_feas  <- compute_ISE(F_smooth_feas_grid, F_true_grid)
  
  ## --- FEASIBLE: KDE-CDF (Gaussian) w/ LSCV h ---
  h_grid2 <- h_candidates_of_n(n)
  cv_h2   <- vapply(h_grid2, cv_score_kde_feasible, numeric(1),
                    Y = Y, w = w_feas, n_total = n, y_grid = y_grid)
  h_hat2  <- h_grid2[which.min(cv_h2)]
  F_kde_feas_grid <- ipw_kde_cdf(Y, w_feas, y_grid, n_total = n, h = h_hat2, kernel = kde_kernel)
  ISE_kde_feas    <- compute_ISE(F_kde_feas_grid, F_true_grid)
  
  list(
    ISE_uns_pseudo = ISE_uns_pseudo,
    ISE_sm_pseudo  = ISE_sm_pseudo,
    ISE_kde_pseudo = ISE_kde_pseudo,
    
    ISE_uns_feas   = ISE_uns_feas,
    ISE_sm_feas    = ISE_sm_feas,
    ISE_kde_feas   = ISE_kde_feas
  )
}

## -----------------------
## TEST: LSCV curves for Bernstein (m) and KDE (h) for all n in n_grid
## Save plots + print argmin markers
## -----------------------

for (idx in seq_along(n_grid)) {
  n_test <- n_grid[idx]
  set.seed(12345 + idx)   # different seed per n for reproducibility
  
  # Simulate one dataset for this n
  dat_test <- simulate_dataset(n_test)
  Yt <- dat_test$Y; Xt <- dat_test$X; pit <- dat_test$pi; deltat <- dat_test$delta
  
  # Weights for pseudo and feasible
  w_pseudo_t <- deltat / pit
  pi_hat_t   <- estimate_pi_hat(Xt, deltat, pi_min = pi_min)
  w_feas_t   <- deltat / pi_hat_t
  
  # Grid and candidate sets
  y_grid_test <- seq(0, 1, length.out = grid_points)
  
  ## -------- Bernstein (m) --------
  m_grid_test <- m_candidates_of_n(n_test)
  
  # LSCV over m (pseudo)
  cv_pseudo_m <- vapply(m_grid_test, cv_score_pseudo, numeric(1),
                        Y = Yt, w = w_pseudo_t, n_total = n_test, y_grid = y_grid_test)
  m_hat_pseudo_test <- m_grid_test[which.min(cv_pseudo_m)]
  
  # LSCV over m (feasible)
  cv_feas_m <- vapply(m_grid_test, cv_score_feasible, numeric(1),
                      Y = Yt, w = w_feas_t, n_total = n_test, y_grid = y_grid_test)
  m_hat_feas_test <- m_grid_test[which.min(cv_feas_m)]
  
  # Simple convexity diagnostics: fraction of nonnegative second differences
  conv_pseudo_m <- mean(diff(diff(cv_pseudo_m)) >= 0)
  conv_feas_m   <- mean(diff(diff(cv_feas_m)) >= 0)
  
  cat(sprintf("[LSCV-m] n=%d | pseudo: m*=%d, frac nonneg 2nd diff=%.3f\n",
              n_test, m_hat_pseudo_test, conv_pseudo_m))
  cat(sprintf("[LSCV-m] n=%d | feasible: m*=%d, frac nonneg 2nd diff=%.3f\n",
              n_test, m_hat_feas_test, conv_feas_m))
  
  # Save Bernstein LSCV plots
  pdf(file.path(FIG_DIR, sprintf("LSCV_Bernstein_pseudo_n%d.pdf", n_test)), width = 6, height = 4)
  plot(m_grid_test, cv_pseudo_m, type = "b",
       xlab = "m (Bernstein degree)", ylab = "LSCV(m)",
       main = sprintf("Bernstein LSCV (pseudo, n=%d)", n_test))
  abline(v = m_hat_pseudo_test, lty = 2)
  legend("topright", legend = sprintf("m* = %d", m_hat_pseudo_test), bty = "n")
  dev.off()
  
  pdf(file.path(FIG_DIR, sprintf("LSCV_Bernstein_feasible_n%d.pdf", n_test)), width = 6, height = 4)
  plot(m_grid_test, cv_feas_m, type = "b",
       xlab = "m (Bernstein degree)", ylab = "LSCV(m)",
       main = sprintf("Bernstein LSCV (feasible, n=%d)", n_test))
  abline(v = m_hat_feas_test, lty = 2)
  legend("topright", legend = sprintf("m* = %d", m_hat_feas_test), bty = "n")
  dev.off()
  
  ## -------- KDE CDF (Gaussian by default) bandwidth h --------
  h_grid_test <- h_candidates_of_n(n_test)
  
  # LSCV over h (pseudo)
  cv_pseudo_h <- vapply(h_grid_test, cv_score_kde_pseudo, numeric(1),
                        Y = Yt, w = w_pseudo_t, n_total = n_test, y_grid = y_grid_test)
  h_hat_pseudo_test <- h_grid_test[which.min(cv_pseudo_h)]
  
  # LSCV over h (feasible)
  cv_feas_h <- vapply(h_grid_test, cv_score_kde_feasible, numeric(1),
                      Y = Yt, w = w_feas_t, n_total = n_test, y_grid = y_grid_test)
  h_hat_feas_test <- h_grid_test[which.min(cv_feas_h)]
  
  conv_pseudo_h <- mean(diff(diff(cv_pseudo_h)) >= 0)
  conv_feas_h   <- mean(diff(diff(cv_feas_h)) >= 0)
  
  cat(sprintf("[LSCV-h] n=%d | pseudo: h*=%0.4f, frac nonneg 2nd diff=%.3f\n",
              n_test, h_hat_pseudo_test, conv_pseudo_h))
  cat(sprintf("[LSCV-h] n=%d | feasible: h*=%0.4f, frac nonneg 2nd diff=%.3f\n",
              n_test, h_hat_feas_test, conv_feas_h))
  
  # Save KDE LSCV plots
  pdf(file.path(FIG_DIR, sprintf("LSCV_KDE_pseudo_n%d.pdf", n_test)), width = 6, height = 4)
  plot(h_grid_test, cv_pseudo_h, type = "b",
       xlab = "h (bandwidth)", ylab = "LSCV(h)",
       main = sprintf("KDE CDF LSCV (pseudo, n=%d, kernel=%s)", n_test, kde_kernel))
  abline(v = h_hat_pseudo_test, lty = 2)
  legend("topright", legend = sprintf("h* = %.4f", h_hat_pseudo_test), bty = "n")
  dev.off()
  
  pdf(file.path(FIG_DIR, sprintf("LSCV_KDE_feasible_n%d.pdf", n_test)), width = 6, height = 4)
  plot(h_grid_test, cv_feas_h, type = "b",
       xlab = "h (bandwidth)", ylab = "LSCV(h)",
       main = sprintf("KDE CDF LSCV (feasible, n=%d, kernel=%s)", n_test, kde_kernel))
  abline(v = h_hat_feas_test, lty = 2)
  legend("topright", legend = sprintf("h* = %.4f", h_hat_feas_test), bty = "n")
  dev.off()
}

## -----------------------
## One-shot CDF plots per n (pseudo & feasible; Unsmoothed vs Bernstein vs KDE)
## Saved in FIG_DIR as CDF_pseudo_n{n}.pdf and CDF_feasible_n{n}.pdf
## -----------------------

if (!exists("FIG_DIR")) {
  FIG_DIR <- file.path(getwd(), "simulation_figures")
  dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
}

for (n in n_grid) {
  set.seed(4242 + n)  # reproducible per n
  
  ## Simulate one dataset and build weights
  dat   <- simulate_dataset(n)
  Y     <- dat$Y
  X     <- dat$X
  pi    <- dat$pi
  delta <- dat$delta
  
  w_pseudo <- delta / pi
  pi_hat   <- estimate_pi_hat(X, delta, pi_min = pi_min)
  w_feas   <- delta / pi_hat
  
  ## Use the global y_grid if it exists; otherwise define a default
  if (!exists("y_grid")) y_grid <- seq(0, 1, length.out = 512)
  
  ## ========= PSEUDO: compute the three CDFs =========
  # Unsmoothed
  F_uns_pseudo <- ipw_ecdf_at(Y, w_pseudo, y_grid, n_total = n)
  F_uns_pseudo <- pmin(pmax(F_uns_pseudo, 0), 1)
  
  # Bernstein (LSCV over m)
  m_grid_here   <- m_candidates_of_n(n)
  cv_pseudo_m   <- vapply(m_grid_here, cv_score_pseudo, numeric(1),
                          Y = Y, w = w_pseudo, n_total = n, y_grid = y_grid)
  m_hat_pseudo  <- m_grid_here[which.min(cv_pseudo_m)]
  t_k_pseudo    <- seq(0, 1, length.out = m_hat_pseudo + 1L)
  F_knots_pseud <- ipw_ecdf_at(Y, w_pseudo, t_k_pseudo, n_total = n)
  F_bern_pseudo <- bernstein_smooth(F_knots_pseud, y_grid, m_hat_pseudo)
  F_bern_pseudo <- pmin(pmax(F_bern_pseudo, 0), 1)
  
  # KDE (LSCV over h)
  h_grid_here   <- h_candidates_of_n(n)
  cv_pseudo_h   <- vapply(h_grid_here, cv_score_kde_pseudo, numeric(1),
                          Y = Y, w = w_pseudo, n_total = n, y_grid = y_grid)
  h_hat_pseudo  <- h_grid_here[which.min(cv_pseudo_h)]
  F_kde_pseudo  <- ipw_kde_cdf(Y, w_pseudo, y_grid, n_total = n, h = h_hat_pseudo)
  F_kde_pseudo  <- pmin(pmax(F_kde_pseudo, 0), 1)
  
  ## ========= FEASIBLE: compute the three CDFs =========
  # Unsmoothed
  F_uns_feas <- ipw_ecdf_at(Y, w_feas, y_grid, n_total = n)
  F_uns_feas <- pmin(pmax(F_uns_feas, 0), 1)
  
  # Bernstein (LSCV over m)
  cv_feas_m   <- vapply(m_grid_here, cv_score_feasible, numeric(1),
                        Y = Y, w = w_feas, n_total = n, y_grid = y_grid)
  m_hat_feas  <- m_grid_here[which.min(cv_feas_m)]
  t_k_feas    <- seq(0, 1, length.out = m_hat_feas + 1L)
  F_knots_feas<- ipw_ecdf_at(Y, w_feas, t_k_feas, n_total = n)
  F_bern_feas <- bernstein_smooth(F_knots_feas, y_grid, m_hat_feas)
  F_bern_feas <- pmin(pmax(F_bern_feas, 0), 1)
  
  # KDE (LSCV over h)
  cv_feas_h   <- vapply(h_grid_here, cv_score_kde_feasible, numeric(1),
                        Y = Y, w = w_feas, n_total = n, y_grid = y_grid)
  h_hat_feas  <- h_grid_here[which.min(cv_feas_h)]
  F_kde_feas  <- ipw_kde_cdf(Y, w_feas, y_grid, n_total = n, h = h_hat_feas)
  F_kde_feas  <- pmin(pmax(F_kde_feas, 0), 1)
  
  ## ========= Plot: PSEUDO =========
  pdf(file.path(FIG_DIR, sprintf("CDF_pseudo_n%d.pdf", n)),
      width = 6.5, height = 4.5)
  plot(y_grid, F_uns_pseudo, type = "s", lwd = 1.2, col = "blue",
       xlab = "y", ylab = "Estimated CDF",
       main = sprintf("CDF (Pseudo weights), n = %d", n),
       ylim = c(0, 1), xaxs = "i", yaxs = "i")
  lines(y_grid, F_bern_pseudo, lwd = 1.2, col = "red", lty = 2)
  lines(y_grid, F_kde_pseudo,  lwd = 1.2, col = "darkgreen", lty = 3)
  legend("bottomright",
         legend = c(sprintf("Unsmoothed (IPW)"),
                    sprintf("Bernstein (m*=%d)", m_hat_pseudo),
                    sprintf("I-IPW KDE (h*=%.4f)", h_hat_pseudo)),
         col = c("blue", "red", "darkgreen"),
         lty = c(1, 2, 3), lwd = 1.2, bty = "n")
  dev.off()
  
  ## ========= Plot: FEASIBLE =========
  pdf(file.path(FIG_DIR, sprintf("CDF_feasible_n%d.pdf", n)),
      width = 6.5, height = 4.5)
  plot(y_grid, F_uns_feas, type = "s", lwd = 1.2, col = "blue",
       xlab = "y", ylab = "Estimated CDF",
       main = sprintf("CDF (Feasible weights), n = %d", n),
       ylim = c(0, 1), xaxs = "i", yaxs = "i")
  lines(y_grid, F_bern_feas, lwd = 1.2, col = "red", lty = 2)
  lines(y_grid, F_kde_feas,  lwd = 1.2, col = "darkgreen", lty = 3)
  legend("bottomright",
         legend = c(sprintf("Unsmoothed (IPW)"),
                    sprintf("Bernstein (m*=%d)", m_hat_feas),
                    sprintf("I-IPW KDE (h*=%.4f)", h_hat_feas)),
         col = c("blue", "red", "darkgreen"),
         lty = c(1, 2, 3), lwd = 1.2, bty = "n")
  dev.off()
}

## -----------------------
## Main Monte Carlo loops (parallel)
## -----------------------

t_global_start <- Sys.time()

y_grid      <- seq(0, 1, length.out = grid_points)
F_true_grid <- F_true(y_grid)

## Create cluster
cl <- parallel::makeCluster(n_cores)
on.exit(parallel::stopCluster(cl), add = TRUE)

## Reproducible RNG across workers
parallel::clusterSetRNGStream(cl, 20250913)

## Export necessary objects/functions to workers
parallel::clusterExport(
  cl,
  varlist = c(
    # hyperparameters / globals used inside workers
    "px", "y_dist", "pi_vals", "pi_min",
    "grid_points", "m_factor", "m_cap", "m_min",
    "kde_kernel",
    
    # grids precomputed on master
    "y_grid", "F_true_grid",
    
    # LSCV helpers and scoring functions
    "bernstein_term1_exact", "bernstein_tail_integrals",
    "cv_score_pseudo", "cv_score_feasible",
    
    # KDE CDF + LSCV bandwidth
    "get_kernel_cdf", "ipw_kde_cdf", "h_candidates_of_n",
    "cv_score_kde_core", "cv_score_kde_pseudo", "cv_score_kde_feasible",
    
    # functions used inside one_replication (and their deps)
    "one_replication", "simulate_dataset", "rX", "rY", "F_true", "pi_fun",
    "estimate_pi_hat", "ipw_ecdf_at", "bernstein_smooth", "m_candidates_of_n",
    "integral_mean", "compute_ISE"
  ),
  envir = environment()
)

## Containers for summary tables
pseudo_rows         <- list()  # (median/IQR — already present)
feasible_rows       <- list()

## NEW: mean/var containers
pseudo_rows_meanvar   <- list()
feasible_rows_meanvar <- list()

## Containers for boxplots (store full ISE distributions by n & method)
ise_pseudo_long   <- list()
ise_feasible_long <- list()

for (n in n_grid) {
  if (verbose) cat(sprintf("\nRunning simulations for n = %d ...\n", n))
  
  ## Run MC replications in parallel for this n
  res_list <- parallel::parLapplyLB(
    cl,
    X = seq_len(MC),
    fun = function(rep_idx, n_here, y_grid_arg, F_true_grid_arg) {
      one_replication(n = n_here, y_grid = y_grid_arg, F_true_grid = F_true_grid_arg)
    },
    n_here = n, y_grid_arg = y_grid, F_true_grid_arg = F_true_grid
  )
  
  ## Extract vectors
  ISE_uns_pseudo_vec <- vapply(res_list, function(z) z$ISE_uns_pseudo,  numeric(1))
  ISE_sm_pseudo_vec  <- vapply(res_list, function(z) z$ISE_sm_pseudo,   numeric(1))
  ISE_kde_pseudo_vec <- vapply(res_list, function(z) z$ISE_kde_pseudo,  numeric(1))
  
  ISE_uns_feas_vec   <- vapply(res_list, function(z) z$ISE_uns_feas,    numeric(1))
  ISE_sm_feas_vec    <- vapply(res_list, function(z) z$ISE_sm_feas,     numeric(1))
  ISE_kde_feas_vec   <- vapply(res_list, function(z) z$ISE_kde_feas,    numeric(1))
  
  ## Accumulate per-rep ISEs for boxplots (long format)
  ise_pseudo_long[[length(ise_pseudo_long) + 1L]] <- rbind(
    data.frame(n = n, estimator = "Unsmoothed",   ISE = ISE_uns_pseudo_vec),
    data.frame(n = n, estimator = "Bernstein",    ISE = ISE_sm_pseudo_vec),
    data.frame(n = n, estimator = "I-IPW KDE",    ISE = ISE_kde_pseudo_vec)
  )
  ise_feasible_long[[length(ise_feasible_long) + 1L]] <- rbind(
    data.frame(n = n, estimator = "Unsmoothed",   ISE = ISE_uns_feas_vec),
    data.frame(n = n, estimator = "Bernstein",    ISE = ISE_sm_feas_vec),
    data.frame(n = n, estimator = "I-IPW KDE",    ISE = ISE_kde_feas_vec)
  )
  
  ## Row for PSEUDO table (per n), metrics with (Unsmoothed, Bernstein, KDE)
  pseudo_rows[[length(pseudo_rows) + 1L]] <- data.frame(
    n = n,
    
    median_ISE.Unsmoothed = median(ISE_uns_pseudo_vec),
    median_ISE.Bernstein  = median(ISE_sm_pseudo_vec),
    median_ISE.KDE        = median(ISE_kde_pseudo_vec),
    
    IQR_ISE.Unsmoothed    = IQR(ISE_uns_pseudo_vec),
    IQR_ISE.Bernstein     = IQR(ISE_sm_pseudo_vec),
    IQR_ISE.KDE           = IQR(ISE_kde_pseudo_vec)
  )
  
  ## Row for FEASIBLE table (per n), metrics with (Unsmoothed, Bernstein, KDE)
  feasible_rows[[length(feasible_rows) + 1L]] <- data.frame(
    n = n,
    
    median_ISE.Unsmoothed = median(ISE_uns_feas_vec),
    median_ISE.Bernstein  = median(ISE_sm_feas_vec),
    median_ISE.KDE        = median(ISE_kde_feas_vec),
    
    IQR_ISE.Unsmoothed    = IQR(ISE_uns_feas_vec),
    IQR_ISE.Bernstein     = IQR(ISE_sm_feas_vec),
    IQR_ISE.KDE           = IQR(ISE_kde_feas_vec)
  )
  
  ## --- NEW: Row for PSEUDO (mean/var) ---
  pseudo_rows_meanvar[[length(pseudo_rows_meanvar) + 1L]] <- data.frame(
    n = n,
    mean_ISE.Unsmoothed = mean(ISE_uns_pseudo_vec),
    mean_ISE.Bernstein  = mean(ISE_sm_pseudo_vec),
    mean_ISE.KDE        = mean(ISE_kde_pseudo_vec),
    var_ISE.Unsmoothed  = var(ISE_uns_pseudo_vec),
    var_ISE.Bernstein   = var(ISE_sm_pseudo_vec),
    var_ISE.KDE         = var(ISE_kde_pseudo_vec)
  )
  
  ## --- NEW: Row for FEASIBLE (mean/var) ---
  feasible_rows_meanvar[[length(feasible_rows_meanvar) + 1L]] <- data.frame(
    n = n,
    mean_ISE.Unsmoothed = mean(ISE_uns_feas_vec),
    mean_ISE.Bernstein  = mean(ISE_sm_feas_vec),
    mean_ISE.KDE        = mean(ISE_kde_feas_vec),
    var_ISE.Unsmoothed  = var(ISE_uns_feas_vec),
    var_ISE.Bernstein   = var(ISE_sm_feas_vec),
    var_ISE.KDE         = var(ISE_kde_feas_vec)
  )
  
}

pseudo_table_grouped       <- do.call(rbind, pseudo_rows)            # median/IQR (already present)
feasible_table_grouped     <- do.call(rbind, feasible_rows)

## NEW:
pseudo_table_meanvar       <- do.call(rbind, pseudo_rows_meanvar)    # mean/var
feasible_table_meanvar     <- do.call(rbind, feasible_rows_meanvar)

## -----------------------
## Display the two grouped tables (console)
## -----------------------

cat("\n===== Pseudo-estimators (known propensity) — grouped columns =====\n")
print(pseudo_table_grouped, row.names = FALSE)

cat("\n===== Feasible estimators (estimated propensity) — grouped columns =====\n")
print(feasible_table_grouped, row.names = FALSE)

## -----------------------
## LaTeX writer: stacked sections (Median/IQR on top, Mean/Var below)
## -----------------------

## Fixed-width decimal formatter with truncation to N places
options(scipen = 999)  # avoid scientific notation globally
n_dig <- 7

.tex_fmt <- function(x, digits = n_dig, mode = c("truncate", "round")) {
  mode <- match.arg(mode)
  if (is.numeric(x)) {
    if (mode == "truncate") {
      # truncate toward 0 to exactly `digits` places
      scaled <- trunc(x * (10^digits)) / (10^digits)
      return(formatC(scaled, format = "f", digits = digits))
    } else {
      # standard rounding to `digits` places
      return(formatC(x, format = "f", digits = digits))
    }
  }
  x
}

write_grouped_table_tex_stacked <- function(df_top, df_bottom, file, caption, label,
                                            digits = 4, use_hlines = TRUE, tol = 1e-12) {
  # Expect column names:
  #   top:    n, median_ISE.{Unsmoothed,Bernstein,I-IPW KDE}, IQR_ISE.{Unsmoothed,Bernstein,I-IPW KDE}
  #   bottom: n, mean_ISE.{Unsmoothed,Bernstein,I-IPW KDE},   var_ISE.{Unsmoothed,Bernstein,I-IPW KDE}
  stopifnot(all(df_top$n == df_bottom$n), ncol(df_top) == 7, ncol(df_bottom) == 7)
  
  con <- file(file, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  
  writeLines("\\begin{table}[H]", con)
  writeLines("\\centering", con)
  writeLines("\\normalsize", con)
  
  colspec <- "|r|ccc|ccc|"
  writeLines(paste0("\\begin{tabular}{", colspec, "}"), con)
  if (use_hlines) writeLines("\\hline", con)
  
  ## ---------- TOP SECTION: Median / IQR ----------
  # Header row 1
  hdr1_top <- c(
    "$n$",
    "\\multicolumn{3}{|c|}{Median ISE}",
    "\\multicolumn{3}{|c|}{Interquartile Range ISE}"
  )
  writeLines(paste(hdr1_top, collapse = " & "), con); writeLines(" \\\\", con)
  if (use_hlines) writeLines("\\hline", con)
  
  # Header row 2 (methods)
  methods <- c("Unsmoothed","Bernstein","I-IPW KDE")
  hdr2 <- c(" ", methods, methods)
  writeLines(paste(hdr2, collapse = " & "), con); writeLines(" \\\\", con)
  if (use_hlines) writeLines("\\hline", con)
  
  # Identify top columns
  median_cols <- grep("^median_ISE\\.", names(df_top))
  iqr_cols    <- grep("^IQR_ISE\\.",    names(df_top))
  
  for (i in seq_len(nrow(df_top))) {
    row <- df_top[i, ]
    row_chars <- character(ncol(df_top))
    row_chars[1] <- as.character(as.integer(row[[1]]))
    row_chars[2:ncol(df_top)] <- vapply(row[2:ncol(df_top)], .tex_fmt, "", digits = digits)
    
    # Bold minima in each group
    vals_med <- as.numeric(row[median_cols]); if (all(is.finite(vals_med))) {
      mmin <- min(vals_med); idx <- median_cols[abs(vals_med - mmin) <= tol]
      for (j in idx) row_chars[j] <- paste0("\\textbf{", row_chars[j], "}")
    }
    vals_iqr <- as.numeric(row[iqr_cols]); if (all(is.finite(vals_iqr))) {
      mmin <- min(vals_iqr); idx <- iqr_cols[abs(vals_iqr - mmin) <= tol]
      for (j in idx) row_chars[j] <- paste0("\\textbf{", row_chars[j], "}")
    }
    
    writeLines(paste(paste(row_chars, collapse = " & "), "\\\\"),
               con)
  }
  
  if (use_hlines) writeLines("\\hline", con)
  if (use_hlines) writeLines("\\hline", con)
  
  ## ---------- BOTTOM SECTION: Mean / Var ----------
  # Reprint header row 1 (mean/var groups)
  hdr1_bot <- c(
    "$n$",
    "\\multicolumn{3}{|c|}{Mean ISE}",
    "\\multicolumn{3}{|c|}{Variance ISE ($\\times 10^3$)}"
  )
  writeLines(paste(hdr1_bot, collapse = " & "), con); writeLines(" \\\\", con)
  if (use_hlines) writeLines("\\hline", con)
  
  # Header row 2 (methods again)
  hdr2 <- c(" ", methods, methods)
  writeLines(paste(hdr2, collapse = " & "), con); writeLines(" \\\\", con)
  if (use_hlines) writeLines("\\hline", con)
  
  mean_cols <- grep("^mean_ISE\\.", names(df_bottom))
  var_cols  <- grep("^var_ISE\\.",  names(df_bottom))
  
  ## SCALE ONLY THE VARIANCE COLUMNS FOR THE LaTeX TABLE
  if (length(var_cols)) df_bottom[var_cols] <- lapply(df_bottom[var_cols], `*`, 1000)
  
  for (i in seq_len(nrow(df_bottom))) {
    row <- df_bottom[i, ]
    row_chars <- character(ncol(df_bottom))
    row_chars[1] <- as.character(as.integer(row[[1]]))
    row_chars[2:ncol(df_bottom)] <- vapply(row[2:ncol(df_bottom)], .tex_fmt, "", digits = digits)
    
    # Bold minima in each group
    vals_mean <- as.numeric(row[mean_cols]); if (all(is.finite(vals_mean))) {
      mmin <- min(vals_mean); idx <- mean_cols[abs(vals_mean - mmin) <= tol]
      for (j in idx) row_chars[j] <- paste0("\\textbf{", row_chars[j], "}")
    }
    vals_var <- as.numeric(row[var_cols]); if (all(is.finite(vals_var))) {
      mmin <- min(vals_var); idx <- var_cols[abs(vals_var - mmin) <= tol]
      for (j in idx) row_chars[j] <- paste0("\\textbf{", row_chars[j], "}")
    }
    
    writeLines(paste(paste(row_chars, collapse = " & "), "\\\\"),
               con)
  }
  
  if (use_hlines) writeLines("\\hline", con)
  writeLines("\\end{tabular}", con)
  writeLines(paste0("\\caption{", caption, "}"), con)
  writeLines(paste0("\\label{", label, "}"), con)
  writeLines("\\end{table}", con)
}

## -----------------------
## Save stacked LaTeX tables to .tex files
## -----------------------

write_grouped_table_tex_stacked(
  df_top    = pseudo_table_grouped,      # Median / IQR
  df_bottom = pseudo_table_meanvar,      # Mean / Var
  file      = file.path(TABLE_DIR, "table_pseudo.tex"),
  caption   = "Pseudo-estimators (known propensity). Top: Median ISE and Interquartile Range ISE. Bottom: Mean ISE and Variance ISE. Columns subdivided into Unsmoothed, Bernstein-smoothed, and Integrated-IPW KDE. Bold indicates the minimum for each $n$ within each group.",
  label     = "tab:pseudo.estimators",
  digits    = n_dig
)

write_grouped_table_tex_stacked(
  df_top    = feasible_table_grouped,    # Median / IQR
  df_bottom = feasible_table_meanvar,    # Mean / Var
  file      = file.path(TABLE_DIR, "table_feasible.tex"),
  caption   = "Feasible estimators (estimated propensity). Top: Median ISE and Interquartile Range ISE. Bottom: Mean ISE and Variance ISE. Columns subdivided into Unsmoothed, Bernstein-smoothed, and Integrated-IPW KDE. Bold indicates the minimum for each $n$ within each group.",
  label     = "tab:feasible.estimators",
  digits    = n_dig
)

cat("\nLaTeX tables saved as:\n  - table_pseudo.tex\n  - table_feasible.tex\n")

## -----------------------
## Boxplots of ISE by estimator and n
## -----------------------

if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(ggplot2)

df_pseudo_all   <- do.call(rbind, ise_pseudo_long)
df_feasible_all <- do.call(rbind, ise_feasible_long)

## Order the methods consistently
method_levels <- c("Unsmoothed", "Bernstein", "I-IPW KDE")
df_pseudo_all$estimator   <- factor(df_pseudo_all$estimator,   levels = method_levels)
df_feasible_all$estimator <- factor(df_feasible_all$estimator, levels = method_levels)

## (Optional) if you want log-scale on y because of outliers, set use_log <- TRUE
use_log <- FALSE

# Pseudo plot
p_pseudo <- ggplot(df_pseudo_all, aes(x = estimator, y = ISE)) +
  geom_boxplot(outlier.size = 0.6) +
  facet_wrap(
    ~ n,
    scales   = "free_y",
    labeller = labeller(n = function(v) paste0("n = ", v))
  ) +
  labs(title = "Pseudo-estimators (known propensity): ISE by method and n",
       x = NULL, y = "ISE") +
  theme_minimal(base_size = 12) +
  theme(
    strip.placement   = "outside",
    panel.spacing     = unit(6, "pt"),
    axis.text.x       = element_text(angle = 0, hjust = 0.5, size = 7)  # <-- smaller labels
  )

if (use_log) p_pseudo <- p_pseudo + scale_y_log10()

ggsave(filename = file.path(FIG_DIR, "boxplot_ISE_pseudo.pdf"),
       plot = p_pseudo, width = 7.5, height = 4.8, device = cairo_pdf)

# Feasible plot
p_feasible <- ggplot(df_feasible_all, aes(x = estimator, y = ISE)) +
  geom_boxplot(outlier.size = 0.6) +
  facet_wrap(
    ~ n,
    scales   = "free_y",
    labeller = labeller(n = function(v) paste0("n = ", v))
  ) +
  labs(title = "Feasible estimators (estimated propensity): ISE by method and n",
       x = NULL, y = "ISE") +
  theme_minimal(base_size = 12) +
  theme(
    strip.placement   = "outside",
    panel.spacing     = unit(6, "pt"),
    axis.text.x       = element_text(angle = 0, hjust = 0.5, size = 7)  # <-- smaller labels
  )

if (use_log) p_feasible <- p_feasible + scale_y_log10()

ggsave(filename = file.path(FIG_DIR, "boxplot_ISE_feasible.pdf"),
       plot = p_feasible, width = 7.5, height = 4.8, device = cairo_pdf)

cat("\nBoxplots saved as:\n  - boxplot_ISE_pseudo.pdf\n  - boxplot_ISE_feasible.pdf\n")

## -----------------------
## End timer
## -----------------------

t_global_end <- Sys.time()
cat(sprintf("\n[Timing] Monte Carlo section ran in %.2f minutes.\n",
            as.numeric(difftime(t_global_end, t_global_start, units = "mins"))))

