## ============================================
## Real-data application: NHANES 2017-2018
## Feasible CDF (unsmoothed vs. Bernstein-smoothed with LSCV)
## Outcome: fasting plasma glucose (LBXGLU, mg/dL)
## Discrete auxiliary X: sex x exam period (RIAGENDR x RIDEXMON)
## Analysis sample = fasting-lab SUBSAMPLE ONLY
## ============================================

## -----------------------
## Packages
## -----------------------

pkgs <- c("nhanesA", "dplyr", "haven", "ggplot2")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(nhanesA)
library(dplyr)
library(ggplot2)

## -----------------------
## Hyperparameters (edit me)
## -----------------------

set.seed(2025)

# Transform Y (mg/dL) to [0,1]
# Option "quantile": use caps at empirical [q_low, q_high] of observed Y; clip outside.
# Option "fixed": use fixed clinical caps [a,b]; clip outside.
y_rescale_mode <- "fixed"          # "quantile" or "fixed"
y_quantile_caps <- c(0.01, 0.99)   # used if mode == "quantile"
y_fixed_caps    <- c(40, 460)      # mg/dL caps if mode == "fixed"

# LSCV grid & numerics
grid_points <- 501
m_min       <- 1L
m_cap       <- 1000L
m_factor    <- 3

# pi-hat floor to avoid division by ~0
pi_min_floor <- 1e-4

# Output directories
out_dir <- "nhanes_2017_2018_realdata"
if (!dir.exists(out_dir)) dir.create(out_dir, showWarnings = FALSE)

## -----------------------
## Set the path, then define/create output directory
## -----------------------

path <- "C://Users//fred1//Dropbox//Gharbi_Jedidi_Khardani_Ouimet_2025_missing_data"
setwd(path)

# Define output directories relative to 'path'
out_dir    <- file.path(getwd(), "nhanes_2017_2018_realdata")
output_dir <- file.path(getwd(), "simulation_tables")

# Create them (recursive=TRUE handles nested paths)
dir.create(out_dir,    showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Optional: sanity checks and write permission check
stopifnot(dir.exists(out_dir), dir.exists(output_dir))
if (file.access(out_dir, 2) != 0) stop("No write permission for out_dir: ", out_dir)
if (file.access(output_dir, 2) != 0) stop("No write permission for output_dir: ", output_dir)

cat("Outputs will be written to:\n  - ", normalizePath(out_dir), "\n  - ", normalizePath(output_dir), "\n", sep = "")

## -----------------------
## Download/load NHANES tables (via nhanesA)
## -----------------------

DEMO <- nhanesA::nhanes("DEMO_J")  # demographics 2017-2018
GLU  <- nhanesA::nhanes("GLU_J")   # fasting plasma glucose 2017-2018

# Keep the variables we need
DEMO_small <- DEMO %>%
  dplyr::select(SEQN, RIAGENDR, RIDAGEYR, RIDEXMON)  # sex, age, exam period
GLU_small  <- GLU %>%
  dplyr::select(SEQN, LBXGLU)                        # fasting plasma glucose (mg/dL)

## -----------------------
## ANALYSIS SAMPLE: fasting-lab SUBSAMPLE ONLY
## (units present in GLU_J define the subsample; LBXGLU may still be NA)
## -----------------------

# Keep only DEMO rows whose SEQN appears in GLU (fasting subsample), then attach LBXGLU
dat_sub <- DEMO_small %>%
  dplyr::semi_join(GLU_small, by = "SEQN") %>%
  dplyr::left_join(GLU_small, by = "SEQN")

# Missingness indicator within the subsample: 1 if glucose observed, 0 otherwise
dat_sub <- dat_sub %>%
  dplyr::mutate(delta = as.integer(!is.na(LBXGLU)))

# Discrete auxiliary X = (exam period) x (sex)
# RIAGENDR: 1=Male, 2=Female; RIDEXMON: 1=Nov–Apr, 2=May–Oct
dat_sub <- dat_sub %>%
  dplyr::mutate(
    X = interaction(RIDEXMON, RIAGENDR, drop = TRUE),
    X = factor(X)
  )

# Analysis frame: subsample units with X observed
dat_use <- dat_sub %>% dplyr::filter(!is.na(X))

# Quick counts so it is explicit we are in the subsample
n_total <- nrow(dat_use)
n_obs   <- sum(dat_use$delta)
n_miss  <- n_total - n_obs
cat("Analysis sample: fasting-lab subsample only\n")
cat(sprintf("n (subsample, X observed) = %d\n", n_total))
cat(sprintf("LBXGLU observed = %d | missing = %d (%.1f%% missing)\n",
            n_obs, n_miss, 100 * n_miss / n_total))

## -----------------------
## Rescale glucose Y to [0,1]
## -----------------------

Y_raw_obs <- dat_use$LBXGLU[dat_use$delta == 1]
if (y_rescale_mode == "quantile") {
  qs <- quantile(Y_raw_obs, probs = y_quantile_caps, na.rm = TRUE)
  a <- as.numeric(qs[1]); b <- as.numeric(qs[2])
} else {
  a <- y_fixed_caps[1]; b <- y_fixed_caps[2]
}
rescale01 <- function(y, a, b) pmin(pmax((y - a) / (b - a), 0), 1)

dat_use <- dat_use %>%
  dplyr::mutate(
    Y_raw = LBXGLU,
    Y01   = ifelse(is.na(Y_raw), NA_real_, rescale01(Y_raw, a, b))
  )

## -----------------------
## Estimate pi-hat(X) from discrete X (within subsample)
## -----------------------

estimate_pi_hat <- function(X, delta, pi_min = 1e-4) {
  Xf <- factor(X)  # ensure factor indexing
  cell_rate <- tapply(delta, Xf, function(d) max(mean(d), pi_min))
  as.numeric(cell_rate[as.character(Xf)])
}

dat_use <- dat_use %>%
  dplyr::mutate(
    pi_hat = estimate_pi_hat(X, delta, pi_min = pi_min_floor),
    w_feas = delta / pi_hat
  )

## -----------------------
## Core estimation helpers (unsmoothed, Bernstein smoothing, LSCV)
## -----------------------

ipw_ecdf_at <- function(Y, w, t_points, n_total) {
  obs <- which(w > 0 & !is.na(Y))
  if (!length(obs)) return(rep(0, length(t_points)))
  Yobs <- Y[obs]; wobs <- w[obs]
  ord  <- order(Yobs); Yso <- Yobs[ord]; wso <- wobs[ord]
  cs   <- cumsum(wso)
  idx  <- findInterval(t_points, Yso, rightmost.closed = TRUE, left.open = FALSE)
  Fv <- numeric(length(idx))
  pos <- idx > 0L
  Fv[pos] <- cs[idx[pos]]
  (1 / n_total) * Fv
}

bernstein_smooth <- function(F_at_knots, y_points, m) {
  ks <- 0:m
  vapply(y_points, function(y) sum(F_at_knots * dbinom(ks, size = m, prob = y)), numeric(1))
}

## --- Exact Term 1 via Beta integrals: ∫ B_m(F)^2 dy ---
bernstein_term1_exact <- function(F_knots, m) {
  ks <- 0:m
  Ikl <- outer(
    ks, ks,
    function(k, l) choose(m, k) * choose(m, l) * beta(k + l + 1, 2 * m - k - l + 1)
  )
  sum(outer(F_knots, F_knots) * Ikl)
}

## --- Tail integral of the basis: ∫_{y0}^1 b_{m,k}(y) dy ---
bernstein_tail_integrals <- function(y0, m) {
  ks <- 0:m
  (1 - pbeta(y0, ks + 1, m - ks + 1)) / (m + 1)
}

## --- LSCV for feasible estimator under Lebesgue ISE ---
cv_score_feasible <- function(m, Y, w, n_total) {
  t_k    <- seq(0, 1, length.out = m + 1L)
  Fk     <- ipw_ecdf_at(Y, w, t_k, n_total = n_total)
  
  # Term 1: exact ∫ B_m(F)^2 dy
  term1  <- bernstein_term1_exact(Fk, m)
  
  # Term 2: (2/n) * sum_i w_i * ∫_{Y_i}^1 B_m(F_{-i})(y) dy
  obs <- which(w > 0 & !is.na(Y))
  if (!length(obs)) return(term1)
  
  sumW <- n_total * Fk
  acc  <- 0
  for (i in obs) {
    ind_i    <- as.numeric(t_k >= Y[i])                  # 1{Y_i ≤ k/m}
    F_minus  <- (sumW - w[i] * ind_i) / (n_total - 1)
    tail_i   <- bernstein_tail_integrals(Y[i], m)
    acc      <- acc + w[i] * sum(F_minus * tail_i)
  }
  term2 <- (2 / n_total) * acc
  
  term1 - term2
}

## Build m-grid
m_candidates_of_n <- function(n) {
  m_max <- min(floor(m_factor * n^(2/3)), m_cap, n)
  if (m_max < m_min) m_max <- m_min
  seq.int(m_min, m_max, by = 1L)
}

## -----------------------
## Fit feasible CDFs (unsmoothed vs. smoothed via LSCV)
## -----------------------

y_grid  <- seq(0, 1, length.out = grid_points)

# Unsmoothed feasible + timing
t0 <- proc.time()[3]
F_uns_grid <- ipw_ecdf_at(Y = dat_use$Y01, w = dat_use$w_feas, t_points = y_grid, n_total = n_total)
time_uns_ms <- 1000 * (proc.time()[3] - t0)

# Smoothed feasible via LSCV over data-driven m-grid + timing
t0 <- proc.time()[3]
m_grid <- m_candidates_of_n(n_total)
cv_vals <- vapply(m_grid, cv_score_feasible, numeric(1), Y = dat_use$Y01, w = dat_use$w_feas, n_total = n_total)
m_hat   <- m_grid[which.min(cv_vals)]
t_k     <- seq(0, 1, length.out = m_hat + 1L)
F_knots <- ipw_ecdf_at(dat_use$Y01, dat_use$w_feas, t_points = t_k, n_total = n_total)
F_sm_grid <- bernstein_smooth(F_knots, y_grid, m_hat)
time_sm_ms <- 1000 * (proc.time()[3] - t0)

## -----------------------
## Diagnostics & outputs
## -----------------------

cat(sprintf("Sample size (subsample, X observed): %d\n", n_total))
cat(sprintf("Overall observed rate (delta=1): %.1f%%\n", 100 * mean(dat_use$delta)))
cat(sprintf("Chosen Bernstein degree by LSCV: m* = %d\n", m_hat))
cat(sprintf("Rescaling caps for glucose (mg/dL): [a, b] = [%.1f, %.1f]\n", a, b))
cat(sprintf("Timing (ms): unsmoothed = %.1f, smoothed+LSCV = %.1f\n", time_uns_ms, time_sm_ms))

# Cell-level table: counts, observed rate, pi-hat (within subsample)
cell_tab <- dat_use %>%
  dplyr::group_by(X) %>%
  dplyr::summarise(
    n = dplyr::n(),
    obs_rate = mean(delta),
    pi_hat   = first(pi_hat),
    .groups = "drop"
  )
print(cell_tab)

## -----------------------
## LSCV curve on the real data
## -----------------------

pdf(file.path(out_dir, "LSCV_feasible_curve.pdf"), width = 6, height = 4)
plot(m_grid, cv_vals, type = "b", xlab = "m (Bernstein degree)", ylab = "LSCV(m)",
     main = "LSCV curve (feasible, NHANES 2017-2018)")
abline(v = m_hat, lty = 2)
legend("topright", legend = sprintf("m* = %d", m_hat), bty = "n")
dev.off()

## -----------------------
## Sanity checks before plotting
## -----------------------

stopifnot(length(y_grid) == length(F_uns_grid),
          length(y_grid) == length(F_sm_grid))

# Monotonicity and endpoint checks (allow tiny numerical jitter)
eps <- 1e-12
if (any(diff(F_uns_grid) < -eps)) {
  bad <- which(diff(F_uns_grid) < -eps)[1]
  warning(sprintf("Unsmoothed CDF decreased at y ≈ %.6f: Δ=%.3e",
                  y_grid[bad+1], diff(F_uns_grid)[bad]))
}
if (any(diff(F_sm_grid) < -eps)) {
  bad <- which(diff(F_sm_grid) < -eps)[1]
  warning(sprintf("Smoothed CDF decreased at y ≈ %.6f: Δ=%.3e",
                  y_grid[bad+1], diff(F_sm_grid)[bad]))
}

cat(sprintf("Tail values (unsmoothed): F(0)=%.6f, F(1)=%.6f\n",
            F_uns_grid[1], F_uns_grid[length(F_uns_grid)]))
cat(sprintf("Tail values (smoothed)  : F(0)=%.6f, F(1)=%.6f\n",
            F_sm_grid[1],  F_sm_grid[length(F_sm_grid)]))

# Keep only finite values to avoid stray segments
keep <- is.finite(y_grid) & is.finite(F_uns_grid) & is.finite(F_sm_grid)
y_plot   <- y_grid[keep]
Funs_plot<- F_uns_grid[keep]
Fsm_plot <- F_sm_grid[keep]

# Rebuild plotting frame in strictly increasing x
ord <- order(y_plot)
df_plot <- data.frame(
  y = y_plot[ord],
  Feasible_unsmoothed = Funs_plot[ord],
  Feasible_smoothed   = Fsm_plot[ord]
)

## -----------------------
## Plot and save: step for ECDF, line for smoothed
## -----------------------

# Narrow x-range (for a better view)
x_min <- 0.05
x_max <- 0.40

pdf(file.path(out_dir, "Feasible_CDFs_NHANES.pdf"), width = 6.5, height = 4.5)

ggplot(df_plot, aes(x = y)) +
  geom_step(
    aes(y = Feasible_unsmoothed,
        color = "Unsmoothed IPW",
        linetype = "Unsmoothed IPW"),
    direction = "hv",
    linewidth = 0.6,
    na.rm = TRUE
  ) +
  geom_line(
    aes(y = Feasible_smoothed,
        color = "Bernstein-smoothed",
        linetype = "Bernstein-smoothed"),
    linewidth = 0.6,
    na.rm = TRUE
  ) +
  scale_color_manual(
    name   = NULL,
    values = c("Unsmoothed IPW" = "blue",
               "Bernstein-smoothed" = "red")
  ) +
  scale_linetype_manual(
    name   = NULL,
    values = c("Unsmoothed IPW" = "solid",
               "Bernstein-smoothed" = "dashed")
  ) +
  guides(
    color = guide_legend(override.aes = list(linewidth = 1.0))
  ) +
  coord_cartesian(xlim = c(x_min, x_max), ylim = c(0, 1)) +
  scale_x_continuous(breaks = seq(x_min, x_max, by = 0.10),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25),
                     expand = expansion(mult = c(0, 0))) +
  labs(
    x = "Transformed glucose",
    y = "Estimated CDF",
    title = "         Feasible CDF estimates of fasting plasma glucose"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position      = c(0.88, 0.50),        # middle-right
    legend.justification = c(1, 0.5),
    legend.key.width     = grid::unit(2.4, "cm")
  )

dev.off()

## -----------------------
## Build a small summary table and write to LaTeX in ./simulation_tables/
## (no timing columns)
## -----------------------

summary_df <- data.frame(
  n = n_total,
  observed_rate_percent = round(100 * mean(dat_use$delta), 1),
  m_star = m_hat
)

print(summary_df)

# Simple LaTeX writer for this one-row summary table (no timing columns)
write_realdata_table_tex <- function(df, file, caption, label, use_hlines = TRUE) {
  con <- file(file, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  
  writeLines("\\begin{table}[!ht]", con)
  writeLines("\\centering", con)
  
  writeLines("\\begin{tabular}{c c c}", con)
  if (use_hlines) writeLines("\\hline", con)
  
  # Header
  hdr <- c("$n$", "Observed rate (\\%)", "$m^*$")
  writeLines(paste(hdr, collapse = " & "), con)
  writeLines(" \\\\", con)
  if (use_hlines) writeLines("\\hline", con)
  
  # Single row
  row <- df[1, ]
  vals <- c(
    as.character(as.integer(row$n)),
    formatC(row$observed_rate_percent, digits = 1, format = "f"),
    as.character(as.integer(row$m_star))
  )
  writeLines(paste(paste(vals, collapse = " & "), "\\\\"),
             con)
  
  if (use_hlines) writeLines("\\hline", con)
  writeLines("\\end{tabular}", con)
  writeLines(paste0("\\caption{", caption, "}"), con)
  writeLines(paste0("\\label{", label, "}"), con)
  writeLines("\\end{table}", con)
}

write_realdata_table_tex(
  df      = summary_df,
  file    = file.path(output_dir, "table_realdata_feasible.tex"),
  caption = "NHANES 2017--2018 real-data application (feasible estimator with estimated propensity) on the fasting-lab subsample. The table reports the subsample size used (with $X$ observed), overall observed rate, and the LSCV-chosen Bernstein degree $m^*$.",
  label   = "tab:realdata.feasible"
)

cat(sprintf("\nSaved files:\n  - %s\n  - %s\n  - %s\n",
            file.path(out_dir, "LSCV_feasible_curve.pdf"),
            file.path(out_dir, "Feasible_CDFs_NHANES.pdf"),
            file.path(output_dir, "table_realdata_feasible.tex")))
