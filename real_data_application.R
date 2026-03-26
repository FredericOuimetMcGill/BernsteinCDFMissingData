## -----------------------
## ANALYSIS SAMPLE: fasting-lab SUBSAMPLE ONLY
## GLU_J defines eligibility for the fasting-lab analysis sample.
## Units outside this file are NOT treated as missing outcomes here.
## Missingness is defined only within this subsample:
##   delta = 1 if LBXGLU is observed, 0 otherwise.
## -----------------------

# Keep only DEMO rows whose SEQN appears in GLU (fasting subsample), then attach LBXGLU
dat_sub <- DEMO_small %>%
  dplyr::semi_join(GLU_small, by = "SEQN") %>%
  dplyr::left_join(GLU_small, by = "SEQN")

# Missingness indicator within the fasting-lab subsample
dat_sub <- dat_sub %>%
  dplyr::mutate(delta = as.integer(!is.na(LBXGLU)))

# Robust recoding helpers: nhanesA may return coded values either as 1/2 or as labels
std_exam_period <- function(z) {
  s <- trimws(tolower(as.character(z)))
  dplyr::case_when(
    s %in% c("1", "1.0") | grepl("nov", s) ~ "November--April",
    s %in% c("2", "2.0") | grepl("may", s) | grepl("oct", s) ~ "May--October",
    TRUE ~ NA_character_
  )
}

std_sex <- function(z) {
  s <- trimws(tolower(as.character(z)))
  dplyr::case_when(
    s %in% c("1", "1.0", "male") ~ "Male",
    s %in% c("2", "2.0", "female") ~ "Female",
    TRUE ~ NA_character_
  )
}

x_levels <- c(
  "November--April, Male",
  "November--April, Female",
  "May--October, Male",
  "May--October, Female"
)

# Discrete auxiliary X = (exam period) x (sex)
dat_sub <- dat_sub %>%
  dplyr::mutate(
    exam_period = std_exam_period(RIDEXMON),
    sex         = std_sex(RIAGENDR),
    X_desc = dplyr::case_when(
      exam_period == "November--April" & sex == "Male"   ~ x_levels[1],
      exam_period == "November--April" & sex == "Female" ~ x_levels[2],
      exam_period == "May--October"    & sex == "Male"   ~ x_levels[3],
      exam_period == "May--October"    & sex == "Female" ~ x_levels[4],
      TRUE ~ NA_character_
    ),
    X_desc = factor(X_desc, levels = x_levels),
    X      = factor(X_desc, levels = x_levels, labels = c("1", "2", "3", "4"))
  )

# Analysis frame: subsample units with X observed
dat_use <- dat_sub %>% dplyr::filter(!is.na(X))

cat("\nCounts of X cells (including NA, before filtering):\n")
print(table(dat_sub$X_desc, useNA = "ifany"))

if (nrow(dat_use) == 0L) {
  stop("No observations left after constructing X. Check how RIAGENDR and RIDEXMON were imported on this machine.")
}

# Quick counts so it is explicit we are in the subsample
n_total <- nrow(dat_use)
n_obs   <- sum(dat_use$delta)
n_miss  <- n_total - n_obs
cat("Analysis sample: fasting-lab subsample only\n")
cat("Missingness mechanism analyzed here: LBXGLU item nonresponse within the subsample\n")
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

# Raw glucose G = LBXGLU; transformed outcome Y01 lives in [0,1]
dat_use <- dat_use %>%
  dplyr::mutate(
    Y_raw = LBXGLU,
    Y01   = ifelse(is.na(Y_raw), NA_real_, rescale01(Y_raw, a, b))
  )

## -----------------------
## Estimate pi-hat(X) from discrete X (within subsample)
## Working MAR model in this application:
##   P(delta = 1 | Y_raw, X) = P(delta = 1 | X) = pi(X),
## with X in {1,2,3,4}. Thus pi(X) is constant within each cell.
## We estimate pi(X = x) by the observed fraction in cell x.
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

# Cell-level table: counts and estimated observation probabilities pi-hat(x) within the subsample
cell_tab <- data.frame(
  X = factor(c("1", "2", "3", "4"), levels = c("1", "2", "3", "4")),
  X_desc = factor(x_levels, levels = x_levels)
) %>%
  dplyr::left_join(
    dat_use %>%
      dplyr::group_by(X, X_desc) %>%
      dplyr::summarise(
        n        = dplyr::n(),
        observed = sum(delta),
        missing  = dplyr::n() - sum(delta),
        obs_rate = mean(delta),
        pi_hat   = first(pi_hat),
        .groups  = "drop"
      ),
    by = c("X", "X_desc")
  ) %>%
  dplyr::mutate(
    n        = dplyr::coalesce(n, 0L),
    observed = dplyr::coalesce(observed, 0L),
    missing  = dplyr::coalesce(missing, 0L),
    obs_rate = ifelse(n > 0, obs_rate, NA_real_),
    pi_hat   = ifelse(n > 0, pi_hat, NA_real_),
    X_num    = as.integer(X),
    observed_rate_percent = 100 * obs_rate
  ) %>%
  dplyr::arrange(X_num)

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
                  y_grid[bad + 1], diff(F_uns_grid)[bad]))
}
if (any(diff(F_sm_grid) < -eps)) {
  bad <- which(diff(F_sm_grid) < -eps)[1]
  warning(sprintf("Smoothed CDF decreased at y ≈ %.6f: Δ=%.3e",
                  y_grid[bad + 1], diff(F_sm_grid)[bad]))
}

cat(sprintf("Tail values (unsmoothed): F(0)=%.6f, F(1)=%.6f\n",
            F_uns_grid[1], F_uns_grid[length(F_uns_grid)]))
cat(sprintf("Tail values (smoothed)  : F(0)=%.6f, F(1)=%.6f\n",
            F_sm_grid[1],  F_sm_grid[length(F_sm_grid)]))

# Keep only finite values to avoid stray segments
keep <- is.finite(y_grid) & is.finite(F_uns_grid) & is.finite(F_sm_grid)
y_plot    <- y_grid[keep]
Funs_plot <- F_uns_grid[keep]
Fsm_plot  <- F_sm_grid[keep]

# Rebuild plotting frame in strictly increasing x
ord_plot <- order(y_plot)
df_plot <- data.frame(
  y = y_plot[ord_plot],
  glucose_mgdl = a + (b - a) * y_plot[ord_plot],
  Feasible_unsmoothed = Funs_plot[ord_plot],
  Feasible_smoothed   = Fsm_plot[ord_plot]
)

## -----------------------
## Plot and save: zoomed view on transformed scale
## -----------------------

# Narrow x-range on transformed scale (for a better view)
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
    y = "Estimated CDF"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position      = c(0.88, 0.50),
    legend.justification = c(1, 0.5),
    legend.key.width     = grid::unit(2.4, "cm")
  )

dev.off()

## -----------------------
## Full-range figure on the original mg/dL scale
## (equivalently, this spans the full transformed range [0,1])
## -----------------------

full_x_label <- if (y_rescale_mode == "fixed") {
  "Fasting plasma glucose (mg/dL)"
} else {
  "Fasting plasma glucose (capped scale, mg/dL)"
}

pdf(file.path(out_dir, "Feasible_CDFs_NHANES_original_scale_full.pdf"), width = 6.5, height = 4.5)

ggplot(df_plot, aes(x = glucose_mgdl)) +
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
  coord_cartesian(xlim = c(a, b), ylim = c(0, 1)) +
  scale_x_continuous(breaks = pretty(c(a, b), n = 6),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25),
                     expand = expansion(mult = c(0, 0))) +
  labs(
    x = full_x_label,
    y = "Estimated CDF"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position      = c(0.88, 0.50),
    legend.justification = c(1, 0.5),
    legend.key.width     = grid::unit(2.4, "cm")
  )

dev.off()

## -----------------------
## Build summary tables and write to LaTeX in ./simulation_tables/
## -----------------------

summary_df <- data.frame(
  n = n_total,
  observed_rate_percent = round(100 * mean(dat_use$delta), 1),
  m_star = m_hat
)

print(summary_df)

cell_table_df <- cell_tab %>%
  dplyr::transmute(
    X = X_num,
    Cell = as.character(X_desc),
    n = n,
    pi_hat = pi_hat
  )

print(cell_table_df)

# Simple LaTeX writer for the one-row summary table (no timing columns)
write_realdata_table_tex <- function(df, outfile, caption, label, use_hlines = TRUE) {
  con <- base::file(outfile, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  
  writeLines("\\begin{table}[!ht]", con)
  writeLines("\\centering", con)
  
  writeLines("\\begin{tabular}{c c c}", con)
  if (use_hlines) writeLines("\\hline", con)
  
  hdr <- c("$n$", "Observed rate (\\%)", "$m^*$")
  writeLines(paste(hdr, collapse = " & "), con)
  writeLines(" \\\\", con)
  if (use_hlines) writeLines("\\hline", con)
  
  row <- df[1, ]
  vals <- c(
    as.character(as.integer(row$n)),
    formatC(row$observed_rate_percent, digits = 1, format = "f"),
    as.character(as.integer(row$m_star))
  )
  writeLines(paste(paste(vals, collapse = " & "), "\\\\"), con)
  
  if (use_hlines) writeLines("\\hline", con)
  writeLines("\\end{tabular}", con)
  writeLines(paste0("\\caption{", caption, "}"), con)
  writeLines(paste0("\\label{", label, "}"), con)
  writeLines("\\end{table}", con)
}

# LaTeX writer for the four-cell propensity table
write_realdata_cell_table_tex <- function(df, outfile, caption, label, use_hlines = TRUE) {
  con <- base::file(outfile, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  
  writeLines("\\begin{table}[!ht]", con)
  writeLines("\\centering", con)
  
  writeLines("\\begin{tabular}{c l c c}", con)
  if (use_hlines) writeLines("\\hline", con)
  
  hdr <- c("$X$", "Cell", "$n_x$", "$\\widehat{\\pi}(x)$")
  writeLines(paste(hdr, collapse = " & "), con)
  writeLines(" \\\\", con)
  if (use_hlines) writeLines("\\hline", con)
  
  for (i in seq_len(nrow(df))) {
    vals <- c(
      as.character(as.integer(df$X[i])),
      df$Cell[i],
      as.character(as.integer(df$n[i])),
      formatC(df$pi_hat[i], digits = 3, format = "f")
    )
    writeLines(paste(paste(vals, collapse = " & "), "\\\\"), con)
  }
  
  if (use_hlines) writeLines("\\hline", con)
  writeLines("\\end{tabular}", con)
  writeLines(paste0("\\caption{", caption, "}"), con)
  writeLines(paste0("\\label{", label, "}"), con)
  writeLines("\\end{table}", con)
}

write_realdata_table_tex(
  df      = summary_df,
  outfile = file.path(output_dir, "table_realdata_feasible.tex"),
  caption = "NHANES 2017--2018 real-data application (feasible estimator with estimated propensity) on the fasting-lab subsample. The table reports the subsample size used (with $X$ observed), overall observed rate, and the LSCV-chosen Bernstein degree $m^*$.",
  label   = "tab:realdata.feasible"
)

write_realdata_cell_table_tex(
  df      = cell_table_df,
  outfile = file.path(output_dir, "table_realdata_cell_propensities.tex"),
  caption = "NHANES 2017--2018 fasting-lab subsample: the four sex-by-exam-period cells used for the discrete auxiliary variable $X$, together with their cell sizes and estimated observation probabilities $\\widehat{\\pi}(x)$.",
  label   = "tab:realdata.cell.propensities"
)

cat(sprintf("\nSaved files:\n  - %s\n  - %s\n  - %s\n  - %s\n  - %s\n",
            file.path(out_dir, "LSCV_feasible_curve.pdf"),
            file.path(out_dir, "Feasible_CDFs_NHANES.pdf"),
            file.path(out_dir, "Feasible_CDFs_NHANES_original_scale_full.pdf"),
            file.path(output_dir, "table_realdata_feasible.tex"),
            file.path(output_dir, "table_realdata_cell_propensities.tex")))

