# BFI: An R Package for Bayesian Federated Inference
# Reproducibility code for the Computer Methods and Programs in Biomedicine manuscript
# Verified by the author on 2026-11-10: 31/31 numerical checks passed.
# R 4.5.1; BFI 3.4.0 (GitHub); survival 3.8-12; Windows 11 (64-bit).
# Run in a clean R session: source(file.choose(), echo = FALSE)
# The script prints analyses and checks to the R console and displays the
# Weibull plot in the active graphics device. It creates no output files.
# All datasets are distributed with BFI or survival.

expected_bfi_version <- package_version("3.4.0")
expected_survival_version <- package_version("3.8-12")

if (!requireNamespace("BFI", quietly = TRUE)) {
  stop("BFI is not installed; install version 3.4.0 before running this script.")
}
if (!requireNamespace("survival", quietly = TRUE)) {
  stop("The survival R package is required.")
}
if (packageVersion("BFI") != expected_bfi_version) {
  stop("Expected BFI 3.4.0; installed version is ",
       as.character(packageVersion("BFI")), ".")
}
if (packageVersion("survival") != expected_survival_version) {
  warning("Manuscript used survival 3.8-12; installed version is ",
          as.character(packageVersion("survival")),
          ". Check the numerical comparison report carefully.")
}

library(BFI)
library(survival)

if (getRversion() != numeric_version("4.5.1")) {
  warning("The manuscript used R 4.5.1; current R version is ",
          as.character(getRversion()), ".")
}

cat("R:", as.character(getRversion()), "| BFI:", as.character(packageVersion("BFI")),
    "| survival:", as.character(packageVersion("survival")), "\n")

# Compare numerical quantities with values rounded for the manuscript.
# A failure is reported, not silently passed; it does not stop later analyses.
checks <- data.frame(section = character(), quantity = character(),
                     observed = numeric(), expected = numeric(),
                     tolerance = numeric(), passed = logical())
check_value <- function(section, quantity, observed, expected, tolerance = 0.001) {
  observed <- as.numeric(observed)
  if (length(observed) != 1L || is.na(observed) || !is.finite(observed)) {
    stop("Invalid observed value for ", section, ": ", quantity)
  }
  ok <- abs(observed - expected) <= tolerance
  checks <<- rbind(checks, data.frame(
    section = section, quantity = quantity, observed = observed,
    expected = expected, tolerance = tolerance, passed = ok
  ))
  cat(if (ok) "[PASS] " else "[MISMATCH] ", section, ": ", quantity,
      " | observed=", signif(observed, 8),
      " expected=", expected, "\n", sep = "")
  invisible(ok)
}






# ============================================================================
# Section 4. A complete BFI workflow: Nurses data
# ============================================================================


data("Nurses", package = "BFI")

local_data <- split(Nurses, Nurses$hospital)
centers <- sort(unique(Nurses$hospital))
L <- length(centers)
L

local_analysis <- function(dat) {
  X <- data.frame(
    age = dat$age,
    gender = factor(dat$gender, levels = c(0, 1)),
    experience = dat$experien,
    wardtype = factor(dat$wardtype, levels = c(0, 1))
  )

  Lambda <- inv.prior.cov(
    X = X,
    lambda = 0.01,
    family = "gaussian"
  )

  fit <- MAP.estimation(
    y = dat$stress,
    X = X,
    family = "gaussian",
    Lambda = Lambda,
    control = list(maxit = 1000)
  )

  stopifnot(fit$convergence == 0)

  list(
    theta_hat = fit$theta_hat,
    A_hat = fit$A_hat,
    Lambda = fit$Lambda
  )
}

local_results <- lapply(local_data, local_analysis)


# Formula interface for the first hospital
d1 <- local_data[[1]]
d1$gender <- factor(d1$gender, levels = c(0, 1))
d1$wardtype <- factor(d1$wardtype, levels = c(0, 1))
d1$experience <- d1$experien

form1 <- stress ~ age + gender + experience + wardtype

Lambda_formula <- inv.prior.cov(
  formula = form1,
  data = d1,
  lambda = 0.01,
  family = "gaussian"
)

fit_one <- MAP.estimation(
  formula = form1,
  data = d1,
  family = "gaussian",
  Lambda = Lambda_formula,
  control = list(maxit = 1000)
)

summary(fit_one)

# Equivalence with the original y/X interface
stopifnot(
  isTRUE(all.equal(local_results[[1]]$Lambda, Lambda_formula)),
  isTRUE(all.equal(local_results[[1]]$theta_hat, fit_one$theta_hat)),
  isTRUE(all.equal(local_results[[1]]$A_hat, fit_one$A_hat))
)

# Section 3.2 example: extended 3.4.0 formula interface, using local data.
# This demonstration is separate from the models used in the numerical tables.
d <- d1
form <- stress ~ age * factor(gender) + I(experience^2)
Lambda_extended <- inv.prior.cov(
  formula = form, data = d, lambda = 0.01, family = "gaussian"
)
fit_extended <- MAP.estimation(
  formula = form, data = d, family = "gaussian", Lambda = Lambda_extended,
  control = list(maxit = 1000)
)
stopifnot(fit_extended$convergence == 0)


theta_hats <- lapply(local_results, `[[`, "theta_hat")
A_hats <- lapply(local_results, `[[`, "A_hat")
Lambda <- local_results[[1]]$Lambda

fit_bfi <- bfi(
  theta_hats = theta_hats,
  A_hats = A_hats,
  Lambda = Lambda,
  family = "gaussian"
)

class(fit_bfi)
summary(fit_bfi)
summary(fit_bfi, cur_mat = TRUE)

# Pooled benchmark
X_pooled <- data.frame(
  age = Nurses$age,
  gender = factor(Nurses$gender, levels = c(0, 1)),
  experience = Nurses$experien,
  wardtype = factor(Nurses$wardtype, levels = c(0, 1))
)

Lambda_pooled <- inv.prior.cov(
  X = X_pooled,
  lambda = 0.01,
  family = "gaussian"
)

fit_pooled <- MAP.estimation(
  y = Nurses$stress,
  X = X_pooled,
  family = "gaussian",
  Lambda = Lambda_pooled,
  control = list(maxit = 1000)
)

stopifnot(fit_pooled$convergence == 0)
#summary(fit_pooled)

bfi_est <- coef(fit_bfi)
pooled_est <- coef(fit_pooled)

sd_bfi <- sqrt(diag(vcov(fit_bfi)))
sd_pooled <- sqrt(diag(vcov(fit_pooled)))

comparison <- data.frame(
  parameter = names(bfi_est),
  BFI = bfi_est,
  pooled = pooled_est,
  difference = bfi_est - pooled_est,
  BFI_sd = sd_bfi,
  pooled_sd = sd_pooled
)

print(comparison, digits = 3, row.names = FALSE)
max(abs(bfi_est - pooled_est))
max(abs(sd_bfi - sd_pooled))


check_value("4", "BFI intercept", coef(fit_bfi)["(Intercept)"], 5.6313)
check_value("4", "pooled intercept", coef(fit_pooled)["(Intercept)"], 5.4593)
check_value("4", "BFI residual variance", coef(fit_bfi)["sigma2"], 0.4352)
check_value("4", "pooled residual variance", coef(fit_pooled)["sigma2"], 0.8701)

# ============================================================================
# Section 5. Handling heterogeneity across centers
# ============================================================================

# Center-specific intercepts
Lambda_intercept <- inv.prior.cov(
  X = X_pooled,
  lambda = 0.01,
  family = "gaussian",
  stratified = TRUE,
  strat_par = 1,
  L = L
)

fit_intercept <- bfi(
  theta_hats = theta_hats,
  A_hats = A_hats,
  Lambda = list(Lambda, Lambda_intercept),
  family = "gaussian",
  stratified = TRUE,
  strat_par = 1
)
summary(fit_intercept)

est <- coef(fit_intercept)
ints <- est[grep("Intercept", names(est))]
round(range(ints), 4)

# Center-specific residual variances
Lambda_sigma <- inv.prior.cov(
  X = X_pooled,
  lambda = 0.01,
  family = "gaussian",
  stratified = TRUE,
  strat_par = 2,
  L = L
)

fit_sigma <- bfi(
  theta_hats = theta_hats,
  A_hats = A_hats,
  Lambda = list(Lambda, Lambda_sigma),
  family = "gaussian",
  stratified = TRUE,
  strat_par = 2
)
summary(fit_sigma)

est <- coef(fit_sigma)
s2 <- est[grep("sigma2_loc", names(est))]
round(range(s2), 4)
round(max(s2) / min(s2), 1)

# Center-specific intercepts and residual variances
Lambda_both <- inv.prior.cov(
  X = X_pooled,
  lambda = 0.01,
  family = "gaussian",
  stratified = TRUE,
  strat_par = c(1, 2),
  L = L
)

fit_both <- bfi(
  theta_hats = theta_hats,
  A_hats = A_hats,
  Lambda = list(Lambda, Lambda_both),
  family = "gaussian",
  stratified = TRUE,
  strat_par = c(1, 2)
)
summary(fit_both)

tb <- drop(fit_both$theta_hat)
max(abs(tb[grep("Intercept", names(tb))] - ints))
max(abs(tb[grep("sigma2_loc", names(tb))] - s2))


# Center-specific variable: hospital size
Hsize <- vapply(
  local_data,
  function(dat) dat$hospsize[1],
  numeric(1)
)
Hsize <- factor(Hsize)

Lambda_hsize <- inv.prior.cov(
  X = X_pooled,
  lambda = 0.01,
  family = "gaussian",
  stratified = TRUE,
  center_spec = Hsize,
  L = L
)

fit_hsize <- bfi(
  theta_hats = theta_hats,
  A_hats = A_hats,
  Lambda = list(Lambda, Lambda_hsize),
  family = "gaussian",
  stratified = TRUE,
  center_spec = Hsize
)
summary(fit_hsize)


hsize_intercepts <- coef(fit_hsize)[grep("Intercept", names(coef(fit_hsize)))]
print(hsize_intercepts)
if (length(hsize_intercepts) == 3L) {
  check_value("5", "small-hospital intercept", hsize_intercepts[1], 5.1141)
  check_value("5", "medium-hospital intercept", hsize_intercepts[2], 5.5953)
  check_value("5", "large-hospital intercept", hsize_intercepts[3], 6.0640)
} else {
  warning("Expected three hospital-size intercepts; found ",
          length(hsize_intercepts))
}

# ============================================================================
# Section 6. Survival models in practice
# ============================================================================

set.seed(112358)
p <- 3
beta <- 1:3
a <- 5
b <- 6

n1 <- 30
X1 <- data.frame(matrix(rnorm(n1 * p), n1, p))
y1 <- surv.simulate(
  Z = list(X1),
  beta = beta,
  a = a,
  b = b,
  u1 = 0.1,
  cen_rate = 0.3,
  gen_data_from = "weibul"
)$D[[1]][, 1:2]

n2 <- 30
X2 <- data.frame(matrix(rnorm(n2 * p), n2, p))
y2 <- surv.simulate(
  Z = list(X2),
  beta = beta,
  a = a,
  b = b,
  u1 = 0.1,
  cen_rate = 0.3,
  gen_data_from = "weibul"
)$D[[1]][, 1:2]

# Weibull baseline hazard
Lambda_weib <- inv.prior.cov(
  X = X1,
  lambda = c(0.1, 1),
  family = "survival",
  basehaz = "weibul"
)

fit_weib1 <- MAP.estimation(
  y = y1,
  X = X1,
  family = "survival",
  Lambda = Lambda_weib,
  basehaz = "weibul"
)

fit_weib2 <- MAP.estimation(
  y = y2,
  X = X2,
  family = "survival",
  Lambda = Lambda_weib,
  basehaz = "weibul"
)

fit_weib_bfi <- bfi(
  theta_hats = list(fit_weib1$theta_hat, fit_weib2$theta_hat),
  A_hats = list(fit_weib1$A_hat, fit_weib2$A_hat),
  Lambda = Lambda_weib,
  family = "survival",
  basehaz = "weibul"
)
summary(fit_weib_bfi)

# Exponentiated polynomial baseline hazard
alpha_poly <- 0.1
max_order <- 2

Lambda_poly <- inv.prior.cov(
  X = X1,
  lambda = c(0.1, 1),
  family = "survival",
  basehaz = "poly",
  max_order = max_order
)

fit_poly1 <- MAP.estimation(
  y = y1,
  X = X1,
  family = "survival",
  Lambda = Lambda_poly,
  basehaz = "poly",
  alpha = alpha_poly,
  max_order = max_order
)

fit_poly2 <- MAP.estimation(
  y = y2,
  X = X2,
  family = "survival",
  Lambda = Lambda_poly,
  basehaz = "poly",
  alpha = alpha_poly,
  max_order = max_order
)

theta_A_polys <- list(fit_poly1$theta_A_poly, fit_poly2$theta_A_poly)
q_ls <- c(fit_poly1$q_l, fit_poly2$q_l)
q_max <- max(q_ls)
print(q_ls)

fit_poly_bfi <- bfi(
  Lambda = Lambda_poly,
  family = "survival",
  basehaz = "poly",
  theta_A_polys = theta_A_polys,
  q_ls = q_ls
)
summary(fit_poly_bfi)

# Piecewise exponential baseline hazard
n_intervals <- 3
local_max_times <- c(max(y1$time), max(y2$time))
min_max_times <- min(local_max_times)
min_max_times

Lambda_pwexp <- inv.prior.cov(
  X = X1,
  lambda = c(0.1, 1),
  family = "survival",
  basehaz = "pwexp",
  n_intervals = n_intervals
)

fit_pwexp1 <- MAP.estimation(
  y = y1,
  X = X1,
  family = "survival",
  Lambda = Lambda_pwexp,
  basehaz = "pwexp",
  n_intervals = n_intervals,
  min_max_times = min_max_times
)

fit_pwexp2 <- MAP.estimation(
  y = y2,
  X = X2,
  family = "survival",
  Lambda = Lambda_pwexp,
  basehaz = "pwexp",
  n_intervals = n_intervals,
  min_max_times = min_max_times
)

fit_pwexp_bfi <- bfi(
  theta_hats = list(fit_pwexp1$theta_hat, fit_pwexp2$theta_hat),
  A_hats = list(fit_pwexp1$A_hat, fit_pwexp2$A_hat),
  Lambda = Lambda_pwexp,
  family = "survival",
  basehaz = "pwexp"
)
summary(fit_pwexp_bfi)

times <- seq(0.05, min_max_times, length.out = 100)
h_bfi <- hazards.fun(time = times, p = 3, theta_hat = fit_weib_bfi$theta_hat, basehaz = "weibul")
#plot(times, h_bfi$bhazard, type = "l", xlab = "Time", ylab = "Baseline hazard")


check_value("6", "Weibull BFI X1", as.numeric(fit_weib_bfi$theta_hat[1]), 1.3653)
check_value("6", "Weibull BFI X2", as.numeric(fit_weib_bfi$theta_hat[2]), 2.1354)
check_value("6", "Weibull BFI X3", as.numeric(fit_weib_bfi$theta_hat[3]), 3.5779)
check_value("6", "poly local order center 1", fit_poly1$q_l, 2, 0)
check_value("6", "poly local order center 2", fit_poly2$q_l, 2, 0)
check_value("6", "piecewise min_max_times", min_max_times, 1.286064, 0.00001)

# ============================================================================
# Section 7. Treatment-effect estimation
# ============================================================================

# ----------------------------------------------------------------------------
# Section 7.2. Survival treatment-effect example: Rotterdam data
# ----------------------------------------------------------------------------

# survival was loaded and checked at the start of the script.

data(cancer, package = "survival")
df <- rotterdam

ignore <- with(
  df,
  recur == 0 & death == 1 & rtime < dtime
)

df$status <- with(
  df,
  ifelse(recur == 1 | ignore, recur, death)
)

df$time <- with(
  df,
  ifelse(recur == 1 | ignore, rtime, dtime)
) / 365.25

df <- df[, c(
  "time", "status", "chemo", "year", "age", "meno",
  "size", "grade", "nodes", "pgr", "er", "hormon"
)]

df$size <- factor(df$size)
df$grade <- factor(df$grade)

df$age <- as.numeric(scale(df$age))
df$pgr <- as.numeric(scale(df$pgr))
df$er <- as.numeric(scale(df$er))
df$year <- as.numeric(scale(df$year))

# Main pseudo-center allocation: random allocation without treatment stratification
L <- 3
set.seed(11235)

groups <- cut(
  sample(nrow(df)),
  breaks = L,
  labels = FALSE
)

datasets <- split(df, groups)
print(table(groups, df$chemo))

X_list <- lapply(
  datasets,
  function(dat) {
    dat[, !(names(dat) %in% c("time", "status")), drop = FALSE]
  }
)

Y_list <- lapply(
  datasets,
  function(dat) {
    dat[, c("time", "status"), drop = FALSE]
  }
)

# First round: federated propensity-score model
Lambda_r1 <- lapply(
  X_list,
  function(X) {
    inv.prior.cov(
      X,
      lambda = 0.01,
      family = "survival",
      treatment = "chemo",
      treat_round = "first"
    )
  }
)

fit_r1 <- lapply(
  seq_len(L),
  function(l) {
    MAP.estimation(
      Y_list[[l]],
      X_list[[l]],
      family = "survival",
      Lambda = Lambda_r1[[l]],
      treatment = "chemo",
      treat_round = "first"
    )
  }
)

fitbfi_r1 <- bfi(
  theta_hats = lapply(fit_r1, function(x) x$theta_hat),
  A_hats = lapply(fit_r1, function(x) x$A_hat),
  Lambda = Lambda_r1[[1]],
  family = "survival",
  treat_round = "first"
)

# Second round: weighted Cox treatment analysis
Lambda_r2 <- lapply(
  X_list,
  function(X) {
    inv.prior.cov(
      X,
      lambda = 0.01,
      family = "survival",
      basehaz = "unspecified",
      treatment = "chemo",
      treat_round = "second"
    )
  }
)

fit_r2 <- lapply(
  seq_len(L),
  function(l) {
    MAP.estimation(
      Y_list[[l]],
      X_list[[l]],
      family = "survival",
      Lambda = Lambda_r2[[l]],
      basehaz = "unspecified",
      treatment = "chemo",
      treat_round = "second",
      gamma_bfi = fitbfi_r1$theta_hat
    )
  }
)

fitbfi_r2 <- bfi(
  theta_hats = lapply(fit_r2, function(x) x$theta_hat),
  A_hats = lapply(fit_r2, function(x) x$A_hat),
  Lambda = Lambda_r2[[1]],
  family = "survival",
  basehaz = "unspecified",
  treat_round = "second"
)

print(fitbfi_r2$theta_hat)
summary(fitbfi_r2)

check_value("7.2", "weighted Cox BFI coefficient",
            as.numeric(coef(fitbfi_r2)["chemo"]), -0.1984)

# Save main-allocation objects for Section 8.3 validation
rotterdam_main <- list(
  df = df,
  groups = groups,
  datasets = datasets,
  X_list = X_list,
  Y_list = Y_list,
  Lambda_r1 = Lambda_r1,
  fit_r1 = fit_r1,
  fitbfi_r1 = fitbfi_r1,
  Lambda_r2 = Lambda_r2,
  fit_r2 = fit_r2,
  fitbfi_r2 = fitbfi_r2
)


# ----------------------------------------------------------------------------
# Section 7.3. Average treatment-effect outputs for binomial models
# ----------------------------------------------------------------------------
# IMPORTANT: This block preserves the exact random-number order used for the
# reported results in the manuscript.

set.seed(112358)

beta_bin <- 1:3
n <- c(100, 200)
sim_center <- function(n_l) {X <- data.frame(
             x1 = rnorm(n_l),
             treatment = factor(sample(1:2, n_l, replace = TRUE))
         )
       eta <- beta_bin[1] + beta_bin[2] * X$x1 +
                beta_bin[3] * as.numeric(X$treatment == "2")
       list(
             X = X,
             y = rbinom(n_l, size = 1, prob = binomial()$linkinv(eta))
         )
 }
centers_bin <- lapply(n, sim_center)

local_fit <- function(d, round, ...) {
  Lambda <- inv.prior.cov(
             d$X,
             lambda = 0.01,
             family = "binomial",
             treatment = "treatment",
             treat_round = round
         )
       MAP.estimation(
             d$y,
             d$X,
             family = "binomial",
             Lambda = Lambda,
             treatment = "treatment",
             treat_round = round,
             ...
         )
 }
central_fit <- function(fits, round, ...) {
       bfi(
             theta_hats = lapply(fits, `[[`, "theta_hat"),
             A_hats = lapply(fits, `[[`, "A_hat"),
             Lambda = fits[[1]]$Lambda,
             family = "binomial",
             treat_round = round,
             ...
         )
}

fit_r1_bin <- lapply(centers_bin, local_fit, round = "first")
fitbfi_r1_bin <- central_fit(fit_r1_bin, round = "first")


fit_r2_bin <- lapply(
       centers_bin,
       local_fit,
       round = "second",
       gamma_bfi = fitbfi_r1_bin$theta_hat
 )
fitbfi_r2_bin <- central_fit(
     fit_r2_bin,
     round = "second",
     for_ATE = lapply(fit_r2_bin, `[[`, "for_ATE")
 )

print(fitbfi_r2_bin$theta_hat)
print(fitbfi_r2_bin$Ave_Treat)
summary(fitbfi_r2_bin)


check_value("7.3", "binomial treatment coefficient",
            as.numeric(coef(fitbfi_r2_bin)["treatment"]), 2.6191)
check_value("7.3", "binomial IPTW", fitbfi_r2_bin$Ave_Treat$IPTW, 0.3365)
check_value("7.3", "binomial wIPTW", fitbfi_r2_bin$Ave_Treat$wIPTW, 0.3364)

# ----------------------------------------------------------------------------
# Section 7.4. Randomized studies with known treatment probabilities
# ----------------------------------------------------------------------------
# Reuse the same simulated binomial data from Section 7.3. Known assignment
# probability is 0.5 for every participant, so no propensity-model round is run.


fit_rct <- lapply(
       centers_bin,
       function(d) local_fit(
             d,
             round = "second",
             RCT_propens = rep(0.5, nrow(d$X))
         )
 )
fitbfi_rct <- central_fit(
       fit_rct,
       round = "second",
       for_ATE = lapply(fit_rct, `[[`, "for_ATE") )

summary(fitbfi_rct)



check_value("7.4", "RCT treatment coefficient",
            as.numeric(coef(fitbfi_rct)["treatment"]), 2.9070)
check_value("7.4", "RCT IPTW", fitbfi_rct$Ave_Treat$IPTW, 0.6000)
check_value("7.4", "RCT wIPTW", fitbfi_rct$Ave_Treat$wIPTW, 0.3714)
check_value("7.4", "RCT number assigned to treatment",
            sum(vapply(centers_bin, function(z)
              sum(as.character(z$X$treatment) == "2"), numeric(1))), 172, 0)

# ============================================================================
# Section 8. Numerical validation against pooled analyses
# ============================================================================

# ----------------------------------------------------------------------------
# Section 8.1. Gaussian regression under a homogeneous partition
# ----------------------------------------------------------------------------

data("Nurses", package = "BFI")

hospital_sizes <- as.integer(table(Nurses$hospital))
set.seed(11235)
random_order <- sample(seq_len(nrow(Nurses)))

Nurses_homo <- Nurses[random_order, ]
Nurses_homo$hospital_random <- rep(
  seq_along(hospital_sizes),
  times = hospital_sizes
)

X_combined_g <- data.frame(
  age = Nurses_homo$age,
  gender = Nurses_homo$gender,
  exp = Nurses_homo$experien
)

y_combined_g <- Nurses_homo$stress

Lambda_combined_g <- inv.prior.cov(
  X = X_combined_g,
  lambda = 0.01,
  family = "gaussian"
)

L_g <- length(hospital_sizes)
X_list_g <- vector("list", L_g)
y_list_g <- vector("list", L_g)
Lambda_list_g <- vector("list", L_g)
fit_local_g <- vector("list", L_g)

for (l in seq_len(L_g)) {
  id <- Nurses_homo$hospital_random == l

  X_list_g[[l]] <- data.frame(
    age = Nurses_homo$age[id],
    gender = Nurses_homo$gender[id],
    exp = Nurses_homo$experien[id]
  )

  y_list_g[[l]] <- Nurses_homo$stress[id]

  Lambda_list_g[[l]] <- inv.prior.cov(
    X = X_list_g[[l]],
    lambda = 0.01,
    family = "gaussian"
  )

  fit_local_g[[l]] <- MAP.estimation(
    y = y_list_g[[l]],
    X = X_list_g[[l]],
    family = "gaussian",
    Lambda = Lambda_list_g[[l]],
    control = list(maxit = 1000)
  )
}

stopifnot(all(vapply(fit_local_g, function(x) x$convergence, numeric(1)) == 0))

fit_bfi_g <- bfi(
  theta_hats = lapply(fit_local_g, function(x) x$theta_hat),
  A_hats = lapply(fit_local_g, function(x) x$A_hat),
  Lambda = append(lapply(fit_local_g, function(x) x$Lambda), list(Lambda_combined_g)),
  family = "gaussian"
)
summary(fit_bfi_g)

fit_pooled_g <- MAP.estimation(
  y = y_combined_g,
  X = X_combined_g,
  family = "gaussian",
  Lambda = Lambda_combined_g,
  control = list(maxit = 1000)
)
stopifnot(fit_pooled_g$convergence == 0)

# Regression coefficients
p_g <- length(fit_bfi_g$theta_hat) - 1

bfi_beta_g <- as.numeric(fit_bfi_g$theta_hat)[seq_len(p_g)]
pooled_beta_g <- as.numeric(fit_pooled_g$theta_hat)[seq_len(p_g)]

comparison_gaussian <- data.frame(
  parameter = colnames(fit_bfi_g$theta_hat)[seq_len(p_g)],
  BFI = bfi_beta_g,
  pooled = pooled_beta_g,
  difference = bfi_beta_g - pooled_beta_g,
  BFI_sd = as.numeric(fit_bfi_g$sd)[seq_len(p_g)],
  pooled_sd = as.numeric(fit_pooled_g$sd)[seq_len(p_g)]
)

print(comparison_gaussian, digits = 6, row.names = FALSE)

cat(
  "Maximum absolute regression-coefficient difference:",
  max(abs(comparison_gaussian$difference)),
  "\n"
)

cat(
  "Maximum absolute posterior-SD difference:",
  max(abs(comparison_gaussian$BFI_sd - comparison_gaussian$pooled_sd)),
  "\n"
)

# Dispersion parameter
sigma2_bfi_g <- tail(as.numeric(fit_bfi_g$theta_hat), 1)
sigma2_pooled_g <- tail(as.numeric(fit_pooled_g$theta_hat), 1)

cat("BFI sigma2:", sigma2_bfi_g, "\n")
cat("Pooled sigma2:", sigma2_pooled_g, "\n")
cat(
  "Absolute sigma2 difference:",
  abs(sigma2_bfi_g - sigma2_pooled_g),
  "\n"
)

c(
  BFI_sd = sqrt(diag(vcov(fit_bfi_g)))["sigma2"],
  pooled_sd = sqrt(diag(vcov(fit_pooled_g)))["sigma2"]
)


# Table 3 from the manuscript: includes the residual variance on the correct scale.
all_bfi_g <- coef(fit_bfi_g)
all_pooled_g <- coef(fit_pooled_g)
all_sd_bfi_g <- sqrt(diag(vcov(fit_bfi_g)))
all_sd_pooled_g <- sqrt(diag(vcov(fit_pooled_g)))
stopifnot(identical(names(all_bfi_g), names(all_pooled_g)))
table_gaussian <- data.frame(
  parameter = names(all_bfi_g), BFI = as.numeric(all_bfi_g),
  pooled = as.numeric(all_pooled_g),
  difference = as.numeric(all_bfi_g - all_pooled_g),
  BFI_sd = as.numeric(all_sd_bfi_g[names(all_bfi_g)]),
  pooled_sd = as.numeric(all_sd_pooled_g[names(all_bfi_g)])
)
print(table_gaussian, digits = 5, row.names = FALSE)
check_value("8.1", "Gaussian BFI residual variance", sigma2_bfi_g, 0.7659)
check_value("8.1", "Gaussian pooled residual variance", sigma2_pooled_g, 0.8718)
check_value("8.1", "max regression coefficient difference",
            max(abs(comparison_gaussian$difference)), 0.0913)

# ----------------------------------------------------------------------------
# Section 8.2. Logistic regression with the trauma data
# ----------------------------------------------------------------------------

data("trauma", package = "BFI")

trauma$age <- as.numeric(scale(trauma$age))
trauma$ISS <- as.numeric(scale(trauma$ISS))
trauma$GCS <- as.numeric(scale(trauma$GCS))
trauma$hospital <- factor(trauma$hospital)

centers <- levels(trauma$hospital)

X_list <- lapply(
  centers,
  function(l) {
    subset(
      trauma,
      hospital == l,
      select = c(sex, age, ISS, GCS)
    )
  }
)

y_list <- lapply(
  centers,
  function(l) trauma$mortality[trauma$hospital == l]
)

Lambda_list <- lapply(
  X_list,
  function(X) {
    inv.prior.cov(
      X,
      lambda = 0.01,
      family = "binomial"
    )
  }
)

fit_local <- lapply(
  seq_along(centers),
  function(l) {
    MAP.estimation(
      y = y_list[[l]],
      X = X_list[[l]],
      family = "binomial",
      Lambda = Lambda_list[[l]]
    )
  }
)

stopifnot(all(vapply(fit_local, function(x) x$convergence, numeric(1)) == 0))

X_combined <- trauma[, c("sex", "age", "ISS", "GCS")]

Lambda_combined <- inv.prior.cov(
  X_combined,
  lambda = 0.01,
  family = "binomial"
)

fit_bfi_bin <- bfi(
  theta_hats = lapply(fit_local, function(x) x$theta_hat),
  A_hats = lapply(fit_local, function(x) x$A_hat),
  Lambda = c(Lambda_list, list(Lambda_combined)),
  family = "binomial"
)

fit_pooled_bin <- MAP.estimation(
  y = trauma$mortality,
  X = X_combined,
  family = "binomial",
  Lambda = Lambda_combined
)
stopifnot(fit_pooled_bin$convergence == 0)

bfi_est_t <- as.numeric(fit_bfi_bin$theta_hat)
names(bfi_est_t) <- colnames(fit_bfi_bin$theta_hat)
pooled_est_t <- as.numeric(fit_pooled_bin$theta_hat)
names(pooled_est_t) <- names(fit_pooled_bin$theta_hat)

comparison_binomial <- data.frame(
  parameter = names(bfi_est_t),
  BFI = bfi_est_t,
  pooled = pooled_est_t,
  difference = bfi_est_t - pooled_est_t
)

print(comparison_binomial, digits = 6, row.names = FALSE)
cat("Maximum absolute estimate difference:", max(abs(bfi_est_t - pooled_est_t)), "\n")
cat(
  "Maximum absolute posterior-SD difference:",
  max(abs(as.numeric(fit_bfi_bin$sd) - as.numeric(fit_pooled_bin$sd))),
  "\n"
)


# Table 4 from the manuscript.
table_logistic <- data.frame(
  parameter = names(bfi_est_t),
  BFI = unname(bfi_est_t), pooled = unname(pooled_est_t),
  difference = unname(bfi_est_t - pooled_est_t),
  BFI_sd = as.numeric(sqrt(diag(vcov(fit_bfi_bin)))),
  pooled_sd = as.numeric(sqrt(diag(vcov(fit_pooled_bin))))
)
print(table_logistic, digits = 5, row.names = FALSE)
check_value("8.2", "max logistic coefficient difference",
            max(abs(bfi_est_t - pooled_est_t)), 0.2606)
check_value("8.2", "max logistic posterior SD difference",
            max(abs(table_logistic$BFI_sd - table_logistic$pooled_sd)), 0.0197)

# ----------------------------------------------------------------------------
# Section 8.3. Survival treatment-effect validation
# ----------------------------------------------------------------------------
# Restore the main Rotterdam allocation objects saved in Section 7.2.

df <- rotterdam_main$df
groups <- rotterdam_main$groups
fitbfi_r1 <- rotterdam_main$fitbfi_r1
fitbfi_r2 <- rotterdam_main$fitbfi_r2

# Pooled Cox model using BFI-derived propensity scores
X_ps <- model.matrix(
  ~ year + age + meno + size + grade + nodes + pgr + er + hormon,
  data = df
)

gamma_bfi <- coef(fitbfi_r1)

stopifnot(identical(colnames(X_ps), names(gamma_bfi)))

ps_bfi <- as.numeric(plogis(X_ps %*% gamma_bfi))
Wprop_bfi <- df$chemo / ps_bfi + (1 - df$chemo) / (1 - ps_bfi)

df_pooled <- df
df_pooled$center <- factor(groups)

fit_pooled_bfi_ps <- coxph(
  Surv(time, status) ~ chemo + strata(center),
  data = df_pooled,
  weights = Wprop_bfi,
  robust = TRUE,
  ties = "breslow"
)

coef(fit_pooled_bfi_ps)["chemo"]

# Pooled Cox model using a separately fitted pooled propensity model
fit_ps_pooled <- glm(
  chemo ~ year + age + meno + size + grade + nodes + pgr + er + hormon,
  data = df,
  family = binomial()
)

ps_pooled <- predict(fit_ps_pooled, type = "response")
Wprop_pooled <- df$chemo / ps_pooled + (1 - df$chemo) / (1 - ps_pooled)

fit_pooled <- coxph(
  Surv(time, status) ~ chemo + strata(center),
  data = df_pooled,
  weights = Wprop_pooled,
  robust = TRUE,
  ties = "breslow"
)

coef(fit_pooled)["chemo"]

summary(fit_pooled_bfi_ps)$coefficients

# Sensitivity allocation: random allocation stratified by treatment group
set.seed(11235)

df_alt <- df
df_alt$orig_row <- seq_len(nrow(df_alt))

treated <- df_alt[df_alt$chemo == 1, ]
control <- df_alt[df_alt$chemo == 0, ]

treated <- treated[sample(nrow(treated)), ]
control <- control[sample(nrow(control)), ]

treated_split <- vector("list", L)
control_split <- vector("list", L)

for (l in seq_len(L)) {
  treated_split[[l]] <- treated[seq(l, nrow(treated), by = L), , drop = FALSE]
  control_split[[l]] <- control[seq(l, nrow(control), by = L), , drop = FALSE]
}

datasets_alt <- lapply(
  seq_len(L),
  function(l) rbind(treated_split[[l]], control_split[[l]])
)

groups_alt <- rep(NA_integer_, nrow(df_alt))
for (l in seq_len(L)) {
  groups_alt[datasets_alt[[l]]$orig_row] <- l
  datasets_alt[[l]]$orig_row <- NULL
}
df_alt$orig_row <- NULL

print(table(groups_alt, df_alt$chemo))

X_list_alt <- lapply(
  datasets_alt,
  function(dat) dat[, !(names(dat) %in% c("time", "status")), drop = FALSE]
)

Y_list_alt <- lapply(
  datasets_alt,
  function(dat) dat[, c("time", "status"), drop = FALSE]
)

# First round under alternative allocation
Lambda_r1_alt <- lapply(
  X_list_alt,
  function(X) {
    inv.prior.cov(
      X,
      lambda = 0.01,
      family = "survival",
      treatment = "chemo",
      treat_round = "first"
    )
  }
)

fit_r1_alt <- lapply(
  seq_len(L),
  function(l) {
    MAP.estimation(
      Y_list_alt[[l]],
      X_list_alt[[l]],
      family = "survival",
      Lambda = Lambda_r1_alt[[l]],
      treatment = "chemo",
      treat_round = "first"
    )
  }
)

fitbfi_r1_alt <- bfi(
  theta_hats = lapply(fit_r1_alt, function(x) x$theta_hat),
  A_hats = lapply(fit_r1_alt, function(x) x$A_hat),
  Lambda = Lambda_r1_alt[[1]],
  family = "survival",
  treat_round = "first"
)

# Second round under alternative allocation
Lambda_r2_alt <- lapply(
  X_list_alt,
  function(X) {
    inv.prior.cov(
      X,
      lambda = 0.01,
      family = "survival",
      basehaz = "unspecified",
      treatment = "chemo",
      treat_round = "second"
    )
  }
)

fit_r2_alt <- lapply(
  seq_len(L),
  function(l) {
    MAP.estimation(
      Y_list_alt[[l]],
      X_list_alt[[l]],
      family = "survival",
      Lambda = Lambda_r2_alt[[l]],
      basehaz = "unspecified",
      treatment = "chemo",
      treat_round = "second",
      gamma_bfi = fitbfi_r1_alt$theta_hat
    )
  }
)

fitbfi_r2_alt <- bfi(
  theta_hats = lapply(fit_r2_alt, function(x) x$theta_hat),
  A_hats = lapply(fit_r2_alt, function(x) x$A_hat),
  Lambda = Lambda_r2_alt[[1]],
  family = "survival",
  basehaz = "unspecified",
  treat_round = "second"
)

# Pooled validation under alternative allocation using BFI propensity scores
X_ps_alt <- model.matrix(
  ~ year + age + meno + size + grade + nodes + pgr + er + hormon,
  data = df_alt
)

gamma_bfi_alt <- as.numeric(fitbfi_r1_alt$theta_hat)
names(gamma_bfi_alt) <- colnames(fitbfi_r1_alt$theta_hat)
stopifnot(identical(colnames(X_ps_alt), names(gamma_bfi_alt)))

ps_bfi_alt <- as.numeric(plogis(X_ps_alt %*% gamma_bfi_alt))
Wprop_bfi_alt <- df_alt$chemo / ps_bfi_alt +
  (1 - df_alt$chemo) / (1 - ps_bfi_alt)

df_pooled_alt <- df_alt
df_pooled_alt$center <- factor(groups_alt)

fit_pooled_bfi_ps_alt <- coxph(
  Surv(time, status) ~ chemo + strata(center),
  data = df_pooled_alt,
  weights = Wprop_bfi_alt,
  robust = TRUE,
  ties = "breslow"
)

# Pooled validation under alternative allocation using pooled propensity scores
fit_ps_pooled_alt <- glm(
  chemo ~ year + age + meno + size + grade + nodes + pgr + er + hormon,
  data = df_alt,
  family = binomial()
)

ps_pooled_alt <- predict(fit_ps_pooled_alt, type = "response")
Wprop_pooled_alt <- df_alt$chemo / ps_pooled_alt +
  (1 - df_alt$chemo) / (1 - ps_pooled_alt)

fit_pooled_ps_alt <- coxph(
  Surv(time, status) ~ chemo + strata(center),
  data = df_pooled_alt,
  weights = Wprop_pooled_alt,
  robust = TRUE,
  ties = "breslow"
)

# Table corresponding to the manuscript survival-validation comparison
validation_survival <- data.frame(
  allocation = c("Unstratified by treatment", "Stratified by treatment"),
  BFI = c(
    as.numeric(fitbfi_r2$theta_hat[1, "chemo"]),
    as.numeric(fitbfi_r2_alt$theta_hat[1, "chemo"])
  ),
  pooled_BFI_propensity = c(
    as.numeric(coef(fit_pooled_bfi_ps)["chemo"]),
    as.numeric(coef(fit_pooled_bfi_ps_alt)["chemo"])
  ),
  pooled_pooled_propensity = c(
    as.numeric(coef(fit_pooled)["chemo"]),
    as.numeric(coef(fit_pooled_ps_alt)["chemo"])
  )
)

print(validation_survival, digits = 7, row.names = FALSE)
abs(validation_survival$BFI - validation_survival$pooled_BFI_propensity)
abs(validation_survival$BFI - validation_survival$pooled_pooled_propensity)




#####################################################
#### Plots
#####################################################
tmin  <- min(c(y1$time, y2$time))
tmax  <- min(max(y1$time), max(y2$time))
times <- seq(tmin, tmax, length.out = 200)
bh <- function(fit) {
  hazards.fun(
    time = times,
    p = p,
    theta_hat = fit$theta_hat,
    basehaz = "weibul"
  )$bhazard
}
h1 <- bh(fit_weib1); h2 <- bh(fit_weib2); hb <- bh(fit_weib_bfi)
# range(c(h1, h2, hb))

cols <- c("#0072B2", "#E69F00", "#000000")

op <- par(mar = c(4, 4, 1, 1), mgp = c(2.4, 0.8, 0))
matplot(times, cbind(h1, h2, hb), type = "l", log = "y",
        lty = c(2, 3, 1), lwd = c(1.8, 1.8, 2.4), col = cols,
        xlab = "Time", ylab = "Baseline hazard (log scale)")
legend("topleft", bty = "n",
       lty = c(2, 3, 1), lwd = c(1.8, 1.8, 2.4), col = cols,
       legend = c("Center 1", "Center 2", "BFI"))
par(op)

#----

# Numerical verification of Section 8.3 pooled/alternative comparisons.
check_value("8.3", "pooled Cox with BFI propensity (main)",
            coef(fit_pooled_bfi_ps)["chemo"], -0.1992)
check_value("8.3", "pooled Cox with pooled propensity (main)",
            coef(fit_pooled)["chemo"], -0.2034)
check_value("8.3", "BFI Cox coefficient (alternative)",
            as.numeric(coef(fitbfi_r2_alt)["chemo"]), -0.2004)
check_value("8.3", "pooled BFI-propensity Cox (alternative)",
            coef(fit_pooled_bfi_ps_alt)["chemo"], -0.2023)
check_value("8.3", "pooled pooled-propensity Cox (alternative)",
            coef(fit_pooled_ps_alt)["chemo"], -0.2068)

# Compare the aggregate Weibull hazard with the local-hazard envelope.
outside <- hb < pmin(h1, h2) | hb > pmax(h1, h2)
cat("Weibull Figure: fraction outside two local curves = ",
    mean(outside), "\n", sep = "")

# Numerical checks are held in memory and printed to the R console.
cat("\nNumerical verification: ", sum(checks$passed), "/",
    nrow(checks), " passed.\n", sep = "")
if (!all(checks$passed)) {
  stop("One or more manuscript values do not match within tolerance. ",
       "Review the [MISMATCH] messages above.")
}

# ============================================================================
# End
# ============================================================================

cat("\nSession information:\n")
print(sessionInfo())
