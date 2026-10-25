rm(list = ls())

set.seed(12345)

out_dir <- "reviewer5_repo_benchmark"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Initial extracellular Legionella counts to test.
initial_doses <- c(1, 5, 10, 50, 100)

# Number of successful threshold-crossing times to collect per model/dose.
# Start with 100 to test. Use 1000+ for final manuscript results.
target_successes <- 10000

# Maximum attempts per model/dose. Low-dose stochastic models may go extinct often.
max_attempts <- 500000

# Model stopping threshold for symptom-associated bacterial load.
LT <- 50661

# Time horizon.
t_max <- 500

# Max Gillespie events.
max_events <- 1e7

# ODE and SDE time steps.
ode_dt <- 0.01
sde_dt <- 0.1   # this matches your repo style in random sde

# Switch states from the current manuscript.
# These are the state-space thresholds for switching.
L_sde_switch <- 49.13
M_sde_switch <- 25.01

L_ode_switch <- 92.14
M_ode_switch <- 46.91

# If you want to reproduce older repo thresholds instead, set these manually,
# e.g. L_sde_switch <- 33.35577; M_sde_switch <- 7.00475
# and L_ode_switch <- 181.9034; M_ode_switch <- 38.19993.

a <- 0.089      # alpha
b <- 1.088      # beta
h <- 0.041      # lambda
G <- 62.792
r <- 32.072

time_it <- function(expr) {
  start <- proc.time()[["elapsed"]]
  value <- eval.parent(substitute(expr))
  end <- proc.time()[["elapsed"]]
  list(value = value, seconds = end - start)
}

safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
safe_median <- function(x) if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)
safe_sd <- function(x) if (sum(!is.na(x)) <= 1) NA_real_ else sd(x, na.rm = TRUE)
safe_mae <- function(x) if (all(is.na(x))) NA_real_ else mean(abs(x), na.rm = TRUE)
safe_rmse <- function(x) if (all(is.na(x))) NA_real_ else sqrt(mean(x^2, na.rm = TRUE))
safe_q <- function(x, p) {
  if (all(is.na(x))) return(NA_real_)
  as.numeric(quantile(x, probs = p, na.rm = TRUE, names = FALSE))
}

  ode_step_rhs <- function(L, M) {
  c(
    dL = -(a + b) * L + h * G * M,
    dM = a * L - h * M
  )
}

ode_simulate_until <- function(
    dose,
    L0 = dose,
    M0 = 0,
    target = c("LT", "SDE_SWITCH", "ODE_SWITCH"),
    t_max = 500,
    dt = 0.01
) {
  target <- match.arg(target)
  
  n_steps <- ceiling(t_max / dt)
  
  L <- L0
  M <- M0
  t <- 0
  
  target_hit <- function(L, M) {
    if (target == "LT") {
      return(L >= LT)
    }
    if (target == "SDE_SWITCH") {
      return(L >= L_sde_switch && M >= M_sde_switch)
    }
    if (target == "ODE_SWITCH") {
      return(L >= L_ode_switch && M >= M_ode_switch)
    }
    FALSE
  }
  
  if (target_hit(L, M)) {
    return(list(status = target, time = 0, L = L, M = M))
  }
  
  # RK4 integration to avoid requiring deSolve.
  for (i in seq_len(n_steps)) {
    k1 <- ode_step_rhs(L, M)
    k2 <- ode_step_rhs(L + 0.5 * dt * k1[1], M + 0.5 * dt * k1[2])
    k3 <- ode_step_rhs(L + 0.5 * dt * k2[1], M + 0.5 * dt * k2[2])
    k4 <- ode_step_rhs(L + dt * k3[1], M + dt * k3[2])
    
    L_new <- L + dt * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6
    M_new <- M + dt * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]) / 6
    
    t <- t + dt
    
    L <- max(L_new, 0)
    M <- max(M_new, 0)
    
    if (target_hit(L, M)) {
      return(list(status = target, time = t, L = L, M = M))
    }
    
    if (L <= 0 && M <= 0) {
      return(list(status = "extinct", time = NA_real_, L = L, M = M))
    }
  }
  
  list(status = "max_time", time = NA_real_, L = L, M = M)
}

run_ode_only_once <- function(dose) {
  ode_simulate_until(dose = dose, target = "LT", t_max = t_max, dt = ode_dt)
}

mtbp_one_step <- function(time, L, M) {
  rates <- c(
    sphago = a * L,
    fphago = b * L,
    burst = h * M
  )
  
  total_rates <- sum(rates)
  
  if (total_rates <= 0) {
    return(list(time = time, L = L, M = M, status = "extinct"))
  }
  
  tau <- rexp(n = 1, rate = total_rates)
  event <- sample.int(n = 3, size = 1, prob = rates / total_rates)
  
  if (event == 1) {
    # successful phagocytosis
    L <- L - 1
    M <- M + 1
  } else if (event == 2) {
    # failed phagocytosis
    L <- L - 1
  } else {
    # burst
    L <- L + rnbinom(n = 1, mu = G, size = r)
    M <- M - 1
  }
  
  L <- max(L, 0)
  M <- max(M, 0)
  
  list(time = time + tau, L = L, M = M, status = "running")
}

mtbp_simulate_until <- function(
    dose,
    target = c("LT", "SDE_SWITCH", "ODE_SWITCH"),
    t_max = 500,
    max_events = 1e7
) {
  target <- match.arg(target)
  
  time <- 0
  L <- dose
  M <- 0
  
  target_hit <- function(L, M) {
    if (target == "LT") {
      return(L >= LT)
    }
    if (target == "SDE_SWITCH") {
      return(L >= L_sde_switch && M >= M_sde_switch)
    }
    if (target == "ODE_SWITCH") {
      return(L >= L_ode_switch && M >= M_ode_switch)
    }
    FALSE
  }
  
  if (target_hit(L, M)) {
    return(list(status = target, time = time, L = L, M = M, events = 0))
  }
  
  for (event_i in seq_len(max_events)) {
    if (L >= LT) {
      return(list(status = "LT", time = time, L = L, M = M, events = event_i))
    }
    
    if (target_hit(L, M)) {
      return(list(status = target, time = time, L = L, M = M, events = event_i))
    }
    
    if (L <= 0 && M <= 0) {
      return(list(status = "extinct", time = NA_real_, L = L, M = M, events = event_i))
    }
    
    if (time >= t_max) {
      return(list(status = "max_time", time = NA_real_, L = L, M = M, events = event_i))
    }
    
    step <- mtbp_one_step(time, L, M)
    
    if (step$status == "extinct") {
      return(list(status = "extinct", time = NA_real_, L = L, M = M, events = event_i))
    }
    
    time <- step$time
    L <- step$L
    M <- step$M
  }
  
  list(status = "max_events", time = NA_real_, L = L, M = M, events = max_events)
}

run_mtbp_reference_once <- function(dose) {
  mtbp_simulate_until(dose = dose, target = "LT", t_max = t_max, max_events = max_events)
}

sde_simulate_until <- function(
    L0,
    M0,
    target = c("LT", "ODE_SWITCH"),
    dt = 0.1,
    t_max = 500
) {
  target <- match.arg(target)
  
  n_steps <- ceiling(t_max / dt)
  
  l <- L0
  m <- M0
  time <- 0
  
  target_hit <- function(l, m) {
    if (target == "LT") {
      return(l >= LT)
    }
    if (target == "ODE_SWITCH") {
      return(l >= L_ode_switch && m >= M_ode_switch)
    }
    FALSE
  }
  
  if (target_hit(l, m)) {
    return(list(status = target, time = 0, L = l, M = m))
  }
  
  for (i in seq_len(n_steps)) {
    if (l >= LT) {
      return(list(status = "LT", time = time, L = l, M = m))
    }
    
    if (target_hit(l, m)) {
      return(list(status = target, time = time, L = l, M = m))
    }
    
    if (l <= 0 && m <= 0) {
      return(list(status = "extinct", time = NA_real_, L = l, M = m))
    }
    
    noise1 <- rnorm(n = 1, mean = 0, sd = sqrt(dt))
    noise2 <- rnorm(n = 1, mean = 0, sd = sqrt(dt))
    
    A <- ((a + b) * l + h * (G + (1 / r + 1) * G^2) * m) *
      (a * l + h * m) -
      (a * l + h * G * m)^2
    
    if (is.na(A) || A < 0) A <- 0
    B <- sqrt(A)
    
    C <- (2 * a + b) * l + h * (G + (1 / r + 1) * G^2 + 1) * m
    D <- sqrt(C + 2 * B)
    
    if (is.na(D) || D <= 0) {
      # drift-only fallback if noise matrix degenerates
      l_new <- l - dt * (a + b) * l + dt * h * G * m
      m_new <- m + dt * a * l - dt * h * m
    } else {
      l_new <- l -
        dt * (a + b) * l +
        dt * h * G * m +
        ((a + b) * l + h * (G + (1 / r + 1) * G^2) * m + B) * noise1 / D -
        (a * l + h * G * m) * noise2 / D
      
      m_new <- m +
        dt * a * l -
        dt * h * m -
        (a * l + h * G * m) * noise1 / D +
        (a * l + h * m + B) * noise2 / D
    }
    
    if (is.na(l_new)) l_new <- 0
    if (is.na(m_new)) m_new <- 0
    
    l <- max(l_new, 0)
    m <- max(m_new, 0)
    
    time <- time + dt
  }
  
  list(status = "max_time", time = NA_real_, L = l, M = m)
}

run_sde_only_once <- function(dose) {
  sde_simulate_until(L0 = dose, M0 = 0, target = "LT", dt = sde_dt, t_max = t_max)
}

run_mtbp_to_ode_once <- function(dose) {
  first <- mtbp_simulate_until(
    dose = dose,
    target = "ODE_SWITCH",
    t_max = t_max,
    max_events = max_events
  )
  
  if (first$status == "LT") {
    return(first)
  }
  
  if (first$status != "ODE_SWITCH") {
    return(first)
  }
  
  second <- ode_simulate_until(
    dose = dose,
    L0 = first$L,
    M0 = first$M,
    target = "LT",
    t_max = t_max - first$time,
    dt = ode_dt
  )
  
  if (second$status == "LT") {
    return(list(
      status = "LT",
      time = first$time + second$time,
      L = second$L,
      M = second$M
    ))
  }
  
  list(status = second$status, time = NA_real_, L = second$L, M = second$M)
}

run_mtbp_to_sde_to_ode_once <- function(dose) {
  first <- mtbp_simulate_until(
    dose = dose,
    target = "SDE_SWITCH",
    t_max = t_max,
    max_events = max_events
  )
  
  if (first$status == "LT") {
    return(first)
  }
  
  if (first$status != "SDE_SWITCH") {
    return(first)
  }
  
  second <- sde_simulate_until(
    L0 = first$L,
    M0 = first$M,
    target = "ODE_SWITCH",
    dt = sde_dt,
    t_max = t_max - first$time
  )
  
  if (second$status == "LT") {
    return(list(
      status = "LT",
      time = first$time + second$time,
      L = second$L,
      M = second$M
    ))
  }
  
  if (second$status != "ODE_SWITCH") {
    return(list(status = second$status, time = NA_real_, L = second$L, M = second$M))
  }
  
  third <- ode_simulate_until(
    dose = dose,
    L0 = second$L,
    M0 = second$M,
    target = "LT",
    t_max = t_max - first$time - second$time,
    dt = ode_dt
  )
  
  if (third$status == "LT") {
    return(list(
      status = "LT",
      time = first$time + second$time + third$time,
      L = third$L,
      M = third$M
    ))
  }
  
  list(status = third$status, time = NA_real_, L = third$L, M = third$M)
}

simulate_success_times <- function(model_name, dose, target_successes, max_attempts) {
  times <- numeric(0)
  statuses <- character(0)
  
  total_seconds <- 0
  attempts <- 0
  
  model_fun <- switch(
    model_name,
    "MTBP reference" = run_mtbp_reference_once,
    "ODE-only" = run_ode_only_once,
    "SDE-only" = run_sde_only_once,
    "MTBP-to-ODE hybrid" = run_mtbp_to_ode_once,
    "MTBP-to-SDE-to-ODE hybrid" = run_mtbp_to_sde_to_ode_once,
    stop("Unknown model name")
  )
  
  while (length(times) < target_successes && attempts < max_attempts) {
    attempts <- attempts + 1
    
    timed <- time_it(model_fun(dose))
    total_seconds <- total_seconds + timed$seconds
    
    res <- timed$value
    statuses <- c(statuses, res$status)
    
    if (res$status == "LT" && !is.na(res$time)) {
      times <- c(times, res$time)
    }
    
    # ODE-only is deterministic, so repeat it without actually rerunning forever.
    if (model_name == "ODE-only" && attempts == 1) {
      if (length(times) == 1) {
        times <- rep(times[1], target_successes)
        statuses <- rep("LT", target_successes)
      }
      break
    }
    
    if (attempts %% 1000 == 0) {
      cat(
        "    ", model_name,
        "| dose", dose,
        "| attempts", attempts,
        "| successes", length(times),
        "\n"
      )
    }
  }
  
  list(
    model = model_name,
    dose = dose,
    times = times,
    attempts = attempts,
    successes = length(times),
    statuses = statuses,
    total_seconds = total_seconds,
    seconds_per_attempt = total_seconds / max(attempts, 1),
    seconds_per_success = total_seconds / max(length(times), 1)
  )
}

models <- c(
  "MTBP reference",
  "ODE-only",
  "SDE-only",
  "MTBP-to-ODE hybrid",
  "MTBP-to-SDE-to-ODE hybrid"
)

all_results <- list()
idx <- 1

for (dose in initial_doses) {
  cat("\n============================\n")
  cat("Initial dose:", dose, "\n")
  cat("============================\n")
  
  for (model_name in models) {
    cat("  Running:", model_name, "\n")
    
    all_results[[idx]] <- simulate_success_times(
      model_name = model_name,
      dose = dose,
      target_successes = target_successes,
      max_attempts = max_attempts
    )
    
    cat(
      "  Finished:", model_name,
      "| successes:", all_results[[idx]]$successes,
      "| attempts:", all_results[[idx]]$attempts,
      "| sec/success:", round(all_results[[idx]]$seconds_per_success, 4),
      "\n"
    )
    
    idx <- idx + 1
  }
}

raw_times <- do.call(
  rbind,
  lapply(all_results, function(x) {
    if (length(x$times) == 0) {
      return(data.frame(
        dose = x$dose,
        model = x$model,
        replicate = integer(0),
        time_hours = numeric(0),
        time_days = numeric(0)
      ))
    }
    
    data.frame(
      dose = x$dose,
      model = x$model,
      replicate = seq_along(x$times),
      time_hours = x$times,
      time_days = x$times / 24
    )
  })
)

write.csv(
  raw_times,
  file.path(out_dir, "raw_times.csv"),
  row.names = FALSE
)

timing_summary <- do.call(
  rbind,
  lapply(all_results, function(x) {
    data.frame(
      dose = x$dose,
      model = x$model,
      attempts = x$attempts,
      successes = x$successes,
      success_probability = x$successes / max(x$attempts, 1),
      total_seconds = x$total_seconds,
      seconds_per_attempt = x$seconds_per_attempt,
      seconds_per_success = x$seconds_per_success
    )
  })
)

write.csv(
  timing_summary,
  file.path(out_dir, "timing_summary.csv"),
  row.names = FALSE
)

get_times <- function(dose, model) {
  raw_times$time_days[raw_times$dose == dose & raw_times$model == model]
}

wasserstein_1d <- function(x, y, n_grid = 1000) {
  if (length(x) == 0 || length(y) == 0) return(NA_real_)
  u <- seq(0.001, 0.999, length.out = n_grid)
  qx <- quantile(x, probs = u, names = FALSE, na.rm = TRUE)
  qy <- quantile(y, probs = u, names = FALSE, na.rm = TRUE)
  mean(abs(qx - qy))
}

dist_summary_rows <- list()
k <- 1

for (dose in initial_doses) {
  ref <- get_times(dose, "MTBP reference")
  
  for (model_name in models[models != "MTBP reference"]) {
    model_times <- get_times(dose, model_name)
    
    if (length(model_times) == 0 || length(ref) == 0) {
      ks_stat <- NA_real_
      ks_p <- NA_real_
    } else {
      ks <- suppressWarnings(ks.test(model_times, ref))
      ks_stat <- unname(ks$statistic)
      ks_p <- ks$p.value
    }
    
    dist_summary_rows[[k]] <- data.frame(
      dose = dose,
      model = model_name,
      
      n_reference = length(ref),
      n_model = length(model_times),
      
      mean_reference_days = safe_mean(ref),
      mean_model_days = safe_mean(model_times),
      mean_error_days = safe_mean(model_times) - safe_mean(ref),
      
      median_reference_days = safe_median(ref),
      median_model_days = safe_median(model_times),
      median_error_days = safe_median(model_times) - safe_median(ref),
      
      q025_reference_days = safe_q(ref, 0.025),
      q025_model_days = safe_q(model_times, 0.025),
      q025_error_days = safe_q(model_times, 0.025) - safe_q(ref, 0.025),
      
      q975_reference_days = safe_q(ref, 0.975),
      q975_model_days = safe_q(model_times, 0.975),
      q975_error_days = safe_q(model_times, 0.975) - safe_q(ref, 0.975),
      
      wasserstein_days = wasserstein_1d(model_times, ref),
      ks_statistic = ks_stat,
      ks_p_value = ks_p
    )
    
    k <- k + 1
  }
}

distribution_error_summary <- do.call(rbind, dist_summary_rows)

distribution_error_summary <- merge(
  distribution_error_summary,
  timing_summary[timing_summary$model != "MTBP reference",
                 c("dose", "model", "attempts", "successes",
                   "success_probability", "seconds_per_success")],
  by = c("dose", "model"),
  all.x = TRUE
)

write.csv(
  distribution_error_summary,
  file.path(out_dir, "distribution_error_summary.csv"),
  row.names = FALSE
)

quantile_probs <- c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)

quantile_rows <- list()
qk <- 1

for (dose in initial_doses) {
  ref <- get_times(dose, "MTBP reference")
  
  for (model_name in models[models != "MTBP reference"]) {
    model_times <- get_times(dose, model_name)
    
    for (p in quantile_probs) {
      quantile_rows[[qk]] <- data.frame(
        dose = dose,
        model = model_name,
        quantile = p,
        reference_days = safe_q(ref, p),
        model_days = safe_q(model_times, p),
        error_days = safe_q(model_times, p) - safe_q(ref, p)
      )
      qk <- qk + 1
    }
  }
}

quantile_error_summary <- do.call(rbind, quantile_rows)

write.csv(
  quantile_error_summary,
  file.path(out_dir, "quantile_error_summary.csv"),
  row.names = FALSE
)

# Boxplot of threshold-crossing time distributions.
pdf(file.path(out_dir, "threshold_time_boxplots.pdf"), width = 12, height = 8)
par(mfrow = c(ceiling(length(initial_doses) / 2), 2), mar = c(8, 4, 3, 1))

for (dose in initial_doses) {
  dd <- raw_times[raw_times$dose == dose, ]
  dd$model <- factor(dd$model, levels = models)
  
  boxplot(
    time_days ~ model,
    data = dd,
    las = 2,
    ylab = "Threshold-crossing time (days)",
    main = paste("Initial dose =", dose),
    outline = TRUE
  )
}

dev.off()

# Distribution-level error plot.
pdf(file.path(out_dir, "distribution_error_plot.pdf"), width = 10, height = 6)

plot(
  distribution_error_summary$seconds_per_success,
  distribution_error_summary$wasserstein_days,
  log = "x",
  pch = 19,
  xlab = "Seconds per successful threshold-crossing sample (log scale)",
  ylab = "Wasserstein distance from MTBP reference (days)",
  main = "Accuracy--runtime trade-off"
)

text(
  distribution_error_summary$seconds_per_success,
  distribution_error_summary$wasserstein_days,
  labels = paste0(
    distribution_error_summary$model,
    ", d=", distribution_error_summary$dose
  ),
  pos = 4,
  cex = 0.6
)

dev.off()

# Density plots by dose.
pdf(file.path(out_dir, "threshold_time_densities.pdf"), width = 12, height = 8)
par(mfrow = c(ceiling(length(initial_doses) / 2), 2), mar = c(4, 4, 3, 1))

for (dose in initial_doses) {
  dd <- raw_times[raw_times$dose == dose, ]
  
  ref <- dd$time_days[dd$model == "MTBP reference"]
  if (length(ref) > 1) {
    d_ref <- density(ref)
    plot(
      d_ref,
      main = paste("Initial dose =", dose),
      xlab = "Threshold-crossing time (days)",
      ylab = "Density",
      lwd = 2
    )
    
    for (model_name in models[models != "MTBP reference"]) {
      xx <- dd$time_days[dd$model == model_name]
      if (length(xx) > 1) {
        lines(density(xx), lwd = 2, lty = which(models == model_name))
      }
    }
    
    legend(
      "topright",
      legend = models,
      lty = seq_along(models),
      lwd = 2,
      cex = 0.7,
      bty = "n"
    )
  }
}

dev.off()

cat("\nTiming summary:\n")
print(timing_summary)

cat("\nDistribution error summary:\n")
print(distribution_error_summary)
