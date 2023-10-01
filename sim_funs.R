gen_data <- function (seed, n, n_scl = 1.5) {
  set.seed(seed)
  n_target <- n # target number of observational units
  n <- n * n_scl
  n_i <- 15  # number of (planned) encounters per unit
  tmax <- 7 # maximum follow-up time (type I censoring)
  remove(n_scl)
  # longitudinal outcome 1/2
  ## parameters true values
  betas <- c("Intercept" = 6.94, "Time1" = 1.30, "Time2" = 1.84, "Time3" = 1.82)
  sigma_y <- 0.6 # measurement error sd
  D <- matrix(0, 4, 4)
  D[lower.tri(D, TRUE)] <- c(0.71, 0.33, 0.07, 1.26, 2.68, 3.81, 4.35, 7.62, 5.4, 8)
  D <- D + t(D)
  diag(D) <- diag(D) * 0.5
  b <- MASS::mvrnorm(n, rep(0, nrow(D)), D)
  Bkn <- c(0, 7)
  kn <- c(1, 3)
  remove(D)
  # terminal outcome
  ## parameters true values
  gammas_t <- c("(Intercept)" = -9, "Group" = 0.5, "Age" = 0.05) # phi = exp(Intercept)
  sigma_t <- 2
  alpha <- 0.5
  ## terminal data
  group <- rep(0:1, each = n/2)
  age <- runif(n, 30, 70)
  W_t <- cbind("(Intercept)" = 1, "Group" = group, "Age" = age)
  eta_t <- as.vector(W_t %*% gammas_t) 
  invS_t <- function(t, u, i) {
    h <- function(s) { 
      NS <- splines::ns(s, knots = kn, Boundary.knots = Bkn)
      X <- cbind(1, NS)
      Z <- cbind(1, NS)
      eta_y <- as.vector(X %*% betas + rowSums(Z * b[rep(i, nrow(Z)), ]))
      exp(log(sigma_t) + (sigma_t - 1) * log(s) + eta_t[i] + eta_y * alpha) 
    }
    integrate(h, lower = 0, upper = t)$value + log(u)
  }
  u_t <- runif(n)
  ter_times <- numeric(n)
  for(i in seq_len(n)) {
    root <- try(uniroot(invS_t, interval = c(1e-05, 250),
                        u = u_t[i], i = i)$root, TRUE)  
    ter_times[i] <- if (!inherits(root, "try-error")) root else NA
  }
  surv_na <- !is.na(ter_times)
  if(sum(surv_na) < n_target) stop("Not enough patients. Increase 'n_scl'.")
  rmv_ids <- sample(which(surv_na), sum(surv_na) - n_target)
  surv_na[rmv_ids] <- FALSE # remove the excess of units
  surv <- data.frame(id    = seq_len(n)[surv_na],
                     time  = ter_times[surv_na],
                     group = group[surv_na],
                     age   = age[surv_na])
  b <- b[surv_na, , drop = FALSE]
  cens_times <- tmax
  surv$Tstatus <- as.numeric(surv$time <= cens_times) # event indicator
  surv$time <- pmin(surv$time, cens_times) # add censoring time
  surv$id <- seq_along(surv$id)
  remove(gammas_t, sigma_t, group, W_t, eta_t, alpha, invS_t, u_t, i, root, 
         rmv_ids, ter_times, cens_times, n, age, surv_na)
  # longitudinal outcome 2/2
  long <- data.frame(id   = rep(surv$id, each = n_i),
                     time = c(replicate(length(surv$id), c(0, sort(runif(n_i - 1, 1, tmax))))))
  X <- model.matrix(~ 1 + splines::ns(time, knots = kn, Boundary.knots = Bkn), 
                    data = long)
  Z <- model.matrix(~ 1 + splines::ns(time, knots = kn, Boundary.knots = Bkn), 
                    data = long)
  eta_y <- as.vector(X %*% betas + rowSums(Z * b[long$id, ]))
  long$y <- rnorm(length(eta_y), eta_y, sigma_y)
  long_cens <- long$time <= rep(surv$time, times = rle(long$id)$lengths) 
  long <- long[long_cens, , drop = FALSE] # drop censored encounters
  long$group <- surv$group[long$id]
  remove(kn, Bkn, X, betas, Z, b, eta_y, sigma_y, n_i, tmax, long_cens)
  # return
  long$seed <- surv$seed <- seed # save seed
  return(list(long = long, surv = surv))
}

slicer <- function(n_slices, id_var, data_Long, data_Surv, seed = 123L) {
  ids_unq <- unique(c(data_Long[[id_var]], data_Surv[[id_var]]))
  set.seed(seed)
  ids_slc <- split(sample(ids_unq), (seq_along(ids_unq) %% n_slices) + 1)
  Long <- lapply(ids_slc, function(ids) 
    data_Long[data_Long[[id_var]] %in% ids, ])
  Surv <- lapply(ids_slc, function(ids) 
    data_Surv[data_Surv[[id_var]] %in% ids, ])
  class(Long) <- class(Surv) <- "sliced_data"
  list(Long = Long, Surv = Surv)
}

lme_sliced_data <- function(fixed, data = sys.frame(sys.parent()), random, correlation = NULL, 
                            weights = NULL, subset, method = c("REML", "ML"), na.action = na.fail, 
                            control = list(), contrasts = NULL, keep.data = TRUE,
                            cores = NULL) {
  args <- as.list(match.call())[-c(1L, 3L)] # remove $"" and $data 
  args["cores"] <- NULL
  if(is.null(cores)) cores <- max(parallel::detectCores() - 1, 1)
  f <- function(i) do.call(nlme::lme.formula, c(list(data = data[[i]]), args))
  cl <- parallel::makeCluster(cores)
  parallel::clusterEvalQ(cl, library(nlme))
  res <- parallel::parLapply(cl, seq_along(data), f)
  parallel::stopCluster(cl)
  structure(res, class = "sliced_lme")
}

coxph_sliced_data <- function (formula, data, weights, subset, na.action, init, control, 
                               ties = c("efron", "breslow", "exact"), singular.ok = TRUE, 
                               robust, model = FALSE, x = FALSE, y = TRUE, tt, method = ties, 
                               id, cluster, istate, statedata, nocenter = c(-1, 0, 1),
                               cores = NULL, ...) {
  args <- as.list(match.call())[-c(1L, 3L, 4L)] # remove $"" and $data
  args["cores"] <- NULL
  if(is.null(cores)) cores <- max(parallel::detectCores() - 1, 1)
  f <- function(i) do.call(survival::coxph, c(list(data = data[[i]]), args))
  cl <- parallel::makeCluster(cores)
  parallel::clusterEvalQ(cl, library(survival))
  res <- parallel::parLapply(cl, seq_along(data), f)
  parallel::stopCluster(cl)
  structure(res, class = "sliced_coxph")
}

par_diag <- function(n_slc, n_chains, ncl_in, ncl_out, cores, diff_time) {
  # slice format
  nrow <- ceiling(n_chains / ncl_in)
  ncol <- ncl_in
  v <- seq_len(n_chains)
  length(v) <- nrow * ncol
  m <- matrix(v, byrow = TRUE, nrow = nrow, ncol = ncol)
  slc <- NULL
  for(r in seq_len(nrow)) {
    slc <- paste0(slc, m[r, !is.na(m[r, ])], collapse = "|")
    if(r < nrow) slc <- paste0(slc, "| \n ", collapse = "")
  }
  slc <- paste0("[", slc, "]")
  # diagram
  n_loops <- ceiling(n_slc / ncl_out)
  if(length(diff_time) != n_slc) {
    stop("The length of 'diff_time' doesn't match the expected size.")
  }
  length(diff_time) <- n_loops * ncl_out
  diff_time <- matrix(diff_time, byrow = TRUE, nrow = n_loops, ncol = ncl_out)
  diff_time <- apply(diff_time, 1, max, na.rm = TRUE)
  for(i in seq_len(n_loops)) {
    times <- min(n_slc - (i-1) * ncl_out, ncl_out)
    cat(rep(slc, times))
    free_cores <- cores - times * ncl_in
    if(free_cores > 0) {
      cat(paste0(" ", paste0(rep("_", free_cores), collapse = "|")))
    }
    cat(paste0(" ", round(diff_time[i], 2), "s"))
    cat("\n\n")
  }
  
}

jm_sliced_data <- function (Surv_object, Mixed_objects, time_var, recurrent = FALSE,
                            functional_forms = NULL, data_Surv = NULL, id_var = NULL,
                            priors = NULL, control = NULL, 
                            print = FALSE, cores = NULL, ...) {
  args <- as.list(match.call())[-c(1L, 2L, 3L)] # remove $"", Surv_object, Mixed_objects
  args["print"] <- NULL
  if(is.null(cores)) cores <- max(parallel::detectCores() - 1, 1L)
  if(!exists("n_chains")) n_chains <- 3L
  n_slc <- length(Surv_object) # number of slices
  ncl_in  <- min(n_chains, cores) # number of inner clusters (per outer cluster)
  ncl_out <- max(floor(cores / ncl_in), 1) # number of outer clusters
  args["cores"] <- ncl_in
  outer <- function(i){
    tic2 <- Sys.time()
    jm_fit <- do.call(JMbayes2::jm, c(list(Surv_object = Surv_object[[i]], 
                                           Mixed_objects = lapply(Mixed_objects, "[[", i)), 
                                      args))
    timer <- difftime(Sys.time(), tic2)
    list(mcmc = jm_fit$mcmc)
  }
  tic1 <- Sys.time()
  cl_out <- parallel::makeCluster(ncl_out)
  invisible(parallel::clusterEvalQ(cl_out, library(JMbayes2)))
  parallel::clusterExport(cl_out, list("args", "Surv_object", "Mixed_objects"), env = environment())
  res <- parallel::parLapply(cl_out, seq_len(n_slc), outer)
  parallel::stopCluster(cl_out)
  if(print) { par_diag(n_slc, dots$n_chains, ncl_in, ncl_out, cores, diff_time) }
  res
}