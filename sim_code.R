pcks <- c("parallel", "JMbayes2")
invisible(lapply(pcks, require, character.only = TRUE))
source("sim_funs.R")
# scenario details
nslices_scn <- c(1L, 2L, 5L, 10L) # number of data slices, include always 1L (no slicing)
n_scn <- c(500L, 1000L, 2500L, 5000L) # number of patients per dataset
n_data <- 200L
# sim details
seed <- 2022
n_cores <- max(parallel::detectCores() - 1, 1L)
#n_cores <- 7
# generate data
tic1 <- Sys.time()
for(n in n_scn) {
    # generate data (in parallel)
    fun <- function(i) {
      source("sim_funs.R")
      gen_data(seed = seed + i, n = n)
    }
    cl <- parallel::makeCluster(n_cores)
    parallel::clusterExport(cl, list("seed", "n")) # load vars in each cl
    dataL <- parallel::parLapply(cl, seq_len(n_data), fun)
    parallel::stopCluster(cl)
    # slice data
    for(nslices in nslices_scn) {
      dataL_surv <- lapply(dataL, "[[", "surv")
      dataL_long <- lapply(dataL, "[[", "long")
      dataL_slc <- mapply(slicer, n_slices = nslices, id_var = "id", 
                          data_Long = dataL_long, data_Surv = dataL_surv, 
                          SIMPLIFY = FALSE)
      saveRDS(dataL_slc, file = paste0("dataL_n", n, "_slc", nslices, ".rds"))
      scn_count <- length(nslices_scn) * (match(n, n_scn) - 1) + match(nslices, nslices_scn)
      cat(paste0("\r Data scenarios generated: ", scn_count, "/", length(n_scn) * length(nslices_scn),
                 " (", Sys.time(),", ", round(difftime(Sys.time(), tic1, units = "mins"), 2), " mins)"))
    }
}
toc1 <- Sys.time()
dur_min1 <- difftime(toc1, tic1, units = "min")
beepr::beep(8)
print(round(difftime(toc1, tic1), 2))
# fit data
tic2 <- Sys.time()
for(n in n_scn) {
    for(nslices in nslices_scn) {
      dataL_slc <- readRDS(file = paste0("dataL_n", n, "_slc", nslices, ".rds"))
      # fit lme
      dataL_slc_long <- lapply(dataL_slc, "[[", "Long")
      lme_all <- lapply(dataL_slc_long, function(data_i) {
        lme_sliced_data(data = data_i,
                        fixed = y ~ splines::ns(time, k =  c(1, 3), B = c(0, 7)),
                        random = list(id = pdDiag(form = ~ splines::ns(time, k = c(1, 3),
                                                                       B = c(0, 7)))),
                        control = list(lmeControl(opt = "optim", niterEM = 45)),
                        cores = n_cores)
      })
      remove(dataL_slc_long)
      # fit coxph
      dataL_slc_surv <- lapply(dataL_slc, "[[", "Surv")
      remove(dataL_slc)
      coxph_all <- lapply(dataL_slc_surv, function(data_i) {
        coxph_sliced_data(data = data_i,
                          formula = Surv(time, Tstatus) ~ group + age,
                          cores = n_cores)
      })
      remove(dataL_slc_surv)
      # fit jmjm_all[[1]]$
      jm_all <- lapply(seq_along(coxph_all), function(i) {
        tic3 <- Sys.time()
        jm_fit <- try(jm_sliced_data(Surv_object = coxph_all[[i]],
                            Mixed_objects = list(lme_all[[i]]),
                            time_var = "time", cores = n_cores), TRUE)
        toc3 <- Sys.time()
        dur_min3 <- difftime(toc3, tic3, units = "min")
        if (inherits(jm_fit, "try-error")) {jm_fit <- dur_min3 <- NA}
        return(list(fit = jm_fit, time = dur_min3))
        })
      remove(coxph_all, lme_all)
      jm_fit  <- lapply(jm_all, "[[", "fit")
      jm_time <- lapply(jm_all, "[[", "time")
      remove(jm_all)
      saveRDS(list(fit = jm_fit, time = jm_time), file = paste0("jm_all_n", n, "_slc", nslices, ".rds"))
      remove(jm_fit, jm_time)
      scn_count <- length(nslices_scn) * (match(n, n_scn) - 1) + match(nslices, nslices_scn)
      cat(paste0("\r Data scenarios fitted: ", scn_count, "/", length(n_scn) * length(nslices_scn),
                 ", n=", n, ", nslices=", nslices, " (", Sys.time(), ", ", 
                 round(difftime(Sys.time(), tic2, units = "mins"), 2), " mins)"))
    }
}
toc2 <- Sys.time()
dur_min2 <- difftime(toc2, tic2, units = "min")
beepr::beep(8)
print(round(difftime(toc2, tic2), 2))