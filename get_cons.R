par_nms <- c("betas1", "gammas", "alphas")
nslices_scn <- c(1L, 2L, 5L, 10L) # number of data slices, include always 1L (no slicing)
n_scn <- c(500L, 1000L, 2500L, 5000L) # number of patients per dataset
cons_scn <- c("union", "eq_wt", "vr_wt")
# make consensus
cons <- list()
dim <- sapply(list(n_scn, nslices_scn), length)
n_cons <- length(cons_scn)
tic1 <- Sys.time()
for(par in par_nms) { # par <- par_nms[1]
  jm_all <- readRDS(file = paste0("jm_all_n", n_scn[1], "_slc", nslices_scn[1], ".rds"))
  if(par == par_nms[1]) {
    n_data <- length(jm_all$fit)
    n_chain <- length(jm_all$fit[[1]][[1]]$mcmc[[par]]) # jm_all$fit[[data_i]][[slice_i]]
    n_iter <- nrow(jm_all$fit[[1]][[1]]$mcmc[[par]][[1]]) * n_chain
    cons[["time"]] <- array(NA, dim = c(dim, n_data))
  }
  n_par <- ncol(jm_all$fit[[1]][[1]]$mcmc[[par]][[1]])
  cons[["mcmc"]][[par]] <- array(NA, 
                                 dim = c(dim, n_par, n_iter * max(nslices_scn), n_data, n_cons)) # array dimensions, n_iter * max(nslices_scn) because for the union we are putting all together
  cons[["lw"]][[par]] <- cons[["up"]][[par]] <- cons[["mn"]][[par]] <- 
    array(NA, dim = c(dim, n_par, n_data, n_cons)) # array dimensions
  cons[["wts"]][[par]] <- array(NA, dim = c(dim, n_par, max(nslices_scn), n_data, n_cons))
  for(n in n_scn) { # n <- n_scn[1]
    row <- match(n, n_scn)
    for(nslices in nslices_scn) { # nslices <- nslices_scn[1]
      col <- match(nslices, nslices_scn)
      if(!(row == 1 & col == 1)) { # when (row=1, col=1) jm_all is already loaded above 
        file_nm <- paste0("jm_all_n", n, "_slc", nslices, ".rds")
        if(file.exists(file_nm)){
          jm_all <- readRDS(file = file_nm)
        } else {
          break; print(paste0("File not found in the dir: ", file_nm))
        }
      }
      if(par == par_nms[1]) { # the time is the same for all parameters in the same model
        cons[["time"]][row, col, ] <- unlist(jm_all$time)
      }
      for(data_i in seq_len(n_data)) { # data_i <- seq_len(n_data)[1]
        mcmc_lst <- lapply(seq_len(nslices), function(slc_i) {
          do.call("rbind", jm_all$fit[[data_i]][[slc_i]]$mcmc[[par]])
        })
        mcmc_arr <- array(unlist(mcmc_lst), dim = c(dim(mcmc_lst[[1]]), nslices),
                          dimnames = list(NULL, colnames(mcmc_lst[[1]]), NULL))
        remove(mcmc_lst)
        # union
        if("union" %in% cons_scn) {
          mcmc <- apply(mcmc_arr, 2, rbind)
          cons[["mcmc"]][[par]][row, col, , seq_len(n_iter * nslices), data_i, 1] <-  t(mcmc)
          cons[["lw"]][[par]][row, col, , data_i, 1] <- apply(mcmc, 2, quantile, probs = c(0.025))
          cons[["up"]][[par]][row, col, , data_i, 1] <- apply(mcmc, 2, quantile, probs = c(0.975))
          cons[["mn"]][[par]][row, col, , data_i, 1] <- apply(mcmc, 2, mean)
          remove(mcmc)
        }
        # equal-weight
        if("eq_wt" %in% cons_scn) {
          eq_w <- rep(1/nslices, nslices)
          cons[["wts"]][[par]][row, col, , seq_len(nslices), data_i, 2] <- rep(eq_w, each = n_par)
          eq_W <- array(rep(eq_w, each = dim(mcmc_arr)[1]), dim = dim(mcmc_arr))
          mcmc <- apply(mcmc_arr * eq_W, c(1:2), sum)
          cons[["mcmc"]][[par]][row, col, , seq_len(n_iter), data_i, 2] <- t(mcmc)
          cons[["lw"]][[par]][row, col, , data_i, 2] <- apply(mcmc, 2, quantile, probs = c(0.025))
          cons[["up"]][[par]][row, col, , data_i, 2] <- apply(mcmc, 2, quantile, probs = c(0.975))
          cons[["mn"]][[par]][row, col, , data_i, 2] <- apply(mcmc, 2, mean)
          remove(mcmc, eq_w, eq_W)
        }
        # var-weight
        if("vr_wt" %in% cons_scn) {
          vr_w <- 1/apply(mcmc_arr, 3, matrixStats::colVars)
          if(is.vector(vr_w)) vr_w <- t(as.matrix(vr_w)) #?? probably there is smarter way to do this
          vr_w <- vr_w/rowSums(vr_w)
          cons[["wts"]][[par]][row, col, , seq_len(nslices), data_i, 3] <- vr_w
          vr_W <- array(rep(vr_w, each = dim(mcmc_arr)[1]), dim = dim(mcmc_arr))
          mcmc <- apply(mcmc_arr * vr_W, c(1:2), sum)
          cons[["mcmc"]][[par]][row, col, , seq_len(n_iter), data_i, 3] <- t(mcmc)
          cons[["lw"]][[par]][row, col, , data_i, 3] <- apply(mcmc, 2, quantile, probs = c(0.025))
          cons[["up"]][[par]][row, col, , data_i, 3] <- apply(mcmc, 2, quantile, probs = c(0.975))
          cons[["mn"]][[par]][row, col, , data_i, 3] <- apply(mcmc, 2, mean)
          remove(mcmc, vr_w, vr_W)
        }
        remove(mcmc_arr)
      }
      remove(jm_all)
      scn_count <- length(n_scn) * (match(nslices, nslices_scn) - 1) + match(n, n_scn) 
      cat(paste0("\r Param: ", match(par, par_nms), "/", length(par_nms), ", Data scenarios analysed: ", scn_count, "/", length(n_scn) * length(nslices_scn),
                 ", n=", n, ", nslices=", nslices, " (", Sys.time(), ", ", 
                 round(difftime(Sys.time(), tic1, units = "mins"), 2), " mins)     "))
    }
  }
}
toc1 <- Sys.time()
dur_min1 <- difftime(toc1, tic1, units = "min")
beepr::beep(8)
print(round(difftime(toc1, tic1), 2))
# save results
saveRDS(cons, file = paste0("cons.rds"))
remove(cons)