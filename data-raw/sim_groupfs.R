# Simulation results used in vignettes/articles/selectiveInference.Rmd.
# Run from the package root after installing the package:
#   Rscript data-raw/sim_groupfs.R [niters] [cores]
library(unbiasedgoodness)
args <- commandArgs(trailingOnly = TRUE)
niters <- if (length(args) >= 1) as.integer(args[1]) else 1000L
cores <- if (length(args) >= 2) as.integer(args[2]) else parallel::detectCores()

settings <- list(n = 100, p = 200, s0 = 5, m = 1.2, b = 10,
                 xtype = "toeplitz", btype = "U[-2,2]")
start <- Sys.time()
output <- do.call(simulate_gof,
                  c(list(niters = niters, method = "groupfs", cores = cores,
                         seed = 20221010), settings))
elapsed <- difftime(Sys.time(), start, units = "mins")
message("groupfs simulation: ", nrow(output), " iterations (",
        attr(output, "n_failed"), " failed) in ", round(elapsed, 1), " minutes")
write.csv(output, "inst/extdata/sim_groupfs.csv", row.names = FALSE)
