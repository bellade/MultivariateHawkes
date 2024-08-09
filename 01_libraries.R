
library("data.table")
library("ggplot2")
library("reshape2")
library("rstan")

# for plots
library("xtable")
library("latex2exp")
options(mc.cores = parallel::detectCores()) # allows stan to use multiple cores
rstan_options(auto_write = TRUE) # avoid unnecessary recompilations 