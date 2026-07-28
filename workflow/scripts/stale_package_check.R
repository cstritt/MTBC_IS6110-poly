

install.packages("funchir")
funchir::stale_package_check("scripts/collect_results.R")

scripts <- list.files("scripts", pattern = "\\.R$", full.names = TRUE)
lapply(scripts, funchir::stale_package_check)