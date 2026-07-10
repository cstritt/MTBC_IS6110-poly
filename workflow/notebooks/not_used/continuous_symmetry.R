# continuous_symmetry.R
# No external packages required (uses base R besselI).
# Weighted von Mises KDE + L2 reflection statistic

# Convert genome positions (bp) to radians
angles_from_positions <- function(positions, genome_length) {
  (2 * pi * (positions %% genome_length) / genome_length) %% (2*pi)
}

# Reflect angles across axis phi
reflect_angles <- function(angles, phi) {
  (2*phi - angles) %% (2*pi)
}

# trapz numerical integral
trapz <- function(x, y) {
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

# von Mises kernel value (vectorized)
vm_kernel <- function(x, mu, kappa) {
  # x and mu are in radians; return matrix-like broadcasting via outer
  # but we'll call it with scalars/vectors via vectorized sum in kde function
  exp(kappa * cos(x - mu)) / (2 * pi * besselI(kappa, 0))
}

# Weighted KDE on grid using von Mises kernels
weighted_vm_kde <- function(angles, weights = NULL, kappa = 10, grid_size = 1024) {
  if (is.null(weights)) weights <- rep(1, length(angles))
  weights <- weights / sum(weights)        # normalize so densities integrate to 1
  grid <- seq(0, 2*pi, length.out = grid_size)
  # For each grid point g, density = sum_i w_i * vm(g, angles[i])
  dens <- vapply(grid, function(g) sum(weights * exp(kappa * cos(g - angles)) / (2*pi*besselI(kappa,0))), numeric(1))
  dens <- dens / trapz(grid, dens)
  list(theta = grid, density = dens)
}

# L2 distance between two density vectors on same grid
l2_distance <- function(f1, f2, grid) {
  dx <- mean(diff(grid))
  sum((f1 - f2)^2) * dx
}

# Compute symmetry statistic for a single axis phi
symmetry_stat_phi <- function(angles, weights = NULL, phi, kappa = 10, grid_size = 1024) {
  kde1 <- weighted_vm_kde(angles, weights, kappa = kappa, grid_size = grid_size)
  reflected <- reflect_angles(angles, phi)
  kde2 <- weighted_vm_kde(reflected, weights, kappa = kappa, grid_size = grid_size)
  l2_distance(kde1$density, kde2$density, kde1$theta)
}

# Scan candidate axes and return full vector of stats + best
find_best_axis_continuous <- function(angles, weights = NULL, n_phi = 360, kappa = 10, grid_size = 1024) {
  phis <- seq(0, 2*pi, length.out = n_phi)
  stats <- vapply(phis, function(p) symmetry_stat_phi(angles, weights, p, kappa = kappa, grid_size = grid_size), numeric(1))
  idx <- which.min(stats)
  list(best_phi = phis[idx], best_stat = stats[idx], phis = phis, stats = stats)
}

plot_kde_and_reflection <- function(angles, weights = NULL, phi = 0, kappa = 10, grid_size = 1024) {
  kde1 <- weighted_vm_kde(angles, weights, kappa = kappa, grid_size = grid_size)
  kde2 <- weighted_vm_kde(reflect_angles(angles, phi), weights, kappa = kappa, grid_size = grid_size)
  plot(kde1$theta, kde1$density, type = "l", xlab = "angle (rad)", ylab = "density", main = paste0("KDE and reflected KDE (phi=", round(phi*180/pi,2),"°)"))
  lines(kde2$theta, kde2$density, col = "red", lty = 2)
  legend("topright", legend = c("original KDE", "reflected KDE"), col = c("black", "red"), lty = c(1,2))
}

n_phi = 360
kappa = 10
grid_size = 1024

obs <- find_best_axis_continuous(angles, weights, n_phi = n_phi, kappa = kappa, grid_size = grid_size)
obs_stat <- obs$best_stat

plot_kde_and_reflection(angles, weights=NULL, phi=obs$best_phi)


