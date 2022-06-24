
E.inv <- function(z, om) {
  1 / sqrt(om * (1 + z)^3 + (1 - om))
}



D <- function(H0, z, s, om) {
  s * integrate(E.inv, lower = 0, upper = z, om = om)$value / H0 / (1 + z)
}



D.diff <- function(H0, z_d, z_s, s, om) {
  s * {integrate(E.inv, lower = 0, upper = z_s, om = om)$value - 
   integrate(E.inv, lower = 0, upper = z_d, om = om)$value} / H0 / (1 + z_s)
}



TD.dist <- function(H0, z_d, z_s, s, om) {
  D_d <- D(H0, z_d, s, om)
  D_s <- D(H0, z_s, s, om)
  D_ds <- D.diff(H0, z_d, z_s, s, om)
  (1 + z_d) * D_d * D_s / D_ds
}



log.post.H0.temp <- function(H0, z.d, z.s, om,
                        fermat.diff.mean, fermat.diff.se, 
                        time.delay.mean, time.delay.se) {

  s <- 299792.458
  arcsec2radian <- pi / 180 / 3600
  Mpc2km <- 3.085677581e19
  day2s <- 24 * 3600

  td <- TD.dist(H0, z.d, z.s, s, om)
  post.mean <- fermat.diff.mean * arcsec2radian^2 * td * Mpc2km / s / day2s
  post.var <- fermat.diff.se^2 * arcsec2radian^4 * td^2 * Mpc2km^2 / s^2 / day2s^2 + time.delay.se^2
  dnorm(time.delay.mean, mean = post.mean, sd = sqrt(post.var), log = TRUE)
}



log.post.H0 <- function(param, z.d, z.s, 
                        fermat.diff.mean, fermat.diff.se,
                        time.delay.mean, time.delay.se) {
  n.obs <- length(z.d)
  sum(sapply(1 : n.obs, function(k) {
      log.post.H0.temp(H0 = param[1], z.d = z.d[k], z.s = z.s[k], om = param[2],
                  fermat.diff.mean[k], fermat.diff.se[k], 
                  time.delay.mean[k], time.delay.se[k])}))

}



h0 <- function(TD.est, TD.se, 
               FPD.est, FPD.se, 
               z.d, z.s,
               h0.bound = c(60, 80),
               h0.grid.size = 400,
               omega.bound = c(0.05, 0.5),
               omega.grid.size = 200,
               sample.size = 1e5, 
               method = "mle") {

  if (method == "bayes") {

    h0.grid <- seq(h0.bound[1], h0.bound[2], length.out = h0.grid.size)
    omega.grid <- seq(omega.bound[1], omega.bound[2], length.out = omega.grid.size)
    exp.ab <- expand.grid(h0.grid, omega.grid)
    n.eval <- nrow(exp.ab)

    log.post.z <- sapply(1 : n.eval, function(a){
      log.post.H0(param = c(exp.ab[a, 1], exp.ab[a, 2]), z.d = z.d, z.s = z.s, 
                  fermat.diff.mean = FPD.est, fermat.diff.se = FPD.se, 
                  time.delay.mean = TD.est, time.delay.se = TD.se)
    })

    post.z <- matrix(exp(log.post.z - max(log.post.z)), nrow = length(h0.grid), ncol = length(omega.grid))

    h0.marg <- rowSums(post.z)
    h0.index <- sample(1 : h0.grid.size, size = sample.size, replace = TRUE, prob = h0.marg)
    h0.sample <- h0.grid[h0.index]
    h0.sample <- h0.sample + runif(sample.size, -(h0.grid[2] - h0.grid[1]) / 2, (h0.grid[2] - h0.grid[1]) / 2)

    omega.sample <- sapply(h0.index, function(i) { sample(omega.grid, size = 1, replace = TRUE, prob = post.z[i, ])})
    omega.sample <- omega.sample + 
                    runif(sample.size, -(omega.grid[2] - omega.grid[1]) / 2, (omega.grid[2] - omega.grid[1]) / 2)    

    output <- list(h0 = h0.sample, omega = omega.sample, contour = post.z, h0.grid = h0.grid, omega.grid = omega.grid)

  } else {

    res <- optim(c(70, 0.3), log.post.H0, control = list(fnscale = -1), hessian = TRUE, method = "L-BFGS-B",
                 lower = c(h0.bound[1], omega.bound[1]), upper = c(h0.bound[2], omega.bound[2]), z.d = z.d, z.s = z.s, 
                 fermat.diff.mean = FPD.est, fermat.diff.se = FPD.se, 
                 time.delay.mean = TD.est, time.delay.se = TD.se)
    mle <- res$par
    se <- sqrt(diag(-solve(res$hessian)))

    output <- list(h0.mle = mle[1], h0.se = se[1], omega.mle = mle[2], omega.se = se[2])

  }

  output

}

