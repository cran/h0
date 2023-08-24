
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



lpost.H0.temp <- function(H0, z.d, z.s, om, kappa,
                        fermat.diff.mean, fermat.diff.se, 
                        time.delay.mean, time.delay.se) {

  s <- 299792.458
  arcsec2radian <- pi / 180 / 3600
  Mpc2km <- 3.085677581e19
  day2s <- 24 * 3600

  td <- TD.dist(H0, z.d, z.s, s, om)
  post.mean <- fermat.diff.mean * arcsec2radian^2 * td * Mpc2km / s / (1 - kappa) / day2s
  post.var <- fermat.diff.se^2 * arcsec2radian^4 * td^2 * Mpc2km^2 / s^2 / (1 - kappa)^2 / day2s^2 + time.delay.se^2
  dt((time.delay.mean - post.mean ) / sqrt(post.var), df = 4, log = TRUE) - log(sqrt(post.var)) 

}



# log posteior density
lpost.H0 <- function(param, z.d, z.s, 
                        fermat.diff.mean, fermat.diff.se,
                        time.delay.mean, time.delay.se) {
  n.obs <- length(z.d)
  kappa <- param[3 : length(param)]
  kappa.rep <- rep(kappa, times = rle(z.s)$length)

  sum(sapply(1 : n.obs, function(k) {
      lpost.H0.temp(H0 = param[1], z.d = z.d[k], z.s = z.s[k], om = param[2],
                       kappa = kappa.rep[k], fermat.diff.mean[k], fermat.diff.se[k], 
                       time.delay.mean[k], time.delay.se[k])})) + 
  sum(dcauchy(kappa, scale = 0.025, log = TRUE))
}



# main function to fit the data
h0 <- function(TD.est, TD.se, 
               FPD.est, FPD.se, 
               z.d, z.s, 
               multimodal = FALSE,
               h0.bound = c(0, 150), h0.scale = 10,
               omega.bound = c(0.05, 0.5),
               initial.param,
               sample.size, burnin.size) {

  epsilon = 1e-308 # for multimodal case
  n.total <- sample.size + burnin.size
  h0.t <- initial.param[1]
  omega.t <- initial.param[2]
  kappa.t <- initial.param[3 : length(initial.param)]
  n.kappa <- length(kappa.t)

  h0.out <- rep(NA, n.total)
  omega.out <- rep(NA, n.total)
  kappa.out <- matrix(NA, nrow = n.total, ncol = length(kappa.t))

  h0.accept <- rep(0, n.total)
  omega.accept <- rep(0, n.total)
  kappa.accept <- matrix(0, nrow = n.total, ncol = length(kappa.t))


  # starting density evaluation
  lpost.t <- lpost.H0(param = c(h0.t, omega.t, kappa.t), 
                              z.d = z.d, z.s = z.s,
                              fermat.diff.mean = FPD.est, fermat.diff.se = FPD.se, 
                              time.delay.mean = TD.est, time.delay.se = TD.se)

  # starting density evaluation for the multimodal sampling
  if (multimodal == TRUE) {
    z.t <- rnorm(1, h0.t, 0.01)
    lpost.z.t <- lpost.H0(param = c(z.t, omega.t, kappa.t), 
                                  z.d = z.d, z.s = z.s,
                                  fermat.diff.mean = FPD.est, fermat.diff.se = FPD.se, 
                                  time.delay.mean = TD.est, time.delay.se = TD.se)
  }

  for (i in 1 : n.total) {
      
    # h0 update
    if (multimodal == TRUE) {
      # initial downhill h0 update
      h0.d <- -1
      while (h0.d < h0.bound[1] | h0.d > h0.bound[2]) {
        h0.d <- rnorm(1, mean = h0.t, sd = h0.scale)
      }
      lpost.h0.d <- lpost.H0(param = c(h0.d, omega.t, kappa.t), 
                                   z.d = z.d, z.s = z.s,
                                   fermat.diff.mean = FPD.est, fermat.diff.se = FPD.se, 
                                   time.delay.mean = TD.est, time.delay.se = TD.se)
      # repeat until the downhill move is accepted
      while (-rexp(1) > log(exp(lpost.t) + epsilon) - log(exp(lpost.h0.d) + epsilon)) {
        h0.d <- -1
        while (h0.d < h0.bound[1] | h0.d > h0.bound[2]) {
          h0.d <- rnorm(1, mean = h0.t, sd = h0.scale)
        }
        lpost.h0.d <- lpost.H0(param = c(h0.d, omega.t, kappa.t), 
                                     z.d = z.d, z.s = z.s,
                                     fermat.diff.mean = FPD.est, fermat.diff.se = FPD.se, 
                                     time.delay.mean = TD.est, time.delay.se = TD.se)
      }

      # initial uphill h0 update
      h0.u <- -1
      while (h0.u < h0.bound[1] | h0.u > h0.bound[2]) {
        h0.u <- rnorm(1, mean = h0.d, sd = h0.scale)
      }
      lpost.h0.u <- lpost.H0(param = c(h0.u, omega.t, kappa.t), 
                                   z.d = z.d, z.s = z.s,
                                   fermat.diff.mean = FPD.est, fermat.diff.se = FPD.se, 
                                   time.delay.mean = TD.est, time.delay.se = TD.se)
      # repeat until the uphill move is accepted
      while (-rexp(1) > log(exp(lpost.h0.u) + epsilon) - log(exp(lpost.h0.d) + epsilon)) {
        h0.u <- -1
        while (h0.u < h0.bound[1] | h0.u > h0.bound[2]) {
          h0.u <- rnorm(1, mean = h0.d, sd = h0.scale)
        }
        lpost.h0.u <- lpost.H0(param = c(h0.u, omega.t, kappa.t), 
                                     z.d = z.d, z.s = z.s,
                                     fermat.diff.mean = FPD.est, fermat.diff.se = FPD.se, 
                                     time.delay.mean = TD.est, time.delay.se = TD.se)
      }

      # initial downhill auxiliary h0 update
      h0.z.d <- -1
      while (h0.z.d < h0.bound[1] | h0.z.d > h0.bound[2]) {
        h0.z.d <- rnorm(1, mean = h0.u, sd = h0.scale)
      }
      lpost.z.d <- lpost.H0(param = c(h0.z.d, omega.t, kappa.t), 
                                  z.d = z.d, z.s = z.s,
                                  fermat.diff.mean = FPD.est, fermat.diff.se = FPD.se, 
                                  time.delay.mean = TD.est, time.delay.se = TD.se)
      # repeat until the downhill auxiliary move is accepted
      while (-rexp(1) > log(exp(lpost.h0.u) + epsilon) - log(exp(lpost.z.d) + epsilon)) {
        h0.z.d <- -1
        while (h0.z.d < h0.bound[1] | h0.z.d > h0.bound[2]) {
          h0.z.d <- rnorm(1, mean = h0.u, sd = h0.scale)
        }
        lpost.z.d <- lpost.H0(param = c(h0.z.d, omega.t, kappa.t), 
                                    z.d = z.d, z.s = z.s,
                                    fermat.diff.mean = FPD.est, fermat.diff.se = FPD.se, 
                                    time.delay.mean = TD.est, time.delay.se = TD.se)
      }

      # accept or reject the proposal
      min.nu <- min(1, (exp(lpost.t) + epsilon) / (exp(lpost.z.t) + epsilon))
      min.de <- min(1, (exp(lpost.h0.u) + epsilon) / (exp(lpost.z.d) + epsilon))
      l.mh <- lpost.h0.u - lpost.t + log(min.nu) - log(min.de)

      if (l.mh > -rexp(1)) {
        h0.t <- h0.u
        z.t <- h0.z.d
        lpost.t <- lpost.h0.u
        lpost.z.t <- lpost.z.d
        h0.accept[i] <- 1
      }

    } else { 
    # for default MCMC (non-multimodal)

      h0.p <- -1
      while (h0.p < h0.bound[1] | h0.p > h0.bound[2]) {
        h0.p <- rnorm(1, mean = h0.t, sd = h0.scale)
      }
      lpost.p <- lpost.H0(param = c(h0.p, omega.t, kappa.t), 
                                z.d = z.d, z.s = z.s,
                                fermat.diff.mean = FPD.est, fermat.diff.se = FPD.se, 
                                time.delay.mean = TD.est, time.delay.se = TD.se)
      log.metropolis <- lpost.p - lpost.t

      if (log.metropolis > -rexp(1)) {
        h0.t <- h0.p
        h0.accept[i] <- 1
        lpost.t <- lpost.p
      }
    }

    h0.out[i] <- h0.t

    # adaptive MCMC 10% acceptence rate for multimodal and 40% for default
    if (multimodal == TRUE) {
      if (i %% 100 == 0) {
        if(mean(h0.accept[(i - 99) : i]) > 0.1) {
          scale.adj <- exp(min(0.1, 1 / sqrt(i / 100)))
        } else if (mean(h0.accept[(i - 99) : i]) < 0.1) {
          scale.adj <- exp(-min(0.1, 1 / sqrt(i / 100)))
        }
        h0.scale <- h0.scale * scale.adj
      }
    } else {
      if (i %% 100 == 0) {
        if(mean(h0.accept[(i - 99) : i]) > 0.4) {
           scale.adj <- exp(min(0.1, 1 / sqrt(i / 100)))
        } else if (mean(h0.accept[(i - 99) : i]) < 0.4) {
           scale.adj <- exp(-min(0.1, 1 / sqrt(i / 100)))
        }
        h0.scale <- h0.scale * scale.adj
      }
    }

    # omega update
    omega.p <- runif(1, min = omega.bound[1], max = omega.bound[2])
    lpost.p <- lpost.H0(param = c(h0.t, omega.p, kappa.t), 
                              z.d = z.d, z.s = z.s,
                              fermat.diff.mean = FPD.est, fermat.diff.se = FPD.se, 
                              time.delay.mean = TD.est, time.delay.se = TD.se)
    log.metropolis <- lpost.p - lpost.t

    if (log.metropolis > -rexp(1)) {
      omega.t <- omega.p
      omega.accept[i] <- 1
      lpost.t <- lpost.p
    }

    omega.out[i] <- omega.t

    # kappa updates
    for (j in 1 : n.kappa) {
      kappa.p <- kappa.t
      kappa.p[j] <- rcauchy(1, scale = 0.025)
      
      lpost.p <- lpost.H0(param = c(h0.t, omega.t, kappa.p), 
                                z.d = z.d, z.s = z.s,
                                fermat.diff.mean = FPD.est, fermat.diff.se = FPD.se, 
                                time.delay.mean = TD.est, time.delay.se = TD.se)
      log.metropolis <- lpost.p - lpost.t
      log.hastings <-  sum(dcauchy(kappa.t, scale = 0.025, log = TRUE)) - 
                       sum(dcauchy(kappa.p, scale = 0.025, log = TRUE))

      if ((log.metropolis + log.hastings) > -rexp(1)) {
        kappa.t <- kappa.p
        kappa.accept[i, j] <- 1
        lpost.t <- lpost.p
      }
       
    } # end of updating kappa
  
    kappa.out[i, ] <- kappa.t

  } # end of iterations

  # saving outcomes into output
  if (length(kappa.t) == 1) {
    output <- list(h0 = h0.out[-c(1 : burnin.size)], 
                   omega = omega.out[-c(1 : burnin.size)], 
                   kappa = kappa.out[-c(1 : burnin.size), ],
                   h0.accept.rate = mean(h0.accept[-c(1 : burnin.size)]),
                   omega.accept.rate = mean(omega.accept[-c(1 : burnin.size)]),
                   kappa.accept.rate = mean(kappa.accept[-c(1 : burnin.size)]),
                   h0.scale = h0.scale)

  } else {
    output <- list(h0 = h0.out[-c(1 : burnin.size)], 
                   omega = omega.out[-c(1 : burnin.size)], 
                   kappa = kappa.out[-c(1 : burnin.size), ],
                   h0.accept.rate = mean(h0.accept[-c(1 : burnin.size)]),
                   omega.accept.rate = mean(omega.accept[-c(1 : burnin.size)]),
                   kappa.accept.rate = colMeans(kappa.accept[-c(1 : burnin.size), ]),
                   h0.scale = h0.scale)
  } # end of saving

  # reporting output
  output
}
