\name{h0}
\alias{h0}
\title{Obtaining Bayesian posterior sample or maximum likelihood estimate of the Hubble constant via a meta analysis}
\usage{
h0(TD.est, TD.se, 
   FPD.est, FPD.se, 
   z.d, z.s,
   h0.bound = c(60, 80),
   h0.grid.size = 400,
   omega.bound = c(0.05, 0.5),
   omega.grid.size = 200,
   sample.size = 1e5, 
   method = "mle")
}
\arguments{
\item{TD.est}{A vector for time delay estimates. For example, if a data set contained three time delay estimates of a quad-lens system and one time delay estimate of a double-lens system, then this vector would be composed of the three time delay estimates from the quad and one time delay estimate from the double.}

\item{TD.se}{A vector for standard errors (one-sigma uncertainties) of the time delay estimates. For example, if a data set contained three time delay estimates of a quad-lens system and one time delay estimate of a double-lens system, then this vector would be composed of the three standard errors corresponding to the first three time delay estimates from the quad and one standard error corresponding to the single time delay estimate from the double.}

\item{FPD.est}{A vector for Fermat potential difference estimates. For example, if a data set contained three Fermat potential difference estimates of a quad-lens system and one Fermat potential difference estimate of a double-lens system, then this vector would be composed of the three Fermat potential difference estimates from the quad and one Fermat potential difference estimate from the double.}

\item{FPD.se}{A vector for standard errors (one-sigma uncertainties) of the Fermat potential difference estimates. For example, if a data set contained three Fermat potential difference estimates of a quad-lens system and one Fermat potential difference estimate of a double-lens system, then this vector would be composed of the three standard errors corresponding to the first three Fermat potential difference estimates from the quad and one standard error corresponding to the single Fermat potential difference estimate from the double.}

\item{z.d}{A vector for redshifts of deflectors (lenses). For example, if a data set contained three time delay estimates of a quad-lens system and one time delay estimate of a double-lens system, then users need to fill out the first three spots of this vector by the same redshift of the quad-lens, and fill out the last spot by the redshift of the double-lens.}

\item{z.s}{A vector for redshifts of sources (quasars). For example, if a data set contained three time delay estimates of a quad-lens system and one time delay estimate of a double-lens system, then users need to fill out the first three spots of this vector by the same redshift of the quadruply lensed quasar, and fill out the last spot by the redshift of the doubly-lensed quasar.}

\item{h0.bound}{It determins the search range for the Hubble constant for both maximum likelihood estimation and posterior sampling. Default is to search possible values of the Hubble constant between 60 and 80. Users can set a wider range, if needed.}

\item{h0.grid.size}{Optional if a Bayesian inference is to be conducted with \code{method = "bayes"}. It determins the number of values on a grid of the Hubble constant. Default is to set it to 400, meaning that the grid sampling method will start on a grid of 400 values of the Hubble constant within the range (set by \code{h0.bound}). For more accurate grid sampling, it is recommended to increase this size when the search region increases.}

\item{omega.bound}{It determins the search range for the current dark matter density for both maximum likelihood estimation and posterior sampling. Default is to search possible values of the current dark matter density between 0.05 and 0.5. Users can set a wider range, if needed.}

\item{omega.grid.size}{Optional if a Bayesian inference is to be conducted with \code{method = "bayes"}. It determins the number of values on a grid of the current dark matter density. Default is to set it to 200, meaning that the grid sampling method will start on a grid of 200 values of the current dark matter density within the range (set by \code{omega.bound}).  For more accurate grid sampling, it is recommended to increase this size when the search region increases.}

\item{sample.size}{Optional if a Bayesian inference is to be conducted with \code{method = "bayes"}. It determines the posterior sample size for each of the two model parameters, the Hubble constant and the current dark matter density, to be drawn by a grid sampling method. Its default value is 100,000.}

\item{method}{Maximum likelihood estimation is conducted if \code{method = "mle"}, while Bayesian posterior sampling is done if \code{method = "bayes"}.}

}
\value{
The outcomes of \code{h0} are composed of:
\describe{
\item{h0.mle}{A constant. The maximum likelihood estimate of the Hubble constant, available when \code{method = "mle"}.}

\item{h0.se}{A constant. An asymptotic one-sigma uncertainty (standard error) of the maximum likelihood estimate of the Hubble constant, available when \code{method = "mle"}. This is computed from the curvature of the likelihood function evaluated at the maximum likelihood estimate.}

\item{omega.mle}{A constant. The maximum likelihood estimate of the current dark matter density, available when \code{method = "mle"}.}

\item{omega.se}{A constant. An asymptotic one-sigma uncertainty (standard error) of the current dark matter density estimate of the Hubble constant, available when \code{method = "mle"}. This is computed from the curvature of the likelihood function evaluated at the maximum likelihood estimate.}

\item{h0}{A vector. The posterior sample of the Hubble constant, available when \code{method = "bayes"}. Its sample size is set by users.}

\item{omega}{A vector. The posterior sample of the current dark matter density, available when \code{method = "bayes"}. Its sample size is set by users.}

\item{contour}{A matrix. All of the likelihood values computed on the two-dimensional grid, available when \code{method = "bayes"}. This matrix can be used to produce a contour plot}

\item{h0.grid}{A vector. It contains the grid of values of the Hubble constant that has been used during the grid sampling procedure, available when \code{method = "bayes"}. This can also be used to draw a contour plot.}

\item{omega.grid}{A vector. It contains the grid of values of the current dark matter density that has been used during the grid sampling procedure, available when \code{method = "bayes"}. This can also be used to draw a contour plot.}
}
}
\description{
The function \code{h0} conducts a meta analysis to infer the current expansion rate of the University (Hubble constant) via maximum likelihood estimation or Bayesian posterior sampling. The analysis is constructed under time delay cosmography. The minimum inputs to implement this function are (1) time delay estimates (\code{TD.est}); (2) their standard errors (\code{TD.se}); (3) Fermat potential difference estimates (\code{FPD.est}); (4) their standard errors (\code{FPD.se}); (5) redshifts of deflectors (lenses) (\code{z.d}); and (6) redshifts of sources (quasars) (\code{z.s}). We note that the length of each input vector must be identical, and the order of input values in each vector must be consistent across all of the vectors. For example, let's say the third input value in one of the input vectors, \code{TD.est}, is a time delay estimate between gravitationally lensed images A and C of a specific quasar. Then the third input value in the other input vectors must be a quantity corresponding to the gravitationally lensed images A and C of that specific quasar.
}
\details{
Time delay cosmography infers the Hubble constant \eqn{H_0} by modeling time delays caused by strong gravitational lensing. The following equation is the basis of time delay cosmography:
  \deqn{c\Delta_{ijk}=D_{\Delta}(H_0, z_k, \Omega) \phi_{ijk},}
where \eqn{c} denotes the speed of light, and \eqn{\Delta_{ijk}} is the time delay in days between lensed images \eqn{i} and \eqn{j} of quasar \eqn{k} (\eqn{k=1, 2, \ldots, n}). The notation \eqn{z_k=\{z_{sk}, z_{dk}\}} indicates a vector for two redshifts of the \eqn{k}-th source (\eqn{z_{sk}}) and deflector (\eqn{z_{dk}}). The next notation \eqn{\Omega=\{\Omega_{m}, \Omega_{\Lambda}\}} indicates a vector for the present-day dark matter density \eqn{\Omega_{m}} and dark energy density \eqn{\Omega_{\Lambda}}. Since their sum becomes one under the flat \eqn{\Lambda}CDM model (standard cosmology),  we treat \eqn{\Omega_{m}} as the only unknown parameter. The notation  \eqn{D_{\Delta}(H_0, z_k, \Omega)} denotes the time delay distance in the unit of megaparsec (Mpc), which is a deterministic function of \eqn{H_0}, \eqn{z_k}, and \eqn{\Omega}.

We emphasize that the physical equation  deterministically relates five quantities, i.e., time delay \eqn{\Delta_{ijk}}, the Hubble constant \eqn{H_0}, redshifts \eqn{z_k}, cosmological parameter \eqn{\Omega}, and Fermat potential difference \eqn{\phi_{ijk}}. That means, if we completely knew  four quantities,  \eqn{\Delta_{ijk}}, \eqn{z_k}, \eqn{\Omega}, and \eqn{\phi_{ijk}}, then \eqn{H_0} can be easily determined to be a single number. The meta analysis takes advantage of this relationship, marginalizes nuisance parameters, and returns estimates of \eqn{H_0} and \eqn{\Omega} via maximum likelihood estimation or Bayesian posterior sampling.
}
\examples{

#####  Suppose we have three time delay estimates of a quad-lens system, 
#####  and one time delay estimate of a double-lens system. 
time.delay.est <- c(-101, -91.7, -104.2, -18.8)

#####  Each time delay estimate accompanies its one-sigma uncertainty (standard error)
time.delay.se <- c(4.4, 4.2, 4.1, 0.9)

#####  The redshift of the deflector in the quad-lens system is 0.868, 
#####  and that in the double-lens system is 0.512.
z.d.obs <- c(0.868, 0.868, 0.868, 0.512)

#####  The redshift of the source in the quad-lens system is 2.621, 
#####  and that in the double-lens system is 3.37.
z.s.obs <- c(2.621, 2.621, 2.621, 3.37)

#####  Three Fermat potential difference estimates are obtained from a quad-lens syatem, 
#####  and one estimate is from the double-lens system.
fermat.diff.est <- c(-0.6395, -0.5734, -0.6431, -0.2625)

#####  Each Fermat potential difference estimate accompanies 
#####  its one-sigma uncertainty (standard error).
fermat.diff.se <- c(0.0265, 0.0254, 0.0299, 0.0148)

#### Next, we make sure that each vector has the same length.
length(z.d.obs); length(z.s.obs); length(fermat.diff.est); 
length(fermat.diff.se); length(time.delay.est); length(time.delay.se)

#####  The following is an example code to obtain 1000 posterior samples of the Hubble constant
#####  and those of the current dark matter density.
#####  In practice, users may want to increase the grid sizes for a more accuracte inference, 
#####  even though it increases computational burden.
res <- h0(TD.est = time.delay.est, TD.se = time.delay.se, 
          FPD.est = fermat.diff.est, FPD.se = fermat.diff.se,
          z.d = z.d.obs, z.s = z.s.obs, 
          h0.bound = c(60, 80), h0.grid.size = 50,
          omega.bound = c(0.05, 0.5), omega.grid.size = 25, 
          sample.size = 1000, method = "bayes")

#####  We can draw a contour plot, and the following two histograms, 
#####  using the output from the Bayesian fit. 
contour(res$h0.grid, res$omega.grid, res$contour, 
        nlevels = 5, drawlabels = FALSE,
        ylim = c(0.05, 0.5), xlim = c(60, 80), levels = c(0.01, 0.25, 0.5, 0.75, 0.99))

hist(res$h0, 20, prob = TRUE, col = 0, xlim = c(60, 80))

hist(res$omega, 30, prob = TRUE, col = 0, xlim = c(0.05, 0.55))

#####  The following is an example code to obtain the maximum likelihood estimate of 
#####  the Hubble constant and its asymptotic uncertainty. 
res.mle <- h0(TD.est = time.delay.est, TD.se = time.delay.se, 
              FPD.est = fermat.diff.est, FPD.se = fermat.diff.se,
              z.d = z.d.obs, z.s = z.s.obs, 
              h0.bound = c(60, 80), omega.bound = c(0.05, 0.5),
              method = "mle")

res.mle$h0.mle; res.mle$h0.se; res.mle$omega.mle; res.mle$omega.se
}
\author{
Hyungsuk Tak
}
\references{
    T. Treu and P. Marshall (2016) <doi:10.1007/s00159-016-0096-8> "Time Delay Cosmography".
}
