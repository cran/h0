\name{h0}
\alias{h0}
\title{Obtaining Bayesian posterior sample of the Hubble constant via the proposed meta-analysis}
\usage{
h0(TD.est, TD.se, 
   FPD.est, FPD.se, 
   z.d, z.s,
   multimodal = FALSE,
   h0.bound = c(0, 150),
   h0.scale = 10,
   omega.bound = c(0.05, 0.5),
   initial.param,
   sample.size, burnin.size)
}
\arguments{
\item{TD.est}{A vector for time delay estimates (in days). For example, if a data set contained three time delay estimates of a quad-lens system and one time delay estimate of a double-lens system, then this vector would be composed of the three time delay estimates from the quad and one time delay estimate from the double.}

\item{TD.se}{A vector for standard errors (one-sigma uncertainties) of the time delay estimates. For example, if a data set contained three time delay estimates of a quad-lens system and one time delay estimate of a double-lens system, then this vector would be composed of the three standard errors corresponding to the first three time delay estimates from the quad and one standard error corresponding to the single time delay estimate from the double.}

\item{FPD.est}{A vector for Fermat potential difference estimates. For example, if a data set contained three Fermat potential difference estimates of a quad-lens system and one Fermat potential difference estimate of a double-lens system, then this vector would be composed of the three Fermat potential difference estimates from the quad and one Fermat potential difference estimate from the double.}

\item{FPD.se}{A vector for standard errors (one-sigma uncertainties) of the Fermat potential difference estimates. For example, if a data set contained three Fermat potential difference estimates of a quad-lens system and one Fermat potential difference estimate of a double-lens system, then this vector would be composed of the three standard errors corresponding to the first three Fermat potential difference estimates from the quad and one standard error corresponding to the single Fermat potential difference estimate from the double.}

\item{z.d}{A vector for redshifts of deflectors (lenses). For example, if a data set contained three time delay estimates of a quad-lens system and one time delay estimate of a double-lens system, then users need to fill out the first three spots of this vector by the same redshift of the quad-lens, and fill out the last spot by the redshift of the double-lens.}

\item{z.s}{A vector for redshifts of sources (quasars). For example, if a data set contained three time delay estimates of a quad-lens system and one time delay estimate of a double-lens system, then users need to fill out the first three spots of this vector by the same redshift of the quadruply lensed quasar, and fill out the last spot by the redshift of the doubly-lensed quasar.}

\item{multimodal}{If "TRUE", the sampling method switches from Metropolis-Hastings within Gibbs to repelling-attracting Metropolis-Hastings within Gibbs to enhance the convergence of Markov chain to a multimodal target distribution. Default is "FALSE".}

\item{h0.bound}{It determins the Uniform prior range for the Hubble constant. Default is to search possible values of the Hubble constant between 0 and 150. Users can set a wider range, if needed.}

\item{h0.scale}{It determins the proposal scale of the Metropolis update for the Hubble constant. Default is 10. Due to the adaptive Markov chain Monte Carlo method adopted in this package, this scale changes at every 100 iterations. If "multimodal = FALSE", then the scale changes to achieve the acceptance rate equal to 0.4. If "multimodal = TRUE", the scale is adapted to achieve the acceptance rate equal to 0.1.}

\item{omega.bound}{It determins the prior range for the current dark matter density. Default is to search possible values of the current dark matter density between 0.05 and 0.5. Users can set a wider range, if needed.}

\item{initial.param}{A vector of (K+2) model parameters, that is, the Hubble constant, the dark matter density, and K external convergences, where K is the number of lens systems in the data set. A reasonable set of values is  c(runif(1, 0, 150), runif(1, 0.05, 0.5), rcaucy(K, scale = 0.025)).}


\item{sample.size}{It determines the posterior sample size for each of the model parameters, the Hubble constant, the dark matter density, and K external convergences.}

\item{burnin.size}{It determines the warming-up posterior sample size to be discarded before obtaining the desired sample size. The total number of iterations is the sum of sample.size and burnin.size.}

}
\value{
The outcomes of \code{h0} are composed of:
\describe{

\item{h0}{A vector. The posterior sample of the Hubble constant.}

\item{omega}{A vector. The posterior sample of the dark matter density.}

\item{kappa}{A matrix. Each row contains one posterior sample of K external convergences.}

\item{h0.accept.rate}{A scalar. The acceptance rate of the posterior sample of the Hubble constant.}

\item{omega.accept.rate}{A scalar. The acceptance rate of the posterior sample of the dark matter density.}

\item{kappa.accept.rate}{A vector. The acceptance rates of the posterior samples of K external convergences.}

\item{h0.scale}{A scalar. The proposal scale of the Metropolis update when the last ieration is made because it changes at every 100 iterations for adaptive Markov chain Monte Carlo.}
}
}
\description{
The function \code{h0} conducts a Bayesian meta-analysis to infer the current expansion rate of the University (Hubble constant). The analysis is constructed under time delay cosmography. The minimum inputs of the meta-analysis are (1) time delay estimates (\code{TD.est}); (2) their standard errors (\code{TD.se}); (3) Fermat potential difference estimates (\code{FPD.est}); (4) their standard errors (\code{FPD.se}); (5) redshifts of deflectors (lenses) (\code{z.d}); and (6) redshifts of sources (quasars) (\code{z.s}). We note that the length of each input vector must be identical, and the order of input values in each vector must be consistent across all of the vectors. For example, let's say the third input value in one of the input vectors, \code{TD.est}, is a time delay estimate between gravitationally lensed images A and C of a specific quasar. Then the third input value in the other input vectors must be a quantity corresponding to the gravitationally lensed images A and C of that specific quasar.
}
\details{
Time delay cosmography infers the Hubble constant \eqn{H_0} by modeling time delays caused by strong gravitational lensing. The following equation is the basis of time delay cosmography:
  \deqn{c\Delta_{ijk}=D_{\Delta}(H_0, z_k, \Omega) \phi_{ijk},}
where \eqn{c} denotes the speed of light, and \eqn{\Delta_{ijk}} is the time delay in days between lensed images \eqn{i} and \eqn{j} of quasar \eqn{k} (\eqn{k=1, 2, \ldots, n}). The notation \eqn{z_k=\{z_{sk}, z_{dk}\}} indicates a vector for two redshifts of the \eqn{k}-th source (\eqn{z_{sk}}) and deflector (\eqn{z_{dk}}). The next notation \eqn{\Omega=\{\Omega_{m}, \Omega_{\Lambda}\}} indicates a vector for the present-day dark matter density \eqn{\Omega_{m}} and dark energy density \eqn{\Omega_{\Lambda}}. Since their sum becomes one under the flat \eqn{\Lambda}CDM model (standard cosmology),  we treat \eqn{\Omega_{m}} as the only unknown parameter. The notation  \eqn{D_{\Delta}(H_0, z_k, \Omega)} denotes the time delay distance in the unit of megaparsec (Mpc), which is a deterministic function of \eqn{H_0}, \eqn{z_k}, and \eqn{\Omega}.

In addition, the meta-analysis accounts for the impact of the mass structure in the line of sight between the observer and the lens via external convergences, \eqn{\kappa_{ext, k}}, where \eqn{k} is the lens index (\eqn{k=1, 2, ..., K}). The Hubble constant after accounting for this impact is obtained by \eqn{H_0^{ext}=(1-\kappa_{ext, k})H_0}, where \eqn{H_0} is the Hubble constant without considering the external convergense.

The meta analysis takes advantage of these physical relationships, and returns posterior samples of the model parameters, \eqn{H_0} \eqn{\Omega}, and \eqn{\kappa_{ext, k}}'s.

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

#####  The following is an example code to obtain 10 posterior samples of the model parameters.
initial.values <- c(runif(1, 0, 150), runif(1, 0.05, 0.5), rcauchy(2, scale = 0.025))
res <- h0(TD.est = time.delay.est, TD.se = time.delay.se, 
          FPD.est = fermat.diff.est, FPD.se = fermat.diff.se,
          z.d = z.d.obs, z.s = z.s.obs, 
          initial.param = initial.values,
          sample.size = 10, burnin.size = 5)
}
\author{
Hyungsuk Tak
}
\references{
    T. Treu and P. Marshall (2016) <doi:10.1007/s00159-016-0096-8> "Time Delay Cosmography".
}
