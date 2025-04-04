\name{BNLPWA}
\alias{BNLPWA}
\title{Bayesian Non-Local Prior-Based Wavelet Analysis}
\description{\code{BNLPWA} is the main function of this package that performs Bayesian wavelet analysis using individual non-local priors as described in Sanyal & Ferreira (2017) and non-local prior mixtures as described in Sanyal (2025). It currently works with one-dimensional data. The usage is described below.}
\usage{
BNLPWA(
  data, 
  func=NULL, 
  method=c("mixture","mom","imom"), 
  mixprob_dist=c("logit","genlogit","hypsec","gennormal"), 
  scale_dist=c("polynom","doubleexp"),
  r=1, 
  nu=1, 
  wave.family="DaubLeAsymm", 
  filter.number=6, 
  bc="periodic"
)
}
\arguments{
  \item{data}{Vector of noisy data.}
  \item{func}{Vector of true functional values. NULL by default. If available, they are used to compute and return mean squared error (MSE) of the estimates.}
  \item{method}{"mixture" for non-local prior mixture-based analysis, "mom" or "imom" for individual non-local prior-based analysis.}
  \item{mixprob_dist}{Specification for the mixture probabilities of the spike-and-slab prior. Equations given in the Details.}
  \item{scale_dist}{Specification for the scale parameters of the non-local priors. Equations given in the Details.}
  \item{r}{Integer specifying a) the order of the MOM prior or the shape parameter of the IMOM prior for individual non-local prior-based analysis, or b) the order of the MOM prior for non-local prior mixture-based analysis.}
  \item{nu}{Integer specifying the shape parameter of the IMOM prior for non-local prior mixture-based analysis. Not used for individual non-local prior-based analysis.}
  \item{wave.family}{The family of wavelets to use - "DaubExPhase" or "DaubLeAsymm". Default is "DaubLeAsymm".}
  \item{filter.number}{The number of vanishing moments of the wavelet. Default is 6.}
  \item{bc}{The boundary condition to use - "periodic" or "symmetric". Default is "periodic".}
}
\value{A list containing the following.
  \item{data }{The data vector.}
  \item{func.post.mean }{Posterior estimate (mean) of the function that generated the data.}
  \item{wavelet.empirical }{Empirical wavelet coefficients obtained by wavelet transformation of the data.}
  \item{wavelet.post.mean }{Posterior estimate (mean) of the true wavelet coefficients obtained by wavelet transformation of the underlying function. }
  \item{hyperparam }{Estimates of the hyperparameters that specify the spike-and-slab prior for the wavelet coefficients.}
  \item{sigma }{Estimate of the standard deviation of the error.}
  \item{MSE.mean }{Mean squared error of the estimates, computable only if true functional values are supplied in the input argument \code{func}. }
  \item{runtime }{System run-time of the function.}
}
\details{
  
  \subsection{Spike-and-slab prior for the wavelet coefficients}{

    For individual MOM prior-based analysis, the spike-and-slab prior for the wavelet coefficient \eqn{d_{lj}} is given by
    \deqn{ d_{lj} \mid \gamma_l, \tau_l, \sigma^2, r \sim \gamma_l \; \text{MOM}(\tau_l, \sigma^2, r) + (1-\gamma_l) \; \delta(0), }

    for individual IMOM prior-based analysis, the spike-and-slab prior for the wavelet coefficient \eqn{d_{lj}} is given by
    \deqn{ d_{lj} \mid \gamma_l, \tau_l, \sigma^2, r \sim \gamma_l \; \text{IMOM}(\tau_l, \sigma^2, r) + (1-\gamma_l) \; \delta(0), }

    and for non-local prior mixture-based analysis, the spike-and-slab prior for the wavelet coefficient \eqn{d_{lj}} is given by
    \deqn{ d_{lj} \mid \gamma_l^{(1)}, \gamma_l^{(2)}, \tau_l^{(1)}, \tau_l^{(2)}, \sigma^2, r, \nu \sim \gamma_l^{(1)} \; \text{MOM}(\tau_l^{(1)}, r, \sigma^2) + (1-\gamma_l^{(1)})\gamma_l^{(2)} \;\text{IMOM}(\tau_l^{(2)}, \nu, \sigma^2) + (1-\gamma_l^{(1)})(1-\gamma_l^{(2)}) \; \delta(0), }

    where the density of the MOM prior is
    \deqn{ mom(d_{lj} | \tau_{l}^{(1)},r,\sigma^{2}) = \widetilde{M}_{r} \left(\tau_{l}^{(1)}\sigma^{2}\right)^{-r-1/2} d_{lj}^{2r} \exp\left(-\frac{d_{lj}^{2}}{2\tau_{l}^{(1)}\sigma^{2}}\right), \quad r>1, \tau_{l}^{(1)}>0, \widetilde{M}_{r} = \frac{(2\pi)^{-1/2}}{(2r-1)!!} }

    and the density of the IMOM prior is
    \deqn{ imom(d_{lj} | \tau_{l}^{(2)},\nu,\sigma^{2}) = \frac{\left(\tau_{l}^{(2)}\sigma^{2}\right)^{\nu/2}}{\Gamma(\nu/2)} |d_{lj}|^{-\nu-1} \exp\left( -\frac{\tau_{l}^{(2)}\sigma^{2}}{d_{lj}^{2}} \right),\quad \nu>1, \tau_{l}^{(2)}>0. } 
  }
  

  \subsection{Hyperparameter specifications}{

    For non-local prior mixture-based analysis, the available specifications for the mixture probabilities are
    \enumerate{
      \item \bold{Logit:} 
      \deqn{
        \gamma_l^{(1)} = \frac{\exp(\theta^{\gamma}_{1} - \theta^{\gamma}_{2}l)}
        {1 + \exp(\theta^{\gamma}_{1} - \theta^{\gamma}_{2}l)}, \quad
        \theta^{\gamma}_{1} \in \mathbb{R}, \; \theta^{\gamma}_{2} > 0
      }{
        gamma_l^(1) = exp(theta^gamma_1 - theta^gamma_2 * l) / 
        (1 + exp(theta^gamma_1 - theta^gamma_2 * l)), 
        where theta^gamma_1 in R, theta^gamma_2 > 0
      }
      \deqn{
        \gamma_l^{(2)} = \frac{\exp(\theta^{\gamma}_{3} - \theta^{\gamma}_{4}l)}
        {1 + \exp(\theta^{\gamma}_{3} - \theta^{\gamma}_{4}l)}, \quad
        \theta^{\gamma}_{3} \in \mathbb{R}, \; \theta^{\gamma}_{4} > 0
      }{
        gamma_l^(2) = exp(theta^gamma_3 - theta^gamma_4 * l) / 
        (1 + exp(theta^gamma_3 - theta^gamma_4 * l)), 
        where theta^gamma_3 in R, theta^gamma_4 > 0
      }

      \item \bold{Generalized logit or Richards:} 
      \deqn{
        \gamma_l^{(1)} = \frac{1}{[1 + \exp\{-(\theta^{\gamma}_{1} - \theta^{\gamma}_{2}l)\}]^{\theta^{\gamma}_{3}}}, \quad 
        \theta^{\gamma}_{1} \in \mathbb{R}, \; \theta^{\gamma}_{2},\theta^{\gamma}_{3} > 0 
      }{
        gamma_l^(1) = 1 / [1 + exp(-(theta^gamma_1 - theta^gamma_2 * l))]^(theta^gamma_3), 
        where theta^gamma_1 is in R, and theta^gamma_2, theta^gamma_3 > 0
      }
      \deqn{
        \gamma_l^{(2)} = \frac{1}{[1 + \exp\{-(\theta^{\gamma}_{4} - \theta^{\gamma}_{5}l)\}]^{\theta^{\gamma}_{6}}}, \quad
        \theta^{\gamma}_{4} \in \mathbb{R}, \; \theta^{\gamma}_{5},\theta^{\gamma}_{6} > 0;
      }{
        gamma_l^(2) = 1 / [1 + exp(-(theta^gamma_4 - theta^gamma_5 * l))]^(theta^gamma_6), 
        where theta^gamma_4 is in R, and theta^gamma_5, theta^gamma_6 > 0
      }

      \item \bold{Hyperbolic secant:}
      \deqn{
        \gamma_l^{(1)} = \frac{2}{\pi} \arctan\left[\exp\left(\frac{\pi}{2} \left(\theta^{\gamma}_{1} - \theta^{\gamma}_{2}l\right)\right)\right], \quad
        \theta^{\gamma}_{1} \in \mathbb{R}, \; \theta^{\gamma}_{2} > 0
      }{
        gamma_l^(1) = (2 / pi) * arctan(exp((pi / 2) * (theta^gamma_1 - theta^gamma_2 * l))),
        where theta^gamma_1 is in R, and theta^gamma_2 > 0
      }
      \deqn{
        \gamma_l^{(2)} = \frac{2}{\pi} \arctan\left[\exp\left(\frac{\pi}{2} \left(\theta^{\gamma}_{3} - \theta^{\gamma}_{4}l\right)\right)\right], \quad
        \theta^{\gamma}_{3} \in \mathbb{R}, \; \theta^{\gamma}_{4} > 0
      }{
        gamma_l^(2) = (2 / pi) * arctan(exp((pi / 2) * (theta^gamma_3 - theta^gamma_4 * l))),
        where theta^gamma_3 is in R, and theta^gamma_4 > 0
      }

      \item \bold{Generalized normal:}
      \deqn{
        \gamma_l^{(1)} = \frac{1}{2} + \text{sign}(\theta^{\gamma}_{1}-l) \frac{1}{2\Gamma(1/\theta^{\gamma}_{2})} 
        \gamma\left(1/\theta^{\gamma}_{2} ,\left|\frac{\theta^{\gamma}_{1}-l}{\theta^{\gamma}_{3}}\right|^{\theta^{\gamma}_{2}}\right), \quad 
        \theta^{\gamma}_{1} \in \mathbb{R}, \; \theta^{\gamma}_{2},\theta^{\gamma}_{3} > 0
      }{
        gamma_l^(1) = (1/2) + sign(theta^gamma_1 - l) * (1 / (2 * Gamma(1 / theta^gamma_2))) * 
        gamma(1 / theta^gamma_2, |(theta^gamma_1 - l) / theta^gamma_3| ^ theta^gamma_2),
        where theta^gamma_1 is in R, and theta^gamma_2, theta^gamma_3 > 0
      }
      \deqn{
        \gamma_l^{(2)} = \frac{1}{2} + \text{sign}(\theta^{\gamma}_{4}-l) \frac{1}{2\Gamma(1/\theta^{\gamma}_{5})} 
        \gamma\left(1/\theta^{\gamma}_{5} ,\left|\frac{\theta^{\gamma}_{4}-l}{\theta^{\gamma}_{6}}\right|^{\theta^{\gamma}_{5}}\right), \quad
        \theta^{\gamma}_{4} \in \mathbb{R}, \; \theta^{\gamma}_{5},\theta^{\gamma}_{6} > 0
      }{
        gamma_l^(2) = (1/2) + sign(theta^gamma_4 - l) * (1 / (2 * Gamma(1 / theta^gamma_5))) * 
        gamma(1 / theta^gamma_5, |(theta^gamma_4 - l) / theta^gamma_6| ^ theta^gamma_5),
        where theta^gamma_4 is in R, and theta^gamma_5, theta^gamma_6 > 0
      }
    }
    For individual non-local prior based analysis, \eqn{gamma_l} is defined likewise. 

    For non-local prior mixture-based analysis, the available specifications for the scale parameters are
    \enumerate{
      \item \bold{Polynomial decay:}
      \deqn{
        \tau_{l}^{(1)} = \theta^{\tau}_{1} l^{-\theta^{\tau}_{2}}, \quad \theta^{\tau}_{1},\theta^{\tau}_{2} > 0
      }{
        tau_l^(1) = theta^tau_1 * l^(-theta^tau_2), where theta^tau_1, theta^tau_2 > 0
      }
      \deqn{
        \tau_{l}^{(2)} = \theta^{\tau}_{3} l^{-\theta^{\tau}_{4}}, \quad \theta^{\tau}_{3},\theta^{\tau}_{4} > 0
      }{
        tau_l^(2) = theta^tau_3 * l^(-theta^tau_4), where theta^tau_3, theta^tau_4 > 0
      }

      \item \bold{Double-exponential decay:}
      \deqn{
        \tau_{l}^{(1)} = \theta^{\tau}_{1} \exp(-\theta^{\tau}_{2} l) + \theta^{\tau}_{3} \exp(-\theta^{\tau}_{4} l), \quad
        \theta^{\tau}_{1},\theta^{\tau}_{2},\theta^{\tau}_{3},\theta^{\tau}_{4} > 0
      }{
        tau_l^(1) = theta^tau_1 * exp(-theta^tau_2 * l) + theta^tau_3 * exp(-theta^tau_4 * l), 
        where theta^tau_1, theta^tau_2, theta^tau_3, theta^tau_4 > 0
      }
      \deqn{
        \tau_{l}^{(2)} = \theta^{\tau}_{5} \exp(-\theta^{\tau}_{6} l) + \theta^{\tau}_{7} \exp(-\theta^{\tau}_{8} l), \quad
        \theta^{\tau}_{5},\theta^{\tau}_{6},\theta^{\tau}_{7},\theta^{\tau}_{8} > 0
      }{
        tau_l^(2) = theta^tau_5 * exp(-theta^tau_6 * l) + theta^tau_7 * exp(-theta^tau_8 * l), 
        where theta^tau_5, theta^tau_6, theta^tau_7, theta^tau_8 > 0
      }
    }
    For individual non-local prior based analysis, \eqn{tau_l} is defined likewise.

  }

  Note: The wavelet computations are performed by using the R package \pkg{wavethresh}.
}
\references{
  Sanyal, Nilotpal. "Nonlocal prior mixture-based Bayesian wavelet regression." arXiv preprint arXiv:2501.18134 (2025).

  Sanyal, Nilotpal, and Marco AR Ferreira. "Bayesian wavelet analysis using nonlocal priors with an application to FMRI analysis." Sankhya B 79.2 (2017): 361-388.
}
\seealso{
  \code{\link[wavethresh]{wd}}, \code{\link[wavethresh]{wr}}
}
\author{Nilotpal Sanyal <nsanyal@utep.edu>

Maintainer: Nilotpal Sanyal <nsanyal@utep.edu>
}
\examples{
  # Using the well-known Doppler function to 
  # illustrate the use of the function BNLPWA

  # set seed for reproducibility
  set.seed(1)

  # Define the doppler function
  doppler <- function(x) { 
    sqrt(x * (1 - x)) * sin((2 * pi * 1.05) / (x + 0.05)) 
  }

  # Generate true values over a grid of length an integer power of 2 
  n <- 128  # Number of points
  x <- seq(0, 1, length.out = n)
  true_signal <- doppler(x) 

  # Add noise to generate data
  sigma <- 0.2  # Noise level
  y <- true_signal + rnorm(n, mean = 0, sd = sigma)

  # BNLPWA analysis based on MOM prior using logit specification
  # for the mixture probabilities and polynomial decay
  # specification for the scale parameter
  fit_mom <- BNLPWA(data=y, func=true_signal, r=1, wave.family=
    "DaubLeAsymm", filter.number=6, bc="periodic", method="mom", 
    mixprob_dist="logit", scale_dist="polynom")

  plot(y,type="l",col="grey") # plot of data
  lines(fit_mom$func.post.mean,col="blue") # plot of posterior 
  # smoothed estimates
  fit_mom$MSE.mean

  \donttest{
    # BNLPWA analysis using non-local prior mixtures using generalized 
    # logit (Richard's) specification for the mixture probabilities and 
    # double exponential decay specification for the scale parameter
    fit_mixture <- BNLPWA(data=y, func=true_signal, r=1, nu=1, wave.family=
      "DaubLeAsymm", filter.number=6, bc="periodic", method="mixture", 
      mixprob_dist="genlogit", scale_dist="doubleexp")

    plot(y,type="l",col="grey") # plot of data
    lines(fit_mixture$func.post.mean,col="blue") # plot of posterior 
    # smoothed estimates
    fit_mixture$MSE.mean
  }
  
  # Compare with other wavelet methods
  library(wavethresh)
  wd <- wd(y, family="DaubLeAsymm", filter.number=6, bc="periodic")  # Wavelet decomposition  
  
  wd_thresh_universal <- threshold(wd, policy="universal", type="hard")
  fit_universal <- wr(wd_thresh_universal)
  MSE_universal <- mean((true_signal-fit_universal)^2)
  MSE_universal

  wd_thresh_sure <- threshold(wd, policy="sure", type="soft")
  fit_sure <- wr(wd_thresh_sure)
  MSE_sure <- mean((true_signal-fit_sure)^2)
  MSE_sure

  wd_thresh_BayesThresh <- threshold(wd, policy="BayesThresh", type="hard")
  fit_BayesThresh <- wr(wd_thresh_BayesThresh)
  MSE_BayesThresh <- mean((true_signal-fit_BayesThresh)^2)
  MSE_BayesThresh

  wd_thresh_cv <- threshold(wd, policy="cv", type="hard")
  fit_cv <- wr(wd_thresh_cv)
  MSE_cv <- mean((true_signal-fit_cv)^2)
  MSE_cv

  wd_thresh_fdr <- threshold(wd, policy="fdr", type="hard")
  fit_fdr <- wr(wd_thresh_fdr)
  MSE_fdr <- mean((true_signal-fit_fdr)^2)
  MSE_fdr

  # Compare with non-wavelet methods
      # Kernel smoothing
  fit_ksmooth <- ksmooth(x, y, kernel="normal", bandwidth=0.05)
  MSE_ksmooth <- mean((true_signal-fit_ksmooth$y)^2)
  MSE_ksmooth
      # LOESS smoothing
  fit_loess <- loess(y ~ x, span=0.1)  # Adjust span for more or less smoothing
  MSE_loess <- mean((true_signal-predict(fit_loess))^2)
  MSE_loess
      # Cubic spline smoothing
  spline_fit <- smooth.spline(x, y, spar=0.5)  # Adjust spar for smoothness
  MSE_spline <- mean((true_signal-spline_fit$y)^2)
  MSE_spline
}










