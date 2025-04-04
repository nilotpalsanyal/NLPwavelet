\name{NLPwavelet-package}

\alias{NLPwavelet-package}

\alias{NLPwavelet}

\docType{package}

\title{
  Bayesian Wavelet Analysis Using Non-local Priors
}

\description{
  Performs Bayesian wavelet analysis using individual non-local priors as described in Sanyal & Ferreira (2017) <DOI:10.1007/s13571-016-0129-3> and non-local prior mixtures as described in Sanyal (2025) <DOI:10.48550/arXiv.2501.18134>.
}

\details{
  The main function is \code{\link{BNLPWA}}, which has arguments for specifying analysis using individual non-local priors or non-local prior mixtures and various hyperparameter specifications for the wavelet coefficients and scale parameters of the non-local priors. See the manual of \code{\link{BNLPWA}} for examples.
}

\author{
Nilotpal Sanyal <nsanyal@utep.edu>

Maintainer: Nilotpal Sanyal <nsanyal@utep.edu>
}

\references{
  Sanyal, Nilotpal. "Nonlocal prior mixture-based Bayesian wavelet regression." arXiv preprint arXiv:2501.18134 (2025).

  Sanyal, Nilotpal, and Marco AR Ferreira. "Bayesian wavelet analysis using nonlocal priors with an application to FMRI analysis." Sankhya B 79.2 (2017): 361-388.
}
% \keyword{ package }
% \seealso{
%   Optional links to other man pages
% }
% \examples{
%   \dontrun{
%      ## Optional simple examples of the most important functions
%      ## These can be in \dontrun{} and \donttest{} blocks.   
%   }
% }