\name{Exponentiated Weibull distribution}
\alias{se_ew}
\alias{re_ew}
\alias{hce_ew}
\alias{ae_ew}
\title{Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the exponentiated Weibull distribution}
\usage{
se_ew(a, beta, zeta)
re_ew(a, beta, zeta, delta)
hce_ew(a, beta, zeta, delta)
ae_ew(a, beta, zeta, delta)
}
\description{
Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the exponentiated Weibull distribution.
}

\arguments{
\item{a}{
The strictly positive shape parameter of the exponentiated Weibull distribution (\eqn{a > 0}).
}
\item{beta}{
The strictly positive scale parameter of the baseline Weibull distribution (\eqn{\beta > 0}).
}
\item{zeta}{
The strictly positive shape parameter of the baseline Weibull distribution (\eqn{\zeta > 0}).
}
\item{delta}{
The strictly positive parameter  (\eqn{\delta > 0}) and (\eqn{\delta \ne 1}).
}
}

\details{The following is the probability density function of the exponentiated Weibull distribution:
\deqn{
f(x)=a\zeta\beta^{-\zeta}x^{\zeta-1}e^{-\left(\frac{x}{\beta}\right)^{\zeta}}\left[1-e^{-\left(\frac{x}{\beta}\right)^{\zeta}}\right]^{a-1},
}
where \eqn{x > 0}, \eqn{a > 0}, \eqn{\beta > 0} and \eqn{\zeta > 0}.
}
\value{
The functions se_ew, re_ew, hce_ew, and ae_ew provide the Shannon entropy, Rényi entropy, Havrda and Charvat entropy, and Arimoto entropy, respectively, depending on the selected parametric values of the exponentiated Weibull distribution and \eqn{\delta}.
}
\references{
Nadarajah, S., Cordeiro, G. M., & Ortega, E. M. (2013). The exponentiated Weibull distribution: a survey. Statistical Papers, 54, 839-877.
}

\author{
Muhammad Imran, Christophe Chesneau and Farrukh Jamal

R implementation and documentation: Muhammad Imran <imranshakoor84@yahoo.com>, Christophe Chesneau <christophe.chesneau@unicaen.fr> and Farrukh Jamal \email{farrukh.jamal@iub.edu.pk}.
}


%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{re_exp}, \link{re_wei}, \link{re_ew}
}
}


\examples{
se_ew(0.8, 0.2, 0.8)
delta <- c(1.5, 2, 3)
re_ew(1.2, 1.2, 1.4, delta)
hce_ew(1.2, 1.2, 1.4, delta)
ae_ew(1.2, 1.2, 1.4, delta)
}


