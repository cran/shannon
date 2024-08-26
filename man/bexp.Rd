\name{Beta exponential distribution}
\alias{se_bexp}
\alias{re_bexp}
\alias{hce_bexp}
\alias{ae_bexp}
\title{Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the beta exponential distribution}
\usage{
se_bexp(lambda, alpha, beta)
re_bexp(lambda, alpha, beta, delta)
hce_bexp(lambda, alpha, beta, delta)
ae_bexp(lambda, alpha, beta, delta)
}
\description{
Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the beta exponential distribution.
}
\arguments{
\item{lambda}{
The strictly positive scale parameter of the exponential  distribution (\eqn{\lambda > 0}).
}
\item{alpha}{
The strictly positive shape parameter of the  beta distribution (\eqn{\alpha > 0}).
}
\item{beta}{
The strictly positive shape parameter of the  beta  distribution (\eqn{\beta > 0}).
}
\item{delta}{
The strictly positive parameter  (\eqn{\delta > 0}) and (\eqn{\delta \ne 1}).
}
}

\details{The following is the probability density function of the  beta exponential distribution:
\deqn{
f(x)=\frac{\lambda e^{-\beta\lambda x}}{B(\alpha,\beta)}\left(1-e^{-\lambda x}\right)^{\alpha-1},
}
where \eqn{x > 0}, \eqn{\alpha > 0}, \eqn{\beta > 0} and \eqn{\lambda > 0}, and \eqn{B(a,b)} denotes the standard beta function.
}
\value{
The functions se_bexp, re_bexp, hce_bexp, and ae_bexp provide the Shannon entropy, Rényi entropy, Havrda and Charvat entropy, and Arimoto entropy, respectively, depending on the selected parametric values of the beta exponential distribution and \eqn{\delta}.
}
\references{
Nadarajah, S., & Kotz, S. (2006). The beta exponential distribution. Reliability Engineering & System Safety, 91(6), 689-697.

}

\author{
Muhammad Imran, Christophe Chesneau and Farrukh Jamal

R implementation and documentation: Muhammad Imran <imranshakoor84@yahoo.com>, Christophe Chesneau <christophe.chesneau@unicaen.fr> and Farrukh Jamal \email{farrukh.jamal@iub.edu.pk}.
}


%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{re_beta}, \link{re_exp}
}
}


\examples{
# Computation of the Shannon entropy
se_bexp(1.2, 0.2, 1.5)
delta <- c(0.2, 0.3, 0.5)
# Computation of the Rényi entropy
re_bexp(1.2, 0.2, 0.5, delta)
# Computation of the Havrda and Charvat entropy
hce_bexp(1.2, 0.2, 1.5, delta)
# Computation of the Arimoto entropy
ae_bexp(1.2, 0.2, 1.5, delta)
}


