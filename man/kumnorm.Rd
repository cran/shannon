\name{Kumaraswamy normal  distribution}
\alias{se_kumnorm}
\alias{re_kumnorm}
\alias{hce_kumnorm}
\alias{ae_kumnorm}
\title{Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the Kumaraswamy normal  distribution}
\usage{
se_kumnorm(mu, sigma, a, b)
re_kumnorm(mu, sigma, a, b, delta)
hce_kumnorm(mu, sigma, a, b, delta)
ae_kumnorm(mu, sigma, a, b, delta)
}
\description{
Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the Kumaraswamy normal  distribution.
}

\arguments{
\item{mu}{
The location parameter of the normal distribution (\eqn{\mu\in\left(-\infty,+\infty\right)).
}
}
\item{sigma}{
The strictly positive scale parameter of the normal distribution (\eqn{\sigma > 0}).
}
\item{a}{
The strictly positive shape parameter of the Kumaraswamy distribution (\eqn{a > 0}).
}
\item{b}{
The strictly positive shape parameter of the Kumaraswamy distribution (\eqn{b > 0}).
}
\item{delta}{
The strictly positive parameter  (\eqn{\delta > 0}) and (\eqn{\delta \ne 1}).
}
}

\details{The following is the probability density function of the Kumaraswamy normal distribution:
\deqn{
f(x)=\frac{ab}{\sigma}\phi\left(\frac{x-\mu}{\sigma}\right)\left[\Phi\left(\frac{x-\mu}{\sigma}\right)\right]^{a-1}\left[1-\Phi\left(\frac{x-\mu}{\sigma}\right)^{a}\right]^{b-1},
}
where \eqn{x\in\left(-\infty,+\infty\right)}, \eqn{\mu\in\left(-\infty,+\infty\right)}, \eqn{\sigma > 0}, \eqn{a > 0} and \eqn{b > 0}, and the functions \eqn{\phi(t)} and \eqn{\Phi(t) }, denote the probability density function and cumulative distribution function of the standard normal distribution, respectively.
}
\value{
The functions se_kumnorm, re_kumnorm, hce_kumnorm, and ae_kumnorm provide the Shannon entropy, Rényi entropy, Havrda and Charvat entropy, and Arimoto entropy, respectively, depending on the selected parametric values of the Kumaraswamy normal  distribution and \eqn{\delta}.
}
\references{
Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized distributions. Journal of Statistical Computation and Simulation, 81(7), 883-898.

}

\author{
Muhammad Imran, Christophe Chesneau and Farrukh Jamal

R implementation and documentation: Muhammad Imran <imranshakoor84@yahoo.com>, Christophe Chesneau <christophe.chesneau@unicaen.fr> and Farrukh Jamal \email{farrukh.jamal@iub.edu.pk}.
}


%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{re_norm}, \link{re_kum}
}
}


\examples{
se_kumnorm(0.2, 1.5, 1, 1)
delta <- c(1.5, 2, 3)
re_kumnorm(1.2, 1, 2, 1.5, delta)
hce_kumnorm(1.2, 1, 2, 1.5, delta)
ae_kumnorm(1.2, 1, 2, 1.5, delta)
}


