\name{Log-normal distribution}
\alias{se_lnorm}
\alias{re_lnorm}
\alias{hce_lnorm}
\alias{ae_lnorm}
\title{Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the log-normal distribution}
\usage{
se_lnorm(mu, sigma)
re_lnorm(mu, sigma, delta)
hce_lnorm(mu, sigma, delta)
ae_lnorm(mu, sigma, delta)
}
\description{
Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the log-normal distribution.
}

\arguments{
\item{mu}{
The location parameter (\eqn{\mu\in\left(-\infty,+\infty\right)}).
}
\item{sigma}{
The strictly positive scale parameter of the log-normal distribution (\eqn{\sigma > 0}).
}
\item{delta}{
The strictly positive parameter  (\eqn{\delta > 0}) and (\eqn{\delta \ne 1}).
}
}
\details{The following is the probability density function of the log-normal distribution:
\deqn{
f(x)=\frac{1}{x\sigma\sqrt{2\pi}}e^{-\frac{\left(\log(x)-\mu\right)^{2}}{2\sigma^{2}}},
}
where \eqn{x > 0}, \eqn{\mu\in\left(-\infty,+\infty\right)} and \eqn{\sigma > 0}.
}
\value{
The functions se_lnorm, re_lnorm, hce_lnorm, and ae_lnorm provide the Shannon entropy, Rényi entropy, Havrda and Charvat entropy, and Arimoto entropy, respectively, depending on the selected parametric values of the log-normal distribution and \eqn{\delta}.
}
\references{
Johnson, N. L., Kotz, S., & Balakrishnan, N. (1995). Continuous univariate distributions, Volume 1, Chapter 14. Wiley, New York.

}
\author{
Muhammad Imran, Christophe Chesneau and Farrukh Jamal

R implementation and documentation: Muhammad Imran <imranshakoor84@yahoo.com>, Christophe Chesneau <christophe.chesneau@unicaen.fr> and Farrukh Jamal \email{farrukh.jamal@iub.edu.pk}.
}


%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{re_wei}, \link{re_norm}
}
}


\examples{
se_lnorm(0.2, 1.4)
delta <- c(2, 3)
re_lnorm(1.2, 0.4, delta)
hce_lnorm(1.2, 0.4, delta)
ae_lnorm(1.2, 0.4, delta)
}


