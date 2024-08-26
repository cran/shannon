\name{Logistic distribution}
\alias{se_logis}
\alias{re_logis}
\alias{hce_logis}
\alias{ae_logis}
\title{Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the logistic distribution}
\usage{
se_logis(mu, sigma)
re_logis(mu, sigma, delta)
hce_logis(mu, sigma, delta)
ae_logis(mu, sigma, delta)
}
\description{
Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the logistic distribution.
}

\arguments{
\item{mu}{
The location parameter of the logistic distribution (\eqn{\mu\in\left(-\infty,+\infty\right)}).
}
\item{sigma}{
The strictly positive scale parameter of the logistic distribution (\eqn{\sigma > 0}).
}
\item{delta}{
The strictly positive parameter  (\eqn{\delta > 0}) and (\eqn{\delta \ne 1}).
}
}
\details{ The following is the probability density function of the logistic distribution:
\deqn{
f(x)=\frac{e^{-\frac{\left(x-\mu\right)}{\sigma}}}{\sigma\left(1+e^{-\frac{\left(x-\mu\right)}{\sigma}}\right)^{2}},
}
where \eqn{x\in\left(-\infty,+\infty\right)}, \eqn{\mu\in\left(-\infty,+\infty\right)} and \eqn{\sigma > 0}.
}
\value{
The functions se_logis, re_logis, hce_logis, and ae_logis provide the Shannon entropy, Rényi entropy, Havrda and Charvat entropy, and Arimoto entropy, respectively, depending on the selected parametric values of the logistic distribution and \eqn{\delta}.
}
\references{
Johnson, N. L., Kotz, S., & Balakrishnan, N. (1995). Continuous univariate distributions, Volume 2 (Vol. 289). John Wiley & Sons.

}
\author{
Muhammad Imran, Christophe Chesneau and Farrukh Jamal

R implementation and documentation: Muhammad Imran <imranshakoor84@yahoo.com>, Christophe Chesneau <christophe.chesneau@unicaen.fr> and Farrukh Jamal \email{farrukh.jamal@iub.edu.pk}.
}


%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{re_gum}, \link{re_norm}
}
}


\examples{
se_logis(0.2, 1.4)
delta <- c(2, 3)
re_logis(1.2, 0.4, delta)
hce_logis(1.2, 0.4, delta)
ae_logis(1.2, 0.4, delta)
}


