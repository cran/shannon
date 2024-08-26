\name{Gamma distribution}
\alias{Se_gamma}
\alias{re_gamma}
\alias{hce_gamma}
\alias{ae_gamma}
\title{Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the gamma distribution}
\usage{
Se_gamma(alpha, beta)
re_gamma(alpha, beta, delta)
hce_gamma(alpha, beta, delta)
ae_gamma(alpha, beta, delta)
}
\description{
Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the gamma distribution.
}

\arguments{
\item{alpha}{
The strictly positive shape parameter of the gamma distribution (\eqn{\alpha > 0}).
}
\item{beta}{
The strictly positive scale parameter of the gamma distribution (\eqn{\beta > 0}).
}
\item{delta}{
The strictly positive parameter  (\eqn{\delta > 0}) and (\eqn{\delta \ne 1}).
}
}

\details{The following is the probability density function of the gamma distribution:
\deqn{
f(x)=\frac{\beta^{\alpha}}{\Gamma(\alpha)}x^{\alpha-1}e^{-\beta x},
}
where \eqn{x > 0}, \eqn{\alpha > 0} and \eqn{\beta > 0}, and \eqn{\Gamma(a)} is the standard gamma function.
}
\value{
The functions Se_gamma, re_gamma, hce_gamma, and ae_gamma provide the Shannon entropy, Rényi entropy, Havrda and Charvat entropy, and Arimoto entropy, respectively, depending on the selected parametric values of the gamma distribution and \eqn{\delta}.
}
\references{
Burgin, T. A. (1975). The gamma distribution and inventory control. Journal of the Operational Research Society, 26(3), 507-525.
}

\author{
Muhammad Imran, Christophe Chesneau and Farrukh Jamal

R implementation and documentation: Muhammad Imran <imranshakoor84@yahoo.com>, Christophe Chesneau <christophe.chesneau@unicaen.fr> and Farrukh Jamal \email{farrukh.jamal@iub.edu.pk}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{re_exp}, \link{re_wei}
}
}


\examples{
Se_gamma(1.2, 1.4)
delta <- c(1.5, 2, 3)
re_gamma(1.2, 1.4, delta)
hce_gamma(1.2, 1.4, delta)
ae_gamma(1.2, 1.4, delta)
}


