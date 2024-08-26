\name{Nakagami distribution}
\alias{se_naka}
\alias{re_naka}
\alias{hce_naka}
\alias{ae_naka}
\title{Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the Nakagami distribution}
\usage{
se_naka(alpha, beta)
re_naka(alpha, beta, delta)
hce_naka(alpha, beta, delta)
ae_naka(alpha, beta, delta)
}
\description{
Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the Nakagami distribution.
}

\arguments{
\item{alpha}{
The strictly positive scale parameter of the Nakagami distribution (\eqn{\alpha > 0}).
}
\item{beta}{
The strictly positive shape parameter of the Nakagami distribution (\eqn{\beta > 0}).
}
\item{delta}{
The strictly positive parameter  (\eqn{\delta > 0}) and (\eqn{\delta \ne 1}).
}
}

\details{The following is the probability density function of the Nakagami  distribution:
\deqn{
f(x)=\frac{2\alpha^{\alpha}}{\Gamma(\alpha)\beta^{\alpha}}x^{2\alpha-1}e^{-\frac{\alpha x^{2}}{\beta}},
}
where \eqn{x > 0}, \eqn{\alpha > 0} and \eqn{\beta > 0}, and \eqn{\Gamma(a)} is the standard gamma function.
}
\value{
The functions se_naka, re_naka, hce_naka, and ae_naka provide the Shannon entropy, Rényi entropy, Havrda and Charvat entropy, and Arimoto entropy, respectively, depending on the selected parametric values of the Nakagami distribution and \eqn{\delta}.
}
\references{
Schwartz, J., Godwin, R. T., & Giles, D. E. (2013). Improved maximum-likelihood estimation of the shape parameter in the Nakagami distribution. Journal of Statistical Computation and Simulation, 83(3), 434-445.


}

\author{
Muhammad Imran, Christophe Chesneau and Farrukh Jamal

R implementation and documentation: Muhammad Imran <imranshakoor84@yahoo.com>, Christophe Chesneau <christophe.chesneau@unicaen.fr> and Farrukh Jamal \email{farrukh.jamal@iub.edu.pk}.
}


%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{re_exp}, \link{re_gamma}, \link{re_wei}
}
}


\examples{
se_naka(1.2, 0.2)
delta <- c(1.5, 2, 3)
re_naka(1.2, 0.2, delta)
hce_naka(1.2, 0.2, delta)
ae_naka(1.2, 0.2, delta)
}


