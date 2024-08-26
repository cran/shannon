\name{Kumaraswamy exponential distribution}
\alias{re_kexp}
\alias{hce_kexp}
\alias{ae_kexp}
\title{Compute the Rényi, Havrda and Charvat, and Arimoto entropies of the Kumaraswamy exponential distribution}
\usage{
re_kexp(lambda, a, b, delta)
hce_kexp(lambda, a, b, delta)
ae_kexp(lambda, a, b, delta)
}
\description{
Compute the Rényi, Havrda and Charvat, and Arimoto entropies of the Kumaraswamy exponential distribution.
}

\arguments{
\item{a}{
The strictly positive shape parameter of the Kumaraswamy distribution (\eqn{a > 0}).
}
\item{b}{
The strictly positive shape parameter of the Kumaraswamy distribution (\eqn{b > 0}).
}
\item{lambda}{
The strictly positive parameter of the exponential distribution (\eqn{\lambda > 0}).
}
\item{delta}{
The strictly positive parameter  (\eqn{\delta > 0}) and (\eqn{\delta \ne 1}).
}
}

\details{The following is the probability density function of the Kumaraswamy exponential distribution:
\deqn{
f(x)=ab\lambda e^{-\lambda x}\left(1-e^{-\lambda x}\right)^{a-1}\left\{ 1-\left(1-e^{-\lambda x}\right)^{a}\right\} ^{b-1},
}
where \eqn{x > 0}, \eqn{a > 0}, \eqn{b > 0} and \eqn{\lambda > 0}.
}
\value{
The functions re_kexp, hce_kexp, and ae_kexp provide the Rényi entropy, Havrda and Charvat entropy, and Arimoto entropy, respectively, depending on the selected parametric values of the Kumaraswamy exponential distribution and \eqn{\delta}.
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
\code{\link{re_exp}, \link{re_kum}
}
}


\examples{
delta <- c(1.5, 2, 3)
re_kexp(1.2, 1.2, 1.4, delta)
hce_kexp(1.2, 1.2, 1.4, delta)
ae_kexp(1.2, 1.2, 1.4, delta)
}


