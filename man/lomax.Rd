\name{Lomax distribution}
\alias{se_lom}
\alias{re_lom}
\alias{hce_lom}
\alias{ae_lom}
\title{Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the Lomax distribution}
\usage{
se_lom(alpha, beta)
re_lom(alpha, beta, delta)
hce_lom(alpha, beta, delta)
ae_lom(alpha, beta, delta)
}
\description{
Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the Lomax distribution.
}

\arguments{
\item{alpha}{
The strictly positive shape parameter of the Lomax distribution (\eqn{\alpha > 0}).
}
\item{beta}{
The strictly positive scale parameter of the Lomax distribution (\eqn{\beta > 0}).
}
\item{delta}{
The strictly positive parameter  (\eqn{\delta > 0}) and (\eqn{\delta \ne 1}).
}
}

\details{The following is the probability density function of the Lomax distribution:
\deqn{
f(x)=\frac{\alpha}{\beta}\left(1+\frac{x}{\beta}\right)^{-\alpha-1},
}
where \eqn{x > 0}, \eqn{\alpha > 0} and \eqn{\beta > 0}.
}
\value{
The functions se_lom, re_lom, hce_lom, and ae_lom provide the Shannon entropy, Rényi entropy, Havrda and Charvat entropy, and Arimoto entropy, respectively, depending on the selected parametric values of the Lomax distribution and \eqn{\delta}.
}
\references{
Abd-Elfattah, A. M., Alaboud, F. M., & Alharby, A. H. (2007). On sample size estimation for Lomax distribution. Australian Journal of Basic and Applied Sciences, 1(4), 373-378.

}

\author{
Muhammad Imran, Christophe Chesneau and Farrukh Jamal

R implementation and documentation: Muhammad Imran <imranshakoor84@yahoo.com>, Christophe Chesneau <christophe.chesneau@unicaen.fr> and Farrukh Jamal \email{farrukh.jamal@iub.edu.pk}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{re_exp}, \link{re_gamma}
}
}


\examples{
se_lom(1.2, 0.2)
delta <- c(1.5, 2, 3)
re_lom(1.2, 0.2, delta)
hce_lom(1.2, 0.2, delta)
ae_lom(1.2, 0.2, delta)
}


