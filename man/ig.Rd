\name{Inverse-gamma distribution}
\alias{se_ig}
\alias{re_ig}
\alias{hce_ig}
\alias{ae_ig}
\title{Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the inverse-gamma distribution}
\usage{
se_ig(alpha, beta)
re_ig(alpha, beta, delta)
hce_ig(alpha, beta, delta)
ae_ig(alpha, beta, delta)
}
\description{
Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the inverse-gamma distribution.
}

\arguments{
\item{alpha}{
The strictly positive shape parameter of the inverse-gamma  distribution (\eqn{\alpha > 0}).
}
\item{beta}{
The strictly positive scale parameter of the inverse-gamma  distribution (\eqn{\beta > 0}).
}
\item{delta}{
The strictly positive parameter  (\eqn{\delta > 0}) and (\eqn{\delta \ne 1}).
}
}

\details{The following is the probability density function of the inverse-gamma distribution:
\deqn{
f(x)=\frac{\beta^{\alpha}}{\Gamma(\alpha)}x^{-\alpha-1}e^{-\frac{\beta}{x}},
}
where \eqn{x > 0}, \eqn{\alpha > 0} and \eqn{\beta > 0}, and \eqn{\Gamma(a)} is the standard gamma function.
}
\value{
The functions se_ig, re_ig, hce_ig, and ae_ig provide the Shannon entropy, Rényi entropy, Havrda and Charvat entropy, and Arimoto entropy, respectively, depending on the selected parametric values of the inverse-gamma distribution and \eqn{\delta}.
}
\references{
Rivera, P. A., Calderín-Ojeda, E., Gallardo, D. I., & Gómez, H. W. (2021). A compound class of the inverse Gamma and power series distributions. Symmetry, 13(8), 1328.

Glen, A. G. (2017). On the inverse gamma as a survival distribution. Computational Probability Applications, 15-30.

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
se_ig(1.2, 0.2)
delta <- c(1.5, 2, 3)
re_ig(1.2, 0.2, delta)
hce_ig(1.2, 0.2, delta)
ae_ig(1.2, 0.2, delta)
}


