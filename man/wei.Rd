\name{Weibull distribution}
\alias{se_wei}
\alias{re_wei}
\alias{hce_wei}
\alias{ae_wei}
\title{Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the Weibull distribution}
\usage{
se_wei(alpha, beta)
re_wei(alpha, beta, delta)
hce_wei(alpha, beta, delta)
ae_wei(alpha, beta, delta)
}
\description{
Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the Weibull distribution.
}

\arguments{
\item{alpha}{
The strictly positive scale parameter of the Weibull distribution (\eqn{\alpha > 0}).
}
\item{beta}{
The strictly positive shape parameter of the Weibull distribution (\eqn{\beta > 0}).
}
\item{delta}{
The strictly positive parameter  (\eqn{\delta > 0}) and (\eqn{\delta \ne 1}).
}
}

\details{The following is the probability density function of the Weibull distribution:
\deqn{
f(x)=\frac{\beta}{\alpha}\left(\frac{x}{\alpha}\right)^{\beta-1}e^{-(\frac{x}{\alpha})^{\beta}},
}
where \eqn{x > 0}, \eqn{\alpha > 0} and \eqn{\beta > 0}.
}
\value{
The functions se_wei, re_wei, hce_wei, and ae_wei provide the Shannon entropy, Rényi entropy, Havrda and Charvat entropy, and Arimoto entropy, respectively, depending on the selected parametric values of the Weibull distribution and \eqn{\delta}.
}
\references{
Weibull, W. (1951). A statistical distribution function of wide applicability. Journal of applied mechanics, 18, 293-297.

}

\author{
Muhammad Imran, Christophe Chesneau and Farrukh Jamal

R implementation and documentation: Muhammad Imran <imranshakoor84@yahoo.com>, Christophe Chesneau <christophe.chesneau@unicaen.fr> and Farrukh Jamal \email{farrukh.jamal@iub.edu.pk}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{re_exp}, \link{re_gamma}, \link{re_ee}
}
}


\examples{
se_wei(1.2, 0.2)
delta <- c(1.5, 2, 3)
re_wei(1.2, 0.2, delta)
hce_wei(1.2, 0.2, delta)
ae_wei(1.2, 0.2, delta)
}


