\name{F distribution}
\alias{se_f}
\alias{re_f}
\alias{hce_f}
\alias{ae_f}
\title{Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the F distribution}
\usage{
se_f(alpha, beta)
re_f(alpha, beta, delta)
hce_f(alpha, beta, delta)
ae_f(alpha, beta, delta)
}
\description{
Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the F distribution.
}

\arguments{
\item{alpha}{
The strictly positive parameter (first degree of freedom) of the F distribution (\eqn{\alpha > 0}).
}
\item{beta}{
The strictly positive parameter (second degree of freedom) of the F distribution (\eqn{\beta > 0}).
}
\item{delta}{
The strictly positive parameter  (\eqn{\delta > 0}) and (\eqn{\delta \ne 1}).
}
}

\details{The following is the probability density function of the F distribution:
\deqn{
f(x)=\frac{1}{B(\frac{\alpha}{2},\frac{\beta}{2})}\left(\frac{\alpha}{\beta}\right)^{\frac{\alpha}{2}}x^{\frac{\alpha}{2}-1}\left(1+\frac{\alpha}{\beta}x\right)^{-\left(\frac{\alpha+\beta}{2}\right)},
}
where \eqn{x > 0}, \eqn{\alpha > 0} and \eqn{\beta > 0}, and \eqn{B(a,b)} is the standard beta function.
}
\value{
The functions se_f, re_f, hce_f, and ae_f provide the Shannon entropy, Rényi entropy, Havrda and Charvat entropy, and Arimoto entropy, respectively, depending on the selected parametric values of the F distribution and \eqn{\delta}.
}
\references{
Johnson, N. L., Kotz, S., & Balakrishnan, N. (1995). Continuous univariate distributions, volume 2 (Vol. 289). John Wiley & Sons.
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
se_f(1.2, 1.4)
delta <- c(2.2, 2.3)
re_f(1.2, 0.4, delta)
hce_f(1.2, 1.4, delta)
ae_f(1.2, 1.4, delta)
}


