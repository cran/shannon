\name{Gompertz distribution}
\alias{se_gomp}
\alias{re_gomp}
\alias{hce_gomp}
\alias{ae_gomp}
\title{Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the Gompertz distribution}
\usage{
se_gomp(alpha, beta)
re_gomp(alpha, beta, delta)
hce_gomp(alpha, beta, delta)
ae_gomp(alpha, beta, delta)
}
\description{
Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the Gompertz distribution.
}

\arguments{
\item{alpha}{
The strictly positive parameter of the Gompertz  distribution (\eqn{\alpha > 0}).
}
\item{beta}{
The strictly positive parameter of the Gompertz  distribution (\eqn{\beta > 0}).
}
\item{delta}{
The strictly positive parameter  (\eqn{\delta > 0}) and (\eqn{\delta \ne 1}).
}
}
\details{The following is the probability density function of the Gompertz distribution:
\deqn{
f(x)=\alpha e^{\beta x-\frac{\alpha}{\beta}\left(e^{\beta x}-1\right)},
}
where \eqn{x > 0}, \eqn{\alpha > 0} and \eqn{\beta > 0}.
}
\value{
The functions se_gomp, re_gomp, hce_gomp, and ae_gomp provide the Shannon entropy, Rényi entropy, Havrda and Charvat entropy, and Arimoto entropy, respectively, depending on the selected parametric values of the Gompertz distribution and \eqn{\delta}.
}
\references{
Soliman, A. A., Abd-Ellah, A. H., Abou-Elheggag, N. A., & Abd-Elmougod, G. A. (2012). Estimation of the parameters of life for Gompertz distribution using progressive first-failure censored data. Computational Statistics & Data Analysis, 56(8), 2471-2485.


}
\author{
Muhammad Imran, Christophe Chesneau and Farrukh Jamal

R implementation and documentation: Muhammad Imran <imranshakoor84@yahoo.com>, Christophe Chesneau <christophe.chesneau@unicaen.fr> and Farrukh Jamal \email{farrukh.jamal@iub.edu.pk}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{re_exp}, \link{re_gamma}, \link{re_ray}
}
}


\examples{
se_gomp(2.4,0.2)
delta <- c(2, 3)
re_gomp(2.4,0.2, delta)
hce_gomp(2.4,0.2, delta)
ae_gomp(2.4,0.2, delta)
}


