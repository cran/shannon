\name{Truncated Gompertz distribution}
\alias{rlse_gomp}
\alias{rlre_gomp}
\alias{rlhce_gomp}
\alias{rlae_gomp}
\title{Relative loss for various entropy measures using the truncated Gompertz distribution}
\usage{
rlse_gomp(p, alpha, beta)
rlre_gomp(p, alpha, beta, delta)
rlhce_gomp(p, alpha, beta, delta)
rlae_gomp(p, alpha, beta, delta)
}
\description{
Compute the relative information loss of the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the truncated Gompertz distribution.
}

\arguments{
\item{alpha}{
The strictly positive parameter of the Gompertz  distribution (\eqn{\alpha > 0}).
}
\item{beta}{
The strictly positive parameter of the Gompertz  distribution (\eqn{\beta > 0}).
}
\item{p}{
The truncation time \eqn{(p>0)}.
}
\item{delta}{
The strictly positive parameter  (\eqn{\delta > 0}) and (\eqn{\delta \ne 1}).
}
}

\value{
The functions rlse_gomp, rlre_gomp, rlhce_gomp, and rlae_gomp provide the relative information loss based on the Shannon entropy, Rényi entropy, Havrda and Charvat entropy, and Arimoto entropy, respectively, depending on the selected parametric values of the truncated Gompertz distribution, \eqn{p} and \eqn{\delta}.
}
\references{
Awad, A. M., & Alawneh, A. J. (1987). Application of entropy to a life-time model. IMA Journal of Mathematical Control and Information, 4(2), 143-148.

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
\code{\link{re_gomp}
}
}


\examples{
p <- c(0.25, 0.50)
rlse_gomp(p, 2.4,0.2)
rlre_gomp(p, 2.4,0.2, 0.5)
rlhce_gomp(p, 2.4,0.2, 0.5)
rlae_gomp(p, 2.4,0.2, 0.5)
}


