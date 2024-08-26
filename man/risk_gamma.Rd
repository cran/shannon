\name{Truncated gamma distribution}
\alias{rlse_gamma}
\alias{rlre_gamma}
\alias{rlhce_gamma}
\alias{rlae_gamma}
\title{Relative loss for various entropy measures using the truncated gamma distribution}
\usage{
rlse_gamma(p, alpha, beta)
rlre_gamma(p, alpha, beta, delta)
rlhce_gamma(p, alpha, beta, delta)
rlae_gamma(p, alpha, beta, delta)
}
\description{
Compute the relative information loss of the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the truncated gamma distribution.
}

\arguments{
\item{alpha}{
The strictly positive shape parameter of the gamma distribution (\eqn{\alpha > 0}).
}
\item{beta}{
The strictly positive scale parameter of the gamma distribution (\eqn{\beta > 0}).
}

\item{p}{
The truncation time \eqn{(p>0)}.
}
\item{delta}{
The strictly positive parameter  (\eqn{\delta > 0}) and (\eqn{\delta \ne 1}).
}
}

\value{
The functions rlse_gamma, rlre_gamma, rlhce_gamma, and rlae_gamma provide the relative information loss based on the Shannon entropy, Rényi entropy, Havrda and Charvat entropy, and Arimoto entropy, respectively, depending on the selected parametric values of the truncated gamma distribution, \eqn{p} and \eqn{\delta}.
}
\references{
Awad, A. M., & Alawneh, A. J. (1987). Application of entropy to a life-time model. IMA Journal of Mathematical Control and Information, 4(2), 143-148.

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
\code{\link{re_gamma}
}
}


\examples{
p <- c(1, 1.50, 1.75)
rlse_gamma(p, 0.2, 1)
rlre_gamma(p, 0.2, 1, 0.5)
rlhce_gamma(p, 0.2, 1, 0.5)
rlae_gamma(p, 0.2, 1, 0.5)
}


