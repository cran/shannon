\name{Truncated Weibull distribution}
\alias{rlse_wei}
\alias{rlre_wei}
\alias{rlhce_wei}
\alias{rlae_wei}
\title{Relative loss for various entropy measures using the truncated Weibull distribution}
\usage{
rlse_wei(p, alpha, beta)
rlre_wei(p, alpha, beta, delta)
rlhce_wei(p, alpha, beta, delta)
rlae_wei(p, alpha, beta, delta)
}
\description{
Compute the relative information loss of the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the truncated Weibull distribution.
}

\arguments{
\item{alpha}{
The strictly positive scale parameter of the Weibull distribution (\eqn{\alpha > 0}).
}
\item{beta}{
The strictly positive shape parameter of the Weibull distribution (\eqn{\beta > 0}).
}

\item{p}{
The truncation time \eqn{(p>0)}.
}
\item{delta}{
The strictly positive parameter  (\eqn{\delta > 0}) and (\eqn{\delta \ne 1}).
}
}

\value{
The functions rlse_wei, rlre_wei, rlhce_wei, and rlae_wei provide the relative information loss based on the Shannon entropy, Rényi entropy, Havrda and Charvat entropy, and Arimoto entropy, respectively, depending on the selected parametric values of the truncated Weibull distribution, \eqn{p} and \eqn{\delta}.
}
\references{
Awad, A. M., & Alawneh, A. J. (1987). Application of entropy to a life-time model. IMA Journal of Mathematical Control and Information, 4(2), 143-148.

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
\code{\link{re_wei}
}
}


\examples{
p <- c(1, 1.7, 3)
rlse_wei(p, 2, 1)
rlre_wei(p, 2, 1, 0.5)
rlhce_wei(p, 2, 1, 0.5)
rlae_wei(p, 2, 1, 0.5)
}


