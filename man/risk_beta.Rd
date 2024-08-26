\name{Truncated beta distribution}
\alias{rlse_beta}
\alias{rlre_beta}
\alias{rlhce_beta}
\alias{rlae_beta}
\title{Relative loss for various entropy measures using the truncated beta distribution}
\usage{
rlse_beta(p, alpha, beta)
rlre_beta(p, alpha, beta, delta)
rlhce_beta(p, alpha, beta, delta)
rlae_beta(p, alpha, beta, delta)
}
\description{
Compute the relative information loss of the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the truncated beta distribution.
}

\arguments{
\item{alpha}{
The strictly positive shape parameter of the beta distribution (\eqn{\alpha > 0}).
}
\item{beta}{
The strictly positive shape parameter of the beta distribution (\eqn{\beta > 0}).
}

\item{p}{
The truncation time \eqn{(p>0)}.
}
\item{delta}{
The strictly positive parameter  (\eqn{\delta > 0}) and (\eqn{\delta \ne 1}).
}
}

\value{
The functions rlse_beta, rlre_beta, rlhce_beta, and rlae_beta provide the relative information loss based on the Shannon entropy, Rényi entropy, Havrda and Charvat entropy, and Arimoto entropy, respectively, depending on the selected parametric values of the truncated beta distribution, \eqn{p} and \eqn{\delta}.
}
\references{
Gupta, A. K., & Nadarajah, S. (2004). Handbook of beta distribution and its applications. CRC Press.

Awad, A. M., & Alawneh, A. J. (1987). Application of entropy to a life-time model. IMA Journal of Mathematical Control and Information, 4(2), 143-148.

}

\author{
Muhammad Imran, Christophe Chesneau and Farrukh Jamal

R implementation and documentation: Muhammad Imran <imranshakoor84@yahoo.com>, Christophe Chesneau <christophe.chesneau@unicaen.fr> and Farrukh Jamal \email{farrukh.jamal@iub.edu.pk}.
}


%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{re_beta}
}
}


\examples{
p <- c(0.25, 0.50, 0.75)
rlse_beta(p, 0.2, 0.4)
rlre_beta(p, 0.2, 0.4, 0.5)
rlhce_beta(p, 0.2, 0.4, 0.5)
rlae_beta(p, 0.2, 0.4, 0.5)
}


