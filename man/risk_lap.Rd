\name{Truncated Laplace  distribution}
\alias{rlse_lap}
\alias{rlre_lap}
\alias{rlhce_lap}
\alias{rlae_lap}
\title{Relative loss for various entropy measures using the truncated Laplace  distribution}
\usage{
rlse_lap(p, alpha, beta)
rlre_lap(p, alpha, beta, delta)
rlhce_lap(p, alpha, beta, delta)
rlae_lap(p, alpha, beta, delta)
}
\description{
Compute the relative information loss of the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the truncated Laplace  distribution.
}

\arguments{
\item{alpha}{
Location parameter of the Laplace distribution (\eqn{\alpha\in\left(-\infty,+\infty\right)}).
}
\item{beta}{
The strictly positive scale parameter of the Laplace  distribution (\eqn{\beta > 0}).
}

\item{p}{
The truncation time \eqn{(p>0)}.
}
\item{delta}{
The strictly positive parameter  (\eqn{\delta > 0}) and (\eqn{\delta \ne 1}).
}
}

\value{
The functions rlse_lap, rlre_lap, rlhce_lap, and rlae_lap provide the relative information loss based on the Shannon entropy, Rényi entropy, Havrda and Charvat entropy, and Arimoto entropy, respectively, depending on the selected parametric values of the truncated Laplace  distribution, \eqn{p} and \eqn{\delta}.
}
\references{
Awad, A. M., & Alawneh, A. J. (1987). Application of entropy to a life-time model. IMA Journal of Mathematical Control and Information, 4(2), 143-148.

Cordeiro, G. M., & Lemonte, A. J. (2011). The beta Laplace distribution. Statistics & Probability Letters, 81(8), 973-982.

}

\author{
Muhammad Imran, Christophe Chesneau and Farrukh Jamal

R implementation and documentation: Muhammad Imran <imranshakoor84@yahoo.com>, Christophe Chesneau <christophe.chesneau@unicaen.fr> and Farrukh Jamal \email{farrukh.jamal@iub.edu.pk}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{re_lap}
}
}


\examples{
p <- c(0.25, 0.50, 0.75)
rlse_lap(p, 0.2, 0.4)
rlre_lap(p, 0.2, 0.4, 0.5)
rlhce_lap(p, 0.2, 0.4, 0.5)
rlae_lap(p, 0.2, 0.4, 0.5)
}


