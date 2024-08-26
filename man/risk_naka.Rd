\name{Truncated Nakagami distribution}
\alias{rlse_naka}
\alias{rlre_naka}
\alias{rlhce_naka}
\alias{rlae_naka}
\title{Relative loss for various entropy measures using the truncated Nakagami distribution}
\usage{
rlse_naka(p, alpha, beta)
rlre_naka(p, alpha, beta, delta)
rlhce_naka(p, alpha, beta, delta)
rlae_naka(p, alpha, beta, delta)
}
\description{
Compute the relative information loss of the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the truncated Nakagami distribution.
}

\arguments{
\item{alpha}{
The strictly positive scale parameter of the Nakagami distribution (\eqn{\alpha > 0}).
}
\item{beta}{
The strictly positive shape parameter of the Nakagami distribution (\eqn{\beta > 0}).
}

\item{p}{
The truncation time \eqn{(p>0)}.
}
\item{delta}{
The strictly positive parameter  (\eqn{\delta > 0}) and (\eqn{\delta \ne 1}).
}
}

\value{
The functions rlse_naka, rlre_naka, rlhce_naka, and rlae_naka provide the relative information loss based on the Shannon entropy, Rényi entropy, Havrda and Charvat entropy, and Arimoto entropy, respectively, depending on the selected parametric values of the truncated Nakagami distribution, \eqn{p} and \eqn{\delta}.
}
\references{
Awad, A. M., & Alawneh, A. J. (1987). Application of entropy to a life-time model. IMA Journal of Mathematical Control and Information, 4(2), 143-148.

Schwartz, J., Godwin, R. T., & Giles, D. E. (2013). Improved maximum-likelihood estimation of the shape parameter in the Nakagami distribution. Journal of Statistical Computation and Simulation, 83(3), 434-445.


}

\author{
Muhammad Imran, Christophe Chesneau and Farrukh Jamal

R implementation and documentation: Muhammad Imran <imranshakoor84@yahoo.com>, Christophe Chesneau <christophe.chesneau@unicaen.fr> and Farrukh Jamal \email{farrukh.jamal@iub.edu.pk}.
}


%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{re_naka}
}
}


\examples{
p <- c(1.25, 1.50, 1.75)
rlse_naka(p, 0.2, 1)
rlre_naka(p, 0.2, 1, 0.5)
rlhce_naka(p, 0.2, 1, 0.5)
rlae_naka(p, 0.2, 1, 0.5)
}


