\name{Truncated exponential distribution}
\alias{rlse_exp}
\alias{rlre_exp}
\alias{rlhce_exp}
\alias{rlae_exp}
\title{Relative loss for various entropy measures using the truncated exponential distribution}
\usage{
rlse_exp(p, alpha)
rlre_exp(p, alpha, delta)
rlhce_exp(p, alpha, delta)
rlae_exp(p, alpha, delta)
}
\description{
Compute the relative information loss of the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the truncated exponential distribution.
}

\arguments{
\item{alpha}{
The strictly positive scale parameter of the exponential distribution (\eqn{\alpha > 0}).
}

\item{p}{
The truncation time \eqn{(p>0)}.
}
\item{delta}{
The strictly positive parameter  (\eqn{\delta > 0}) and (\eqn{\delta \ne 1}).
}
}

\value{
The functions rlse_exp, rlre_exp, rlhce_exp, and rlae_exp provide the relative information loss based on the Shannon entropy, Rényi entropy, Havrda and Charvat entropy, and Arimoto entropy, respectively, depending on the selected parametric values of the truncated exponential distribution, \eqn{p} and \eqn{\delta}.
}
\references{
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
\code{\link{re_exp}
}
}


\examples{
p <- c(1, 1.7, 3)
rlse_exp(p, 2)
rlre_exp(p, 2, 0.5)
rlhce_exp(p, 2, 0.5)
rlae_exp(p, 2, 0.5)
}


