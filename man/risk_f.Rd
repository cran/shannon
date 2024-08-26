\name{Truncated F distribution}
\alias{rlse_f}
\alias{rlre_f}
\alias{rlhce_f}
\alias{rlae_f}
\title{Relative loss for various entropy measures using the truncated F distribution}
\usage{
rlse_f(p, alpha, beta)
rlre_f(p, alpha, beta, delta)
rlhce_f(p, alpha, beta, delta)
rlae_f(p, alpha, beta, delta)
}
\description{
Compute the relative information loss of the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the truncated F distribution.
}

\arguments{
\item{alpha}{
The strictly positive parameter (first degree of freedom) of the F distribution (\eqn{\alpha > 0}).
}
\item{beta}{
The strictly positive parameter (second degree of freedom) of the F distribution (\eqn{\beta > 0}).
}

\item{p}{
The truncation time \eqn{(p>0)}.
}
\item{delta}{
The strictly positive parameter  (\eqn{\delta > 0}) and (\eqn{\delta \ne 1}).
}
}

\value{
The functions rlse_f, rlre_f, rlhce_f, and rlae_f provide the relative information loss based on the Shannon entropy, Rényi entropy, Havrda and Charvat entropy, and Arimoto entropy, respectively, depending on the selected parametric values of the truncated F distribution, \eqn{p} and \eqn{\delta}.
}
\references{
Awad, A. M., & Alawneh, A. J. (1987). Application of entropy to a life-time model. IMA Journal of Mathematical Control and Information, 4(2), 143-148.
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
\code{\link{re_f}
}
}


\examples{
p <- c(1.25, 1.50, 1.75)
rlse_f(p, 4, 6)
rlre_f(p, 4, 6, 0.5)
rlhce_f(p, 4, 6, 0.5)
rlae_f(p, 4, 6, 0.5)
}


