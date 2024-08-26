\name{Truncated Birnbaum-Saunders distribution}
\alias{rlse_bs}
\alias{rlre_bs}
\alias{rlhce_bs}
\alias{rlae_bs}
\title{Relative loss for various entropy measures using the truncated Birnbaum-Saunders distribution}
\usage{
rlse_bs(p, v)
rlre_bs(p, v, delta)
rlhce_bs(p, v, delta)
rlae_bs(p, v, delta)
}
\description{
Compute the relative information loss of the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the truncated Birnbaum-Saunders distribution.
}

\arguments{
\item{v}{
The strictly positive scale parameter of the Birnbaum-Saunders distribution (\eqn{v > 0}).
}
\item{p}{
The truncation time \eqn{(p>0)}.
}
\item{delta}{
The strictly positive parameter  (\eqn{\delta > 0}) and (\eqn{\delta \ne 1}).
}
}

\value{
The functions rlse_bs, rlre_bs, rlhce_bs, and rlae_bs provide the relative information loss based on the Shannon entropy, Rényi entropy, Havrda and Charvat entropy, and Arimoto entropy, respectively, depending on the selected parametric values of the truncated Birnbaum-Saunders distribution, \eqn{p} and \eqn{\delta}.
}
\references{
Awad, A. M., & Alawneh, A. J. (1987). Application of entropy to a life-time model. IMA Journal of Mathematical Control and Information, 4(2), 143-148.

Chan, S., Nadarajah, S., & Afuecheta, E. (2016). An R package for value at risk and expected shortfall. Communications in Statistics Simulation and Computation, 45(9), 3416-3434.

Arimoto, S. (1971). Information-theoretical considerations on estimation problems. Inf. Control, 19, 181–194.

}

\author{
Muhammad Imran, Christophe Chesneau and Farrukh Jamal

R implementation and documentation: Muhammad Imran <imranshakoor84@yahoo.com>, Christophe Chesneau <christophe.chesneau@unicaen.fr> and Farrukh Jamal \email{farrukh.jamal@iub.edu.pk}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{re_bs}
}
}


\examples{
p <- c(1, 1.7, 3)
rlse_bs(p, 0.2)
rlre_bs(p, 0.2, 0.5)
rlhce_bs(p, 0.2, 0.5)
rlae_bs(p, 0.2, 0.5)
}


