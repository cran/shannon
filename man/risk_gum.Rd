\name{Truncated Gumbel distribution}
\alias{rlse_gum}
\alias{rlre_gum}
\alias{rlhce_gum}
\alias{rlae_gum}
\title{Relative loss for various entropy measures using the truncated Gumbel distribution}
\usage{
rlse_gum(p, alpha, beta)
rlre_gum(p, alpha, beta, delta)
rlhce_gum(p, alpha, beta, delta)
rlae_gum(p, alpha, beta, delta)
}
\description{
Compute the relative information loss of the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the truncated Gumbel distribution.
}

\arguments{
\item{alpha}{
The location parameter of the Gumbel distribution (\eqn{\alpha\in\left(-\infty,+\infty\right)}).
}
\item{beta}{
The strictly positive scale parameter of the Gumbel distribution (\eqn{\beta > 0}).
}

\item{p}{
The truncation time \eqn{(p>0)}.
}
\item{delta}{
The strictly positive parameter  (\eqn{\delta > 0}) and (\eqn{\delta \ne 1}).
}
}

\value{
The functions rlse_gum, rlre_gum, rlhce_gum, and rlae_gum provide the relative information loss based on the Shannon entropy, Rényi entropy, Havrda and Charvat entropy, and Arimoto entropy, respectively, depending on the selected parametric values of the truncated Gumbel distribution, \eqn{p} and \eqn{\delta}.
}
\references{
Awad, A. M., & Alawneh, A. J. (1987). Application of entropy to a life-time model. IMA Journal of Mathematical Control and Information, 4(2), 143-148.

Gomez, Y. M., Bolfarine, H., & Gomez, H. W. (2019). Gumbel distribution with heavy tails and applications to environmental data. Mathematics and Computers in Simulation, 157, 115-129.

}

\author{
Muhammad Imran, Christophe Chesneau and Farrukh Jamal

R implementation and documentation: Muhammad Imran <imranshakoor84@yahoo.com>, Christophe Chesneau <christophe.chesneau@unicaen.fr> and Farrukh Jamal \email{farrukh.jamal@iub.edu.pk}.
}


%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{re_gum}
}
}


\examples{
p <- c(1.8,2.2)
rlse_gum(p, 4, 2)
rlre_gum(p, 4, 2, 2)
rlhce_gum(p, 4, 2, 2)
rlae_gum(p, 4, 2, 2)
}


