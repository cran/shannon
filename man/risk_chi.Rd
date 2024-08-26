\name{Truncated Chi-squared distribution}
\alias{rlse_chi}
\alias{rlre_chi}
\alias{rlhce_chi}
\alias{rlae_chi}
\title{Relative loss for various entropy measures using the truncated Chi-squared distribution}
\usage{
rlse_chi(p, n)
rlre_chi(p, n, delta)
rlhce_chi(p, n, delta)
rlae_chi(p, n, delta)
}
\description{
Compute the relative information loss of the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the truncated Chi-squared distribution.
}

\arguments{
\item{n}{
The degree of freedom and positive parameter of the Chi-squared distribution (\eqn{n > 0}).
}

\item{p}{
The truncation time \eqn{(p>0)}.
}
\item{delta}{
The strictly positive parameter  (\eqn{\delta > 0}) and (\eqn{\delta \ne 1}).
}
}

\value{
The functions rlse_chi, rlre_chi, rlhce_chi, and rlae_chi provide the relative information loss based on the Shannon entropy, Rényi entropy, Havrda and Charvat entropy, and Arimoto entropy, respectively, depending on the selected parametric values of the truncated Chi-squared distribution, \eqn{p} and \eqn{\delta}.
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
\code{\link{re_chi}
}
}


\examples{
p <- c(1, 1.7, 3)
rlse_chi(p, 2)
rlre_chi(p, 2, 0.5)
rlhce_chi(p, 2, 0.5)
rlae_chi(p, 2, 0.5)
}


