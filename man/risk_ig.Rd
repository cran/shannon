\name{Truncated inverse-gamma distribution}
\alias{rlse_ig}
\alias{rlre_ig}
\alias{rlhce_ig}
\alias{rlae_ig}
\title{Relative loss for various entropy measures using the truncated inverse-gamma distribution}
\usage{
rlse_ig(p, alpha, beta)
rlre_ig(p, alpha, beta, delta)
rlhce_ig(p, alpha, beta, delta)
rlae_ig(p, alpha, beta, delta)
}
\description{
Compute the relative information loss of the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the truncated inverse-gamma distribution.
}

\arguments{
\item{alpha}{
The strictly positive shape parameter of the inverse-gamma  distribution (\eqn{\alpha > 0}).
}
\item{beta}{
The strictly positive scale parameter of the inverse-gamma  distribution (\eqn{\beta > 0}).
}
\item{p}{
The truncation time \eqn{(p>0)}.
}
\item{delta}{
The strictly positive parameter  (\eqn{\delta > 0}) and (\eqn{\delta \ne 1}).
}
}

\value{
The functions rlse_ig, rlre_ig, rlhce_ig, and rlae_ig provide the relative information loss based on the Shannon entropy, Rényi entropy, Havrda and Charvat entropy, and Arimoto entropy, respectively, depending on the selected parametric values of the truncated inverse-gamma distribution, \eqn{p} and \eqn{\delta}.
}
\references{
Awad, A. M., & Alawneh, A. J. (1987). Application of entropy to a life-time model. IMA Journal of Mathematical Control and Information, 4(2), 143-148.

Rivera, P. A., Calderín-Ojeda, E., Gallardo, D. I., & Gómez, H. W. (2021). A compound class of the inverse Gamma and power series distributions. Symmetry, 13(8), 1328.

}

\author{
Muhammad Imran, Christophe Chesneau and Farrukh Jamal

R implementation and documentation: Muhammad Imran <imranshakoor84@yahoo.com>, Christophe Chesneau <christophe.chesneau@unicaen.fr> and Farrukh Jamal \email{farrukh.jamal@iub.edu.pk}.
}


%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{re_ig}
}
}


\examples{
p <- c(1.25, 1.50)
rlse_ig(p, 1.2, 0.2)
rlre_ig(p, 1.2, 0.2, 0.5)
rlhce_ig(p, 1.2, 0.2, 0.5)
rlae_ig(p, 1.2, 0.2, 0.5)
}


