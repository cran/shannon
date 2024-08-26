\name{Truncated Kumaraswamy distribution}
\alias{rlse_kum}
\alias{rlre_kum}
\alias{rlhce_kum}
\alias{rlae_kum}
\title{Relative loss for various entropy measures using the truncated Kumaraswamy distribution}
\usage{
rlse_kum(p, alpha, beta)
rlre_kum(p, alpha, beta, delta)
rlhce_kum(p, alpha, beta, delta)
rlae_kum(p, alpha, beta, delta)
}
\description{
Compute the relative information loss of the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the truncated Kumaraswamy distribution.
}

\arguments{
\item{alpha}{
The strictly positive shape parameter of the Kumaraswamy distribution (\eqn{\alpha > 0}).
}
\item{beta}{
The strictly positive shape parameter of the Kumaraswamy distribution (\eqn{\beta > 0}).
}

\item{p}{
The truncation time \eqn{(p>0)}.
}
\item{delta}{
The strictly positive parameter  (\eqn{\delta > 0}) and (\eqn{\delta \ne 1}).
}
}

\value{
The functions rlse_kum, rlre_kum, rlhce_kum, and rlae_kum provide the relative information loss based on the Shannon entropy, Rényi entropy, Havrda and Charvat entropy, and Arimoto entropy, respectively, depending on the selected parametric values of the truncated Kumaraswamy distribution, \eqn{p} and \eqn{\delta}.
}
\references{
Awad, A. M., & Alawneh, A. J. (1987). Application of entropy to a life-time model. IMA Journal of Mathematical Control and Information, 4(2), 143-148.

El-Sherpieny, E. S. A., & Ahmed, M. A. (2014). On the Kumaraswamy distribution. International Journal of Basic and Applied Sciences, 3(4), 372.

Al-Babtain, A. A., Elbatal, I., Chesneau, C., & Elgarhy, M. (2021). Estimation of different types of entropies for the Kumaraswamy distribution. PLoS One, 16(3), e0249027.

}

\author{
Muhammad Imran, Christophe Chesneau and Farrukh Jamal

R implementation and documentation: Muhammad Imran <imranshakoor84@yahoo.com>, Christophe Chesneau <christophe.chesneau@unicaen.fr> and Farrukh Jamal \email{farrukh.jamal@iub.edu.pk}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{re_kum}
}
}


\examples{
p <- c(0.25, 0.50, 0.75)
rlse_kum(p, 0.2, 0.4)
rlre_kum(p, 0.2, 0.4, 0.5)
rlhce_kum(p, 0.2, 0.4, 0.5)
rlae_kum(p, 0.2, 0.4, 0.5)
}


