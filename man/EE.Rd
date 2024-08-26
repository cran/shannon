\name{Exponentiated exponential distribution}
\alias{se_ee}
\alias{re_ee}
\alias{hce_ee}
\alias{ae_ee}
\title{Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the exponentiated exponential distribution}
\usage{
se_ee(alpha, beta)
re_ee(alpha, beta, delta)
hce_ee(alpha, beta, delta)
ae_ee(alpha, beta, delta)
}
\description{
Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the exponentiated exponential distribution.
}

\arguments{
\item{alpha}{
The strictly positive scale parameter of the exponentiated exponential distribution (\eqn{\alpha > 0}).
}
\item{beta}{
The strictly positive shape parameter of the exponentiated exponential distribution (\eqn{\beta > 0}).
}
\item{delta}{
The strictly positive parameter  (\eqn{\delta > 0}) and (\eqn{\delta \ne 1}).
}
}

\details{The following is the probability density function of the  exponentiated exponential distribution:
\deqn{
f(x)=\alpha\beta e^{-\alpha x}\left(1-e^{-\alpha x}\right)^{\beta-1},
}
where \eqn{x > 0}, \eqn{\alpha > 0} and \eqn{\beta > 0}.
}
\value{
The functions se_ee, re_ee, hce_ee, and ae_ee provide the Shannon entropy, Rényi entropy, Havrda and Charvat entropy, and Arimoto entropy, respectively, depending on the selected parametric values of the exponentiated exponential distribution and \eqn{\delta}.
}
\references{
Nadarajah, S. (2011). The exponentiated exponential distribution: a survey. AStA Advances in Statistical Analysis, 95, 219-251.

Gupta, R. D., & Kundu, D. (2007). Generalized exponential distribution: Existing results and some recent developments. Journal of Statistical Planning and Inference, 137(11), 3537-3547.

}
\author{
Muhammad Imran, Christophe Chesneau and Farrukh Jamal

R implementation and documentation: Muhammad Imran <imranshakoor84@yahoo.com>, Christophe Chesneau <christophe.chesneau@unicaen.fr> and Farrukh Jamal \email{farrukh.jamal@iub.edu.pk}.
}


%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{re_exp}, \link{re_wei}, \link{re_nh}
}
}


\examples{
se_ee(0.2, 1.4)
delta <- c(1.5, 2, 3)
re_ee(0.2, 1.4, delta)
hce_ee(0.2, 1.4, delta)
ae_ee(0.2, 1.4, delta)
}


