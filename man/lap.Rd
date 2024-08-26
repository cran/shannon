\name{Laplace  distribution}
\alias{Se_lap}
\alias{re_lap}
\alias{hce_lap}
\alias{ae_lap}
\title{Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the Laplace or the double exponential distributiondistribution}
\usage{
Se_lap(alpha, beta)
re_lap(alpha, beta, delta)
hce_lap(alpha, beta, delta)
ae_lap(alpha, beta, delta)
}
\description{
Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the Laplace  distribution.
}

\arguments{
\item{alpha}{
The location parameter of the Laplace distribution (\eqn{\alpha\in\left(-\infty,+\infty\right)}).
}
\item{beta}{
The strictly positive scale parameter of the Laplace  distribution (\eqn{\beta > 0}).
}
\item{delta}{
The strictly positive parameter  (\eqn{\delta > 0}) and (\eqn{\delta \ne 1}).
}
}
\details{The following is the probability density function of the Laplace  distribution:
\deqn{
f(x)=\frac{1}{2\beta}e^{\frac{-|x-\alpha|}{\beta}},
}
where \eqn{x\in\left(-\infty,+\infty\right)}, \eqn{\alpha\in\left(-\infty,+\infty\right)} and \eqn{\beta > 0}.
}
\value{
The functions Se_lap, re_lap, hce_lap, and ae_lap provide the Shannon entropy, Rényi entropy, Havrda and Charvat entropy, and Arimoto entropy, respectively, depending on the selected parametric values of the Laplace  distribution and \eqn{\delta}.
}
\references{
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
\code{\link{re_gum}, \link{re_norm}
}
}


\examples{
Se_lap(0.2, 1.4)
delta <- c(2, 3)
re_lap(1.2, 0.4, delta)
hce_lap(1.2, 0.4, delta)
ae_lap(1.2, 0.4, delta)
}


