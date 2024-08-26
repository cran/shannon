\name{Normal distribution}
\alias{se_norm}
\alias{re_norm}
\alias{hce_norm}
\alias{ae_norm}
\title{Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the normal distribution}
\usage{
se_norm(alpha, beta)
re_norm(alpha, beta, delta)
hce_norm(alpha, beta, delta)
ae_norm(alpha, beta, delta)
}
\description{
Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the normal distribution.
}

\arguments{
\item{alpha}{
The location parameter of the normal distribution (\eqn{\alpha\in\left(-\infty,+\infty\right)}).
}
\item{beta}{
The strictly positive scale parameter of the normal distribution (\eqn{\beta > 0}).
}
\item{delta}{
The strictly positive parameter  (\eqn{\delta > 0}) and (\eqn{\delta \ne 1}).
}
}

\details{The following is the probability density function of the normal distribution:
\deqn{
f(x)=\frac{1}{\beta\sqrt{2\pi}}e^{-0.5\left(\frac{x-\alpha}{\beta}\right)^{2}},
}
where \eqn{x\in\left(-\infty,+\infty\right)}, \eqn{\alpha\in\left(-\infty,+\infty\right)} and \eqn{\beta > 0}. The parameters \eqn{\alpha} and \eqn{\beta} represent the mean and standard deviation, respectively.
}
\value{
The functions se_norm, re_norm, hce_norm, and ae_norm provide the Shannon entropy, Rényi entropy, Havrda and Charvat entropy, and Arimoto entropy, respectively, depending on the selected parametric values of the Normal distribution and \eqn{\delta}.
}
\references{
Patel, J. K., & Read, C. B. (1996). Handbook of the normal distribution (Vol. 150). CRC Press.

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
se_norm(0.2, 1.4)
delta <- c(1.5, 2, 3)
re_norm(0.2, 1.4, delta)
hce_norm(0.2, 1.4, delta)
ae_norm(0.2, 1.4, delta)
}


