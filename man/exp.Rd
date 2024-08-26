\name{Exponential distribution}
\alias{Se_exp}
\alias{re_exp}
\alias{hce_exp}
\alias{ae_exp}
\title{ Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the exponential distribution}
\usage{
Se_exp(alpha)
re_exp(alpha, delta)
hce_exp(alpha, delta)
ae_exp(alpha, delta)
}
\description{
Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the exponential distribution.
}

\arguments{
\item{alpha}{
The strictly positive scale parameter of the exponential distribution (\eqn{\alpha > 0}).
}
\item{delta}{
The strictly positive parameter  (\eqn{\delta > 0}) and (\eqn{\delta \ne 1}).
}
}

\details{The following is the probability density function of the exponential distribution:
\deqn{
f(x)=\alpha e^{-\alpha x},
}
where \eqn{x > 0} and \eqn{\alpha > 0}.
}
\value{
The functions Se_exp, re_exp, hce_exp, and ae_exp provide the Shannon entropy, Rényi entropy, Havrda and Charvat entropy, and Arimoto entropy, respectively, depending on the selected parametric values of the exponential distribution and \eqn{\delta}.
}
\references{
Balakrishnan, K. (2019). Exponential distribution: theory, methods and applications. Routledge.

Singh, A. K. (1997). The exponential distribution-theory, methods and applications, Technometrics,
39(3), 341-341.

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
\code{\link{re_chi}, \link{re_gamma}, \link{re_wei}
}
}


\examples{
Se_exp(0.2)
delta <- c(1.5, 2, 3)
re_exp(0.2, delta)
hce_exp(0.2, delta)
ae_exp(0.2, delta)
}


