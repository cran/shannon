\name{Exponential extension distribution}
\alias{se_nh}
\alias{re_nh}
\alias{hce_nh}
\alias{ae_nh}
\title{Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the exponential extension distribution}
\usage{
se_nh(alpha, beta)
re_nh(alpha, beta, delta)
hce_nh(alpha, beta, delta)
ae_nh(alpha, beta, delta)
}
\description{
Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the exponential extension distribution.
}

\arguments{
\item{alpha}{
The strictly positive parameter of the exponential extension distribution (\eqn{\alpha > 0}).
}
\item{beta}{
The strictly positive parameter of the exponential extension distribution (\eqn{\beta > 0}).
}
\item{delta}{
The strictly positive parameter  (\eqn{\delta > 0}) and (\eqn{\delta \ne 1}).
}
}

\details{The following is the probability density function of the exponential extension distribution:
\deqn{
f(x)=\alpha\beta(1+\alpha x)^{\beta-1}e^{1-(1+\alpha x)^{\beta}},
}
where \eqn{x > 0}, \eqn{\alpha > 0} and \eqn{\beta > 0}.
}
\value{
The functions se_nh, re_nh, hce_nh, and ae_nh provide the Shannon entropy, Rényi entropy, Havrda and Charvat entropy, and Arimoto entropy, respectively, depending on the selected parametric values of the exponential extension distribution and \eqn{\delta}.
}
\references{
Nadarajah, S., & Haghighi, F. (2011). An extension of the exponential distribution. Statistics, 45(6), 543-558.

}

\author{
Muhammad Imran, Christophe Chesneau and Farrukh Jamal

R implementation and documentation: Muhammad Imran <imranshakoor84@yahoo.com>, Christophe Chesneau <christophe.chesneau@unicaen.fr> and Farrukh Jamal \email{farrukh.jamal@iub.edu.pk}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{re_exp}, \link{re_gamma}, \link{re_ee}, \link{re_wei}
}
}


\examples{
se_nh(1.2, 0.2)
delta <- c(1.5, 2, 3)
re_nh(1.2, 0.2, delta)
hce_nh(1.2, 0.2, delta)
ae_nh(1.2, 0.2, delta)
}


