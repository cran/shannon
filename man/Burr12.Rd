\name{Burr XII distribution}
\alias{se_burr}
\alias{re_burr}
\alias{hce_burr}
\alias{ae_burr}
\title{Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the Burr XII distribution}
\usage{
se_burr(k, c)
re_burr(k, c, delta)
hce_burr(k, c, delta)
ae_burr(k, c, delta)
}
\description{
Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the Burr XII distribution.
}

\arguments{
\item{k}{
The strictly positive shape parameter of the Burr XII distribution (\eqn{k > 0}).
}
\item{c}{
The strictly positive shape parameter of the Burr XII distribution (\eqn{c > 0}).
}
\item{delta}{
The strictly positive parameter  (\eqn{\delta > 0}) and (\eqn{\delta \ne 1}).
}
}

\details{The following is the probability density function of the Burr XII distribution:
\deqn{
f(x)=kcx^{c-1}\left(1+x^{c}\right)^{-k-1},
}
where \eqn{x > 0}, \eqn{c > 0} and \eqn{k > 0}.
}
\value{
The functions se_burr, re_burr, hce_burr, and ae_burr provide the Shannon entropy, Rényi entropy, Havrda and Charvat entropy, and Arimoto entropy, respectively, depending on the selected parametric values of the Burr XII distribution and \eqn{\delta}.
}
\references{
Rodriguez, R. N. (1977). A guide to the Burr type XII distributions. Biometrika, 64(1), 129-134.

Zimmer, W. J., Keats, J. B., & Wang, F. K. (1998). The Burr XII distribution in reliability analysis. Journal of Quality Technology, 30(4), 386-394.
}

\author{
Muhammad Imran, Christophe Chesneau and Farrukh Jamal

R implementation and documentation: Muhammad Imran <imranshakoor84@yahoo.com>, Christophe Chesneau <christophe.chesneau@unicaen.fr> and Farrukh Jamal \email{farrukh.jamal@iub.edu.pk}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{re_gamma}, \link{re_wei}
}
}


\examples{
se_burr(0.2, 1.4)
delta <- c(2, 3)
re_burr(1.2, 1.4, delta)
hce_burr(1.2, 1.4, delta)
ae_burr(1.2, 1.4, delta)
}


