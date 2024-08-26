\name{Kumaraswamy distribution}
\alias{se_kum}
\alias{re_kum}
\alias{hce_kum}
\alias{ae_kum}
\title{Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the Kumaraswamy distribution}
\usage{
se_kum(alpha, beta)
re_kum(alpha, beta, delta)
hce_kum(alpha, beta, delta)
ae_kum(alpha, beta, delta)
}
\description{
Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the Kumaraswamy distribution.
}

\arguments{
\item{alpha}{
The strictly positive shape parameter of the Kumaraswamy distribution (\eqn{\alpha > 0}).
}
\item{beta}{
The strictly positive shape parameter of the Kumaraswamy distribution (\eqn{\beta > 0}).
}
\item{delta}{
The strictly positive scale parameter  (\eqn{\delta > 0}).
}
}

\details{ The following is the probability density function of the Kumaraswamy distribution:
\deqn{
f(x)=\alpha\beta x^{\alpha-1}\left(1-x^{\alpha}\right)^{\beta-1},
}
where \eqn{0\leq x\leq1}, \eqn{\alpha > 0} and \eqn{\beta > 0}.
}
\value{
The functions se_kum, re_kum, hce_kum, and ae_kum provide the Shannon entropy, Rényi entropy, Havrda and Charvat entropy, and Arimoto entropy, respectively, depending on the selected parametric values of the Kumaraswamy distribution and \eqn{\delta}.
}
\references{
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
\code{\link{re_beta}
}
}


\examples{
se_kum(1.2, 1.4)
delta <- c(1.5, 2, 3)
re_kum(1.2, 1.4, delta)
hce_kum(1.2, 1.4, delta)
ae_kum(1.2, 1.4, delta)
}


