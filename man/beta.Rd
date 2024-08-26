\name{Beta distribution}
\alias{Se_beta}
\alias{re_beta}
\alias{hce_beta}
\alias{ae_beta}
\title{Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the beta distribution}
\usage{
Se_beta(alpha, beta)
re_beta(alpha, beta, delta)
hce_beta(alpha, beta, delta)
ae_beta(alpha, beta, delta)
}
\description{
Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the beta distribution.
}
\arguments{
\item{alpha}{
The strictly positive shape parameter of the beta distribution (\eqn{\alpha > 0}).
}
\item{beta}{
The strictly positive shape parameter of the beta distribution (\eqn{\beta > 0}).
}
\item{delta}{
The strictly positive parameter  (\eqn{\delta > 0}) and (\eqn{\delta \ne 1}).
}
}

\details{The following is the probability density function of the beta distribution:
\deqn{
f(x)=\frac{\Gamma\left(\alpha+\beta\right)}{\Gamma\left(\alpha\right)\Gamma\left(\beta\right)}x^{\alpha-1}\left(1-x\right)^{\beta-1},
}
where \eqn{0\leq x\leq1}, \eqn{\alpha > 0} and \eqn{\beta > 0}, and \eqn{\Gamma(a)} denotes the standard gamma function.

}
\value{
The functions Se_beta, re_beta, hce_beta, and ae_beta provide the Shannon entropy, Rényi entropy, Havrda and Charvat entropy, and Arimoto entropy, respectively, depending on the selected parametric values of the beta distribution and \eqn{\delta}.
}
\references{
Gupta, A. K., & Nadarajah, S. (2004). Handbook of beta distribution and its applications. CRC Press.

Johnson, N. L., Kotz, S., & Balakrishnan, N. (1994). Beta distributions. Continuous univariate distributions. 2nd ed. New York, NY: John Wiley and Sons, 221-235.
}
\author{
Muhammad Imran, Christophe Chesneau and Farrukh Jamal

R implementation and documentation: Muhammad Imran <imranshakoor84@yahoo.com>, Christophe Chesneau <christophe.chesneau@unicaen.fr> and Farrukh Jamal \email{farrukh.jamal@iub.edu.pk}.
}


%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{se_kum}, \link{re_kum}, \link{hce_kum}, \link{ae_kum}
}
}


\examples{
# Computation of the Shannon entropy
Se_beta(2, 4)
delta <- c(1.2, 3)
# Computation of the Rényi  entropy
re_beta(2, 4, delta)
# Computation of the Havrda and Charvat entropy
hce_beta(2, 4, delta)
# Computation of the Arimoto entropy
ae_beta(2, 4, delta)
# A graphic presentation of the Havrda and Charvat entropy (HCE)
library(ggplot2)
delta <- c(0.2, 0.3, 0.5, 0.8, 1.2, 1.5, 2.5, 3, 3.5)
hce_beta(2, 1.2, delta)
z <- hce_beta(2, 1.2, delta)
dat <- data.frame(x = delta , HCE = z)
p_hce <- ggplot(dat, aes(x = delta, y = HCE)) + geom_line()
plot <- p_hce + ggtitle(expression(alpha == 2~~beta == 1.2))
}


