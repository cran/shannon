\name{Chi-squared distribution}
\alias{se_chi}
\alias{re_chi}
\alias{hce_chi}
\alias{ae_chi}
\title{Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the Chi-squared distribution}
\usage{
se_chi(n)
re_chi(n, delta)
hce_chi(n, delta)
ae_chi(n, delta)
}
\description{
Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the chi-squared distribution.
}

\arguments{
\item{n}{
The degree of freedom and the positive parameter of the Chi-squared distribution (\eqn{n > 0}).
}
\item{delta}{
The strictly positive parameter  (\eqn{\delta > 0}) and (\eqn{\delta \ne 1}).
}
}
\details{The following is the probability density function of the (non-central) Chi-squared distribution:
\deqn{
f(x)=\frac{1}{2^{\frac{n}{2}}\Gamma(\frac{n}{2})}x^{\frac{n}{2}-1}e^{-\frac{x}{2}},
}
where \eqn{x > 0} and \eqn{n > 0}, and \eqn{\Gamma(a)} denotes the standard gamma function.
}
\value{
The functions se_chi, re_chi, hce_chi, and ae_chi provide the Shannon entropy, Rényi entropy, Havrda and Charvat entropy, and Arimoto entropy, respectively, depending on the selected parametric values of the Chi-squared distribution and \eqn{\delta}.
}
\references{

Johnson, N. L., Kotz, S., & Balakrishnan, N. (1995). Continuous univariate distributions, volume 2 (Vol. 289). John Wiley & Sons.

}

\author{
Muhammad Imran, Christophe Chesneau and Farrukh Jamal

R implementation and documentation: Muhammad Imran <imranshakoor84@yahoo.com>, Christophe Chesneau <christophe.chesneau@unicaen.fr> and Farrukh Jamal \email{farrukh.jamal@iub.edu.pk}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{re_exp}, \link{re_gamma}, \link{re_bs}
}
}


\examples{
se_chi(1.2)
delta <- c(0.2, 0.3)
re_chi(1.2, delta)
hce_chi(1.2, delta)
ae_chi(1.2, delta)
}


