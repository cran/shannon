\name{Frechet distribution}
\alias{se_fre}
\alias{re_fre}
\alias{hce_fre}
\alias{ae_fre}
\title{Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the Fréchet distribution}
\usage{
se_fre(alpha, beta, zeta)
re_fre(alpha, beta, zeta, delta)
hce_fre(alpha, beta, zeta, delta)
ae_fre(alpha, beta, zeta, delta)
}
\description{
Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the Fréchet distribution.
}

\arguments{
\item{alpha}{
The parameter of the Fréchet  distribution (\eqn{\alpha>0}).
}
\item{beta}{
The parameter of the Fréchet  distribution (\eqn{\beta\in\left(-\infty,+\infty\right)}).
}
\item{zeta}{
The parameter of the Fréchet  distribution (\eqn{\zeta>0}).
}
\item{delta}{
The strictly positive parameter  (\eqn{\delta > 0}) and (\eqn{\delta \ne 1}).
}
}
\details{The following is the probability density function of the Fréchet distribution:

\deqn{
f(x)=\frac{\alpha}{\zeta}\left(\frac{x-\beta}{\zeta}\right)^{-1-\alpha}e^{-(\frac{x-\beta}{\zeta})^{-\alpha},}
}
where \eqn{x>\beta}, \eqn{\alpha>0}, \eqn{\zeta>0} and \eqn{\beta\in\left(-\infty,+\infty\right)}. The Fréchet distribution is also known as inverse Weibull distribution and special case of the generalized extreme value distribution.
}
\value{
The functions se_fre, re_fre, hce_fre, and ae_fre provide the Shannon entropy, Rényi entropy, Havrda and Charvat entropy, and Arimoto entropy, respectively, depending on the selected parametric values of the Fréchet distribution distribution and \eqn{\delta}.
}
\references{
Abbas, K., & Tang, Y. (2015). Analysis of Fréchet distribution using reference priors. Communications in Statistics-Theory and Methods, 44(14), 2945-2956.
}
\author{
Muhammad Imran, Christophe Chesneau and Farrukh Jamal

R implementation and documentation: Muhammad Imran <imranshakoor84@yahoo.com>, Christophe Chesneau <christophe.chesneau@unicaen.fr> and Farrukh Jamal \email{farrukh.jamal@iub.edu.pk}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{re_exp}, \link{re_gum}
}
}


\examples{
se_fre(0.2, 1.4, 1.2)
delta <- c(2, 3)
re_fre(1.2, 0.4, 1.2, delta)
hce_fre(1.2, 0.4, 1.2, delta)
ae_fre(1.2, 0.4, 1.2, delta)
}


