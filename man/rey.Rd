\name{Rayleigh  distribution}
\alias{se_ray}
\alias{re_ray}
\alias{hce_ray}
\alias{ae_ray}
\title{ Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the Rayleigh  distribution}
\usage{
se_ray(alpha)
re_ray(alpha, delta)
hce_ray(alpha, delta)
ae_ray(alpha, delta)
}
\description{
Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the Rayleigh  distribution.
}

\arguments{
\item{alpha}{
The strictly positive parameter of the Rayleigh distribution (\eqn{\alpha > 0}).
}

\item{delta}{
The strictly positive parameter  (\eqn{\delta > 0}) and (\eqn{\delta \ne 1}).
}
}

\details{The following is the probability density function of the Rayleigh  distribution:
\deqn{
f(x)=2\alpha xe^{-\alpha x^{2}},
}
where \eqn{x > 0} and \eqn{\alpha > 0}.
}
\value{
The functions se_ray, re_ray, hce_ray, and ae_ray provide the Shannon entropy, Rényi entropy, Havrda and Charvat entropy, and Arimoto entropy, respectively, depending on the selected parametric values of the Rayleigh  distribution and \eqn{\delta}.
}
\references{


Dey, S., Maiti, S. S., & Ahmad, M. (2016). Comparison of different entropy measures. Pak. J. Statist, 32(2), 97-108.

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
\code{\link{re_exp}, \link{re_gamma}, \link{re_wei}
}
}


\examples{
se_ray(0.2)
delta <- c(1.5, 2, 3)
re_ray(0.2, delta)
hce_ray(0.2, delta)
ae_ray(0.2, delta)
# A graphic representation of the Rényi entropy (RE)
library(ggplot2)
delta <- c(1.5, 2, 3)
z <- re_ray(0.2, delta)
dat <- data.frame(x = delta , RE = z)
p_re <- ggplot(dat, aes(x = delta, y = RE)) + geom_line()
plot <- p_re + ggtitle(expression(alpha == 0.2))
# A graphic presentation of the Havrda and Charvat entropy (HCE)
delta <- c(1.5, 2, 3)
z <- hce_ray(0.2, delta)
dat <- data.frame(x = delta , HCE = z)
p_hce <- ggplot(dat, aes(x = delta, y = HCE)) + geom_line()
plot <- p_hce + ggtitle(expression(alpha == 0.2))
}


