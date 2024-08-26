\name{Truncated Rayleigh  distribution}
\alias{rlse_ray}
\alias{rlre_ray}
\alias{rlhce_ray}
\alias{rlae_ray}
\title{Relative loss for various entropy measures using the truncated Rayleigh  distribution}
\usage{
rlse_ray(p, alpha)
rlre_ray(p, alpha, delta)
rlhce_ray(p, alpha, delta)
rlae_ray(p, alpha, delta)
}
\description{
Compute the relative information loss of the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the truncated Rayleigh  distribution.
}

\arguments{
\item{alpha}{
The strictly positive scale parameter of the Rayleigh distribution (\eqn{\alpha > 0}).
}

\item{p}{
The truncation time \eqn{(p>0)}.
}
\item{delta}{
The strictly positive parameter  (\eqn{\delta > 0}) and (\eqn{\delta \ne 1}).
}
}

\value{
The functions rlse_ray, rlre_ray, rlhce_ray, and rlae_ray provide the relative information loss based on the Shannon entropy, Rényi entropy, Havrda and Charvat entropy, and Arimoto entropy, respectively, depending on the selected parametric values of the truncated Rayleigh  distribution, \eqn{p} and \eqn{\delta}.
}
\references{
Dey, S., Maiti, S. S., & Ahmad, M. (2016). Comparison of different entropy measures. Pak. J. Statist, 32(2), 97-108.

Awad, A. M., & Alawneh, A. J. (1987). Application of entropy to a life-time model. IMA Journal of Mathematical Control and Information, 4(2), 143-148.
}

\author{
Muhammad Imran, Christophe Chesneau and Farrukh Jamal

R implementation and documentation: Muhammad Imran <imranshakoor84@yahoo.com>, Christophe Chesneau <christophe.chesneau@unicaen.fr> and Farrukh Jamal \email{farrukh.jamal@iub.edu.pk}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{re_ray}
}
}


\examples{
p <- seq(0.25, 2, by=0.25)
rlse_ray(p, 2)
rlre_ray(p, 2, 0.5)
rlhce_ray(p, 2, 0.5)
rlae_ray(p, 2, 0.5)

# A graphic representation of relative loss (RL)
library(ggplot2)
# p is a truncation time vector
p <- seq(0.25, 2, by = 0.25)
# RL based on the Rényi entropy
z1 <- rlre_ray(p, 0.1, 0.5)
# RL based on the Havrda and Charvat entropy
z2 <- rlhce_ray(p, 0.1, 0.5)
# RL based on the Arimoto entropy
z3 <- rlae_ray(p, 0.1, 0.5)
# RL based on the Shannon entropy
z4 <- rlse_ray(p, 0.1)
df <- data.frame(x = p, RL = z1, z2, z3, z4)
head(df)
p1 <- ggplot(df, aes(x = p, y = RL, color = Entropy))
p1 + geom_line(aes(colour = "RE"), size = 1) + geom_line(aes(x,
    y = z2, colour = "HCE"), size = 1) + geom_line(aes(x, y = z3,
    colour = "AR"), size = 1) + geom_line(aes(x, y = z4, colour = "SE"),
    size = 1) + ggtitle(expression(delta == 0.5 ~ ~alpha == 0.1))
}


