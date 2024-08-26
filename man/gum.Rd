\name{Gumbel distribution}
\alias{Se_gum}
\alias{re_gum}
\alias{hce_gum}
\alias{ae_gum}
\title{Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the Gumbel distribution}
\usage{
Se_gum(alpha, beta)
re_gum(alpha, beta, delta)
hce_gum(alpha, beta, delta)
ae_gum(alpha, beta, delta)
}
\description{
Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the Gumbel distribution.
}

\arguments{
\item{alpha}{
The location parameter of the Gumbel distribution (\eqn{\alpha\in\left(-\infty,+\infty\right)}).
}
\item{beta}{
The strictly positive scale parameter of the Gumbel distribution (\eqn{\beta > 0}).
}
\item{delta}{
The strictly positive parameter  (\eqn{\delta > 0}) and (\eqn{\delta \ne 1}).
}
}
\details{The following is the probability density function of the Gumbel distribution:
\deqn{f(x)=\frac{1}{\beta}e^{-(z+e^{-z})},}
where \eqn{z=\frac{x-\alpha}{\beta}}, \eqn{x\in\left(-\infty,+\infty\right)}, \eqn{\alpha\in\left(-\infty,+\infty\right)} and \eqn{\beta > 0}.
}
\value{
The functions Se_gum, re_gum, hce_gum, and ae_gum provide the Shannon entropy, Rényi entropy, Havrda and Charvat entropy, and Arimoto entropy, respectively, depending on the selected parametric values of the Gumbel distribution and \eqn{\delta}.
}
\references{
Gomez, Y. M., Bolfarine, H., & Gomez, H. W. (2019). Gumbel distribution with heavy tails and applications to environmental data. Mathematics and Computers in Simulation, 157, 115-129.
}
\author{
Muhammad Imran, Christophe Chesneau and Farrukh Jamal

R implementation and documentation: Muhammad Imran <imranshakoor84@yahoo.com>, Christophe Chesneau <christophe.chesneau@unicaen.fr> and Farrukh Jamal \email{farrukh.jamal@iub.edu.pk}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{re_norm}
}
}


\examples{
Se_gum(1.2, 1.4)
delta <- c(2, 3)
re_gum(1.2, 0.4, delta)
hce_gum(1.2, 0.4, delta)
ae_gum(1.2, 0.4, delta)
}


