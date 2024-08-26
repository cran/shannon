\name{Truncated Student's t distribution}
\alias{rlse_st}
\alias{rlre_st}
\alias{rlhce_st}
\alias{rlae_st}
\title{Relative loss for various entropy measures using the truncated Student's t distribution}
\usage{
rlse_st(p, v)
rlre_st(p, v, delta)
rlhce_st(p, v, delta)
rlae_st(p, v, delta)
}
\description{
Compute the relative information loss of the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the truncated Student's t distribution.
}

\arguments{
\item{v}{
The strictly positive parameter of the Student distribution (\eqn{v > 0}), also called a degree of freedom.
}

\item{p}{
The truncation time \eqn{(p>0)}.
}
\item{delta}{
The strictly positive parameter  (\eqn{\delta > 0}) and (\eqn{\delta \ne 1}).
}
}

\value{
The functions rlse_st, rlre_st, rlhce_st, and rlae_st provide the relative information loss based on the Shannon entropy, Rényi entropy, Havrda and Charvat entropy, and Arimoto entropy, respectively, depending on the selected parametric values of the truncated Student's t distribution, \eqn{p} and \eqn{\delta}.
}
\references{
Yang, Z., Fang, K. T., & Kotz, S. (2007). On the Student's t-distribution and the t-statistic. Journal of Multivariate Analysis, 98(6), 1293-1304.

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
\code{\link{re_st}
}
}


\examples{
p <- c(1, 1.7, 3)
rlse_st(p, 4)
rlre_st(p, 4, 0.5)
rlhce_st(p, 4, 0.5)
rlae_st(p, 4, 0.5)
}


