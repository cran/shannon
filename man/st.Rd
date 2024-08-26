\name{Student's t distribution}
\alias{se_st}
\alias{re_st}
\alias{hce_st}
\alias{ae_st}
\title{ Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the Student's t distribution}
\usage{
se_st(v)
re_st(v, delta)
hce_st(v, delta)
ae_st(v, delta)
}
\description{
Compute the Shannon, Rényi, Havrda and Charvat, and Arimoto entropies of the Student's t distribution.
}

\arguments{
\item{v}{
The strictly positive parameter of the Student's t distribution (\eqn{v > 0}), also called a degree of freedom.
}
\item{delta}{
The strictly positive parameter  (\eqn{\delta > 0}) and (\eqn{\delta \ne 1}).
}
}

\details{The following is the probability density function of the Student t distribution:
\deqn{
f(x)=\frac{\Gamma(\frac{v+1}{2})}{\sqrt{v\pi}\Gamma(\frac{v}{2})}\left(1+\frac{x^{2}}{v}\right)^{-(v+1)/2},
}
where \eqn{x\in\left(-\infty,+\infty\right)} and \eqn{v > 0}, and \eqn{\Gamma(a)} is the standard gamma function.
}
\value{
The functions se_st, re_st, hce_st, and ae_st provide the Shannon entropy, Rényi entropy, Havrda and Charvat entropy, and Arimoto entropy, respectively, depending on the selected parametric values of the Student's t distribution and \eqn{\delta}.
}
\references{
Yang, Z., Fang, K. T., & Kotz, S. (2007). On the Student's t-distribution and the t-statistic. Journal of Multivariate Analysis, 98(6), 1293-1304.

Ahsanullah, M., Kibria, B. G., & Shakil, M. (2014). Normal and Student's t distributions and their applications (Vol. 4). Paris, France: Atlantis Press.


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
\code{\link{re_exp}, \link{re_gamma}
}
}


\examples{
se_st(4)
delta <- c(1.5, 2, 3)
re_st(4, delta)
hce_st(4, delta)
ae_st(4, delta)
}


