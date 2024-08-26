\name{shannon-package}
\alias{shannon-package}
\docType{package}
\title{
Computation of Entropy Measures and Relative Loss
}

\description{

The functions allow for the numerical evaluation of some commonly used entropy measures, such as Shannon entropy, Rényi entropy, Havrda and Charvat entropy, and Arimoto entropy, at selected parametric values from several well-known and widely used probability distributions. Moreover, the functions also compute the relative loss of these entropies using the truncated distributions. Let \eqn{X} be an absolutely continuous random variable having the probability density function \eqn{f(x)}. Then, the Shahnon entropy is as follows:
\deqn{
H(X)=-\intop_{-\infty}^{+\infty}f(x)\log f(x)dx.
}
The  Rényi entropy is as follows:
\deqn{
H_{\delta}(X)=\frac{1}{1-\delta}\log\intop_{-\infty}^{+\infty}f(x)^{\delta}dx;\qquad\delta>0,\delta\ne1.
}
The Havrda and Charvat entropy is as follows:
\deqn{
H_{\delta}(X)=\frac{1}{2^{1-\delta}-1}\left(\intop_{-\infty}^{+\infty}f(x)^{\delta}dx-1\right);\qquad\delta>0,\delta\ne1.
}
The Arimoto entropy is as follows:
\deqn{
H_{\delta}(X)=\frac{\delta}{1-\delta}\left[\left(\intop_{-\infty}^{+\infty}f(x)^{\delta}dx\right)^{\frac{1}{\delta}}-1\right];\qquad\delta>0,\delta\ne1.
}

Let \eqn{D(X)} be an entropy, and \eqn{D_p(X)} be its truncated integral version at \eqn{p}, i.e., defined with the truncated version of \eqn{f(x)} over the interval \eqn{(-\infty,p)}. Then we define the corresponding relative loss entropy is defined by

\deqn{
S_D(p)= \frac{D(X)-D_p(X)}{D(X)}.
}

}


\details{
\tabular{ll}{
Package: \tab shannon\cr
Type: \tab Package\cr
Version: \tab 0.2.0 \cr
Date: \tab 2024-08-21\cr
License: \tab GPL-2\cr
}
}

\author{
Muhammad Imran \email{imranshakoor84@yahoo.com},
Christophe Chesneau
\email{christophe.chesneau@unicaen.fr}
and Farrukh Jamal \email{farrukh.jamal@iub.edu.pk}.
}


\section{Maintainers}{
Muhammad Imran <imranshakoor84@yahoo.com>
}

\keyword{Distribution theory}
\keyword{Entropy measures}
\keyword{Shannon entropy}
\keyword{Rényi entropy}
\keyword{Havrda and Charvat entropy}
\keyword{Arimoto entropy}
\keyword{Relative loss}
\keyword{Truncation time}
\keyword{Truncated model}

\references{
Shannon, C. E. (1948). A mathematical theory of communication. The Bell system technical journal, 27(3), 379-423.

Rényi, A. (1961). On measures of entropy and information,
Hungarian Academy of Sciences, Budapest, Hungary, 547-
561.

Havrda, J., & Charvat, F. (1967). Quantification method of classification processes. Concept of structural \eqn{\alpha}-entropy. Kybernetika, 3(1), 30-35.

Arimoto, S. (1971). Information-theoretical considerations on estimation problems. Information and control, 19(3), 181-194.

Awad, A. M., & Alawneh, A. J. (1987). Application of entropy to a life-time model. IMA Journal of Mathematical Control and Information, 4(2), 143-148.
}

