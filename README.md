# semipartail
Semiparametric tail estimation for heavy-tailed data. This package implements the method proposed by Fithian and Wager (2015).

To install this package in R, run the following commands:

```R
library(devtools) 
install_github("swager/semipartail")
```

Example usage:

```R
library(semipartail)
X = exp(rgamma(n=1000, shape=4, scale=0.45))
Z = exp(rgamma(n=1000000, shape=3, scale=0.45))
semipar.mean(main.sample=X, background.sample=Z, threshold=quantile(X, 0.8))
```

#### References
William Fithian and Stefan Wager.
<b>Semiparametric Exponential Families for Heavy-Tailed Data.</b>
<i>Biometrika</i>, 102(2):486–493, 2015.
[<a href="http://arxiv.org/abs/1307.7830">arxiv</a>]
