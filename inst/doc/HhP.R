## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  library(HhP)
#  library(Matrix)
#  library(MASS)
#  library(fmrs)
#  data(example.data.reg)
#  n   = example.data.reg$n
#  q   = example.data.reg$q
#  p   = example.data.reg$p
#  # ------------ Necessary parameters to support algorithm implementation --------
#  beta.init.list  =  gen_int_beta(n, p, q, example.data.reg)
#  beta.init  =  beta.init.list$beta.init
#  lambda  =  genelambda.obo()
#  result  =  HhP.reg(lambda, example.data.reg, n, q, p, beta.init)
#  index.list  =  evaluation.sum(n,q,p, result$admmres, result$abic.n, result$admmres2, example.data.reg$Beta0, result$bic.var)
#  index.list$err.s
#  

