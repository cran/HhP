ADMM = function(n, q, p, data.total, beta.init, lambda1, lambda2, p_y=1,
                iter.max=50, epsi=0.1, gam.mcp=3, penal.para=1, merge.all=F,
                sample.index, a.matrix, one.matrix, one.matrix.main, one.matrix.GE)
{

  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Description:
  ##            The key function of hierarchical heterogeneity analysis:
  ##            The implementation of ADMM.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ n: The sample size.
  ## @ q: The dimension of type 1 features.
  ## @ p: The dimension of type 2 features.
  ## @ data.total: The input data analyzed (a list including the response and design matrix).
  ## @ beta.init: The Initial values of regression coefficients.
  ## @ lambda1: The tuning parameter controlling the number of refined subgroup.
  ## @ lambda2: the tuning parameter controlling the number of rough subgroup.
  ## @ p_y: The dimension of y.
  ## @ iter.max: int, Maximum number of cycles of the ADMM algorithm, the default setting is 5.
  ## @ epsi: a float value, algorithm termination threshold.
  ## @ gam.mcp: The regularization parameter in MCP, the default setting is 3.
  ## @ penal.para: The penalty parameter in ADMM algorithm, the default setting is 1.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Output:
  ## A list including estimated regression coefficients, subgroups, BIC value, and so on
  ## ------------------------------------------------------------------------------------------------------------------------------------------

  pp=q+p
  main.n=c(1:q)
  GE.n=setdiff(1:pp,main.n)
  a.main.n = sort(rep((0:(n-1))*pp,q))+main.n
  a.GE.n = sort(rep((0:(n-1))*pp,p))+GE.n
  d.main.n =sort(rep((0:(n*(n-1)/2-1))*pp,q)) + main.n
  d.GE.n = sort(rep((0:(n*(n-1)/2-1))*pp,p)) + GE.n

  cova = data.total$data_x
  data.response = data.total$data_y

  row_vec = vector(mode="list",length=n)
  for(i5 in 1:n){
    row_vec[[i5]] = matrix(cova[i5,],1,pp)
  }
  cova_diag = bdiag(row_vec)

  #################################### Initialize beta #########################################
  # group
  eta.init = a.matrix%*%beta.init
  gamma.init = Matrix(0,nrow=p_y,ncol=pp*n*(n-1)/2,sparse = T)

  # iteration
  iter = 1
  # coe
  beta.est.list = vector(mode="list",length=iter.max);beta.est.list[[iter]] = beta.init
  # differences of coe
  eta.est.list = vector(mode="list",length=iter.max);eta.est.list[[iter]] = eta.init
  # the dual variables
  gamma.est.list = vector(mode="list",length=iter.max);gamma.est.list[[iter]] = gamma.init

  while(iter<=iter.max){
    iter = iter+1
    if(iter>iter.max){
      break
    }

    beta.est.list[[iter]] = solve(t(cova_diag)%*%cova_diag+penal.para*t(a.matrix)%*%a.matrix)%*%
      (t(cova_diag)%*%data.response-t(a.matrix)%*%t(gamma.est.list[[iter-1]])+penal.para*(t(a.matrix)%*%eta.est.list[[iter-1]]))

    eta.i = a.matrix%*%beta.est.list[[iter]]
    gamma.i = gamma.est.list[[iter-1]]
    eta.tem1 = eta.i+t(gamma.i)/penal.para
    eta.tem1.main = as.matrix(eta.tem1[d.main.n,],sparse=T)
    eta.tem1.GE = as.matrix(eta.tem1[d.GE.n,],sparse=T)
    eta.i1 = eta.i

    # 2-norm
    eta.tem2norm.main = sqrt(t(one.matrix.main)%*%(eta.tem1.main^2))
    eta.tem2norm.GE = sqrt(t(one.matrix.GE)%*%(eta.tem1.GE^2))
    eta.tem2norm = sqrt(eta.tem2norm.GE^2 + eta.tem2norm.main^2)

    # all:non-shrinkage & main:non-shrinkage
    num_n_n = which(eta.tem2norm > gam.mcp*lambda1 & eta.tem2norm.main > gam.mcp*lambda2)

    # all:shrinkage & main:non-shrinkage
    vu.coe = apply(1-lambda1/penal.para/eta.tem2norm,1,function(x) max(x,0)) / (1-1/(gam.mcp*penal.para))
    vv =  vu.coe * eta.tem2norm.main
    num_s_n = which(eta.tem2norm <= gam.mcp*lambda1 & vv > gam.mcp*lambda2)

    # all:non-shrinkage & main:shrinkage
    vv.coe = apply(1-lambda2/penal.para/eta.tem2norm.main,1,function(x) max(x,0)) / (1-1/(gam.mcp*penal.para))
    ww = sqrt(eta.tem2norm.GE^2 + ( vv.coe * eta.tem2norm.main)^2)
    num_n_s = which(ww > gam.mcp*lambda1 & eta.tem2norm.main <= gam.mcp*lambda2)

    # all:shrinkage & main:shrinkage
    num_s_s = setdiff(1:(n*(n-1)/2),Reduce(union,list(num_n_n,num_s_n,num_n_s)))

    if(length(num_n_n)<2){
      num_n_n2 = which(rowSums(as.matrix(one.matrix[,num_n_n]))!=0)
    }else{num_n_n2 = which(rowSums(one.matrix[,num_n_n])!=0)}
    if(length(num_s_n)<2){
      num_s_n2 = which(rowSums(as.matrix(one.matrix[,num_s_n]))!=0)
    }else{num_s_n2 = which(rowSums(one.matrix[,num_s_n])!=0)}
    if(length(num_n_s)<2){
      num_n_s2 = which(rowSums(as.matrix(one.matrix[,num_n_s]))!=0)
    }else{num_n_s2 = which(rowSums(one.matrix[,num_n_s])!=0)}
    if(length(num_s_s)<2){
      num_s_s2 = which(rowSums(as.matrix(one.matrix[,num_s_s]))!=0)
    }else{num_s_s2 = which(rowSums(one.matrix[,num_s_s])!=0)}

    # all:non-shrinkage & main:non-shrinkage
    if(length(num_n_n2) > 0){eta.i1[num_n_n2,] = eta.i[num_n_n2,]}
    length(num_n_n2)
    # all:shrinkage & main:non-shrinkage
    num = num_s_n; num2 = num_s_n2;
    if(length(num) > 0){
      eta.tem3 = as.matrix(vu.coe[num])
      eta.tem4 = as.vector(apply(eta.tem3, 1, function(x) rep(x,pp)))
      eta.i1[num2,] = eta.tem4*eta.i[num2,]
    }

    # all:non-shrinkage & main:shrinkage
    num = num_n_s; num2 = num_n_s2;
    if(length(num) > 0){
      num.main = num2[sort(rep((1:length(num)-1)*pp,q)) + main.n]
      num.GE = num2[sort(rep((1:length(num)-1)*pp,p)) + GE.n]
      eta.tem3 = as.matrix(vv.coe[num])
      eta.tem4 = as.vector(apply(eta.tem3, 1, function(x) rep(x,q)))
      eta.i1[num.main,] = eta.tem4*eta.i[num.main,]
      eta.i1[num.GE,] = eta.i[num.GE,]
    }

    # all:shrinkage & main:shrinkage
    num = num_s_s; num2 = num_s_s2;
    if(length(num) > 0){
      num.main = num2[sort(rep((1:length(num)-1)*pp,q)) + main.n]
      num.GE = num2[sort(rep((1:length(num)-1)*pp,p)) + GE.n]
      eta.i0norm = sqrt(t(one.matrix)%*%(eta.est.list[[iter-1]]^2))
      eta.i0norm.main = sqrt(t(one.matrix.main)%*%(eta.est.list[[iter-1]][d.main.n,]^2))
      mcp.u = mcp_d(eta.i0norm[num],lambda1,gam.mcp)/penal.para/(eta.i0norm[num]+10^(-7))
      mcp.v = mcp_d(eta.i0norm.main[num],lambda2,gam.mcp)/penal.para/(eta.i0norm.main[num]+10^(-7))
      eta.tem3.main = as.matrix(mcp.u+mcp.v)
      eta.tem3.GE = as.matrix(mcp.u)
      eta.tem4.main = as.vector(apply(eta.tem3.main, 1, function(x) rep(x,q)))
      eta.tem4.GE = as.vector(apply(eta.tem3.GE, 1, function(x) rep(x,p)))
      eta.i1[num.main,] = eta.i[num.main,]/(1+eta.tem4.main)
      eta.i1[num.GE,] = eta.i[num.GE,]/(1+eta.tem4.GE)
      eta.i1[abs(eta.i1) < 10^(-9)]=0
    }

    eta.est.list[[iter]] = eta.i1
    gamma.est.list[[iter]] = gamma.est.list[[iter-1]]+penal.para*t(a.matrix%*%beta.est.list[[iter]]-eta.est.list[[iter]])
    eps.group = sqrt(sum((a.matrix%*%beta.est.list[[iter]]-eta.est.list[[iter]])^2))

    if(eps.group<epsi){
      break
    }

    beta.est.list[iter-1] = list(NULL)
    eta.est.list[iter-1] = list(NULL)
    gamma.est.list[iter-1] = list(NULL)
  }

  if(iter>iter.max){
    beta.gfl = beta.est.list[[iter-1]]
    eta.gfl = eta.est.list[[iter-1]]
  }else{
    beta.gfl = beta.est.list[[iter]]
    eta.gfl = eta.est.list[[iter]]
  }

  b.main = as.matrix(beta.gfl[a.main.n,])
  b.GE = as.matrix(beta.gfl[a.GE.n,])
  eta.main = as.matrix(eta.gfl[d.main.n,],sparse=T)
  eta.GE = as.matrix(eta.gfl[d.GE.n,],sparse=T)

  rm(beta.est.list)
  rm(eta.est.list)
  rm(gamma.est.list)
  gc()
  ###################
  num.main  =  get.gr.num(n, b.main,eta.main, q, sample.index, merge.all=merge.all)
  num.GE  =  get.gr.num(n, b.GE,eta.GE, p, sample.index, merge.all=merge.all)
  gr.main.est = num.main$gr.num
  gr.GE.est = num.GE$gr.num
  residual  =  log(sum((data.response-cova_diag%*%beta.gfl)^2/n))
  BIC.var = residual+log(n*pp)*log(n)*(gr.main.est*(q)+gr.GE.est*(p))/n
  BIC.o = residual+log(n)*(gr.main.est*q+gr.GE.est*p)/n

  return(list(beta.gfl=beta.gfl, b.main=b.main, b.GE=b.GE, eta.gfl=eta.gfl, eta.main=eta.main, eta.GE=eta.GE,
              residual=residual,BIC.var=BIC.var,BIC.o=BIC.o,num.main=num.main,num.GE=num.GE))
}

InitialBeta  =  function(data.total, q, size, lambda.min=0.001, gr.init=10, a.matrix){

  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: InitialBeta
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description:
  ##            Generating the initial values using the ridge fusion method.
  ## -----------------------------------------------------------------------------------------------------------------

  cova = data.total$data_x
  data.response = data.total$data_y

  # equal beta
  row_vec = vector(mode="list",length=size)
  for(i5 in 1:size){
    row_vec[[i5]] = matrix(cova[i5,],1,q)
  }
  cova_diag = bdiag(row_vec)
  beta.ridge = solve(t(cova_diag)%*%cova_diag+lambda.min*(t(a.matrix)%*%a.matrix))%*%t(cova_diag)%*%data.response

  beta.ridge.mat = matrix(beta.ridge,nrow=q,ncol=size)
  beta.ridge.med = apply(beta.ridge.mat,2,stats::median)

  group.size = size/gr.init

  # lasso
  beta.order = order(beta.ridge.med)
  cova_group = vector(mode="list",length=gr.init)
  respon_group = vector(mode="list",length=gr.init)
  beta.init.list = vector(mode="list",length=size)
  for(i in 1:gr.init){
    ind = beta.order[c(((i-1)*group.size+1):(i*group.size))]
    cova_group[[i]] = cova[ind,]
    respon_group[[i]] = data.response[ind,]
    beta.ols = unname(solve(t(cova_group[[i]])%*%cova_group[[i]])%*%t(cova_group[[i]])%*%respon_group[[i]])
    for(j in 1:length(ind)){
      beta.init.list[[ind[j]]] = beta.ols
    }
  }
  beta.init  =  NULL
  for (i in 1:size) {
    bi  =  beta.init.list[[i]]
    if(length(bi) > 0){
      beta.init  =  c(beta.init,bi)
    } else {
      bi  =  rep(0,q)
      beta.init  =  c(beta.init,bi)
    }
  }
  beta.init  =  as.matrix(beta.init)
  return(list(beta.init=beta.init))
}
