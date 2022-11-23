evaluation  =  function(n, q, p, beta.gfl,eta.gfl,Beta0,
                       sample.index, a.matrix.main, a.matrix.GE)
{

  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: evaluation
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description:
  ##            Evaluating performances of proposed methods.
  ## -----------------------------------------------------------------------------------------------------------------


  pp = q+p
  main.n = c(1:q)
  GE.n = setdiff(1:pp,main.n)
  a.main.n = sort(rep((0:(n-1))*pp,q))+main.n
  a.GE.n = sort(rep((0:(n-1))*pp,p))+GE.n
  d.main.n =sort(rep((0:(n*(n-1)/2-1))*pp,q)) + main.n
  d.GE.n = sort(rep((0:(n*(n-1)/2-1))*pp,p)) + GE.n

  b.main = as.matrix(beta.gfl[a.main.n,])
  b.GE = as.matrix(beta.gfl[a.GE.n,])
  eta.main = as.matrix(eta.gfl[d.main.n,],sparse=T)
  eta.GE = as.matrix(eta.gfl[d.GE.n,],sparse=T)

  beta.true.list = Beta0$beta.true.list
  beta.true.main.list = Beta0$beta.true.main.list0
  beta.true.GE.list = Beta0$beta.true.GE.list0

  beta.main.esi = t(matrix(b.main,nrow = q))
  beta.main.true = t(do.call(cbind,beta.true.main.list))
  error.main  =  sqrt(sum((beta.main.esi - beta.main.true)^2)) / sqrt(n*q)
  beta.GE.esi = t(matrix(b.GE,nrow = p))
  beta.GE.true = t(do.call(cbind,beta.true.GE.list))
  error.GE  =  sqrt(sum((beta.GE.esi - beta.GE.true)^2)) / sqrt(n*p)

  RI.main = getErrorRate(n, b.main,eta.main, beta.true.main.list, q, a.matrix.main, sample.index)
  gr.main.est = RI.main$gr.num
  RI.GE = getErrorRate(n, b.GE,eta.GE, beta.true.GE.list, p, a.matrix.GE, sample.index)
  gr.GE.est = RI.GE$gr.num

  return(list(RI.main=RI.main, RI.GE=RI.GE, error.main=error.main, error.GE=error.GE))
}


getErrorRate  =  function(size, beta.gfl, eta.gfl, beta.true.list.aa,
                         q, a.matrix.aa, sample.index)
{

  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: getErrorRate
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description:
  ##            Calculating the indicators of grouping accuracy.
  ## -----------------------------------------------------------------------------------------------------------------

  # true group matrix
  beta.true.mat = do.call(rbind,beta.true.list.aa)
  diffbeta = a.matrix.aa%*%as.matrix(beta.true.mat)
  dim(diffbeta) = c(q,size*(size-1)/2)
  epsi.betatrue = apply(diffbeta,2,sum)

  captrue.matrix = matrix(0,nrow=size,ncol=size)
  for(i in 1:length(epsi.betatrue)){
    if(epsi.betatrue[i]==0){
      captrue.matrix[sample.index[1,i],sample.index[2,i]] = 1
    }
  }

  diff.gfl = apply(matrix(eta.gfl,nrow=q,ncol=size*(size-1)/2),2,sum)

  capgfl.matrix = matrix(0,nrow=size,ncol=size)
  if(length(which(diff.gfl==0))==0){
    captrue.matrix2 = captrue.matrix+t(captrue.matrix); diag(captrue.matrix2)=1
    capgfl.matrix2 = capgfl.matrix+t(capgfl.matrix); diag(capgfl.matrix2)=1
    group_error_rate = sum(abs(captrue.matrix-capgfl.matrix))/(size*(size-1)/2)
    group_tp = length(intersect(which(as.vector(captrue.matrix)==1),which(as.vector(capgfl.matrix)==1)))/sum(captrue.matrix)
    group_fp = length(intersect(which(as.vector(captrue.matrix)==0),which(as.vector(capgfl.matrix)==1)))/(size*(size-1)/2-sum(captrue.matrix))
    return(list(group.error=c(group_error_rate,group_tp,group_fp), gr.num=size, matrix=list(captrue.matrix2=captrue.matrix2,capgfl.matrix2=capgfl.matrix2)))
  }
  if(length(which(diff.gfl==0))==size*(size-1)/2){
    capgfl.matrix[upper.tri(capgfl.matrix)]=1
    captrue.matrix2 = captrue.matrix+t(captrue.matrix); diag(captrue.matrix2)=1
    capgfl.matrix2 = capgfl.matrix+t(capgfl.matrix); diag(capgfl.matrix2)=1
    group_error_rate = sum(abs(captrue.matrix-capgfl.matrix))/(size*(size-1)/2)
    group_tp = length(intersect(which(as.vector(captrue.matrix)==1),which(as.vector(capgfl.matrix)==1)))/sum(captrue.matrix)
    group_fp = length(intersect(which(as.vector(captrue.matrix)==0),which(as.vector(capgfl.matrix)==1)))/(size*(size-1)/2-sum(captrue.matrix))
    return(list(group.error=c(group_error_rate,group_tp,group_fp), gr.num=1, matrix=list(captrue.matrix2=captrue.matrix2,capgfl.matrix2=capgfl.matrix2)))
  }else{
    sample.index.gfl = sample.index[,which(diff.gfl==0)]
    if(length(which(diff.gfl==0))==1){
      capgfl.matrix[sample.index.gfl[1],sample.index.gfl[2]] = 1
    }else{
      for(i in 1:length(which(diff.gfl==0))){
        capgfl.matrix[sample.index.gfl[1,i],sample.index.gfl[2,i]] = 1
      }
    }

    captrue.matrix2 = captrue.matrix+t(captrue.matrix); diag(captrue.matrix2)=1
    capgfl.matrix2 = capgfl.matrix+t(capgfl.matrix); diag(capgfl.matrix2)=1
    group.num.gf = nrow(unique(capgfl.matrix2))

    group_error_rate = sum(abs(captrue.matrix-capgfl.matrix))/(size*(size-1)/2)
    group_tp = length(intersect(which(as.vector(captrue.matrix)==1),which(as.vector(capgfl.matrix)==1)))/sum(captrue.matrix)
    group_fp = length(intersect(which(as.vector(captrue.matrix)==0),which(as.vector(capgfl.matrix)==1)))/(size*(size-1)/2-sum(captrue.matrix))
    group_fn = length(intersect(which(as.vector(captrue.matrix)==1),which(as.vector(capgfl.matrix)==0)))/sum(captrue.matrix)
    return(list(group.error=c(group_error_rate,group_tp,group_fp), gr.num=group.num.gf, matrix=list(captrue.matrix2=captrue.matrix2,capgfl.matrix2=capgfl.matrix2)))
  }
}




