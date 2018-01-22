library(dscr)
library(flashr)
sessionInfo()

dsc_flashsim=new_dsc("dsc-flashsim","../output/dsc-flashsim-files")

sim_r1_paper =function(n=200,p=300,pi0=0.9,s=sqrt(c(0.25,.5,1,2,4)),tau=NULL,missing=FALSE){

  #set tau as in paper for defaults
  if(is.null(tau)){
    if(pi0==0.9){tau=1}
    else if(pi0==0.3){tau=1/16}
    else {stop("must specify tau")}
  }

  f = rnorm(p)

  z = sample(1:5,n,replace=TRUE)
  l = rnorm(n,mean=0,sd=s[z])

  zero = rbinom(n,1,pi0)
  l[zero==1] = 0

  LF = l %*% t(f)

  Y = LF + sqrt(1/tau) * matrix(rnorm(n*p),nrow=n,ncol=p)
  if(missing){
    miss = matrix(rbinom(n*p,1,0.5),nrow=n)
    Y[miss==1] = NA
  }

  return(list(Y=Y,LF=LF))
}

sim_r1_paper.wrapper = function(args){
  res = do.call(sim_r1_paper,args = args)
  return(list(meta = list(LF=res$LF),input = list(Y=res$Y)))
}


flash_r1.wrapper = function(input,args){
  f=flashr::flash_r1(input$Y)
  return(list(LFhat = flash_get_lf(f)))
}

flash_r1_pn.wrapper = function(input,args){
  f=flashr::flash_r1(input$Y,ebnm_fn=ebnm_pn)
  return(list(LFhat = flash_get_lf(f)))
}

ssvd.wrapper = function(input,args){
  res = ssvd::ssvd(input$Y,method = "theory")
  return(list(LFhat = res$d * res$u %*% t(res$v)))
}

ssvd.method.wrapper = function(input,args){
  res = ssvd::ssvd(input$Y,method = "method")
  return(list(LFhat = res$d * res$u %*% t(res$v)))
}

svd.wrapper = function(input,args){
  res = svd(input$Y)
  return(list(LFhat = res$d[1] * res$u[,1] %*% t(res$v[,1])))
}

pmd.wrapper = function(input,args){
  res = PMA::PMD(input$Y,center=FALSE)
  return(list(LFhat = res$d[1] * res$u[,1] %*% t(res$v[,1])))
}

pmd.cv.wrapper = function(input,args){
  cv.out <- PMA::PMD.cv(input$Y, type="standard", sumabss=seq(0.1, 1, len=20),center=FALSE)
  res <-  PMA::PMD(input$Y, type="standard", sumabs=cv.out$bestsumabs, K=1, v=cv.out$v.init,center=FALSE)
  return(list(LFhat = res$d[1] * res$u[,1] %*% t(res$v[,1])))
}


# this a helper function by Wei to do CV
CVPMD_softImpute=function(Y,c_s,K,fold = 10, method = "PMD"){
  N = dim(Y)[1]
  P = dim(Y)[2]
  colindex = matrix(sample(P,P),ncol = fold)
  rowindex = matrix(sample(N,N),ncol = fold)

  missing= array(0,dim = c(fold,N,P))
  foldindex = array(0,dim = c(fold,fold,2))
  for(i in 1:fold){
    for(j in 1:fold){
      foldindex[i,j,1] = i
      foldindex[i,j,2] = (i+j) %% fold
    }
  }
  foldindex[which(foldindex == 0)] = fold
  for(i in 1:fold){
    missing[i, , ] = Y
    for(j in 1:fold){
      missing[i,rowindex[,foldindex[j,i,1]],colindex[,foldindex[j,i,2]]] = NA
    }
    missing[i,,which(colSums(missing[i,,],na.rm = T) ==0)] = Y[,which(colSums(missing[i,,],na.rm = T) ==0)]
  }
  # c_s is the candicate of shrinkage parameter
  n_s = length(c_s)
  # rmse for each grids
  CVRMSE = rep(0,n_s)
  minrmse = Inf
  opt_s = 0
  # for each candidate, we run it N_sim times
  for(t_s in 1:n_s){
    # for each grid
    # each time we set the rmse to zeros
    rmse = rep(0,fold)
    for(i in 1:fold){
      if(method == "PMD"){
        res_log = capture.output({out = PMD(missing[i,,], sumabs = c_s[t_s], sumabsv = NULL, sumabsu = NULL,K = K)})
      }else{
        out = softImpute::softImpute(missing[i,,], rank.max = K,lambda = c_s[t_s])
      }
      if(length(out$d)==1){
        misshat = (out$d) * out$u %*% t(out$v)
      }else{
        misshat = out$u %*%  diag(out$d) %*% t(out$v)
      }
      for(j in 1:fold){
        # for each fold j
        rmse[i] = rmse[i] + sum((Y[rowindex[,foldindex[j,i,1]],colindex[,foldindex[j,i,2]]] -
                                   misshat[rowindex[,foldindex[j,i,1]],colindex[,foldindex[j,i,2]]])^2,na.rm = TRUE)
      }
    } #get the result for one run
    CVRMSE[t_s] = CVRMSE[t_s] + sqrt(sum(rmse)/(N*P))
    if(CVRMSE[t_s] < minrmse){
      minrmse = CVRMSE[t_s]
      opt_s = c_s[t_s]
    }
  }
  return(list(opt_s = opt_s, output = CVRMSE))
}


# PMA.wrapper = function(Y,ngrids = 10,K,fold = 10){
#   library(PMA)
#   N = dim(Y)[1]
#   P = dim(Y)[2]
#   c_s = seq(0.3,0.8,len=ngrids)
#   cvout = CVPMD_softImpute(Y,c_s,K ,fold , method = "PMD")
#   res_log = capture.output({out = PMD(Y,sumabsu = NULL, sumabsv = NULL, sumabs = cvout$opt_s ,K = K)})
#   return(list(d = out$d, u = out$u, v = out$v))
# }

softImpute.cv.wrapper = function(input,args){
  Y= input$Y
  K=1
  ngrids = args$ngrids
  fold = args$fold
  max = args$max

  N = dim(Y)[1]
  P = dim(Y)[2]
  c_s = seq(0,max,len=ngrids)
  cvout = CVPMD_softImpute(Y,c_s,K ,fold , method = "softImpute")
  res = softImpute::softImpute(Y, rank.max = K,lambda = cvout$opt_s)
  return(list(LFhat = res$d[1] * res$u %*% t(res$v)))
}



### Scores


rmse = function(a,b){return(sqrt(mean((a-b)^2)))}
rrmse.wrapper = function(data,output){
  return(list(rmse=rmse(data$meta$LF,output$LFhat), rrmse=rmse(data$meta$LF,output$LFhat)/sqrt(mean(data$meta$LF^2))))
}

add_scenario(dsc_flashsim,"strong-sparsity-r1",sim_r1_paper.wrapper,args=list(),seed=1:20)
add_scenario(dsc_flashsim,"intermediate-sparsity-r1",sim_r1_paper.wrapper,args=list(pi0=0.3),seed=1:20)
add_scenario(dsc_flashsim,"no-sparsity-r1",sim_r1_paper.wrapper,args=list(pi0=0,tau=1/16),seed=1:20)


add_method(dsc_flashsim,"flash_r1",flash_r1.wrapper)
add_method(dsc_flashsim,"flash_r1_pn",flash_r1_pn.wrapper)
add_method(dsc_flashsim,"ssvd",ssvd.wrapper)
add_method(dsc_flashsim,"ssvd.method",ssvd.method.wrapper)
add_method(dsc_flashsim,"svd",svd.wrapper)
add_method(dsc_flashsim,"pmd",pmd.wrapper)
add_method(dsc_flashsim,"pmd.cv",pmd.cv.wrapper)
add_method(dsc_flashsim,"softimpute.cv",softImpute.cv.wrapper,args=list(ngrids=10,fold=5,max=100))
add_method(dsc_flashsim,"softimpute2.cv",softImpute.cv.wrapper,args=list(ngrids=10,fold=5,max=10))


add_score(dsc_flashsim,rrmse.wrapper,"rmse")


res=run_dsc(dsc_flashsim)
save(res,dsc_flashsim,file="../output/dsc-flashsim-files/res.RData")


