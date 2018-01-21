library(dscr)
library(flashr)
sessionInfo()

dsc_flashsim=new_dsc("dsc-flashsim","../output/dsc-flashsim-files")

sim_r1_paper =function(args){
  n=200;p=300;pi0=0.9;s=c(0.25,.5,1,2,4);tau=NULL;missing=FALSE

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
  return(list(input=list(Y=Y),meta=list(LF=LF)))
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
  res = ssvd::ssvd(input$Y)
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
  cv.out <- PMA::PMD.cv(input$Y, type="standard", sumabss=seq(0.1, 0.6, len=20),center=FALSE)
  res <-  PMA::PMD(input$Y, type="standard", sumabs=cv.out$bestsumabs, K=1, v=cv.out$v.init,center=FALSE)
  return(list(LFhat = res$d[1] * res$u[,1] %*% t(res$v[,1])))
}


rmse = function(a,b){return(sqrt(mean((a-b)^2)))}
rmse.wrapper = function(data,output){
  return(list(rmse=rmse(data$meta$LF,output$LFhat)))
}

add_scenario(dsc_flashsim,"strong-sparsity-r1",sim_r1_paper,seed=1:20)

add_method(dsc_flashsim,"flash_r1",flash_r1.wrapper)
add_method(dsc_flashsim,"flash_r1_pn",flash_r1_pn.wrapper)
add_method(dsc_flashsim,"ssvd",ssvd.wrapper)
add_method(dsc_flashsim,"svd",svd.wrapper)
add_method(dsc_flashsim,"pmd",pmd.wrapper)
add_method(dsc_flashsim,"pmd.cv",pmd.cv.wrapper)

add_score(dsc_flashsim,rmse.wrapper,"rmse")


res=run_dsc(dsc_flashsim)
save(res,dsc_flashsim,file="../output/dsc-flashsim-files/res.RData")

