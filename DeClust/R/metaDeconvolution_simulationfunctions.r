

##################################################################################################
####
#' Simple version of DeClust with linear regression and  without iterative optimization 
#'
#' This function is a simple version of DeClust. It uses linear regression instead of nonlinear optimization assuming the gene expression follows normal distribution instead of log-normal distribution. In addition, unlike the function deClustFromMarker(), it doesn't use iteractive optimization.
#'
#' @param exprM a gene by sample expression matrix (in the original expression scale, not log-transformed,has to be non-negative)
#' @param k number of clusters
#' @param seed the random seed set for kmeans clustering
#' @return A named list with three components: subtype, subtypeprofileM and subtypefractionM. Subtype is a vector with length equal to the sample size, and it stores the sample clustering results;subtypeprofileM is a gene by compartment matrix, and it stores the inferred expression profile for each cancer subtype, immune and stromal compartment;subtypefractionM is a sample by compartment matrix, and it stores the estimated fraction of each compartment for each sample.
#' @examples  data(exprM);
#' r<-simpleFromMarker(exprM,3);
#' @export
simpleFromMarker<-function(exprM,k,seed=1)
  {
data(SI_geneset)
stromal<-as.matrix(SI_geneset)[1,-1];
immune<-as.matrix(SI_geneset)[2,-1];
if(sum(rownames(exprM)%in%immune)<10|sum(rownames(exprM)%in%stromal)<10)stop("too few markers for immune and stromal")
stromalMean<-exp(apply(log(exprM)[rownames(exprM)%in%stromal,],2,median))
immuneMean<-exp(apply(log(exprM)[rownames(exprM)%in%immune,],2,median))
markerM<-cbind(stromal=stromalMean/quantile(stromalMean,0.99),immune=immuneMean/quantile(immuneMean,0.99))


totalrM<-foreach( rates = seq(0.1,1,by=0.1))%dopar%
{
  print(rates)
  rM<-c()
  for(ratei in seq(0.1,1,by=0.1))
    {      
print(ratei)
fractionM<-cbind(stromal=markerM[,"stromal"]*rates,immune=markerM[,"immune"]*ratei)
    fractionM<-cbind(fractionM,cancer=1-fractionM[,"stromal"]-fractionM[,"immune"])
    if(mean(fractionM[,"cancer"]<0)>0.05)
    {
      ############we stop when there are more than 5% of samples with purity<0
      next();
    }else{
    fractionM[fractionM>1]<-1;
    fractionM[fractionM<0]<-0;
    r<-simpleFromCC(exprM,fractionM,k)
    reM<-r$subtypeprofileM%*%t(r$subtypefractionM);
    reM[reM<1]<-1           
    rM<-rbind(rM,c(rates=rates,ratei=ratei,MSE=mean(as.vector((log(reM)-log(exprM))^2))))
  }
}
rM;
}

rM<-do.call("rbind",totalrM)  
rM<-rM[which.min(rM[,"MSE"]),]
rates<-rM["rates"]
ratei<-rM["ratei"]
fractionM<-cbind(stromal=markerM[,"stromal"]*rates,immune=markerM[,"immune"]*ratei)
fractionM<-cbind(fractionM,cancer=1-fractionM[,"stromal"]-fractionM[,"immune"])
fractionM[fractionM>1]<-1;
fractionM[fractionM<0]<-0;

r<-simpleFromCC(exprM,fractionM,k,seed=seed)
r;
}


            


#######################################################################################
########################################################################################
####Use marker to infer fraction
#' Simple version of DeClust with linear regression 
#'
#' This function is a simple version of DeClust. It uses linear regression instead of nonlinear optimization assuming the gene expression follows normal distribution instead of log-normal distribution. Unlike simpleFromMarker, it uses iteractive optimization.
#'
#' @param exprM a gene by sample expression matrix (in the original expression scale, not log-transformed, has to be non-negative)
#' @param k number of clusters
#'  @param maxsubtypestep the maximum number of iterations in the inner lay of optimization for sample clustering
#' @param maxstep the maximum number of iterations in the outer lay of optimization for the fraction of each compartment
#' @param seed the random seed set for kmeans clustering
#' @return A named list with three components: subtype, subtypeprofileM and subtypefractionM. Subtype is a vector with length equal to the sample size, and it stores the sample clustering results;subtypeprofileM is a gene by compartment matrix, and it stores the inferred expression profile for each cancer subtype, immune and stromal compartment;subtypefractionM is a sample by compartment matrix, and it stores the estimated fraction of each compartment for each sample.
#' @examples  data(exprM);
#' r<-deClustFromMarker(exprM,3);
#' @export
deClustFromMarker<-function(exprM,k,maxsubtypestep=100,maxstep=100,seed=1)
  {
data(SI_geneset)
stromal<-as.matrix(SI_geneset)[1,-1];
immune<-as.matrix(SI_geneset)[2,-1];
if(sum(rownames(exprM)%in%immune)<10|sum(rownames(exprM)%in%stromal)<10)stop("too few markers for immune and stromal")
stromalMean<-exp(apply(log(exprM)[rownames(exprM)%in%stromal,],2,median))
immuneMean<-exp(apply(log(exprM)[rownames(exprM)%in%immune,],2,median))
markerM<-cbind(stromal=stromalMean/quantile(stromalMean,0.99),immune=immuneMean/quantile(immuneMean,0.99))


totalrM<-foreach( rates = seq(0.1,1,by=0.1))%dopar%
{
  print(rates)
  rM<-c()
  for(ratei in seq(0.1,1,by=0.1))
    {      
print(ratei)
fractionM<-cbind(stromal=markerM[,"stromal"]*rates,immune=markerM[,"immune"]*ratei)
    fractionM<-cbind(fractionM,cancer=1-fractionM[,"stromal"]-fractionM[,"immune"])
    if(mean(fractionM[,"cancer"]<0)>0.05)
    {
      ############we stop when there are more than 5% of samples with purity<0
      next();
    }else{
    fractionM[fractionM>1]<-1;
    fractionM[fractionM<0]<-0;
    r<-deClustFromCC(exprM,fractionM,k=k,maxsubtypestep=maxsubtypestep,seed=seed)
    reM<-r$subtypeprofileM%*%t(r$subtypefractionM);
    reM[reM<1]<-1           
    rM<-rbind(rM,c(rates=rates,ratei=ratei,MSE=mean(as.vector((log(reM)-log(exprM))^2))))
  }
}
rM;
}

rM<-do.call("rbind",totalrM)  
rM<-rM[which.min(rM[,"MSE"]),]
rates<-rM["rates"]
ratei<-rM["ratei"]
fractionM<-cbind(stromal=markerM[,"stromal"]*rates,immune=markerM[,"immune"]*ratei)
fractionM<-cbind(fractionM,cancer=1-fractionM[,"stromal"]-fractionM[,"immune"])
fractionM[fractionM>1]<-1;
fractionM[fractionM<0]<-0;

#r0<-simpleDeconv(exprM,fractionM,k=k,maxsubtypestep=maxsubtypestep)
r<-deClustFromCCiterative(exprM,fractionM,k=k,maxstep=maxstep,maxsubtypestep=maxsubtypestep,seed=seed)


r;
}





##############################################################################################################################from marker: exprM in the original scale, it has to be a matrix
###################

#' Full version of DeClust 
#'
#' This function is a full version of DeClust. It uses nonlinear optimization assuming the gene expression follows log-normal distribution. It also contains iteractive optimization procedures.
#'
#' @param exprM a gene by sample expression matrix (in the original expression scale, not log-transformed, has to be non-negative)
#' @param k number of clusters
#'  @param maxsubtypestep the maximum number of iterations in the inner lay of optimization for sample clustering
#' @param maxstep the maximum number of iterations in the outer lay of optimization for the fraction of each compartment
#' @param seed the random seed set for kmeans clustering
#' @return A named list with three components: subtype, subtypeprofileM and subtypefratioM. Subtype is a vector with length equal to the sample size, and it stores the sample clustering results;subtypeprofileM is a gene by compartment matrix, and it stores the inferred expression profile for each cancer subtype, immune and stromal compartment;subtypefractionM is a sample by compartment matrix, and it stores the estimated fraction of each compartment for each sample.
#' @examples  data(exprM);
#' library("doParallel");
#' cl<-makeCluster(5,type="FORK",outfile="");
#' registerDoParallel(cl);
#' r<-deClustFromMarkerlognormal(exprM,3);
#' @export
deClustFromMarkerlognormal<-function(exprM,k,maxsubtypestep=20,maxstep=20,seed=1)
  {
   options(warn=-1)
   data(SI_geneset)
stromal<-as.matrix(SI_geneset)[1,-1];
immune<-as.matrix(SI_geneset)[2,-1];
if(sum(rownames(exprM)%in%immune)<10|sum(rownames(exprM)%in%stromal)<10)stop("too few markers for immune and stromal")
stromalMean<-exp(apply(log(exprM)[rownames(exprM)%in%stromal,],2,median))
immuneMean<-exp(apply(log(exprM)[rownames(exprM)%in%immune,],2,median))
markerM<-cbind(stromal=stromalMean/quantile(stromalMean,0.99),immune=immuneMean/quantile(immuneMean,0.99))

#registerDoMC(10)

   rates<-c()
   ratei<-c()
totalrM<-foreach( rates = seq(0.1,1,by=0.1))%dopar%
{
 # source("~/projects/Li/TCGAtools/Rcode/metaDeconvolution_simulationfunctions.r")
  print(rates)
  rM<-c()
  for(ratei in seq(0.1,1,by=0.1))
    {      
print(ratei)
fractionM<-cbind(stromal=markerM[,"stromal"]*rates,immune=markerM[,"immune"]*ratei)
    fractionM<-cbind(fractionM,cancer=1-fractionM[,"stromal"]-fractionM[,"immune"])
    if(mean(fractionM[,"cancer"]<0)>0.05)
    {
      ############we stop when there are more than 5% of samples with purity<0
      next();
    }else{
    fractionM[fractionM>1]<-1;
    fractionM[fractionM<0]<-0;
    r<-deClustFromCC(exprM+1,fractionM,k=k,maxsubtypestep=maxsubtypestep,seed=seed)
    reM<-r$subtypeprofileM%*%t(r$subtypefractionM);
    reM[reM<1]<-1           
    rM<-rbind(rM,c(rates=rates,ratei=ratei,MSE=mean(as.vector((log(reM)-log(exprM+1))^2))))
  }
}
rM;
}

rM<-do.call("rbind",totalrM)  
rM<-rM[which.min(rM[,"MSE"]),]
rates<-rM["rates"]
ratei<-rM["ratei"]


fractionM<-cbind(stromal=markerM[,"stromal"]*rates,immune=markerM[,"immune"]*ratei)
fractionM<-cbind(fractionM,cancer=1-fractionM[,"stromal"]-fractionM[,"immune"])
fractionM[fractionM>1]<-1;
fractionM[fractionM<0]<-0;

# save(rates,ratei,fractionM,file=paste(outputDir,label,"_bestrate.rda",sep=""))

r<-deClustFromCCiterativelognormal(exprM,fractionM,k,maxstep=maxstep,maxsubtypestep=maxsubtypestep,seed=seed)
r;
}



#########################################################################################
 
#' calculate BIC from the output of DeClust 
#'
#'
#'
#' @param exprM the bulky expression data as input for deClust functions
#' @param deClustoutput the output of deClust functions
#'  
#' @return BIC value
#' @examples  data(exprM);
#' library("doParallel");
#' cl<-makeCluster(5,type="FORK",outfile="");
#' registerDoParallel(cl);
#' r<-deClustFromMarkerlognormal(exprM,3);
#' BIC<-calculateBIC(exprM,r)
#' @export

calculateBIC<-function(exprM,deClustoutput)
 {
    r<-deClustoutput
    mse<-mean(as.vector((log(r$subtypeprofileM%*%t(r$subtypefractionM)+1)-log(exprM+1))^2))
   sampleN<-length(r$subtype)
   subtypeN<-length(unique(r$subtype))
  return(log(mse)*sampleN+log(sampleN)*subtypeN);
}



##################################################################################################internal function
############################################################################################
#######################################################################################################optimizefractionM

optimizefractionM<-function(exprM,profileM,subtype,initfractionM)
  {
    ###########for samples with unknown initfractionM, assign to the subtype with the most similar profile of top 10% cancer genes, assign initfractionM randomly
    if(any(is.na(subtype)))
      {
        FC<-apply(profileM[,grep("subtype",colnames(profileM)),drop=FALSE],1,min)-apply(profileM[,c("stromal","immune")],1,max)
        cancergenes<-which(FC>quantile(FC,0.9))
      ms<-names(subtype)[is.na(subtype)]
      for(ss in ms)
        {
          subtypenames<-setdiff(unique(subtype),NA)
        subtype[ss]<-subtypenames[which.max(cor(exprM[cancergenes,ss],profileM[cancergenes,subtypenames]))[1]]
         initfractionM[ss,]<-0;
         initfractionM[ss,c("stromal","immune",subtype[ss])]<-1/3
       }
      }
    #####when geneV is too small(or even zero), it will cause problem, so remove extreme genes
    geneV<-apply((log((exp(profileM)-1)%*%t(initfractionM)+1)-exprM)^2,1,mean)
   # gg<-which(geneV<quantile(geneV,0.95)&geneV>quantile(geneV,0.05))
     
    r<-foreach(i = 1:ncol(exprM))%dopar%{
     # print(i)
      initfractionV<-initfractionM[i,c("stromal","immune",subtype[i])];
       ######when the initfractionM has close to zero fraction,the optimization may have issues
      initfractionV[initfractionV<0.0001]<-0.0001;
      t<-optimizefractionV(exprM[,i],profileM[,c("stromal","immune",subtype[i])],geneV,initfractionV)
      names(t)[3]<-"cancer"
      t;
              }
      #  proc.time()-pt
fractionM<-do.call("rbind",r)
fractionM<-as.matrix(fractionM[,1:3])
rownames(fractionM)<-rownames(initfractionM);
colnames(fractionM)<-c("stromal","immune","cancer");
fractionM<-getsubtypefractionM(fractionM,subtype)
fractionM;    
  }

##########################################################################################################
optimizefractionV<-function(exprV,profileM,geneV,initfractionV)
  {
fn1<-function(fractionV, profileM,exprV,geneV)
  {
    mean((log((exp(profileM)-1)%*%fractionV+1)-exprV)^2/geneV)
  }
gn1<-function(fractionV, profileM,exprV,geneV)
  {
    t(exp(profileM)-1)%*%((log((exp(profileM)-1)%*%fractionV+1)-exprV)*2/geneV/((exp(profileM)-1)%*%fractionV+1))/nrow(profileM)
  }

 optimx(par=initfractionV, fn=fn1,gr=gn1,lower=0, upper=Inf, itnmax=NULL, hessian=FALSE, method="L-BFGS-B",profileM=profileM,exprV=exprV,geneV=geneV)
}



########################################################################################################getsubtypefractionM
getsubtypefractionM<-function(fractionM,subtype)
  {

    subtypeM<-matrix(0,length(subtype),length(unique(subtype)),dimnames=list(names(subtype),unique(subtype)))
for(i in colnames(subtypeM))subtypeM[names(subtype)[subtype==i],i]<-1
if(any(sort(rownames(subtypeM))!=sort(rownames(fractionM))))stop("subtype samples don't match")
subtypeM<-subtypeM[rownames(fractionM),,drop=FALSE]
subtypeM<-sweep(subtypeM,1,fractionM[,"cancer"],"*")
fractionM<-cbind(fractionM[,c("stromal","immune")],subtypeM)
fractionM[fractionM<0]<-0;
fractionM
  }



#######################################################################################################optimizeprofileM
optimizeprofileM<-function(exprM,fractionM)
  {
    #registerDoMC()


fn1<-function(profilei,expri,fractionM)
  {
    mean((log(fractionM%*%(exp(profilei)-1)+1)-expri)^2)
  }

    
gr1<-function(profilei,expri,fractionM)
  {
    t<-2*(log(fractionM%*%(exp(profilei)-1)+1)-expri)/(fractionM%*%(exp(profilei)-1)+1)
    t(t)%*%fractionM/nrow(fractionM)*exp(profilei)
  }

print("start optimize profileM")
#pt<-proc.time()

r<-foreach(k=1:nrow(exprM))%dopar%{
 # print(k)
  options(warn=-1)
  p<-optimx(par=rep(mean(exprM[k,]),ncol(fractionM)), fn=fn1, gr=gr1, lower=0, upper=Inf, itnmax=NULL, hessian=FALSE, method="L-BFGS-B",expri=exprM[k,],fractionM=fractionM)
  ####sometimes optimx outputs NA for unclear reason, it might be fixed by changing initial value.
  if(any(is.na(p[1:ncol(fractionM)])))p<-optimx(par=rep(mean(exprM[k,])+1,ncol(fractionM)), fn=fn1, gr=gr1, lower=0, upper=Inf, itnmax=NULL, hessian=FALSE, method="L-BFGS-B",expri=exprM[k,],fractionM=fractionM)
  p;    
}#
     

                                        #  proc.time()-pt
r<-do.call("rbind",r)
curprofileM<-as.matrix(r[,1:ncol(fractionM)])
rownames(curprofileM)<-rownames(exprM);
colnames(curprofileM)<-colnames(fractionM);
curprofileM[curprofileM<0]<-0;
curprofileM;
}


########################################################################################
#############################################       
###calculate mse
calculateMSE<-function(fractionM,exprM,profileM)
{
  if(any(colnames(fractionM)!=colnames(profileM)))stop("fraction names need to match")
   mean(as.vector((log((exp(profileM)-1)%*%t(fractionM)+1)-exprM)^2))
}





################################################################################################from CC
#########################################################################################
#####################################################################################################interatively updating fractionM

deClustFromCCiterativelognormal<-function(exprM,fractionM,subtypeN,maxstep=maxstep,maxsubtypestep=maxsubtypestep,seed=1)
  {
   
   
    if(sum(is.na(exprM))>0)stop("NA is not allowed in exprM")
    if(any(colnames(exprM)!=rownames(fractionM)))stop("sample Name doesn't match")
   
    stepi<-1
  
    while(stepi<=maxstep)
      {
        print(paste("step",stepi))

       r<-deClustFromCClognormal(exprM=exprM,fractionM=fractionM,subtypeNum=subtypeN,maxsubtypestep=maxsubtypestep,seed=seed)
              
     newsubtypefractionM<-optimizefractionM(log(exprM+1),log(r$subtypeprofileM+1),r$subtype,r$subtypefractionM) 
     newfractionM<-cbind(newsubtypefractionM[,c("stromal","immune")],cancer=apply(newsubtypefractionM[,grep("subtype",colnames(newsubtypefractionM)),drop=FALSE],1,sum))
        print(max(abs(newfractionM-fractionM)))
     if(all(abs(newfractionM-fractionM)<0.01))break;
     fractionM<-newfractionM;
     stepi<-stepi+1;
     # save(r,stepi,file=paste(outputDir,label,"_",stepi,".rda",sep=""))
      }

return(r)
    
  }



##################################################################################################################################################################################
################## deconvolution from fractionM lognormal
####################fractionM has 3 columns: stromal,immune and cancer
      deClustFromCClognormal<-function(exprM,fractionM,subtypeNum,maxsubtypestep=100,seed=1)
{
  exprM<-log(exprM+1)
      
  #####initial ressubtype using logresidual without subtypes
profileM<-optimizeprofileM(exprM,fractionM)
if(subtypeNum==1){
  mse<-calculateMSE(fractionM,exprM,profileM)
  subtypefractionM<-fractionM;
  colnames(subtypefractionM)[colnames(subtypefractionM)=="cancer"]<-"subtype1"
  subtypeprofileM<-profileM;
  colnames(subtypeprofileM)[colnames(subtypeprofileM)=="cancer"]<-"subtype1"
  subtype<-rep("subtype1",nrow(subtypefractionM))
  names(subtype)<-rownames(subtypefractionM)
  return(list(mse=mse,subtypefractionM=subtypefractionM,subtypeprofileM=exp(subtypeprofileM)-1,subtype=subtype))
}

if(subtypeNum>1)
  {
logresM<-exprM-log((exp(profileM)-1)%*%t(fractionM)+1)
set.seed(seed)
r<-kmeans(t(logresM),subtypeNum,nstart=5)
ressubtype<-r$cluster;
ressubtype<-paste("subtype",ressubtype,sep="")
names(ressubtype)<-names(r$cluster)
ressubtype<-ressubtype[rownames(fractionM)]

##################
oldsubtype<-ressubtype
       step<-1
while(step<maxsubtypestep)
  {
    print(paste("subtypestep",step))
    subtypefractionM<-getsubtypefractionM(fractionM,oldsubtype)
    subtypeprofileM<-optimizeprofileM(exprM,subtypefractionM)
    mse<-calculateMSE(subtypefractionM,exprM,subtypeprofileM)
    print(mse)
####reassign cluster
   subtypeNames<-unique(oldsubtype)
   msesubtype<-sapply(subtypeNames,function(subtype)apply((exprM-log((exp(subtypeprofileM[,c("stromal","immune",subtype)])-1)%*%t(fractionM)+1))^2,2,mean))
   newsubtype<-apply(msesubtype,1,function(x) subtypeNames[which.min(x)])
   if(all(newsubtype==oldsubtype))break;
   print(sum(newsubtype!=oldsubtype))
   oldsubtype<-newsubtype
    step<-step+1;
  }   
#####
subtype<-oldsubtype
return(list(subtypefractionM=subtypefractionM,subtypeprofileM=exp(subtypeprofileM)-1,subtype=subtype))
}
}







########################################################################################################################################################################
#########################
simpleFromCC<-function(exprM,fractionM,k,seed=1)
  {
if(any(colnames(exprM)!=rownames(fractionM)))stop("sample doesn't match")
 model<-lm.fit(x=fractionM,y=t(exprM))
  reM<-t(model$coefficients)%*%t(fractionM)
  reM[reM<1]<-1
  adjexprM<-log(exprM)-log(reM)
set.seed(seed)
r<-kmeans(t(adjexprM),k,nstart=5)

ressubtype<-r$cluster;
ressubtype<-paste("subtype",ressubtype,sep="")
names(ressubtype)<-names(r$cluster)
ressubtype<-ressubtype[rownames(fractionM)]

subtypefractionM<-getsubtypefractionM(fractionM,ressubtype)
     model<-lm.fit(x=subtypefractionM,y=t(exprM))
     subtypeprofileM<-t(model$coefficients)
     subtypeprofileM[subtypeprofileM<1]<-1

return(list(subtype=ressubtype,subtypeprofileM=subtypeprofileM,subtypefractionM=subtypefractionM))
}



########################################################################################
#############################

deClustFromCC<-function(exprM,fractionM,k,maxsubtypestep=100,seed=1)
  {
 model<-lm.fit(x=fractionM,y=t(exprM))
  reM<-t(model$coefficients)%*%t(fractionM)
  reM[reM<1]<-1
  adjexprM<-log(exprM)-log(reM)
set.seed(seed)
r<-kmeans(t(adjexprM),k,nstart=5)

ressubtype<-r$cluster;
ressubtype<-paste("subtype",ressubtype,sep="")
names(ressubtype)<-names(r$cluster)
ressubtype<-ressubtype[rownames(fractionM)]

oldsubtype<-ressubtype
       step<-1
while(step<maxsubtypestep)
  {
    print(paste("subtypestep",step))
    subtypefractionM<-getsubtypefractionM(fractionM,oldsubtype)
     model<-lm.fit(x=subtypefractionM,y=t(exprM))
     subtypeprofileM<-t(model$coefficients)
     subtypeprofileM[subtypeprofileM<1]<-1;
####reassign cluster
   subtypeNames<-unique(oldsubtype)
   msesubtype<-sapply(subtypeNames,function(subtype)
                      {
                        reM<-subtypeprofileM[,c("stromal","immune",subtype)]%*%t(fractionM[,c("stromal","immune","cancer")])
                        reM[reM<1]<-1;
                        apply((log(exprM)-log(reM))^2,2,mean)
                      })
                    
   newsubtype<-apply(msesubtype,1,function(x) subtypeNames[which.min(x)])
   if(all(newsubtype==oldsubtype))break;
   print(sum(newsubtype!=oldsubtype))
   oldsubtype<-newsubtype
    step<-step+1;
  }

#####
 subtypeprofileM[subtypeprofileM<1]<-1;
 
list(subtype=newsubtype,subtypeprofileM=subtypeprofileM,subtypefractionM=subtypefractionM)
}


#####################################################################################
########################################################################################

deClustFromCCiterative<-function(exprM,fractionM,k=2,weight=TRUE,maxstep=100,maxsubtypestep=100,seed=1)
  {
    stepi<-1
    while(stepi<=maxstep)
      {
        print(paste("step",stepi))
        #if(stepi==20)browser()
        r<-deClustFromCC(exprM,fractionM,maxsubtypestep=maxsubtypestep,k=k,seed=seed)
        subtype<-r$subtype
  
        geneV<-apply((r$subtypeprofileM%*%t(r$subtypefractionM)-exprM)^2,1,mean)
        subtypeName<-unique(subtype)
    
     subtypefractionML<-lapply(subtypeName,function(st)
        {
      if(weight) model<-lm.wfit(x=r$subtypeprofileM[,c("stromal","immune",st)],y=exprM[,which(subtype==st)],w=1/geneV)
     if(!weight) model<-lm.fit(x=r$subtypeprofileM[,c("stromal","immune",st)],y=exprM[,which(subtype==st)])
     subtypefractionM<-t(model$coefficients)
      ###when there is only one sample there, the rownames is lost
     rownames(subtypefractionM)<-colnames(exprM)[which(subtype==st)] 
     subtypefractionM[subtypefractionM<0]<-0;
     subtypefractionM<-cbind(subtypefractionM,matrix(0,nrow(subtypefractionM),length(subtypeName)-1,dimnames=list(rownames(subtypefractionM),setdiff(subtypeName,st))))
     subtypefractionM[,colnames(r$subtypeprofileM),drop=FALSE];
    })

    subtypefractionM<-do.call("rbind",subtypefractionML)
    subtypefractionM<-subtypefractionM[colnames(exprM),]
     newfractionM<-cbind(subtypefractionM[,c("stromal","immune")],cancer=apply(subtypefractionM[,grep("subtype",colnames(subtypefractionM)),drop=FALSE],1,sum))

   
     if(all(abs(newfractionM-fractionM)<0.01))break;
      print(max(abs(newfractionM-fractionM)))
     fractionM<-newfractionM;
     stepi<-stepi+1;
      }

    r;
  }






