#!/usr/bin/R
#liusiyang 20210315
#---
#perform cox regression on biomarkers for stroke reoccurrence or death;
#perform logistic regression on biomarkers for MRS status
#---
#install.packages(c("survival", "survminer","metafor"))
library("survival")
library("survminer")
library("metafor")
library("dplyr")
library("magrittr")
library("broom")
library ("plyr")

znorm<-function(x){
  k<-which(!is.na(x))
  y<-(x[k]-mean(x[k]))/sd(x[k])
  x[k]<-y
  x
}

NA_counts<-function(x){
  k<-which(!is.na(x))
  y<-length(k)
  y
}

nsnp_function<-function(x,df){
    nsnp<-sort(df[df$Biomarker==x,"Num_SNP"])
    if(length(nsnp)==1){
        return(nsnp)
    }else{
        return(nsnp[1])
    }
}

args<-commandArgs(T)
qtransdata<-args[1]
sample_size<-args[2]
tag<-args[3] #what type of data it is? use the directory name
stroke_type<-args[4]
outcome_type<-args[5]
IV_type<-args[6]
covariates<-args[7]
flag<-args[8] #whether the PRS was built based on different type or by the general population
nsnpdata<-args[9] #add
groupdata<-args[10]

for (ccs_type in c(1,3,12345)){
      prefix<-paste0(sample_size,"_",stroke_type,"_",outcome_type,"_",IV_type,"_",covariates,"_",flag,"_ccs",ccs_type)
      prs<-read.table("biomarker.55.h2_r2_rmout.txt",header=T)  #corresponding to sampleinfo.9739.PRS10X3.txt

if(outcome_type=="m3_stroke"){
    outcome_dd="m3_stroke_dd"
}else if(outcome_type=="y1_stroke"){
    outcome_dd="y1_stroke_dd"
}else if(outcome_type=="m3_death"){
    outcome_dd="m3_death_dd"
}else if(outcome_type=="y1_death"){
    outcome_dd="y1_death_dd"
}

if(flag=="All" | flag=="HST" | flag=="HSF"){

if(flag=="HST" && ccs_type!=12345){
        next;
}
    infile<-qtransdata
    data<-read.table(qtransdata,header=T,stringsAsFactors=F) #base
    nsnp_file<-nsnpdata
    nsnp_df<-read.table(nsnp_file,header=T,stringsAsFactors=F)
    nsnp_df2<-unique(nsnp_df[,c("Biomarker","Num_SNP")])
}else if(flag=="CCS" && ccs_type!="12345"){
    infile<-paste0(qtransdata,"/HSF_CCS",ccs_type,"/1e-6","/HSF_CCS",ccs_type,".1e-6.sampleinfo.PRS.qtrans.xls")
    data<-read.table(infile,header=T,stringsAsFactors=F)
    nsnp_file<-paste0(nsnpdata,"/HSF_CCS",ccs_type,"/1e-6","/HSF_CCS",ccs_type,".qtrans.summary.bygroup.xls")
    nsnp_df<-read.table(nsnp_file,header=T,stringsAsFactors=F)
    nsnp_df2<-unique(nsnp_df[,c("Biomarker","Num_SNP")])
}else if(flag=="CCS" && ccs_type=="12345"){
        next
}

if(flag=="All" | flag=="HST" | flag=="HSF"){
    if(ccs_type==12345){
      data<-data[data$IMG_ccs==1 | data$IMG_ccs==2 | data$IMG_ccs==3 | data$IMG_ccs==4 | data$IMG_ccs==5,]
    }else if(ccs_type==45){
      data<-data[data$IMG_ccs==4 | data$IMG_ccs==5,]
    }else{
      data<-data[data$IMG_ccs==ccs_type,]
    }
}
  
#riskfactors <- colnames(data)

if(T){
    if(flag=="All"){
    riskfactors<-c("BSL_Fib_PRS","BSL_DD_PRS")
    riskfactors1<-c("BSL_Fib","BSL_DD")
    colors<-c("#FF0000","#FF0000")
    Biomarkers<-c("Fib","D-dimer")
    }else if(flag=="HST"){
    riskfactors<-c("BSL_Fib_PRS","BSL_DD_PRS")
    riskfactors1<-c("BSL_Fib","BSL_DD")
    colors<-c("#FF0000","#FF0000")
    Biomarkers<-c("Fib","D-dimer")
    }else if(flag=="HSF"){
    #No.2 "BSL_AGEs_PRS",
    riskfactors<-c("BSL_Fib_PRS","BSL_DD_PRS")
    riskfactors1<-c("BSL_Fib","BSL_DD")
    colors<-c("#FF0000","#FF0000")
    Biomarkers<-c("Fib","D-dimer")
    }else if(flag=="CCS"){
        if(ccs_type==1){
            riskfactors<-c("BSL_Fib_PRS","BSL_DD_PRS")
            riskfactors1<-c("BSL_Fib","BSL_DD")
            colors<-c("#FF0000","#FF0000")
            Biomarkers<-c("Fib","D-dimer")
            groups<-c("BSL_DD_Group","BSL_Fib_Group")
         }else if(ccs_type==2){
            riskfactors<-c("BSL_Fib_PRS","BSL_DD_PRS")
            riskfactors1<-c("BSL_Fib","BSL_DD")
            colors<-c("#FF0000","#FF0000")
            Biomarkers<-c("Fib","D-dimer")
            groups<-c("BSL_DD_Group","BSL_Fib_Group")
        }else if(ccs_type==3){
            riskfactors<-c("BSL_Fib_PRS","BSL_DD_PRS")
            riskfactors1<-c("BSL_Fib","BSL_DD")
            colors<-c("#FF0000","#FF0000")
            Biomarkers<-c("Fib","D-dimer")
            groups<-c("BSL_DD_Group","BSL_Fib_Group")
        }else if(ccs_type==5){
            riskfactors<-c("BSL_Fib_PRS","BSL_DD_PRS")
            riskfactors1<-c("BSL_Fib","BSL_DD")
            colors<-c("#FF0000","#FF0000")
            Biomarkers<-c("Fib","D-dimer")
            groups<-c("BSL_DD_Group","BSL_Fib_Group")
        }
    }
}

if(T){  
  prs_t<-as.matrix(t(prs[,3]))
  colnames(prs_t)<-prs[,1]
  prs_df<-as.data.frame(prs_t)
  h2_t<-as.matrix(t(prs[,2]))
  colnames(h2_t)<-prs[,1]
  h2_df<-as.data.frame(h2_t)
}  
  
#"BMI","HEIGHT"
print(paste0("check data:",infile))
str(data)
data[,c("diff_MRS_state","m3_diff_MRS_state")][data[,c("diff_MRS_state","m3_diff_MRS_state")]==0]<-"NA"
data[,c("diff_MRS_state","m3_diff_MRS_state")][data[,c("diff_MRS_state","m3_diff_MRS_state")]==1]<-0
data[,c("diff_MRS_state","m3_diff_MRS_state")][data[,c("diff_MRS_state","m3_diff_MRS_state")]==2]<-1
data$m3_diff_MRS_state<-as.factor(data$m3_diff_MRS_state)
data$diff_MRS_state<-as.factor(data$diff_MRS_state)

if(T){
data[,c("m3_stroke","y1_stroke")][data[,c("m3_stroke","y1_stroke")]==0]<-"NA"
data[,c("m3_stroke","y1_stroke")][data[,c("m3_stroke","y1_stroke")]==1]<-0
data[,c("m3_stroke","y1_stroke")][data[,c("m3_stroke","y1_stroke")]==2]<-1
}

if(T){
data[,c("m3_death","y1_death")][data[,c("m3_death","y1_death")]==0]<-"NA"
data[,c("m3_death","y1_death")][data[,c("m3_death","y1_death")]==1]<-0
data[,c("m3_death","y1_death")][data[,c("m3_death","y1_death")]==2]<-1
}


if(T){
data[,c("PE","VTE","GIBLD","EP","INFECT")][data[,c("PE","VTE","GIBLD","EP","INFECT")]==0]<-"NA"
data[,c("PE","VTE","GIBLD","EP","INFECT")][data[,c("PE","VTE","GIBLD","EP","INFECT")]==1]<-0
data[,c("PE","VTE","GIBLD","EP","INFECT")][data[,c("PE","VTE","GIBLD","EP","INFECT")]==2]<-1
}

print(table(data[,"diff_MRS_state"]))
print(table(data[,"m3_diff_MRS_state"]))

#----
#compute R2
#----

Pheno_df<-data[,c(riskfactors1)]
PRS_df<-data[,c(riskfactors)]
R_square<-diag(outer(Pheno_df, PRS_df, function(X, Y){
mapply(function(...) cor.test(..., na.action = "na.exclude")$estimate^2,
           X, Y)
}))

R_square_P<-diag(outer(Pheno_df, PRS_df, function(X, Y){
mapply(function(...) cor.test(..., na.action = "na.exclude")$p.value,
           X, Y)
}))


data[,c(riskfactors)]=apply(data[,c(riskfactors)],2,znorm)
print(paste0("ccs",ccs_type,"check data after extracting riskfactors"))
str(data$m3_diff_MRS_state)
effective_size=apply(data[,riskfactors],2,NA_counts)
nsnp=sapply(riskfactors1,nsnp_function,nsnp_df2)
case_size=length(data[ ,outcome_type][data[,outcome_type]==1])
control_size=length(data[ ,outcome_type][data[,outcome_type]==0])
prs_value<-prs_df[,c(riskfactors1)]
h2_value<-h2_df[,c(riskfactors1)]



if(outcome_type=="m3_diff_MRS_state" | outcome_type == "diff_MRS_state" | outcome_type == "PE" | outcome_type == "VTE" | outcome_type == "GIBLD" | outcome_type == "EP" | outcome_type == "INFECT"){
     group_function<-function(biomarker){
                            index<-match(biomarker,riskfactors)
                            group_col<-groups[index]
                            data<-data[,c(biomarker,group_col,outcome_type,"GENDER","AGE","PCA1","PCA2","PCA3","PCA4","PCA5","H_STROKE01","IMG_ccs","onset_multic_day")]
                            dots <- lapply(group_col, as.symbol)
                            return(data %>% group_by(.dots=dots) %>% group_map(~ glm(as.formula(paste0(outcome_type,"~GENDER+AGE+PCA1+PCA2+PCA3+PCA4+PCA5+H_STROKE01+IMG_ccs+onset_multic_day+",biomarker)),family=binomial(link='logit'),data=.x)))
                          }

    univ_models <- lapply(riskfactors,group_function)

    univ_results <- lapply(univ_models,function(x){
                            return(lapply(x,function(x){
                                 x <- summary(x)
                                 biomarker<-dimnames(x$coefficients)[[1]][12]
                                 p.value<-signif(rev(x$coefficients[,"Pr(>|z|)"])[1], digits=3)
                                 wald.test<-signif(rev(x$coefficients[,"z value"][1]), digits=3)
                                 BETA<-signif(rev(x$coefficients[,"Estimate"])[1], digits=3);#coeficient beta
                                 SE<-signif(rev(x$coefficients[,"Std. Error"])[1], digits=3);#coeficient beta
                                 HR <-signif(exp(rev(x$coefficients[,"Estimate"])[1]), digits=3);#exp(beta)
                                 BETA.confint.lower <- signif(BETA-2*SE,3)
                                 BETA.confint.upper <- signif(BETA+2*SE,3)
                                 HR.confint.lower <- signif(exp(BETA.confint.lower),3)
                                 HR.confint.upper <- signif(exp(BETA.confint.upper),3)
                                 res<-c(biomarker,BETA, SE,BETA.confint.lower,BETA.confint.upper,HR, HR.confint.lower,HR.confint.upper,wald.test, p.value)
                                 names(res)<-c("Biomarker","BETA","SE","BETA.confint.lower","BETA.confint.upper", "HR", "HR.confint.lower","HR.confint.upper", "wald.test","p.value")
                                 return(res)
                                }))
                       })
}else{
    group_function<-function(biomarker){
                            index<-match(biomarker,riskfactors)
                            group_col<-groups[index]
                            data<-data[,c(biomarker,group_col,outcome_dd,outcome_type,"GENDER","AGE","PCA1","PCA2","PCA3","PCA4","PCA5","H_STROKE01","IMG_ccs","onset_multic_day")]
                            dots <- lapply(group_col, as.symbol)
                            return(data %>% group_by(.dots=dots) %>% group_map(~ coxph(as.formula(paste0("Surv(",outcome_dd, ",",outcome_type, ")~GENDER+AGE+PCA1+PCA2+PCA3+PCA4+PCA5+H_STROKE01+IMG_ccs+onset_multic_day+",biomarker)),data=.x)))
#                           group_map(~ coxph(as.formula(paste0("Surv(",outcome_dd, ",",outcome_type, ")~GENDER+AGE+PCA1+PCA2+PCA3+PCA4+PCA5+H_STROKE01+IMG_ccs+",biomarker)),data=.x))
                          }

     univ_models <- lapply(riskfactors,group_function)
#     str(univ_models)

    # Extract data
    univ_results <- lapply(univ_models,function(x){
                            return(lapply(x,function(x){
                                 print("first x")
                                 str(x)
                                 x <- summary(x)
                                 print("first summary(x)")
                                 str(x)
                                 print(dimnames(x$coefficients))
                                 biomarker<-dimnames(x$coefficients)[[1]][12]
                                 print(biomarker)
                                 p.value<-signif(rev(x$coefficients[,"Pr(>|z|)"])[1], digits=3)
                                 wald.test<-signif(rev(x$coefficients[,"z"][1]), digits=3)
                                 BETA<-signif(rev(x$coefficients[,"coef"])[1], digits=3);#coeficient beta
                                 SE<-signif(rev(x$coefficients[,"se(coef)"])[1], digits=3);#coeficient beta
                                 HR <-signif(rev(x$coefficients[,"exp(coef)"])[1], digits=3);#exp(beta)
                                 HR.confint.lower <- signif(rev(x$conf.int[,"lower .95"])[1], 3)
                                 HR.confint.upper <- signif(rev(x$conf.int[,"upper .95"])[1],3)
                                 BETA.confint.lower <- signif(rev(log(x$conf.int[,"lower .95"])[1]),3)
                                 BETA.confint.upper <- signif(rev(log(x$conf.int[,"upper .95"])[1]),3)

                                 res<-c(biomarker,BETA, SE,BETA.confint.lower,BETA.confint.upper,HR, HR.confint.lower,HR.confint.upper,wald.test, p.value)
                                 names(res)<-c("Biomarker","BETA","SE","BETA.confint.lower","BETA.confint.upper", "HR", "HR.confint.lower","HR.confint.upper", "wald.test","p.value")
                                 return(res)
                                }))
                            
                       })    
}

write.table(prefix,paste0(prefix,".meta.detail.txt"),quote=F,sep="\t",row.names=F,append=F) #change 2
all_meta_list<-lapply(univ_results,function(x){
    res <- t(as.data.frame(x, check.names = FALSE))
    df<-as.data.frame(res,stringsAsFactors = F)
#   df<-df[order(df$p.value),]
    print(paste(prefix,"df"))
    print(df)
    df[is.na(df)]<-0
    df[df==0]<-0.00001
    meta_out<-rma.uni(yi=as.numeric(df$BETA),sei=as.numeric(df$SE),control=list(maxiter=100000,stepadj=0.5))
    print(paste(prefix,"str meta_out"))
    print(str(meta_out))
if(T){
     meta_df<-data.frame("Biomarker"=df$Biomarker[1],"BETA"=signif(as.numeric(meta_out$beta[1,1]),digits=3),
                         "SE"=signif(meta_out$se,digits=3),
                         "BETA.confint.lower"=signif(meta_out$ci.lb,digits=3),
                         "BETA.confint.upper"=signif(meta_out$ci.ub,digits=3),
                         "HR"=signif(exp(as.numeric(meta_out$beta[1,1])),digits=3),
                         "HR.confint.lower"=signif(exp(meta_out$ci.lb),digits=3),
                         "HR.confint.upper"=signif(exp(meta_out$ci.ub),digits=3),
                         "wald.test"=signif(meta_out$zval,3),"p.value"=signif(meta_out$pval,3))
     df<-rbind(df,meta_df)

#    meta_df<-ldply (meta_out, data.frame)
#    print(paste(prefix,"<F3>meta_df"))<F4><F3>
#    print(meta_df)
}

    print(paste(prefix,"str meta_df"))
    str(meta_df)
    write.table(df,paste0(prefix,".meta.detail.txt"),quote=F,sep="\t",row.names=F,append=T) #change 2
    return(meta_df)
})

print(paste(prefix,"str all_meta_df"))
str(all_meta_list)


all_meta_df = do.call(what = rbind, args = all_meta_list)
print(paste(prefix,"str all_meta_df as.data.frame"))
str(all_meta_df)
str(effective_size)

    all_meta_df$size<-effective_size
    all_meta_df$nsnp<-nsnp
    all_meta_df$case<-case_size
    all_meta_df$control<-control_size
    all_meta_df$prs<-as.numeric(prs_value)
    all_meta_df$h2<-as.numeric(h2_value)
    all_meta_df$p.category<-ifelse(all_meta_df$p.value<0.05,paste0("\"","#55b0f6","\""),paste0("\"","#132b43","\""))
    all_meta_df$xlabel<-paste0("\"",colors,"\"")
    all_meta_df$Biomarker<-Biomarkers
    all_meta_df$R2<-R_square
    all_meta_df$R2_P<-R_square_P
print(paste(prefix,"str all_meta_df as.data.frame add effective size"))
str(all_meta_df)
    all_meta_df$CCStype<-paste0(flag,"_CCS",rep(ccs_type,length(Biomarkers)))

print(paste(prefix,"str all_meta_df as.data.frame add info"))
str(all_meta_df)
write.table(all_meta_df,paste0(prefix,".meta.txt"),quote=F,sep="\t",row.names=F,append=F)


if(T){
#making forest plot
forest(
  x=as.numeric(df$HR),
  ci.lb=as.numeric(df$HR.confint.lower),
  ci.ub=as.numeric(df$HR.confint.upper),
  
  #x=as.numeric(df$BETA),
  #ci.lb=as.numeric(df$BETA.confint.lower),
  #ci.ub=as.numeric(df$BETA.confint.upper),
  
  slab = sprintf("    %s", rownames(df)),
  xlab = "HR for m3_death",
  annotate = FALSE,
  ilab = data.frame(
    sprintf("%.0f", as.numeric(df$size)),
    sprintf("%.3f", as.numeric(df$h2)),
    sprintf("%.3f", as.numeric(df$prs)),
    sprintf("%.2f", as.numeric(df$HR)),
    sprintf("(%.2f, %.2f)", as.numeric(df$HR.confint.lower), as.numeric(df$HR.confint.upper)),
    #sprintf("%.2f", as.numeric(df$BETA)),
    #sprintf("(%.2f, %.2f)", as.numeric(df$BETA.confint.lower), as.numeric(df$BETA.confint.upper)),
    sprintf("%.2e", as.numeric(df$p.value))),
  ilab.xpos = c(4, 5.5, 7, 8.5, 10,12),
  pch = 16,
  #atransf = exp,
  #at = log(c(0.5, 1, 2, 4, 8, 16)),
  #at = c(0,1,2,3,4,5,6),
  at=c(0,0.5,1,1.5,2,2.5,3),
  #rows = rev(c(0:2, 4:10, 12:13, 15:18, 20:23, 25:28,30:32)),
  xlim = c(-5,14),
  ylim = c(0, length(riskfactors)+3), #number of biomarkers plus 3
  refline=1
)

title(main=paste0("CCS",ccs_type,": ",stroke_type,": Logit(",outcome_type,")~GENDER+AGE+onset_multic_day+PCA1+PCA2+PCA3+PCA4+PCA5+H_STROKE01+IMG_ccs")) #change 4
text(c(4, 5.5, 7, 8.5, 10,12),length(riskfactors)+2.5,pos=1,c("Size","H2(GCTA)","PRS_R2","HR","95% CI", "P"))
}
dev.off()

}
}
