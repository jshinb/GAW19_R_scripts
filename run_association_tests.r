# --- Loading necessary packages and functions ---- #
require(stringr,lib.loc="/home/bulllab/jshin/R/x86_64-unknown-linux-gnu-library/3.0",quietly=T)

require(SKAT,lib.loc="/home/bulllab/jshin/R/x86_64-unknown-linux-gnu-library/3.0",quietly=T)
# source SKAT-functions modified in order to extract extra information
# (such as kurtosis or df estimates)
for (i in 1:length(dir('/home/bulllab/gaw18/gaw19/jshin/scripts/SKAT_R/',full.name=TRUE))) {
  source(dir('/home/bulllab/gaw18/gaw19/jshin/scripts/SKAT_R/',full.name=TRUE)[i])
}
assignInNamespace("SKAT_PValue_Logistic_VarMatching",
                  SKAT_PValue_Logistic_VarMatching_JS,ns="SKAT")
assignInNamespace("KMTest.logistic.Linear.VarMatching",
                  KMTest.logistic.Linear.VarMatching_JS,ns="SKAT")
assignInNamespace("SKAT_Get_DF_Sim",
                  SKAT_Get_DF_Sim_JS,ns="SKAT")
assignInNamespace("SKAT_GET_kurtosis",
                  SKAT_GET_kurtosis_JS,ns="SKAT")
assignInNamespace("SKAT_Logistic_VarMatching_GetParam1_QuantileAdj",
                  SKAT_Logistic_VarMatching_GetParam1_QuantileAdj_JS,ns="SKAT")
assignInNamespace("SKAT_Logistic_VarMatching_GetParam",
                  SKAT_Logistic_VarMatching_GetParam_JS,ns="SKAT")

# source penalized-likelihood approach implemented by a previous trainee(NAME?) of SB
for (i in 1:length(dir('/home/bulllab/jshin/pmlr11/R/',full.name=TRUE))){
source(dir('/home/bulllab/jshin/pmlr11/R/',full.name=TRUE)[i])
}
rm(i)

out.file = "/home/bulllab/gaw18/gaw19/results/chr3_MAP4_res_no_imputation_qsub.out"
# ------------------ read in data ------------------#
# (phenotype + genotype - created by 'Make_analysible_data.Rnw'
# in JS's desktop-should mv file in the future)
# each genetic marker column codes for the number of index allele
data <- read.csv('/home/bulllab/gaw18/gaw19/data/chr3_pheno_MAP4_var_sites_MAF.csv',
                 header=T, stringsAsFactors=F)
print(head(data)) 
print(dim(data))

marker = names(data)[-c(1:9)]
rs_dist = read.csv('/home/bulllab/gaw18/gaw19/data/chr3_MAP4_positions_from_snpnexus_30132.csv',header=T,stringsAsFactors=F)
rs_dist = rs_dist[,c("SNP","chromPosition")]

dist = str_replace(marker,"var_3_","")

for(i in 1:nrow(rs_dist)){
  dist[which(dist==rs_dist$SNP[i])] <- rs_dist$chromPosition[i]  
}
dist <- as.numeric(dist) #34 of them do not have distances

map.info <- cbind.data.frame(marker,dist,stringsAsFactors=FALSE)
rm(marker,dist)

# -------------- analysis begins here --------------#
# column numbers where marker data are included
first.G.col = 10
last.G.col = ncol(data) #87 polymorphic markers

#25 output variables
                  
# convert the numeric 'hypt' column into factor for pmlr() function
data$y <- as.factor(data$hypt)
org.data <- data 
rm(data)

for(i in first.G.col:last.G.col){
  marker = names(org.data)[i]
  pos = map.info$dist[(i-9)]
  cat((i-9),'-th marker ', marker, "\n", sep="")

  #missing genotype rates - in the data set with complet info on hypt.
  missing_geno_rate = (sum(is.na(org.data[!is.na(org.data$y),marker]))/sum(!is.na(org.data$y))) 

  # creating a dataset with complete information 
  # to prevent SKAT from imputing missing genotypes
  no.missing.ind <- !is.na(org.data$y) & !is.na(org.data[,marker])
  data = org.data[no.missing.ind,]

  n00 = sum(data$hypt==0 & data[,i]==0,na.rm=T)
  n01 = sum(data$hypt==0 & data[,i]==1,na.rm=T)
  n02 = sum(data$hypt==0 & data[,i]==2,na.rm=T)
  n10 = sum(data$hypt==1 & data[,i]==0,na.rm=T)
  n11 = sum(data$hypt==1 & data[,i]==1,na.rm=T)
  n12 = sum(data$hypt==1 & data[,i]==2,na.rm=T)
  n = n00+n01+n02+n10+n11+n12
  allele.freq = (n01+n11+2*(n02+n12))/n

  if(allele.freq == 0){
    beta_MLE <- SE.beta_MLE <- beta_PMLE <- SE.beta_PMLE <- NA
    stat_lrt <- stat_plrt <- NA
    stat_score <- stat_score_var_adj <- stat_score_var_kurt_adj <- NA
    df_score_var_kurt_adj <- NA
    pval_lrt <- pval_plrt <- pval_score <- NA
    pval_score_var_adj <- pval_score_var_kurt_adj <- NA
  }

  if(allele.freq > 0 ){
  # applying MLE
  # fitting an additive model
  reg.model = as.formula(paste('y~',marker))
  print(reg.model)

  # standard likelihood ratio test
  lrt.res = pmlr(reg.model, data=data, method="likelihood", penalized=F)
  stat_lrt = lrt.res$stat[,marker,TRUE]
  pval_lrt = lrt.res$pval[,marker,TRUE]
  beta_MLE = lrt.res$coef[,,TRUE][marker]
  SE.beta_MLE = sqrt(lrt.res$var[marker,marker,TRUE])
    
  # penalized likelihood ratio test
  plrt.res = pmlr(reg.model, data=data, method="likelihood", penalized=T)
  stat_plrt = plrt.res$stat[,marker,TRUE]
  pval_plrt = plrt.res$pval[,marker,TRUE]
  beta_PMLE = plrt.res$coef[,,TRUE][marker]
  SE.beta_PMLE = sqrt(plrt.res$var[marker,marker,TRUE])

  # standard score test
  score.res = pmlr(reg.model, data=data, method="score", penalized=F)
  stat_score = score.res$stat[,marker,TRUE]  
  pval_score = score.res$pval[,marker,TRUE]

  do.not.run <- function(){
  #comparison.score
    lm0 = glm(y~1,data=data[!is.na(data[,i]),],family=binomial())
    lm1 = glm(reg.model,data=data[!is.na(data[,i]),],family=binomial())
    anova(lm0,lm1,test="Rao")
    anova(lm0,lm1,test="Chi")
  }
  
  #SKAT - small-sample-adjustments to var or to var and kurtosis
  Z <- as.matrix(data[,i],ncol=1)
  # var-adj test
  obj.kurtosis.adj <- SKAT_Null_Model_MomentAdjust(hypt~1, data=data)
  obj.no.kurtosis.adj <- SKAT_Null_Model_MomentAdjust(hypt~1, data=data,is_kurtosis_adj=FALSE)

  #cutoff of the missing rates of the SNPs - default is 15%
  #not sure if I remove SNPs with missing rate >= 15%
  missing_rate_threshold = 1 # not filtering anything
  
  #score test with small-sample-adjusted variance  
  skat.no.kurtosis.adj <- SKAT(Z,obj.no.kurtosis.adj,weights=1,
                               missing_cutoff=missing_rate_threshold,estimate_MAF=2)
  stat_score_var_adj <- sqrt(2*skat.no.kurtosis.adj$param$df)*(skat.no.kurtosis.adj$Q-skat.no.kurtosis.adj$param$muQ)/sqrt(skat.no.kurtosis.adj$param$varQ)+skat.no.kurtosis.adj$param$df
  pval_score_var_adj <- skat.no.kurtosis.adj$p.value


  #score test with small-sample-adjusted variance and kurtosis
  skat.kurtosis.adj <- SKAT(Z,obj.kurtosis.adj,weights=1,
                            missing_cutoff=missing_rate_threshold,estimate_MAF=2)
  df_score_var_kurt_adj <- skat.kurtosis.adj$param$df
  stat_score_var_kurt_adj <- sqrt(2*df_score_var_kurt_adj)*(skat.kurtosis.adj$Q-skat.kurtosis.adj$param$muQ)/sqrt(skat.kurtosis.adj$param$varQ)+df_score_var_kurt_adj
  pval_score_var_kurt_adj <- skat.kurtosis.adj$p.value

  }

  res <- cbind.data.frame(marker,pos,allele.freq,
    n00,n01,n02,n10,n11,n12,n,missing_geno_rate,
    beta_MLE,SE.beta_MLE,beta_PMLE,SE.beta_PMLE,
    stat_lrt,stat_plrt,
    stat_score,stat_score_var_adj,stat_score_var_kurt_adj,
    df_score_var_kurt_adj,
    pval_lrt,pval_plrt,pval_score,pval_score_var_adj,pval_score_var_kurt_adj)

  if(i==first.G.col){
    #printing column names
    #do this once - using the default sep=" " 
    write.table(rbind(names(res)),file=out.file,
    quote=F,col.names=F,row.names=F,append=TRUE,sep=" ")
  }

  write.table(rbind(res),file=out.file,
    quote=F,col.names=F,row.names=F,append=TRUE,sep=" ")

}

write.table(map.info, file = "/home/bulllab/gaw18/gaw19/results/chr3_MAP4.map",
  quote=F,col.names=T,row.names=F,sep=" ")

