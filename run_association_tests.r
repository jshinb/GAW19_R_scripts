# --- Loading necessary packages and functions ---- #
#install.packages("SKAT")

library(SKAT)

# source SKAT-functions modified in order to extract extra information
# (such as kurtosis or df estimates)
for (i in 1:length(dir('/home/bulllab/gaw18/gaw19/jshin/scripts/R/',full.name=TRUE))) {
  source(dir('/home/bulllab/gaw18/gaw19/jshin/scripts/R/',full.name=TRUE)[i])
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

# ------------------ read in data ------------------#
# (phenotype + genotype - created by 'Make_analysible_data.Rnw'
# in JS's desktop-should mv file in the future)
# each genetic marker column codes for the number of index allele
data <- read.csv('/home/bulllab/gaw18/gaw19/data/chr3_pheno_MAP4_var_sites_MAF_0.01.csv',
                 header=T, stringsAsFactors=F)
print(head(data)) #with 10 variants (with MAF>=0.01

# -------------- analysis begins here --------------#
# column numbers where marker data are included
first.G.col = 10
last.G.col = 19

#25 output variables
                  
pos <- NA #not-available yet
# convert the numeric 'hypt' column into factor for pmlr() function
data$y <- as.factor(data$hypt)
i=12
for(i in first.G.col:last.G.col){
  
  n00 = sum(data$hypt==0 & data[,i]==0,na.rm=T)
  n01 = sum(data$hypt==0 & data[,i]==1,na.rm=T)
  n02 = sum(data$hypt==0 & data[,i]==2,na.rm=T)
  n10 = sum(data$hypt==1 & data[,i]==0,na.rm=T)
  n11 = sum(data$hypt==1 & data[,i]==1,na.rm=T)
  n12 = sum(data$hypt==1 & data[,i]==2,na.rm=T)
  n = n00+n01+n02+n10+n11+n12
  missing_rate = 1-n/nrow(data)
  allele.freq = (n01+n11+2*(n02+n12))/n

  # applying MLE
  # fitting an additive model
  marker = names(data)[i]
  reg.model = as.formula(paste('y~',marker))
  print(reg.model)

  # standard likelihood ratio test
  lrt_MLE.res = pmlr(reg.model, data=data, method="likelihood", penalized=F)
  beta_MLE = lrt_MLE.res$coef[,,TRUE][marker]
  SE.beta_MLE = sqrt(lrt_MLE.res$var[marker,marker,TRUE])
  pval_beta_MLE = lrt_MLE.res$pval[,marker,TRUE]
    
  # penalized likelihood ratio test
  plrt_MLE.res = pmlr(reg.model, data=data, method="likelihood", penalized=T)
  beta_PMLE = plrt_MLE.res$coef[,,TRUE][marker]
  SE.beta_PMLE = sqrt(plrt_MLE.res$var[marker,marker,TRUE])
  pval_beta_PMLE = plrt_MLE.res$pval[,marker,TRUE]

  # standard score test
  score_MLE.res = pmlr(reg.model, data=data, method="score", penalized=F)
  stat_score_MLE = score_MLE.res$test$statistic[2]
  pval_score_MLE = score_MLE.res$test$pvalue[2]

  do.not.run <- function(){
  #comparison.score
    lm0 = glm(y~1,data=data[!is.na(data[,i]),],family=binomial())
    lm1 = glm(reg.model,data=data[!is.na(data[,i]),],family=binomial())
    anova(lm0,lm1,test="Rao")
  }
  
  #SKAT 
  Z <- as.matrix(data[,i],ncol=1)
  # var-adj test
  obj.kurtosis.adj <- SKAT_Null_Model_MomentAdjust(hypt~1, data=data)
  obj.no.kurtosis.adj <- SKAT_Null_Model_MomentAdjust(hypt~1, data=data,is_kurtosis_adj=FALSE)

  #cutoff of the missing rates of the SNPs - default is 15%
  #not sure if I remove SNPs with missing rate >= 15%
  missing_rate_threshold = 0.5
  
  skat.kurtosis.adj <- SKAT(Z,obj.kurtosis.adj,weights=1,
                            missing_cutoff=missing_rate_threshold,estimate_MAF=2)
  skat.no.kurtosis.adj <- SKAT(Z,obj.no.kurtosis.adj,weights=1,
                               missing_cutoff=missing_rate_threshold,estimate_MAF=2)
  
}

res_colnames <- c("marker","pos","allele.freq",
                  "n00","n01","n02","n10","n11","n12","n", #cell-numbers
                  "beta_MLE","SE.beta_MLE","beta_PMLE","SE.beta_PMLE",#genetic effect estimates
                  "stat_lrt_MLE","stat_plrt_MLE",#LRT test statistics
                  "stat_score_MLE",
                  "stat_SKAT_score_var_adj","stat_SKAT_score_var_kurt_adj",#score statistics
                  "df_SKAT_score_var_kurt_adj",#df estimates from SKAG-var-kurt-adj
                  "pval_lrt_MLE","pval_plrt_MLE", #LRT-pvalues
                  "pval_score_MLE","pval_SKAT_score_var_adj","pval_SKAT_score_var_kurt_adj")#score statistics pvalues
                  
