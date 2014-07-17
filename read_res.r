# The file is used for summarizing association tests between 
# hypertension status genetic markers in MAP4 region on chromosome 3
# Do this on Desktop

data = read.csv('~/Desktop/GAW19/data/proc/chr3_pheno_MAP4_var_sites_MAF.csv',
	header=T,stringsAsFactors=F)
res = read.table('~/Desktop/GAW19/results/chr3_MAP4_res_no_imputation_qsub.out',header=T,stringsAsFactors=F,sep=" ")

# LDplot
require(LDheatmap)
require(genetics)
marker.dat = data[data$hypt==0,-c(1:9)]
gdat = NULL
i=1
gdat <- as.genotype.allele.count(marker.dat[,i],alleles=c("0","1"))
for(i in 2:ncol(marker.dat)){
gdat = cbind.data.frame(gdat,as.genotype.allele.count(marker.dat[,i],alleles=c("0","1")))	
}
names(gdat) <- names(marker.dat)
gdat.dist <- as.vector(res$pos)
gdat.dist[73] <- mean(res$pos[72:74],na.rm=T)

#which one are monomorphic in controls?
#var_3_47894706
#var_3_47894808
#var_3_47898875
#var_3_47898934
#rs150980381
#var_3_47913606
#var_3_47913624
#var_3_47951185
#var_3_47957635
#var_3_47957709
#rs76333690
#var_3_47969700
#rs76535647
#var_3_48016938
mono.controls = res$n01+res$n02 == 0
my.heatmap <- LDheatmap(gdat[,!mono.controls],genetic.dist=gdat.dist[!mono.controls])

pdf('chr3_MAP4_pvals.pdf',width=11,height=4)
par(mfrow=c(1,5))
with(res,plot(pos/10^6,-log10(pval_lrt),pch=20,type="p",
	ylab="-log10 (P-value) ",xlab="",ylim=c(0,4)))
title('Standard likelihood Ratio test')

with(res,plot(pos/10^6,-log10(pval_plrt),pch=20,type="p",
	ylab=" ",xlab="",ylim=c(0,4)))
title('Penalized likliehood test')

with(res,plot(pos/10^6,-log10(pval_score),pch=20,type="p",
	xlab="chr3: Position (Mb)",ylim=c(0,4)))
title('Standard score test')

with(res,plot(pos/10^6,-log10(pval_score_var_adj),pch=20,type="p",
	ylab=" ",xlab="",ylim=c(0,4)))
title('Score test with small-sample-adjusted variance')

with(res,plot(pos/10^6,-log10(pval_score_var_kurt_adj),pch=20,type="p",
	ylab=" ",xlab="",ylim=c(0,4)))
title('Score test with small-sample-adjusted variance and kurtosis')
dev.off()


pdf('chr3_MAP4_pvals_pairwise.pdf',width=8,height=8)
par(mfrow=c(2,2))
with(res,plot(-log10(pval_score_var_kurt_adj),-log10(pval_lrt),pch=20,type="p",
	ylab="standard log likelihood",
	xlab="score test with adjusted variance and kurtosis",ylim=c(0,4)));
abline(0,1,col="red")


with(res,plot(-log10(pval_score_var_kurt_adj),-log10(pval_plrt),pch=20,type="p",
	ylab="penalized log likelihood",
	xlab="score test with adjusted variance and kurtosis",ylim=c(0,4)));
abline(0,1,col="red")

with(res,plot(-log10(pval_score_var_kurt_adj),-log10(pval_score),pch=20,type="p",
	ylab="score test",
	xlab="score test with adjusted variance and kurtosis",ylim=c(0,4)));
abline(0,1,col="red")

with(res,plot(-log10(pval_score_var_kurt_adj),-log10(pval_score_var_adj),pch=20,type="p",
	ylab="score test with adjusted variance",
	xlab="score test with adjusted variance and kurtosis",ylim=c(0,4)));
abline(0,1,col="red")
dev.off()
mean(data$hypt)#0.26
