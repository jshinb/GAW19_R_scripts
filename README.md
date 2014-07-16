GAW19_R_scripts
===============
* These are the scripts created to analyses GAW19 data comparing sparse-data and standard likelihood methods.

=== To-Do ===
* Need to make '.pbs' file to submit the job on cluster
* Need to 'cat' the column names to 
  '/home/bulllab/gaw18/gaw19/results/chr3_MAP4_res.out'

=== To-Do (july 16) ===
* Need to think about the 'imputation' in SKAT - leave as is since it is the default option for SKAT??

According to SKAT manual - 
'impute.method' 
a method to impute missing genotypes (default= "fixed"). 
"random" imputes missing genotypes by generating 
binomial(2,p) random variables (p is the MAF), 
and "fixed" imputes missing genotypes by assigning 
the mean genotype value (2p). 
If you use "random", you will have different p-values for 
different runs because imputed values are randomly assigned.
