marker_dist = read.table('/home/bulllab/gaw18/gaw19/data/chr3-dose_markers.out',header=T,sep="_",stringsAsFactors=FALSE)

gene.name="MAP4"
gene.start = 47892000
gene.end = 48131000

ind = marker_dist$POS >= gene.start & marker_dist$POS <= gene.end #409 markers
range(which(ind))#50220 50628


# after running the above lines
# sed -n ' 1,1p;50220,50628p;50628q' chr3-dose.csv > chr3_MAP4-dose.csv #410 lines including headers
