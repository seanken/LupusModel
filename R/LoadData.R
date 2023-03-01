library(DropletUtils)


loadDrops<-function(lst,numCells)
{
print("Get QC")
fils_h5=sub("$","/molecule_info.h5",lst)
fils_qc=sub("$","/metrics_summary.csv",lst)
fils_cbc=sub("$","/filtered_feature_bc_matrix/barcodes.tsv.gz",lst)
fils_genes=sub("$","/filtered_feature_bc_matrix/features.tsv.gz",lst)
out=lapply(fils_qc,read.csv)
out=lapply(out,function(x){as.numeric(gsub(",","",x[1,2]))})
out=as.numeric(out)
rat=min(out)/out
names(rat)=names(lst)

print("Downsample")
out=lapply(1:length(lst),function(x){
print(x)
ds=rat[x]
h5=fils_h5[x]
nam=names(lst)[x]
print("Read genes!")
genes=read.table(fils_genes[x],sep="\t",stringsAsFactors=F)
#print("Get cells")
#cells=scan(fils_cbc[x],"")
print("DS")
dat=downsampleReads(h5,prop=ds)
print("Only cells")
tab=emptyDropsCellRanger(dat,n.expected.cells=numCells[x])
cells=rownames(tab)[tab$FDR <= .01 & !is.na(tab[,"FDR"]) ]
dat=dat[,cells]
colnames(dat)=sub("^",paste(nam,"_",sep=""),colnames(dat))
print("Rename genes")
rownames(genes)=genes[,1]
genes["genes"]=make.unique(genes[,2])
rownames(dat)=genes[rownames(dat),"genes"]
print(dim(dat))
return(dat)
})
print("Combine!")
dat=do.call(cbind,out)
print("Return!")
return(dat)
}
