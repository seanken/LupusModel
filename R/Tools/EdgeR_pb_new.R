library(Seurat)
library(Matrix)
library(edgeR)
library(fdrtool)

##modified from seurat code
combData<-function(object,genes=c(),assay="RNA",DSTo=0)
{
genes.use=rownames(object[[assay]]@counts)
#genes.use=rownames(object@assays$RNA@data)

if(length(genes)>0){genes.use=genes}
data.all=data.frame(row.names = genes.use)
levs=c()
for(i in levels(object@active.ident)) {
temp.cells=WhichCells(object,idents=i)
if(DSTo>0)
{
temp.cells=sample(temp.cells,DSTo)
}
levs=c(levs,i)
if (length(temp.cells)==1) data.temp=(object[[assay]]@counts[genes.use,temp.cells])
#if (length(temp.cells)==1) data.temp=(object@assays$RNA@counts[genes.use,temp.cells])
if (length(temp.cells)>1) data.temp=apply(object[[assay]]@counts[genes.use,temp.cells],1,sum)
#if (length(temp.cells)>1) data.temp=apply(object@assays$RNA@counts[genes.use,temp.cells],1,sum)
data.all=cbind(data.all,data.temp)
colnames(data.all)[ncol(data.all)]=i

}
colnames(data.all)=levs
return(data.all)
}




##modified from seurat code
combData_mean<-function(object,genes=c(),assay="RNA",DSTo=1)
{
genes.use=rownames(object[[assay]]@counts)
#genes.use=rownames(object@assays$RNA@data)

if(length(genes)>0){genes.use=genes}
data.all=data.frame(row.names = genes.use)
levs=c()
for(i in levels(object@active.ident)) {
temp.cells=WhichCells(object,idents=i)
levs=c(levs,i)
if (length(temp.cells)==1) data.temp=(object[[assay]]@counts[genes.use,temp.cells])
#if (length(temp.cells)==1) data.temp=(object@assays$RNA@counts[genes.use,temp.cells])
if (length(temp.cells)>1) data.temp=apply(object[[assay]]@counts[genes.use,temp.cells],1,function(x){floor(DSTo*mean(x))})
#if (length(temp.cells)>1) data.temp=apply(object@assays$RNA@counts[genes.use,temp.cells],1,sum)
data.all=cbind(data.all,data.temp)
colnames(data.all)[ncol(data.all)]=i

}
colnames(data.all)=levs
return(data.all)
}




######
##Performs DE by combining cells in same batch
######
##Input:
##seur: seurat object, with count data stored in seur@raw.data, condition information (ko_vs_wt) is stored in seur@data.info["ko_vs_wt"]. Assume two conditions, one ko, one for wt
##id: the cluster (or list of clusters) of interest
##kowt: the name of the column in the meta.data table containing ko vs wt information
##combineOn: gives the name of the column in seur@meta.data with batch information (so one factor per batch)
##method: method used for DE. Options: EdgeR, DESeq,
##form: usually leave blank. Can use to include extra covariate infomation. NOT IMPLEMENTED YET FOR edgeR.
##minCells: min number of cells per batch for batch to be included.
##minBatches: minimum number of batches per condition for analysis to run (if not met returns NULL).
##minReads: minimum number of reads in genes (requires this number of reads in at least 2 organoids)
##
##Note: ignore genes that do not have at least 10 reads mapped to them in at least two batches.
###########################
##Output:
##res: a data frame with DE information (pvalue, FDR corrected pvalue, logFC, etc). Exact formatting depends on method selected.
##########################################
combDE_edgeR<-function(seur,cell_id=c(),kowt="ko_vs_wt",combineOn="orig.ident",form=~condition,minCells=20,minBatches=2,minReads=10,contrast=NULL,coef=2,assay="RNA",combmethod="sum",useLRT=T,performDS=F,trend.method="locfit",tagwise=T,getPseudo=F,prior.df=NULL,robust=F)
{
if(length(cell_id)>0)
{
print("subsample")
seur<-subset(seur,cells=names(seur@active.ident)[seur@active.ident %in% cell_id])
}


seur<-SetIdent(seur,value=combineOn)
DSTo=min(as.numeric(table(seur@active.ident)))
print("combine data")
dat=c()
if(combmethod=="sum")
{
if(!performDS){DSTo=0}
dat<-combData(seur,assay=assay,DSTo=DSTo)
}
if(combmethod=="mean")
{
dat<-combData_mean(seur,DSTo=DSTo,assay=assay)
}

print("Filter 10X Lanes")
keep=names(summary(seur@active.ident))[summary(seur@active.ident)>minCells]
numOrg=length(keep)
dat=dat[keep]
print(dim(dat))
print("get ko vs wt data")
print(head(dat))
cols=colnames(dat)
print(cols)
lst=c()
lst=cols


keep<-strsplit(as.character(form)[2],"+",fixed=T)[[1]]
keep=trimws(keep)
keep=setdiff(keep,kowt)
keep=setdiff(keep,"condition")

print(c(kowt,combineOn,keep))
print(head(seur@meta.data))

#keep<-strsplit(form,"+",fixed=T)[[1]]
val=seur@meta.data[,c(kowt,combineOn,keep)]
val=val[!duplicated(val[,2]),]
rownames(val)=val[,2]
lst_backup=lst
#lst=as.character(val[lst,1])
lst=val[lst,1]


print(lst)
print(summary(lst))
if(min(summary(factor(lst)))<minBatches)
{
print("Not enough batches!")
return(NULL)
}



if(length(summary(factor(lst)))<2)
{
print("Not enough batches!")
return(NULL)
}


print("Make col data")
print(head(lst))
colDat=data.frame(lst)
#colDat=data.frame(colDat)
colnames(colDat)="condition"
rownames(colDat)=colnames(dat)
print(head(colDat))
#colDat[1]=factor(colDat[,1])
print(class(colDat[,1]))
if(length(keep)>0)
{
colDat[keep]=val[lst_backup,keep]
}


print(head(colDat))

print("Filter")


print("Filter")
#lst=list("meta"=colDat,"dat"=dat)
#return(lst)
group=colDat[,1]
dat=dat[rowSums(dat > minReads)>2,]
print("run")
y <- DGEList(counts=dat,group=group)
keep_genes <- filterByExpr(y)
y <- y[keep_genes,,keep.lib.sizes=FALSE]

y <- calcNormFactors(y)
print(head(colDat))
design <- model.matrix(form,colDat)
if(is.character(contrast))
{
contrast <- makeContrasts(contrasts=c(contrast),levels=design)[,1]
}
print(contrast)
print(head(design))
y <- estimateDisp(y,design,trend.method=trend.method,tagwise=tagwise,prior.df=prior.df,robust=robust)
print("Filter")
keep_edge <- filterByExpr(y)
y <- y[keep_edge, , keep.lib.sizes=FALSE]

if(getPseudo)
{
ret=list()
ret[[1]]=y
ret[[2]]=design
return(ret)
}


if(useLRT)
{
fit <- glmFit(y,design)
qlf=NULL
if(is.null(contrast))
 {
  qlf <- glmLRT(fit,coef=coef)
  }
  else{
  qlf <- glmLRT(fit,contrast=contrast)
  }
}
else
{
fit <- glmQLFit(y,design)
qlf=NULL
if(is.null(contrast))
 {
  qlf <- glmQLFTest(fit,coef=coef)
  }
  else{
  qlf <- glmQLFTest(fit,contrast=contrast)
  }
}

#return(qlf)
 res=data.frame(qlf$table)
  res["padj"]=p.adjust(res$PValue,"fdr")
res=res[order(res$PValue),]
res["Gene"]=rownames(res)

numCol=dim(res)[2]
print(head(res))
res=res[c(numCol,1:(numCol-1))]
print(head(res))
norm=-sign(res$logFC)*qnorm(res$PValue/2)
pval=fdrtool(norm,plot=F)$pval
res["pval_fdrtools"]=pval
res["padj_fdrtools"]=p.adjust(pval,"fdr")
return(res)

}








