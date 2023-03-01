library(speckle)
library(SingleCellExperiment)
#library(CellBench)
library(limma)
library(ggplot2)
library(scater)
library(patchwork)
library(edgeR)

##example form: ~0+condition+batch
##example contrasts: conditionko-conditionwt
runPropel<-function(seur,form,contrasts,samples,clusters,transform="asin")
{
print("Get Prop")
props <- getTransformedProps(seur@meta.data[,clusters],seur@meta.data[,samples], transform=transform)
print("Get metadata")
tokeep=trimws(strsplit(as.character(form)[[2]],"+",fixed=T)[[1]])
tokeep=intersect(tokeep,colnames(seur@meta.data))
print(tokeep)
meta=seur@meta.data[,c(samples,tokeep)]
print(head(meta))
tab<-dplyr::distinct(seur@meta.data[,c(samples,tokeep)])
rownames(tab)=tab[,samples]
tab=tab[colnames(props$TransformedProps),]
print("Set up")
design <- model.matrix(form,tab)
print(head(design))
print(contrasts)
mycontr <- makeContrasts(contrasts=contrasts,levels=design)
print("Run!")
out=propeller.ttest(props, design, contrasts = mycontr, robust=TRUE, trend=FALSE,sort=TRUE)
return(out)
}
