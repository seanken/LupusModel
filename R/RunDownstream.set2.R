library(Seurat)
library(openxlsx)
source("Tools/RunPropel.R")
source("Tools/EdgeR_pb_new.R")
source("Tools/Enrich_FGSEA.R")

print("Load data")
seur=readRDS("seur.set2.micro.only.mat.RDS")
system("mkdir output_set2")

print("Test for differences in cell type composition")
system("mkdir output_set2/CellComposition")
mrk1=runPropel(seur,~0+Condition,"ConditionD-ConditionC","orig.ident","CellType")
mrk1["Test"]="ConditionD-ConditionC"
mrk2=runPropel(seur,~0+Condition,"ConditionP-ConditionC","orig.ident","CellType")
mrk2["Test"]="ConditionP-ConditionC"
mrk3=runPropel(seur,~0+Condition,"ConditionD-ConditionP","orig.ident","CellType")
mrk3["Test"]="ConditionD-ConditionP"
colnames(mrk3)[1]="PropMean.Condition1"
colnames(mrk3)[2]="PropMean.Condition2"
mrk3["CellType"]=rownames(mrk3)

colnames(mrk1)[1]="PropMean.Condition1"
colnames(mrk1)[2]="PropMean.Condition2"
mrk1["CellType"]=rownames(mrk1)


colnames(mrk2)[1]="PropMean.Condition1"
colnames(mrk2)[2]="PropMean.Condition2"
mrk2["CellType"]=rownames(mrk2)
mrk=rbind(mrk1,mrk2,mrk3)
mrk=mrk[order(mrk$P.Value),]
mrk["padj_global"]=p.adjust(mrk[,"P.Value"],"fdr")

mrk=mrk[,c(8,7,5,9,3,1,2,4)]

mrk["log2Ratio"]=log(mrk[,"PropRatio"],2)

write.table(mrk,"output_set2/CellComposition/Propel.txt",row.names=F,quote=F,sep="\t")

print(head(mrk))

print("Test for DE")
system("mkdir output_set2/DE")
refs=c("C","C","P")
kos=c("D","P","D")

celltypes=table(seur@meta.data[,"CellType"])
celltypes=names(celltypes)[as.numeric(celltypes)>500]

DE=list()
for(num in 1:3){
ko=kos[num]
ref=refs[num]
#for(ko in kos){
for(celltype in celltypes){
nam=paste(ref,ko,celltype,sep="_")
print(nam)
cells=names(seur@active.ident)[seur@meta.data[,"CellType"]==celltype & seur@meta.data[,"Condition"] %in% c(ref,ko)]
tmp=subset(seur,cells=cells)
tmp@meta.data["Condition"]=factor(tmp@meta.data[,"Condition"],c(ref,ko))
de=combDE_edgeR(tmp,kowt="Condition")
DE[[nam]]=de
print(head(de))
print(" ")
}}

for(i in names(DE)){DE[[i]]["Test"]=i}

de_comb=do.call(rbind,DE)

de_comb["padj_global"]=p.adjust(de_comb[,"PValue"],"fdr")
rownames(de_comb)=NULL
de_comb=de_comb[,grep("fdrtools",colnames(de_comb),invert=T)]
#DE=split(de_comb,Test)
DE=split(de_comb,de_comb[,"Test"])
write.xlsx(DE,"output_set2/DE/DE.results.xlsx")
saveRDS(DE,"output_set2/DE/DE.results.RDS")


print("Test for Enrichment")
system("mkdir output_set2/FGSEA")
Kegg=lapply(DE,function(mrk){mrk["score"]=-sign(mrk[,"logFC"])*log(mrk[,"PValue"]);RunEnrichment_KEGG(df=mrk,genes="Gene",Score="score")})
Go=lapply(DE,function(mrk){mrk["score"]=-sign(mrk[,"logFC"])*log(mrk[,"PValue"]);RunEnrichment_GO(df=mrk,genes="Gene",Score="score")})
Kegg=lapply(Kegg,function(x){x[8]=as.character(lapply(x[,8],function(i){paste(i,collapse="/")}));return(x)})
Go=lapply(Go,function(x){x[8]=as.character(lapply(x[,8],function(i){paste(i,collapse="/")}));return(x)})

write.xlsx(Kegg,"output_set2/FGSEA/Kegg.results.xlsx")
write.xlsx(Go,"output_set2/FGSEA/Go.results.xlsx")
saveRDS(Kegg,"output_set2/FGSEA/Kegg.results.RDS")
saveRDS(Go,"output_set2/FGSEA/Go.results.RDS")

