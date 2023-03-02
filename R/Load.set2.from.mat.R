source("Tools/Azimuth.R")
source("Tools/RemoveDoublet.R")
library(stringr)
library(dplyr)
library(tidyr)
library(Matrix)

print("Load data!")
mat=readRDS("mat.DS2.Aggregated.RDS")

print("Make Seurat!")
minGenes=300
colnames(mat)=sub("-1","",colnames(mat))
seur<-CreateSeuratObject(mat,"Seurat",min.features=minGenes)
seur<-NormalizeData(seur,normalization.method="LogNormalize",scale.factor=1000000)

print("Get Shared Variable Genes")
seur<-FindVariableFeatures(seur)
regress=c("nFeature_RNA")
seur<-ScaleData(seur,features=seur@assays$RNA@var.features,vars.to.regress=regress)

seur<-RunPCA(seur,npcs=60)
seur<-RunUMAP(seur,dims=1:15)
seur<-FindNeighbors(seur,dims=1:15)
seur<-FindClusters(seur)


bot=Matrix::colSums(seur@assays$RNA@counts)

print("Get Mito")
mito=grep("^mt-",rownames(seur@assays$RNA@counts))
top=Matrix::colSums(seur@assays$RNA@counts[mito,])
mn=top/bot
seur@meta.data["mito"]=mn[names(seur@active.ident)]

print("Get ribo")
mito=grep("^Rp[s,l]",rownames(seur@assays$RNA@counts))
top=Matrix::colSums(seur@assays$RNA@counts[mito,])
mn=top/bot
seur@meta.data["ribo"]=mn[names(seur@active.ident)]

print("Get doublet scores")
seur=DoubletScores(seur)


print("Azimuth")
ref=readRDS("Azimuth.ref.RDS")
meta=RunAzimuth(seur,ref)
for(i in colnames(meta))
{seur@meta.data[i]=meta[,i]}
seur@meta.data["Condition"]=str_sub(seur@meta.data[,"orig.ident"],1,1)
seur@meta.data["Cluster"]=as.numeric(as.character(seur@meta.data[,"RNA_snn_res.0.8"]))+1

print("Assign Cell Types to Clusters")
lst=rep("Homeostatic",21)

lst[21]="S100a4+"
lst[17]="Cycling"
lst[c(7,9)]="Ms4a7+"
lst[15]="Ms4a7+ low activation"
lst[13]="Homeostatic low activation"
lst[c(18,19,20)]="Non-microglia"
lst[16]="Tmem119- Myeloid"
lst[14]="Interferon-responsive"


seur@meta.data["CellType_ref"]=lst[seur@meta.data[,"Cluster"]]
lst[15]="Ms4a7+"
lst[13]="Homeostatic"
seur@meta.data["CellType"]=lst[seur@meta.data[,"Cluster"]]

saveRDS(seur,"seur.set2.with.doublets.mat.RDS")

seur=subset(seur,CellType!="Non-microglia")

saveRDS(seur,"seur.set2.micro.only.mat.RDS")

