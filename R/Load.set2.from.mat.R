source("Tools/Azimuth.R")
source("Tools/RemoveDoublet.R")
library(stringr)
library(dplyr)
library(tidyr)
library(Matrix)

print("Load data!")
#lst=readRDS("CellRanger.DS1.RDS")
lst=readRDS("CellRanger.DS2.RDS")
lst=sub("/filtered_feature_bc_matrix","",lst)
#mat=Read10X(lst)
#mat=loadDrops(lst)
mat=readRDS("mat.DS2.Aggregated.RDS")
print("Make Seurat!")
minGenes=300
colnames(mat)=sub("-1","",colnames(mat))
seur<-CreateSeuratObject(mat,"Seurat",min.features=minGenes)
seur<-NormalizeData(seur,normalization.method="LogNormalize",scale.factor=1000000)
print("Get Shared Variable Genes")
seur<-FindVariableFeatures(seur)
#seur=SharedVariable(seur)
regress=c("nFeature_RNA")
seur<-ScaleData(seur,features=seur@assays$RNA@var.features,vars.to.regress=regress)

seur<-RunPCA(seur,npcs=60)
seur<-RunUMAP(seur,dims=1:15)
seur<-FindNeighbors(seur,dims=1:15)
seur<-FindClusters(seur)

#seur@meta.data["demux_type"]=dat[names(seur@active.ident),"demux_type"]
#seur@meta.data["assignment"]=dat[names(seur@active.ident),"assignment"]

#seur@meta.data["orig.ident"]=sub("^","Run",as.character(lapply(names(seur@active.ident),function(x){strsplit(x,"_")[[1]][2]})))

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

#print("Cells to Remove")
#seur@meta.data["Rem"]=F
#seur@meta.data[seur@meta.data$demux_type=="doublet","Rem"]=T
#seur@meta.data[seur@meta.data$demux_type=="unknown","Rem"]=T
#seur@meta.data[seur@meta.data$scds>1.2,"Rem"]=T

#tab<-seur@meta.data %>% group_by(seurat_clusters) %>% summarise(PercRem=mean(Rem)) %>% as.data.frame()
#rem=tab[tab$PercRem>.5,2]
#seur@meta.data[seur@meta.data$seurat_clusters %in% rem,"Rem"]=T


#print("Save!")
#saveRDS(seur,"seur.set1.with.doublets.RDS")

#print("Remove doublets")
#seur=subset(seur,Rem==F)

#print("Recluster")
#seur<-RunUMAP(seur,dims=1:20)
#seur<-FindNeighbors(seur,dims=1:20)
#seur<-FindClusters(seur)

print("Azimuth")
ref=readRDS("Azimuth.ref.RDS")
meta=RunAzimuth(seur,ref)
for(i in colnames(meta))
{seur@meta.data[i]=meta[,i]}
seur@meta.data["Condition"]=str_sub(seur@meta.data[,"orig.ident"],1,1)
#saveRDS(seur,"seur.set2.with.doublets.RDS")
seur@meta.data["Cluster"]=as.numeric(as.character(seur@meta.data[,"RNA_snn_res.0.8"]))+1

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

#seur@meta.data["decom"]=seur@active.ident %in% c(15,17,18)
#seur=subset(seur,decom==F)

#seur<-FindVariableFeatures(seur)
#seur=SharedVariable(seur)
#regress=c("nFeature_RNA")
#seur<-ScaleData(seur,features=seur@assays$RNA@var.features,vars.to.regress=regress)

#seur<-RunPCA(seur,npcs=60)
#seur<-RunUMAP(seur,dims=1:10)
#seur<-FindNeighbors(seur,dims=1:10)
#seur<-FindClusters(seur)
saveRDS(seur,"seur.set2.micro.only.mat.RDS")

