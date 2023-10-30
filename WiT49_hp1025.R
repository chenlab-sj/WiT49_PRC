
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(readxl)
library(devtools)
library(pheatmap)
library(ggplotify)
library(dendsort)
library(edgeR)

args = commandArgs(trailingOnly=TRUE)
Usage <-"Rscript WiT49_hp1025.R genes.xls Wit49_count.csv NC_count.csv NC_samplemeta.csv ouput.path" 
if (length(args)==0) {
  stop(Usage, call.=FALSE)
} else {
  genes.xls <- args[1]
  Wit49_count.csv <- args[2]
  NC_count.csv <- args[3]
  NC_samplemeta.csv <- args[4]
  res.folder <-args[5]
  
}

genes <- read_excel(genes.xls)
WiT49_counts <- read.csv(Wit49_count.csv,sep = "\t")
NC_counts.all <- read.csv(NC_count.csv, sep = "\t")
NCsample_index <- read.csv(NC_samplemeta.csv, sep = "\t")



genes <-as.data.frame(genes) %>%
  filter(`Marker of` != "Stroma")

WiT49_meta<-data.frame(
  id = c(colnames(WiT49_counts)[2:7]),
  Type = c(rep("WiT49-Parental",3),rep("WiT49-PRC",3)),
  Blastema.level=rep('NA',6)
)

NCsample_index%>%
  mutate(Blastema.level = 
           case_when(Percent.Blastema >= 80 ~ 'high_blastema',
                     Percent.Blastema <= 20 ~ 'low_blastema',
                     TRUE ~ "Null"))%>%
  filter(`Blastema.level`!="Null")%>%
  filter(!is.na(id))%>% 
  separate(id2, into = c("Type", "id2.label"), sep = "-")%>%
  filter(Type=='PT' | Type=='KT')->NCsample_index2

NC_meta <- NCsample_index2[,c("id","Type","Blastema.level")]
NC_counts <- NC_counts.all[,c('GeneName',NCsample_index2$id)]

rownames(WiT49_counts)<-WiT49_counts$gene_name
rownames(NC_counts)<-NC_counts$GeneName
WiT49_counts <- WiT49_counts[,-1]
NC_counts <- NC_counts[,-1]
WiT49_counts.cpm <- cpm(WiT49_counts, log=TRUE,normalized.lib.sizes=TRUE)
NC_counts.cpm<- cpm(NC_counts, log=TRUE,normalized.lib.sizes=TRUE)
combined.cpm <- merge(WiT49_counts.cpm,NC_counts.cpm,by = 'row.names', all = TRUE)

## meta
combined.meta <-rbind(WiT49_meta,NC_meta)
#combined.meta <- combined.meta %>% arrange(factor(id, levels=sample_name))
combined.meta <- combined.meta %>%
  mutate(Type = 
           case_when(Type == 'PT' ~ 'Patient sample',
                     Type == 'KT'~ 'PDX',
                     TRUE ~ combined.meta$Type),
         Blastema.level =
           case_when(Blastema.level == 'low_blastema'~ 'Low-level blastema',
                     Blastema.level == 'high_blastema'~ 'High-level blastema',
                     TRUE ~ Blastema.level))
print("work")
##################################################################
## combine WiT and NC data 
genes %>% merge(combined.cpm, by.x = 'Gene list', by.y ='Row.names',all.x= TRUE )%>%
  drop_na()->combined_genesig.cpm

rownames(combined_genesig.cpm)<-combined_genesig.cpm$`Gene list`
combined_genesig.cpm <- combined_genesig.cpm[,-c(1,2)]

gene_name<-rownames(combined_genesig.cpm)
gene.meta<-genes[genes$`Gene list` %in% gene_name,]
gene_meta <- gene.meta %>% arrange(factor(`Gene list`, levels=gene_name))
row.names(gene.meta)<-gene.meta$`Gene list`
gene.meta<-gene.meta %>% select(-`Gene list`)


##################################################################

## 2, limma
# dds <- DESeq(dds)
# normalized_counts <- counts(dds, normalized = TRUE)
# vsd <- vst(dds, blind = FALSE)
# mat <- assay(vsd)

WiT49.id <- WiT49_meta$id
PT.id <-NC_meta$id[NC_meta$Type=='PT']
KT.id <-NC_meta$id[NC_meta$Type=='KT']
combined.batch<-case_when(
  colnames(combined.cpm) %in% WiT49.id ~ 1,
  colnames(combined.cpm) %in% PT.id ~2,
  colnames(combined.cpm)%in% KT.id ~3,
  TRUE ~0
)
combined.batch<-combined.batch[-1]
rownames(combined.cpm)<-combined.cpm$Row.names
combined.cpm <- combined.cpm[,-1]


mat <- limma::removeBatchEffect(combined.cpm, batch = combined.batch)
as.data.frame(mat)%>%
  mutate(`Row.names`=rownames(mat))-> combined.cpm.nobatch

genes %>% merge(combined.cpm.nobatch, by.x = 'Gene list', by.y ='Row.names',all.x= TRUE )%>%
  drop_na()->combined_genesig.cpm

rownames(combined_genesig.cpm)<-combined_genesig.cpm$`Gene list`
combined_genesig.cpm <- combined_genesig.cpm[,-c(1,2)]

colnames(combined.meta)<-c("id","Type","Blastema level")
## 1. ipsc and PT samples

sample_name.sub1 <- combined.meta$id[combined.meta$Type!='PDX']
combined_genesig.cpm1 <- combined_genesig.cpm[,sample_name.sub1]
combined.meta1<-combined.meta[combined.meta$Type!='PDX',]
row.names(combined.meta1) <- combined.meta1$id
combined.meta1 <-  combined.meta1[, -1]
combined_genesig.std1 <- t(scale(t(combined_genesig.cpm1),center = TRUE, scale = TRUE))
#cols <- colorRampPalette(rev(brewer.pal(10, "RdYlBu")))(256)
cols = colorRampPalette( c("green", "black", "red"), space="rgb")(32)

anno.col1<-list(Type=c(`WiT49-Parental`="#d95f02", `WiT49-PRC`="#1b9e77", `Patient sample`="#7570b3"),
               `Blastema level`=c(`High-level blastema`="#e41a1c",`Low-level blastema`="#377eb8",`NA`="gray"),
               `Marker of`=c(Blastemal = "#fb8072",Epithelial="#b3de69"))

marker_levels<-c("Blastemal","Stroma","Epithelia")
gene_meta.order <- gene_meta[order(factor(gene_meta$`Marker of`,levels = marker_levels)),]
combined_genesig.std1<-combined_genesig.std1[gene_meta.order$`Gene list`,]

callback = function(hc, ...){dendsort(hc,isReverse=TRUE)}
pheatmap(combined_genesig.std1,
         #distfun = function(x) as.dist((1-cor(t(x)))/2),
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "complete",
         clustering_callback = callback,
         color=cols,
         cluster_rows = F,
         annotation_col=combined.meta1,
         annotation_colors = anno.col1,
         annotation_row = gene.meta,
         border_color=NA,
         annotation_legend=TRUE,
         fontsize = 5,
         cellheight=5, cellwidth = 5)->p1

# group = kmeans(t(mat), centers = 3)$cluster
# Heatmap(mat, name = "mat", cluster_columns = cluster_within_group(mat, group))

## 2. ipsc and KT samples

sample_name.sub2 <- combined.meta$id[combined.meta$Type!='Patient sample']
combined_genesig.cpm2 <- combined_genesig.cpm[,sample_name.sub2]
combined.meta2<-combined.meta[combined.meta$Type!='Patient sample',]
row.names(combined.meta2) <- combined.meta2$id
combined.meta2 <-  combined.meta2[, -1]
#combined_counts.std2 <-scale(combined_counts.mtx2, center = TRUE, scale = TRUE)
combined_genesig.std2 <- t(scale(t(combined_genesig.cpm2),center = TRUE, scale = TRUE))
combined_genesig.std2<-combined_genesig.std2[gene_meta.order$`Gene list`,]

anno.col2<-list(Type=c(`WiT49-Parental`="#d95f02", `WiT49-PRC`="#1b9e77",  PDX="#e7298a"),
                `Blastema level`=c(`High-level blastema`="#e41a1c",`Low-level blastema`="#377eb8",`NA`="gray"),
                `Marker of`=c(Blastemal = "#fb8072",Epithelial="#b3de69"))
pheatmap(combined_genesig.std2,
         #distfun = function(x) as.dist((1-cor(t(x)))/2),
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "complete",
         color=cols,
         cluster_rows = F,
         annotation_col=combined.meta2,
         annotation_row = gene.meta,
         annotation_colors = anno.col2,
         border_color=NA,
         fontsize = 5,
         cellheight=5, cellwidth = 5)->p2

if (!file.exists(res.folder)) {
  dir.create(res.folder)
}

png(paste(res.folder, "Wilm_1_genesig_1026.png", sep=""), res =300, width=2200, height=1500)
as.ggplot(p1)
dev.off()

png(paste(res.folder, "Wilm_2_genesig_1026.png", sep=""), res =300, width=2200, height=1500)
as.ggplot(p2)
dev.off()

