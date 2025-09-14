library(rms);library(survival);library(survminer)
library(ggplot2);library(ggsci);library(cowplot);library(plotly)
library(tidyverse);library(PredictABEL);library(ggpubr)
library(openxlsx);library(tableone);library(reshape2);library(caret)
library(pROC);library(timeROC);library(dplyr);library(plyr)
library(ComplexHeatmap); library(gridExtra);library(ggrepel)
library(gplots);library(data.table);library(pheatmap)
library(DESeq2);library(clusterProfiler);library(org.Hs.eg.db)
library(enrichplot);library(TCGAbiolinks);library(stringr)
library(SummarizedExperiment);library(limma);library(IOBR)
library(estimate);library(RColorBrewer);library(STRINGdb)
library(igraph);library(ggraph);library(ggsignif)
library(circlize)

setwd("path/to/data")

train_set<-read.xlsx("/train_set_data.xlsx")
inval_set<-read.xlsx("/inval_set_data.xlsx")
test_set<-read.xlsx("/test_set_data.xlsx")
TCGA_TCIAset<-read.xlsx("TCGA_TCIA_set_data.xlsx")


# PCA analysis
data  <- train_set[, c("ITH_proba", "GTR_proba")]
pca_result <- prcomp(data, center = TRUE, scale. = TRUE)
summary(pca_result)
print(pca_result$rotation)
train_set$IntegratedIndex <- -pca_result$x[,1]

newdata <- inval_set[, c("ITH_proba", "GTR_proba")]
inval_set$IntegratedIndex <- -predict(pca_result, newdata = newdata)[, 1]

newdata <- test_set[, c("ITH_proba", "GTR_proba")]
test_set$IntegratedIndex <- -predict(pca_result, newdata = newdata)[, 1]

newdata <- TCGA_TCIAset[, c("ITH_proba", "GTR_proba")]
TCGA_TCIAset$IntegratedIndex <- -predict(pca_result, newdata = newdata)[, 1]

# IntegratedIndex risk value
summary(train_set$IntegratedIndex)
cutoffs_value_Inte<-median(train_set$IntegratedIndex)
train_set$Inte_label<-ifelse(train_set$IntegratedIndex<=cutoffs_value_Inte,0,1)
inval_set$Inte_label<-ifelse(inval_set$IntegratedIndex<=cutoffs_value_Inte,0,1)
test_set$Inte_label<-ifelse(test_set$IntegratedIndex<=cutoffs_value_Inte,0,1)
TCGA_TCIAset$Inte_label<-ifelse(TCGA_TCIAset$IntegratedIndex<=cutoffs_value_Inte,0,1)

#------baseline table-----
combined_df <- bind_rows(
  train_set %>% mutate(Group = "1"),
  inval_set %>% mutate(Group = "2"),
  test_set %>% mutate(Group = "3"),
  TCGA_TCIAset %>% mutate(Group = "4"))
table1 <- CreateTableOne (vars =c("Var1","Var2","Var3"),
                          factorVars = c("Var1","Var2","Var3"), 
                          data=combined_df, strata ='Group',addOverall = TRUE )
table11<-print(table1,nonnormal = c("Var1"), quote = FALSE, noSpaces = TRUE, printToggle = FALSE, showAllLevels = TRUE, smd = TRUE);table11
write.csv(table11, file = "results/baseline_table.csv")


#------------Logist/Cox regression------------
fit_log <- glm(ORR ~ X1 + X2 + X3, data = train_set, family = binomial())
summary(fit_log)

fit_cox <- coxph(Surv(OS, dead) ~ X1 + X2 + X3, data = train_set)
summary(fit_cox)

#--------Survival analysis----------
#replaced per different datasets
fit <- survfit(Surv(OS,dead)~Inte_label,data =train_set);fit  
survdiff(Surv(OS,dead)~Inte_label,data =train_set,rho=0)
ggsurvplot(fit,data =train_set, 
                pval = TRUE, pval.size = 4.5, pval.coord = c(2.5,0.1), risk.table = TRUE, 
                xlim = c(0,40), break.x.by = 5, ylim = c(0,1), break.y.by = 0.1,  
                xlab = "Months", ylab = "Overall Survival",
                legend = c(0.75,0.9),  legend.title = "", theme = theme(legend.background = element_rect(fill = "transparent", colour = NA)), 
                legend.labs = c( "Low risk","High risk"))
 


#------------ROC curve------------
model_train_ITH <- roc(train_set$ORR,train_set$ITH_proba)
best_coords_train_ITH <- coords(model_train_ITH, x = "best", best.method = "youden",
                                 ret = c("threshold","sensitivity","specificity","accuracy","ppv","npv"),
                                 transpose = FALSE);best_coords_train_ITH
best_threshold_train_ITH <- as.numeric(best_coords_train_ITH["threshold"])
  
model_train_GTR <- roc(train_set$ORR,train_set$GTR_proba)
best_coords_train_GTR <- coords(model_train_GTR, x = "best", best.method = "youden",
                                 ret = c("threshold","sensitivity","specificity","accuracy","ppv","npv"),
                                 transpose = FALSE);best_coords_train_GTR
best_threshold_train_GTR <- as.numeric(best_coords_train_GTR["threshold"])
  
model_train_Inte <- roc(train_set$ORR,train_set$IntegratedIndex)
best_coords_train_Inte <- coords(model_train_Inte, x = "best", best.method = "youden",
                                 ret = c("threshold","sensitivity","specificity","accuracy","ppv","npv"),
                                 transpose = FALSE);best_coords_train_Inte
best_threshold_train_Inte <- as.numeric(best_coords_train_Inte["threshold"])
  
f1 <- glm(ORR~x1+x2+x3+x4+x5, data = train_set, family = binomial())
train_set$pre<-predict(f1,type='response')
model_train_clin <- roc(train_set$ORR,train_set$pre)

plot(model_train_ITH,col="#086FA1", legacy.axes=TRUE,reuse.auc=TRUE,axes=TRUE,
       xlim=c(1, 0), ylim=c(0, 1),xlab="1 - Specificity", ylab="Sensitivity",
       lty=par("lty")) 
lines(model_train_GTR, type="l",lty=1,lwd=2,col="#009E73FF")
lines(model_train_Inte, type="l",lty=1,lwd=2,col="#BC3C29FF")
lines(model_train_clin, type="l",lty=1,lwd=2,col="#FF7F00")
title(main="", col.main="black", font.main=2)


#------------Calibration Curve------------
par(mfrow=c(1, 1))
par(oma=c(2, 2, 2, 2))
groups<-10
rangeaxis<-c(0,1)
cOutcome<-3 

rc1<-glm(ORR ~ IntegratedIndex, family=binomial(), data=train_set)
rc2<-glm(ORR ~ IntegratedIndex, family=binomial(), data=inval_set)
rc3<-glm(ORR ~ IntegratedIndex, family=binomial(), data=test_set)

prisk1<-rc1$fitted.values
prisk2<-rc2$fitted.values
prisk3<-rc3$fitted.values

a<-plotCalibration(data=train_set,cOutcome = cOutcome, predRisk=prisk1,groups=groups,rangeaxis = rangeaxis)
b<-plotCalibration(data=inval_set,cOutcome = cOutcome, predRisk=prisk2,groups=groups,rangeaxis = rangeaxis)
c<-plotCalibration(data=test_set,cOutcome = cOutcome, predRisk=prisk3,groups=groups,rangeaxis = rangeaxis)


plot_calibration_curve <- function(data, color, pch, lty) {
  num_points <- 4  
  n <- nrow(data)  
  indices <- seq(from = 1, to = n, length.out = num_points) 
  indices1 <- round(indices+n) 
  indices2 <- round(indices1+n) 
  lines(data[indices1], data[indices2], type="b", lwd=2, col=color, pch=pch, lty=lty, xlim=c(0,1), ylim=c(0,1))}


plot(NA, NA, xlim=c(0, 1), ylim=c(0, 1), xlab="Predicted probability", ylab="Actual probability", type="n")
plot_calibration_curve(a[[1]], "#BC3C29FF", 16, 1)
plot_calibration_curve(b[[1]], "#0072B5FF", 16, 1)
plot_calibration_curve(c[[1]], "#E18727FF", 16, 1)
abline(0,1,lty=2,lwd=2,col="#708090")

legend("bottomright", bty="n",c("Ideal","Training set", "Internal validation set","Test set")  
       ,lty=c(2,1,1,1),col=c("#708090","#BC3C29FF", "#0072B5FF","#E18727FF"))  
  

#-------DCA curve-------

library(ggDCA)


m1 <- glm(ORR ~ITH_proba,data = train_set)
m2 <- glm(ORR ~GTR_proba,data = train_set)
m3 <- glm(ORR ~IntegratedIndex,data = train_set)
m4 <- glm(ORR~x1+x2+x3+x4,data = train_set, family = binomial())

dd <- datadist(train_set)
options(datadist = "dd")
d_train <- dca(m1,m2,m3,m4,model.names = c('ITH','GTR','GTR-ITH','Clinical'))
ggplot(d_train,xlim(0,1),linetype = FALSE, smooth = TRUE,
       color = c("#20854EFF", "#0072B5FF","#BC3C29FF","#E18727FF","#708090","#000000"))+
  scale_x_continuous(breaks = seq(0, 1, by = 0.1))


#-----Model comparison-Delong-NRI-IDI-----
library(PredictABEL)
groups<-10
rangeaxis<-c(0,1)
cOutcome<-3 

reclassification(data=train_set,cOutcome=cOutcome,predrisk1=train_set$GTR_proba,predrisk2=train_set$ITH_proba,cutoff=c(0, 0.2, 0.4, 0.6, 0.8, 1))
reclassification(data=inval_set,cOutcome=cOutcome,predrisk1=inval_set$GTR_proba,predrisk2=inval_set$ITH_proba,cutoff=c(0, 0.2, 0.4, 0.6, 0.8, 1))
reclassification(data=test_set,cOutcome=cOutcome,predrisk1=test_set$GTR_proba,predrisk2=test_set$ITH_proba,cutoff=c(0, 0.2, 0.4, 0.6, 0.8, 1))

#-----Delong test-----
roc1 <- roc(train_set$ORR, train_set$ITH_proba)
roc2 <- roc(train_set$ORR, train_set$GTR_proba)
roc.test(roc1, roc2, method = "delong")

roc1 <- roc(inval_set$ORR, inval_set$ITH_proba)
roc2 <- roc(inval_set$ORR, inval_set$GTR_proba)
roc.test(roc1, roc2, method = "delong")

roc1 <- roc(test_set$ORR, test_set$ITH_proba)
roc2 <- roc(test_set$ORR, test_set$GTR_proba)
roc.test(roc1, roc2, method = "delong")





#----Bioinformatic analysis----
options(timeout = 3000)

# Download data
# query <- GDCquery(project = "TCGA-LIHC",
#                   data.category = "Transcriptome Profiling",  
#                   data.type = "Gene Expression Quantification",  
#                   workflow.type = "STAR - Counts")  
# GDCdownload(query, method = "api", files.per.chunk = 100) 
# GDCprepare(query, save = TRUE, save.filename = "TCGA-LIHC_mRNA.Rdata")
# data = GDCprepare(query) 
# exp = assay(data) 

#----Load data----
load(file="path/to/TCGA-LIHC_data.Rdata") 
expr_counts_mrna <- assay(data,"unstranded") 
expr_tpm_mrna <- assay(data,"tpm_unstrand") 
expr_fpkm_mrna <- assay(data,"fpkm_unstrand") 

sample_types <- sapply(strsplit(colnames(expr_counts_mrna), "-"), function(x) x[4])
table(sample_types)  
filtered_expr_counts_mrna <- expr_counts_mrna[, sample_types %in% c("01A")] 
filtered_colnames <- colnames(filtered_expr_counts_mrna)
simplified_colnames <- sapply(strsplit(filtered_colnames, "-"), function(x) paste(x[1:3], collapse = "-"))
colnames(filtered_expr_counts_mrna) <- simplified_colnames

bio_label=TCGA_TCIAset[,c("case_submitter_id","ITH_label","GTR_label","Inte_label")] 

#-----------Expression matrix-----------
filtered_expr_counts_mrna <- filtered_expr_counts_mrna[, colnames(filtered_expr_counts_mrna) %in% bio_label$case_submitter_id]

RNA_symbol <-rowData(data)$gene_name
filtered_expr_counts_mrna<-cbind(as.data.frame(RNA_symbol),filtered_expr_counts_mrna)
dim(filtered_expr_counts_mrna)

qc= as.matrix(filtered_expr_counts_mrna)
rownames(qc)=qc[,1] 
exp=qc[,2:ncol(qc)]  
dimnames=list(rownames(exp),colnames(exp))
data_qc=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
dim(data_qc) 

data_qc_test = as.data.frame(data_qc[rowSums(data_qc)>0,]) 
nrow(data_qc_test)
data_qc_test=data_qc_test[rowMeans(data_qc_test)>1,] 
nrow(data_qc_test) 

data_qc_test$genes<-rownames(data_qc_test)
data_qc_test=aggregate(.~genes,mean,data=data_qc_test)
nrow(data_qc_test) 

original_colnames <- colnames(data_qc_test[,-1])
original_rownames <- data_qc_test[,1]

data_qc_test_integer <- as.matrix(data_qc_test[,-1])

rownames(data_qc_test_integer) <- original_rownames
colnames(data_qc_test_integer) <- original_colnames
exprSet <- data_qc_test_integer
write.csv(exprSet, file=".path/to/exprSet.csv")

#----label input-----
ref_df <- data.frame(case_submitter_id = colnames(filtered_expr_counts_mrna)[-1])
bio_label_adjust <- left_join(ref_df, bio_label, by = "case_submitter_id")
write.csv(bio_label_adjust, file=".path/to/bio_label_adjust.csv")

Inte_group_list=factor(ifelse(bio_label_adjust$Inte_label==0,'GTR-ITH high-risk','GTR-ITH low-risk'),
                       levels = c("GTR-ITH high-risk", "GTR-ITH low-risk"))
table(Inte_group_list)

#----ABRS score-----
vsd <- varianceStabilizingTransformation(exprSet)

selected_genes <- c("CXCR2P1", "ICOS", "TIMD4", "CTLA4", "PAX5", "KLRC3", "FCRL3", "AIM2", "GBP5", "CCL4")
vsd_selected_genes <- as.data.frame(t(vsd[rownames(vsd) %in% selected_genes, ]))

vsd_selected_genes$ABRS_score  <-rowMeans(vsd_selected_genes)
original_rownames <- rownames(vsd_selected_genes)
vsd_selected_genes$case_submitter_id <- original_rownames

vsd_selected_genes <- left_join(vsd_selected_genes, TCGA_TCIAset_bioinf[,c("case_submitter_id","ITH_proba","ITH_label","GTR_proba","GTR_label","IntegratedIndex","Inte_label")] , by = "case_submitter_id")
rownames(vsd_selected_genes) <- original_rownames
vsd_selected_genes$ABRS_score_HL<-ifelse(vsd_selected_genes$ABRS_score<=median(vsd_selected_genes$ABRS_score),0,1)

#Intel
correlation_result <- cor.test(vsd_selected_genes$IntegratedIndex, vsd_selected_genes$ABRS_score, method = "pearson");correlation_result
t.test(ABRS_score ~ Inte_label, data = vsd_selected_genes)

p3<-ggplot(vsd_selected_genes, aes(x = as.factor(Inte_label), y = ABRS_score, 
                                   fill = as.factor(Inte_label), color = as.factor(Inte_label))) +
  geom_boxplot(width = 0.5) +   
  labs(title = "", x = "GTR-ITH index", y = "") +
  theme_minimal() +  
  scale_fill_manual(name = NULL,labels = c("High risk", "Low risk"),values=c("#B22222","#446BB3"))+ 
  scale_color_manual(name = NULL,labels = c("High risk", "Low risk"),values=c("black","black"))+
  geom_signif(comparisons = list(c("0", "1")),   
              annotations = "** P = 0.002",             
              y_position = 10,                          
              tip_length = 0.05,                        
              color = "black",textsize = 4)+       
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        axis.title.x  = element_text(size = 16, face = "bold"),
        axis.title.y  = element_text(size = 16, face = "bold"),
        axis.text.y   = element_text(size = 14, face = "bold")
  );p3


#----Differentially expressed genes----
cut_off_pvalue = 0.05 
cut_off_logFC = 1       

colData_Inte <- data.frame(row.names = colnames(exprSet),Inte_group_list= Inte_group_list)
dds_Inte <- DESeqDataSetFromMatrix(countData = exprSet,colData = colData_Inte,design = ~ Inte_group_list)
dds_Inte <- DESeq(dds_Inte)
resultsNames(dds_Inte)
res_Inte <-  results(dds_Inte, contrast=c("Inte_group_list","GTR-ITH high-risk","GTR-ITH low-risk"))
resOrdered_Inte <- res_Inte[order(res_Inte$padj),]
resOrdered_Inte=as.data.frame(resOrdered_Inte)
resOrdered_Inte$change = ifelse(resOrdered_Inte$padj < cut_off_pvalue & abs(resOrdered_Inte$log2FoldChange) >= cut_off_logFC, 
                                ifelse(resOrdered_Inte$log2FoldChange>= cut_off_logFC ,'Up','Down'),'Stable') 

resOrdered_Inte <- dplyr::filter(resOrdered_Inte,  !is.na(change))
write.csv(resOrdered_Inte, file="/path/to/DGE order_Inte.csv")


#-----volcano plot------
ggplot(resOrdered_Inte, 
        aes(x = log2FoldChange, y = -log10(padj), colour=change)) +
    geom_point(alpha=0.5, size=3.5)+  
    scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+ 
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) + 
    geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +
    labs(x="log2 (FoldChange)", y="-log10 (adjusted P value)")  

table(resOrdered_Inte$change)


#-----Heatmap------
gene_up_Inte = resOrdered_Inte[resOrdered_Inte$change == 'Up',] 
gene_down_Inte = resOrdered_Inte[resOrdered_Inte$change == 'Down',] 

deg_opt <- resOrdered_Inte %>%
  rownames_to_column(var = "rownames") %>%  
  filter(change != "Stable") %>%         
  arrange(log2FoldChange) %>%             
  column_to_rownames(var = "rownames")    

gene_Inte=rownames(deg_opt)
mrna_expr_count1=log(exprSet+1)
exp_heatmap <- mrna_expr_count1[gene_Inte,]
annotation_col <- data.frame(Group = Inte_group_list)
rownames(annotation_col) <- colnames(exp_heatmap) 
annotation_col <- annotation_col %>%
  rownames_to_column(var = "rownames") %>%
  arrange(Group) %>%                     
  column_to_rownames(var = "rownames") 

annotation_col2 <- annotation_col %>%
  rownames_to_column("case_submitter_id") %>%
  left_join(TCGA_TCIAset[,c("case_submitter_id","sex","dead")], by="case_submitter_id")%>%
  column_to_rownames("case_submitter_id") %>%
  dplyr::select(Dead=dead, Sex= sex, Group)   

my_GTRot<-pheatmap(exp_heatmap[,rownames(annotation_col2)], 
                  show_colnames = F, show_rownames =F,  
                  scale = "row", 
                  cluster_cols = F,  
                  annotation_col = annotation_col2,
                  breaks = seq(-1, 1, length.out = 100),  
                  color = colorRampPalette(c("#0099CC", "white", "#CC3333"))(100),
                  fontsize = 12, fontsize_row = 10, fontsize_col = 10,
                  main = "",
                  gaps_col = cumsum(table(annotation_col2$Group)),
                  treeheight_row = 50, treeheight_col = 50, 
                  annotation_colors=list(
                    Group=c("GTR-ITH low-risk"="#446BB3","GTR-ITH high-risk"="#B22222"),
                    Sex  = c("Male"="#66CC00", "Female"="#FFCC33"),
                    Dead = c("Dead"="#3399CC", "Alive"="#99CC33")
                  ));my_GTRot



#-----GO------
gene_up_Inte_entrez <- as.character(na.omit(bitr(rownames(gene_up_Inte), 
                                                 fromType="SYMBOL", toType="ENTREZID", 
                                                 OrgDb="org.Hs.eg.db",drop = T)[,2])) 
gene_down_Inte_entrez <- as.character(na.omit(bitr(rownames(gene_down_Inte), 
                                                   fromType="SYMBOL", toType="ENTREZID", 
                                                   OrgDb="org.Hs.eg.db",drop = T)[,2])) 
gene_diff_Inte_entrez <- unique(c(gene_up_Inte_entrez ,gene_down_Inte_entrez ))

gene_diff_Inte = c(rownames(gene_up_Inte),rownames(gene_down_Inte)) 
df_id_Inte<-bitr(gene_diff_Inte,fromType = "SYMBOL", toType = "ENTREZID",
                 OrgDb = "org.Hs.eg.db",drop=TRUE)  

GO_Inte <- enrichGO(gene = df_id_Inte$ENTREZID,keyType = "ENTREZID",
                    OrgDb= org.Hs.eg.db,  ont = "ALL", pAdjustMethod = "fdr",
                    pvalueCutoff  = 0.05, qvalueCutoff  = 0.05,
                    minGSSize = 10, maxGSSize = 500)

my_GTRot<-dotplot(GO_Inte,showCategory=20,font.size = 12,label_format=50,title = "")
my_GTRot


#-----KEGG------
KEGG_Inte<-enrichKEGG(gene=df_id_Inte$ENTREZID,
                      keyType = "kegg", organism  = "hsa", 
                      pvalueCutoff = 0.05, qvalueCutoff = 0.2)

my_GTRot<-dotplot(KEGG_Inte,title = "");my_GTRot


#----PPI network----
string_db <- STRINGdb$new( version="12", species=9606,  score_threshold=400, 
                           input_directory="")
data_mapped <- df_id_Inte %>% string_db$map(my_data_frame_id_col_names = "ENTREZID", 
                                            removeUnmappedRows = TRUE)
string_db$plot_network(data_mapped$STRING_id)

hit<-data_mapped$STRING_id
info <- string_db$get_interactions(hit)

links <- info %>%
  mutate(from = data_mapped[match(from, data_mapped$STRING_id), "SYMBOL"]) %>% 
  mutate(to = data_mapped[match(to, data_mapped$STRING_id), "SYMBOL"]) %>%  
  dplyr::select(from, to , last_col()) %>% 
  dplyr::rename(weight = combined_score)

nodes <- links %>% { data.frame(gene = c(.$from, .$to)) } %>% distinct()

net <- igraph::graph_from_data_frame(d=links,vertices=nodes,directed = F)

igraph::V(net)$deg <- igraph::degree(net)
igraph::V(net)$size <- igraph::degree(net)/5 #
igraph::E(net)$width <- igraph::E(net)$weight/10

ggraph(net,layout = "kk")+
  geom_edge_fan(aes(edge_width=width), color = "lightblue", show.legend = F)+
  geom_node_point(aes(size=size), color="orange", alpha=0.7)+
  geom_node_text(aes(filter=deg>5, label=name), size = 5, repel = T)+
  scale_edge_width(range = c(0.2,1))+
  scale_size_continuous(range = c(1,10) )+
  guides(size=F)+
  theme_graph()

ggraph(net,layout = "stress")+ 
  geom_edge_fan(aes(edge_width=width), color = "lightblue", show.legend = F)+
  geom_node_point(aes(size=size), color="orange", alpha=0.7)+
  geom_node_text(aes(filter=deg>5, label=name), size = 5, repel = T)+
  scale_edge_width(range = c(0.2,1))+
  scale_size_continuous(range = c(1,10) )+
  guides(size=F)+
  theme_graph()

ggraph(net,layout = "linear", circular = TRUE)+
  geom_edge_fan(aes(edge_width=width), color = "lightblue", show.legend = F)+
  geom_node_point(aes(size=size), color="orange", alpha=0.7)+
  geom_node_text(aes(filter=deg>5, label=name), size = 5, repel = F)+
  scale_edge_width(range = c(0.2,1))+
  scale_size_continuous(range = c(1,10) )+
  guides(size=F)+
  theme_graph()


links_2 <- links %>%  add_count(from, name = "from_c") %>%  
  add_count(to, name = "to_c") %>%
  filter(!(from_c == 1 & to_c == 1)) %>% 
  dplyr::select(1, 2, 3)  


nodes_2 <- links_2 %>% { data.frame(gene = c(.$from, .$to)) } %>% distinct()

net_2 <- igraph::graph_from_data_frame(d=links_2,vertices=nodes_2,directed = F)

igraph::V(net_2)$deg <- igraph::degree(net_2)
igraph::V(net_2)$size <- igraph::degree(net_2)/5
igraph::E(net_2)$width <- igraph::E(net_2)$weight/10


#------Immune infiltration cibersort--------
exprSet_cibersort<-exprSet

cibersort_ori <- deconvo_tme(eset = exprSet_cibersort, method = "cibersort", 
                             arrays = TRUE,tumor = TRUE,
                             perm = 100)

write_tsv(cibersort_ori, ".path/to/immune/cibersort.txt")

cibersort<-cibersort_ori
cibersort$ID <- gsub("\\.", "-", sub("^([^.]+\\.[^.]+\\.[^.]+)\\..*$", "\\1", cibersort$ID))

cell_bar_GTRot(input = cibersort, title = "Cell Fraction", legend.position = "top", 
              pattern ="CIBERSORT", features = colnames(cibersort)[2:23], 
              coord_filp = T, palette = 2) 

# #epic
epic_res <- deconvo_tme(eset = exprSet_cibersort, method = "epic", arrays = TRUE)
epic_GTRot <- cell_bar_GTRot(input = epic_res, features = colnames(epic_res)[2:9],
                           title = "EPIC Cell Fraction")

# #quantiseq
quantiseq_res <- deconvo_tme(eset = exprSet_cibersort, tumor = F, arrays = TRUE,
                             scale_mrna = TRUE, method = "quantiseq")
quantiseq_GTRot <- cell_bar_GTRot(input = quantiseq_res, features = colnames(quantiseq_res)[2:12],
                                title = "Quantiseq Cell Fraction")

ips_res <- deconvo_tme(eset = exprSet_cibersort, method = "ips", plot= F);ips_res

##------pheatmap------
data_total <- cbind(cibersort[,-c(24:26)],epic_res[,-1],quantiseq_res[,-1])
data_total <- as.data.frame(t(data_total))
colnames(data_total) <- data_total[1,]
data_total <- dplyr::slice(data_total,-1)

data_total <- data_total %>% mutate_all(as.numeric)


stand_fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {
  outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata > halfwidth] = halfwidth
    outdata[outdata < (-halfwidth)] = -halfwidth
  }
  return(outdata)
}
data_total <- stand_fun(data_total,halfwidth = 2)
data_total <- data_total[!apply(is.na(data_total), 1, all), ]

## Inte
bio_label_Inte=bio_label[,c("case_submitter_id","Inte_label")]
rownames(bio_label_Inte)<-bio_label_Inte$case_submitter_id
bio_label_Inte <- subset(bio_label_Inte, select = -1)

sorted_index <- order(bio_label_Inte$Inte_label, decreasing = FALSE)
bio_label_Inte  <- bio_label_Inte[sorted_index, , drop = FALSE]
sorted_names <- rownames(bio_label_Inte)
data_total_Inte <- data_total[, sorted_names]
bio_label_Inte$Inte_label<-as.factor(bio_label_Inte$Inte_label)
colnames(bio_label_Inte)[colnames(bio_label_Inte) == "Inte_label"]  <- "Group"
bio_label_Inte$Group <- factor(ifelse(bio_label_Inte$Group==0,'GTR-ITH high-risk','GTR-ITH low-risk'),
                               levels = c('GTR-ITH high-risk','GTR-ITH low-risk'))


annCol <- data.frame(Group = bio_label_Inte,     
                     row.names =rownames(bio_label_Inte),stringsAsFactors = FALSE)

methods <- sub('.*_', '', rownames(data_total_Inte)) 

annRow <- data.frame(Methods = factor(methods, levels = unique(methods)),
                     row.names = rownames(data_total_Inte),
                     stringsAsFactors = FALSE)

breaksList = seq(-2, 2, by = 0.1)  
colors <- colorRampPalette(c("#336699", "white", "red"))(length(breaksList))


my_GTRot<-pheatmap(data_total_Inte,  
                  annotation_col = annCol,  annotation_row = annRow,
                  color = colors, breaks = breaksList,
                  cluster_rows = FALSE,  cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = FALSE,
                  gaps_col = cumsum(table(annCol$Group)),  
                  gaps_row = cumsum(table(annRow$Methods)),
                  annotation_colors=list(Group=c("GTR-ITH low-risk"="#446BB3","GTR-ITH high-risk"="#B22222")),
                  fontsize_row = 6, fontsize_col = 6,
                  annotation_names_row = FALSE);my_GTRot

##------boxplot------
cibersort2<- dplyr::select(cibersort_ori, -"P-value_CIBERSORT", -"Correlation_CIBERSORT", -"RMSE_CIBERSORT")
cibersort2 <- cibersort2 %>%  left_join(bio_label, by = c("ID" = "case_submitter_id")) 

colnames(cibersort2) <- gsub("_CIBERSORT", "", colnames(cibersort2))
cibersort2 <- cibersort2 %>% column_to_rownames(var = "ID")
cibersort2 <- cibersort2[, c((ncol(cibersort2)-2):ncol(cibersort2), 1:(ncol(cibersort2)-3))]
cibersort2 <-cibersort2[, !names(cibersort2) %in% c("T_cells_CD4_naive", "T_cells_gamma_delta")]

cibersort_Inte <- cibersort2[ , -c(1:2)]
x=colnames(cibersort_Inte)[1]
colnames(cibersort_Inte)[1]="Type"

data<- reshape2::melt(cibersort_Inte,id.vars=c("Type"))
colnames(data)=c("Type","Gene","Expression")
data$Type <- factor(ifelse(data$Type==0,'GTR-ITH high-risk','GTR-ITH low-risk'),
                    levels = c('GTR-ITH high-risk','GTR-ITH low-risk'))

data$Gene <- gsub("_", " ", data$Gene)

p=ggboxplot(data, x="Gene", y="Expression", fill = "Type", 
            ylab="Infiltration levels",
            xlab="",
            legend.title=x,
            palette = c("cornflowerblue","orange1"),
            add = "none")
p=p+rotate_x_text(60)
p3=p+stat_compare_means(aes(group=Type),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif", size = 8 )+
  theme(legend.text    = element_text(size = 14, face="bold"),  
        legend.key.size= unit(1.5, "lines"),                 
        legend.title = element_blank(),
        axis.title.y  = element_text(size = 14, face = "bold"),
        axis.text.y   = element_text(size = 12, face = "bold"));p3


#----Immune score-----

file_dir="path/to/exprSet.txt"  
data <- read.csv(".path/to/exprSet.csv", header = TRUE, sep = ",") 
write.table(data, file_dir, sep = "\t", row.names =FALSE, col.names =  TRUE,quote = FALSE)

filterCommonGenes(input.f=file_dir, output.f=".path/to/immunegene.gct",  
                  id="GeneSymbol") 

estimateScore(".path/to/immune/gene.gct",".path/to/immune/estimat_score.gct",platform="affymetr")

score_table=read.table(".path/to/immune/estimat_score.gct",skip=2,header=1) 
rownames(score_table)=score_table[,1] 
score=t(score_table)  
colnames(score)=score[1,] 
score_matrix=score[-1,]
write.table(score,".path/to/immune/tumor_purityscore.txt",sep="\t")

purityscore<-as.data.frame(score_matrix[-1,])
purityscore <- rownames_to_column(purityscore, var = "ID")

purityscore$ID <- gsub("\\.", "-", sub("^([^.]+\\.[^.]+\\.[^.]+)\\..*$", "\\1", purityscore$ID))

bio_label<-read.csv('./TCGA_TCIAset_bioinf.csv')


#Inte
bio_label_Inte=bio_label[,c("case_submitter_id","Inte_label")] 

purityscore_Inte <- purityscore %>%  left_join(bio_label_Inte, by = c("ID" = "case_submitter_id")) 
purityscore_Inte[,c("StromalScore","ImmuneScore","ESTIMATEScore","TumorPurity")]<-
  lapply(purityscore_Inte[,c("StromalScore","ImmuneScore","ESTIMATEScore","TumorPurity")], as.numeric)


purityscore_Inte$Inte_label<-as.factor(ifelse(purityscore_Inte$Inte_label==0,"GTR-ITH high-risk","GTR-ITH low-risk"))
purityscore_Inte <- column_to_rownames(purityscore_Inte, var = "ID")


scores2 <- pivot_longer(data = purityscore_Inte,cols = colnames(purityscore_Inte)[1:4],names_to = "purityscore_Inte",values_to = 'value')
scores2$Inte_label <- as.character(scores2$Inte_label)
p3<-ggplot(scores2, aes(x = Inte_label, y = value, fill = Inte_label)) +
  geom_boxplot(position = position_dodge(0.8)) +
  theme_bw(base_size = 16) + 
  theme(panel.grid=element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
        
  ) +
  facet_wrap(. ~ purityscore_Inte,scales = 'free',ncol = 4)+
  scale_fill_nejm() +
  stat_compare_means(comparisons = combn(unique(scores2$Inte_label), 2, 
                                         simplify =FALSE), method = 'wilcox.test')+
  theme(axis.text.x = element_blank(),
        legend.text    = element_text(size = 14, face="bold"), 
        legend.key.size= unit(1.5, "lines"),                   
        legend.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x  = element_blank(),
        axis.title.y  = element_text(size = 14, face = "bold"),
        axis.text.y   = element_text(size = 12, face = "bold"));p3  
