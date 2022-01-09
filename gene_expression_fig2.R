library(AnnotationDbi)
library(tidyverse)
library(hgu133a.db)
library(affy)
library(ggplot2)
library(ggrepel)
library(limma)
library(GEOquery)

### Example of gene expression analysis for volcano plots in figure 2 ####

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 2)
raw <- getGEO(GEO = 'GSE63067', AnnotGPL = TRUE)
raw <- raw[[1]]
status <- raw@phenoData@data$title
status <- sub(',.+', '', status)
disease <- status
status <- ifelse(status == 'Steatosis', 'NAFL',ifelse(status=='Healthy','Healthy','NASH'))
tissue <- rep('Liver', length(status))
accession <- raw@phenoData@data$geo_accession %>% as.character()

exprs_rma <- exprs(raw)
anno <- AnnotationDbi::select(hgu133a.db,
                              keys = (rownames(exprs_rma)),
                              columns = "ENTREZID",
                              keytype = "PROBEID")

### Insert this annotation to our feature data of the expressionset, after removing NAs.
anno <- anno %>% group_by(PROBEID) %>% 
  summarise(ENTREZID = paste(ENTREZID, collapse = " /// ")) %>% 
  as.data.frame(.) %>% column_to_rownames("PROBEID")


raw@featureData@data <- anno


exprs_annotated <- exprs(raw) %>% 
  aggregate(
    by = list(raw@featureData@data$ENTREZID), FUN = mean
  ) %>% 
  rename_("ENTREZID" = "Group.1") %>% 
  filter(!grepl("///", ENTREZID)) %>% 
  column_to_rownames("ENTREZID")
head(exprs_annotated)

condition = factor(status)

design_matrix <- model.matrix(~0+condition)
design_matrix

fit <- lmFit(object = exprs_annotated ,design = design_matrix, method = "ls")


contrasting <- c(paste0("conditionNAFL","-","conditionHealthy"),paste0("conditionNASH","-","conditionHealthy"))
title <- c('NAFL vs Healthy','NASH vs Healthy')
# 'Obese vs Healthy'
experiment <- "GSE63067"

for (i in 1:(length(contrasting))){
  contr <- makeContrasts(contrasts = contrasting[i], levels = colnames(coef(fit)))
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  top.table <- topTable(tmp,number = 20000)
  top.table <- top.table %>% mutate(ENTREZID=rownames(top.table))
  
  anno <- AnnotationDbi::select(hgu133a.db,
                                keys = (rownames(top.table)),
                                columns = "SYMBOL",
                                keytype = "ENTREZID")
  top.table <- left_join(top.table,anno,by="ENTREZID")
  names <- rownames(top.table)
  top.table <- top.table %>% mutate(label_pval="Not differentiated")
  top.table$label_pval[which((top.table$P.Value<0.05)& (top.table$logFC>1))] <- 'Up-regulated'
  top.table$label_pval[which((top.table$P.Value<0.05)& (top.table$logFC< -1))] <- 'Down-regulated'
  tt <- top.table %>% filter(P.Value<0.05)
  names <- tt$SYMBOL[order(abs(tt$logFC),decreasing = TRUE)][1:17]
  top.table <- top.table %>% mutate(Label=NA)
  top.table$Label[which(top.table$SYMBOL %in% names)] <- top.table$SYMBOL[which(top.table$SYMBOL %in% names)]
  top.table$label_pval <- factor(top.table$label_pval, levels = c('Up-regulated','Down-regulated','Not differentiated'))
  p1 <- ggplot(top.table,aes(x=logFC,y=-log10(P.Value))) +geom_point(aes(color = label_pval),size=1.2) +xlim(c(-4,4)) + ylim(c(0,4.75))+
    xlab('Log2(Fold Change)') + ylab('-log10(P-value)') +
    #scale_color_manual(values=c("#628FCA","#E66252","#BFBFBF"),name = "Mode of Regulation")+
    scale_color_manual(values=c("#628FCA","#E66252","#BFBFBF"))+
    theme_bw() +  
    geom_hline(yintercept=-log10(0.05), linetype="dashed",color = "black", size=0.5) +
    geom_vline(xintercept=-1, linetype="dashed",color = "#E66252", size=1) +
    geom_vline(xintercept=1, linetype="dashed",color = "#628FCA", size=1)+
    annotate("text",x=-1.8,y=1.37,label="Cut-off P-Value=0.05",size=4.5)+
    #annotate("curve", curvature = 0.2,x=-2.7, xend = -2.5, y= 1.2 , yend= -log10(0.05), arrow = arrow(length = unit(2, "mm"),type = "closed"))+
    annotate("text",x=-2,y=0.2,label="Log2(FC) = -1",size=4.5)+
    annotate("curve",x=-2.7, xend = -1, curvature=0,y= 0 , yend= 0, arrow = arrow(length = unit(2, "mm"),type = "closed"),color="#E66252")+
    annotate("text",x=2,y=0.2,label="Log2(FC) = 1",size=4.5)+
    annotate("curve",x=2.8, xend = 1, curvature=0,y= 0 , yend= 0, arrow = arrow(length = unit(2, "mm"),type = "closed"),color="#628FCA")+
    theme(text = element_text(size=13),legend.position = "none",plot.title = element_text(hjust = 0.5)) +ggtitle(title[i])+
    geom_label_repel(aes(label = Label),
                     box.padding   = 0.5, 
                     point.padding = 0.5,
                     segment.color = 'grey50')
  print(p1)
  png(paste0('supplementary_',experiment,"_",title[i],'.png'),width=9,height=12,units = "in",res=600)
  print(p1)
  dev.off()
  setEPS()
  postscript(paste0('supplementary_',experiment,"_",title[i],'.eps'))
  print(p1)
  dev.off()
}
