## Introduction
#The aim of the present study is to do Pathway Analysis (PA) over NonAlcoholic Fat Liver Disease (NAFLD) data. The procedure we will follow is:
# Download GSE data from GEO.
# Use the LIMMA algorithm to analyze the gene expression data.
# Use the Piano algorithm to do pathway analysis.

rm(list = ls())
wd <- file.path(getwd(), '..')
library(GEOquery)
library(affy)
library(tidyverse)
library(gridExtra)
library(limma)
library(viridis)
library(piano)
library(ggfortify)
library(sva)
wd <- file.path(getwd(), '..')

pianoPathWays <- function(pathFile){
  library(piano)
  library(GSEABase)
  
  #Load the KnowledgeBase
  genesets = getGmt(con =
                      "/Users/Nafsika/Dropbox/Diplwmatiki/piano/c2.cp.v6.0.symbols.gmt")
  genesets = unlist(genesets)
  
  gene_to_term <- data.frame(NA,NA)
  names(gene_to_term) <- c("gene","term")
  for (geneset in genesets){
    temp <- geneIds(geneset)
    temp2 <- setName(geneset)
    temp3 <- as.data.frame(cbind(temp,rep(temp2,length(temp))))
    names(temp3) <- c("gene","term")
    gene_to_term <- rbind(gene_to_term,temp3)
  }
  
  geneSet <- loadGSC(gene_to_term)
  
  #Load the LimmaResult
  gene_list <- read.csv2(pathFile, sep = ';')
  
  myFC <- gene_list$logFC
  names(myFC) <- toupper(gene_list$X)
  
  myPval <- gene_list$adj.P.Val
  names(myPval) <- toupper(gene_list$X)
  
  myTval <- gene_list$t
  names(myTval) <- toupper(gene_list$X)
  
  ###Run the GSA
  gsaRes1 <- runGSA(geneLevelStats = myPval, directions = myFC, gsc =
                      geneSet, geneSetStat = 'fisher')
  gsaRes2 <- runGSA(geneLevelStats = myPval, directions = myFC, gsc =
                      geneSet, geneSetStat = 'stouffer')
  gsaRes3 <- runGSA(geneLevelStats = myPval, directions = myFC, gsc =
                      geneSet, geneSetStat = 'reporter')
  gsaRes4 <- runGSA(geneLevelStats = myPval, directions = myFC, gsc =
                      geneSet, geneSetStat = 'tailStrength')
  gsaRes5 <- runGSA(geneLevelStats = myTval, gsc = geneSet, geneSetStat
                    = 'sum')
  gsaRes6 <- runGSA(geneLevelStats = myTval, gsc = geneSet, geneSetStat
                    = 'maxmean')
  gsaRes7 <- runGSA(geneLevelStats = myTval, gsc = geneSet, geneSetStat
                    = 'page')
  gsaRes8 <- runGSA(geneLevelStats = myFC, gsc = geneSet, geneSetStat =
                      'mean')
  gsaRes9 <- runGSA(geneLevelStats = myFC, gsc = geneSet, geneSetStat =
                      'median')
  
  resList <-
    list(gsaRes1,gsaRes2,gsaRes3,gsaRes4,gsaRes5,gsaRes6,gsaRes7,gsaRes8,gsaRes9)
  names(resList) <-
    c('fisher','stoufer','reporter','tailStrength','sum','maxmean','page','mean','media
n')
  
  return(resList)
}

pianoAnalysis <- function(wd, experiments, i){
  
  resList <- readRDS(file.path(wd, 'Results', 'AnalysisByTissue',
                               'PathwayAnalysis',
                               paste0('resList_', experiments, '.rds')))
  
  resList <- resList[i]
  names <- names(resList)
  
  summary_results <- unlist(resList, recursive = FALSE)
  names(summary_results) <- sub('.', '_', names(summary_results), fixed =
                                  TRUE)
  
  summary_df <- lapply(names(summary_results), function(x)
    GSAsummaryTable(summary_results[[x]], save = TRUE,
                    file = paste0(wd,'/Results/AnalysisByTissue/','/PathwayAnalysis/',
                                  '/IndividualSummary/', i, '_', x,'_', experiments, '.txt')))
  
  names(summary_df) <- names(summary_results)
  
  df_all <- bind_rows(summary_df, .id = "uid") %>%
    separate(col = uid, into = c('Contrast','Method'), sep = '_')
  
  ConsensusPValsAndRanks <- function(df_list,class_type,direction_type){
    consensus_scores <- consensusScores(summary_results, class =
                                          class_type, direction = direction_type,
                                        n = 10e6, adjusted = TRUE, method = 'mean', plot =
                                          FALSE)
    
    pVals_pathways <- consensus_scores$pMat %>%
      -log10(.) %>%
      apply(., 2, function(x) replace(x, is.infinite(x),
                                      mean(as.numeric(x[which(!is.infinite(x))])))) %>%
      rowMeans(TRUE) %>%
      as.data.frame
    colnames(pVals_pathways) <- paste0(class_type,'_',direction_type)
    pVals_pathways$Pathway <- rownames(pVals_pathways)
    pVals_pathways <- pVals_pathways[,c(2,1)]
    
    rank_pathways <- consensus_scores$rankMat
    
    write.csv2(pVals_pathways, file =
                 paste0(wd,'/Results/AnalysisByTissue/','/PathwayAnalysis/',
                        
                        '/PValsPathways/',class_type,'_',direction_type,'PVals_Pathways_',
                        experiments,'_',i, '.csv'))
    
    write.csv2(pVals_pathways, file =
                 paste0(wd,'/Results/AnalysisByTissue/','/PathwayAnalysis/',
                        
                        '/RanksPathways/',class_type,'_',direction_type,'Ranks_Pathways_',
                        experiments,'_',i, '.csv'))
    
    return(list(pVals = pVals_pathways, ranks = rank_pathways))
  }
  
  distinct_up <- ConsensusPValsAndRanks(summary_results,'distinct','up')
  distinct_down <-
    ConsensusPValsAndRanks(summary_results,'distinct','down')
  mixed_up <- ConsensusPValsAndRanks(summary_results,'mixed','up')
  mixed_down <- ConsensusPValsAndRanks(summary_results,'mixed','down')
  nond <- ConsensusPValsAndRanks(summary_results,'non','')
  
  pVals_all <- Reduce(function(x, y) merge(x, y, all=TRUE),
                      list(distinct_up$pVals,mixed_up$pVals,
                           nond$pVals,mixed_down$pVals,distinct_down$pVals))
  
  write.csv2(pVals_all, file =
               paste0(wd,'/Results/AnalysisByTissue/','/PathwayAnalysis/',
                      '/','_','_PVals_Pathways_', experiments,'_', i , '.csv'))
  
  range01 <- function(x){(x - min(x,na.rm=TRUE))/(max(x,na.rm=TRUE)-
                                                    min(x,na.rm=TRUE))}
  
  pVals_all_norm <- pVals_all
  pVals_all_norm$distinct_up <- range01(pVals_all$distinct_up)
  pVals_all_norm$mixed_up <- range01(pVals_all$mixed_up)
  pVals_all_norm$non_ <- range01(pVals_all$non_)
  pVals_all_norm$mixed_down <- range01(pVals_all$mixed_down)
  pVals_all_norm$distinct_down <- range01(pVals_all$distinct_down)
  
  pVals_df <- bind_rows(arrange(pVals_all_norm,-distinct_up)[1:12,],
                        arrange(pVals_all_norm,-mixed_up)[1:12,],
                        arrange(pVals_all_norm,-non_)[1:12,],
                        arrange(pVals_all_norm,-mixed_down)[1:12,],
                        arrange(pVals_all_norm,-distinct_down)[1:12,]) %>%
    melt()
  pVals_df <- pVals_df %>% mutate(Pathway = factor(Pathway, levels =
                                                     unique(pVals_df$Pathway)))
  
  hm1 <- ggplot(pVals_df, aes(x=variable, y=Pathway, fill=value)) +
    geom_tile(color="white", size=0.1) +
    scale_fill_viridis(name="St(-log10(PVal))", begin = 0, end = 1) +
    coord_equal() +
    labs(x=NULL, y=NULL, title="Affected Pathways NASHVsHealthy") +
    scale_y_discrete(label=function(x) abbreviate(x, minlength=40)) +
    theme_tufte(base_family="Helvetica") +
    theme(axis.ticks=element_blank(), axis.text=element_text(size=6),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(paste0(wd,'/Results/AnalysisByTissue/','/PathwayAnalysis/',
                '/','_','_HeatMap_Pathways_', experiments,'_', i, '.jpeg'), width = 6,
         height = 9, units = 'in')
  
  ranks_all <- list(distinct_up = distinct_up$ranks,
                    mixed_up = mixed_up$ranks,
                    nond = nond$ranks,
                    mixed_down = mixed_down$ranks,
                    distinct_down = distinct_down$ranks)
  
  return(list(pVals = pVals_all, ranks = ranks_all, HeatMap = hm1))
}

plotGenes <- function(wd, gene_to_term,limma_data,pathway){
  
  genes <- filter(gene_to_term,term == pathway)$gene
  limma_data_path <- filter(limma_data,TF %in% genes)
  if (experiments == "GSE63067") {
    limma_data_path$logP <- -log10(limma_data_path$P.Value)
  } else {
    limma_data_path$logP <- -log10(limma_data_path$adj.P.Val)
  }
  limma_data_path$Activity <-
    (limma_data_path$logFC/abs(limma_data_path$logFC)) *
    sqrt(limma_data_path$logFC^2 + limma_data_path$logP^2)
  
  vp1 <- ggplot(limma_data_path, aes(logFC, logP)) +
    geom_point(aes(size = abs(Activity), colour = logFC)) +
    scale_colour_gradient2(low = 'firebrick', mid = 'gray90', high = 'forestgreen',
                           limits = c(min(limma_data_path$logFC),
                                      max(limma_data_path$logFC))) +
    geom_text(data = subset(limma_data_path, logP > -log10(0.05) &
                              abs(logFC) > 0.5),
              hjust = 0, nudge_x = 0.05, aes(label = TF), size = 4) +
    labs(title = paste0(pathway), x = 'log2(FC)', y = '-log10(Pvalue)') +
    theme(plot.title = element_text(hjust = 0.7))
  
  ggsave(file.path(wd, 'Results', 'AnalysisByTissue', 'PathwayAnalysis',
                   'VolcanoPlots',
                   paste0(pathway, '_ByContrast_', experiments,'_',i, '.pdf')), width =
           12, height = 9, units = 'in')
  
  ldTF <- limma_data_path$TF
  ldFC <- limma_data_path$logFC
  limma_data_p <- data.frame(ldTF, ldFC)
  colnames(limma_data_p) <- c('TF','logFC')
  
  return(list('GeneValues' = limma_data_p, 'VolcanoPlotContrast' = vp1))
} 

## Download GSE data
#Two different data sets (GSEs) were identified: GSE63067 and GSE89632.
# GSE63067 ----------------------------------------------------------------
exp <- 'GSE63067'
filePaths <- getGEOSuppFiles(exp, baseDir = paste0(wd,'/Data/CELFiles'))
#to download the cel files

if(file.exists(file.path('Data','CELFiles',exp,paste0(exp,'_RAW')))) {
  untar(rownames(filePaths)[1], exdir = file.path("Data", "CELFiles", exp,
                                                  paste0(exp,'_RAW')))
} else {
  dir.create(file.path('Data', 'CELFiles', exp,paste0(exp, '_RAW')))
  untar(rownames(filePaths)[1], exdir = file.path("Data", "CELFiles", exp,
                                                  paste0(exp,'_RAW')))
}
cel_data <- ReadAffy(celfile.path =
                       file.path('Data','CELFiles',exp,paste0(exp,'_RAW')))
eset <- rma(cel_data)
mat2 <- exprs(eset) #Expression Data
gse2 <- getGEO(GEO = 'GSE63067')
gse2 <- gse2[[1]]
status <- gse2@phenoData@data$title
status <- sub(',.+', '', status)
disease <- status
status <- ifelse(status == 'Healthy', 'Healthy', 'Disease')
tissue <- rep('Liver', length(status))
accession <- gse2@phenoData@data$geo_accession %>% as.character()
GeneSymbol2 <- list(gse2@featureData@data$`Gene Symbol`)
#Corresponding Gene
mat2 <- as.data.frame(mat2)
mat2 <- aggregate(mat2, by = GeneSymbol2, FUN = mean)
colnames(mat2) <- gsub('_.+','',colnames(mat2))
colnames(mat2)[1] <- 'GeneSymbol'
mat2 <- filter(mat2, !grepl('///', mat2$GeneSymbol))
write.csv2(mat, file.path(wd, 'Data', 'GSEData', 'ExpressionData', paste0(exp,
                                                                          '_data.txt')),
           row.names = FALSE)
platform <- as.list(levels(phenoData(gse)$platform))
platform <- platform[[1]]
experiment <- rep(exp,ncol(mat)-1)
platform <- rep(platform,ncol(mat)-1)

target <- data.frame(accession, experiment, status, disease, tissue, platform)
write.csv2(target, file.path(wd, 'Data', 'GSEData', 'TargetData', paste0(exp,
                                                                         '_target.txt')),
           row.names = FALSE)
# GSE89632 ----------------------------------------------------------------
exp <- 'GSE89632'
gse <- getGEO(GEO = exp)
gse <- gse[[1]]
mat <- exprs(gse)
pgse <- pData(gse)
status <- gse@phenoData@data$characteristics_ch1.1
status <- sub('diagnosis: ', '', status)
disease <- status
status <- ifelse(status == 'HC', 'Healthy', 'Disease')
tissue <- rep('Liver', length(status))
accession <- gse@phenoData@data$geo_accession %>% as.character()
GeneSymbol <- list(gse@featureData@data$`Symbol`)
mat <- as.data.frame(mat)
mat <- aggregate(mat, by = GeneSymbol, FUN = mean)
mat <- as.data.frame(mat) %>% mutate(PROBEID = rownames(.))
colnames(mat) <- gsub('_.+','',colnames(mat))
colnames(mat)[1] <- 'GeneSymbol'
mat <- filter(mat, !grepl('///', mat$GeneSymbol))
write.csv2(mat, file.path(wd, 'Data', 'GSEData', 'ExpressionData', paste0(exp,
                                                                          '_data.txt')),
           row.names = FALSE)
platform <- as.list(levels(phenoData(gse)$platform))
platform <- platform[[1]]
experiment <- rep(exp,ncol(mat)-2)
platform <- rep(platform,ncol(mat)-2)
target <- data.frame(accession, experiment, status, disease, tissue, platform)
write.csv2(target, file.path(wd, 'Data', 'GSEData', 'TargetData', paste0(exp,
                                                                         '_target.txt')),
           row.names = FALSE)

## Load and Filter Samples
wd <- file.path(getwd(), '..')
tissue <- 'Liver'
expr <- ifelse(tissue == 'Liver', 'Liver')
#experiments <- "GSE63067"
experiments <- "GSE89632"
target <- read_delim(file = file.path(wd, 'Data', 'GSEData','TargetData',
                                      paste0(experiments, '_target.txt')), delim = ';')
Gtarget <- target[grep(expr,target$tissue),]
if (experiments == "GSE63067") {
  RefDiseases <- c('Steatosis', 'Healthy', 'NASH') #changes
  Gtarget <- Gtarget %>% mutate(disease = str_replace_all(disease, "Nonalcoholic steatohepatitis" ,"NASH"))
} else {
  RefDiseases <- c('SS', 'HC', 'NASH') #changes
}
Gtarget <- Gtarget[which(Gtarget$disease %in% RefDiseases),] %>%
  droplevels()
write.table(table(Gtarget[c("tissue", "disease")]),
            file.path(wd, 'Data', 'GSEData', paste0('Table_TissueVsDisease_',
                                                    experiments,'_', tissue, '.txt')), sep=',')
table(Gtarget[c('tissue','disease')])
write.csv2(Gtarget, file.path(wd, 'Data', 'GSEData', paste0(tissue,'_',
                                                            experiments, '_target.txt')), sep = ';')
### Patient Expression Data
experiments <- levels(as.factor(as.character(Gtarget$experiment)))
files <- lapply(experiments, function(x) list.files(path = file.path(wd, 'Data',
                                                                     'GSEData', 'ExpressionData'), pattern = x)) %>% unlist #changes
gse <- read.csv2(file = file.path(wd, 'Data', 'GSEData','ExpressionData',
                                  paste0(experiments, '_data.txt')), dec = ',', sep = ";")
rownames(gse) <- gse$GeneSymbol
gse <- gse[, which(colnames(gse) %in%
                     levels(as.factor(as.character(Gtarget$accession))))]
df <- apply(gse, 1, function(x) replace(x, is.na(x),
                                        mean(as.numeric(x[which(!is.na(x))]))))
df <- data.frame(t(df))
write.csv2(df, file.path(wd, 'Data', 'GSEData', paste0('Liver_data_',
                                                       experiments, '.txt')), sep = ';')
## Gene-level statistics
#Two different subtypes of the disease were included in this study (plus
#                                                                   control). For each subtype, a contrast was built comparing the samples
#belonging to the specific subtype with the healthy samples. For each contrast,
#Gene-Level Statistics (GLSs) were calculated, i.e., fold changes (fc), p-values
#(PVal), and t statistics (t-stat) or Gene Set Analysis (GSA).
#LIMMA algorithm was used to calculate such GLSs.
#LIMMA first fits a linear model for each gene (using the different subtypes of
#                                              the disease as variables). Then, given the linear models, computes the
#estimated coefficients and standard errors for a given set of contrasts. Finally,
#computes GLSs by empirical Bayes.

f1 <- factor(Gtarget$disease)
design <- model.matrix(~0 + f1)
colnames(design) <- levels(f1)
fit <- lmFit(df, design)
# set up a contrast matrix to compare between states
if (experiments == "GSE63067") {
  contrast.matrix <- makeContrasts(NASHvHealthy = NASH - Healthy,
                                   StvHealthy = Steatosis - Healthy, levels = design) #changed
} else {
  contrast.matrix <- makeContrasts(NASHvHealthy = NASH - HC, StvHealthy
                                   = SS - HC, levels = design) #changed
}


contrast.matrix
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
gene_list <- lapply(1:2, function(x) topTable(fit2, coef = x, number = 1000000,
                                              sort.by = "logFC")) #changed
names(gene_list) <- colnames(contrast.matrix)
for(i in names(gene_list)){
  print(i)
  print(gene_list[[i]][1:10,])
  write.csv2(gene_list[[i]], file = file.path(wd, 'Results', 'AnalysisByTissue',
                                              'LIMMA',
                                              paste0('limma_', experiments, '_', i, '.txt')))
}
i <- names(gene_list)[2]
contrast <- i
if (experiments == "GSE63067") {
  tmp_genes <- gene_list[[contrast]] %>%
    mutate(GeneID = rownames(.)) %>%
    mutate(logP = -log10(P.Value)) %>%
    mutate(Status = ifelse(abs(logFC) < 1 | logP < -log10(0.05),
                           'Undifferentiated', 'Differentiated'))
} else {
  tmp_genes <- gene_list[[contrast]] %>%
    mutate(GeneID = rownames(.)) %>%
    mutate(logP = -log10(adj.P.Val)) %>%
    mutate(Status = ifelse(abs(logFC) < 1 | logP < -log10(0.05),
                           'Undifferentiated', 'Differentiated'))
}
most_expressed <- tmp_genes %>%
  arrange(-abs(logFC)) %>%
  arrange(Status) %>%
  head(n = 15)
ggplot(tmp_genes, aes(logFC, logP, colour = Status)) +
  geom_point() +
  geom_label_repel(data = most_expressed, aes(logFC, logP, label = GeneID,
                                              hjust = 1)) + ggtitle(contrast) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(file.path(wd, 'Results', 'AnalysisByTissue', 'PathwayAnalysis',
                 paste0('Volcano_plot_GLS_', experiments,'_',i, '.jpeg')), width = 6,
       height = 5, units = 'in')
## Gene Set Analysis (GSA)
#Diseases are rarely originated due to alterations in individual genes/proteins.
#What its true for TFs is not necessarily true for proteins.
#We are not only interested in big changes in individual TFs, but also in 'small'
#changes in functionally correlated TFs.
#Gain functional insights about the disease and its mechanisms.
#GSA helps understanding the underlying mechanism of the disease.

library(piano)
library(GSEABase)
genesets = getGmt(con =
                    file.path("/Users/Nafsika/Dropbox/Diplwmatiki/piano/c2.cp.v6.0.symbols.gmt"))
genesets = unlist(genesets)
gene_to_term <- data.frame(NA, NA)
names(gene_to_term) <- c("gene", "term")
for (geneset in genesets){
  temp <- geneIds(geneset)
  temp2 <- setName(geneset)
  temp3 <- as.data.frame(cbind(temp, rep(temp2, length(temp))))
  names(temp3) <- c("gene", "term")
  gene_to_term <- rbind(gene_to_term, temp3)
}
head(gene_to_term, n = 15)

#We ran [piano]
#(https://bioconductor.org/packages/release/bioc/html/piano.html), a package
#from R that performs several GSA methods and then consensuates the
#different results.
#The justification for this is that there is not a GSA method clearly better than
#the rest.

#Piano needs a template where the individual genes are mapped to sets. The
#template selected for this study was the [MSigDB-CuratedPathwaysCanonicalPathways]
#(http://software.broadinstitute.org/gsea/msigdb/collections.jsp#C2) knowledge
#  base. This knowledge base is represented by a bipartite graph mapping
#  individual genes to functional sets (terms).
#  With the GLSs and the knowledge base, we already have the two inputs
#  necessary to run any piano GSA method. For each contrast, 9 GSA methods
#  were run:
    
#  * Fisher (PVal)
#  * Stouffer (PVal)
#  * Reporter (PVal)
#  * Tail Strength (PVal)
#  * Sum (TVal)
#  * MaxMean (TVal)
#  * PAGE (TVal)
#  * Mean (FC)
#  * Median (FC)

source(file.path(wd, 'PIANO_mysigDB_cp.R'))
resList <- lapply(colnames(contrast.matrix),function(x) pianoPathWays(file.path(wd, 'Results',
                                                      'AnalysisByTissue','LIMMA', paste0('limma_', experiments,'_', x,'.txt'))))
names(resList) <- colnames(contrast.matrix)
saveRDS(resList,file = file.path(wd, 'Results', 'AnalysisByTissue','PathwayAnalysis',
                                 paste0('resList_', experiments, '.rds')))

i <- names(gene_list)[2]
source(file.path(wd,'PIANO_analysis_results.R'))
gsa_result <- pianoAnalysis(wd, experiments, i)

#For each GSA run (contrast-method pair) and pathway, 5 p-value classes are
#obtained:

#* Distinct-up: Up-regulation (up and down are cancelled out).
#* Mixed-up: Up-regulation (up and down are not cancelled out, a pathway can
#                           be mixed-up and mixed-down regulated at the same time).
#* Non-Directional: Differentially regulated (regardless the sign).
#* Mixed-down: Down-regulation (up and down are not cancelled out).
#* Distinct-down: Down-regulation (up and down are cancelled out).
#For each run, 5 ranked list are generated (one for each p-value).
#Eventually, for each pathway a group of 5 p-values is obtained (one of each
#                                                                class). For the different p-values the -log10 is computed.
#Once the most affected pathways have been identified, we can go one step
#backwards and take a look to the GLSs of all the individual genes participating
#in such pathways.

limma_data <- read.csv2(file.path(wd, 'Results', 'AnalysisByTissue', 'LIMMA',
                                  paste0('limma_', experiments,'_', i, '.txt')), sep = ';',
                        row.names = 1) %>%
  mutate(TF = rownames(.))
pathway <- arrange(gsa_result$pVals, -distinct_up)[1,]$Pathway
source(file.path(wd,'PIANO_analysis_results.R'))
individual_genes <- plotGenes(wd, gene_to_term, limma_data, pathway)