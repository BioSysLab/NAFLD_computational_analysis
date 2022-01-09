library(tidyverse)
library(FactoMineR)
library(factoextra)
library(animation)
library(ggbiplot)


proteomics <- read.csv("C:/Users/PC/Desktop/Manuscript/Results/PCA/ProteomicsBio.csv", row.names=1)
head(proteomics)
summary(proteomics)

str(proteomics)
summary(proteomics)
colSums(is.na(proteomics))

#PCA
pca_proteomics2 <- PCA(proteomics,
                       scale.unit = FALSE,
                       graph = FALSE)
summary(pca_proteomics2)


plot.PCA(pca_proteomics2,
         choix = c("ind","var","varcor"))

ggplot(as.data.frame(pca_proteomics2$ind$coord[,1:2]), aes(x=pca_proteomics2$ind$coord[,1], y=pca_proteomics2$ind$coord[,2])) +
  geom_point(size=3) +  
  geom_text(label=rownames(proteomics), nudge_x = 0.05, nudge_y = 0.05, size=3) +
  ggtitle("Principal Component Analysis") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Principal Component 1") + 
  ylab("Principal Component 2")
ggsave("PCA.png", dpi=600, units = "cm", width = 20, height = 20)

#K-MEANS
RNGkind(sample.kind = "Rounding") #to get the set.seed numbers not changed everytime executed
kmeansTunning <- function(data, maxK) {
  withinall <- NULL
  total_k <- NULL
  for (i in 2:maxK) {
    set.seed(100)
    temp <- kmeans(data,i)$tot.withinss
    withinall <- append(withinall, temp)
    total_k <- append(total_k,i)
  }
  plot(x = total_k, y = withinall, type = "o", col="blue", 
       xlab = "Number of Cluster", 
       ylab = "Total within", 
       main="Elbow Method for Optimal K",
       frame.plot = axes)
}


# kmeansTunning
kmeansTunning(proteomics, maxK = 10)

set.seed(100)
proteomic_cluster <- kmeans(proteomics, centers = 4)
fviz_cluster(proteomic_cluster, data = proteomics, ellipse = FALSE)

#Plot Kmeans on PCA

ggplot(as.data.frame(pca_proteomics2$ind$coord[,1:2]), 
       aes(x=pca_proteomics2$ind$coord[,1], y=pca_proteomics2$ind$coord[,2], color=as.factor(proteomic_cluster$cluster))) +
  geom_point(size=3) +  
  geom_text(label=rownames(proteomics), nudge_x = 0.1, nudge_y = 0.1, size=3) +
  ggtitle("Principal Component Analysis & k-means clustering") + 
  guides(col=guide_legend("Cluster")) +
  theme(plot.title = element_text(hjust = 0.5), legend.title = element_text(hjust = 0.5)) +
  xlab("Principal Component 1") + 
  ylab("Principal Component 2")
ggsave("PCA_kmeans.png", dpi=600, units = "cm", width = 20, height = 20)
ggsave("PCA_kmeans.eps", dpi=600, units = "cm", width = 20, height = 20)



