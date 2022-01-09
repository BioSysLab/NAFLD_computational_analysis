############ Packages ############
library(reshape2)
library(AICcmodavg)
library(car)
library(multcomp)
library (tidyverse)
library(datasets)
library(gapminder)
library(purrr)
library(broom)
library(cowplot)
library(FactoMineR)
library(factoextra)
library(devtools)
library(rstatix)
library(ggpubr)
library(ggsignif)

############ Import files ############

pho_raw <- read.csv("D:/BIOL NTUA/NAFLD/Manuscript/Results/phosphos.csv")
head(pho_raw)

cyt_raw <- read.csv("D:/BIOL NTUA/NAFLD/Manuscript/Results/cytos.csv")
head(cyt_raw)

se <- function(x) sqrt(var(x)/length(x))                                    ##SEM function 


############ Phosphos ############

pho <- pho_raw %>% gather (key="Protein", value ="FD", SMAD3:RS6) %>% drop_na()
head(pho)
summary(pho)

#################   Graphs per Drug 

############# Treatment tables

dmem <- pho %>% filter(Sample == 'DMEM')
etoh <- pho %>% filter(Sample == 'ETOH')
et_dm <- pho %>% filter(Sample == 'ETOH_DMSO')
ffa <- pho %>% filter(Sample == 'FFA')
ffa_dm <- pho %>% filter(Sample == 'FFA_DMSO')
rsv_ffa <- pho %>% filter(Sample == 'RSV_FFA')
sir_ffa <- pho %>% filter(Sample == 'SIR_FFA')
dif_ffa <- pho %>% filter(Sample == 'DIF_FFA')
pra_ffa <- pho %>% filter(Sample == 'PRA_FFA')
fen_ffa <- pho %>% filter(Sample == 'FEN_FFA')
gal_ffa <- pho %>% filter(Sample == 'GAL_FFA')
rsv <- pho %>% filter(Sample == 'RSV')
sir <- pho %>% filter(Sample == 'SIR')
dif <- pho %>% filter(Sample == 'DIF')
pra <- pho %>% filter(Sample == 'PRA')
fen <- pho %>% filter(Sample == 'FEN')
gal <- pho %>% filter(Sample == 'GAL')

####### Sirolimus

sir_all <- full_join(dmem, etoh) %>%  full_join(et_dm) %>% 
  full_join(ffa) %>%  full_join(ffa_dm) %>%
  full_join(sir) %>%  full_join(sir_ffa) 
sir_all[sir_all==""] <- NA
transform(sir_all, as.numeric(FD))

head(sir_all)
summary(sir_all)

sir_all <- sir_all %>% group_by(Sample,Protein) %>%
    mutate(med=median(FD),
            iqr=IQR(FD)) %>%
  filter((FD<=med+1*iqr) & (FD>=med-1*iqr)) %>%
  ungroup()

sir_stats <- sir_all %>% group_by(Sample, Protein) %>%
  summarise(mean=median(FD),
            sd=sd(FD),
            cv=sd/mean*100,
            sem=se(FD)) %>%
  ungroup()

max_sir <- max(sir_all$FD)

sir_final <- left_join(sir_all,sir_stats)

sir_final$Sample <- factor(sir_final$Sample, levels=c("DMEM", "ETOH","ETOH_DMSO", 
                                                    "FFA","FFA_DMSO", "SIR_FFA", "SIR"))

bar_sir <- sir_final %>% group_by(Protein) %>%
  ggplot(aes(x=Sample, y=mean, fill=Sample)) + geom_col(position = position_dodge()) + ylim(0,2) +
  geom_point(aes(y=FD)) +
  theme_minimal() +
  scale_fill_manual(values=c("#5C9EFF", "#5C9EFF", "#5C9EFF", "#F23B4F", "#F23B4F", "#5E4DFF", "#B4D194")) +
  geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.8) +
  facet_wrap(~ Protein, scales = "free") +
  labs(title="Sirolimus", x ="Protein", y = "Fold Change") +
  theme(plot.title = element_text(size=20, hjust=0.5),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        text = element_text(size=10),
        axis.text.x = element_text(angle=45, hjust=1),
        strip.text.x = element_text(size = 14, color="black"),
        legend.position = c(0.9, 0.1),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14),
        legend.background = element_rect(fill = "transparent", colour = NA)) +
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=.2, position=position_dodge(9))
bar_sir
  
#CREB

rsv_creb <- rsv_final  %>% filter(Protein=='CREB1')

rsv_creb$Sample <- factor(rsv_creb$Sample, levels=c("DMEM", "ETOH","ETOH_DMSO", 
                                                      "FFA","FFA_DMSO", "RSV_FFA", "RSV"))
anov_rsv_creb <- aov(FD ~ Sample, data=rsv_creb)
res.tukey <- TukeyHSD((anov_rsv_creb),"Sample",ordered=TRUE)
my_comparisons <- c('DMEM-ETOH', 'ETOH_DMSO-DMEM', 'FFA-DMEM', 'FFA_DMSO-DMEM',
                    'FFA-ETOH','ETOH_DMSO-FFA_DMSO','FFA_DMSO-RSV_FFA',
                    'ETOH_DMSO-RSV')

df <- as.data.frame(res.tukey[[1]]) %>% rownames_to_column('comparisons') %>% filter(comparisons %in% my_comparisons)
colnames(df)[5] <- 'p.adj'
df_rsv_creb <- df %>% separate(comparisons,c("group1","group2"),"-") %>%
  mutate(Treatment=group1) %>% 
  mutate(p.sig=ifelse(p.adj<0.0001,'<0.0001',format(round(p.adj, digits = 4),scientific = F)))

#ERK

rsv_erk <- rsv_final  %>% filter(Protein=='ERK')

rsv_erk$Sample <- factor(rsv_erk$Sample, levels=c("DMEM", "ETOH","ETOH_DMSO", 
                                                    "FFA","FFA_DMSO", "RSV_FFA", "RSV"))
anov_rsv_erk <- aov(FD ~ Sample, data=rsv_erk)
res.tukey <- TukeyHSD((anov_rsv_creb),"Sample",ordered=TRUE)
my_comparisons <- c('DMEM-ETOH', 'ETOH_DMSO-DMEM', 'FFA-DMEM', 'FFA_DMSO-DMEM',
                    'FFA-ETOH','ETOH_DMSO-FFA_DMSO','FFA_DMSO-RSV_FFA',
                    'ETOH_DMSO-RSV')

df <- as.data.frame(res.tukey[[1]]) %>% rownames_to_column('comparisons') %>% filter(comparisons %in% my_comparisons)
colnames(df)[5] <- 'p.adj'
df_rsv_erk <- df %>% separate(comparisons,c("group1","group2"),"-") %>%
  mutate(Treatment=group1) %>% 
  mutate(p.sig=ifelse(p.adj<0.0001,'<0.0001',format(round(p.adj, digits = 4),scientific = F)))

#MEK

rsv_mek <- rsv_final  %>% filter(Protein=='MEK')

rsv_mek$Sample <- factor(rsv_mek$Sample, levels=c("DMEM", "ETOH","ETOH_DMSO", 
                                                    "FFA","FFA_DMSO", "RSV_FFA", "RSV"))
anov_rsv_mek <- aov(FD ~ Sample, data=rsv_mek)
res.tukey <- TukeyHSD((anov_rsv_mek),"Sample",ordered=TRUE)
my_comparisons <- c('DMEM-ETOH', 'ETOH_DMSO-DMEM', 'FFA-DMEM', 'FFA_DMSO-DMEM',
                    'FFA-ETOH','FFA_DMSO-ETOH_DMSO','FFA_DMSO-RSV_FFA',
                    'ETOH_DMSO-RSV')

df <- as.data.frame(res.tukey[[1]]) %>% rownames_to_column('comparisons') %>% filter(comparisons %in% my_comparisons)
colnames(df)[5] <- 'p.adj'
df_rsv_mek <- df %>% separate(comparisons,c("group1","group2"),"-") %>%
  mutate(Treatment=group1) %>% 
  mutate(p.sig=ifelse(p.adj<0.0001,'<0.0001',format(round(p.adj, digits = 4),scientific = F)))

#PTN11

rsv_ptn <- rsv_final  %>% filter(Protein=='PTN11')

rsv_ptn$Sample <- factor(rsv_ptn$Sample, levels=c("DMEM", "ETOH","ETOH_DMSO", 
                                                  "FFA","FFA_DMSO", "RSV_FFA", "RSV"))
anov_rsv_ptn <- aov(FD ~ Sample, data=rsv_ptn)
res.tukey <- TukeyHSD((anov_rsv_ptn),"Sample",ordered=TRUE)
my_comparisons <- c('DMEM-ETOH', 'ETOH_DMSO-DMEM', 'FFA-DMEM', 'FFA_DMSO-DMEM',
                    'FFA-ETOH','FFA_DMSO-ETOH_DMSO','FFA_DMSO-RSV_FFA',
                    'ETOH_DMSO-RSV')

df <- as.data.frame(res.tukey[[1]]) %>% rownames_to_column('comparisons') %>% filter(comparisons %in% my_comparisons)
colnames(df)[5] <- 'p.adj'
df_rsv_ptn <- df %>% separate(comparisons,c("group1","group2"),"-") %>%
  mutate(Treatment=group1) %>% 
  mutate(p.sig=ifelse(p.adj<0.0001,'<0.0001',format(round(p.adj, digits = 4),scientific = F)))

#STAT3

rsv_stat <- rsv_final  %>% filter(Protein=='STAT3')

rsv_stat$Sample <- factor(rsv_stat$Sample, levels=c("DMEM", "ETOH","ETOH_DMSO", 
                                                  "FFA","FFA_DMSO", "RSV_FFA", "RSV"))
anov_rsv_stat <- aov(FD ~ Sample, data=rsv_stat)
res.tukey <- TukeyHSD((anov_rsv_stat),"Sample",ordered=TRUE)
my_comparisons <- c('DMEM-ETOH', 'ETOH_DMSO-DMEM', 'FFA-DMEM', 'FFA_DMSO-DMEM',
                    'FFA-ETOH','FFA_DMSO-ETOH_DMSO','FFA_DMSO-RSV_FFA',
                    'ETOH_DMSO-RSV')

df <- as.data.frame(res.tukey[[1]]) %>% rownames_to_column('comparisons') %>% filter(comparisons %in% my_comparisons)
colnames(df)[5] <- 'p.adj'
df_rsv_stat <- df %>% separate(comparisons,c("group1","group2"),"-") %>%
  mutate(Treatment=group1) %>% 
  mutate(p.sig=ifelse(p.adj<0.0001,'<0.0001',format(round(p.adj, digits = 4),scientific = F)))

#RS6

rsv_rs6 <- rsv_final  %>% filter(Protein=='RS6')

rsv_rs6$Sample <- factor(rsv_rs6$Sample, levels=c("DMEM", "ETOH","ETOH_DMSO", 
                                                    "FFA","FFA_DMSO", "RSV_FFA", "RSV"))
anov_rsv_rs6 <- aov(FD ~ Sample, data=rsv_rs6)
res.tukey <- TukeyHSD((anov_rsv_rs6),"Sample",ordered=TRUE)
my_comparisons <- c('DMEM-ETOH', 'ETOH_DMSO-DMEM', 'FFA-DMEM', 'DMEM-FFA_DMSO',
                    'FFA-ETOH','ETOH_DMSO-FFA_DMSO','RSV_FFA-FFA_DMSO',
                    'ETOH_DMSO-RSV')

df <- as.data.frame(res.tukey[[1]]) %>% rownames_to_column('comparisons') %>% filter(comparisons %in% my_comparisons)
colnames(df)[5] <- 'p.adj'
df_rsv_rs6 <- df %>% separate(comparisons,c("group1","group2"),"-") %>%
  mutate(Treatment=group1) %>% 
  mutate(p.sig=ifelse(p.adj<0.0001,'<0.0001',format(round(p.adj, digits = 4),scientific = F)))
