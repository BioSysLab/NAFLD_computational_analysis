library(AICcmodavg)
library(car)
library(multcomp)
library(tidyverse)
library(devtools)
library(rstatix)
library(ggpubr)
library(ggsignif)

data <- read.csv("C:/Users/PC/Desktop/Manuscript/Results/GraphPad_csv/ROS/ROS_Steat_2.csv", row.names = c("N1", "N2", "N3") ) %>%
  as_tibble() 
head (data)
summary (data)

raw <-data %>% gather (key = Treatment, value = FC)
head(raw)


se <- function(x) sqrt(var(x)/length(x))
stat <- raw %>% group_by(Treatment) %>% summarise(
  mean = mean(FC, na.rm = F),
  stdev = sd(FC, na.rm = F),
  sem = se(FC)
)

final <- left_join(raw, stat)
head(final)

### Tuckey test ###
final$Treatment <- factor(final$Treatment, levels=c("DMEM","DMSO", "H2O2",
                                                    "ETOH","FFA","CLO", 
                                                    "MEF","PIM"))

anov <- aov(FC ~ Treatment, data=final)

res.tukey <- TukeyHSD((anov),"Treatment",ordered=TRUE)

my_comparisons <- c('H2O2-DMEM',
                    'DMEM-DMSO','DMEM-ETOH','FFA-ETOH','CLO-DMSO',
                    'MEF-DMSO','PIM-DMSO')


df <- as.data.frame(res.tukey[[1]]) %>% rownames_to_column('comparisons') %>% filter(comparisons %in% my_comparisons)
colnames(df)[5] <- 'p.adj'
df <- df %>% separate(comparisons,c("group1","group2"),"-") %>%
  mutate(Treatment=group1) %>% 
  mutate(p.sig=ifelse(p.adj<0.0001,'<0.0001',format(round(p.adj, digits = 4),scientific = F)))

ggplot(final, aes(x=Treatment, y=mean, fill=Treatment)) + geom_col(position = position_dodge()) +
  geom_point(aes(y=final$FC)) +
  theme_classic() +
  scale_fill_manual(values = c("#628FCA", "#628FCA", "#628FCA", "#628FCA", "#E76152", "#F69780", "#F69780", "#F69780")) +
  scale_y_continuous(expand = c(0,0), limits = c(0,3.8)) + 
  ggtitle("Quantification of Intracellular ROS Producttion") + 
  xlab ("Treatment") + ylab("Fold-Change of Intracellular ROS Producttion \n (RFU/ug protein)") +
  theme(text = element_text(size = 16, family = "Sans"),
        plot.title = element_text(size=24, hjust = 0.5),
        legend.position = "none",
        axis.text.x=element_text(size=14, colour="black", angle = 45, hjust = 1),
        axis.text.y=element_text(size=14, colour="black")) +
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=.2, position=position_dodge(9)) +
  stat_pvalue_manual(df,label='p.sig',c(2.6,3.4,3.6,3.2,3,3,2.8),tip.length = 0.02, bracket.size=0.5, label.size=5,
                     step.increase = 0, vjust = 0.2)
ggsave("ROS_steat_DR.eps", scale=2, units = c("cm"), dpi=300, width=13, height=10)
ggsave("ROS_steat_DR.png", scale=2, units = c("cm"), dpi=300, width=13, height=10)
