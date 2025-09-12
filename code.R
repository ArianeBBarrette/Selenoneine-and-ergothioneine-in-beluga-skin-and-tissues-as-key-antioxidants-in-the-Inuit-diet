#Research Article for	Scientific Reports 
#Title: 	Distribution of Selenoneine and Ergothioneine in Beluga Skin and Tissues: Key Antioxidants in Inuit diet
#Authors: Ariane B. Barrette, Philippe Archambault, Mélanie Lemire, Corinne Zinflou, Nathalie Ouellet, Pierre Dumas, Adel Achouba, Jean-Éric Tremblay, Matthew Little and Pierre Ayotte
#Contact: ariane.barrette.2@ulaval.ca



# #Figure 1, data on adult beluga whales (1-14) concentrations of  --------
library(readr)
library(dplyr)
library(tidyverse)
#name file : SEN_EGT_tissues_adult_beluga_ABB 
beluga_tissus <- read_delim(file.choose(), delim = ";", trim_ws = TRUE, show_col_types = FALSE)
beluga_tissus <- mutate(beluga_tissus, Concentration = as.numeric(Concentration))

# Data table - descriptive stats  ----------------------------------------------------
std_err <- function(x, n) sd(x) / sqrt(n)
library(inTextSummaryTable)
library(data.table)
library(gt)
beluga_data_summ <- beluga_tissus %>% mutate(Concentration = as.numeric(Concentration)) %>%
setDT()
donnee_descriptive <- beluga_data_summ[, .(mean = mean(Concentration), geom_mean = geomMean(Concentration), min = min(Concentration),
                                             max = max(Concentration),
                                             median = median(Concentration),
                                             percentile25 = quantile(Concentration, probs = 0.25),
                                             percentile75 = quantile(Concentration, probs = 0.75),
                                             geom_se = geomSE(Concentration),
                                             se = std_err(Concentration, .N), n = .N ),
                                         by=.(Tissue, Analyte)]

gt(donnee_descriptive)

# #comparaison of tissues for all 4 analytes in adult beluga whales,
#one-way analysis of variance (ANOVA) followed by the post-hoc Tukey pairwise comparison analysi----------------------------------------------------
library(dplyr)
library(ggplot2)
library(Rmisc)
library(rstatix)
library(multcompView)

# Helper: replace non-positive with NA for log, and LOD/2 if you prefer
to_log <- function(x) ifelse(x > 0, x, NA_real_)
halfLOD <- function(x, lod) ifelse(is.na(x) | x <= 0, lod/2, x)

# ---- SELENONEINE (SEN) ----
SENto <- beluga_tissus %>%
  filter(Analyte == "Selenoneine")

# one way ANOVA 
SENtoresult <- aov(Concentration ~ Tissue, data = SENto)
summary(SENtoresult)
#Tukey
TukeySEN <- TukeyHSD(SENtoresult, ordered = TRUE, conf.level = 0.95)
TukeySEN.cld <- multcompLetters4(SENtoresult, TukeySEN); print(TukeySEN.cld)

# Assumptions
layout(matrix(c(1:5, 5), 2, 3))
plot(SENtoresult, which = c(1,3:5), pch = 20, cex = 1.5)
qqnorm(residuals(SENtoresult)); qqline(residuals(SENtoresult))
layout(1)
shapiro.test(SENto$Concentration)
#not respected

# log-transformation
SENto <- SENto %>% mutate(conclog = log(Concentration)) 
SENtoresult <- aov(conclog ~ Tissue, data = SENto)
summary(SENtoresult)
#Tukey
TukeySEN <- TukeyHSD(SENtoresult, ordered = TRUE, conf.level = 0.95)
TukeySEN.cld <- multcompLetters4(SENtoresult, TukeySEN); print(TukeySEN.cld)

layout(matrix(c(1:5, 5), 2, 3))
plot(SENtoresult, which = c(1,3:5), pch = 20, cex = 1.5)
qqnorm(residuals(SENtoresult)); qqline(residuals(SENtoresult))
layout(1)
shapiro.test(SENto$conclog)
#much better



## ---- ERGOTHIONEINE (EGT) ----
EGTto <- beluga_tissus %>% 
  filter(Analyte == "Ergothioneine")

# one-way ANOVA
EGTtoresult <- aov(Concentration ~ Tissue, data = EGTto)
summary(EGTtoresult)

# Tukey
TukeyEGT <- TukeyHSD(EGTtoresult, ordered = TRUE, conf.level = 0.95)
TukeyEGT.cld <- multcompLetters4(EGTtoresult, TukeyEGT); print(TukeyEGT.cld)

# Assumptions
layout(matrix(c(1:5, 5), 2, 3))
plot(EGTtoresult, which = c(1,3:5), pch = 20, cex = 1.5)
qqnorm(residuals(EGTtoresult)); qqline(residuals(EGTtoresult))
layout(1)
shapiro.test(EGTto$Concentration)

# log-transform
EGTto <- EGTto %>% mutate(conclog = log(Concentration))
EGTtoresult <- aov(conclog ~ Tissue, data = EGTto)
summary(EGTtoresult)

TukeyEGT <- TukeyHSD(EGTtoresult, ordered = TRUE, conf.level = 0.95)
TukeyEGT.cld <- multcompLetters4(EGTtoresult, TukeyEGT); print(TukeyEGT.cld)

layout(matrix(c(1:5, 5), 2, 3))
plot(EGTtoresult, which = c(1,3:5), pch = 20, cex = 1.5)
qqnorm(residuals(EGTtoresult)); qqline(residuals(EGTtoresult))
layout(1)
shapiro.test(EGTto$conclog)
#Much better

## ---- Se-methyl-selenoneine (Se_methyl_selenoneine) ----
SeMethylSEN <- beluga_tissus %>% 
  filter(Analyte == "Se_methyl_selenoneine")

# one-way ANOVA
SeMethylSEN_res <- aov(Concentration ~ Tissue, data = SeMethylSEN)
summary(SeMethylSEN_res)

# Tukey
Tukey_SeMSEN <- TukeyHSD(SeMethylSEN_res, ordered = TRUE, conf.level = 0.95)
Tukey_SeMSEN.cld <- multcompLetters4(SeMethylSEN_res, Tukey_SeMSEN); print(Tukey_SeMSEN.cld)

# Assumptions
layout(matrix(c(1:5, 5), 2, 3))
plot(SeMethylSEN_res, which = c(1,3:5), pch = 20, cex = 1.5)
qqnorm(residuals(SeMethylSEN_res)); qqline(residuals(SeMethylSEN_res))
layout(1)
shapiro.test(SeMethylSEN$Concentration)

# log-transform (use LOD/2 for non-positive values)
SeMethylSEN <- SeMethylSEN %>% 
  mutate(conc_adj = ifelse(Concentration > 0, Concentration, Limit_of_detection/2),
         conclog  = log(conc_adj))

SeMethylSEN_res <- aov(conclog ~ Tissue, data = SeMethylSEN)
summary(SeMethylSEN_res)

Tukey_SeMSEN <- TukeyHSD(SeMethylSEN_res, ordered = TRUE, conf.level = 0.95)
Tukey_SeMSEN.cld <- multcompLetters4(SeMethylSEN_res, Tukey_SeMSEN); print(Tukey_SeMSEN.cld)

layout(matrix(c(1:5, 5), 2, 3))
plot(SeMethylSEN_res, which = c(1,3:5), pch = 20, cex = 1.5)
qqnorm(residuals(SeMethylSEN_res)); qqline(residuals(SeMethylSEN_res))
layout(1)
shapiro.test(SeMethylSEN$conclog)

###For methylated compounds, as many were under the LOD, results are not statistically powerful
## Results were presented only as descriptive analyses


## ---- S-methyl-ergothioneine (S_methyl_ergothioneine) ----
SMethylEGT <- beluga_tissus %>% 
  filter(Analyte == "S_methyl_ergothioneine")

# one-way ANOVA
SMethylEGT_res <- aov(Concentration ~ Tissue, data = SMethylEGT)
summary(SMethylEGT_res)

# Tukey
Tukey_SMEGT <- TukeyHSD(SMethylEGT_res, ordered = TRUE, conf.level = 0.95)
Tukey_SMEGT.cld <- multcompLetters4(SMethylEGT_res, Tukey_SMEGT); print(Tukey_SMEGT.cld)

# Assumptions
layout(matrix(c(1:5, 5), 2, 3))
plot(SMethylEGT_res, which = c(1,3:5), pch = 20, cex = 1.5)
qqnorm(residuals(SMethylEGT_res)); qqline(residuals(SMethylEGT_res))
layout(1)
shapiro.test(SMethylEGT$Concentration)

# log-transform (use LOD/2 for non-positive values)
SMethylEGT <- SMethylEGT %>% 
  mutate(conc_adj = ifelse(Concentration > 0, Concentration, Limit_of_detection/2),
         conclog  = log(conc_adj))

SMethylEGT_res <- aov(conclog ~ Tissue, data = SMethylEGT)
summary(SMethylEGT_res)

Tukey_SMEGT <- TukeyHSD(SMethylEGT_res, ordered = TRUE, conf.level = 0.95)
Tukey_SMEGT.cld <- multcompLetters4(SMethylEGT_res, Tukey_SMEGT); print(Tukey_SMEGT.cld)

layout(matrix(c(1:5, 5), 2, 3))
plot(SMethylEGT_res, which = c(1,3:5), pch = 20, cex = 1.5)
qqnorm(residuals(SMethylEGT_res)); qqline(residuals(SMethylEGT_res))
layout(1)
shapiro.test(SMethylEGT$conclog)

###For methylated compounds, as many were under the LOD, results are not statistically powerful
## Results were presented only as descriptive analyses


# Graphs for Figure 2  ----------------------------------------------------------------
library(ggbreak)
library(wesanderson)
library(ggplot2)
library(ggthemes)
library(gridExtra)
#general theme
theme_1 <-    theme_foundation(base_size=20, base_family="")+
  theme(axis.title.x =  element_text(size=30), 
        title = element_text(size =30),
        plot.title = element_text(hjust = -0.05, size = 30, face = "bold"),
        axis.title.y =  element_text(size=30), 
        axis.text.x = element_text(size=25),
        axis.text.y = element_text(size=25),
        panel.background = element_rect(colour = "black"),
        plot.background = element_rect(colour = NA, size =1),
        panel.border = element_rect(colour = NA, size = 1),
        axis.line = element_line(colour="black", size = 0.5, 
                                 arrow = arrow(type = "closed", length = unit(0.05, "inches"), angle = 90)),
        axis.ticks = element_line(),
        panel.grid.major =  element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text( size = 30, face = "bold"), 
        strip.background = element_rect(fill = "grey"))


##boxplots 
'Selenoneine figure 1a'
print(TukeySEN.cld)
cld = data.frame(tukey=c("d", "b,c","c,d","b","e","f","a"), 
                 Tissue = c("Blood", "Brain", "Intestine", "Kidney", "Liver", "Muscle", "Skin"))
SENmax <- SENto %>% aggregate(Concentration ~ Analyte + Tissue, FUN = max)
cld = cbind(SENmax$Concentration, cld)
SENtotukey <- merge(cld, SENto,  by = "Tissue")

SENtographbox<- ggplot(SENtotukey, aes(x=Tissue, y=Concentration)) + 
  geom_boxplot( outlier.colour="black", outlier.shape= 16,
                outlier.size=2, notch = FALSE, fill = wes_palette(n=1, name="Zissou1"))+
  labs(x = "Tissues",  y = "Selenoneine concentration (µg/g)") +
  facet_wrap(~Analyte, 
             labeller = labeller(Analyte = c("seleno" = "Selenoneine")))+
  theme(axis.title.x = element_text(size = 20), axis.text.x = element_text(size=15),  
        axis.text.y = element_text(size=20), face = "bold")+
  
  scale_y_break(breaks = c(4,7), scale = 1, expand = c(0,0), space =0.5)+
  geom_text(data = filter(SENtotukey, Whale_ID == 13), 
            aes(y = (SENmax$Concentration),label = tukey), vjust=-0.5, size = 10)+
  ggtitle("a")+
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0,35)) +
  theme_foundation(base_size=20, base_family="")+
  theme_1+
  theme(axis.text.x.top  = element_blank(),
        axis.ticks.x.top = element_blank(), 
        axis.line.x.top = element_blank(),
        axis.text.y.right  = element_blank(),
        axis.ticks.y.right = element_blank(), 
        axis.line.y.right = element_blank() )
SENtographbox

'EGT tissues figure 1b'
print(TukeyEGT.cld)
cld = data.frame(tukey=c("b,c","c","b","d","a","e", "d"), 
                 Tissue = c("Brain", "Intestine", "Kidney", "Liver", "Skin", "Muscle", "Blood"))
EGTmax <- EGTto %>% aggregate(Concentration ~ Analyte + Tissue, FUN = max)
cld = cbind(EGTmax$Concentration, cld)
ERGOtotukey <- merge(cld, EGTto, by = "Tissue")

ERGOtographbox <- ggplot(ERGOtotukey, aes(x=Tissue, y=Concentration)) + 
  geom_boxplot( outlier.colour="black", outlier.shape= 16,
                outlier.size=2, notch = FALSE, fill = wes_palette(n=1, name="FantasticFox1"))+
  labs(x = "Tissues",  y = "Ergothioneine concentration (µg/g)") +
  facet_wrap(~Analyte, 
             labeller = labeller(Analyte = c("ergo" = "Ergothioneine")))+
  theme(axis.title.x = element_text(size = 20), axis.text.x = element_text(size=15),  
        axis.text.y = element_text(size=20), face = "bold")+
  geom_text(data = filter(ERGOtotukey, Whale_ID == 13), 
            aes(y = (EGTmax$Concentration),label = tukey), vjust=-0.5, size = 10)+
  scale_y_break(breaks = c(15,35), scale = 1, expand = c(0,0), space =0.5)+
  ggtitle("b")+
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0,200)) +
  theme_foundation(base_size=20, base_family="")+
  theme_1+
  theme(axis.text.x.top  = element_blank(),
        axis.ticks.x.top = element_blank(), 
        axis.line.x.top = element_blank(),
        axis.text.y.right  = element_blank(),
        axis.ticks.y.right = element_blank(), 
        axis.line.y.right = element_blank())
ERGOtographbox 


library(aplot)
figure2 <- plot_list(gglist=list(SENtographbox , ERGOtographbox ))
figure2
#save as PDF 20X25
ggsave(figure2, 
       filename = "figure2.jpg",
       device = "jpg",
       height = 10, width = 25, units = "in",
       path = "downloads")



# Graphs for Supplementary Figure 1  ----------------------------------------------------------------

'Se-methyl-SEN tissues supp. figure 1a'
SeMethylSEN_plot <- summarySE(SeMethylSEN, measurevar="Concentration", groupvars=c("Tissue", "Analyte"))
SeMethylSEN_plot_all <- merge(SeMethylSEN_plot, SeMethylSEN, by = c("Tissue", "Analyte"))


MSENtographf <- ggplot(data=SeMethylSEN_plot_all, aes(x=Tissue, y=Concentration.x)) +
  geom_bar(stat = "identity", width = 0.7, fill = wes_palette(n=1, name="Moonrise3"), color = "black", position = "identity") +
  labs(x = "Tissues",  y = "Concentration (ng/g)") +
  geom_point(aes(y = Concentration.y), position = position_jitter(width = 0.2), color = "black", alpha = 1/2, size =3) +
  scale_y_continuous(expand = c(0, 0), 
                     limits = c(0,90))+
  geom_errorbar(aes(ymin=Concentration.x-se, ymax=Concentration.x+se), width=.2 ) +
  theme_1+
  facet_wrap(~Analyte, 
             labeller = labeller(Analyte = c("Se_methyl_selenoneine" = "Se-Methyl-Selenoneine")))+
  scale_y_break(breaks = c(64,83), scale = 0.05, expand = c(0,0), space =0.5)+
  geom_abline(slope = 0)+
  geom_abline(slope = 0, intercept = 64, size = 1)+
  geom_abline(slope = 0, intercept = 83, size = 1)+
  ggtitle("a")+
  annotate("segment", x = 0.1, xend = 7.5, y = 2, yend = 2,
           colour = "red", alpha = 1/2, size=2)+
  annotate("text", label = expression(bold("ND")), color = "red", x = 0.4, y = 4, size = 8)+
  coord_cartesian(xlim = c(0.5, 7),
                  clip = 'off') +
  theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(), 
        axis.line.y.right = element_blank(),
        plot.title = element_text(hjust=-0.03, vjust =-2, size = 40, face = "bold"),
        plot.margin = unit(c(1,1,1,3), "lines"))
MSENtographf 


'S-methyl-ergothioneine tissues supp. figure 1b'
SMethylEGT_plot <- summarySE(SMethylEGT, measurevar="Concentration", groupvars=c("Tissue", "Analyte"))
SMethylEGT_plot_all <- merge(SMethylEGT_plot, SMethylEGT, by = c("Tissue", "Analyte"))


MEGTtographf <- ggplot(data=SMethylEGT_plot_all, aes(x=Tissue, y=Concentration.x)) +
  geom_bar(stat = "identity", width = 0.7, fill = wes_palette(n=1, name="Moonrise3"), color = "black", position = "identity") +
  labs(x = "Tissues",  y = "Concentration (ng/g)") +
  geom_point(aes(y = Concentration.y), position = position_jitter(width = 0.2), color = "black", alpha = 1/2, size =3) +
  scale_y_continuous(expand = c(0, 0), 
                     limits = c(0,60)) +
  geom_errorbar(aes(ymin=Concentration.x-se, ymax=Concentration.x+se), width=.2 ) +
  theme_1+
  facet_wrap(~Analyte, 
             labeller = labeller(Analyte = c("S_methyl_ergothioneine" = "S-Methyl-Ergothioneine")))+
  ggtitle("b")+
  annotate("segment", x = 0.5, xend = 7.5, y = 6.32, yend = 6.32,
           colour = "red", alpha = 1/2, size =2)+
  annotate("text", x = 2, y =2 , label = "ND", size = 8)+
  annotate("text", x = 7, y =58 , label = "NS", size = 10)+
  coord_cartesian(xlim = c(0.5, 7),
                  clip = 'off') +
  theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(), 
        axis.line.y.right = element_blank(),
        plot.title = element_text(hjust=-0.03, vjust =-2, size = 40, face = "bold"),
        plot.margin = unit(c(1,1,1,3), "lines"))
MEGTtographf 


'S-methyl-EGT organs/tissues figure 1d'
brain <- beluga %>% filter( tissue == "brain", mesure == "methergo") %>%
  aggregate(conc ~  whale + tissue + mesure + LOD, FUN = mean)

MERGOto <- beluga %>% filter(mesure == "methergo", tissue != "skin", peau == "dorsal" | is.na(peau)) %>%
  aggregate(conc ~ whale + tissue + mesure + LOD, FUN = mean) %>%
  rbind(brain)

print(TukeyMERGO.cld)
cld = data.frame(tukey=c("", "a","a","a","a","a", "a"), tissue = c("brain", "int", "kidney", "liver", "mattaaq", "muscle", "sang"))
MERGOto <- mutate(MERGOto, concng = conc*1000)
MERGOto <- mutate(MERGOto, LODng = LOD*1000)
MERGOtograph <- MERGOto %>% summarySE(measurevar="concng", groupvars=c("tissue", "mesure")) 
MERGOtotukey <- merge(cld, MERGOto, by = "tissue")
MERGOtotukey <- merge(MERGOtotukey, MERGOtograph, by = c("tissue", "mesure"))

MERGOtographf <- ggplot(data=MERGOtotukey, aes(x=tissue, y=concng.y)) +
  geom_bar(stat = "identity", width = 0.7, fill = wes_palette(n=1, name="Moonrise1"), color = "black", position = "identity") +
  labs(x = "Organs/Tissues",  y = "Concentration (ng/g)") +
  geom_point(aes(y = concng.x), position = position_jitter(width = 0.2), color = "black", alpha = 1/2, size=3) +
  scale_y_continuous(expand = c(0, 0), 
                     limits = c(0,60)) +
  geom_errorbar(aes(ymin=concng.y-se, ymax=concng.y+se), width=.2 ) +
  theme_1+
  ggtitle("b")+
  annotate("segment", x = 0.5, xend = 6.5, y = 6.32, yend = 6.32,
           colour = "red", alpha = 1/2, size =2)+
  annotate("segment", x = 6.5, xend = 7.5, y = 1.137121, yend = 1.137121,
           colour = "red", alpha = 1/2, size =2)+
  annotate("text", x = 1, y =2 , label = "ND", size = 8)+
  annotate("text", x = 7, y =58 , label = "NS", size = 10)+
  scale_x_discrete(labels = c("Brain", "Intestine", "Kidney", "Liver", "Skin", "Muscle", "Blood"))+
  facet_wrap(~mesure, 
             labeller = labeller(mesure = c("methergo" = "S-Methyl-Ergothioneine")))+
  geom_abline(slope = 0)+
  theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(), 
        axis.line.y.right = element_blank(), 
        axis.title.y.left = element_blank(),
        plot.title = element_text(hjust=-0.11, vjust =-2, size = 40, face = "bold"),
        plot.margin = unit(c(1,1,1,3), "lines"))
MERGOtographf 

library(aplot)
figureext2 <- plot_list(gglist=list(MSENtographf ,MERGOtographf))
figureext2
#save as PDF 20X25
ggsave(figureext2, 
       filename = "figureext2.jpg",
       device = "jpg",
       height = 10, width = 25, units = "in",
       path = "downloads")








# DATA reorganisation for correlations ------------------------------------
 
"correlations"
seleno <- beluga_tissus %>%  aggregate(Concentration ~ Analyte + Whale_ID + Tissue, FUN = mean) %>%
  filter(Analyte == "Selenoneine") %>% 
  subset(select = -c(Analyte)) %>% drop_na(Concentration)
colnames(seleno) <- c( "whale", "tissue","seleno")

ergo <- beluga_tissus %>%  aggregate(Concentration ~ Analyte + Whale_ID + Tissue, FUN = mean) %>%
  filter(Analyte == "Ergothioneine") %>% 
  subset(select = -c(Analyte)) %>% drop_na(Concentration)
colnames(ergo) <- c( "whale", "tissue","ergo")

methergo <- beluga_tissus %>%  aggregate(Concentration ~ Analyte + Whale_ID + Tissue, FUN = mean) %>%
  filter(Analyte == "S_methyl_ergothioneine") %>% 
  subset(select = -c(Analyte)) %>% drop_na(Concentration)
colnames(methergo) <- c( "whale", "tissue","methergo")

methseleno <- beluga_tissus %>%  aggregate(Concentration ~ Analyte + Whale_ID + Tissue, FUN = mean) %>%
  filter(Analyte == "Se_methyl_selenoneine") %>% 
  subset(select = -c(Analyte)) %>% drop_na(Concentration)
colnames(methseleno) <- c( "whale", "tissue","methseleno")

# #correlation S-methyl-EGT vs EGT ** non statistically significant ** -------------------------------------
corr_MERGO_ERGO <- join(methergo, ergo)

corr_MERGO_ERGO_rein <- corr_MERGO_ERGO %>% filter( tissue == "Kidney") 
lm <- lm(ergo ~ methergo, data = corr_MERGO_ERGO_rein)
summary(lm)
#p>0.05 so no significant correlations
layout(matrix(c(1:5, 5), 2, 3)) +
  plot(lm, c(1, 3:5), pch = 20, cex = 1.5) 
qqnorm(lm$residuals)  
layout(1)
MERGO_ERGO_rein <- ggplot(corr_MERGO_ERGO_rein , mapping = aes(y= ergo, x = methergo))+
  geom_point( shape = 21, fill = 'black', color = 'white', size = 5)+
  geom_smooth(method = "lm", se=T, color="black", fill = wes_palette(n=1, name = "Royal1"), formula = y ~ x)+
  facet_wrap(~tissue,
             labeller = labeller(tissue = c("kidney" = "Kidney")))+
  labs(x = "S-Methyl-Ergothioneine (nmol/g)", y = "Ergothioneine (nmol/g)")+
  theme_1
MERGO_ERGO_rein

corr_MERGO_ERGO_liver <- corr_MERGO_ERGO %>% filter( tissue == "Liver") 
lm <- lm(ergo ~ methergo, data = corr_MERGO_ERGO_liver)
summary(lm)
#p>0.05 so no significant correlations
layout(matrix(c(1:5, 5), 2, 3)) +
  plot(lm, c(1, 3:5), pch = 20, cex = 1.5) +
  qqnorm(lm$residuals)  +
  layout(1)
MERGO_ERGO_liver <- ggplot(corr_MERGO_ERGO_liver , mapping = aes(y= ergo, x = methergo))+
  geom_point( shape = 21, fill = 'black', color = 'white', size = 5)+
  geom_smooth(method = "lm", se=T, color="black", fill = wes_palette(n=1, name = "Royal1"), formula = y ~ x)+
  facet_wrap(~tissue,
             labeller = labeller(tissue = c("liver" = "Liver")))+
  labs(x = "S-Methyl-Ergothioneine (nmol/g)", y = "Ergothioneine (nmol/g)")+
  theme_1
MERGO_ERGO_liver

corr_MERGO_ERGO_gut <- corr_MERGO_ERGO %>% filter( tissue == "Intestine") 
lm <- lm(ergo ~ methergo, data = corr_MERGO_ERGO_gut)
summary(lm)
#p>0.05 so no significant correlations
layout(matrix(c(1:5, 5), 2, 3)) +
  plot(lm, c(1, 3:5), pch = 20, cex = 1.5) +
  qqnorm(lm$residuals)  +
  layout(1)
MERGO_ERGO_gut <- ggplot(corr_MERGO_ERGO_gut , mapping = aes(y= ergo, x = methergo))+
  geom_point( shape = 21, fill = 'black', color = 'white', size = 5)+
  geom_smooth(method = "lm", se=T, color="black", fill = wes_palette(n=1, name = "Royal1"), formula = y ~ x)+
  facet_wrap(~tissue,
             labeller = labeller(tissue = c("int" = "Gut")))+
  labs(x = "S-Methyl-Ergothioneine (nmol/g)", y = "Ergothioneine (nmol/g)")+
  theme_1
MERGO_ERGO_gut

corr_MERGO_ERGO_muscle <- corr_MERGO_ERGO %>% filter( tissue == "Muscle") 
lm <- lm(ergo ~ methergo, data = corr_MERGO_ERGO_muscle)
summary(lm)
#p>0.05 so no significant correlations
layout(matrix(c(1:5, 5), 2, 3)) +
  plot(lm, c(1, 3:5), pch = 20, cex = 1.5) +
  qqnorm(lm$residuals)  +
  layout(1)
MERGO_ERGO_muscle <- ggplot(corr_MERGO_ERGO_muscle , mapping = aes(y= ergo, x = methergo))+
  geom_point( shape = 21, fill = 'black', color = 'white', size = 5)+
  geom_smooth(method = "lm", se=T, color="black", fill = wes_palette(n=1, name = "Royal1"), formula = y ~ x)+
  facet_wrap(~tissue,
             labeller = labeller(tissue = c("muscle" = "Muscle")))+
  labs(x = "S-Methyl-Ergothioneine (nmol/g)", y = "Ergothioneine (nmol/g)")+
  theme_1
MERGO_ERGO_muscle

corr_MERGO_ERGO_blood <- corr_MERGO_ERGO %>% filter( tissue == "Blood") %>% drop_na(ergo)
lm <- lm(ergo ~ methergo, data = corr_MERGO_ERGO_blood)
summary(lm)
#p>0.05 so no significant correlations
layout(matrix(c(1:5, 5), 2, 3)) +
  plot(lm, c(1, 3:5), pch = 20, cex = 1.5) +
  qqnorm(lm$residuals)  +
  layout(1)
MERGO_ERGO_blood <- ggplot(corr_MERGO_ERGO_blood , mapping = aes(y= ergo, x = methergo))+
  geom_point( shape = 21, fill = 'black', color = 'white', size = 5)+
  geom_smooth(method = "lm", se=T, color="black", fill = wes_palette(n=1, name = "Royal1"), formula = y ~ x)+
  facet_wrap(~tissue,
             labeller = labeller(tissue = c("sang" = "Whole Blood")))+
  labs(x = "S-Methyl-Ergothioneine (nmol/g)", y = "Ergothioneine (nmol/g)")+
  theme_1
MERGO_ERGO_blood

#not enough data above LOD for brain
corr_MERGO_ERGO_brain <- corr_MERGO_ERGO %>% filter( tissue == "Brain") 
corr_MERGO_ERGO_brain

corr_MERGO_ERGO_peau <- corr_MERGO_ERGO %>% filter( tissue == "Skin") 
lm <- lm(ergo ~ methergo, data = corr_MERGO_ERGO_peau)
summary(lm)
#p>0.05 so no significant correlations
layout(matrix(c(1:5, 5), 2, 3)) +
  plot(lm, c(1, 3:5), pch = 20, cex = 1.5) +
  qqnorm(lm$residuals)  +
  layout(1)
MERGO_ERGO_peau <- ggplot(corr_MERGO_ERGO_peau , mapping = aes(y= ergo, x = methergo))+
  geom_point( shape = 21, fill = 'black', color = 'white', size = 5)+
  geom_smooth(method = "lm", se=T, color="black", fill = wes_palette(n=1, name = "Royal1"), formula = y ~ x)+
  facet_wrap(~tissue,
             labeller = labeller(tissue = c("Skin" = "Skin")))+
  labs(x = "S-Methyl-Ergothioneine (nmol/g)", y = "Ergothioneine (nmol/g)")+
  theme_1
MERGO_ERGO_peau

"no correlations are significative"

grid.arrange(MERGO_ERGO_peau, MERGO_ERGO_liver, MERGO_ERGO_rein, MERGO_ERGO_gut,
             MERGO_ERGO_muscle, MERGO_ERGO_blood, 
             ncol =2)
#save as 45x40

# #correlation S-methyl-SEN vs SEN ** non statistically significant ** -------------------------------------
corr_MSEN_SEN <- join(methseleno, seleno)

corr_MSEN_SEN_rein <- corr_MSEN_SEN %>% filter( tissue == "Kidney") 
lm <- lm(seleno ~ methseleno, data = corr_MSEN_SEN_rein)
summary(lm)
#p>0.05 so no significant correlations
layout(matrix(c(1:5, 5), 2, 3)) +
  plot(lm, c(1, 3:5), pch = 20, cex = 1.5) +
  qqnorm(lm$residuals)  +
  layout(1)
MSEN_SEN_rein <- ggplot(corr_MSEN_SEN_rein , mapping = aes(y= seleno, x = methseleno))+
  geom_point( shape = 21, fill = 'black', color = 'white', size = 5)+
  geom_smooth(method = "lm", se=T, color="black", fill = wes_palette(n=1, name = "Royal1"), formula = y ~ x)+
  facet_wrap(~tissue,
             labeller = labeller(tissue = c("Kidney" = "Kidney")))+
  labs(x = "Se-Methyl-Selenoneine (nmol/g)", y = "Selenoneine (nmol/g)")+
  theme_1
MSEN_SEN_rein

corr_MSEN_SEN_int <- corr_MSEN_SEN %>% filter( tissue == "Intestine") 
lm <- lm(seleno ~ methseleno, data = corr_MSEN_SEN_int)
summary(lm)
#p>0.05 so no significant correlations
layout(matrix(c(1:5, 5), 2, 3)) +
  plot(lm, c(1, 3:5), pch = 20, cex = 1.5) +
  qqnorm(lm$residuals)  +
  layout(1)
MSEN_SEN_int <- ggplot(corr_MSEN_SEN_int , mapping = aes(y= seleno, x = methseleno))+
  geom_point( shape = 21, fill = 'black', color = 'white', size = 5)+
  geom_smooth(method = "lm", se=T, color="black", fill = wes_palette(n=1, name = "Royal1"), formula = y ~ x)+
  theme_1
MSEN_SEN_int

corr_MSEN_SEN_mattaaq <- corr_MSEN_SEN %>% filter( tissue == "Skin") 
lm <- lm(seleno ~ methseleno, data = corr_MSEN_SEN_mattaaq)
summary(lm)
#p>0.05 so no significant correlations
layout(matrix(c(1:5, 5), 2, 3)) +
  plot(lm, c(1, 3:5), pch = 20, cex = 1.5) +
  qqnorm(lm$residuals)  +
  layout(1)
MSEN_SEN_mattaaq <- ggplot(corr_MSEN_SEN_mattaaq , mapping = aes(y= seleno, x = methseleno))+
  geom_point( shape = 21, fill = 'black', color = 'white', size = 5)+
  geom_smooth(method = "lm", se=T, color="black", fill = wes_palette(n=1, name = "Royal1"), formula = y ~ x)+
  facet_wrap(~tissue,
             labeller = labeller(tissue = c("Skin" = "Skin (dorsal)")))+
  labs(x = "Se-Methyl-Selenoneine (nmol/g)", y = "Selenoneine (nmol/g)")+
  theme_1
MSEN_SEN_mattaaq

corr_MSEN_SEN_blood <- corr_MSEN_SEN %>% filter( tissue == "Blood")
lm <- lm(seleno ~ methseleno, data = corr_MSEN_SEN_blood)
summary(lm)
#significant stat
layout(matrix(c(1:5, 5), 2, 3)) +
  plot(lm, c(1, 3:5), pch = 20, cex = 1.5) 
qqnorm(lm$residuals)  
layout(1)
shapiro.test(corr_MSEN_SEN$seleno)
shapiro.test(corr_MSEN_SEN$methseleno)
#accepted postulat

MSEN_SEN_blood <- ggplot(corr_MSEN_SEN_blood , mapping = aes(y= seleno, x = methseleno))+
  geom_point( shape = 21, fill = 'black', color = 'white', size = 5)+
  geom_smooth(method = "lm", se=T, color="black", fill = wes_palette(n=1, name = "Royal1"), formula = y ~ x)+
  facet_wrap(~tissue,
             labeller = labeller(tissue = c("Blood" = "Blood")))+
  labs(x = "Se-Methyl-Selenoneine (nmol/g)", y = "Selenoneine (nmol/g)")+
  theme_1
MSEN_SEN_blood

#heavily influenced by two data (cook's distance above 2)-both from different hunting season -remove to see

corr_MSEN_SEN_blood <- corr_MSEN_SEN %>% filter( tissue == "Blood")
corr_MSEN_SEN_blood <- corr_MSEN_SEN_blood[-c(1,2), ]
lm <- lm(seleno ~ methseleno, data = corr_MSEN_SEN_blood)
summary(lm)
#not statistically significant anymore 

#not enough data above LOD
corr_MSEN_SEN_muscle <- corr_MSEN_SEN %>% filter( tissue == "Muscle") 
corr_MSEN_SEN_brain<- corr_MSEN_SEN %>% filter( tissue == "Brain") 
corr_MSEN_SEN_liver <- corr_MSEN_SEN %>% filter( tissue == "Liver") 




# Correlation ergo/seleno -- figure 3 -------------------------------------

corr_ERGO_SEN <- join(ergo, seleno)

corr_ERGO_SEN_rein <- corr_ERGO_SEN %>% filter( tissue == "Kidney") 
lm <- lm(ergo ~ seleno, data = corr_ERGO_SEN_rein)
summary(lm)
layout(matrix(c(1:5, 5), 2, 3)) 
plot(lm, c(1, 3:5), pch = 20, cex = 1.5)
qqnorm(lm$residuals)  
layout(1)
shapiro.test(corr_ERGO_SEN_rein$ergo)
shapiro.test(corr_ERGO_SEN_rein$seleno)

#correlation
library(smplot2)
N_rein <- nrow(corr_ERGO_SEN_rein)
cor.test(corr_ERGO_SEN_rein$seleno, corr_ERGO_SEN_rein$ergo)
ERGO_SEN_rein <- ggplot(corr_ERGO_SEN_rein , mapping = aes(y= ergo, x = seleno))+
  geom_point( shape = 21, fill = 'black', color = 'white', size = 5)+
  sm_statCorr(color = "grey", text_size	= 10)+
  annotate("text", x = min(corr_ERGO_SEN_rein$seleno), 
           y = max(corr_ERGO_SEN_rein$ergo), 
           label = paste0("n = ", N_rein),
           hjust = 0, vjust = 3, size = 10, color = "black") + 
  facet_wrap(~tissue,
             labeller = labeller(tissue = c("Kidney" = "Kidney")))+
  labs(x = "Selenoneine (µg/g)", y = "Ergothioneine (µg/g)")+
  ggtitle(expression(bold("b")))+
  theme_1+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())
ERGO_SEN_rein


corr_ERGO_SEN_brain <- corr_ERGO_SEN %>% filter( tissue == "Brain") 
lm <- lm(ergo ~ seleno, data = corr_ERGO_SEN_brain)
summary(lm)
layout(matrix(c(1:5, 5), 2, 3)) +
  plot(lm, c(1, 3:5), pch = 20, cex = 1.5) 
  qqnorm(lm$residuals)  
  layout(1)

#correlation de pearson
library(smplot2)
N_brain <- nrow(corr_ERGO_SEN_brain)
cor.test(corr_ERGO_SEN_brain$seleno, corr_ERGO_SEN_brain$ergo)
ERGO_SEN_brain <- ggplot(corr_ERGO_SEN_brain , mapping = aes(y= ergo, x = seleno))+
  geom_point( shape = 21, fill = 'black', color = 'white', size = 5)+
  sm_statCorr(color = "grey", text_size	= 10)+
  facet_wrap(~tissue,
             labeller = labeller(tissue = c("Brain" = "Brain")))+
  annotate("text", x = min(corr_ERGO_SEN_brain$seleno), 
           y = max(corr_ERGO_SEN_brain$ergo), 
           label = paste0("n = ", N_brain),
           hjust = 0, vjust = 3, size = 10, color = "black") + 
  labs(x = "Selenoneine (µg/g)", y = "Ergothioneine (µg/g)")+
  ggtitle(expression(bold("d")))+
  theme_1+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())
ERGO_SEN_brain


corr_ERGO_SEN_int <- corr_ERGO_SEN %>% filter( tissue == "Intestine") 
lm <- lm(ergo ~ seleno, data = corr_ERGO_SEN_int)
summary(lm)
layout(matrix(c(1:5, 5), 2, 3)) +
  plot(lm, c(1, 3:5), pch = 20, cex = 1.5) 
  qqnorm(lm$residuals)  
  layout(1)


#correlation de pearson
N_int <- nrow(corr_ERGO_SEN_int)
cor.test(corr_ERGO_SEN_int$seleno, corr_ERGO_SEN_int$ergo)
ERGO_SEN_int <- ggplot(corr_ERGO_SEN_int , mapping = aes(y= ergo, x = seleno))+
  geom_point( shape = 21, fill = 'black', color = 'white', size = 5)+
  sm_statCorr(color = "grey", text_size	= 10)+
  facet_wrap(~tissue,
             labeller = labeller(tissue = c("Intestine" = "Intestine")))+
  annotate("text", x = min(corr_ERGO_SEN_int$seleno), 
           y = max(corr_ERGO_SEN_int$ergo), 
           label = paste0("n = ", N_int),
           hjust = 0, vjust = 3, size = 10, color = "black") + 
  labs(x = "Selenoneine (µg/g)", y = "Ergothioneine (µg/g)")+
  ggtitle(expression(bold("e")))+
  theme_1+
  theme(axis.title.x = element_blank())
ERGO_SEN_int

corr_ERGO_SEN_muscle <- corr_ERGO_SEN %>% filter( tissue == "Muscle") 
lm <- lm(ergo ~ seleno, data = corr_ERGO_SEN_muscle)
summary(lm)
layout(matrix(c(1:5, 5), 2, 3)) +
  plot(lm, c(1, 3:5), pch = 20, cex = 1.5) 
  qqnorm(lm$residuals)  
  layout(1)

#correlation de pearson
N_muscle <- nrow(corr_ERGO_SEN_muscle)
cor.test(corr_ERGO_SEN_muscle$seleno, corr_ERGO_SEN_muscle$ergo)
ERGO_SEN_muscle <- ggplot(corr_ERGO_SEN_muscle , mapping = aes(y= ergo, x = seleno))+
  geom_point( shape = 21, fill = 'black', color = 'white', size = 5)+
  sm_statCorr(color = "grey", text_size	= 10)+
  facet_wrap(~tissue,
             labeller = labeller(tissue = c("Muscle" = "Muscle")))+
  labs(x = "Selenoneine (µg/g)", y = "Ergothioneine (µg/g)")+
  annotate("text", x = min(corr_ERGO_SEN_muscle$seleno), 
           y = max(corr_ERGO_SEN_muscle$ergo), 
           label = paste0("n = ", N_muscle),
           hjust = 0, vjust = 3, size = 10, color = "black") + 
  ggtitle(expression(bold("f")))+
  theme_1+
  theme(axis.title.y = element_blank())
ERGO_SEN_muscle

corr_ERGO_SEN_liver <- corr_ERGO_SEN %>% filter( tissue == "Liver") 
lm <- lm(ergo ~ seleno, data = corr_ERGO_SEN_liver)
summary(lm)
layout(matrix(c(1:5, 5), 2, 3)) +
  plot(lm, c(1, 3:5), pch = 20, cex = 1.5) 
  qqnorm(lm$residuals)  
  layout(1)


#correlation de pearson
N_liver <- nrow(corr_ERGO_SEN_liver)
cor.test(corr_ERGO_SEN_liver$seleno, corr_ERGO_SEN_liver$ergo)
ERGO_SEN_liver <- ggplot(corr_ERGO_SEN_liver , mapping = aes(y= ergo, x = seleno))+
  geom_point( shape = 21, fill = 'black', color = 'white', size = 5)+
  sm_statCorr(color = "grey", text_size	= 10 )+
  facet_wrap(~tissue,
             labeller = labeller(tissue = c("Liver" = "Liver")))+
  labs(x = "Selenoneine (µg/g)", y = "Ergothioneine (µg/g)")+
  annotate("text", x = min(corr_ERGO_SEN_liver$seleno), 
           y = max(corr_ERGO_SEN_liver$ergo), 
           label = paste0("n = ", N_liver),
           hjust = 0, vjust = 3, size = 10, color = "black") + 
  ggtitle(expression(bold("c")))+
  theme_1+
  theme(axis.title.x = element_blank())
ERGO_SEN_liver

corr_ERGO_SEN_mattaaq <- corr_ERGO_SEN %>% filter( tissue == "Skin") 
lm <- lm(ergo ~ seleno, data = corr_ERGO_SEN_mattaaq)
summary(lm)
layout(matrix(c(1:5, 5), 2, 3)) +
  plot(lm, c(1, 3:5), pch = 20, cex = 1.5) +
  qqnorm(lm$residuals)  +
  layout(1)

#correlation de pearson
N_mattaaq <- nrow(corr_ERGO_SEN_mattaaq)
cor.test(corr_ERGO_SEN_mattaaq$seleno, corr_ERGO_SEN_mattaaq$ergo)
ERGO_SEN_mattaaq <- ggplot(corr_ERGO_SEN_mattaaq , mapping = aes(y= ergo, x = seleno))+
  geom_point( shape = 21, fill = 'black', color = 'white', size = 5)+
  sm_statCorr(color = "grey", text_size	= 10)+
  facet_wrap(~tissue,
             labeller = labeller(tissue = c("Skin" = "Skin")))+
  labs(x = "Selenoneine (µg/g)", y = "Ergothioneine (µg/g)")+
  annotate("text", x = min(corr_ERGO_SEN_mattaaq$seleno), 
           y = max(corr_ERGO_SEN_mattaaq$ergo), 
           label = paste0("n = ", N_mattaaq),
           hjust = 0, vjust = 3, size = 10, color = "black") + 
  ggtitle(expression(bold("a")))+
  theme_1+
  theme(axis.title.x = element_blank())
ERGO_SEN_mattaaq

corr_ERGO_SEN_blood <- corr_ERGO_SEN %>% filter( tissue == "Blood") 
lm <- lm(ergo ~ seleno, data = corr_ERGO_SEN_blood)
summary(lm)
layout(matrix(c(1:5, 5), 2, 3)) +
  plot(lm, c(1, 3:5), pch = 20, cex = 1.5) +
  qqnorm(lm$residuals)  +
  layout(1)


#correlation de pearson
N_blood <- nrow(corr_ERGO_SEN_blood)
cor.test(corr_ERGO_SEN_blood$seleno, corr_ERGO_SEN_blood$ergo)
ERGO_SEN_blood <- ggplot(corr_ERGO_SEN_blood , mapping = aes(y= ergo, x = seleno))+
  geom_point( shape = 21, fill = 'black', color = 'white', size = 5)+
  sm_statCorr(color = "grey", text_size	= 10)+
  facet_wrap(~tissue,
             labeller = labeller(tissue = c("Blood" = "Whole Blood")))+
  labs(x = "Selenoneine (µg/g)", y = "Ergothioneine (µg/g)")+
  annotate("text", x = min(corr_ERGO_SEN_blood$seleno), 
           y = max(corr_ERGO_SEN_blood$ergo), 
           label = paste0("n = ", N_blood),
           hjust = 0, vjust = 3, size = 10, color = "black") + 
  ggtitle(expression(bold("g")))+
  theme_1
ERGO_SEN_blood


#figure 4
figure3 <- grid.arrange(ERGO_SEN_mattaaq, ERGO_SEN_rein, ERGO_SEN_liver, ERGO_SEN_brain,
                        ERGO_SEN_int, ERGO_SEN_muscle, ERGO_SEN_blood,
                        ncol =2)
figure3
#save as 20x20
#save as PDF 20X25
ggsave(figure4, 
       filename = "figure3.jpg",
       device = "jpg",
       height = 25, width = 25, units = "in",
       path = "downloads")



# Figure 4, Descriptive data on skin  -------------------------------------
#naame SEN_EGT_skin_leyers_beluga_ABB
beluga_skin <- read_delim(file.choose(), delim = ";", trim_ws = TRUE, show_col_types = FALSE)
donnee_descriptive_peau <- beluga_skin  %>% drop_na(Concentration) %>%
  setDT()
donnee_descriptive_peau <- donnee_descriptive_peau[, .(mean = mean(Concentration), geom_mean = geomMean(Concentration), 
                                                       min = min(Concentration),
                                                       max = max(Concentration),
                                                       median = median(Concentration),
                                                       percentile25 = quantile(Concentration, probs = 0.25),
                                                       percentile75 = quantile(Concentration, probs = 0.75),
                                                       geom_se = geomSE(Concentration),
                                                       se = std_err(Concentration, .N), n = .N ),
                                                   by=.(Skin_Segment, Analyte)]
gt(donnee_descriptive_peau)


# Supp. Table 3 Foetus concentrations -------------------------------------
#name data : SEN_EGT_beluga_foetus_tissues_ABB
beluga_foetus <- read_delim(file.choose(), delim = ";", trim_ws = TRUE, show_col_types = FALSE)
donnee_descriptive_foetus <- beluga_foetus  %>% drop_na(Concentration) %>%
  setDT()
donnee_descriptive_foetus <- donnee_descriptive_foetus[, .(mean = mean(Concentration), geom_mean = geomMean(Concentration), 
                                                       min = min(Concentration),
                                                       max = max(Concentration),
                                                       median = median(Concentration),
                                                       percentile25 = quantile(Concentration, probs = 0.25),
                                                       percentile75 = quantile(Concentration, probs = 0.75),
                                                       geom_se = geomSE(Concentration),
                                                       se = std_err(Concentration, .N), n = .N ),
                                                   by=.(Tissue, Analyte)]
gt(donnee_descriptive_foetus)

