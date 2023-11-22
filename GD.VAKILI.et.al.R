#Goldeneye long term anti-predator behaviour data analyses: 

GD <- GD1984_2021_Vakili_etal

# 0) Correlations: 
install.packages("Hmisc")

GD.COR <- as.data.frame(lapply(GD[, c("ClutchSize", "hatching.success", "Age", "Body.condition", "ParasitizedNest")], as.numeric))
COR <- cor(GD.COR, method = "spearman", use = "complete.obs")
summary(COR)
## Initialize empty matrices to store p-values and correlation coefficients
COR.P <- matrix(NA, nrow = ncol(COR), ncol = ncol(COR))
COR.coefficients <- matrix(NA, nrow = ncol(COR), ncol = ncol(COR))

COR <- rcorr(as.matrix(GD.COR), type = "spearman")
overall_Rs <- COR$r[1,2]
p_value <- COR$P[1,2]
cat("Overall Spearman correlation coefficient (Rs):", overall_Rs, "\n")
cat("P-value:", p_value, "\n")

# Loop through each variable pair to calculate p-values and correlation coefficients
for (i in 1:ncol(GD.COR)) {
  for (j in 1:ncol(GD.COR)) {
    if (i != j) {
      cor_test <- cor.test(GD.COR[, i], GD.COR[, j], method = "spearman", na.action = "na.omit")
      COR.P[i, j] <- cor_test$p.value
      COR.coefficients[i, j] <- cor_test$estimate
    }
  }
}

print(COR.P)
print(COR.coefficients)


# 1) Section 2. "Materials and mthods"
## Age identification of individuals
### Plotting the number of individuals per behavioural levels in each age group: 

ggplot(GD, aes(x=Age, fill=Behavior))+
  geom_bar(position="stack")+
  labs(x= "Age", y="Number of individuals")+
  scale_x_continuous(breaks=1:17)+
  scale_fill_manual(values=c("turquoise", "lightsteelblue2", "royalblue4"),
                    labels=c("Shy", "Average", "Bold"))+ 
  theme(
    panel.background=element_blank(),
    plot.background = element_blank(),
    axis.line = element_line("black"),
    panel.grid.major= element_blank(),
    panel.grid.minor= element_blank(),
    plot.title= element_text(family = "Times New Roman", size= 15))
#2) Section 3. "Results" 
## Population information at the beginning of result + table 1: 
length(unique(GD$FemaleID))

range(GD$Age)
Age.mean <- mean(GD$Age)
Age.sd <- sd(GD$Age)/ sqrt(length(GD$Age))
print(paste("Mean:",Age.mean))
print(paste("SE:",Age.sd))

range(GD$Wing, na.rm=TRUE)
mean.wing <- mean(GD$Wing, na.rm= TRUE)
sd.wing <- sd(GD$Wing, na.rm= TRUE)/ sqrt(length(GD$Wing))
print(paste("Mean:", mean.wing))
print(paste("SE:", sd.wing))

range(GD$Weight, na.rm=TRUE)
mean.weight <- mean(GD$Weight, na.rm= TRUE)
sd.weight <- sd(GD$Weight, na.rm= TRUE)/ sqrt(length(GD$Weight))
print(paste("Mean:", mean.weight))
print(paste("SE:", sd.weight))

range(GD$NestingNumber)
mean.nesting <- mean(GD$NestingNumber)
sd.nesting <- sd(GD$NestingNumber)/ sqrt(length(GD$NestingNumber))
print(paste("Mean:", mean.nesting))
print(paste("SE:", sd.nesting))

Age.percentage <- sum(GD$Age >1 & GD$Age <4) / length(GD$Age) *100
Age.percentage

range(GD$Body.condition, na.rm=TRUE)
body.con.mean <- mean(GD$Body.condition, na.rm=TRUE)
print(paste("Mean:", body.con.mean))
body.con.sd <- sd(GD$Body.condition, na.rm=TRUE)/ sqrt(sum(!is.na(GD$Body.condition)))
body.con.sd

Clutchsize.mean <- mean(GD$ClutchSize)
Clutchsize.mean
Clutchsize.sd <- sd(GD$ClutchSize)/ sqrt(length(GD$ClutchSize))
print(paste("SE:",Clutchsize.sd))
Clutchsize.percentage <- sum(GD$ClutchSize >12) / length(GD$ClutchSize) *100
Clutchsize.percentage
sum(GD$ClutchSize> 12)

more.than.one.breeding.attempt <- GD %>%
  filter(NestingNumber > 1 & NestingNumber < 17) %>%
  distinct(FemaleID)
count_more.than.one.breeding.attempt <- nrow(more.than.one.breeding.attempt)
count_more.than.one.breeding.attempt

# 2.1) Initial distribution of behavioural types among the population, below table 1:
Initial_Behavioural_type <- table(GD$Behavior)
print(Initial_Behavioural_type)

Average_count <- 925
bold_count <- 397
contingency_table <- matrix(c(Average_count, bold_count), nrow=2)
chi_square <- chisq.test(contingency_table)
phenotype.comparison <- table(GD$Behavior %in% c("Bold", "Average"))
chi_square1 <- chisq.test(phenotype.comparison)
chi_square1

print(chi_square)

#########################################################################################################

#3) Section 3.1 "Results"
## Relationship between behavior and age of the breeding females:
library(emmeans)
library(ordinal)
library(lsmeans)

GD$Behavior<- as.factor(GD$Behavior)
GD$Behavior <- ordered(GD$Behavior, levels = c("Shy", "Average", "Bold"))

M.Behavior <- clmm(Behavior ~ Age * ClutchSize + Body.condition + hatching.success + ParasitizedNest + I(ClutchSize^2) + (1 | FemaleID), data = GD)
summary(M.Behavior)
# 3.1)Testing the model fitness: 
## presented in table 2
Null.clmm <- clmm(Behavior~ 1+ (1|FemaleID), data=GD)
AIC_null <- AIC(Null.clmm)
AIC_M.Behavior <- AIC(M.Behavior)
delta_AIC <- AIC_null - AIC_M.Behavior
print(delta_AIC)

loglikeNull <- logLik(Null.clmm)
loglikeModel <- logLik(M.Behavior)
chisq <- -2 * (loglikeNull - loglikeModel)
df <- df.residual(M.Behavior)-df.residual(Null.clmm)
p_value <- 1- pchisq(chisq, df)

LRT_statistic <- -2 * (logLik(Null.clmm) - logLik(M.Behavior))
df <- attr(LRT_statistic, "df")
p_value <- pchisq(LRT_statistic, df, lower.tail = FALSE)
cat("Chi-square statistics:", LRT_statistic, "\n")
cat("Degrees of freedom:", df,"\n")
cat("p-value associated with chi-square:", p_value, "\n")

# 3.2) Emmeans: Behaviour- age, table 3: 
PA2<-data.frame(summary(emmeans(M.Behavior, pairwise~ Behavior | Age,at = list(Age = 2), mode = "prob")))
PA3<-data.frame(summary(emmeans(M.Behavior, pairwise~ Behavior | Age,at = list(Age = 3), mode = "prob")))
PA4<-data.frame(summary(emmeans(M.Behavior, pairwise~ Behavior | Age,at = list(Age = 4), mode = "prob")))
PA5<-data.frame(summary(emmeans(M.Behavior, pairwise~ Behavior | Age,at = list(Age = 5), mode = "prob")))
PA6<-data.frame(summary(emmeans(M.Behavior, pairwise~ Behavior | Age,at = list(Age = 6), mode = "prob")))
PA7<-data.frame(summary(emmeans(M.Behavior, pairwise~ Behavior | Age,at = list(Age = 7), mode = "prob")))
PA8<-data.frame(summary(emmeans(M.Behavior, pairwise~ Behavior | Age,at = list(Age = 8), mode = "prob")))
PA9<-data.frame(summary(emmeans(M.Behavior, pairwise~ Behavior | Age,at = list(Age = 9), mode = "prob")))
PA10<-data.frame(summary(emmeans(M.Behavior, pairwise~ Behavior | Age,at = list(Age = 10), mode = "prob")))
PA11<-data.frame(summary(emmeans(M.Behavior, pairwise~ Behavior | Age,at = list(Age = 11), mode = "prob")))
PA12<-data.frame(summary(emmeans(M.Behavior, pairwise~ Behavior | Age,at = list(Age = 12), mode = "prob")))
PA13<-data.frame(summary(emmeans(M.Behavior, pairwise~ Behavior | Age,at = list(Age = 13), mode = "prob")))
PA14<-data.frame(summary(emmeans(M.Behavior, pairwise~ Behavior | Age,at = list(Age = 14), mode = "prob")))
PA15<-data.frame(summary(emmeans(M.Behavior, pairwise~ Behavior | Age,at = list(Age = 15), mode = "prob")))
PA16<-data.frame(summary(emmeans(M.Behavior, pairwise~ Behavior | Age,at = list(Age = 16), mode = "prob")))
PA17<-data.frame(summary(emmeans(M.Behavior, pairwise~ Behavior | Age,at = list(Age = 17), mode = "prob")))

# 3.3) Plot: Behaviour and age, Fig 3A: 
library(ggplot2)

PA_plot<-rbind(PA2,PA3,PA4,PA5,PA6,PA7,PA8,PA9,PA10,PA11,PA12,PA13,PA14,PA15,PA16,PA17)
names(PA_plot)
PA_plot<- PA_plot[, c(1:7)]

head(PA_plot)
names(PA_plot) <- c("Behavior","Age","probability", "SE","df","LCI","UCI")

PA_plot$Behavior <- ordered(PA_plot$Behavior, Behavior= c("Shy", "Average", "Bold"))


M1.P <- ggplot()+
  geom_smooth(data = PA_plot,
              aes(x = Age,
                  y = probability,
                  color = Behavior,
                  group = Behavior),
              size=1.5, se = F)+
  geom_ribbon(data = PA_plot,
              aes(x = Age,
                  ymax = probability + SE,
                  ymin = probability - SE ,
                  fill = Behavior,
                  group = Behavior),
              alpha = 0.3)

M1.P  <- M1.P  + xlab("Age") + ylab("probability to exhibit the Behavior") + ggtitle("(A)")
M1.P <- M1.P  + theme(text = element_text(size=16), panel.background = element_rect(fill = "transparent",colour = NA),
                      plot.background = element_rect(fill = "transparent",colour = NA),legend.position = c(0.9, 0.7),
                      plot.title= element_text(family = "Times New Roman", size= 16), axis.line = element_line("black"))
M1.P  <- M1.P  + scale_color_manual(values=c("turquoise", "lightsteelblue2", "royalblue4"),labels=c("Shy", "Average", "Bold"))+
  scale_fill_manual(values=c("turquoise", "lightsteelblue2", "royalblue4"),
                    labels=c("Shy", "Average", "Bold"))

M1.P

#4) Section 3.2 "Results"
## Relationship between behaviour and clutch size and nest parasitism: 
M.Behavior <- clmm(Behavior ~ Age * ClutchSize + Body.condition + hatching.success + ParasitizedNest + I(ClutchSize^2) + (1 | FemaleID), data = GD)
summary(M.Behavior)

# 4.1) Emmeans: Behavior and clutch size, table 4: 
P1 <- data.frame(summary(emmeans(M.Behavior, pairwise~ Behavior | ClutchSize,at = list(ClutchSize = 2), mode = "prob"))) 
P2 <- data.frame(summary(emmeans(M.Behavior, pairwise~ Behavior | ClutchSize,at = list(ClutchSize = 8), mode = "prob"))) 
P3 <- data.frame(summary(emmeans(M.Behavior, pairwise~ Behavior | ClutchSize,at = list(ClutchSize = 13), mode = "prob")))

# 4.2) Plot: Behaviour and clutch size, Fig 3B: 

PA_plot<-rbind(P1,P2,P3)
names(PA_plot)
PA_plot<- PA_plot[, c(1:7)]

head(PA_plot)
names(PA_plot) <- c("Behavior","ClutchSize","probability", "SE","df","LCI","UCI")

M1.P <- ggplot()+
  geom_smooth(data = PA_plot,
              aes(x = ClutchSize,
                  y = probability,
                  color = Behavior,
                  group = Behavior),
              size=1.5, se = F)+
  geom_ribbon(data = PA_plot,
              aes( x= ClutchSize,
                   ymax = probability + SE,
                   ymin = probability - SE,
                   fill = Behavior,
                   group = Behavior),
              alpha = 0.3)

M1.P  <- M1.P  + xlab("Clutch size") + ylab("Probability to exhibit the behavior") + ggtitle("(B)")
M1.P <- M1.P  + theme(text = element_text(size=15), panel.background = element_rect(fill = "transparent",colour = NA),
                      plot.background = element_rect(fill = "transparent",colour = NA), legend.position = c(0.9, 0.9),
                      plot.title= element_text(family = "Times New Roman", size= 15), axis.line = element_line("black"))
M1.P  <- M1.P  + scale_color_manual(values=c("turquoise", "lightsteelblue2", "royalblue4"),labels=c("Shy", "Average","Bold"))+
  scale_fill_manual(values=c("turquoise", "lightsteelblue2", "royalblue4"),
                    labels=c("Shy", "Average","Bold"))
M1.P 

#5) Section 3.3. "Results"
## Relationship between behavior and interaction effect of age and clutch size:

M.Behavior <- clmm(Behavior ~ Age * ClutchSize + Body.condition + hatching.success + ParasitizedNest + I(ClutchSize^2) + (1 | FemaleID), data = GD)
summary(M.Behavior)

## 5.1) Emmeans, table 5: 
## for each coparison, replace the "Age" and "Clutchszie" values with the ones presented in table 5:

P1 <- data.frame(summary(emmeans(M.Behavior, pairwise~ Behavior | Age*ClutchSize,at = list(Age = 15, ClutchSize =19), mode = "prob")))

# 5.2) Plotting the interaction between age and clutch size in association with behavior, Fig 4:
library(effects)
library(MASS)

behavior_levels <- ordered(GD$Behavior, levels=c("Shy", "Average", "Bold"))

dev.off()
Effect <- effect("Age:ClutchSize", M.Behavior)
behavior_colors <- c("turquoise", "lightsteelblue2", "royalblue4")
colorsGD <- behavior_colors[behavior_levels]

par(family="Times New Roman")
plot(Effect, main="", lines=TRUE,
     ylab="Behavior (probability of occurance)", colours= colorsGD,
     cex.axis=2, cex.lab=3, cex.main=2)

#6) Section 3.4. "Results"
## Relationship between hatching success and behavior:  
###Since there are a lot of 1s and some 0s in the datset, betaregression can not be the suitable method to use! 
#### The alternative approach turns the hatching success into a binary variable = 0 failure and 1 success:

library(lme4)

M2.Behavior <- glmer(hatching.success ~ Behavior + Age + ClutchSize+ Body.condition+ ParasitizedNest +(1|FemaleID), data= GD, family=binomial)
summary(M2.Behavior)

# 6.1) To check the between-levels effect in behavior:
GD$Behavior <- factor(GD$Behavior, ordered= FALSE)
GD$Behavior <- relevel(GD$Behavior, ref= 2) 

# 6.2) Goodness of fit:
## presented in table 2: 

library(DHARMa)
simulationOutput <- simulateResiduals(fittedModel= M2.Behavior, plot= T) 

null_model <- glmer(hatching.success ~ 1 + (1|FemaleID), data = GD, family = binomial)

AIC_null <- AIC(null_model)
AIC_M3.Behavior <- AIC(M2.Behavior)
delta_AIC <- AIC_null - AIC_M3.Behavior
print(delta_AIC)

loglikeNull <- logLik(null_model)
loglikeModel <- logLik(M2.Behavior)
chisq <- -2 * (loglikeNull - loglikeModel)
df <- df.residual(M2.Behavior)-df.residual(null_model)
p_value <- 1- pchisq(chisq, df)

LRT_statistic <- -2 * (logLik(null_model) - logLik(M2.Behavior))
df <- attr(LRT_statistic, "df")
p_value <- pchisq(LRT_statistic, df, lower.tail = FALSE)
cat("Chi-square statistics:", LRT_statistic, "\n")
cat("Degrees of freedom:", df,"\n")
cat("p-value associated with chi-square:", p_value, "\n")

## 6.3) Plotting, Fig 5: 
library(ggplot2)
library(effects)
dev.off()

GD$Behavior <- ordered(GD$Behavior, levels = c("Shy", "Average", "Bold"))
effect2 <- as.data.frame(effect(M2.Behavior, term= "Behavior")) 
effect2_df <- data.frame(Behavior= effect2$Behavior, Probability= effect2$fit)
effect2_df$SE <- effect2$se
effect2_df$CI_lower <- effect2$lower
effect2_df$CI_upper <- effect2$upper
behavior_colors <- c("turquoise", "lightsteelblue2", "royalblue4")
P <- ggplot(effect2_df, aes(x= Behavior, y= Probability, fill= Behavior))+
  geom_bar(stat="identity", width= 0.5)+
  geom_errorbar(aes(ymin=CI_lower, ymax = CI_upper), width = 0.2)+
  scale_fill_manual(values=behavior_colors)+
  labs(x= "Behavior", y="Hatching success probability")
theme_minimal() x
P <- P + theme(legend.position="none")
P <- P + theme(
  panel.background=element_blank(),
  plot.background = element_blank(),
  axis.line = element_line("black"),
  panel.grid.major= element_blank(),
  panel.grid.minor= element_blank()
)
P

#################### Created by :
#Farshad S. Vakili; 2023, Turku, Finland 
#E-mail: Farshad.S.Vakili@gmail.com 
################
     


