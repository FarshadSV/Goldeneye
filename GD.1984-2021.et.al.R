#Goldeneye long term anti-predator behaviour data analyses: 

library(Hmisc)
library(ggplot2)
library(ggeffects)
library(lsmeans)
library(emmeans)
library(ordinal)
library(lme4)
library(DHARMa)
library(effects)
library(marginaleffects)
library(MASS)
library(dplyr)
library(rptR)
library(tidyr)
library(boot)
library(parallel)


#Uploading the dataset
GD <- GD1984_2021

# 0) Correlations:
GD.COR <- as.data.frame(lapply(GD[, c("ClutchSize", "hatchSuccess","hatchFailure", "Age", "Body.condition","ParasitizedNest")], as.numeric))
GD.COR <- as.data.frame(lapply(GD[, c("ClutchSize", "Age", "Body.condition")], as.numeric))
COR <- cor(GD.COR, method = "spearman", use = "complete.obs")
summary(COR)
## Initialize empty matrices to store p-values and correlation coefficients
COR.P <- matrix(NA, nrow = ncol(COR), ncol = ncol(COR))
COR.coefficients <- matrix(NA, nrow = ncol(COR), ncol = ncol(COR))

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
  labs(x= "Age", y="Number of individuals", fill="Behaviour")+
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

### ALSO, Plotting the number of individuals per behavioural levels in each Clutch size group: 
ggplot(GD, aes(x=ClutchSize, fill=Behavior))+
  geom_bar(position="stack")+
  labs(x= "Clutch size", y="Number of individuals")+
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

#########################################################################################################

#3) Section 3.1 "Results"
## Relationship between behavior and age of the breeding females:

GD$Behavior<- as.factor(GD$Behavior)
GD$Behavior <- ordered(GD$Behavior, levels = c("Shy", "Average", "Bold"))

M.Behavior <- clmm(Behavior ~ Age * ClutchSize + Body.condition + I(ClutchSize^2) + (1| FemaleID), data = GD)
summary(M.Behavior)


### Supplementary material: Testing M.Behavior fit with alternative models even though the only model which meets our prior hypotheses is M.Behavior!
M.Behavior2 <- clmm(Behavior ~ Age + ClutchSize +  Body.condition+ (1| FemaleID), data = GD)
M.Behavior3 <- clmm(Behavior ~ Age * ClutchSize +  Body.condition+ (1| FemaleID), data = GD)

  AIC_BIC <- data.frame(
    Model = c("M.Behavior", "M.Behavior2", "M.Behavior3","Null.clmm"),
    AIC = c(AIC(M.Behavior), AIC(M.Behavior2), AIC(M.Behavior3), AIC(Null.clmm)),
    BIC = c(BIC(M.Behavior), BIC(M.Behavior2), BIC(M.Behavior3), BIC(Null.clmm)),
    logLik = c(logLik(M.Behavior), logLik(M.Behavior2), logLik(M.Behavior3), logLik(Null.clmm))
  )
print(AIC_BIC)

#3.1)Testing the model fitness, compared with null model (presented in table 2): 
# 3.1.1) Fitting the null model
#creating a common dataset without NA. Since the M.Behavior model excludes NAs:
common_GD <- na.omit(GD)
M.Behavior <- clmm(Behavior ~ Age * ClutchSize + Body.condition + I(ClutchSize^2) + (1| FemaleID), data = common_GD)
Null.clmm <- clmm(Behavior~ 1+ (1|FemaleID), data=common_GD)

# 3.1.2) Likehood Rtio Test (LRT)
LRTbehavior <- anova(M.Behavior, Null.clmm)
LRTbehavior

# 3.2) Given the proportion of behavioral categories per age (Fig.2), we investigate the interpretation of the output of M.Behavior through several stages as follows:
#3.2.1) Checking the proportional change in the number of individuals per categories with age:
Age_compare <- subset(GD, Age %in% c(3, 8) & GD$Behavior %in% c("Shy", "Average", "Bold"))
contingency_table <- table(Age_compare$Age, Age_compare$Behavior)
print(contingency_table)

Shy_3 <- contingency_table[1, "Shy"]
Shy_8 <- contingency_table[2, "Shy"]
Bold_3 <- contingency_table[1, "Bold"]
Bold_8 <- contingency_table[2, "Bold"]
Average_3 <- contingency_table[1, "Average"]
Average_8 <- contingency_table[2, "Average"]
# Calculating proportions
prop_bold_3 <- Bold_3 / (Bold_3 + Shy_3+ Average_3)
prop_bold_8 <- Bold_8 / (Bold_8 + Shy_8 + Average_8)
prop_shy_3 <- Shy_3 / (Bold_3 + Shy_3 + Average_3)
prop_shy_8 <- Shy_8 / (Bold_8 + Shy_8 + Average_8)
prop_average_3 <- Average_3 / (Bold_3+ Shy_3 + Average_3)
prop_average_8 <- Average_8 / (Bold_8 + Shy_8 + Average_8)
#difference between behavioral proportions
chis_test <- chisq.test(contingency_table)
print(chis_test) # Significancy of the changes in proportions from age 3 to 8 

prop_sig <- prop.test(c(Bold_3, Bold_8), c(Shy_3 + Bold_3 + Average_3, Shy_8 + Bold_8 + Average_8), correct = FALSE)
print(prop_sig)

#3.2.2) Checking the proportional change in the number of individuals per categories with clutch size:
Clutch_compare <- subset(GD, ClutchSize %in% c(8, 12) & GD$Behavior %in% c("Shy","Average", "Bold"))
contingency_table <- table(Clutch_compare$ClutchSize, Clutch_compare$Behavior)
print(contingency_table)

Shy_8 <- contingency_table[1, "Shy"]
Shy_12 <- contingency_table[2, "Shy"]
Bold_8 <- contingency_table[1, "Bold"]
Bold_12 <- contingency_table[2, "Shy"]
Average_8 <- contingency_table[1, "Average"]
Average_12 <- contingency_table[2, "Average"]
# Calculating proportions
prop_bold_8 <- Bold_8 / (Bold_8 + Shy_8 + Average_8)
prop_bold_12 <- Bold_12 / (Bold_12 + Shy_12 + Average_12)
prop_shy_8 <- Shy_8 / (Bold_8 + Shy_8 + Average_8)
prop_shy_12 <- Shy_12 / (Bold_12 + Shy_12 + Average_12)
prop_average_8 <- Average_8 / (Bold_8+ Shy_8 + Average_8)
prop_average_12 <- Average_12 / (Bold_12 + Shy_12 + Average_12)
#difference between behavioral proportions
chis_test <- chisq.test(contingency_table)
print(chis_test) # Whether the change in the proportion of clutch size among behavioral categories is significant or not 

prop_sig <- prop.test(c(Shy_6, Shy_12), c(Shy_6 + Bold_6 + Average_12, Shy_12 + Bold_12 + Average_12), correct = FALSE)
print(prop_sig)


# 3.3) Assessment of repeatability estimates and behavioral variance
#First: Creating a sub-set of data, including only individuals with more than one breeding attempt 
GD_Sub <- GD %>% 
  group_by(FemaleID) %>%
  filter(n() >1) %>%
  ungroup()

#Then to test how much of total variance is attributed to between-individual differences in behavior. 
#Running the model with the subset of data
M.BehaviorSUB <- clmm(Behavior ~ Age * ClutchSize + Body.condition + I(ClutchSize^2) + (1| FemaleID), data = GD_Sub)
summary(M.BehaviorSUB)
#We have two main choices for repeatability testing. Either converting the Behavior to a numerical variable and run a lmm as below:
GD_Sub <- GD_Sub %>%
  mutate(Behavior_numeric = as.integer(Behavior))
## 3.3.1) Among individual behavioral consistency, via mixed-effect model method, "lme4" package,
# Fitting a linear mixed-effects model to calculate the variance of random effect "FemaleID", while weighting for the fixed effects
##representing the variance of the random intercept for each individual, which represents the within-individual variance
M.Linear <- lmer(Behavior_numeric ~ Age * ClutchSize + Body.condition + ParasitizedNest + I(ClutchSize^2) + (1 | FemaleID), data = GD_Sub, REM=TRUE)
# Extracting the variance component of the model 
var_components <- VarCorr(M.Linear)
var_components
# Extracting the variance of the random intercept for each individual
var_random <- var_components$FemaleID[1]
var_random 
#Residual variance 
residual_var <- attr(VarCorr(M.Linear), "sc")^2
residual_var
# Calculating the total variance (sum of residual variance (variability not explained by fixed or random effect) + variance component for each female)
total_variance <- residual_var + var_random
# Calculating the Coefficient of relavtive plasticity (CRP)/ICC as the proportion of variance in Behavior that can be attributed to differences among individuals
CRP <- var_random / total_variance
# Printing the CRP (Showing the ratio of variance component of the random effect for each individual in behavior to the overall behavioral variance in population)
print(CRP) #CRP = Ï^2_individual / (Ï^2_individual + Ï^2_residual) =X% of total variance in Behavior is due of among-individual difference 


#BUT, 3.3.2) SINCE Behavior is a categorical ordinal variable, it's best to use the "ICC/CRP" method without converting the data into numeric
### Extracting the variance for the random effect (FemaleID) and the residual variance
var_random <- VarCorr(M.BehaviorSUB)$FemaleID[1]
residual_var <- pi^2 / 3 # fixed for logit link in clmm modelig 
ICC <- var_random / (var_random + residual_var)
cat("Intra-class Correlation Coefficient (ICC):", ICC, "\n")

#3.3.3) bootstrap method: calculating repeatability by generating  a distribution of estimates by re-sampling the data with replacement, helping to assess the consistency and variability of the estimates
#Creating 1000 bootstrap samples (runs 1000 times) from "GD_Sub", each a random sample of the original data with replacement
#bootsrap_repeatability calculates the repeatability for each of these samples by fitting a clmm and extracting the variance components

fit_clmm <- function(data) {
  tryCatch({
    clmm(Behavior ~ Age * ClutchSize + Body.condition + I(ClutchSize^2) + (1| FemaleID), 
         data = data)
  }, error = function(e) return(NULL))
}

bootstrap_repeatability <- function(data, indices) {
  data_boot <- data[indices, ]
  model_boot <- fit_clmm(data_boot)
  if (is.null(model_boot)) {
    return(NA)  # Return NA if the model did not converge
  } else {
    var_comp_boot <- VarCorr(model_boot)
    var_random <- var_comp_boot$FemaleID[1] # betwee-individual variance: estimating the variance component of "FemaleID"
    residual_var = pi^2 / 3   # Fixed for logit link in CLMM
    repeatability = var_random / (var_random + residual_var)
    return(repeatability)
  }
}

# Running the bootstrapping with parallel processing
set.seed(123)
boot_results <- boot(data = GD_Sub, statistic = bootstrap_repeatability, R = 1000)
print(boot_results) 
# The output is giving the proportion of total variance in Behavior that is attributed to differences among individuals (between-individual consistency)
# In this context, it's suggesting how individuals maintain a consistent distinct behavior over time (could also be an indication of within-subject consistency. However, it does not directly calculates this one.)
# 95% confidence intervals for the bootstrap results
conf_intervals <- boot.ci(boot_results, type = "perc", conf = 0.95)
print(conf_intervals)
# Bias in bootstrap is the average difference between the bootstrap estimate and the original estimate. Positive bias suggests that the original data might have underestimated the "true repeatability", maybe due to outliers or data characteristics 
# Bootstrap original repeatability estimate (without re-sampling from data) matches with ICC coefficient, confirming consistent behavior among individuals relative to the total variance 
## Based on bootstrap and ICC: A significant portion of behavior (variability in behavior) is consistent across individuals (suggesting existence of personality trait in goldeneye)


# 3.4) To further investigate the interpretation of results:
#calculating the within and between-individual coefficients of age and clutchsize effect (significant in the main model) on behaviour 
# Mean-centering approach 
#Calculating among-individual components (individual means). Intercept of the reaction norm 
individual_means <- GD_Sub %>%
  group_by(FemaleID) %>%
  summarize(mean_Age = mean(Age), 
            mean_ClutchSize = mean(ClutchSize),
            .groups = 'drop')

GD_Sub <- GD_Sub %>%
  left_join(individual_means, by = "FemaleID")
# Calculating the within-individual components (deviation from individual means). The residual. 
GD_Sub <- GD_Sub %>% 
  mutate(
    Age_within = Age - mean_Age,  # Within-individual component for Age
    ClutchSize_within = ClutchSize - mean_ClutchSize  # Within-individual component for ClutchSize
  )

# Fitting the CLMM with both within- and among-individual components
GD_Sub$Behavior<- as.factor(GD_Sub$Behavior)
GD_Sub$Behavior <- ordered(GD_Sub$Behavior, levels = c("Shy", "Average", "Bold"))
clmm_model <- clmm(
  Behavior ~ Age_within  + mean_Age + ClutchSize_within + mean_ClutchSize + Body.condition + I(ClutchSize^2) + (1 | FemaleID),
  data = GD_Sub
)
summary(clmm_model)


# 3.3) Emmeans: Behaviour- age, table 3: 
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

# 3.4) Plot: Behaviour and age, Fig 3A: 
M_pred_age <- ggpredict(M.Behavior, terms="Age[all]")

ggplot(M_pred_age, aes(x = x, y = predicted, group = interaction(group, response.level), color = as.factor(response.level))) +
  geom_line(size=1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = as.factor(response.level)), alpha = 0.4) +
  labs(x = "Age", y = "Probability of exhibiting the behaviour") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 16))+
  scale_color_manual(values = c("turquoise", "lightsteelblue2", "royalblue4"), 
                     labels = c("Shy", "Average", "Bold"), 
                     breaks = c("1", "2", "3")) + 
  scale_fill_manual(values = c("turquoise", "lightsteelblue2", "royalblue4"), 
                    labels = c("Shy", "Average", "Bold"), 
                    breaks = c("1", "2", "3")) + 
  theme(panel.background = element_rect(fill = "transparent", colour = NA), legend.position = c(0.8, 0.65),
        axis.line = element_line(color = "black"),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  guides(color = guide_legend(title = "Behaviour Level"), fill = FALSE) + 
  annotate("text", x= -Inf, y=Inf, label="(a)", hjust= -0.1, vjust=1.1, size=6)

##Second approach for the plot ##
PA_plot<-rbind(PA2$`Differences of cumulative probabilities`,PA3$`Differences of cumulative probabilities`,PA4$`Differences of cumulative probabilities`,PA5$Differences.of.cumulative.probabilities,PA6$Differences.of.cumulative.probabilities,PA7$Differences.of.cumulative.probabilities,PA8$Differences.of.cumulative.probabilities,PA9$Differences.of.cumulative.probabilities,PA10$Differences.of.cumulative.probabilities,PA11$Differences.of.cumulative.probabilities,PA12$Differences.of.cumulative.probabilities,PA13$Differences.of.cumulative.probabilities,PA14$Differences.of.cumulative.probabilities,PA15$Differences.of.cumulative.probabilities,PA16$Differences.of.cumulative.probabilities,PA17$Differences.of.cumulative.probabilities)
names(PA_plot)
PA_plot<- PA_plot[, 1:7]

head(PA_plot)
names(PA_plot) <- c("Behavior","Age","probability", "SE","df","LCI","UCI")


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
              alpha = 0.3) + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  xlab("Age") + ylab("probability to exhibit the Behavior") + ggtitle("(A)") + 
  theme(text = element_text(size=16), panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),legend.position = c(0.9, 0.7),
        plot.title= element_text(family = "Times New Roman", size= 16), axis.line = element_line("black"))+
  scale_color_manual(values=c("turquoise", "lightsteelblue2", "royalblue4"),labels=c("Shy", "Average", "Bold"))+
  scale_fill_manual(values=c("turquoise", "lightsteelblue2", "royalblue4"),
                    labels=c("Shy", "Average", "Bold"))

M1.P

#4) Section 3.2 "Results"
## Relationship between behaviour and clutch size and nest parasitism: 
#The same model as 3.1 
M.Behavior <- clmm(Behavior ~ Age * ClutchSize + Body.condition + I(ClutchSize^2) + (1| FemaleID), data = GD)
summary(M.Behavior)

# 4.1) Emmeans: Behavior and clutch size, table 4: 
P1 <- data.frame(summary(emmeans(M.Behavior, pairwise~ Behavior | ClutchSize,at = list(ClutchSize = 2), mode = "prob"))) 
P2 <- data.frame(summary(emmeans(M.Behavior, pairwise~ Behavior | ClutchSize,at = list(ClutchSize = 8), mode = "prob"))) 
P3 <- data.frame(summary(emmeans(M.Behavior, pairwise~ Behavior | ClutchSize,at = list(ClutchSize = 13), mode = "prob")))
P4 <- data.frame(summary(emmeans(M.Behavior, pairwise~ Behavior | ClutchSize,at = list(ClutchSize = 15), mode = "prob")))

# 4.2) Plot: Behaviour and clutch size, Fig 3B: 

M_pred_clutch <- ggpredict(M.Behavior, terms="ClutchSize[all]")

ggplot(M_pred_clutch, aes(x = x, y = predicted, group = interaction(group, response.level), color = as.factor(response.level))) +
  geom_line(size=1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = as.factor(response.level)), alpha = 0.4) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 14)) +
  labs(x = "Clutch size", y = "Probability of exhibiting the behaviour") +
  scale_color_manual(values = c("turquoise", "lightsteelblue2", "royalblue4"), 
                     labels = c("Shy", "Average", "Bold"), 
                     breaks = c("1", "2", "3")) + 
  scale_fill_manual(values = c("turquoise", "lightsteelblue2", "royalblue4"), 
                    labels = c("Shy", "Average", "Bold"), 
                    breaks = c("1", "2", "3")) + 
  theme(panel.background = element_rect(fill = "transparent", colour = NA), legend.position = c(0.9, 0.9),
        axis.line = element_line(color = "black"),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  guides(color = guide_legend(title = "Behaviour"), fill = FALSE)+
  annotate("text", x= -Inf, y=Inf, label="(b)", hjust= -0.9, vjust=1.1, size=6)


## Alternative plot approach##
PA_plot<-rbind(P1$Differences.of.cumulative.probabilities,P2$Differences.of.cumulative.probabilities,P3$Differences.of.cumulative.probabilities)
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
              alpha = 0.3)+
  xlab("Clutch size") + ylab("Probability to exhibit the behavior") + ggtitle("(B)")+ 
  theme(text = element_text(size=15), panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA), legend.position = c(0.9, 0.9),
        plot.title= element_text(family = "Times New Roman", size= 15), axis.line = element_line("black"))+ 
  scale_color_manual(values=c("turquoise", "lightsteelblue2", "royalblue4"),labels=c("Shy", "Average","Bold"))+ 
  scale_fill_manual(values=c("turquoise", "lightsteelblue2", "royalblue4"),
                    labels=c("Shy", "Average","Bold"))
M1.P 

#5) Section 3.5. "Results"
## Relationship between behavior and interaction effect of age and clutch size:

M.Behavior <- clmm(Behaviour ~ Age * ClutchSize + Body.condition + I(ClutchSize^2) + (1| FemaleID), data = GD)
summary(M.Behavior)

## 5.1) Emmeans, table 5: 
## for each coparison, replace the "Age" and "Clutchszie" values with the ones presented in table 5:

P1 <- data.frame(summary(emmeans(M.Behavior, pairwise~ Behavior | Age*ClutchSize,at = list(Age = 15, ClutchSize =13), mode = "prob")))

# 5.2) Plotting the interaction between age and clutch size in association with behavior, Fig 4:
dev.off()

colnames(GD)[colnames(GD)=="Behavior"] <- "Behaviour"
behaviour_levels <- ordered(GD$Behaviour, levels=c("Shy", "Average", "Bold"))

Effect <- effect("Age:ClutchSize", M.Behavior)
Effect_df <- as.data.frame(Effect)

behavior_colors <- c("turquoise", "lightsteelblue2", "royalblue4")
plot(Effect, multiline = TRUE, ci.style = "bars", main="", colors =behavior_colors)

### Another approach:

colorsGD <- behavior_colors[behavior_levels]

par(family="Times New Roman")

plot(Effect, main="", lines=TRUE,
     ylab="Behavior (probability of occurance)", colours= behavior_colors,
     cex.axis=2, cex.lab=3, cex.main=2)


## OR: 
preds <- ggpredict(M.Behavior, terms = c("Age [all]", "ClutchSize [all]"))

ggplot(preds, aes(x = x, y = predicted, color = group)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  labs(x = "Age", y = "Predicted Behavior", color = "Clutch Size") +
  theme_minimal()

#6) Section 3.4. "Results"
## Relationship between hatching success and behavior:  
###Beta-regression, since there are a lot of extreme values in the data:
#### Transforming the hatching success into hatch-success=number of hatched eggs and hatch-failure=number of unhatched eggs
sum(GD$hatchPercent==0)
GD$hatchSuccess <- round(GD$hatchPercent / 100 * GD$ClutchSize)
GD$hatchFailure <- GD$ClutchSize - GD$hatchSuccess

GD$Behavior<- as.factor(GD$Behavior)
GD$Behavior <- ordered(GD$Behavior, levels = c("Shy", "Average", "Bold"))
#Conventional Beta-regression
M1.Behavior <- glmmTMB(cbind(hatchSuccess, hatchFailure) ~ Behavior + Age + ClutchSize + Body.condition +
                         (1 | FemaleID), family = betabinomial(), data = GD)

#Zero-inflated Beta-regression:
M2.Behavior<- glmmTMB(cbind(hatchSuccess, hatchFailure) ~ Behavior + poly(Age,2) + poly(ClutchSize,2) + Body.condition + 
                             (1 | FemaleID), ziformula = ~1, family = betabinomial(), data = GD)
summary(M2.Behavior)

GD$Behavior <- factor(GD$Behavior, ordered= FALSE)
GD$Behavior <- relevel(GD$Behavior, ref= "Shy") 

#Model # 6.2) Goodness of fit ( table 2):
simulationOutput <- simulateResiduals(M2.Behavior) 
plot(simulationOutput)
##Dispersion
testDispersion(M2.Behavior)
# Outliers using the bootstrap method
testOutliers(M2.Behavior, type = "bootstrap")
#Alternative models:
#Null models:
common_GD <- get_all_vars(cbind(hatchSuccess, hatchFailure) ~ Behavior + poly(Age, 2) + poly(ClutchSize, 2) + Body.condition +(1 | FemaleID), data = GD)
M2.null <- glmmTMB(cbind(hatchSuccess, hatchFailure) ~ 1 + (1 | FemaleID), ziformula = ~1, family = betabinomial(), data = common_GD)
#M2.Behavior without zero inflation 
M2.NoZI <- glmmTMB(cbind(hatchSuccess, hatchFailure) ~ Behavior + poly(Age, 2) + poly(ClutchSize, 2) + Body.condition + 
                     (1 | FemaleID), family = betabinomial(), data = GD)

#Model comparison: 
AIC_BIC <- data.frame(
  Model = c("M2.Behavior", "M1.Behavior", "M2.null", "M2.NoZI"),
  AIC = c(AIC(M2.Behavior), AIC(M1.Behavior), AIC(M2.null), AIC(M2.NoZI)),
  BIC = c(BIC(M2.Behavior), BIC(M1.Behavior), BIC(M2.null), BIC(M2.NoZI)),
  logLik = c(logLik(M2.Behavior), logLik(M1.Behavior), logLik(M2.null), logLik(M2.NoZI))
)
print(AIC_BIC)

## 6.3) Plotting, Fig 5: 
dev.off()
#Ensuring that Behavior is in order
GD$Behavior <- ordered(GD$Behavior, levels = c("Shy", "Average", "Bold"))

# Computing the estimated marginal means for Behavior using emmeans, on the response scale
marginal_means <- emmeans(M2.Behavior, ~ Behavior, type = "response")
marginaleffect <- as.data.frame(marginal_means)
str(marginaleffect)

# Creating a data frame for plotting with correct column names
marginaleffect <- data.frame(
  Behavior = marginaleffect$Behavior,
  Probability = marginaleffect$prob,
  SE = marginaleffect$SE,
  CI_lower = marginaleffect$asymp.LCL,
  CI_upper = marginaleffect$asymp.UCL
)
# Defining colors for each Behavior level
behavior_colors <- c("Shy" = "turquoise", "Average" = "lightsteelblue2", "Bold" = "royalblue4")

# Ploting the estimated marginal means
P <- ggplot(marginaleffect, aes(x = Behavior, y = Probability, fill = Behavior)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +
  scale_fill_manual(values = behavior_colors) +
  labs(
    x = "Behavior", 
    y = "Hatching Success Probability", 
    title = "Estimated Marginal Means for Behavior"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

P


###################################################################
#Supplementary materials: 

# A sample dataset from females with more than 2 behavioral records
# 150 unique individuals
#Individuals with 5 behavioral record 
#1) Age and Behavior
sampled_data <- GD_Sub %>%
  group_by(FemaleID) %>%
  filter(n() >= 5) %>%
  ungroup()

set.seed(123)  
sample_size <- 150
sampled_data <- sampled_data %>%
  filter(FemaleID %in% sample(unique(FemaleID), sample_size))

# Creating a plot with separate squares for each individual, connecting with lines
ggplot(sampled_data, aes(x = Age, y = Behaviour, group = FemaleID)) +
  geom_point() +  # Plot each observation as a dot
  geom_line(aes(x = Age, y = Behaviour), color = "black", linetype = "solid", size = 0.5) +  
  theme_minimal() +
  theme(
    panel.grid = element_blank(),  # Remove the background grid
    axis.line = element_line(color = "black"),  # Keep black axis lines
    strip.background = element_blank(),  # No facet label background
    strip.text = element_text(size = 8)
  ) +
  labs(
    title = "Sample of behaviour trend Over Time",
    x = "Age",
    y = "Behaviour"
  ) +
  facet_wrap(~FemaleID, scales = "free")  # Separate squares for each individual

#2) Clutch Size and Behavior 
ggplot(sampled_data, aes(x = ClutchSize, y = Behaviour, group = FemaleID)) +
  geom_point() +  # Plot each observation as a dot
  geom_line(aes(x = ClutchSize, y = Behaviour), color = "black", linetype = "solid", size = 0.5) +  
  theme_minimal() +
  theme(
    panel.grid = element_blank(),  # Remove the background grid
    axis.line = element_line(color = "black"),  # Keep black axis lines
    strip.background = element_blank(),  # No facet label background
    strip.text = element_text(size = 8)
  ) +
  labs(
    title = "Sample of behaviour trend with different clutch size values",
    x = "Clutch size",
    y = "Behaviour"
  ) +
  
  facet_wrap(~FemaleID, scales = "free")  # Separate squares for each individual

#3) Behavior and clutchsize*age 

GD_summary <- GD %>%
  group_by(Behavior, Age) %>%
  summarize(MeanClutchSize = mean(ClutchSize), .groups = "drop")

# Plot with bar plots and jittered individual points
ggplot(GD_summary, aes(x=Behavior, y=MeanClutchSize, fill=Behavior)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.6) +
  geom_jitter(data=GD, aes(x=Behavior, y=ClutchSize, color=Behavior), 
              position=position_jitter(width=0.2), size=2) +
  facet_wrap(~ Age, nrow=1) +
  scale_fill_discrete(name="Behavior") +
  scale_color_discrete(name="Behavior") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  scale_color_manual(values=c("turquoise", "lightsteelblue2", "royalblue4"),labels=c("Shy", "Average","Bold"))+ 
  scale_fill_manual(values=c("turquoise", "lightsteelblue2", "royalblue4"))+
  labs(x="Behavior", y="Average Clutch Size", 
       title="Average Clutch Size by Behavior and Age Group")


# Saving the dataframe 
write.csv(GD, file= "GD1984_2021.csv", row.names=FALSE)  
#################### Created by :
#Farshad S. Vakili; 2023, Turku, Finland 
#E-mail: Farshad.S.Vakili@gmail.com 
#modified: 4-5.2024
################
     





