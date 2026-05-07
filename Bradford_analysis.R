
bradford <- read.csv("Bradford assay results.csv")

# Ensure factors
bradford$Genotype   <- factor(bradford$Genotype)
bradford$Infection1 <- factor(bradford$Infection1, 
                              levels = c("Control", "3 Tap", "DTI"))
bradford$Cohort     <- factor(bradford$Cohort)

# Set reference levels
bradford$Genotype   <- relevel(bradford$Genotype,   ref = "WT")
bradford$Infection1 <- relevel(bradford$Infection1, ref = "Control")

# Exclude the accidental DTI second infection vial
bradford <- bradford %>% filter(Vial_ID != "1C-02")

# -----------------------------------------------------------------------------
# ANOVA / linear model
# Genotype * Infection1 as fixed effects, Cohort as blocking factor


brad_lm <- lm(
  `Protein.Conc..ul.ml.` ~ Cohort + Genotype * Infection1,
  data = bradford
)

# ANOVA table
anova(brad_lm)

# -----------------------------------------------------------------------------
# Checking normality of residuals 
# -----------------------------------------------------------------------------

# 1. QQ plot of residuals
qqnorm(residuals(brad_lm), main = "QQ Plot — Bradford residuals")
qqline(residuals(brad_lm), col = "red")

# 2. Histogram of residuals
hist(residuals(brad_lm), 
     breaks = 15, 
     main = "Histogram of residuals",
     xlab = "Residuals")

# 3. Formal normality test
shapiro.test(residuals(brad_lm))

# 4. Residuals vs fitted plot (check homogeneity of variance)
plot(fitted(brad_lm), residuals(brad_lm),
     xlab = "Fitted values",
     ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red")

# Shift all values above zero then log transform
bradford$log_shifted <- log(bradford$`Protein.Conc..ul.ml.` + 73)

# Refit ANOVA
brad_lm_log <- lm(
  log_shifted ~ Cohort + Genotype * Infection1,
  data = bradford
)

# Check normality
qqnorm(residuals(brad_lm_log), main = "QQ Plot — log shifted Bradford residuals")
qqline(residuals(brad_lm_log), col = "red")
shapiro.test(residuals(brad_lm_log))
hist(residuals(brad_lm_log), breaks = 15, main = "Histogram — log shifted residuals")


# ANOVA table
anova(brad_lm_log)

# Emmeans for infection effect
library(emmeans)
emm_brad_inf <- emmeans(brad_lm_log, ~ Infection1)
print(emm_brad_inf)
pairs(emm_brad_inf, adjust = "tukey")

# Emmeans for genotype effect
emm_brad_geno <- emmeans(brad_lm_log, ~ Genotype)
print(emm_brad_geno)
pairs(emm_brad_geno, adjust = "tukey")

# Interaction
emm_brad_int <- emmeans(brad_lm_log, ~ Infection1 | Genotype)
print(emm_brad_int)
pairs(emm_brad_int, adjust = "tukey")

library(lme4)

# LMM version with cohort as random effect
brad_lmm <- lmer(
  log_shifted ~ Genotype * Infection1 + (1 | Cohort),
  data = bradford
)

# AIC
AIC(brad_lm_log, brad_lmm)

qqnorm(residuals(brad_lmm), main = "QQ Plot — LMM Bradford residuals")
qqline(residuals(brad_lmm), col = "red")
shapiro.test(residuals(brad_lmm))

# ANOVA table for LMM
install.packages("lmerTest")
library(lmerTest) 
brad_lmm <- lmer(
  log_shifted ~ Genotype * Infection1 + (1 | Cohort),
  data = bradford
)

anova(brad_lmm)
shapiro.test(residuals(brad_lmm))

# Emmeans
emm_brad <- emmeans(brad_lmm, ~ Infection1)
print(emm_brad)
pairs(emm_brad, adjust = "tukey")

emm_brad_geno <- emmeans(brad_lmm, ~ Genotype)
print(emm_brad_geno)

# Additive LMM for emmeans only
brad_lmm_add <- lmer(
  log_shifted ~ Genotype + Infection1 + (1 | Cohort),
  data = bradford
)

# Checking it has no nonEst values
emm_brad_inf_add <- emmeans(brad_lmm_add, ~ Infection1)
print(emm_brad_inf_add)
pairs(emm_brad_inf_add, adjust = "tukey")
emm_brad_geno_add <- emmeans(brad_lmm_add, ~ Genotype)
print(emm_brad_geno_add)

pairs(emm_brad_geno_add, adjust = "tukey")
