# Fecundity analysis — egg count, GLMMs with negative binomial distribution
# Stage 1: effect of first infection only (Infection2 = "Control")
# Stage 2: combined effect of first and second infection
# Viability Analysis 
# Individual Analysis comparing Phag and WT

library(glmmTMB)
library(emmeans)
library(dplyr)


data <- read.csv("Cleaned GL Egg Data.csv")

data$Genotype   <- factor(data$Genotype)
data$Infection1 <- factor(data$Infection1, levels = c("Control", "3 Tap", "DTI"))
data$Infection2 <- factor(data$Infection2, levels = c("Control", "DTI"))
data$Cohort     <- factor(data$Cohort)

data$Genotype   <- relevel(data$Genotype,   ref = "WT")
data$Infection1 <- relevel(data$Infection1, ref = "Control")
data$Infection2 <- relevel(data$Infection2, ref = "Control")

# STAGE 1 — Effect of first infection only
# Subset: Infection2 == "Control", drop NA egg counts

data_s1 <- data %>%
  filter(Infection2 == "Control", !is.na(Eggs))

cat("Stage 1 sample size:", nrow(data_s1), "vials\n")
cat("Vials per Genotype x Infection1:\n")
print(table(data_s1$Genotype, data_s1$Infection1))

data_all <- data %>% 
  filter(!is.na(Eggs))
data_all$Eggs <- as.numeric(data_all$Eggs)

cat("Total vials:", nrow(data), "\n")
cat("Vials per Genotype x Infection1:\n")
print(table(data$Genotype, data$Infection1))

#GLMM model with negative binomial distribution comparing 1st infection dose and genotype interactions
data_s1$Eggs <- as.numeric(data_s1$Eggs)
fec_s1 <- glmmTMB(
  Eggs ~ Genotype * Infection1 + (1 | Cohort),
  family = nbinom2,
  data   = data_s1)

summary(fec_s1)


# Identify cells with potential separation (SE > 50 is a warning sign)
s1_coefs <- as.data.frame(summary(fec_s1)$coefficients$cond)
separation_flags <- s1_coefs[s1_coefs$`Std. Error` > 50, ]
if (nrow(separation_flags) > 0) {
  cat("\nWARNING — potential separation in these terms (SEs very large):\n")
  print(separation_flags)
}


# Stage 1 emmeans — infection effect within each genotype
# Back-transformed to count scale with type = "response"

library(DHARMa)
sim_s1 <- simulateResiduals(fec_s1, n = 1000)
plot(sim_s1)
testDispersion(sim_s1)
testZeroInflation(sim_s1)

fec_s1_zi <- glmmTMB(
  Eggs ~ Genotype * Infection1 + (1 | Cohort),
  ziformula = ~1,
  family = nbinom2,
  data = data_s1
)

AIC(fec_s1, fec_s1_zi)

sim_s1_zi <- simulateResiduals(fec_s1_zi, n = 1000)
plot(sim_s1_zi)
testDispersion(sim_s1_zi)
testZeroInflation(sim_s1_zi)
# Predicted counts per Infection1 level within each genotype
emm_s1 <- emmeans(fec_s1, ~ Infection1 | Genotype, type = "response")
print(emm_s1)

# Pairwise contrasts: all infection levels compared within each genotype
# Bonferroni correction for 3 pairwise comparisons
pairs_s1 <- pairs(emm_s1, adjust = "Bonferroni")
print(pairs_s1)

# Ratio-scale contrasts
pairs_s1_ratio <- pairs(emm_s1, adjust = "Bonferroni", type = "response")
print(pairs_s1_ratio)

# testing dose-response direction:
# Control vs 3-Tap, Control vs DTI, 3-Tap vs DTI
contrast(emm_s1,
  list(
    "Control vs 3-Tap" = c( 1, -1,  0),
    "Control vs DTI"   = c( 1,  0, -1),
    "3-Tap vs DTI"     = c( 0,  1, -1)
  ),
  adjust = "bonferroni",
  by     = "Genotype")

# STAGE 2 — Combined effect of first and second infection
# Full dataset, dropping NA egg counts


data_s2 <- data %>%
  filter(!is.na(Eggs))

cat("\nStage 2 sample size:", nrow(data_s2), "vials\n")
cat("Vials per Infection1 x Infection2:\n")
print(table(data_s2$Infection1, data_s2$Infection2))

# Stage 2 model
# Full three-way interaction: Genotype x Infection1 x Infection2 


data_s2$Eggs <- as.numeric(data_s2$Eggs)
fec_s2 <- glmmTMB(
  Eggs ~ Genotype * Infection1 * Infection2 + (1 | Cohort),
  family = nbinom2,
  data   = data_s2
)

summary(fec_s2)

# checking for separation
s2_coefs <- as.data.frame(summary(fec_s2)$coefficients$cond)
separation_flags_s2 <- s2_coefs[s2_coefs$`Std. Error` > 50, ]
if (nrow(separation_flags_s2) > 0) {
  cat("\nWARNING — potential separation in Stage 2 (SEs very large):\n")
  print(separation_flags_s2)
}

# Stage 2 emmeans

# Predicted counts for each Infection1 x Infection2 combination, per genotype
emm_s2 <- emmeans(fec_s2, ~ Infection1 * Infection2 | Genotype, type = "response")
print(emm_s2)

# Pairwise comparisons of Infection1 within each Genotype,
# averaged over Infection2 (overall first-infection effect)
pairs_s2_inf1 <- pairs(
  emmeans(fec_s2, ~ Infection1 | Genotype, type = "response"),
  adjust = "Bonferroni"
)
print(pairs_s2_inf1)

# Effect of second infection within each Genotype x Infection1 combination
# (does re-infection change egg count relative to no re-infection?)
pairs_s2_inf2 <- pairs(
  emmeans(fec_s2, ~ Infection2 | Genotype * Infection1, type = "response"),
  adjust = "bonferroni"
)
print(pairs_s2_inf2)

# Key terminal investment contrast:
# Does DTI first infection + DTI second infection increase reproduction relative to DTI first infection + no second infection?
emm_s2_ti <- emmeans(fec_s2, ~ Infection2 | Genotype,
                     at     = list(Infection1 = "DTI"),
                     type   = "response")
pairs(emm_s2_ti, adjust = "bonferroni")



library(DHARMa)
# Stage 1
sim_s1 <- simulateResiduals(fec_s1, n = 1000)
plot(sim_s1)                
testDispersion(sim_s1)      
testZeroInflation(sim_s1)   

# Stage 2
sim_s2 <- simulateResiduals(fec_s2, n = 1000)
plot(sim_s2)
testDispersion(sim_s2)
testZeroInflation(sim_s2)

#Viability analysis
#Getting rid of 0s and miscounts 
data$Eggs <- as.numeric(data$Eggs)
data$Hatched <- as.numeric(data$Hatched)

exclude_vials <- c("1A-08", "1I-03", "2A-06", "2D-06", "1F-06", 
                   "3B-09", "3E-09", "4B-09", "2A-09", "2B-04", 
                   "2B-09", "3B-01", "3B-03", "3D-01", "3D-03")

data_v1 <- data %>%
  filter(
    Infection2 == "Control",
    !is.na(Eggs),
    !is.na(Hatched),
    Eggs > 0,
    Hatched <= Eggs,
    !(Vial_ID %in% exclude_vials)
  )

cat("Stage 1 viability sample size:", nrow(data_v1), "\n")
print(table(data_v1$Genotype, data_v1$Infection1))
# this should be negative binomial or beta-binomial
via_s1 <- glmmTMB(
  cbind(Hatched, Eggs - Hatched) ~ Genotype * Infection1 + (1 | Cohort),
  family = binomial,
  data   = data_v1
)

summary(via_s1)

via_s1_add <- glmmTMB(
  cbind(Hatched, Eggs - Hatched) ~ Genotype + Infection1 + (1 | Cohort),
  family = binomial,
  data   = data_v1
)

summary(via_s1_add)

sim_v1 <- simulateResiduals(via_s1_add, n = 1000)
plot(sim_v1)
testDispersion(sim_v1)
testZeroInflation(sim_v1)

via_s1_zbb <- glmmTMB(
  cbind(Hatched, Eggs - Hatched) ~ Genotype + Infection1 + (1 | Cohort),
  ziformula = ~1,
  family = betabinomial,
  data = data_v1
)

summary(via_s1_zbb)

# Diagnostics
sim_v1_zbb <- simulateResiduals(via_s1_zbb, n = 1000)
plot(sim_v1_zbb)
testDispersion(sim_v1_zbb)
testZeroInflation(sim_v1_zbb)

# Comparing AIC with previous model
AIC(via_s1_add, via_s1_zbb)

emm_v1 <- emmeans(via_s1_zbb, ~ Infection1, type = "response")
print(emm_v1)
pairs(emm_v1, adjust = "tukey")

emm_v1_geno <- emmeans(via_s1_zbb, ~ Genotype, type = "response")
print(emm_v1_geno)
pairs(emm_v1_geno, adjust = "tukey")

# Focusing to WT and Phag only
data_phag_wt <- data_s2 %>%
  filter(Genotype %in% c("WT", "Phag")) %>%
  mutate(Genotype = factor(Genotype, levels = c("WT", "Phag")))

cat("Sample size:", nrow(data_phag_wt), "\n")
print(table(data_phag_wt$Genotype, 
            data_phag_wt$Infection1, 
            data_phag_wt$Infection2))

# Full three-way interaction model
fec_phag_wt <- glmmTMB(
  Eggs ~ Genotype * Infection1 * Infection2 + (1 | Cohort),
  family = nbinom2,
  data   = data_phag_wt
)

summary(fec_phag_wt)

# Diagnostics
library(DHARMa)
sim_phag_wt <- simulateResiduals(fec_phag_wt, n = 1000)
plot(sim_phag_wt)
testDispersion(sim_phag_wt)
testZeroInflation(sim_phag_wt)

# Emmeans — predicted counts for all combinations
emm_phag_wt <- emmeans(fec_phag_wt, 
                       ~ Infection1 * Infection2 | Genotype,
                       type = "response")
print(emm_phag_wt)

# effect of re-infection at each first infection dose
pairs(
  emmeans(fec_phag_wt, ~ Infection2 | Genotype * Infection1,
          type = "response"),
  adjust = "bonferroni"
)
pairs(
  emmeans(fec_phag_wt, ~ Genotype | Infection1 * Infection2,
          type = "response"),
  adjust = "bonferroni"
)

# Raw predicted counts for Phag vs WT from focused model
emm_phag_plot <- as.data.frame(summary(emm_phag_wt)) %>%
  mutate(
    Infection_combo = paste(Infection1, Infection2, sep = "\n"),
    Infection_combo = factor(Infection_combo, levels = c(
      "Control\nControl", "3 Tap\nControl", "DTI\nControl",
      "Control\nDTI", "3 Tap\nDTI", "DTI\nDTI"
    )),
    line_type = ifelse(Genotype == "WT", "WT", "Phag")
  )
library(ggplot2)

fig_phag_raw <- ggplot(emm_phag_plot,
                       aes(x = Infection_combo, y = response,
                           colour = Genotype, group = Genotype,
                           linetype = line_type)) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                width = 0.2, linewidth = 0.7) +
  geom_point(size = 3) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = c("Phag" = "#E63946", "WT" = "#999999")) +
  geom_vline(xintercept = 3.5, linetype = "dotted",
             colour = "grey60", linewidth = 0.5) +
  guides(linetype = "none") + 
  annotate("text", x = 2, y = max(emm_phag_plot$asymp.UCL, na.rm = TRUE) * 0.95,
           label = "No re-infection", size = 3, colour = "grey40") +
  annotate("text", x = 5, y = max(emm_phag_plot$asymp.UCL, na.rm = TRUE) * 0.95,
           label = "Re-infected", size = 3, colour = "grey40") +
  annotate("text", x = 6, y = 108,
           label = "*", size = 8, colour = "#E63946")+
  labs(
    x       = "Infection combination\n(First infection / Second infection)",
    y       = "Predicted egg count",
    title   = "ΔPhag shows dose-specific increase in fecundity following re-infection",
    caption  = "Predicted values from focused negative binomial GLMM (WT and ΔPhag only).\n* DTI/DTI re-infection effect in ΔPhag: p = 0.0071 (Bonferroni corrected)"
  ) +
  theme_classic() +
  theme(
    legend.position  = "right",
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(size = 9, colour = "grey40"),
    axis.text.x      = element_text(size = 8)
  )

print(fig_phag_raw)
ggsave("Figure_Phag_raw_counts.png", fig_phag_raw,
       width = 8, height = 5, dpi = 300)

# Comparing each genotype to WT at each infection dose
emm_s1_geno <- emmeans(fec_s1, ~ Genotype | Infection1, type = "response")

# Contrasts of each genotype vs WT within each infection level
contrast(emm_s1_geno,
         method = "trt.vs.ctrl",  # compares all levels to reference (WT)
         ref    = "WT",
         adjust = "bonferroni")
# Stage 2 — compare each genotype to WT at each 
# Infection1 x Infection2 combination
emm_s2_geno <- emmeans(fec_s2, 
                       ~ Genotype | Infection1 * Infection2, 
                       type = "response")

contrast(emm_s2_geno,
         method = "trt.vs.ctrl",
         ref    = "WT",
         adjust = "bonferroni")
