data<-read.csv("Cleaned GL Egg Data.csv")
data$Genotype   <- factor(data$Genotype)
data$Infection2 <- factor(data$Infection2, levels = c("Control", "DTI"))
data$Cohort     <- factor(data$Cohort)

 # Setting WT / Control / Control as the reference level
  data$Genotype   <- relevel(data$Genotype,   ref = "WT")
  data$Infection1 <- relevel(data$Infection1, ref = "Control")
  data$Infection2 <- relevel(data$Infection2, ref = "Control")
    
    library(glmmTMB)
    library(emmeans)
    library(dplyr)
    

    # Null model
    null_fec <- glmmTMB(Eggs ~ 1 + (1 | Cohort),
                        family = nbinom2,
                        data = data_s1)
    
    # Extracting components following Nakagawa et al. 2017
    cohort_var <- VarCorr(null_fec)$cond$Cohort[1]
    sigma <- sigma(null_fec)  # overdispersion parameter
    
    # Mean lambda from fixed intercept
    lambda <- exp(fixef(null_fec)$cond[1] + cohort_var/2)
    
    # Residual variance - precise negative binomial formula
    residual_var <- log(1/lambda + 1/sigma + 1)
    
    # ICC
    ICC <- cohort_var / (cohort_var + residual_var)
    
    # Design effect
    avg_cluster_size <- nrow(data_s1) / length(unique(data_s1$Cohort))
    DEFF <- 1 + (avg_cluster_size - 1) * ICC
    
    cat("Cohort variance:", round(cohort_var, 3), "\n")
    cat("Lambda:", round(lambda, 3), "\n")
    cat("Residual variance:", round(residual_var, 3), "\n")
    cat("ICC:", round(ICC, 3), "\n")
    cat("Average cluster size:", round(avg_cluster_size, 1), "\n")
    cat("Design effect:", round(DEFF, 3), "\n")
    
    
    # POISSON MODEL
    
    data_all <- data %>% 
      filter(!is.na(Eggs))
    data_all$Eggs <- as.numeric(data_all$Eggs)
    
    data_s1$Eggs <- as.numeric(data_s1$Eggs)
    fec_s1 <- glmmTMB(
      Eggs ~ Genotype * Infection1 + (1 | Cohort),
      family = nbinom2,
      data   = data_s1
    )
    
    summary(fec_s1)
    # Fit Poisson version
    fec_s1_pois <- glmmTMB(
      Eggs ~ Genotype * Infection1 + (1 | Cohort),
      family = poisson,
      data = data_s1
    )
    summary(fec_s1_pois)
    # Comparing with negative binomial using AIC
    AIC(fec_s1_pois, fec_s1)
    
    # Checking Poisson diagnostics
    install.packages("DHARMa")
    library(DHARMa)
    
    sim_pois <- simulateResiduals(fec_s1_pois, n = 1000)
    plot(sim_pois)
    testDispersion(sim_pois)
    
    
    # Zero-inflation varies by genotype
    fec_s1_zi_geno <- glmmTMB(
      Eggs ~ Genotype * Infection1 + (1 | Cohort),
      ziformula = ~ Genotype,
      family = nbinom2,
      data = data_s1
    )
    
    # Comparing AIC
    AIC(fec_s1, fec_s1_zi_geno)
    
    
    sim_zi_geno <- simulateResiduals(fec_s1_zi_geno, n = 1000)
    plot(sim_zi_geno)
    testDispersion(sim_zi_geno)
    testZeroInflation(sim_zi_geno)
    
    data_s2 <- data %>%
      filter(!is.na(Eggs))
    
    cat("\nStage 2 sample size:", nrow(data_s2), "vials\n")
    cat("Vials per Infection1 x Infection2:\n")
    print(table(data_s2$Infection1, data_s2$Infection2))
    
    data_s2$Eggs <- as.numeric(data_s2$Eggs)
    fec_s2 <- glmmTMB(
      Eggs ~ Genotype * Infection1 * Infection2 + (1 | Cohort),
      family = nbinom2,
      data   = data_s2
    )
    
    summary(fec_s2)
    
    # Poisson version of Stage 2
    fec_s2_pois <- glmmTMB(
      Eggs ~ Genotype * Infection1 * Infection2 + (1 | Cohort),
      family = poisson,
      data = data_s2
    )
    
    summary(fec_s2_pois)
    
    # AIC comparison
    AIC(fec_s2_pois, fec_s2)
    
    # Zero-inflation tests for Stage 2
    fec_s2_zi <- glmmTMB(
      Eggs ~ Genotype * Infection1 * Infection2 + (1 | Cohort),
      ziformula = ~1,
      family = nbinom2,
      data = data_s2
    )
    
    AIC(fec_s2, fec_s2_zi)
    
    # DHARMa diagnostics on existing Stage 2 model
    sim_s2 <- simulateResiduals(fec_s2, n = 1000)
    plot(sim_s2)
    testDispersion(sim_s2)
    testZeroInflation(sim_s2)
    
    # --- excluding vials with 0 eggs etc. 
    
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
    
    # BB with ZI
    via_s1_zbb <- glmmTMB(
      cbind(Hatched, Eggs - Hatched) ~ Genotype + Infection1 + (1 | Cohort),
      ziformula = ~1,
      family = betabinomial,
      data = data_v1
    )
    
    # Beta-binomial without zero-inflation
    via_s1_bb <- glmmTMB(
      cbind(Hatched, Eggs - Hatched) ~ Genotype + Infection1 + (1 | Cohort),
      family = betabinomial,
      data = data_v1
    )
    
    # Full model selection comparison
    AIC(via_s1_bb, via_s1_zbb)
    