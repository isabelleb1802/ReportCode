surv <- read.csv("GLSurvival_data_final_ 2.csv")
library(survival)
library(coxme)
library(emmeans)
library(dplyr)
library(survminer)



#correcting censoring where it was mixed up in the spreadsheet
surv <- surv %>%
  mutate(
    # Correct coding:
    # status = 1 means died (event)
    # status = 0 means censored (escaped OR survived to day 32)
    status = case_when(
      Censor == 1 ~ 0,          # escaped = censored
      Lifespan == 32 ~ 0,       # survived to end = censored
      TRUE ~ 1                   # everyone else died
    )
  )
# Keeping only main columns
surv <- surv %>%
  select(Vial_ID, Cohort, Genotype, Infection1, Infection2, 
         Lifespan, Censor, Notes) %>%
  # Recoding Infection1 to match fecundity analysis
  mutate(
    Infection1 = recode(Infection1,
                        "0tp" = "Control",
                        "3tp" = "3 Tap",
                        "DTI" = "DTI"),
    Infection1 = factor(Infection1, 
                        levels = c("Control", "3 Tap", "DTI")),
    Infection2 = factor(Infection2,
                        levels = c("NON", "DTI"),
                        labels = c("Control", "DTI")),
    Genotype   = factor(Genotype),
    Genotype   = relevel(Genotype, ref = "WT"),
    Cohort     = factor(Cohort)
  )

# Removing rows with missing Lifespan
surv_clean <- surv %>%
  filter(!is.na(Lifespan))

cat("\nClean sample size:", nrow(surv_clean), "\n")
cat("Censored:", sum(surv_clean$Censor == 1), "\n")
cat("Events (deaths):", sum(surv_clean$Censor == 0), "\n")
surv_clean <- surv_clean %>%
  mutate(
    status = case_when(
      Censor == 1 ~ 0,       # escaped = censored
      Lifespan == 32 ~ 0,    # survived to end = censored  
      TRUE ~ 1               # died = event
    )
  )

cat("Deaths:", sum(surv_clean$status == 1), "\n")
cat("Censored:", sum(surv_clean$status == 0), "\n")

# Checking sample sizes
cat("\nSample size per Genotype x Infection1:\n")
print(table(surv_clean$Genotype, surv_clean$Infection1))

km_inf1 <- survfit(Surv(Lifespan, status) ~ Infection1, 
                   data = surv_clean)
ggsurvplot(km_inf1,
           data        = surv_clean,
           pval        = TRUE,
           conf.int    = TRUE,
           palette     = c("#1D3557", "#457B9D", "#E63946"),
           title       = "Survival by first infection dose",
           xlab        = "Days post first infection",
           ylab        = "Survival probability",
           legend.labs = c("Control", "3-Tap", "DTI"))

# By Genotype
km_geno <- survfit(Surv(Lifespan, status) ~ Genotype,
                   data = surv_clean)
ggsurvplot(km_geno,
           data    = surv_clean,
           pval    = TRUE,
           title   = "Survival by genotype",
           xlab    = "Days post first infection",
           ylab    = "Survival probability")

# Phag vs WT specifically
km_phag <- survfit(Surv(Lifespan, status) ~ Genotype,
                   data = surv_clean %>% 
                     filter(Genotype %in% c("WT", "Phag")))
ggsurvplot(km_phag,
           data      = surv_clean %>% filter(Genotype %in% c("WT", "Phag")),
           pval      = TRUE,
           conf.int  = TRUE,
           palette   = c("#E63946", "#999999"),
           title     = "Survival: ΔPhag vs WT",
           xlab      = "Days post first infection",
           ylab      = "Survival probability",
           legend.labs = c("ΔPhag", "WT"))


# 4. Cox proportional hazards mixed model
# Random effect for Cohort to account for batch variation
# Full model: Genotype * Infection1 * Infection2

cox_full <- coxph(Surv(Lifespan, status) ~ 
                    Genotype * Infection1 * Infection2+ 
                    frailty(Cohort),
                  data = surv_clean)

summary(cox_full)

# Stage 1 model — first infection effects only
cox_s1 <- coxph(Surv(Lifespan, status) ~ 
                  Genotype * Infection1 + 
                  frailty(Cohort),
                data = surv_clean %>% 
                  filter(Infection2 == "Control"))

summary(cox_s1)
# Additive model — genotype and first infection as main effects only
cox_s1_add <- coxph(Surv(Lifespan, status) ~ 
                      Genotype + Infection1 + 
                      frailty(Cohort),
                    data = surv_clean %>% 
                      filter(Infection2 == "Control"))

summary(cox_s1_add)


# Checking what's in the Control second infection subset
surv_ctrl <- surv_clean %>% filter(Infection2 == "Control")

cat("Rows in Control subset:", nrow(surv_ctrl), "\n")
cat("Deaths:", sum(surv_ctrl$status == 1), "\n")
cat("Censored:", sum(surv_ctrl$status == 0), "\n")

cat("\nStatus values present:", unique(surv_ctrl$status), "\n")
cat("Lifespan range:", min(surv_ctrl$Lifespan), "to", max(surv_ctrl$Lifespan), "\n")

cat("\nDeaths per genotype:\n")
print(table(surv_ctrl$Genotype, surv_ctrl$status))

cat("\nInfection2 levels present:\n")
print(table(surv_ctrl$Infection2))


# Reload raw data fresh
surv_raw <- read.csv("GLSurvival_data_final_ 2.csv")

# Building surv_clean in one single pipe
surv_clean <- surv_raw %>%
  select(Vial_ID, Cohort, Genotype, Infection1, Infection2,
         Lifespan, Censor, Notes) %>%
  filter(!is.na(Lifespan)) %>%
  mutate(
    Infection1 = case_when(
      Infection1 == "0tp" ~ "Control",
      Infection1 == "3tp" ~ "3 Tap",
      Infection1 == "DTI" ~ "DTI"
    ),
    Infection1 = factor(Infection1, levels = c("Control", "3 Tap", "DTI")),
    Infection2 = case_when(
      Infection2 == "NON" ~ "Control",
      Infection2 == "DTI" ~ "DTI"
    ),
    Infection2 = factor(Infection2, levels = c("Control", "DTI")),
    Genotype   = factor(Genotype),
    Genotype   = relevel(Genotype, ref = "WT"),
    Cohort     = factor(Cohort),
    status = case_when(
      Censor == 1    ~ 0,
      Lifespan == 32 ~ 0,
      TRUE           ~ 1
    )
  )

# Verify
cat("Infection1:\n"); print(table(surv_clean$Infection1))
cat("Infection2:\n"); print(table(surv_clean$Infection2))
cat("Status:\n");     print(table(surv_clean$status))
cat("Total rows:", nrow(surv_clean), "\n")


# Stage 1 — first infection effects only
# Filter to no second infection group
surv_s1 <- surv_clean %>% filter(Infection2 == "Control")
cat("Stage 1 sample size:", nrow(surv_s1), "\n")
cat("Deaths:", sum(surv_s1$status == 1), "\n")

cox_s1 <- coxph(Surv(Lifespan, status) ~
                  Genotype + Infection1 +
                  frailty(Cohort),
                data = surv_s1)

summary(cox_s1)

# Checking proportional hazards assumption
ph_test <- cox.zph(cox_s1)
print(ph_test)

# Kaplan-Meier by genotype for Stage 1
km_s1_geno <- survfit(Surv(Lifespan, status) ~ Genotype,
                      data = surv_s1)
ggsurvplot(km_s1_geno,
           data        = surv_s1,
           pval        = TRUE,
           conf.int    = FALSE,
           title       = "Stage 1: Survival by genotype (no re-infection)",
           xlab        = "Days post first infection",
           ylab        = "Survival probability")

# Kaplan-Meier by infection dose for Stage 1
km_s1_inf <- survfit(Surv(Lifespan, status) ~ Infection1,
                     data = surv_s1)
ggsurvplot(km_s1_inf,
           data        = surv_s1,
           pval        = TRUE,
           conf.int    = TRUE,
           palette     = c("#1D3557", "#457B9D", "#E63946"),
           title       = "Stage 1: Survival by infection dose (no re-infection)",
           xlab        = "Days post first infection",
           ylab        = "Survival probability",
           legend.labs = c("Control", "3-Tap", "DTI"))



# Stage 2 — full model with both infections
cox_s2 <- coxph(Surv(Lifespan, status) ~
                  Genotype + Infection1 + Infection2 +
                  frailty(Cohort),
                data = surv_clean)

summary(cox_s2)

# Running Phag vs WT focused model
surv_phag_wt <- surv_clean %>%
  filter(Genotype %in% c("WT", "Phag")) %>%
  mutate(Genotype = factor(Genotype, levels = c("WT", "Phag")))


cox_phag <- coxph(Surv(Lifespan, status) ~
                    Genotype * Infection1 * Infection2 +
                    frailty(Cohort),
                  data = surv_phag_wt)

summary(cox_phag)

# Focused Phag vs WT
surv_phag_wt <- surv_clean %>%
  filter(Genotype %in% c("WT", "Phag")) %>%
  mutate(Genotype = factor(Genotype, levels = c("WT", "Phag")))

cox_phag <- coxph(Surv(Lifespan, status) ~
                    Genotype * Infection1 * Infection2 +
                    frailty(Cohort),
                  data = surv_phag_wt)

summary(cox_phag)

# KM plot for Phag vs WT
km_phag <- survfit(Surv(Lifespan, status) ~ Genotype,
                   data = surv_phag_wt)

ggsurvplot(km_phag,
           data        = surv_phag_wt,
           pval        = TRUE,
           conf.int    = TRUE,
           palette     = c("#999999", "#E63946"),
           title       = "Survival: ΔPhag vs WT",
           xlab        = "Days post first infection",
           ylab        = "Survival probability",
           legend.labs = c("WT", "ΔPhag"))

# Infection1 x Infection2 interaction averaged across genotypes
cox_priming_simple <- coxph(Surv(Lifespan, status) ~
                              Genotype + Infection1 * Infection2 +
                              frailty(Cohort),
                            data = surv_clean)

summary(cox_priming_simple)

# Does genotype modify the effect of re-infection?
cox_geno_reinf <- coxph(Surv(Lifespan, status) ~
                          Genotype * Infection2 + Infection1 +
                          frailty(Cohort),
                        data = surv_clean)

summary(cox_geno_reinf)

# KM by genotype and re-infection status
surv_dif_wt <- surv_clean %>%
  filter(Genotype %in% c("WT", "Dif", "Phag")) %>%
  mutate(Genotype = factor(Genotype, levels = c("WT", "Dif", "Phag")))

km_geno_reinf <- survfit(Surv(Lifespan, status) ~ Genotype + Infection2,
                         data = surv_dif_wt)

ggsurvplot(km_geno_reinf,
           data     = surv_dif_wt,
           pval     = FALSE,
           conf.int = FALSE,
           facet.by = "Genotype",
           title    = "Survival by genotype and re-infection status",
           xlab     = "Days post first infection",
           ylab     = "Survival probability",
           palette  = c("#1D3557", "#E63946"),
           legend.labs = c("No re-infection", "Re-infected"))

# Comprehensive KM for all genotypes x re-infection
km_all_reinf <- survfit(Surv(Lifespan, status) ~ Genotype + Infection2,
                        data = surv_clean)

ggsurvplot_facet(km_all_reinf,
                 data     = surv_clean,
                 facet.by = "Genotype",
                 pval     = TRUE,
                 conf.int = FALSE,
                 palette  = c("#1D3557", "#E63946"),
                 xlab     = "Days post first infection",
                 ylab     = "Survival probability")


# Figure — 3-Tap first infection only
surv_3tap <- surv_clean %>% filter(Infection1 == "3 Tap")

km_3tap <- survfit(Surv(Lifespan, status) ~ Genotype + Infection2,
                   data = surv_3tap)

ggsurvplot_facet(km_3tap,
                 data        = surv_3tap,
                 facet.by    = "Genotype",
                 pval        = TRUE,
                 conf.int    = FALSE,
                 palette     = c("#1D3557", "#E63946"),
                 legend.labs = c("No re-infection", "Re-infected"),
                 title       = "Survival following 3-Tap first infection",
                 xlab        = "Days post first infection",
                 ylab        = "Survival probability")

# Figure — DTI first infection only
surv_dti <- surv_clean %>% filter(Infection1 == "DTI")

km_dti <- survfit(Surv(Lifespan, status) ~ Genotype + Infection2,
                  data = surv_dti)

ggsurvplot_facet(km_dti,
                 data        = surv_dti,
                 facet.by    = "Genotype",
                 pval        = TRUE,
                 conf.int    = FALSE,
                 palette     = c("#1D3557", "#E63946"),
                 legend.labs = c("No re-infection", "Re-infected"),
                 title       = "Survival following DTI first infection",
                 xlab        = "Days post first infection",
                 ylab        = "Survival probability")


# WT priming model
surv_wt_priming <- surv_clean %>%
  filter(Genotype == "WT",
         Infection1 %in% c("Control", "3 Tap")) %>%
  mutate(Infection1 = factor(Infection1, 
                             levels = c("Control", "3 Tap")))

cox_wt_priming <- coxph(Surv(Lifespan, status) ~
                          Infection1 * Infection2 +
                          frailty(Cohort),
                        data = surv_wt_priming)
summary(cox_wt_priming)

# Phag priming model
surv_phag_priming <- surv_clean %>%
  filter(Genotype == "Phag",
         Infection1 %in% c("Control", "3 Tap")) %>%
  mutate(Infection1 = factor(Infection1, 
                             levels = c("Control", "3 Tap")))

cox_phag_priming <- coxph(Surv(Lifespan, status) ~
                            Infection1 * Infection2 +
                            frailty(Cohort),
                          data = surv_phag_priming)
summary(cox_phag_priming)

# Phag vs WT priming comparison
surv_phag_wt_priming <- surv_clean %>%
  filter(Genotype %in% c("WT", "Phag"),
         Infection1 %in% c("Control", "3 Tap")) %>%
  mutate(
    Genotype   = factor(Genotype, levels = c("WT", "Phag")),
    Infection1 = factor(Infection1, levels = c("Control", "3 Tap"))
  )

cox_phag_wt_priming <- coxph(Surv(Lifespan, status) ~
                               Genotype * Infection1 * Infection2 +
                               frailty(Cohort),
                             data = surv_phag_wt_priming)
summary(cox_phag_wt_priming)

library(survminer)

# Filtering to WT and Phag, Control and 3-Tap first infection only
surv_priming_plot <- surv_clean %>%
  filter(Genotype %in% c("WT", "Phag"),
         Infection1 %in% c("Control", "3 Tap")) %>%
  mutate(
    Genotype = factor(Genotype, levels = c("WT", "Phag")),
    Infection1 = factor(Infection1, levels = c("Control", "3 Tap"),
                        labels = c("Unprimed", "3-Tap primed"))
  )

surv_priming_plot <- surv_clean %>%
  filter(Genotype %in% c("WT", "Phag"),
         Infection1 %in% c("Control", "3 Tap"),
         as.character(Infection2) == "DTI") %>%
  mutate(
    Genotype = factor(Genotype, 
                      levels = c("WT", "Phag"),
                      labels = c("Wildtype (WT)", "ΔPhag")),
    Infection1 = factor(Infection1, 
                        levels = c("Control", "3 Tap"),
                        labels = c("Unprimed", "3-Tap primed"))
  )
  )
# KM fit
km_priming <- survfit(Surv(Lifespan, status) ~ Infection1,
                      data = surv_priming_plot)


# Faceted plot
ggsurvplot_facet(
  survfit(Surv(Lifespan, status) ~ Infection1,
          data = surv_priming_plot),
  data         = surv_priming_plot,
  facet.by     = "Genotype",
  censor = FALSE, 
  pval         = TRUE,
  conf.int     = FALSE,
  palette      = c("#999999", "#E63946"),
  legend.labs  = c("Unprimed", "3-Tap primed"),
  legend.title = "Prior infection",
  xlab         = "Days post re-infection",
  ylab         = "Survival probability",
  ggtheme      = theme_classic(),
  short.panel.labs = TRUE  # removes the "Genotype:" prefix
) +
  theme(
    plot.title       = element_text(face = "bold"),
    strip.text       = element_text(face = "bold", size = 12),
    strip.background = element_blank()
  )

ggsave("Figure_priming_KM.png",
       width = 10, height = 5, dpi = 300)
table(surv_priming_plot$Genotype, 
      surv_priming_plot$Infection1, 
      surv_priming_plot$Infection2)
