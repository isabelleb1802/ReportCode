# Figure code — fecundity, viability, survival and Bradford analyses
# Uses ggplot2 and emmeans predicted values 

library(ggplot2)
library(emmeans)
library(dplyr)

geno_colours <- c(
  "WT"   = "#1D3557",  
  "Dif"  = "#457B9D", 
  "Mel"  = "#A8DADC", 
  "Phag" = "#E63946", 
  "BOM"  = "#88CCEE",  
  "Toll" = "#6D6875" )  # dusty mauve

  
# =============================================================================
# FIGURE 1 — Stage 1 fecundity
# Predicted egg counts per infection dose for significant/marginal genotypes
# Compared to WT reference line
# =============================================================================
# Extracting emmeans for significant genotypes AND WT
sig_genos <- c("WT", "Dif", "Mel", "Phag", "BOM", "Toll")

emm_s1_sig <- as.data.frame(summary(emm_s1)) %>%
  filter(Genotype %in% sig_genos) %>%
  mutate(
    Genotype = factor(Genotype, levels = sig_genos),
    Infection1 = factor(Infection1, levels = c("Control", "3 Tap", "DTI")),
    line_type = ifelse(Genotype == "WT", "WT", "Genotype")
  )

# Creating a version with WT repeated for each genotype panel
wt_data <- emm_s1_sig %>%
  filter(Genotype == "WT") %>%
  select(-Genotype)

# binding the genotype data with WT data
emm_s1_plot <- emm_s1_sig %>%
  filter(Genotype != "WT")

# Adding WT as a separate repeated dataset for each panel
wt_repeated <- bind_rows(
  lapply(c("Dif", "Mel", "Phag", "BOM", "Toll"), function(g) {
    wt_data %>% mutate(Genotype = g, line_type = "WT")
  })
) %>%
  mutate(Genotype = factor(Genotype, levels = c("Dif", "Mel", "Phag", "BOM", "Toll")))

# Combining
plot_data <- bind_rows(emm_s1_plot, wt_repeated) %>%
  mutate(Genotype = factor(Genotype, levels = c("Dif", "Mel", "Phag", "BOM", "Toll")))


# Significance annotations — genotype vs WT at DTI only
sig_ann_s1 <- data.frame(
  Genotype = c("Mel", "Dif", "Toll", "Phag"),
  Infection1 = c("DTI", "DTI", "DTI", "DTI"),
  label = c("*", "*", "*", "*")
) %>%
  left_join(
    plot_data_s1_final %>% 
      select(Genotype, Infection1, response, SE),
    by = c("Genotype", "Infection1")
  ) %>%
  mutate(Genotype = factor(Genotype, 
                           levels = c("# Significance annotations — genotype vs WT at DTI only 
                           
sig_ann_s1 <- data.frame(
  Genotype = c("Mel", "Dif", "Toll", "Phag"),
  Infection1 = c("DTI", "DTI", "DTI", "DTI"),
  label = c("*", "*", "*", "*")
) %>%
  left_join(
    plot_data_s1_final %>% 
      select(Genotype, Infection1, response, SE),
    by = c("Genotype", "Infection1")
  ) %>%
  mutate(Genotype = factor(Genotype, 
                           levels = c("WT", "Dif", "Mel", "Phag", "BOM", "Toll")))


# Significance annotations — genotype vs WT at DTI only
sig_ann_s1 <- data.frame(
  Genotype = c("Mel", "Dif", "Toll", "Phag"),
  Infection1 = c("DTI", "DTI", "DTI", "DTI"),
  label = c("*", "*", "*", "*")
) %>%
  left_join(
    plot_data_s1_final %>% 
      select(Genotype, Infection1, response, SE),
    by = c("Genotype", "Infection1")
  ) %>%
  mutate(Genotype = factor(Genotype, 
                           levels = c("WT", "Dif", "Mel", "Phag", "BOM", "Toll")))

fig1_bars <- ggplot(plot_data_s1_final,
                    aes(x = Infection1, y = response, fill = Genotype)) +
  geom_col(width = 0.6, alpha = 0.85) +
  geom_errorbar(aes(ymin = response - SE, ymax = response + SE),
                width = 0.2, linewidth = 0.7) +
  geom_text(data = sig_ann_s1,
            aes(x = Infection1, y = response + SE + 1.5, label = label),
            inherit.aes = FALSE,
            size = 6,
            colour = "black") +
  scale_fill_manual(
    values = c(
      "WT"   = "#999999",
                                      "Dif"  = "#457B9D",
                                      "Mel"  = "#A8DADC",
                                      "Phag" = "#E63946",
                                      "BOM"  = "#6D6875",
                                      "Toll" = "#1D3557"
                           ),
                           name = "Genotype"
  ) +
    scale_x_discrete(labels = c("Control", "3-Tap", "DTI")) +
    coord_cartesian(ylim = c(0, 40)) +
    facet_wrap(~ Genotype, nrow = 1) +
    labs(
      x       = "First infection treatment",
      y       = "Predicted egg count",
      title   = "Stage 1: Fecundity response to first infection",
      subtitle = "Predicted egg counts per genotype and infection dose",
      caption = "Predicted values from negative binomial GLMM. Error bars show ±1 SE.\nSome values extend beyond plot range.\n* Significantly lower than WT at DTI (Bonferroni corrected)"
    ) +
    theme_classic() +
    theme(
      legend.position  = "none",
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold", size = 13),
      axis.text.x      = element_text(angle = 45, hjust = 1),
      plot.title       = element_text(face = "bold"),
      plot.subtitle    = element_text(size = 9, colour = "grey40"),
      panel.spacing    = unit(1.5, "lines")
    )
  
  print(fig1_bars)
  ggsave("Figure1_stage1_final.png", fig1_bars,
         width = 14, height = 5, dpi = 300)
  
  
  
  
  
  
# =============================================================================
# FIGURE 2 UPDATED — Stage 2: change in egg count due to re-infection
# across first infection doses
# Y = predicted eggs (re-infected) - predicted eggs (no re-infection)
# X = first infection dose
# Two lines per panel: genotype vs WT
# =============================================================================

# Extract full Stage 2 emmeans for all genotypes x Infection1 x Infection2
emm_s2_all <- as.data.frame(summary(
  emmeans(fec_s2, ~ Infection1 * Infection2 | Genotype, type = "response")
)) %>%
  mutate(Infection1 = factor(Infection1, levels = c("Control", "3 Tap", "DTI")))

# Calculate re-infection change (DTI2 - Control2) per genotype x Infection1
emm_s2_no  <- emm_s2_all %>%
  filter(Infection2 == "Control") %>%
  select(Genotype, Infection1, resp_no = response)

emm_s2_yes <- emm_s2_all %>%
  filter(Infection2 == "DTI") %>%
  select(Genotype, Infection1, resp_yes = response)

emm_s2_delta <- left_join(emm_s2_no, emm_s2_yes,
                          by = c("Genotype", "Infection1")) %>%
  mutate(delta = resp_yes - resp_no)

# Separate WT delta
wt_delta <- emm_s2_delta %>%
  filter(Genotype == "WT") %>%
  select(Infection1, wt_delta = delta)

# Significant genotypes
sig_genos_s2 <- c("Dif", "Mel", "Phag", "BOM", "Toll")

# Genotype delta data
geno_delta <- emm_s2_delta %>%
  filter(Genotype %in% sig_genos_s2) %>%
  mutate(
    Genotype  = factor(Genotype, levels = sig_genos_s2),
    line_type = "Genotype"
  )

# Repeat WT delta for each genotype panel
wt_delta_repeated <- bind_rows(
  lapply(sig_genos_s2, function(g) {
    wt_delta %>%
      mutate(Genotype = g, delta = wt_delta, line_type = "WT")
  })
) %>%
  mutate(Genotype = factor(Genotype, levels = sig_genos_s2))

# Combine
plot_data_s2 <- bind_rows(
  geno_delta %>% select(Genotype, Infection1, delta, line_type),
  wt_delta_repeated %>% select(Genotype, Infection1, delta, line_type)
) %>%
  mutate(Genotype = factor(Genotype, levels = sig_genos_s2))

# Significance annotation — Phag at DTI
sig_ann <- plot_data_s2 %>%
  filter(Genotype == "Phag", line_type == "Genotype",
         Infection1 == "DTI") %>%
  mutate(label = "*")

# Extract Stage 2 emmeans for significant genotypes at all infection combinations
sig_genos_s2 <- c("Dif", "Mel", "Phag", "BOM", "Toll")

emm_s2_plot <- as.data.frame(summary(
  emmeans(fec_s2, ~ Infection1 * Infection2 | Genotype, type = "response")
)) %>%
  filter(Genotype %in% sig_genos_s2) %>%
  mutate(
    Infection_combo = paste(Infection1, Infection2, sep = "\n"),
    Infection_combo = factor(Infection_combo, levels = c(
      "Control\nControl", "3 Tap\nControl", "DTI\nControl",
      "Control\nDTI", "3 Tap\nDTI", "DTI\nDTI"
    )),
    Genotype = factor(Genotype, levels = sig_genos_s2)
  )

# Significance annotation for Phag DTI/DTI
sig_ann_s2 <- emm_s2_plot %>%
  filter(Genotype == "Phag", 
         Infection_combo == "DTI\nDTI") %>%
  mutate(label = "*")

fig2_raw <- ggplot(emm_s2_plot,
                   aes(x = Infection_combo, y = response, fill = Genotype)) +
  geom_col(width = 0.7, alpha = 0.85) +
  geom_errorbar(aes(ymin = response - SE, ymax = response + SE),
                width = 0.2, linewidth = 0.7) +
  geom_text(data = sig_ann_s2,
            aes(x = Infection_combo, y = response + SE, label = label),
            inherit.aes = FALSE,
            colour = "#E63946",
            size = 7,
            vjust = -0.3) +
  geom_vline(xintercept = 3.5, linetype = "dotted",
             colour = "grey60", linewidth = 0.5) +
  annotate("text", x = 2, y = max(emm_s2_plot$response + emm_s2_plot$SE, 
                                  na.rm = TRUE) * 0.95,
           label = "No re-infection", size = 3, colour = "grey40") +
  annotate("text", x = 5, y = max(emm_s2_plot$response + emm_s2_plot$SE, 
                                  na.rm = TRUE) * 0.95,
           label = "Re-infected", size = 3, colour = "grey40") +
  scale_fill_manual(
    values = c(
      "Dif"  = "#457B9D",
      "Mel"  = "#A8DADC",
      "Phag" = "#E63946",
      "BOM"  = "#6D6875",
      "Toll" = "#1D3557"
    ),
    name = "Genotype"
  ) +
  facet_wrap(~ Genotype, nrow = 1) +
  coord_cartesian(ylim = c(0, 45)) +
  labs(
    x       = "Infection combination\n(First infection / Second infection)",
    y       = "Predicted egg count",
    title   = "Stage 2: Fecundity response to first and second infection",
    subtitle = "Positive values = more eggs with re-infection; negative = fewer",
    caption = "Predicted values from negative binomial GLMM. Error bars show ±1 SE.\nSome values extend beyond plot range.\n* Phag DTI/DTI re-infection effect: p = 0.0047 (Bonferroni corrected)"
  ) +
  theme_classic() +
  theme(
    legend.position  = "none",
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 11),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 7),
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(size = 9, colour = "grey40"),
    panel.spacing    = unit(0.8, "lines")
  )

# Calculate WT delta values
wt_delta <- plot_data_s2 %>%
  filter(line_type == "WT") %>%
  select(Infection1, wt_delta = delta)

# Subtract WT delta from each genotype
plot_data_s2_relative <- plot_data_s2 %>%
  filter(line_type == "Genotype") %>%
  left_join(wt_delta, by = "Infection1") %>%
  mutate(delta_vs_wt = delta - wt_delta)

# Calculate WT delta values
wt_delta <- plot_data_s2 %>%
  filter(line_type == "WT") %>%
  select(Infection1, wt_delta = delta)

# Subtract WT delta from each genotype
plot_data_s2_relative <- plot_data_s2 %>%
  filter(line_type == "Genotype") %>%
  left_join(wt_delta, by = "Infection1") %>%
  mutate(delta_vs_wt = delta - wt_delta)

# Significance annotation — Phag DTI only
sig_ann_relative <- plot_data_s2_relative %>%
  filter(Genotype == "Phag", Infection1 == "DTI") %>%
  mutate(label = "*")

fig2_relative <- ggplot(plot_data_s2_relative,
                        aes(x = Infection1, y = delta_vs_wt,
                            fill = Genotype)) +
  geom_hline(yintercept = 0, colour = "grey40",
             linewidth = 0.6, linetype = "dashed") +
  geom_col(width = 0.6, alpha = 0.85) +
  geom_text(data = sig_ann_relative,
            aes(x = Infection1, 
                y = 150,
                label = label),
            inherit.aes = FALSE,
            colour = "#E63946",
            size = 8,
            vjust = 0) +
  scale_fill_manual(
    values = c(
      "Dif"  = "#457B9D",
      "Mel"  = "#A8DADC",
      "Phag" = "#E63946",
      "BOM"  = "#6D6875",
      "Toll" = "#1D3557"
    ),
    name = "Genotype"
  ) +
  scale_x_discrete(labels = c("Control", "3-Tap", "DTI")) +
  facet_wrap(~ Genotype, nrow = 1) +
  coord_cartesian(ylim = c(-120, 165)) +
  labs(
    x        = "First infection treatment",
    y        = "Re-infection response relative to WT\n(egg count change − WT egg count change)",
    title    = "Stage 2: Genotype-specific re-infection responses relative to WT",
    subtitle = "Positive = greater egg count increase with re-infection than WT\nNegative = lesser increase than WT\nDashed line = WT response",
    caption  = "Predicted values from negative binomial GLMM.\n* Phag DTI significantly different from WT re-infection response: p = 0.0047 (Bonferroni corrected)"
  ) +
  theme_classic() +
  theme(
    legend.position  = "none",
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 11),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(size = 9, colour = "grey40"),
    panel.spacing    = unit(0.8, "lines"))

print(fig2_relative)
ggsave("Figure2_relative_to_WT.png", fig2_relative,
       width = 14, height = 5, dpi = 300)

# Extract genotype emmeans
emm_via_plot <- as.data.frame(summary(emm_v1_geno)) %>%
  filter(!(Genotype %in% c("DifCon", "IMD130"))) %>%  # remove separation-affected
  arrange(prob) %>%
  mutate(
    Genotype = factor(Genotype, levels = Genotype),
    category = case_when(
      Genotype == "WT"   ~ "WT",
      Genotype == "Phag" ~ "Significant",
      TRUE               ~ "Other"
    )
  )

# Step 1 — creating plot data with new labels 
plot_data_s1_final <- emm_s1_sig %>%
  filter(Genotype %in% c("WT", "Dif", "Mel", "Phag", "BOM", "Toll")) %>%
  mutate(Genotype = factor(Genotype,
                           levels = c("WT", "Dif", "Mel", "Phag", "BOM", "Toll"),
                           labels = c("WT", "ΔDif", "ΔMel", "ΔPhag", "ΔBOM", "ΔToll")))

# Step 2 — creating sig_ann, using new labels
sig_ann_s1 <- data.frame(
  Genotype = c("ΔMel", "ΔDif", "ΔToll", "ΔPhag"),
  Infection1 = c("DTI", "DTI", "DTI", "DTI"),
  label = c("*", "*", "*", "*")
) %>%
  left_join(
    plot_data_s1_final %>% 
      select(Genotype, Infection1, response, SE),
    by = c("Genotype", "Infection1")
  ) %>%
  mutate(Genotype = factor(Genotype,
                           levels = c("WT", "ΔDif", "ΔMel", "ΔPhag", "ΔBOM", "ΔToll")))

# Step 3 — plot with updated fill names
fig1_bars <- ggplot(plot_data_s1_final,
                    aes(x = Infection1, y = response, fill = Genotype)) +
  geom_col(width = 0.6, alpha = 0.85) +
  geom_errorbar(aes(ymin = response - SE, ymax = response + SE),
                width = 0.2, linewidth = 0.7) +
  geom_text(data = sig_ann_s1,
            aes(x = Infection1, y = response + SE + 1.5, label = label),
            inherit.aes = FALSE,
            size = 6,
            colour = "black") +
  scale_fill_manual(
    values = c(
      "WT"   = "#999999",
      "ΔDif" = "#457B9D",
      "ΔMel" = "#A8DADC",
      "ΔPhag" = "#E63946",
      "ΔBOM" = "#6D6875",
      "ΔToll" = "#1D3557"
    ),
    name = "Genotype"
  ) +
  scale_x_discrete(labels = c("Control", "3-Tap", "DTI")) +
  coord_cartesian(ylim = c(0, 50)) +
  facet_wrap(~ Genotype, nrow = 1) +
  labs(
    x       = "First infection treatment",
    y       = "Predicted egg count",
    caption = "Predicted values from negative binomial GLMM. Error bars show ±1 SE.\nSome values extend beyond plot range.\n* Significantly lower than WT at DTI (Bonferroni corrected)"
  ) +
  theme_classic() +
  theme(
    legend.position  = "none",
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 13),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(size = 9, colour = "grey40"),
    panel.spacing    = unit(1.5, "lines")
  )

print(fig1_bars)
ggsave("Figure1_stage1_final.png", fig1_bars,
       width = 14, height = 5, dpi = 300)
# =============================================================================
# Combined Stage 2 figure showing all infection dose combinations

# =============================================================================

# Calculate WT delta values
wt_delta <- plot_data_s2 %>%
  filter(line_type == "WT") %>%
  select(Infection1, wt_delta = delta)

# Subtract WT delta from each genotype
plot_data_s2_relative <- plot_data_s2 %>%
  filter(line_type == "Genotype") %>%
  left_join(wt_delta, by = "Infection1") %>%
  mutate(delta_vs_wt = delta - wt_delta)
# Calculate WT delta values
wt_delta <- plot_data_s2 %>%
  filter(line_type == "WT") %>%
  select(Infection1, wt_delta = delta)

# Subtract WT delta from each genotype and recode labels
plot_data_s2_relative <- plot_data_s2 %>%
  filter(line_type == "Genotype") %>%
  left_join(wt_delta, by = "Infection1") %>%
  mutate(
    delta_vs_wt = delta - wt_delta,
    Genotype = factor(Genotype,
                      levels = c("Dif", "Mel", "Phag", "BOM", "Toll"),
                      labels = c("ΔDif", "ΔMel", "ΔPhag", "ΔBOM", "ΔToll"))
  )

# Significance annotation — use new label
sig_ann_relative <- plot_data_s2_relative %>%
  filter(Genotype == "ΔPhag", Infection1 == "DTI") %>%
  mutate(label = "*")

# Significance annotation — Phag DTI only
sig_ann_relative <- plot_data_s2_relative %>%
  filter(Genotype == "ΔPhag", Infection1 == "DTI") %>%
  mutate(label = "*")

fig2_relative <- ggplot(plot_data_s2_relative,
                        aes(x = Infection1, y = delta_vs_wt,
                            fill = Genotype)) +
  geom_hline(yintercept = 0, colour = "grey40",
             linewidth = 0.6, linetype = "dashed") +
  geom_col(width = 0.6, alpha = 0.85) +
  geom_text(data = sig_ann_relative,
            aes(x = Infection1, 
                y = 150,  # fixed position near top of visible plot
                label = label),
            inherit.aes = FALSE,
            colour = "#E63946",
            size = 8,
            vjust = 0) +
  coord_cartesian(ylim = c(-120, 165)) +
  scale_fill_manual(
    values = c(
      "ΔDif"  = "#457B9D",
      "ΔMel"  = "#A8DADC",
      "ΔPhag" = "#E63946",
      "ΔBOM"  = "#6D6875",
      "ΔToll" = "#1D3557"
    ),
    name = "Genotype"
  ) +
  scale_x_discrete(labels = c("Control", "3-Tap", "DTI")) +
  facet_wrap(~ Genotype, nrow = 1) +
  labs(
    x        = "First infection treatment",
    y        = "Re-infection response relative to WT\n(egg count change − WT egg count change)",
    subtitle = "Positive = greater egg count increase with re-infection than WT\nNegative = lesser increase than WT\nDashed line = WT response",
    caption  = "Predicted values from negative binomial GLMM.\n* ΔPhag DTI significantly different from WT re-infection response: p = 0.0047 (Bonferroni corrected)"
  ) +
  theme_classic() +
  theme(
    legend.position  = "none",
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 11),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(size = 9, colour = "grey40"),
    panel.spacing    = unit(0.8, "lines")
  )

print(fig2_relative)
ggsave("Figure2_relative_to_WT.png", fig2_relative,
       width = 14, height = 5, dpi = 300)


#Phag vs WT Stage 2 Fecundity Figure
emm_phag <- emmeans(fec_phag_wt, ~ Infection1 * Infection2 * Genotype, 
                    type = "response")

emm_phag_plot <- as.data.frame(summary(emm_phag)) %>%
  mutate(
    Infection_combo = factor(paste(Infection1, Infection2, sep = "\n"),
                             levels = c("Control\nControl", "3 Tap\nControl", 
                                        "DTI\nControl", "Control\nDTI", 
                                        "3 Tap\nDTI", "DTI\nDTI")),
    line_type = ifelse(Genotype == "WT", "WT", "Genotype"),
    Genotype = factor(Genotype, 
                      levels = c("WT", "Phag"),
                      labels = c("WT", "ΔPhag")))

fig_phag_raw <- ggplot(emm_phag_plot,
                       aes(x = Infection_combo, y = response,
                           fill = Genotype, group = Genotype)) +
  geom_col(position = position_dodge(0.7), width = 0.6, alpha = 0.85) +
  geom_errorbar(aes(ymin = response - SE, ymax = response + SE),
                position = position_dodge(0.7), 
                width = 0.2, linewidth = 0.7) +
  scale_fill_manual(values = c("ΔPhag" = "#E63946", "WT" = "#999999")) +
  geom_vline(xintercept = 3.5, linetype = "dotted",
             colour = "grey60", linewidth = 0.5) +
  annotate("text", x = 2, y = max(emm_phag_plot$response + emm_phag_plot$SE, 
                                  na.rm = TRUE) * 0.95,
           label = "No re-infection", size = 3, colour = "grey40") +
  annotate("text", x = 5, y = max(emm_phag_plot$response + emm_phag_plot$SE, 
                                  na.rm = TRUE) * 0.95,
           label = "Re-infected", size = 3, colour = "grey40") +
  annotate("text", x = 6.15, y = max(emm_phag_plot$response + emm_phag_plot$SE, 
                                    na.rm = TRUE) * 1.15,
           label = "*", size = 8, colour = "#E63946") +
  coord_cartesian(ylim = c(0, 50)) +
  labs(
    x       = "Infection combination\n(First infection / Second infection)",
    y       = "Predicted egg count",
    caption = "Predicted values from focused negative binomial GLMM (WT and ΔPhag only). Error bars show ±1 SE.\n* DTI/DTI re-infection effect in ΔPhag: p = 0.0071 (Bonferroni corrected)"
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    plot.title      = element_text(face = "bold"),
    axis.text.x     = element_text(size = 8))

print(fig_phag_raw)

#--------------------------------------------------------------------
#--------------------------------------------------------------------
#stage 3 viability figure
# Extract genotype emmeans
emm_via_plot <- as.data.frame(summary(emm_v1_geno)) %>%
  filter(!(Genotype %in% c("DifCon", "IMD130"))) %>%
  arrange(prob) %>%
  mutate(
    Genotype = recode(Genotype,
                      "Phag"  = "ΔPhag",
                      "Mel"   = "ΔMel",
                      "Dif"   = "ΔDif",
                      "Toll"  = "ΔToll",
                      "Toll2" = "ΔToll2",
                      "BOM"   = "ΔBOM",
                      "IMD"   = "ΔIMD",
                      "AMP14" = "ΔAMP14",
                      "ITPM"  = "ΔITPM",
                      "Phag2" = "ΔPhag2"),
    Genotype = factor(Genotype, levels = Genotype),
    category = case_when(
      Genotype == "WT"    ~ "WT",
      Genotype == "ΔPhag" ~ "Significant",
      TRUE                ~ "Other"
    )
  )

fig_via <- ggplot(emm_via_plot,
                  aes(x = Genotype, y = prob, colour = category)) +
  # WT reference line
  geom_hline(
    yintercept = emm_via_plot$prob[emm_via_plot$Genotype == "WT"],
    linetype = "dashed", colour = "#1D3557", linewidth = 0.6
  ) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                width = 0.3, linewidth = 0.7) +
  geom_point(size = 4) +
  # Significance label above Phag
  geom_text(
    data = emm_via_plot %>% filter(Genotype == "ΔPhag"),
    aes(label = "*"), 
    colour = "#E63946", 
    size = 7,
    vjust = -4
  ) +
  scale_colour_manual(
    values = c(
      "WT"          = "#1D3557",
      "Significant" = "#E63946",
      "Other"       = "#999999"
    ),
    guide = "none"
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  labs(
    x       = "Genotype",
    y       = "Predicted hatching probability",
    subtitle = "Dashed line = WT reference. * ΔPhag significantly lower than WT (p = 0.033)",
    caption  = "Predicted values from zero-inflated beta-binomial GLMM, averaged across infection treatments.\nError bars show 95% CIs. DifCon and ΔIMD130 excluded (complete separation)."
  ) +
  theme_classic() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    plot.title   = element_text(face = "bold"),
    plot.subtitle = element_text(size = 9, colour = "grey40"))

print(fig_via)

ggsave("Figure_viability_genotype.png", fig_via,
       width = 8, height = 5, dpi = 300)