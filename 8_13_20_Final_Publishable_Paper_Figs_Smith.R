#### Load packages and data ####
install.packages("lme4")
install.packages("lmerTest")
install.packages("MuMIn")
install.packages("plot3D")

library(cowplot)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(Hmisc)
library(lme4)
library(lmerTest)
library(MuMIn)
library(plot3D)

# Abbreviate a binomial e.g. Balaenoptera musculus -> B. musculus
abbr_binom <- function(binom) {
  paste(str_sub(binom, 1, 1), 
        str_extract(binom, " .*"), 
        sep = ".")
}

#### Shirel's Allometric Eqs ####
# creating fucntions from Shirel's paper for MW (in kg) for engulfment capacity in liters for each species where we have a known length
Mass_SKR <- tribble(
  ~Species, ~slope,   ~intercept,
  "Minke",     1.8454,  2.1399,
  "Fin",     2.595,  1.28388,
  "Blue",     3.3346, 0.3129,
  "Humpback",     2.3373,  1.8535,
  "Bryde's",      2.7429, 1.0765,
  "Sei",    2.4238, 1.4772
)

#### Flukebeat and Morphometric Data #### 
morphometrics <- read_csv("Finalized Data Sheet For Hayden.csv") %>% 
  filter(DeployID != "bs190322-49") %>% # Removed because mean speed is below 1 m/s
  filter(DeployID != "mn180302-47") %>% # Removed because we don't know when it is lunging ("Unknown" for MaxOrNormal)
  filter(DeployID != "bw170816-41") # Removed because drag was too high

#All Swimming Flukebeat Info
d_all_swimming <- read_csv("AllDronedFlukebeatsFinalized.csv") %>%
  left_join(Mass_SKR, by = "Species") %>% 
  left_join(select(morphometrics, DeployID, c(FinenessRatio, SurfArea, ChordLength)), by = "DeployID") %>% 
  filter(DeployID != "bs190322-49") %>% 
  filter(DeployID != "mn180302-47") %>% 
  filter(DeployID != "bw170816-41") %>% 
  mutate(Mass = (TotLength^slope)*10^intercept,
         TPM = Thrust/Mass,
         Speed_BL = AvgSpeeds/TotLength,
         FA_L = FlukeArea/TotLength,
         FA_FR = FlukeArea/FinenessRatio)

d_all_swimming_summarized <- d_all_swimming %>%  
  group_by(DeployID) %>% 
  summarise(mean_freq = mean(OsFreq),
            sd_freq = sd(OsFreq),
            se_freq = sd_freq / sqrt(n()),
            mean_thrust = mean(Thrust),
            mean_TPM = mean(TPM), 
            sd_TPM = sd(TPM),
            se_TPM = sd_TPM / sqrt(n()),
            mean_drag = mean(`DragCoeff`),
            sd_drag = sd(`DragCoeff`),
            se_drag = sd_drag / sqrt(n()),
            mean_Re = mean(`Reynolds`),
            sd_Re = sd(`Reynolds`),
            se_Re = sd_Re / sqrt(n()),
            mean_E = mean(Eff),
            sd_E = sd(Eff),
            se_E = sd_E / sqrt(n()),
            mean_speed = mean(AvgSpeeds),
            sd_speed = sd(AvgSpeeds),
            se_speed = sd_speed / sqrt(n()),
            Species = first(Species),
            Length = first(TotLength),
            FinenessRatio = first(FinenessRatio),
            Mass = first(Mass),
            SurfArea = first(SurfArea),
            ChordLength = first(ChordLength),
            FA = first(FlukeArea),
            FA_L = first(FA_L),
            FA_FR = first(FA_FR)) %>% 
  mutate(ReynoldsModel = morphometrics$`Reynolds Number Model`,
         DragCoeffModel = morphometrics$`drag coefficient C_D Model`)

d_max_swimming <- d_all_swimming %>% 
  filter(MaxOrNormal == "Lunge-Associated")

d_max_swimming_summarized <- d_max_swimming %>%  
  group_by(DeployID) %>% 
  summarise(mean_freq = mean(OsFreq),
            sd_freq = sd(OsFreq),
            se_freq = sd_freq / sqrt(n()),
            mean_TPM = mean(TPM), 
            sd_TPM = sd(TPM),
            se_TPM = sd_TPM / sqrt(n()),
            mean_drag = mean(`DragCoeff`),
            sd_drag = sd(`DragCoeff`),
            se_drag = sd_drag / sqrt(n()),
            mean_Re = mean(`Reynolds`),
            sd_Re = sd(`Reynolds`),
            se_Re = sd_Re / sqrt(n()),
            mean_E = mean(Eff),
            sd_E = sd(Eff),
            se_E = sd_E / sqrt(n()),
            mean_speed = mean(AvgSpeeds),
            sd_speed = sd(AvgSpeeds),
            se_speed = sd_speed / sqrt(n()),
            Species = first(Species),
            Length = first(TotLength),
            Effort = first(MaxOrNormal),
            FinenessRatio = first(FinenessRatio),
            Mass = first(Mass),
            SurfArea = first(SurfArea),
            ChordLength = first(ChordLength),
            FA = first(FlukeArea),
            FA_L = first(FA_L),
            FA_FR = first(FA_FR)) %>%  
  mutate(ReynoldsModel = morphometrics$`Reynolds Number Model`, 
         DragCoeffModel = morphometrics$`drag coefficient C_D Model`)

d_routine_swimming <- d_all_swimming %>%
  filter(MaxOrNormal == "Routine")

d_routine_swimming_summarized <- d_routine_swimming %>%  
  group_by(DeployID) %>% 
  summarise(mean_freq = mean(OsFreq),
            sd_freq = sd(OsFreq),
            se_freq = sd_freq / sqrt(n()),
            mean_TPM = mean(TPM), 
            sd_TPM = sd(TPM),
            se_TPM = sd_TPM / sqrt(n()),
            mean_drag = mean(`DragCoeff`),
            sd_drag = sd(`DragCoeff`),
            se_drag = sd_drag / sqrt(n()),
            mean_Re = mean(`Reynolds`),
            sd_Re = sd(`Reynolds`),
            se_Re = sd_Re / sqrt(n()),
            mean_E = mean(Eff),
            sd_E = sd(Eff),
            se_E = sd_E / sqrt(n()),
            mean_speed = mean(AvgSpeeds),
            sd_speed = sd(AvgSpeeds),
            se_speed = sd_speed / sqrt(n()),
            Species = first(Species),
            Length = first(TotLength),
            Effort = first(MaxOrNormal),
            FinenessRatio = first(FinenessRatio),
            Mass = first(Mass),
            SurfArea = first(SurfArea),
            ChordLength = first(ChordLength),
            FA = first(FlukeArea),
            FA_L = first(FA_L),
            FA_FR = first(FA_FR)) %>% 
  mutate(ReynoldsModel = morphometrics$`Reynolds Number Model`,
         DragCoeffModel = morphometrics$`drag coefficient C_D Model`)

Deltas <- select(d_routine_swimming_summarized, DeployID, mean_speed_routine = mean_speed, mean_TPM_routine = mean_TPM, mean_drag_routine = mean_drag)  %>% 
  left_join(select(d_max_swimming_summarized, DeployID, mean_speed_max = mean_speed, mean_TPM_max = mean_TPM, mean_drag_max = mean_drag),
            by = "DeployID") %>% 
  mutate(DeltaU = mean_speed_max - mean_speed_routine,
         DeltaDrag = mean_drag_max - mean_drag_routine,
         DeltaTPM = mean_TPM_max - mean_TPM_routine,
         DeltaTPM2 = DeltaTPM^2,
         DeltaTPM3 = DeltaTPM^3)

d_all_swimming_summarized <- d_all_swimming_summarized %>% 
  left_join(select(Deltas, DeployID, DeltaU, DeltaDrag, DeltaTPM, DeltaTPM2, DeltaTPM3),
            by = "DeployID") %>% 
  left_join(select(d_routine_swimming_summarized, DeployID, mean_Re_routine = mean_Re, mean_drag_routine = mean_drag),
            by = "DeployID") %>% 
  left_join(select(d_max_swimming_summarized, DeployID, mean_Re_max = mean_Re, mean_drag_max = mean_drag),
            by = "DeployID") %>% 
  mutate(DeltaDragOverU = DeltaDrag/DeltaU,
         DiffDragDragModelRoutine = mean_drag_routine - DragCoeffModel,
         DiffDragDragModelMax = mean_drag_max - DragCoeffModel)

d_routine_Sp_Sum <- d_routine_swimming_summarized %>% 
  group_by(Species) %>% 
  summarise(sum_freq = mean(mean_freq),
            sumsd_freq = sd(mean_freq),
            sumse_freq = sumsd_freq / sqrt(n()),
            sum_TPM = mean(mean_TPM), 
            sumsd_TPM = sd(mean_TPM),
            sumse_TPM = sumsd_TPM / sqrt(n()),
            sum_drag = mean(mean_drag),
            sumsd_drag = sd(mean_drag),
            sumse_drag = sumsd_drag / sqrt(n()),
            sum_Re = mean(mean_Re),
            sumsd_Re = sd(mean_Re),
            sumse_Re = sumsd_Re / sqrt(n()),
            sum_E = mean(mean_E),
            sumsd_E = sd(mean_E),
            sumse_E = sumsd_E / sqrt(n()),
            sum_speed = mean(mean_speed),
            sumsd_speed = sd(mean_speed),
            sumse_speed = sumsd_speed / sqrt(n()),
            sum_Length = mean(Length),
            sumsd_Length = sd(Length),
            sumse_Length = sumsd_Length / sqrt(n()),
            sum_Fineness = mean(FinenessRatio, na.rm=TRUE),
            sumsd_Fineness = sd(FinenessRatio, na.rm=TRUE),
            sumse_Fineness = sumsd_Fineness / sqrt(n()),
            sum_Mass = mean(Mass),
            sumsd_Mass = sd(Mass),
            sumse_Mass = sumsd_Mass / sqrt(n()),
            sum_SA = mean(SurfArea),
            sumsd_SA = sd(SurfArea),
            sumse_SA = sumsd_SA / sqrt(n()),
            sum_ChordLength = mean(ChordLength),
            sumsd_ChordLength = sd(ChordLength),
            sumse_ChordLength = sumsd_ChordLength / sqrt(n()),
            sum_FA = mean(FA),
            sumsd_FA = sd(FA),
            sumse_FA = sumsd_FA / sqrt(n()),
            sum_FA_L = mean(FA_L),
            sumsd_FA_L = sd(FA_L),
            sumse_FA_L = sumsd_FA_L / sqrt(n()),
            Effort = first(Effort))

d_max_Sp_Sum <- d_max_swimming_summarized %>% 
  group_by(Species) %>% 
  summarise(sum_freq = mean(mean_freq),
            sumsd_freq = sd(mean_freq),
            sumse_freq = sumsd_freq / sqrt(n()),
            sum_TPM = mean(mean_TPM), 
            sumsd_TPM = sd(mean_TPM),
            sumse_TPM = sumsd_TPM / sqrt(n()),
            sum_drag = mean(mean_drag),
            sumsd_drag = sd(mean_drag),
            sumse_drag = sumsd_drag / sqrt(n()),
            sum_Re = mean(mean_Re),
            sumsd_Re = sd(mean_Re),
            sumse_Re = sumsd_Re / sqrt(n()),
            sum_E = mean(mean_E),
            sumsd_E = sd(mean_E),
            sumse_E = sumsd_E / sqrt(n()),
            sum_speed = mean(mean_speed),
            sumsd_speed = sd(mean_speed),
            sumse_speed = sumsd_speed / sqrt(n()),
            sum_Length = mean(Length),
            sumsd_Length = sd(Length),
            sumse_Length = sumsd_Length / sqrt(n()),
            sum_Fineness = mean(FinenessRatio, na.rm=TRUE),
            sumsd_Fineness = sd(FinenessRatio, na.rm=TRUE),
            sumse_Fineness = sumsd_Fineness / sqrt(n()),
            sum_Mass = mean(Mass),
            sumsd_Mass = sd(Mass),
            sumse_Mass = sumsd_Mass / sqrt(n()),
            sum_SA = mean(SurfArea),
            sumsd_SA = sd(SurfArea),
            sumse_SA = sumsd_SA / sqrt(n()),
            sum_ChordLength = mean(ChordLength),
            sumsd_ChordLength = sd(ChordLength),
            sumse_ChordLength = sumsd_ChordLength / sqrt(n()),
            sum_FA = mean(FA),
            sumsd_FA = sd(FA),
            sumse_FA = sumsd_FA / sqrt(n()),
            sum_FA_L = mean(FA_L),
            sumsd_FA_L = sd(FA_L),
            sumse_FA_L = sumsd_FA_L / sqrt(n()),
            Effort = first(Effort))

d_all_binded_summarized <- d_routine_swimming_summarized %>%
  bind_rows(d_max_swimming_summarized)

d_all_spec_binded_summarized <- d_routine_Sp_Sum %>%
  bind_rows(d_max_Sp_Sum)

d_routine_nums <- count(d_routine_swimming_summarized, Species)
d_routine_nums
d_max_nums <- count(d_max_swimming_summarized, Species)
d_max_nums

#### Color Palette #### - look for color blind pallete 
pal <- c("Minke" = "#F28E2B",  "Humpback" = "#59A14F",  "Blue" = "#4E79A7", "Sei" = "#E15759", "Fin" = "#B07AA1", "Bryde's" = "#17BECF", 'Normal' = "Black", 'Lunge-Associated' = "Black")
pal2 <- c("Human" = "#59A14F", "Fish" = "#E15759", "Pinniped" = "#79706E", "Sirenian" = "#B6992D", "Odontocete" = "#B07AA1", "Mysticete" = "#4E79A7", "Rodent" = "#9467BD")



#### Graphs Start Here ####



#### Length ~ Freq ####
fig3Freq <- ggplot(d_routine_swimming_summarized, aes(log(Length), log(mean_freq))) +
  geom_smooth(method = "lm", color = "black", size = 3) +
  geom_smooth(data = d_max_swimming_summarized, aes(log(Length), log(mean_freq)), method = "lm", color = "black", size = 3, linetype = "longdash") +
  geom_point(aes(fill = Species), size = 10, shape = 21) +
  geom_point(data = d_max_swimming_summarized, aes(log(Length), log(mean_freq), fill = Species), size = 10, shape = 24) +
  scale_fill_manual(values = pal) +
  ylim(-2,-0.5) +
  labs(x = bquote('log10 Total Length (m)'),
       y = bquote('log10 Oscillatory Frequency (Hz)')) +
  theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 40),
        axis.title = element_text(size = 48),
        legend.position = "none",
        panel.grid.minor = element_blank())
fig3Freq

# Generalized linear mixed models
GLMM3Freqmax_mean <- lmer(log(mean_freq) ~ log(Length) + (1|Species), 
                      data = d_max_swimming_summarized)
summary(GLMM3Freqmax_mean)
r.squaredGLMM(GLMM3Freqmax_mean)

GLMM3Freqnormal_mean <- lmer(log(mean_freq) ~ log(Length) + (1|Species), 
                         data = d_routine_swimming_summarized)
summary(GLMM3Freqnormal_mean)
r.squaredGLMM(GLMM3Freqnormal_mean)

#### Length ~ Speed (Routine) ####
fig3U <- ggplot(d_routine_swimming_summarized, aes(log(Length), log(mean_speed))) +
  geom_smooth(method = "lm", color = "black", size = 3) +
  geom_smooth(data = d_max_swimming_summarized, aes(log(Length), log(mean_speed)), method = "lm", color = "black", size = 3, linetype = "longdash") +
  geom_point(aes(fill = Species), size = 10, shape = 21) +
  geom_point(data = d_max_swimming_summarized, aes(log(Length), log(mean_speed), fill = Species), size = 10, shape = 24) +
  scale_fill_manual(values = pal) +
  labs(x = bquote('log10 Total Length (m)'),
       y = bquote('log10 Swim Speed'~(m~s^-1))) +
  theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 40),
        axis.title = element_text(size = 48),
        legend.position = "none",
        panel.grid.minor = element_blank())
fig3U

# Generalized linear mixed models
GLMM3Umax_mean <- lmer(log(mean_speed) ~ log(Length) + (1|Species), 
                          data = d_max_swimming_summarized)
summary(GLMM3Umax_mean)
r.squaredGLMM(GLMM3Umax_mean)

GLMM3Unormal_mean <- lmer(log(mean_speed) ~ log(Length) + (1|Species), 
                             data = d_routine_swimming_summarized)
summary(GLMM3Unormal_mean)
r.squaredGLMM(GLMM3Unormal_mean)

# Combine two figues into one
fig3 <- plot_grid(fig3Freq, fig3U,
                  nrow = 1,
                  align = "h",
                  labels = NULL)
ggsave("Figures/fig3.pdf", height = 480, width = 960, units = "mm", dpi = 300)
fig3



#### MST ~ U (Routine) ####
fig4U <- ggplot(d_routine_swimming_summarized, aes(mean_speed, mean_TPM)) +
  geom_smooth(method = "lm", color = "black", size = 3) +
  geom_smooth(data = d_max_swimming_summarized, aes(mean_speed, mean_TPM), method = "lm", color = "black", size = 3, linetype = "longdash") +
  geom_point(aes(fill = Species), size = 10, shape = 21) +
  geom_point(data = d_max_swimming_summarized, aes(mean_speed, mean_TPM, fill = Species), size = 10, shape = 24) +
  scale_fill_manual(values = pal) +
  ylim(0,2) +
  labs(x = bquote('Swim Speed'~(m~s^-1)),
       y = bquote('Mass-Specific Thrust Power'~(W~kg^-1))) +
  theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 40),
        axis.title = element_text(size = 48),
        legend.position = "none",
        panel.grid.minor = element_blank())
fig4U

# Generalized linear mixed models
GLMM4Umax_mean <- lmer(log(mean_TPM) ~ mean_speed + (1|Species), 
                       data = d_max_swimming_summarized)
summary(GLMM4Umax_mean)
r.squaredGLMM(GLMM4Umax_mean)

GLMM4Unormal_mean <- lmer(log(mean_TPM) ~ mean_speed + (1|Species), 
                          data = d_routine_swimming_summarized)
summary(GLMM4Unormal_mean)
r.squaredGLMM(GLMM4Unormal_mean)

#### MST ~ L (Routine vs. Lunge-Associated) ####

fig4TL <- ggplot(d_routine_swimming_summarized, aes(Length, mean_TPM)) +
  geom_smooth(method = "lm", color = "black", size = 3) +
  geom_smooth(data = d_max_swimming_summarized, aes(Length, mean_TPM), method = "lm", color = "black", size = 3, linetype = "longdash") +
  geom_point(aes(fill = Species), size = 10, shape = 21) +
  geom_point(data = d_max_swimming_summarized, aes(Length, mean_TPM, fill = Species), size = 10, shape = 24) +
  scale_fill_manual(values = pal) +
  ylim(0,2) +
  labs(x = bquote('Total Length (m)'),
       y = bquote('Mass-Specific Thrust Power'~(W~kg^-1))) +
  theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 40),
        axis.title = element_text(size = 48),
        legend.position = "none",
        panel.grid.minor = element_blank())
fig4TL

# Generalized linear mixed models
GLMM4TLmax_mean <- lmer(log(mean_TPM) ~ Length + (1|Species), 
                       data = d_max_swimming_summarized)
summary(GLMM4TLmax_mean)
r.squaredGLMM(GLMM4TLmax_mean)

GLMM4TLnormal_mean <- lmer(log(mean_TPM) ~ Length + (1|Species), 
                          data = d_routine_swimming_summarized)
summary(GLMM4TLnormal_mean)
r.squaredGLMM(GLMM4TLnormal_mean)

# Combine two figues into one
fig4 <- plot_grid(fig4U, fig4TL,
                  nrow = 1,
                  align = "h",
                  labels = NULL)
ggsave("Figures/fig4.pdf", height = 480, width = 960, units = "mm", dpi = 300)
fig4


#### Drag ~ U (Routine) ####
fig5U <- ggplot(d_routine_swimming_summarized, aes(mean_speed, mean_drag)) +
  geom_smooth(method = "lm", color = "black", size = 3) +
  geom_smooth(data = d_max_swimming_summarized, aes(mean_speed, mean_drag), method = "lm", color = "black", size = 3, linetype = "longdash") +
  geom_point(aes(fill = Species), size = 10, shape = 21) +
  geom_point(data = d_max_swimming_summarized, aes(mean_speed, mean_drag, fill = Species), size = 10, shape = 24) +
  scale_fill_manual(values = pal) +
  ylim(0, 0.075) +
  labs(x = bquote('Swim Speed'~(m~s^-1)),
       y = bquote('Drag Coefficient')) +
  theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 40),
        axis.title = element_text(size = 48),
        legend.position = "none",
        panel.grid.minor = element_blank())
fig5U

# Generalized linear mixed models
GLMM5Umax_mean <- lmer(log(mean_drag) ~ mean_speed + (1|Species), 
                        data = d_max_swimming_summarized)
summary(GLMM5Umax_mean)
r.squaredGLMM(GLMM5Umax_mean)

GLMM5Unormal_mean <- lmer(log(mean_drag) ~ mean_speed + (1|Species), 
                           data = d_routine_swimming_summarized)
summary(GLMM5Unormal_mean)
r.squaredGLMM(GLMM5Unormal_mean)

#### Drag ~ Length (Routine) ####
fig5TL <- ggplot(d_routine_swimming_summarized, aes(Length, mean_drag)) +
  geom_smooth(method = "lm", color = "black", size = 3) +
  geom_smooth(data = d_max_swimming_summarized, aes(Length, mean_drag), method = "lm", color = "black", size = 3, linetype = "longdash") +
  geom_point(aes(fill = Species), size = 10, shape = 21) +
  geom_point(data = d_max_swimming_summarized, aes(Length, mean_drag, fill = Species), size = 10, shape = 24) +
  scale_fill_manual(values = pal) +
  ylim(0, 0.075) +
  labs(x = bquote('Total Length (m)'),
       y = bquote('Drag Coefficient')) +
  theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 40),
        axis.title.x = element_text(size = 48),
        axis.title.y = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_blank())
fig5TL

# Generalized linear mixed models
GLMM5TLmax_mean <- lmer(log(mean_drag) ~ Length + (1|Species), 
                       data = d_max_swimming_summarized)
summary(GLMM5TLmax_mean)
r.squaredGLMM(GLMM5TLmax_mean)

GLMM5TLnormal_mean <- lmer(log(mean_drag) ~ Length + (1|Species), 
                          data = d_routine_swimming_summarized)
summary(GLMM5TLnormal_mean)
r.squaredGLMM(GLMM5TLnormal_mean)

#### Drag ~ Re (Species Level + Hoerner Models) ####
fig5Re <- ggplot(d_all_spec_binded_summarized, aes(sum_Re, sum_drag)) +
  geom_errorbar(aes(ymax = sum_drag+sumsd_drag, ymin = sum_drag-sumsd_drag)) +
  geom_errorbarh(aes(xmax = sum_Re+sumsd_Re, xmin = sum_Re-sumsd_Re)) +
  geom_point(aes(shape = Effort, fill = Species), size = 10) +
  geom_smooth(data = d_routine_swimming_summarized, aes(mean_Re, DragCoeffModel), method = "lm", color = "black", linetype = 3, size = 3) +
  geom_smooth(data = d_max_swimming_summarized, aes(mean_Re, DragCoeffModel), method = "lm", color = "black", linetype = 4, size = 3) +
  scale_fill_manual(values = pal) +
  scale_shape_manual(values=c(24,21)) +
  ylim(0, 0.075) +
  labs(x = bquote('Reynolds Number'),
       y = bquote('Drag Coefficient')) +
  theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 40),
        axis.title = element_text(size = 48),
        legend.position = "none",
        panel.grid.minor = element_blank())
fig5Re

# Combine three figues into one
fig5 <- plot_grid(fig5U, fig5TL, fig5Re,
                  nrow = 2,
                  align = "h",
                  labels = NULL)
ggsave("Figures/fig5.pdf", height = 960, width = 960, units = "mm", dpi = 300)
fig5


#### Prop Efficiency ~ U, L ####
# Propulsive Efficiency Vs Speed
fig6U <- ggplot(d_routine_swimming_summarized, aes(mean_speed, mean_E)) +
  geom_smooth(method = "lm", color = "black", size = 3) +
  geom_smooth(data = d_max_swimming_summarized, aes(mean_speed, mean_E), method = "lm", color = "black", size = 3, linetype = "longdash") +
  geom_point(aes(fill = Species), size = 10, shape = 21) +
  geom_point(data = d_max_swimming_summarized, aes(mean_speed, mean_E, fill = Species), size = 10, shape = 24) +
  scale_fill_manual(values = pal) +
  labs(x = bquote('Swim Speed'~(m~s^-1)),
       y = "Propulsive Efficiency") +
  expand_limits(y = c(0.75, 1)) +
  theme_bw(base_size = 20, base_family = "Times") +
  theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 40),
        axis.title = element_text(size = 48),
        legend.position = "none",
        panel.grid.minor = element_blank())
fig6U

# Generalized linear mixed models
GLMM6Umax_mean <- lmer(log(mean_E) ~ mean_speed + (1|Species), 
                       data = d_max_swimming_summarized)
summary(GLMM6Umax_mean)
r.squaredGLMM(GLMM6Umax_mean)

GLMM6Unormal_mean <- lmer(log(mean_E) ~ mean_speed + (1|Species), 
                          data = d_routine_swimming_summarized)
summary(GLMM6Unormal_mean)
r.squaredGLMM(GLMM6Unormal_mean)

# Propulsive Efficiency Vs Length
fig6TL <- ggplot(d_routine_swimming_summarized, aes(Length, mean_E)) +
  geom_smooth(method = "lm", color = "black", size = 3) +
  geom_smooth(data = d_max_swimming_summarized, aes(Length, mean_E), method = "lm", color = "black", size = 3, linetype = "longdash") +
  geom_point(aes(fill = Species), size = 10, shape = 21) +
  geom_point(data = d_max_swimming_summarized, aes(Length, mean_E, fill = Species), size = 10, shape = 24) +
  scale_fill_manual(values = pal) +
  labs(x = bquote('Total Length'~(m)),
       y = bquote('Propulsive Efficiency')) +
  expand_limits(y = c(0.75, 1)) +
  theme_bw(base_size = 20, base_family = "Times") +
  theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 40),
        axis.title = element_text(size = 48),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank())
fig6TL

# Generalized linear mixed models
GLMM6TLmax_mean <- lm(log(mean_E) ~ Length, 
                       data = d_max_swimming_summarized)
summary(GLMM6TLmax_mean)
r.squaredGLMM(GLMM6TLmax_mean)

GLMM6TLnormal_mean <- lm(log(mean_E) ~ Length, 
                          data = d_routine_swimming_summarized)
summary(GLMM6TLnormal_mean)
r.squaredGLMM(GLMM6TLnormal_mean)

# Combine two figues into one
fig6 <- plot_grid(fig6U, fig6TL,
                  nrow = 1,
                  align = "h",
                  labels = NULL)
ggsave("Figures/fig6.pdf", height = 480, width = 960, units = "mm", dpi = 300)
fig6


#### Prop Eff ~ L w/ Other Species ####
prop_effs <- read_csv("Propulsive Eff All Species.csv")
fig7 <- ggplot(prop_effs, aes(`Total Length (m)`, `Prop Eff (Max)`)) +
  geom_point(aes(fill = Group, shape = `Type of Swimming`), size = 10) +
  scale_fill_manual(values = pal2) +
  scale_shape_manual(values = c(21, 24, 22)) +
  expand_limits(y = c(0, 1)) +
  labs(x = "Total Length (m)",
       y = "Propulsive Efficiency") +
  theme_classic() +
  theme(axis.text = element_text(size = 40),
        axis.title = element_text(size = 48),
        legend.position = "none",
        panel.grid.minor = element_blank())
ggsave("Figures/fig7.pdf", height = 480, width = 480, units = "mm", dpi = 300)
fig7



#### Extra Figures ####



#### Prop Eff Cetaceans Only ####
cet_prop_effs <- filter(prop_effs, Group %in% c("Odontocete", "Mysticete"))
fig8 <- ggplot(cet_prop_effs, aes(log10(`Total Length (m)`), `Prop Eff (Max)`)) +
  geom_smooth(aes(linetype = Group), method = "lm", size = 3, color = "black") +
  geom_point(aes(color = Group), size = 10) +
  scale_color_manual(values = pal2) +
  labs(x = "Log10 Total Length (m)",
       y = "Propulsive Efficiency") +
  theme_classic() +
  theme(axis.text = element_text(size = 40),
        axis.title = element_text(size = 48),
        legend.position = "none",
        panel.grid.minor = element_blank())
ggsave("Figures/fig8.pdf", height = 480, width = 480, units = "mm", dpi = 300)
fig8

odont_prop_effs <- filter(prop_effs, Group %in% "Odontocete")
mysti_prop_effs <- filter(prop_effs, Group %in% "Mysticete")

#### Delta U ~ Delta MST ####
figDeltas <-ggplot(d_all_swimming_summarized, aes(DeltaU, DeltaTPM)) +
  geom_point(aes(color = Species), size = 10) +
  geom_smooth(method = "lm", color = "black", size = 3) +
  scale_color_manual(values = pal) +
  expand_limits(y = 0) +
  labs(x = bquote('Delta Swim Speed'~(m~s^-1)),
       y = bquote('Delta Mass-Specific Thrust Power'~(Watts~kg^-1))) +
  theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 40),
        axis.title = element_text(size = 48),
        legend.position = "none",
        panel.grid.minor = element_blank())
ggsave("Figures/figDeltas.pdf", height = 480, width = 480, units = "mm", dpi = 300)
figDeltas

####  Fineness Ratio ~ TPM ####
figFR <- ggplot(d_routine_swimming_summarized, aes(FinenessRatio, mean_TPM)) +
  geom_point(aes(color = Species), size = 10) +
  geom_smooth(method = "lm", size = 3, color = "Black") +
  scale_color_manual(values = pal) +
  labs(x = bquote('Fineness Ratio'),
       y = bquote('Mean Mass-Specific Thrust Power'~(Watts~kg^-1))) +
  theme_minimal(base_size = 8) +
  theme_bw(base_size = 20, base_family = "Times") +
  theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 40),
        axis.title = element_text(size = 48),
        legend.position = "right",
        panel.grid.minor = element_blank())
ggsave("Figures/figFR.pdf", height = 480, width = 480, units = "mm", dpi = 300)
figFR

#### Fluke Area / Fineness Ratio ~ TPM ####
figFAFR <- ggplot(d_routine_swimming_summarized, aes(FA_FR, mean_TPM)) +
  geom_point(aes(color = Species), size = 10) +
  geom_smooth(method = "lm", size = 3, color = "Black") +
  scale_color_manual(values = pal) +
  labs(x = bquote('Fluke Area / Fineness Ratio'),
       y = bquote('Mean Mass-Specific Thrust Power'~(Watts~kg^-1))) +
  theme_minimal(base_size = 8) +
  theme_bw(base_size = 20, base_family = "Times") +
  theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 40),
        axis.title = element_text(size = 48),
        legend.position = "none",
        panel.grid.minor = element_blank())
ggsave("Figures/figFAFR.pdf", height = 480, width = 480, units = "mm", dpi = 300)
figFAFR

#### MST ~ Fluke Area (No Length) ####
figFA <- ggplot(d_routine_swimming_summarized, aes(FA, mean_TPM)) +
  geom_point(aes(color = Species), size = 10) +
  geom_smooth(method = "lm", size = 3, color = "Black") +
  scale_color_manual(values = pal) +
  labs(x = bquote('Fluke Area'),
       y = bquote('Mean Mass-Specific Thrust Power'~(Watts~kg^-1))) +
  theme_minimal(base_size = 8) +
  theme_bw(base_size = 20, base_family = "Times") +
  theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 40),
        axis.title = element_text(size = 48),
        legend.position = "none",
        panel.grid.minor = element_blank())
ggsave("Figures/figFA.pdf", height = 480, width = 480, units = "mm", dpi = 300)
figFA

#### MST ~ Freq (Routine) ####
figOsFreq <- ggplot(d_routine_swimming_summarized, aes(mean_freq, mean_TPM)) +
  geom_smooth(method = "lm", color = "black", size = 3) +
  geom_point(aes(color = Species), size = 10) +
  scale_color_manual(values = pal) +
  expand_limits(y = 0) +
  labs(x = bquote('Oscillatory Frequency (Hz)'),
       y = bquote('Mass-Specific Thrust Power'~(Watts~kg^-1))) +
  theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 40),
        axis.title = element_text(size = 48),
        legend.position = "right",
        panel.grid.minor = element_blank())
ggsave("Figures/figOsFreq.pdf", height = 480, width = 480, units = "mm", dpi = 300)
figOsFreq

#### U ~ Freq (Routine) ####
figOsFreqU <- ggplot(d_routine_swimming_summarized, aes(mean_speed, mean_freq)) +
  geom_smooth(method = "lm", color = "black", size = 3) +
  geom_point(aes(color = Species), size = 10) +
  scale_color_manual(values = pal) +
  expand_limits(y = 0) +
  labs(x = bquote('Swim Speed'~(m~s^-1)),
       y = bquote('Oscillatory Frequency (Hz)')) +
  theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 40),
        axis.title = element_text(size = 48),
        legend.position = "right",
        panel.grid.minor = element_blank())
ggsave("Figures/figOsFreqU.pdf", height = 480, width = 480, units = "mm", dpi = 300)
figOsFreqU

#### MST ~ Fluke Area ####
figFluTL <- ggplot(d_routine_swimming_summarized, aes(FA_L, mean_TPM)) +
  geom_point(aes(color = Species), size = 10) +
  geom_smooth(method = "lm", size = 3, color = "Black") +
  scale_color_manual(values = pal) +
  labs(x = bquote('Fluke Area / Total Length (m)'),
       y = bquote('Log(Mean Mass-Specific Thrust Power)'~(Watts~kg^-1))) +
  theme_minimal(base_size = 8) +
  theme_bw(base_size = 20, base_family = "Times") +
  theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 40),
        axis.title = element_text(size = 48),
        legend.position = "none",
        panel.grid.minor = element_blank())
ggsave("Figures/figFluTL.pdf", height = 480, width = 480, units = "mm", dpi = 300)
figFluTL

#### Drag ~ Re (Individuals + Hoerner Models) ####
fig5Re <- ggplot(d_all_swimming_summarized) +
  geom_point(aes(mean_Re, DiffDragDragModelRoutine, color = Species, size = 10)) +
  geom_smooth(aes(mean_Re, DiffDragDragModelRoutine), method = "lm")
fig5Re

# Generalized linear mixed models
GLMMTest <- lmer(log(mean_Re) ~ DiffDragDragModelRoutine + (1|Species), 
                        data = d_all_swimming_summarized)
summary(GLMMTest)
r.squaredGLMM(GLMMTest)
        

#### Drag ~ Re (+ Hoerner Models) ####
fig5Re <- ggplot(d_all_binded_summarized) +
  geom_smooth(method = "lm", aes(mean_Re, mean_drag), color = "black", size = 3) +
  geom_point(aes(mean_Re, mean_drag, color = Species), size = 10) +
  scale_color_manual(values = pal) +
  ylim(0, 0.075) +
  labs(x = bquote('Reynolds Number'),
       y = bquote('Drag Coefficient')) +
  theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 40),
        axis.title = element_text(size = 48),
        legend.position = "none",
        panel.grid.minor = element_blank())
fig5Re

#### Drag ~ Re (+ Hoerner Models) Bad Version ####
fig5Re <- ggplot(d_routine_swimming_summarized) +
  #geom_smooth(method = "lm", aes(mean_Re, mean_drag), color = "black", size = 3) +
  #geom_smooth(data = d_max_swimming_summarized, aes(mean_Re, mean_drag), method = "lm", color = "black", linetype = "longdash", size = 3) +
  geom_smooth(method = "lm", aes(mean_Re, DragCoeffModel), color = "black", linetype = 3, size = 3) +
  geom_smooth(data = d_max_swimming_summarized, aes(mean_Re, DragCoeffModel), method = "lm", color = "black", linetype = 4, size = 3) +
  geom_point(aes(mean_Re, mean_drag, color = Species), size = 10) +
  geom_point(data = d_max_swimming_summarized, aes(mean_Re, mean_drag, color = Species), size = 10, shape = 17) +
  scale_color_manual(values = pal) +
  ylim(0, 0.075) +
  labs(x = bquote('Reynolds Number'),
       y = bquote('Drag Coefficient')) +
  theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 40),
        axis.title = element_text(size = 48),
        legend.position = "none",
        panel.grid.minor = element_blank())
fig5Re
