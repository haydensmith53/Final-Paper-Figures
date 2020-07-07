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
  filter(DeployID != "bs190322-49")

morphometricsFilt <- morphometrics %>% # filter out mn180302-47 (no lunge-associated tailbeats)
  filter(DeployID != "mn180302-47")

#All Swimming Flukebeat Info
d_all_swimming <- read_csv("AllDronedFlukebeatsFinalized.csv") %>%
  left_join(Mass_SKR, by = "Species") %>% 
  left_join(select(morphometrics, DeployID, c(FinenessRatio, SurfArea, ChordLength)), by = "DeployID") %>% 
  filter(DeployID != "bs190322-49") %>% 
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
            FinenessRatio = first(FinenessRatio),
            Mass = first(Mass),
            SurfArea = first(SurfArea),
            ChordLength = first(ChordLength),
            FA = first(FlukeArea),
            FA_L = first(FA_L),
            FA_FR = first(FA_FR)) %>%  
  mutate(ReynoldsModel = morphometricsFilt$`Reynolds Number Model`, 
         DragCoeffModel = morphometricsFilt$`drag coefficient C_D Model`)

d_routine_swimming <- d_all_swimming %>%
  filter(MaxOrNormal == "Routine" || "Unknown")

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
            FinenessRatio = first(FinenessRatio),
            Mass = first(Mass),
            SurfArea = first(SurfArea),
            ChordLength = first(ChordLength),
            FA = first(FlukeArea),
            FA_L = first(FA_L),
            FA_FR = first(FA_FR)) %>% 
  mutate(ReynoldsModel = morphometrics$`Reynolds Number Model`,
         DragCoeffModel = morphometrics$`drag coefficient C_D Model`)

Deltas <- select(d_routine_swimming_summarized, DeployID, mean_speed_routine = mean_speed, mean_TPM_routine = mean_TPM)  %>% 
  left_join(select(d_max_swimming_summarized, DeployID, mean_speed_max = mean_speed, mean_TPM_max = mean_TPM),
            by = "DeployID") %>% 
  mutate(DeltaU = mean_speed_max - mean_speed_routine,
         DeltaTPM = mean_TPM_max - mean_TPM_routine,
         DeltaTPM2 = DeltaTPM^2,
         DeltaTPM3 = DeltaTPM^3)

d_all_swimming_summarized <- d_all_swimming_summarized %>% 
  left_join(select(Deltas, DeployID, DeltaU, DeltaTPM, DeltaTPM2, DeltaTPM3),
            by = "DeployID")

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
            sumse_FA_L = sumsd_FA_L / sqrt(n()))

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
            sumse_FA_L = sumsd_FA_L / sqrt(n()))

d_routine_nums <- count(d_routine_swimming_summarized, Species)
d_routine_nums
d_max_nums <- count(d_max_swimming_summarized, Species)
d_max_nums

#### Color Palette #### - look for color blind pallete 
pal <- c("Minke" = "#F28E2B",  "Humpback" = "#59A14F",  "Blue" = "#4E79A7", "Sei" = "#E15759", "Fin" = "#499894", "Bryde's" = "#808080", 'Normal' = "Black", 'Lunge-Associated' = "Black")
pal2 <- c("Human" = "#59A14F", "Fish" = "#E15759", "Pinniped" = "#79706E", "Sirenian" = "#B6992D", "Odontocete" = "#B07AA1", "Mysticete" = "#4E79A7", "Rodent" = "#9467BD")

#### Graphs Start Here ####

#### MST ~ U (Routine) ####
fig3 <- ggplot(d_routine_swimming_summarized, aes(mean_speed, mean_TPM)) +
  geom_smooth(method = "lm", color = "black", size = 3) +
  geom_point(aes(color = Species), size = 10) +
  scale_color_manual(values = pal) +
  expand_limits(y = 0) +
  labs(x = bquote('Swim Speed'~(m~s^-1)),
       y = bquote('Mass-Specific Thrust Power'~(Watts~kg^-1))) +
  theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 40),
        axis.title = element_text(size = 48),
        legend.position = "right",
        panel.grid.minor = element_blank())
ggsave("Figures/fig3.pdf", height = 480, width = 480, units = "mm", dpi = 300)
fig3

# Generalized linear mixed model
GLMMfig3_mean <- lmer(log(mean_TPM) ~ mean_speed + (1|Species), 
                       data = d_routine_swimming_summarized)
summary(GLMMfig3_mean)
r.squaredGLMM(GLMMfig3_mean)

#### MST ~ L (Routine vs. Lunge-Associated) ####

fig4 <- ggplot(d_routine_swimming_summarized, aes(Length, mean_TPM)) +
  geom_point(aes(color = Species), size = 10) +
  geom_point(data = d_max_swimming_summarized, aes(Length, mean_TPM, color = Species), size = 10, shape = 17) +
  geom_smooth(method = "lm", size = 3) +
  geom_smooth(data = d_max_swimming_summarized, aes(Length, mean_TPM), method = "lm", size = 3, linetype = "longdash") +
  scale_color_manual(values = pal) +
  ylim(0,2) +
  labs(x = bquote('Total Length (m)'),
       y = bquote('Mass-Specific Thrust Power'~(Watts~kg^-1))) +
  theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 40),
        axis.title = element_text(size = 48),
        legend.position = "none",
        panel.grid.minor = element_blank())
ggsave("Figures/fig4.pdf", height = 480, width = 480, units = "mm", dpi = 300)
fig4

# Generalized linear mixed models
GLMM4max_mean <- lmer(log(mean_TPM) ~ Length + (1|Species), 
                      data = d_max_swimming_summarized)
summary(GLMM4max_mean)
r.squaredGLMM(GLMM4max_mean)

GLMM4normal_mean <- lmer(log(mean_TPM) ~ Length + (1|Species), 
                         data = d_routine_swimming_summarized)
summary(GLMM4normal_mean)
r.squaredGLMM(GLMM4normal_mean)


#### MST ~ Fluke Area ####
fig5 <- ggplot(d_routine_swimming_summarized, aes(FA_L, mean_TPM)) +
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
ggsave("Figures/fig5.pdf", height = 480, width = 480, units = "mm", dpi = 300)
fig5


# Generalized linear mixed model
GLMMfig5_mean <- lmer(log(mean_TPM) ~ FA_L + (1|Species), 
                      data = d_routine_swimming_summarized)
summary(GLMMfig5_mean)
r.squaredGLMM(GLMMfig5_mean)

LMfig5test <- lm(log(mean_TPM) ~ FA_L,
                 data = d_routine_swimming_summarized)
summary(LMfig5test)

AIC(GLMMfig5_mean, LMfig5test)

GLMMfig5_test <- lmer(mean_TPM ~ FA_L + (1|Species),
                       family = "quasipoisson",
                      data = d_routine_swimming_summarized)
summary(GLMMfig5_test)
r.squaredGLMM(GLMMfig5_mean)


#### Hoerner Model Comparison ####
fig6 <- ggplot(d_routine_swimming_summarized) +
  geom_point(aes(x = mean_Re, y = mean_drag, color = Species), size = 10) +
  geom_smooth(method = "lm", aes(x = ReynoldsModel, y = DragCoeffModel), color = "black", linetype = 2, size = 3) +
  geom_smooth(method = "lm", aes(x = mean_Re, y = mean_drag), color = "black", size = 3) +
  scale_color_manual(values = pal) +
  ylim(0, 0.05) +
  labs(x = bquote('Reynolds Number'),
       y = bquote('Drag Coefficient')) +
  theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 40),
        axis.title = element_text(size = 48),
        legend.position = "none",
        panel.grid.minor = element_blank())
ggsave("Figures/fig6.pdf", height = 480, width = 480, units = "mm", dpi = 300)
fig6

# Generalized linear mixed model
GLMMfig6Emp_mean <- lmer(log(mean_drag) ~ mean_Re + (1|Species), 
                      data = d_routine_swimming_summarized)
summary(GLMMfig6Emp_mean)
r.squaredGLMM(GLMMfig6Emp_mean)

GLMMfig6Mod_mean <- lmer(log(DragCoeffModel) ~ ReynoldsModel + (1|Species), 
                         data = d_routine_swimming_summarized)
summary(GLMMfig6Mod_mean)
r.squaredGLMM(GLMMfig6Mod_mean)

#### Prop Efficiency ~ U, L ####
# Propulsive Efficiency Vs Speed
fig7U <- ggplot(d_routine_swimming_summarized, aes(mean_speed, mean_E)) +
  geom_point(aes(color = Species), size = 10) +
  geom_smooth(method = "lm", size = 3, color = "black") +
  scale_color_manual(values = pal) +
  labs(x = bquote('Speed'~(m~s^-1)),
       y = "Propulsive Efficiency") +
  expand_limits(y = c(0.75, 1)) +
  theme_bw(base_size = 20, base_family = "Times") +
  theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 40),
        axis.title = element_text(size = 48),
        legend.position = "none",
        panel.grid.minor = element_blank())
fig7U

# Generalized linear mixed model
GLMMfig7U_mean <- lmer(mean_E ~ mean_speed + (1|Species), 
                        data = d_routine_swimming_summarized)
summary(GLMMfig7U_mean)
r.squaredGLMM(GLMMfig7U_mean)

# Propulsive Efficiency Vs Length
fig7L <- ggplot(d_routine_swimming_summarized, aes(Length, mean_E)) +
  geom_point(aes(color = Species), size = 10) +
  geom_smooth(method = "lm", size = 3, color = "black") +
  scale_color_manual(values = pal) +
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
fig7L

# Generalized linear mixed model
GLMMfig7L_mean <- lmer(mean_E ~ Length + (1|Species), 
                        data = d_routine_swimming_summarized)
summary(GLMMfig7L_mean)
r.squaredGLMM(GLMMfig7L_mean)

# Combine two figues into one
fig7 <- plot_grid(fig7U, fig7L,
                  nrow = 1,
                  align = "h",
                  labels = NULL)
ggsave("Figures/fig7.pdf", height = 480, width = 960, units = "mm", dpi = 300)
fig7

#### Prop Eff ~ L w/ Other Species ####
prop_effs <- read_csv("Propulsive Eff All Species.csv")
fig8 <- ggplot(prop_effs, aes(`Total Length (m)`, `Prop Eff (Max)`)) +
  geom_point(aes(color = Group, shape = `Type of Swimming`), size = 10) +
  scale_color_manual(values = pal2) +
  expand_limits(y = c(0, 1)) +
  labs(x = "Total Length (m)",
       y = "Propulsive Efficiency") +
  theme_classic() +
  theme(axis.text = element_text(size = 40),
        axis.title = element_text(size = 48),
        legend.position = "none",
        panel.grid.minor = element_blank())
ggsave("Figures/fig8.pdf", height = 480, width = 480, units = "mm", dpi = 300)
fig8

#### Prop Eff Cetaceans Only ####
cet_prop_effs <- filter(prop_effs, Group %in% c("Odontocete", "Mysticete"))
fig9 <- ggplot(cet_prop_effs, aes(log10(`Total Length (m)`), `Prop Eff (Max)`)) +
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
ggsave("Figures/fig9.pdf", height = 480, width = 480, units = "mm", dpi = 300)
fig9

odont_prop_effs <- filter(prop_effs, Group %in% "Odontocete")
mysti_prop_effs <- filter(prop_effs, Group %in% "Mysticete")

# Linear Regression
GLMMfig9O_mean <- lm(`Prop Eff (Max)` ~ `Total Length (m)`, 
                       data = odont_prop_effs)
summary(GLMMfig9O_mean)
r.squaredGLMM(GLMMfig9O_mean)

GLMMfig9M_mean <- lm(`Prop Eff (Max)` ~ `Total Length (m)`, 
                       data = mysti_prop_effs)
summary(GLMMfig9M_mean)
r.squaredGLMM(GLMMfig9M_mean)

#### Extra Figures ####

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

#### Drag ~ U (Routine) ####
figDragU <- ggplot(d_routine_swimming_summarized, aes(mean_speed, mean_drag)) +
  geom_smooth(method = "lm", color = "black", size = 3) +
  geom_point(aes(color = Species), size = 10) +
  scale_color_manual(values = pal) +
  expand_limits(y = 0) +
  labs(x = bquote('Swim Speed'~(m~s^-1)),
       y = bquote('Drag Coefficient')) +
  theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 40),
        axis.title = element_text(size = 48),
        legend.position = "right",
        panel.grid.minor = element_blank())
ggsave("Figures/figDragU.pdf", height = 480, width = 480, units = "mm", dpi = 300)
figDragU

