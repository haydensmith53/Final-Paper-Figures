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
morphometrics <- read_csv("Finalized Data Sheet For Hayden.csv")

morphometricsFilt <- morphometrics %>% # filter out mn180302-47 (no lunge-associated tailbeats)
  filter(DeployID == "mn180302-47")

#All Swimming Flukebeat Info
d_all_swimming <- read_csv("AllDronedFlukebeatsFinalized.csv") %>%
  left_join(Mass_SKR, by = "Species") %>% 
  left_join(select(morphometrics, DeployID, FinenessRatio), by = "DeployID") %>% 
  mutate(Mass = (TotLength^slope)*10^intercept,
         FinenessRatio = FinenessRatio,
         Freq = OsFreq,
         Thust = Thrust,
         TPM = Thrust/Mass,
         Speed_BL = AvgSpeeds/TotLength,
         FA = FlukeArea,
         FA_L = FlukeArea/TotLength,
         FA_FR = FlukeArea/FinenessRatio)

d_all_swimming_summarized <- d_all_swimming %>%  
  group_by(DeployID) %>% 
  summarise(mean_freq = mean(Freq),
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
            Speed = first(AvgSpeeds),
            Effort = first(MaxOrNormal),
            FinenessRatio = first(FinenessRatio),
            FA = first(FA),
            FA_L = first(FA_L),
            FA_FR = first(FA_FR)) %>% 
  mutate(ReynoldsModel = morphometrics$`Reynolds Number Model`,
         DragCoeffModel = morphometrics$`drag coefficient C_D Model`)

d_max_swimming <- d_all_swimming %>% 
  filter(MaxOrNormal == "Lunge-Associated")

d_max_swimming_summarized <- d_max_swimming %>%  
  group_by(DeployID) %>% 
  summarise(mean_freq = mean(Freq),
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
            Speed = first(AvgSpeeds),
            FinenessRatio = first(FinenessRatio),
            FA = first(FA),
            FA_L = first(FA_L),
            FA_FR = first(FA_FR)) %>%  
  mutate(ReynoldsModel = morphometricsFilt$`Reynolds Number Model`, 
         DragCoeffModel = morphometricsFilt$`drag coefficient C_D Model`)

d_routine_swimming <- d_all_swimming %>%
  filter(MaxOrNormal == "Routine" || "Unknown")

d_routine_swimming_summarized <- d_routine_swimming %>%  
  group_by(DeployID) %>% 
  summarise(mean_freq = mean(Freq),
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
            Speed = first(AvgSpeeds),
            FinenessRatio = first(FinenessRatio),
            FA = first(FA),
            FA_L = first(FA_L),
            FA_FR = first(FA_FR)) %>% 
  mutate(ReynoldsModel = morphometrics$`Reynolds Number Model`,
         DragCoeffModel = morphometrics$`drag coefficient C_D Model`)

Deltas <- select(d_routine_swimming_summarized, DeployID, mean_speed_routine = mean_speed, mean_TPM_routine = mean_TPM)  %>% 
  left_join(select(d_max_swimming_summarized, DeployID, mean_speed_max = mean_speed, mean_TPM_max = mean_TPM),
            by = "DeployID") %>% 
  mutate(DeltaU = mean_speed_max - mean_speed_routine,
         DeltaTPM = mean_TPM_max - mean_TPM_routine)

d_all_swimming_summarized <- d_all_swimming_summarized %>% 
  left_join(select(Deltas, DeployID, DeltaU, DeltaTPM),
            by = "DeployID")

#### Color Palette #### - look for color blind pallete 
pal <- c("Minke" = "#4E79A7",  "Humpback" = "#F28E2B",  "Blue" = "#59A14F", "Sei" = "#E15759", "Fin" = "#499894", "Bryde's" = "Black", 'Normal' = "Black", 'Lunge-Associated' = "Black")
pal2 <- c("Human" = "#59A14F", "Fish" = "#E15759", "Pinniped" = "#79706E", "Sirenian" = "#B6992D", "Odontocete" = "#B07AA1", "Mysticete" = "#4E79A7", "Rodent" = "#9467BD")

#### Graphs Start Here ####

#### MST ~ U (Routine) ####
fig3 <- ggplot(d_routine_swimming_summarized, aes(mean_speed, mean_TPM)) +
  geom_point(aes(color = Species), size = 10) +
  geom_smooth(method = "lm", color = "black", size = 3) +
  scale_color_manual(values = pal) +
  expand_limits(y = 0) +
  labs(x = bquote('Swim Speed'~(m~s^-1)),
       y = bquote('Mass-Specific Thrust Power'~(Watts~kg^-1))) +
  theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 40),
        axis.title = element_text(size = 48),
        legend.position = "none",
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
GLMMfig5_mean <- lmer(log(mean_TPM) ~ `FA_L` + (1|Species), 
                      data = d_routine_swimming_summarized)
summary(GLMMfig5_mean)
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

#### Extra Figures ####

# Delta U ~ Delta MST
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

#  Fineness Ratio ~ TPM
figFR <- ggplot(d_routine_swimming_summarized, aes(FinenessRatio, mean_TPM)) +
  geom_point(aes(color = Species), size = 10) +
  geom_smooth(method = "lm", size = 3, color = "Black") +
  scale_color_manual(values = pal) +
  labs(x = bquote('Fineness Ratio'),
       y = bquote('Log(Mean Mass-Specific Thrust Power)'~(Watts~kg^-1))) +
  theme_minimal(base_size = 8) +
  theme_bw(base_size = 20, base_family = "Times") +
  theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 40),
        axis.title = element_text(size = 48),
        legend.position = "right",
        panel.grid.minor = element_blank())
ggsave("Figures/figFR.pdf", height = 480, width = 480, units = "mm", dpi = 300)
figFR

#  Fluke Area / Fineness Ratio ~ TPM
figFAFR <- ggplot(d_routine_swimming_summarized, aes(FA_FR, mean_TPM)) +
  geom_point(aes(color = Species), size = 10) +
  geom_smooth(method = "lm", size = 3, color = "Black") +
  scale_color_manual(values = pal) +
  labs(x = bquote('Fluke Area / Fineness Ratio'),
       y = bquote('Log(Mean Mass-Specific Thrust Power)'~(Watts~kg^-1))) +
  theme_minimal(base_size = 8) +
  theme_bw(base_size = 20, base_family = "Times") +
  theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 40),
        axis.title = element_text(size = 48),
        legend.position = "none",
        panel.grid.minor = element_blank())
ggsave("Figures/figFAFR.pdf", height = 480, width = 480, units = "mm", dpi = 300)
figFAFR

# MST ~ Fluke Area (No Length)
figFA <- ggplot(d_routine_swimming_summarized, aes(FA, mean_TPM)) +
  geom_point(aes(color = Species), size = 10) +
  geom_smooth(method = "lm", size = 3, color = "Black") +
  scale_color_manual(values = pal) +
  labs(x = bquote('Fluke Area'),
       y = bquote('Log(Mean Mass-Specific Thrust Power)'~(Watts~kg^-1))) +
  theme_minimal(base_size = 8) +
  theme_bw(base_size = 20, base_family = "Times") +
  theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 40),
        axis.title = element_text(size = 48),
        legend.position = "none",
        panel.grid.minor = element_blank())
ggsave("Figures/figFA.pdf", height = 480, width = 480, units = "mm", dpi = 300)
figFA


# Custom functions ----
# Standard error function
SE = function(x){sd(x)/sqrt(sum(!is.na(x)))}

#PROBLEM!!: when I add the data points for each species data set, I am ~7,000 entries short?
# to find SE for individual species, make new data frame for each with a filter 
d_all_minke <- d_combine_swimming %>%
  filter(`Common name` == "Minke")
SE(d_all_minke$'TPM')
SE(d_all_minke$`Drag Coefficient`)
SE(d_all_minke$`Reynolds Number`)
SE(d_all_minke$Efficiency)
SE(d_all_minke$`Fluke Area (m)`)
SE(d_all_minke$`Chord Length (m)`)
SE(d_all_minke$`Total Length (m)`)

d_norm_minke <- d_reg_swimming %>% 
  filter(`Common name` == "Minke")
mean(d_norm_minke$Speed)
SE(d_norm_minke$Speed)
mean(d_norm_minke$Frequency)
SE(d_norm_minke$Frequency)
d_max_minke <- d_max_swimming %>%
  filter(`Common name` == "Minke")
mean(d_max_minke$Speed)
SE(d_max_minke$Speed)
mean(d_max_minke$Frequency)
SE(d_max_minke$Frequency)

d_all_humpback <- d_combine_swimming %>%
  filter(`Common name` == "Humpback")
SE(d_all_humpback$'TPM')
SE(d_all_humpback$`Drag Coefficient`)
SE(d_all_humpback$`Reynolds Number`)
SE(d_all_humpback$Efficiency)
SE(d_all_humpback$`Fluke Area (m)`)
SE(d_all_humpback$`Chord Length (m)`)
SE(d_all_humpback$`Total Length (m)`)

d_norm_humpback <- d_reg_swimming %>% 
  filter(`Common name` == "Humpback")
mean(d_norm_humpback$Speed)
SE(d_norm_humpback$Speed)
mean(d_norm_humpback$Frequency)
SE(d_norm_humpback$Frequency)
d_max_humpback <- d_max_swimming %>%
  filter(`Common name` == "Humpback")
mean(d_max_humpback$Speed)
SE(d_max_humpback$Speed)
mean(d_max_humpback$Frequency)
SE(d_max_humpback$Frequency)

d_all_fin <- d_combine_swimming %>%
  filter(`Common name` == 'Fin')
SE(d_all_fin$'TPM')
SE(d_all_fin$`Drag Coefficient`)
SE(d_all_fin$`Reynolds Number`)
SE(d_all_fin$Efficiency)
SE(d_all_fin$`Fluke Area (m)`)
SE(d_all_fin$`Chord Length (m)`)
SE(d_all_fin$`Total Length (m)`)

d_norm_fin <- d_reg_swimming %>% 
  filter(`Common name` == "Fin")
mean(d_norm_fin$Speed)
SE(d_norm_fin$Speed)
mean(d_norm_fin$Frequency)
SE(d_norm_fin$Frequency)
d_max_fin <- d_max_swimming %>%
  filter(`Common name` == "Fin")
mean(d_max_fin$Speed)
SE(d_max_fin$Speed)
mean(d_max_fin$Frequency)
SE(d_max_fin$Frequency)

d_all_sei <- d_combine_swimming %>%
  filter(`Common name` == 'Sei')
SE(d_all_sei$'TPM')
SE(d_all_sei$`Drag Coefficient`)
SE(d_all_sei$`Reynolds Number`)
SE(d_all_sei$Efficiency)
SE(d_all_sei$`Fluke Area (m)`)
SE(d_all_sei$`Chord Length (m)`)
SE(d_all_sei$`Total Length (m)`)

d_norm_sei <- d_reg_swimming %>% 
  filter(`Common name` == "Sei")
mean(d_norm_sei$Speed)
SE(d_norm_sei$Speed)
mean(d_norm_sei$Frequency)
SE(d_norm_sei$Frequency)
d_max_sei <- d_max_swimming %>%
  filter(`Common name` == "Sei")
mean(d_max_sei$Speed)
SE(d_max_sei$Speed)
mean(d_max_sei$Frequency)
SE(d_max_sei$Frequency)

d_all_brydes <- d_combine_swimming %>%
  filter(`Common name` == 'Brydes')
SE(d_all_brydes$'TPM')
SE(d_all_brydes$`Drag Coefficient`)
SE(d_all_brydes$`Reynolds Number`)
SE(d_all_brydes$Efficiency)
SE(d_all_brydes$`Fluke Area (m)`)
SE(d_all_brydes$`Chord Length (m)`)
SE(d_all_brydes$`Total Length (m)`)

d_norm_brydes <- d_reg_swimming %>% 
  filter(`Common name` == "Brydes")
mean(d_norm_brydes$Speed)
SE(d_norm_brydes$Speed)
mean(d_norm_brydes$Frequency)
SE(d_norm_brydes$Frequency)
d_max_brydes <- d_max_swimming %>%
  filter(`Common name` == "Brydes")
mean(d_max_brydes$Speed)
SE(d_max_brydes$Speed)
mean(d_max_brydes$Frequency)
SE(d_max_brydes$Frequency)

d_all_blue <- d_combine_swimming %>%
  filter(`Common name` == 'Blue')
SE(d_all_blue$'TPM')
SE(d_all_blue$`Drag Coefficient`)
SE(d_all_blue$`Reynolds Number`)
SE(d_all_blue$Efficiency)
SE(d_all_blue$`Fluke Area (m)`)
SE(d_all_blue$`Chord Length (m)`)
SE(d_all_blue$`Total Length (m)`)

d_norm_blue <- d_reg_swimming %>% 
  filter(`Common name` == "Blue")
mean(d_norm_blue$Speed)
SE(d_norm_blue$Speed)
mean(d_norm_blue$Frequency)
SE(d_norm_blue$Frequency)
d_max_blue <- d_max_swimming %>%
  filter(`Common name` == "Blue")
mean(d_max_blue$Speed)
SE(d_max_blue$Speed)
mean(d_max_blue$Frequency)
SE(d_max_blue$Frequency)


