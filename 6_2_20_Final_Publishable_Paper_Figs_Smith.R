#### Load packages and data ####
library(cowplot)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(Hmisc)

install.packages("lme4")
install.packages("lmerTest")
install.packages("MuMIn")
library(lme4)
library(lmerTest)
library(MuMIn)

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
  "Balaenoptera bonaerensis",     1.8454,  2.1399,
  "Balaenoptera physalus",     2.595,  1.28388,
  "Balaenoptera musculus",     3.3346, 0.3129,
  "Megaptera novaeangliae",     2.3373,  1.8535,
  "Balaenoptera brydei",      2.7429, 1.0765,
  "Balaenoptera borealis",    2.4238, 1.4772
)

#### Flukebeat and Morphometric Data #### 
morphometrics <- read_csv("Finalized Data Sheet For Hayden.csv")

#All Swimming Flukebeat Info
d_all_swimming <- read_csv("AllDronedFlukebeatsFinalized.csv") %>% 
left_join(Mass_SKR, by = "Species") %>% 
  mutate(Mass = (TotLength^slope)*10^intercept, 
         TPM = Thrust/Mass,
         Speed_BL = AvgSpeeds/TotLength,
         FA_L = FlukeArea/TotLength)

d_all_swimming_summarized <- d_all_swimming %>%  
  group_by(DeployID) %>% 
  summarise(mean_TPM = mean(TPM), 
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
            Effort = first(MaxOrNormal))

d_max_swimming <- d_all_swimming %>% 
  filter(MaxOrNormal == "Lunge-Associated")
  left_join(Mass_SKR, by = "Species") %>% 
  mutate(Mass = (TotLength^slope)*10^intercept, 
         TPM = Thrust/Mass,
         Speed_BL = AvgSpeeds/TotLength)

d_max_swimming_summarized <- d_max_swimming %>%  
  group_by(DeployID) %>% 
  summarise(mean_TPM = mean(TPM), 
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
            Speed = first(AvgSpeeds))

# Separating normal fluekbeats from all data ----
d_routine_swimming <- d_all_swimming %>%
  filter(MaxOrNormal == "Routine")
left_join(Mass_SKR, by = "Species") %>% 
  mutate(Mass = (TotLength^slope)*10^intercept, 
         TPM = Thrust/Mass,
         Speed_BL = AvgSpeeds/TotLength,
         FA_L = FlukeArea/TotLength)

d_routine_swimming_summarized <- d_routine_swimming %>%  
  group_by(DeployID) %>% 
  summarise(mean_TPM = mean(TPM), 
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
            FA_L = first(FA_L))

d_combine_swimming_summarized <- d_all_swimming %>%  
  group_by(DeployID, MaxOrNormal) %>% 
  summarise(mean_TPM = mean(TPM), 
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
            Speed = first(AvgSpeeds)) 

#### Color Palette #### - look for color blind pallete 
pal <- c("B. bonaerensis" = "CC79A7",  "M. novaeangliae" = "009E73",  "B. musculus" = "#0072B2", "B. borealis" = "E69F00", "B. physalus" = "F0E442", "B. brydei" = "black")

#### Graphs Start Here ####

#### MST ~ U, A/L, L ####
# (kinematics and morphology) 
# (Just routine, no max effort)
routine_effort <- d_routine_swimming_summarized %>% 
  mutate(sp_abbr = abbr_binom(Species))
x_breaks <- seq(1.0, 3.0, by = 0.25)
y_breaks <- seq(0.0, 0.8, by = 0.1)
x_lbls <- ifelse(x_breaks %% 0.5 == 0, as.character(x_breaks), "")
y_lbls <- ifelse(y_breaks %% 0.2 == 0, as.character(y_breaks), "")
x_ticksz <- ifelse(x_breaks %% 0.5 == 0, 1, 0.5)
y_ticksz <- ifelse(y_breaks %% 0.2 == 0, 1, 0.5)
fig3_U <- ggplot(routine_effort, aes(mean_speed, mean_TPM)) +
  geom_point(aes(color = sp_abbr), size = 1) +
  geom_smooth(method = "lm", color = "black", size = 1) +
  scale_color_manual(values = pal) +
  scale_x_continuous(breaks = x_breaks, labels = x_lbls) +
  scale_y_continuous(breaks = y_breaks, labels = y_lbls) +
  expand_limits(y = 0) +
  labs(x = bquote('Speed'~(m~s^-1)),
       y = bquote('Log(Mean Mass-Specific Thrust)'~(N~kg^-1))) +
  theme_bw(base_size = 20, base_family = "Times") +
  theme(legend.position = "none",
        panel.grid = element_blank(), 
        axis.ticks.x = element_line(size = x_ticksz),
        axis.ticks.y = element_line(size = y_ticksz),
        axis.text = element_text(color = "black")) 
fig3_U

# Stats
# Generalized linear mixed model
GLMMfig3U_mean <- lmer(log(mean_TPM) ~ mean_speed + (1|Species), 
                       data = d_combine_swimming_summarized)
summary(GLMMfig3U_mean)
r.squaredGLMM(GLMMfig3U_mean)

x_breaks <- seq(0.1, 0.4, by = 0.05)
y_breaks <- seq(0.0, 0.8, by = 0.1)
x_lbls <- ifelse(x_breaks %% 0.1 == 0, as.character(x_breaks), "")
y_lbls <- ifelse(y_breaks %% 0.2 == 0, as.character(y_breaks), "")
x_ticksz <- ifelse(x_breaks %% 0.1 == 0, 1, 0.4)
y_ticksz <- ifelse(y_breaks %% 0.2 == 0, 1, 0.4)
fig3_AL <- ggplot(routine_effort, aes(FA_L, mean_TPM)) +
  geom_point(aes(color = sp_abbr), size = 1) +
  geom_smooth(method = "lm", color = "black", size = 1) +
  scale_color_manual(values = pal) +
  scale_x_continuous(breaks = x_breaks, labels = x_lbls) +
  scale_y_continuous(breaks = y_breaks, labels = y_lbls) +
  expand_limits(y = 0) +
  labs(x = "Fluke Area / Length (m)") +
  theme_minimal(base_size = 8) +
  theme_bw(base_size = 20, base_family = "Times") +
  theme(axis.title.y = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(),
        axis.ticks.x = element_line(size = x_ticksz),
        axis.ticks.y = element_line(size = y_ticksz),
        axis.text = element_text(color = "black"))
fig3_AL

# Stats
# Generalized linear mixed model
GLMMfig3AL_mean <- lmer(log(mean_TPM) ~ `FA/L` + (1|Species), 
                        data = d_combine_swimming_summarized)
summary(GLMMfig3AL_mean)
r.squaredGLMM(GLMMfig3AL_mean)

## TODO change text and point sizes, fit linear mixed effects model, get prediction intervals from lme
y_breaks <- seq(0.0, 1.4, by = 0.2)
x_lbls <- ifelse(x_breaks %% 5 == 0, as.character(x_breaks), "")
y_lbls <- ifelse(y_breaks %% 0.4 == 0, as.character(y_breaks), "")
x_ticksz <- ifelse(x_breaks %% 5 == 0, 1, 0.5)
y_ticksz <- ifelse(y_breaks %% 0.4 == 0, 0.5, 1)
#### MST ~ L, max vs normal ####
combined_effort <- d_combine_swimming_summarized %>%
  mutate(sp_abbr = abbr_binom(Species))
# Remove outliers from figure
fig3_L <- combined_effort %>% 
  filter(mean_TPM < 1.6) %>% 
  ggplot(aes(Length, mean_TPM, color = sp_abbr)) +
  geom_point(aes(shape = MaxOrNormal)) +
  geom_smooth(method = "lm", aes(color = MaxOrNormal, linetype = MaxOrNormal)) +
  scale_color_manual(values = c(`B. bonaerensis` = "#009E73",
                                `M. novaeangliae` = "#D55E00",
                                `B. musculus` = "#0072B2",
                                `B. borealis` = "E69F00", 
                                `B. physalus` = "F0E442", 
                                `B. brydei` = "black",
                                Routine = "black", Lunge-Associated = "black")) +
  scale_x_continuous(breaks = x_breaks, labels = x_lbls) +
  scale_y_continuous(breaks = y_breaks, labels = y_lbls) +
  scale_shape_manual(values = c(Routine = 16, Lunge-Associated = 1)) +
  scale_linetype_manual(values = c(Routine = "solid", Lunge-Associated = "longdash")) +
  labs(x = bquote('Total Length (m)'),
       y = bquote('Log(Mean Mass-Specific Thrust)'~(N~kg^-1))) +
  theme_minimal(base_size = 8) +
  theme_bw(base_size = 20, base_family = "Times") +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.ticks.x = element_line(size = x_ticksz),
        axis.ticks.y = element_line(size = y_ticksz),
        axis.text = element_text(color = "black"))
fig3_L

# Stats
# Generalized linear mixed model
GLMM4max_mean <- lmer(log(mean_TPM) ~ `Total Length (m)` + (1|Species), 
                      data = filter(d_combine_swimming_summarized, effort_type == "Max"))
summary(GLMM4max_mean)
r.squaredGLMM(GLMM4max_mean)

# Stats
# Generalized linear mixed model
GLMM4normal_mean <- lmer(log(mean_TPM) ~ `Total Length (m)` + (1|Species), 
                         data = filter(d_combine_swimming_summarized, effort_type == "Normal"))
summary(GLMM4normal_mean)
r.squaredGLMM(GLMM4normal_mean)

# Combine into one figure
fig3 <- plot_grid(fig3_U, fig3_AL, fig3_L,
                  nrow = 2,
                  align = "h",
                  rel_widths = c(1.1, 0.9),
                  labels = NULL)
fig3

#### Prop Efficiency ~ U, L ####
x_breaks <- seq(1.0, 3.0, by = .25)
y_breaks <- seq(0.775, 1.00, by = .025)
x_lbls <- ifelse(x_breaks %% 0.5 == 0, as.character(x_breaks), "")
y_lbls <- ifelse(y_breaks %% 0.05 == 0, as.character(y_breaks), "")
x_ticksz <- ifelse(x_breaks %% 0.5 == 0, 1, 0.5)
y_ticksz <- ifelse(y_breaks %% 0.05 == 0, 1, 0.5)
fig5_U <- routine_effort %>% 
  ggplot(aes(mean_speed, mean_E)) +
  geom_point(aes(color = sp_abbr)) +
  geom_smooth(method = "lm", color = "black") +
  scale_color_manual(values = pal) +
  scale_x_continuous(breaks = x_breaks, labels = x_lbls) +
  scale_y_continuous(breaks = y_breaks, labels = y_lbls) +
  labs(x = bquote('Speed'~(m~s^-1)),
       y = "Propulsive Efficiency") +
  expand_limits(y = c(0.75, 1)) +
  theme_bw(base_size = 20, base_family = "Times") +
  theme(legend.position = "none",
        panel.grid = element_blank(), 
        axis.ticks.x = element_line(size = x_ticksz),
        axis.ticks.y = element_line(size = y_ticksz),
        axis.text = element_text(color = "black")) 
fig5_U

# Generalized linear mixed model
GLMM4fig5U_mean <- lmer(mean_E ~ mean_speed + (1|Species), 
                        data = normal_effort)
summary(GLMM4fig5U_mean)
r.squaredGLMM(GLMM4fig5U_mean)

x_breaks <- seq(5, 25, by = 2.5)
y_breaks <- seq(0.775, 1.0, by = 0.025)
x_lbls <- ifelse(x_breaks %% 5 == 0, as.character(x_breaks), "")
y_lbls <- ifelse(y_breaks %% 0.05 == 0, as.character(y_breaks), "")
x_ticksz <- ifelse(x_breaks %% 5 == 0, 1, 0.5)
y_ticksz <- ifelse(y_breaks %% 0.05 == 0, 1, 0.5)
fig5_L <- routine_effort %>% 
  ggplot(aes(Length, mean_E)) +
  geom_point(aes(color = sp_abbr)) +
  geom_smooth(method = "lm", color = "black") +
  scale_color_manual(values = pal) +
  scale_x_continuous(breaks = x_breaks, labels = x_lbls) +
  scale_y_continuous(breaks = y_breaks, labels = y_lbls) +
  expand_limits(y = c(0.75, 1)) +
  labs(x = bquote('Total Length'~(m)),
       y = bquote('Propulsive Efficiency')) +
  theme_minimal(base_size = 8) +
  theme_bw(base_size = 20, base_family = "Times") +
  theme(legend.position = "none",
        panel.grid = element_blank(), 
        axis.ticks.x = element_line(size = x_ticksz),
        axis.ticks.y = element_line(size = y_ticksz),
        axis.text = element_text(color = "black")) 
fig5_L

# Stats
# Generalized linear mixed model
GLMM4fig5L_mean <- lmer(mean_E ~ `Total Length (m)` + (1|Species), 
                        data = normal_effort)
summary(GLMM4fig5L_mean)
r.squaredGLMM(GLMM4fig5L_mean)


fig5 <- plot_grid(fig5_U, fig5_L,
                  nrow = 1,
                  align = "h",
                  labels = NULL)
fig5

#### Prop Eff ~ U (Smoothed Lines + Fish 1998 Data) ####
x_breaks <- seq(2.0, 8.0, by = 1.0)
y_breaks <- seq(0.775, 1.00, by = .025)
x_lbls <- ifelse(x_breaks %% 2 == 0, as.character(x_breaks), "")
y_lbls <- ifelse(y_breaks %% 0.05 == 0, as.character(y_breaks), "")
x_ticksz <- ifelse(x_breaks %% 2.0 == 0, 1, 0.5)
y_ticksz <- ifelse(y_breaks %% 0.05 == 0, 1, 0.5)
fish_prop_eff <- read_csv("fish_prop_eff.csv") 
flukebeats_no_fins <- d_all_swimming_summarized %>%
  mutate(sp_abbr = abbr_binom(Species))
fig6 <- ggplot(flukebeats_no_fins, aes(speed, mean_E)) +
  geom_smooth(aes(color = Species), se = FALSE) +
  geom_smooth(aes(color = Species), fish_prop_eff, se = FALSE, linetype = 2) +
  scale_x_continuous(breaks = x_breaks, labels = x_lbls) +
  scale_y_continuous(breaks = y_breaks, labels = y_lbls) +
  labs(x = bquote('Speed'~(m~s^-1)),
       y = "Propulsive Efficiency") +
  theme_bw(base_size = 20, base_family = "Times") +
  theme(legend.position = "none",
        panel.grid = element_blank(), 
        axis.ticks.x = element_line(size = x_ticksz),
        axis.ticks.y = element_line(size = y_ticksz),
        axis.text = element_text(color = "black")) 
#theme_minimal() +
#theme_classic()
#theme(legend.position = "none",
#     panel.grid.minor = element_blank())
fig6



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


