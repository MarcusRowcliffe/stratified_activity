library(cowplot)
source("simulation_functions.r")

# Define parameters ####
# Habitat 1 is active; Habitat 2 is resting
N <- 50 # population size
nstr <- 2 # number of strata
# number of deployments
deps <- 100
# stratum areas
area <- list(skew = c(9, 1),
             even = c(5, 5))
# observation sample sizes
ssize <- list(low=100, med=200, high=500) 
# Population distribution parameters
pprm <- list(strong=6, weak=1)
# Activity pattern parameters
prms_S1 <- make_act_prms(times = c(6, 10, 15, 18) * pi/12,
                         slopes = c(4, 1, 1, 4),
                         heights = c(dawn=6, day=2, dusk=4, night=0),
                         nper = c(4, 8, 4, 8),
                         mx = 1)
prms_S2 <- make_act_prms(times = c(6, 10, 15, 18) * pi/12,
                         slopes = c(4, 1, 1, 4),
                         heights = c(dawn=6, day=1, dusk=2, night=0),
                         nper = c(4, 8, 4, 8),
                         mx = 1)
prms_S3 <- make_act_prms(times = c(6, 10, 15, 18) * pi/12,
                         slopes = c(4, 1, 1, 4),
                         heights = c(dawn=6, day=1, dusk=2, night=0),
                         nper = c(4, 8, 4, 8),
                         mx = 0.5)
prms_S4 <- make_act_prms(times = c(8, 10, 14, 16) * pi/12,
                         slopes = c(3, 1, 1, 3),
                         heights = c(dawn=1, day=1, dusk=2, night=0),
                         nper = c(4, 8, 4, 8),
                         mx = 0.5)
prms_W1 <- make_act_prms(times = c(6, 10, 15, 18) * pi/12,
                         slopes = c(1, 1, 1, 1),
                         heights = c(dawn=6, day=2, dusk=4, night=2),
                         nper = c(4, 8, 4, 8),
                         mx = 1)
prms_W2 <- make_act_prms(times = c(6, 10, 15, 18) * pi/12,
                         slopes = c(1, 1, 1, 1),
                         heights = c(dawn=6, day=1, dusk=2, night=1),
                         nper = c(4, 8, 4, 8),
                         mx = 1)
prms_W3 <- make_act_prms(times = c(6, 10, 15, 18) * pi/12,
                         slopes = c(1, 1, 1, 1),
                         heights = c(dawn=6, day=1, dusk=2, night=1),
                         nper = c(4, 8, 4, 8),
                         mx = 0.5)
prms_W4 <- make_act_prms(times = c(6, 10, 15, 18) * pi/12,
                         slopes = c(1, 1, 1, 1),
                         heights = c(dawn=2, day=1, dusk=4, night=1),
                         nper = c(4, 8, 4, 8),
                         mx = 0.5)
aprm <- list(strong=list(prms_S1, prms_S2, prms_S3, prms_S4),
             weak=list(prms_W1, prms_W2, prms_W3, prms_W4))

# Plot scenario patterns ####
scenario_plots <- lapply(1:8, function(s) plot_scenario_patterns(s))
lgnd <- plot_grid(get_legend(scenario_plots[[1]]), scale=0.1)
scenario_plots <- scenario_plots %>%
  map(\(p)  p + theme(legend.position = "none"))
scenario_plots <- plot_grid(plotlist=scenario_plots, nrow=4, byrow=F)
void <- ggplot() + theme_void()
xaxis <- ggplot(mapping=aes(x=0, y=0, label="Time")) + geom_text() + theme_void()
yaxis <- ggplot(mapping=aes(x=0, y=0, label="Activity level / Population distribution")) + 
  geom_text(angle=90) + theme_void()
s <- 0.5
sidebar <- ggplot() + 
  lims(x=c(0,2)) +
  annotate("segment", x=1, y=4, xend=1, yend=0, arrow=arrow(type="closed")) +
  geom_text(aes(x=1+s, y=2, label="Decreasing alignment between activity patterns"), angle=270) +
  geom_text(aes(x=1-s, y=0.5+(0:3), label=paste0("A", 1:4))) +
  scale_y_continuous(expand=c(0,0)) +
  theme_void()
topbar <- ggplot() + 
  lims(y=c(0,2)) +
  annotate("segment", x=0, y=1, xend=2, yend=1, arrow=arrow(type="closed")) +
  geom_text(aes(x=1, y=1+s, label="Decreasing strength of shift in patterns")) +
  geom_text(aes(x=0.5+(0:1), y=1-s, label=c("Strong", "Weak"))) +
  scale_x_continuous(expand=c(0,0)) +
  theme_void()
top <- plot_grid(void, topbar, void, nrow=1, rel_widths = c(1,20,10))
mid <- plot_grid(yaxis, scenario_plots, sidebar, lgnd, nrow=1, rel_widths=c(1,20,4,6))
bot <- plot_grid(void, xaxis, void, nrow=1, rel_widths = c(1,20,10))
scenario_patterns <- plot_grid(top, mid, bot, ncol=1, rel_heights=c(3, 20, 1))
ggsave("scenario_patterns.png", scenario_patterns)

# Run simulations ####
scenarios <- expand.grid(alignment=1:4,
                         pattern=c("strong","weak"),
                         config=c("repEven", "repSkew", "strSkew"),
                         effort=c("low", "med", "high"),
                         stringsAsFactors = FALSE) %>%
  mutate(depdat = case_when(config=="strSkew" ~ "strat", .default = "rep"),
         strdat = case_when(config=="repEven" ~ "even", .default = "skew"))

res <- lapply(1:nrow(scenarios), run_scenario_s, reps=4)
errors <- calc_errors(res) %>%
  mutate(effort = fct_relevel(effort, "low", "med", "high"))
sum(errors$n1==0)
sum(errors$n2==0)
View(errors)

# sample sizes
errors %>%
  select(n1, n2, alignment, pattern, config, effort) %>%
  pivot_longer(starts_with("n"), names_to="stratum") %>%
  rename(Observations = value) %>%
  mutate(stratum = ifelse(stratum=="n1", 1, 2),
         Effort_Stratum = paste(effort, stratum, sep="_")) %>%
  mutate(Effort_Stratum = fct_relevel(Effort_Stratum, 
                                      "low_1", "low_2", 
                                      "med_1", "med_2", 
                                      "high_1", "high_2")) %>%
  ggplot(aes(x=config, y=Observations, fill=Effort_Stratum)) + 
  geom_boxplot(coef=100) +
  theme_minimal_hgrid() +
  facet_grid(vars(alignment), vars(pattern))

# Population distribution errors
ggplot(errors, aes(x=config, y=popDist, fill=effort)) + 
  geom_boxplot(coef=50) +
  facet_grid(vars(alignment), vars(pattern))

# Activity level errors
errors %>%
  select(act_strat, act_naive, alignment, pattern, config, effort) %>%
  pivot_longer(starts_with("act_"), names_to="method") %>%
  rename(Error = value) %>%
  mutate(method = ifelse(method=="act_strat", "stratified", "unstratified"),
         Effort_Method = paste(effort, method, sep="_")) %>%
  mutate(Effort_Method = fct_relevel(Effort_Method, 
                                     "low_unstratified", "med_unstratified", "high_unstratified",
                                     "low_stratified", "med_stratified", "high_stratified")) %>%
  ggplot(aes(x=config, y=Error, fill=Effort_Method)) + 
  geom_boxplot(coef=100) +
  geom_hline(yintercept=0) +
  theme_minimal_hgrid() +
  facet_grid(vars(alignment), vars(pattern))

# Stratum-specific activity level errors
errors %>%
  select(act_str1, act_str2, alignment, pattern, config, effort) %>%
  pivot_longer(starts_with("act_"), names_to="stratum") %>%
  rename(Error = value) %>%
  mutate(stratum = ifelse(stratum=="act_str1", 1, 2),
         Effort_Stratum = paste(effort, stratum, sep="_")) %>%
  mutate(Effort_Stratum = fct_relevel(Effort_Stratum, 
                                      "low_1", "med_1", "high_1",
                                      "low_2", "med_2", "high_2")) %>%
  ggplot(aes(x=config, y=Error, fill=Effort_Stratum)) + 
  geom_boxplot(coef=100) +
  geom_hline(yintercept=0) +
  theme_minimal_hgrid() +
  facet_grid(vars(alignment), vars(pattern))

# Activity pattern errors
errors %>%
  select(actPat_strat, actPat_naive, alignment, pattern, config, effort) %>%
  pivot_longer(starts_with("act"), names_to="method") %>%
  rename(Error = value) %>%
  mutate(method = ifelse(method=="actPat_strat", "stratified", "unstratified"),
         Effort_Method = paste(effort, method, sep="_")) %>%
  mutate(Effort_Method = fct_relevel(Effort_Method, 
                                     "low_unstratified", "med_unstratified", "high_unstratified",
                                     "low_stratified", "med_stratified", "high_stratified")) %>%
  ggplot(aes(x=config, y=Error, fill=Effort_Method)) + 
  geom_boxplot(coef=100) +
  geom_hline(yintercept=0) +
  theme_minimal_hgrid() +
  facet_grid(vars(alignment), vars(pattern))

# Population distributions
popdistPlots <- lapply(res, plot_pop_distribution, labels=FALSE)
iStrong <- which(scenarios$pattern=="strong")
iWeak <- which(scenarios$pattern=="weak")
pdStrong_grid <- plot_grid(plotlist=popdistPlots[iStrong], nrow=4, byrow=FALSE)
pdWeak_grid <- plot_grid(plotlist=popdistPlots[iWeak], nrow=4, byrow=FALSE)
effort_lab <- plot_grid(plotlist=list(
  ggplot(mapping=aes(x=0, y=0, label="Low effort")) + geom_text() + theme_void(),
  ggplot(mapping=aes(x=0, y=0, label="Medium effort")) + geom_text() + theme_void(),
  ggplot(mapping=aes(x=0, y=0, label="High effort")) + geom_text() + theme_void()),
  nrow=1)
config_lab <- list(
  ggplot(mapping=aes(x=0, y=0, label="Rep/Even")) + geom_text() + theme_void(),
  ggplot(mapping=aes(x=0, y=0, label="Rep/Skew")) + geom_text() + theme_void(),
  ggplot(mapping=aes(x=0, y=0, label="Strat/Skew")) + geom_text() + theme_void())
config_lab <- plot_grid(plotlist = c(config_lab,config_lab, config_lab), 
                        nrow=1)
x_lab <- ggplot(mapping=aes(x=0, y=0, label="Time")) + geom_text() + theme_void()
y_lab_pd <- ggplot(mapping=aes(x=0, y=0, label="Population proportion", angle=90)) + 
  geom_text() + theme_void()
pdStrong <- plot_grid(effort_lab, config_lab, pdStrong_grid, ncol=1, rel_heights = c(1,1,20))
pdWeak <- plot_grid(effort_lab, config_lab, pdWeak_grid, ncol=1, rel_heights = c(1,1,20))
plot_grid(y_lab_pd, pdStrong, NULL, x_lab, nrow=2, rel_heights=c(30,1), rel_widths=c(1,30))
plot_grid(y_lab_pd, pdWeak, NULL, x_lab, nrow=2, rel_heights=c(30,1), rel_widths=c(1,30))


# Activity patterns
actpatPlots <- lapply(res, plot_activity_patterns, labels=FALSE)
apStrong_grid <- plot_grid(plotlist=actpatPlots[iStrong], nrow=4, byrow=FALSE)
apWeak_grid <- plot_grid(plotlist=actpatPlots[iWeak], nrow=4, byrow=FALSE)
y_lab_ap <- ggplot(mapping=aes(x=0, y=0, label="Activity level", angle=90)) + 
  geom_text() + theme_void()
apStrong <- plot_grid(effort_lab, config_lab, apStrong_grid, ncol=1, rel_heights = c(1,1,20))
apWeak <- plot_grid(effort_lab, config_lab, apWeak_grid, ncol=1, rel_heights = c(1,1,20))
plot_grid(y_lab_ap, apStrong, NULL, x_lab, nrow=2, rel_heights=c(30,1), rel_widths=c(1,30))
plot_grid(y_lab_ap, apWeak, NULL, x_lab, nrow=2, rel_heights=c(30,1), rel_widths=c(1,30))

