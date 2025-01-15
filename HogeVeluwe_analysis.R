# ANALYSIS: ACTIVITY/HABITAT USE ####
# Read deployments and observations 
# See HogeVeluwe_DataProcessing.r for derivation of these files
dir <- "./HogeVeluweData"
deployments <- read.csv(file.path(dir, "./HogeVeluweData/deployments.csv")) %>%
  mutate(start = lubridate::ymd_hms(start, truncated = 3),
         end = lubridate::ymd_hms(end, truncated = 3))

observations <- read.csv(file.path(dir, "./HogeVeluweData/observations.csv")) %>%
  mutate(timestamp = lubridate::ymd_hms(timestamp, truncated = 3))

## Read stratum data ####
# Raw stratum data
full_strata <- read.csv(file.path(dir, "Habitat Availability Hoge Veluwe 2012-2023.csv"))
# Stratum data for density analysis: 
# - filter out water
# - merge habitats 
# - pivot to row per stratum
# - add stratumID
den_strata <- full_strata %>%
  filter(habitat != "water") %>%
  mutate(areaOB = ifelse(year<2018, 0, areaOB),
         habitat = case_when(
           habitat=="driftsand" ~ "sand",
           habitat %in% c("forest_culture", "pine_woodland") ~ "forest",
           habitat %in% c("dry_heath", "wet_heath", "other") ~ "heath",
           habitat=="meadow" ~ "meadow",)) %>%
  pivot_longer(starts_with("area"),
               names_to = "type",
               values_to = "area") %>%
  group_by(species, year, habitat, type) %>%
  summarise(area = sum(area)) %>%
  mutate(type = sub("area", "", type),
         stratumID = paste(habitat, type))
# Stratum data for activity analysis:
# - merge OB with public
# - average areas over years
act_strata <- den_strata %>%
  mutate(type = ifelse(type=="OB", "public", type),
         stratumID = paste(habitat, type),
         species = ifelse(grepl("deer", species),
                          paste0(species, "_badger_red_fox_brown_hare"),
                          species)) %>%
  group_by(species, year, habitat, type, stratumID) %>%
  summarise(area = sum(area)) %>%
  group_by(species, habitat, type, stratumID) %>%
  summarise(area = mean(area))

## Read REM parameters ####
spd_dat <- read.csv(file.path(dir, "speed.csv")) %>%
  group_by(species, time, habitat2) %>%
  summarise(est = mean(speed),
            se = se_from_ses(speed, se)) %>%
  rename(habitat = habitat2)
#  slice(rep(1:n(), each=2)) %>%
#  mutate(timeID = ifelse(time=="day", 1, 2),
#         type = rep(c("rest", "public"), n()/2),
#         stratumID = paste(habitat, type))
#View(spd_dat)

rad_dat <- read.csv(file.path(dir, "radius.csv")) %>%
  group_by(species, habitat2) %>%
  summarise(est = mean(radius),
            se = se_from_ses(radius, se)) %>%
  rename(habitat = habitat2)
#  slice(rep(1:n(), each=2)) %>%
#  mutate(type = rep(c("rest", "public"), n()/2),
#         stratumID = paste(habitat, type))
#View(rad_dat)

ang_dat <- read.csv(file.path(dir, "angle.csv")) %>%
  group_by(time, habitat2) %>%
  summarise(est = mean(angle),
            se = se_from_ses(angle, se)) %>%
  rename(habitat = habitat2)
#  slice(rep(1:n(), each=2)) %>%
#  mutate(timeID = ifelse(time=="day", 1, 2),
#         type = rep(c("rest", "public"), n()/2),
#         stratumID = paste(habitat, type))
#View(ang_dat)


## Stratified activity models ####
lat <- 52.084025
lon <- 5.837495
st <- activity::get_suntimes("2020-04-15", lat, lon, 0)
suntimes <-  c(st$sunrise, st$sunset) * pi / 12
# Function generates a double sigmoid response to time, for 
# a smooth transition between binary day/night parameters
dblsig <- function(time, e1, e2, slp=15){
  e2 + (e1 - e2) *
    (1/(1+exp(-slp*(time-suntimes[1]))) - 
       1/(1+exp(-slp*(time-suntimes[2]))))
}
# Deployment table for activity analyses:
# - merge OB with public
# - add effort field
dep <- deployments %>%
  mutate(stratumID = sub("OB", "public", stratumID),
         effort = as.numeric(difftime(end, start, units="days")))
# function fits stratified activity model for a given species
fitfunc <- function(sp, reps){
  if(sp %in% species[5:7]){
    spd <- ang <- rad <- NULL
  } else{
    spdFunc <- function(time, stratum){
      habs <- strsplit(stratum, " ") %>%
        lapply(function(x) x[1]) %>%
        unlist()
      est <- spd_dat %>%
        filter(species==sp) %>%
        select(-se) %>%
        pivot_wider(names_from = time, values_from = est) %>%
        slice(match(habs, habitat))
      dblsig(time, est$day, est$night)
    }
    angFunc <- function(time, stratum){
      habs <- strsplit(stratum, " ") %>%
        lapply(function(x) x[1]) %>%
        unlist()
      est <- ang_dat %>%
        select(-se) %>%
        pivot_wider(names_from = time, values_from = est) %>%
        slice(match(habs, habitat))
      dblsig(time, est$day, est$night)
    }
    radFunc <- function(stratum){
      habs <- strsplit(stratum, " ") %>%
        lapply(function(x) x[1]) %>%
        unlist()
      rad_dat %>%
        filter(species==sp) %>%
        select(-se) %>%
        slice(match(habs, habitat)) %>%
        .$est
    }
    spd <- pmat_from_func(spdFunc, str)
    ang <- pmat_from_func(angFunc, str)
    rad <- pmat_from_func(radFunc, str)
  }
  str <- act_strata %>%
    filter(grepl(sp, species,))
  obs <- subset(observations, commonName==sp)
  fitact_strat(obs, dep, str, spd, rad, ang, reps=reps)
}
actmods <- lapply(species, fitfunc, 5)
names(actmods) <- species

# Activity / Population distribution / Combo plots
act_plot <- function(sp, lim=0.15){
  dat <- data.frame(x = tm,
                    y = actmods[[sp]]$est$pdf * pi/12,
                    l = actmods[[sp]]$ci$pdf[,1] * pi/12,
                    u = actmods[[sp]]$ci$pdf[,2] * pi/12)
  ggplot(dat, aes(x=x, y=y)) +
    geom_line() + 
    geom_ribbon(aes(ymin=l, ymax=u), alpha=0.2) +
    scale_x_continuous(limits=c(0, 24), breaks=seq(0,24,6)) +
    scale_y_continuous(limits=c(0, lim), breaks=c(0,0.1)) +
    labs(x = NULL, y = NULL) + 
    annotate("text", x=0, y=lim, hjust=0, vjust=1, size=3,
             label = sub("_", " ", sp)) +
    theme_classic()
}
popdist_plot <- function(sp, mods, cols){
  popdist <- mods[[sp]]$est$popdist
  pdat <- data.frame(p = as.vector(popdist),
                     stratum = rep(colnames(popdist), each=nrow(popdist)),
                     Time = rep(tm, ncol(popdist))) %>%
    filter(!grepl("sandasd", stratum))
  ggplot(pdat, aes(x=Time, y=p, fill=stratum)) + 
    geom_area() +
    scale_x_continuous(breaks=seq(0,24,6)) +
    scale_y_continuous(breaks=c(0,1)) +
    scale_fill_manual(values=cols) +
    labs(x = NULL, y = NULL) + 
    theme_classic()
}
popact_plot <- function(sp, mods, cols){
  popdist <- mods[[sp]]$est$popdist
  act <- mods[[sp]]$est$pdf
  f <- function(x) x/sum(x)
  pa <- popdist * act
  
  pdat <- data.frame(pa = as.vector(pa),
                     habitat = rep(colnames(popdist), each=nrow(popdist)),
                     Time = rep(tm, ncol(popdist)))
  plt <- ggplot(pdat, aes(x=Time, y=pa, fill=habitat)) + 
    geom_area() +
    scale_x_continuous(breaks=seq(0,24,6)) +
    scale_y_continuous(breaks=seq(0,0.3,0.1)) +
    scale_fill_manual(values=cols) +
    labs(x = NULL, y = NULL) +
    theme_classic()
  #    theme(panel.background = element_blank())
  plt + theme(legend.position = "none")
}

tm <- seq(0, 24, len=513)
pal <- scales::hue_pal()(8)

act_plots <- lapply(species, act_plot)
act_plots <- plot_grid(plotlist = act_plots, ncol=1)

popdist_plots <- lapply(species, popdist_plot, actmods, pal)
legend <- get_legend(popdist_plots[[1]])
popdist_plots <- popdist_plots %>%
  lapply(function(plt) plt + theme(legend.position = "none"))
popdist_plots <- plot_grid(plotlist = popdist_plots, ncol=1)

#popact_plots <- lapply(species, popact_plot, actmods, pal)
#popact_plots <- plot_grid(plotlist = popact_plots, ncol=1)

xlab <- ggplot(mapping=aes(x=0, y=0, label="Time (h)")) +
  geom_text() +
  theme_void()
ylab1 <- ggplot(mapping=aes(x=0, y=0, label="Activity PDF")) +
  geom_text(angle=90) +
  theme_void()
ylab2 <- ggplot(mapping=aes(x=0, y=0, label="Population proportion")) +
  geom_text(angle=90) +
  theme_void()
#ylab3 <- ggplot(mapping=aes(x=0, y=0, label="Population activity")) +
#  geom_text(angle=90) +
#  theme_void()
lgnd <- plot_grid(legend)

str <- act_strata %>%
  filter(grepl("deer", species)) %>%
  mutate(x=1)
y <- c(0, cumsum(str$area/sum(str$area)))
y <- 1 - (head(y, -1) + tail(y, -1)) / 2 + c(0,0,0,0.01,0.01,-0.01,0,0)
lgnd <- ggplot(str, aes(fill=stratumID, y=area, x=x)) + 
  geom_bar(position="fill", stat="identity", width=0.2, show.legend=FALSE) +
  annotate("text", x=1.12, y=y, label=str$stratumID, hjust=0) +
  xlim(c(0.9,3)) +
  theme_void() 
  
plots <- plot_grid(ylab1, act_plots,
                   ylab2, popdist_plots,
#                   ylab3, popact_plots, 
                   lgnd,
#                   nrow=1, rel_widths = c(rep(c(1,7),3),4))
                    nrow=1, rel_widths = c(rep(c(1,7),2),4))
print(plot_grid(plots, xlab, ncol=1, rel_heights=c(20,1)))


## Habitat selection ####
# Create plot for species sp and habitat hb
sel_plot <- function(sp, hb){
  
  electivity <- function(use, avail){
    t((use-avail) / (use+avail)) %>%
      as.data.frame() %>%
      mutate(time=tm) %>%
      pivot_longer(strata)
  }
  
  pd <- t(actmods[[sp]]$est$popdist)
  lc <- t(actmods[[sp]]$ci$popdist$lcl)
  uc <- t(actmods[[sp]]$ci$popdist$ucl)
  strata <- rownames(pd)
  str <- act_strata %>%
    filter(grepl(sp, species)) %>%
    group_by(species) %>%
    mutate(p.area = area / sum(area))
  pa <- with(str, p.area[match(strata, stratumID)])
  lcelec <- electivity(lc, pa)
  ucelec <- electivity(uc, pa)
  elec <- electivity(pd, pa) %>%
    mutate(lcl = lcelec$value,
           ucl = ucelec$value,
           hab_typ = strsplit(name, " "),
           habitat = unlist(lapply(hab_typ, function(x) x[1])),
           type = unlist(lapply(hab_typ, function(x) x[2])))
  
  cols <- switch(hb,
                 forest = pal[1:2],
                 heath = pal[3:4],
                 meadow = pal[5:6],
                 sand = pal[7:8])
  plt <- elec %>%
    filter(habitat==hb) %>%
    ggplot(aes(x = time, y = value)) +
    geom_line(aes(color = name, linetype = type)) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl, fill=name), alpha = 0.2) +
    geom_hline(yintercept = 0, color="grey") +
    scale_color_manual(values = cols) + 
    scale_fill_manual(values = cols) +
    scale_x_continuous(limits=c(0, 24), breaks=seq(0,24,6)) +
    scale_y_continuous(limits=c(-1, 1), breaks=seq(-1,1,1)) +
    labs(x = NULL, y = NULL) + 
    theme_classic()
  plt + theme(legend.position = "none")
}

# Create plot list
spp_hab <- expand.grid(sp=species, hb=unique(act_strata$habitat),
                       stringsAsFactors=FALSE)
sel_plots <- lapply(1:nrow(spp_hab), function(i) 
  sel_plot(spp_hab$sp[i], spp_hab$hb[i]))
sel_plots <- plot_grid(plotlist=sel_plots, nrow=7, byrow=FALSE)

# Create plots for marginal labels
sp_labs <- lapply(sub("_", " ", species), function(sp)
  ggplot(mapping=aes(x=0, y=0, label=sp)) +
    geom_text(angle=270) +
    theme_void())
sp_labs <- plot_grid(plotlist=sp_labs, ncol=1)
hb_labs <- lapply(unique(act_strata$habitat), function(hb)
  ggplot(mapping=aes(x=0, y=0, label=hb)) +
    geom_text() +
    theme_void())
hb_labs <- plot_grid(plotlist=hb_labs, nrow=1)
x_lab <- ggplot(mapping=aes(x=0, y=0, label="Time (h)")) +
  geom_text() +
  theme_void()
y_lab <- ggplot(mapping=aes(x=0, y=0, label="Selection index")) +
  geom_text(angle=90) +
  theme_void()
blank <- ggplot() + theme_nothing()

# Stich plots together and plot
row1 <- plot_grid(blank, hb_labs, blank, nrow=1,
                  rel_widths = c(1,14,1))
row2 <- plot_grid(y_lab, sel_plots, sp_labs, nrow=1,
                  rel_widths = c(1,16,1))
row3 <- plot_grid(blank, x_lab, blank, nrow=1,
                  rel_widths = c(1,16,1))
plot_grid(row1, row2, row3, ncol=1,
          rel_heights = c(1, 16, 1))


## Stratum activity levels ####
alevel_plot <- function(sp){
  str_actdat <- with(actmods[[sp]], 
                     data.frame(est = est$act_stratum,
                                lcl = ci$act_stratum[1,],
                                ucl = ci$act_stratum[2,])) %>%
    rownames_to_column("stratum")
  
  ggplot(str_actdat, aes(x=stratum, y=est, fill=stratum)) + 
    geom_bar(stat="identity", color=pal, 
             position=position_dodge()) +
    geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.2,
                  position=position_dodge(.9)) +
    labs(x = NULL, y = NULL, title = sp) +
    coord_cartesian(ylim=c(0.2, 0.7)) +
    theme_classic()
}
alevel_plots <- lapply(species, alevel_plot)
alevel_plots <- alevel_plots %>%
  lapply(function(plt) 
    plt + theme(legend.position = "none",
                axis.text.x = element_blank()))
alevel_plots <- plot_grid(plotlist=alevel_plots)
y_lab <- ggplot(mapping=aes(x=0, y=0, label="Activity level")) +
  geom_text(angle=90) +
  theme_void()
plot_grid(y_lab, alevel_plots, lgnd, nrow=1,
          rel_widths = c(1, 15, 3))
