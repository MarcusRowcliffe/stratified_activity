# ANALYSIS: ACTIVITY/HABITAT USE ####
## Read Agouti datapackage ####
dir <- "./Hoge Veluwe data/hoge-veluwe-wildlife-monitoring-project-20240816202112"
ctdpdat <- read_camtrap_dp(file.path(dir, "datapackage.json"), FALSE)

# Check out parsing problems
#probs <- frictionless::problems(ctdpdat$data$observations)
#obs <- read.csv(file.path(dir, "observations.csv"))
#View(ctdpdat$data$observations[probs$row,])
#sum(obs[probs$row,]$observationType=="animal", na.rm=T)
# some malformed dates, but not enough to worry about - only 4 animal observations

# Remove some deployments:
# - missing_dep has no animal observations and is missing from deployments
# - locationName codes V## (not in study design)
missing_dep <- "51f1d2cf-403b-49ff-8cf1-3911db9dcb2c"
ctdpdat$data$media <- data.frame(deploymentID=0) # need this to make next line run
ctdpdat <- ctdpdat %>%
  subset_deployments(deploymentID != missing_dep &
                       nchar(locationName) > 3)

# Standardise location names; add stratumID
fix_locationNames <- function(x){
  ifelse(grepl("OB", x), 
         paste0(substr(x, 4, 5),
                substr(x,1,2), "-", 
                substr(x, nchar(x)-1, nchar(x))),
         ifelse(nchar(x) == 11, 
                substr(x, 1, 7),
                x))
}
ctdpdat$data$deployments <- ctdpdat$data$deployments %>%
  mutate(locationName = fix_locationNames(locationName),
         habitatNum = substr(locationName, 2, 2),
         typeNum = substr(locationName, 4, 4),
         habitat = case_when(
           habitatNum == "1" ~ "sand",
           habitatNum %in% c("2", "3") ~ "heath",
           habitatNum %in% c("4", "5") ~ "forest",
           habitatNum == "6" ~ "meadow"),
         type = case_when(
           typeNum == "1" ~ "rest",
           typeNum == "0" ~ "public",
           typeNum == "B" ~ "OB"),
         stratumID = paste(habitat, type))

# Filter species, add common names and time of day
table(ctdpdat$data$observations$scientificName) %>%
  sort()
spp <- c("Capreolus capreolus", "Cervus elaphus", "Ovis ammon", "Sus scrofa",
         "Vulpes vulpes", "Lepus europaeus", "Meles meles")
species <- c("roe_deer", "red_deer", "mouflon", "wild_boar",
             "red_fox", "brown_hare", "badger")
ctdpdat$data$observations <- ctdpdat$data$observations %>%
  mutate(scientificName = ifelse(scientificName=="Leporidae", 
                                 "Lepus europaeus", scientificName)) %>%
  filter(scientificName %in% spp) %>%
  mutate(commonName = species[match(scientificName, spp)],
         time = time_of_day(timestamp))
#ctdpdat$data$observations %>%
#  filter(is.na(time)) %>%
#  View()
# A handful of malformed dates, hopefully negligible

# slice up the data into annual chunks
starts <- paste0(2014:2020, "-03-01 00:00:00")
ends <- paste0(2014:2020, "-05-31 00:00:00")
ctdpdat$data <- ctdpdat$data[-2] # remove media table to make next step run
ctdpslices <- lapply(1:length(starts), function(i) 
  slice_camtrap_dp(ctdpdat, starts[i], ends[i]))
#lapply(ctdpslices, plot_deployment_schedule)

# bind sliced ctdp list in a single datapackage
ctdpbound <- ctdpdat
ctdpbound$data$deployments <- ctdpslices %>%
  lapply(function(d) d$data$deployments) %>%
  bind_rows() %>%
  mutate(year = year(start),
         locationYear = paste(locationName, year, sep="-"))
ctdpbound$data$observations <- ctdpslices %>%
  lapply(function(d) d$data$observations) %>%
  bind_rows() %>%
  mutate(year = year(timestamp))
table(ctdpbound$data$observations$commonName) %>%
  sort()

## Read stratum data ####
dir <- "./Hoge Veluwe data"
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
         stratumID = paste(habitat, type)) %>%
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
dep <- ctdpbound$data$deployments %>%
  mutate(stratumID = sub("OB", "public", stratumID),
         effort = as.numeric(difftime(end, start, units="days")))
# function fits stratified activity model for a given species
fitfunc <- function(sp, reps){
  if(sp %in% species[5:7]){
    spd <- ang <- rad <- NULL
    str <- act_strata %>%
      filter(grepl("roe_deer", species,))
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
    str <- act_strata %>%
      filter(grepl(sp, species,))
    spd <- pmat_from_func(spdFunc, str)
    ang <- pmat_from_func(angFunc, str)
    rad <- pmat_from_func(radFunc, str)
  }
  obs <- subset(ctdpbound$data$observations, commonName==sp)
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
sp <- "roe_deer"
popdist <- actmods[[sp]]$est$popdist
act <- actmods[[sp]]$est$pdf

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

popact_plots <- lapply(species, popact_plot, actmods, pal)
popact_plots <- plot_grid(plotlist = popact_plots, ncol=1)

xlab <- ggplot(mapping=aes(x=0, y=0, label="Time (h)")) +
  geom_text() +
  theme_void()
ylab1 <- ggplot(mapping=aes(x=0, y=0, label="Activity PDF")) +
  geom_text(angle=90) +
  theme_void()
ylab2 <- ggplot(mapping=aes(x=0, y=0, label="Population proportion")) +
  geom_text(angle=90) +
  theme_void()
ylab3 <- ggplot(mapping=aes(x=0, y=0, label="Population activity")) +
  geom_text(angle=90) +
  theme_void()
lgnd <- plot_grid(legend)

plots <- plot_grid(ylab1, act_plots,
                   ylab2, popdist_plots,
                   ylab3, popact_plots, 
                   lgnd,
                   nrow=1, rel_widths = c(rep(c(1,7),3),4))
#                    nrow=1, rel_widths = c(rep(c(1,7),2),5))
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
sel_plots <- plot_grid(plotlist=sel_plots, nrow=4)

# Create plots for marginal labels
sp_labs <- lapply(sub("_", " ", species), function(sp)
  ggplot(mapping=aes(x=0, y=0, label=sp)) +
    geom_text() +
    theme_void())
sp_labs <- plot_grid(plotlist=sp_labs, nrow=1)
hb_labs <- lapply(unique(act_strata$habitat), function(hb)
  ggplot(mapping=aes(x=0, y=0, label=hb)) +
    geom_text(angle=270) +
    theme_void())
hb_labs <- plot_grid(plotlist=hb_labs, ncol=1)
x_lab <- ggplot(mapping=aes(x=0, y=0, label="Time (h)")) +
  geom_text() +
  theme_void()
y_lab <- ggplot(mapping=aes(x=0, y=0, label="Selection index")) +
  geom_text(angle=90) +
  theme_void()
blank <- ggplot() + theme_nothing()

# Stich plots together and plot
row1 <- plot_grid(blank, sp_labs, blank, nrow=1,
                  rel_widths = c(1,14,1))
row2 <- plot_grid(y_lab, sel_plots, hb_labs, nrow=1,
                  rel_widths = c(1,16,1))
row3 <- plot_grid(blank, x_lab, blank, nrow=1,
                  rel_widths = c(1,16,1))
dev.off()
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
