## Read Agouti datapackage ####
dir <- "C:/Users/rowcliffe.m/OneDrive - Zoological Society of London/CameraTrapping/Activity/Stratified activity/Hoge Veluwe data/hoge-veluwe-wildlife-monitoring-project-20240816202112"
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

dir.create("HogeVeluweData")
ctdpbound$data$deployments %>%
  write.csv("./HogeVeluweData/deployments.csv", row.names = FALSE)
ctdpbound$data$observations %>%
  write.csv("./HogeVeluweData/observations.csv", row.names = FALSE)
