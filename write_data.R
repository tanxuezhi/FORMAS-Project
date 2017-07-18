library(rio)
library(dplyr)

lc_class <- rio::import(file = "../Landcover/clc_legend.xls", which = 1L)

###################
### butterflies ###
###################

# sites data
sites_NL <- read.csv(file = "../Data/Butterflies - Netherlands/Sites_NL_ETRS89_landcover.csv")
sites_FIN <- read.csv(file = "../Data/Butterflies - Finland/Sites_FIN_ETRS89_landcover.csv")

# fragmentation data
frag_data_NL <- read.csv("../Connectivity/Fragmentation/NL/Frag_indices.csv")
frag_data_FIN <- read.csv("../Connectivity/Fragmentation/FIN/Frag_indices.csv")

sites_NL <- merge(sites_NL, frag_data_NL, by.x = "Site", by.y = "LID")
sites_NL <- merge(sites_NL, lc_class, by.x = "Landcover", by.y = "CLC_CODE")
sites_NL <- sites_NL[,-which(names(sites_NL) %in% "GRID_CODE")]

sites_FIN <- merge(sites_FIN, frag_data_FIN, by.x = "Site", by.y = "LID")
sites_FIN <- merge(sites_FIN, lc_class, by.x = "Landcover", by.y = "GRID_CODE")
sites_FIN <- sites_FIN[,-which(names(sites_FIN) %in% "CLC_CODE")]

# cti data
# Abundance-weighted #
cti_AB_NL <- read.csv2("../Data/Butterflies - Netherlands/cti_site_year_abundance_2017-05-19.csv", dec = ".")[,-3]
cti_AB_NL <- merge(cti_AB_NL, sites_NL, by.x = "Site", by.y = "Site")

cti_AB_FIN <- read.csv("../Data/Butterflies - Finland/CTI_Abundance_FINLAND_1999-2016.csv", dec = ".")
cti_AB_FIN <- merge(cti_AB_FIN, sites_FIN, by.x = "Site", by.y = "Site")

cti_AB_butterflies <- as.tbl(bind_rows(cbind.data.frame(cti_AB_NL, country = "NL"), cbind.data.frame(cti_AB_FIN, country = "FIN")))
cti_AB_butterflies$Site <- paste0(cti_AB_butterflies$Site, "_", cti_AB_butterflies$country)

# P/A #
cti_PA_NL <- read.csv2("../Data/Butterflies - Netherlands/cti_site_year_presence_2017-05-19.csv", dec = ".")[,-3]
cti_PA_NL <- merge(cti_PA_NL, sites_NL, by.x = "Site", by.y = "Site")

cti_PA_FIN <- read.csv("../Data/Butterflies - Finland/CTI_presence_FINLAND_1999-2016.csv", dec = ".")
cti_PA_FIN <- merge(cti_PA_FIN, sites_FIN, by.x = "Site", by.y = "Site")

cti_PA_butterflies <- as.tbl(bind_rows(cbind.data.frame(cti_PA_NL, country = "NL"), cbind.data.frame(cti_PA_FIN, country = "FIN")))
cti_PA_butterflies$Site <- paste0(cti_PA_butterflies$Site, "_", cti_PA_butterflies$country)

# merge and write
cti_butterflies <- as.tbl(bind_rows(cbind.data.frame(type = "Presence", cti_PA_butterflies), cbind.data.frame(type = "Abundance", cti_AB_butterflies)))
write.csv(cti_butterflies, "../Data/cti_butterflies_data.csv", row.names = F)

#############
### birds ###
#############

# sites data
sites_SWE <- read.csv(file = "../Data/Birds - Sweden/Sites_SWE_ETRS89_landcover.csv")

# fragmentation data
frag_data_SWE <- read.csv("../Connectivity/Fragmentation/SWE/Frag_indices.csv")

sites_SWE <- merge(sites_SWE, frag_data_SWE, by.x = "Site", by.y = "LID")
sites_SWE <- merge(sites_SWE, lc_class, by.x = "Landcover", by.y = "GRID_CODE")
sites_SWE <- sites_SWE[,-which(names(sites_SWE) %in% "CLC_CODE")]

# cti data
# Abundance-weighted #
cti_AB_SWE <- read.csv("../Data/Birds - Sweden/CTI_abundance_Sweden_1996-2016.csv", dec = ".")
cti_AB_SWE <- as.tbl(inner_join(cti_AB_SWE, sites_SWE, by = "Site"))

# P/A #
cti_PA_SWE <- read.csv("../Data/Birds - Sweden/CTI_presence_Sweden_1996-2016.csv", dec = ".")
cti_PA_SWE <- as.tbl(inner_join(cti_PA_SWE, sites_SWE, by = "Site"))

# merge and write
cti_birds <- as.tbl(bind_rows(cbind.data.frame(type = "Presence", cti_PA_SWE), cbind.data.frame(type = "Abundance", cti_AB_SWE)))
write.csv(cti_birds, "../Data/cti_birds_data.csv", row.names = F)


