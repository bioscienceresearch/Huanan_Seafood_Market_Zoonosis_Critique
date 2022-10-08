rm(list=ls(all=TRUE))

library(spatstat)

# We use UTM Zone 50 reprojected data from Worobey et al. (2022), and digitized data from the Joint WHO-China study (2021) and Gao et al. (2022)
setwd("/mnt/1TB_0/Data/Code/code/Huanan_Market_Analysis/R")

gdf = st_read("/mnt/1TB_0/Data/GIS/Wuhan_Early_Cases/geojson/UTM_Zone50N/huanan-market-internal_poly_UTM_Z50N.geojson")
gdf_data = st_read("/mnt/1TB_0/Data/GIS/Wuhan_Early_Cases/geojson/UTM_Zone50N/huanan-market-internal_data_UTM_Z50N.geojson")
hsm_west_rect = st_read("/mnt/1TB_0/Data/GIS/Wuhan_Early_Cases/geojson/UTM_Zone50N/HSM_west_internal_border_UTM_Z50N.geojson")
hsm_east_rect = st_read("/mnt/1TB_0/Data/GIS/Wuhan_Early_Cases/geojson/UTM_Zone50N/HSM_east_internal_border_UTM_Z50N.geojson")
wildlife_pts= st_read("/mnt/1TB_0/Data/GIS/Wuhan_Early_Cases/geojson/UTM_Zone50N/Wildlife_pt_UTM_Z50N.geojson")
wildlife_poly=st_read("/mnt/1TB_0/Data/GIS/Wuhan_Early_Cases/geojson/Gao_wildlife.geojson")
# convert wildlife polygons to UTM Zone50N
wildlife_poly = st_transform(wildlife_poly, 32650)

tail(gdf_data)
cases = gdf_data[gdf_data$group=='HumanCase',]
market_boundary = gdf[gdf$group=='Boundary',]

# West and East
plot(market_boundary['geometry'])
plot(cases, add=T)
# Just West
#plot(hsm_west_rect, col = alpha('black', 0.1))
plot(market_boundary[1,'geometry'])

cases_coords <- matrix(unlist(cases$geometry), ncol = 2, byrow = T)

# create window
window <- as.owin(market_boundary)
summary(window)

dim(cases_coords)
# add random jitter to case locations so that they are not located on top of eachother
jitterx <- matrix(unlist(runif(dim(cases_coords)[1], min=-0.5, max=0.5)), ncol = 1, byrow = T)
jittery <- matrix(unlist(runif(dim(cases_coords)[1], min=-0.5, max=0.5)), ncol = 1, byrow = T)
xj<- cases_coords[,1]+jitterx
yj<- cases_coords[,2]+jittery

cases_ppp <- ppp(x = xj, y = yj, window = window)
plot(cases_ppp)

# G(r) 
# Values G(r) > G pois (r) suggest that nearest neighbour distances in the point pattern are shorter than for a Poisson process, 
# suggesting a clustered pattern; while values G(r) < G pois (r) suggest a regular (inhibited) pattern

par(pty = "s")
plot(Gest(cases_ppp))

fisher <- function(x) {asin(sqrt(x))}
plot(Gest(cases_ppp), fisher(.) ~ fisher(theo))

par(pty = "s")
plot(Kest(cases_ppp))

L <- Lest(cases_ppp)
plot(L, main = "L function")

# pair correlation function
plot(pcf(cases_ppp))

plot(allstats(cases_ppp))

#Matern cluster model for poisson process
fitM <- kppm(cases_ppp, ~1, "MatClust")
plot(fitM)
plot(envelope(fitM, Lest, nsim = 39))
plot(envelope(fitM, Fest, nsim = 39))
plot(envelope(fitM, Gest, nsim = 39))
plot(envelope(fitM, Jest, nsim = 39))
plot(envelope(fitM, Kest, nsim = 39))

#Thomas model
fitT <- kppm(cases_ppp, ~1, "Thomas")
plot(envelope(fitT, Lest, nsim = 39))
plot(envelope(fitT, Fest, nsim = 39))
plot(envelope(fitT, Gest, nsim = 39))
plot(envelope(fitT, Jest, nsim = 39))
plot(envelope(fitT, Kest, nsim = 39))

# Fitting cluster models using the pair correlation
fitp <- kppm(cases_ppp, ~1, "Thomas", statistic = "pcf")
plot(fitp)

plot(density(cases_ppp, 10))
plot(cases_ppp, add=TRUE)

D <- density(cases_ppp)
plot(D)
persp(D)


# observed and expected quadrant count, may be better to rotate 15.2 degrees
M <- quadrat.test(cases_ppp, nx = 6, ny = 6)
M
plot(cases_ppp)
plot(M, add = TRUE, cex = .5)
M$p.value


# Second order properties
# QC wildlife
plot(market_boundary['geometry'])
plot(wildlife_pts, add=T)

wildlife_coords <- matrix(unlist(wildlife_pts$geometry), ncol = 2, byrow = T)
wildlife_ppp = ppp(x = wildlife_coords[,1], y = wildlife_coords[,2],
                   window = window, check = T)


plot(density(wildlife_ppp, 10))
D <- density(wildlife_ppp)
plot(D)
persp(D)


# create a second order dataset
sup_Cw<-superimpose(A=cases_ppp, B=wildlife_ppp)

sup_Cw
plot(sup_Cw)
plot(density(split(sup_Cw),10), ribbon = FALSE)


# G function A, B Matrix
plot(alltypes(sup_Cw, "G"))

# Distance methods 
plot(Gcross(sup_Cw, "A", "B"))

# One type to any type 
plot(Gdot(sup_Cw, "A"))

# Mark connection function
markconnect(sup_Cw, "A", "B")
plot(alltypes(sup_Cw, markconnect))

# Mark equality function
plot(markcorr(sup_Cw))

# Poisson null 
#pdf('Poinsson_model_Kcross_wildlife_cases.pdf')
E <- envelope(sup_Cw, Kcross, nsim = 39, i = "A", j = "B")
plot(E, main = "test of marked Poisson model")

E <- envelope(sup_Cw, Gcross, nsim = 39, i = "A", j = "B")
plot(E, main = "test of marked Poisson model")

E <- envelope(sup_Cw, Lcross, nsim = 39, i = "A", j = "B")
plot(E, main = "test of marked Poisson model")

E <- envelope(sup_Cw, pcfcross, nsim = 39, i = "A", j = "B")
plot(E, main = "pcf test of marked Poisson model")

E <- envelope(sup_Cw, Jcross, nsim = 139, i = "A", j = "B")
plot(E, main = "J i-to-j test of marked Poisson model")


E <- envelope(sup_Cw, Kcross, nsim = 39, i = "B", j = "A")
plot(E, main = "Kcross test of marked Poisson model")

E <- envelope(sup_Cw, Gcross, nsim = 139, i = "B", j = "A")
plot(E, main = "Gcross test of marked Poisson model")

E <- envelope(sup_Cw, Lcross, nsim = 39, i = "B", j = "A")
plot(E, main = "Lcross test of marked Poisson model")

E <- envelope(sup_Cw, pcfcross, nsim = 39, i = "B", j = "A")
plot(E, main = "pcf test of marked Poisson model")

E <- envelope(sup_Cw, Jcross, nsim = 139, i = "B", j = "A")
plot(E, main = "Jcross test of marked Poisson model")


# Arrays of envelopes
aE <- alltypes(sup_Cw, Kcross, nsim = 39, envelope = TRUE)
plot(aE, sqrt(./pi) - r ~ r, ylab = "L(r)-r")

# Repeat for western section only

cases_west_gdf = st_read("/mnt/1TB_0/Data/GIS/Wuhan_Early_Cases/geojson/UTM_Zone50N/HSM_cases_West_UTM_Z50N.geojson")
tail(cases_west_gdf)
cases = cases_west_gdf[cases_west_gdf$group=='HumanCase',]

plot(market_boundary[1,'geometry'])
window <- as.owin(market_boundary[1,'geometry'])
summary(window)
cases

cases_coords <- matrix(unlist(cases$geometry), ncol = 2, byrow = T)

# add random jitter to case locations so that they are not located on top of eachother
jitterx <- matrix(unlist(runif(dim(cases_coords)[1], min=-0.5, max=0.5)), ncol = 1, byrow = T)
jittery <- matrix(unlist(runif(dim(cases_coords)[1], min=-0.5, max=0.5)), ncol = 1, byrow = T)
xj<- cases_coords[,1]+jitterx
yj<- cases_coords[,2]+jittery
cases_ppp <- ppp(x = xj, y = yj, window = window)
plot(cases_ppp)

par(pty = "s")
plot(Gest(cases_ppp))

fisher <- function(x) {asin(sqrt(x))}
plot(Gest(cases_ppp), fisher(.) ~ fisher(theo))

par(pty = "s")
plot(Kest(cases_ppp))

L <- Lest(cases_ppp)
plot(L, main = "L function")

# pair correlation function
plot(pcf(cases_ppp))

plot(allstats(cases_ppp))

fitM <- kppm(cases_ppp, ~1, "MatClust")
plot(fitM)
plot(envelope(fitM, Lest, nsim = 39))
plot(envelope(fitM, Fest, nsim = 39))
plot(envelope(fitM, Gest, nsim = 39))
plot(envelope(fitM, Jest, nsim = 39))
plot(envelope(fitM, Kest, nsim = 39))

fit <- kppm(cases_ppp, ~1, "Thomas")
plot(envelope(fit, Lest, nsim = 39))
plot(envelope(fitM, Fest, nsim = 39))
plot(envelope(fitM, Gest, nsim = 39))
plot(envelope(fitM, Jest, nsim = 39))
plot(envelope(fitM, Kest, nsim = 39))

# Fitting cluster models using the pair correlation
fitp <- kppm(cases_ppp, ~1, "Thomas", statistic = "pcf")
plot(fitp)

# Second order properties

wildlife_ppp = ppp(x = wildlife_coords[,1], y = wildlife_coords[,2],
                   window = window, check = T)
plot(wildlife_ppp)
# create a second order dataset
sup_Cw<-superimpose(A=cases_ppp, B=wildlife_ppp)
sup_Cw
plot(sup_Cw)

plot(alltypes(sup_Cw, "G"))

# Distance methods 
plot(Gcross(sup_Cw, "A", "B"))
plot(Gcross(sup_Cw, "B", "A"))

# One type to any type 
plot(Gdot(sup_Cw, "A"))

# Mark connection function
markconnect(sup_Cw, "A", "B")
plot(alltypes(sup_Cw, markconnect))

# Mark equality function
plot(markcorr(sup_Cw))

# Poisson null 

E <- envelope(sup_Cw, Kcross, nsim = 39, i = "B", j = "A")
plot(E, main = "Kcross test of marked Poisson model")

E <- envelope(sup_Cw, Gcross, nsim = 139, i = "B", j = "A")
plot(E, main = "Gcross test of marked Poisson model")

E <- envelope(sup_Cw, Lcross, nsim = 39, i = "B", j = "A")
plot(E, main = "Lcross test of marked Poisson model")

E <- envelope(sup_Cw, pcfcross, nsim = 39, i = "B", j = "A")
plot(E, main = "pcf test of marked Poisson model")

E <- envelope(sup_Cw, Jcross, nsim = 139, i = "B", j = "A")
plot(E, main = "Jcross test of marked Poisson model")

# Arrays of envelopes
aE <- alltypes(sup_Cw, Kcross, nsim = 39, envelope = TRUE)
plot(aE, sqrt(./pi) - r ~ r, ylab = "L(r)-r")

#### Repeat for cases from 13th-20th only (which only occurr on the West side)

cases_to_20th_gdf = st_read("/mnt/1TB_0/Data/GIS/Wuhan_Early_Cases/geojson/UTM_Zone50N/WHO_HSM_vendor_cases_1_2019-12-13th-20th_UTM_Z50N.geojson")
cases_to_20th_gdf

plot(market_boundary[1,'geometry'])
window <- as.owin(market_boundary[1,'geometry'])
summary(window)
plot(cases_to_20th_gdf, add=T)

cases_coords <- matrix(unlist(cases_to_20th_gdf$geometry), ncol = 2, byrow = T)

# add random jitter to case locations so that they are not located on top of eachother
jitterx <- matrix(unlist(runif(dim(cases_coords)[1], min=-0.5, max=0.5)), ncol = 1, byrow = T)
jittery <- matrix(unlist(runif(dim(cases_coords)[1], min=-0.5, max=0.5)), ncol = 1, byrow = T)
xj<- cases_coords[,1]+jitterx
yj<- cases_coords[,2]+jittery
cases_ppp <- ppp(x = xj, y = yj, window = window)
plot(market_boundary[1,'geometry'])
plot(cases_ppp, add=T)

par(pty = "s")
plot(Gest(cases_ppp))

fisher <- function(x) {asin(sqrt(x))}
plot(Gest(cases_ppp), fisher(.) ~ fisher(theo))

par(pty = "s")
plot(Kest(cases_ppp))

L <- Lest(cases_ppp)
plot(L, main = "L function")

# pair correlation function
plot(pcf(cases_ppp))

plot(allstats(cases_ppp))

fitM <- kppm(cases_ppp, ~1, "MatClust")
plot(fitM)
plot(envelope(fitM, Lest, nsim = 39))
plot(envelope(fitM, Fest, nsim = 39))
plot(envelope(fitM, Gest, nsim = 39))
plot(envelope(fitM, Jest, nsim = 39))
plot(envelope(fitM, Kest, nsim = 39))

fit <- kppm(cases_ppp, ~1, "Thomas")
plot(envelope(fit, Lest, nsim = 39))
plot(envelope(fitM, Fest, nsim = 39))
plot(envelope(fitM, Gest, nsim = 39))
plot(envelope(fitM, Jest, nsim = 39))
plot(envelope(fitM, Kest, nsim = 39))

# Fitting cluster models using the pair correlation
fitp <- kppm(cases_ppp, ~1, "Thomas", statistic = "pcf")
plot(fitp)

# Second order properties

wildlife_ppp = ppp(x = wildlife_coords[,1], y = wildlife_coords[,2],
                   window = window, check = T)
plot(wildlife_ppp)
# create a second order dataset
sup_Cw<-superimpose(A=cases_ppp, B=wildlife_ppp)
sup_Cw
plot(sup_Cw)

plot(alltypes(sup_Cw, "G"))

# Distance methods 
plot(Gcross(sup_Cw, "A", "B"))
plot(Gcross(sup_Cw, "B", "A"))

# One type to any type 
plot(Gdot(sup_Cw, "A"))

# Mark connection function
markconnect(sup_Cw, "A", "B")
plot(alltypes(sup_Cw, markconnect))

# Mark equality function
plot(markcorr(sup_Cw))

# Poisson null 

E <- envelope(sup_Cw, Kcross, nsim = 39, i = "B", j = "A")
plot(E, main = "Kcross test of marked Poisson model")

E <- envelope(sup_Cw, Gcross, nsim = 139, i = "B", j = "A")
plot(E, main = "Gcross test of marked Poisson model")

E <- envelope(sup_Cw, Lcross, nsim = 39, i = "B", j = "A")
plot(E, main = "Lcross test of marked Poisson model")

E <- envelope(sup_Cw, pcfcross, nsim = 39, i = "B", j = "A")
plot(E, main = "pcf test of marked Poisson model")

E <- envelope(sup_Cw, Jcross, nsim = 139, i = "B", j = "A")
plot(E, main = "Jcross test of marked Poisson model")



