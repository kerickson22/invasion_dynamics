library(rgbif)
library(sp)
library(maps)
library(maptools)
library(rgdal)
library(omnibus)
library(enmSdm)
library(bayesLopod)
library(calibrate)
library(scrubr)
library(plyr)
library(purrr)
library(magick)
library(classInt)
library(gridExtra)

setwd("C:/Users/kerickson/Documents/GitHub/invasion_dynamics")
#load(".RData")

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
memory.limit(memory.limit() *2^30)

#Turn off CAR
mLopodShape_1930_basic = modelLopod(LopodData = ld_Shape_1930,
                              varP = T,
                              q = NULL,
                              pmin = 0,
                              CAR = F,
                              nChains = 3,
                              warmup = 2000,
                              sampling = 1000,
                              nCores =3)
lopodTrace( mLopodShape_1930_basic, inc_warmup = T)

mLopodShape_1960 = modelLopod(LopodData = ld_Shape_1960,
                              varP = T,
                              q = NULL,
                              pmin = 0,
                              CAR = T,
                              nChains = 3,
                              warmup = 500,
                              sampling = 3000,
                              nCores =3)
lopodTrace(mLopodShape_1960, inc_warmup=T)


#Moving to 1970, where there are 68 detections (i = 10)

mLopodShape_1970_basic = modelLopod(LopodData = ld_Shape_1970,
                              varP = F,
                              q = NULL,
                              pmin = 0,
                              CAR =F,
                              nChains = 3,
                              warmup = 2000,
                              sampling = 1000,
                              nCores =4)
lopodTrace(mLopodShape_1970_basic, inc_warmup=T)

sampler_params <- get_sampler_params(mLopodShape_1970@StanFit, inc_warmup=FALSE)
divergent<-sapply(sampler_params, function(x) sum(x[, "divergent__"]))
print(divergent)
#no divergent transitions

mLopodShape_1980 = modelLopod(LopodData = ld_Shape_1980,
                              varP = T,
                              q = NULL,
                              pmin = 0,
                              CAR = T,
                              nChains = 3,
                              warmup = 500,
                              sampling = 1000,
                              nCores =3)
lopodTrace(mLopodShape_1980, inc_warmup=T)
sampler_params <- get_sampler_params(mLopodShape_1980@StanFit, inc_warmup=FALSE)
divergent<-sapply(sampler_params, function(x) sum(x[, "divergent__"]))
print(divergent)
#no divergent transitions



mLopodShape_1990 = modelLopod(LopodData = ld_Shape_1990,
                              varP = T,
                              q = NULL,
                              pmin = 0,
                              CAR = T,
                              nChains = 3,
                              warmup = 4000,
                              sampling = 2000,
                              nCores =3)
lopodTrace(mLopodShape_1990, inc_warmup=T)
sampler_params <- get_sampler_params(mLopodShape_1990@StanFit, inc_warmup=FALSE)
divergent<-sapply(sampler_params, function(x) sum(x[, "divergent__"]))
print(divergent)
#no divergent transitions


mLopodShape_2000 = modelLopod(LopodData = ld_Shape_2000,
                              varP = T,
                              q = NULL,
                              pmin = 0,
                              CAR = T,
                              nChains = 3,
                              warmup = 2000,
                              sampling = 2000,
                              nCores =3)
lopodTrace(mLopodShape_2000, inc_warmup=T)
#There were 10 divergent transitions after warmup

sampler_params <- get_sampler_params(mLopodShape_2000@StanFit, inc_warmup=FALSE)
divergent<-sapply(sampler_params, function(x) sum(x[, "divergent__"]))
print(divergent)

mLopodShape_2010 = modelLopod(LopodData = ld_Shape_2010,
                              varP = T,
                              q = NULL,
                              pmin = 0,
                              CAR = T,
                              nChains = 3,
                              warmup = 500,
                              sampling = 1000,
                              nCores =3)
lopodTrace(mLopodShape_2010, inc_warmup=T)
sampler_params <- get_sampler_params(mLopodShape_2010@StanFit, inc_warmup=FALSE)
divergent<-sapply(sampler_params, function(x) sum(x[, "divergent__"]))
print(divergent)
#No divergent transitions


#Have the same issue that I ran into before...so forced sampling effort in
# Escambia county to be zero
mLopodShape_2020 = modelLopod(LopodData = ld_Shape_2020,
                              varP = T,
                              q = NULL,
                              pmin = 0,
                              CAR = T,
                              nChains = 3,
                              warmup = 1000,
                              sampling = 1000,
                              nCores =3)
lopodTrace(mLopodShape_2020, inc_warmup=T)
sampler_params <- get_sampler_params(mLopodShape_2020@StanFit, inc_warmup=FALSE)
divergent<-sapply(sampler_params, function(x) sum(x[, "divergent__"]))
print(divergent)
#No divergent transitions



SchinusShape_1970_basic = mLopodShape_1970_basic@LopodData@geoDataObject
psiShape_1970_basic = lopodShape(mLopodShape_1970_basic, "psi_Sampled", extrapolate = F,  quant = 0.05)
SchinusShape_1970_basic@data[,"psi_05"] = psiShape_1970_basic@data[,"psi_Sampled"]
psiShape_1970_basic = lopodShape(mLopodShape_1970_basic, "psi_Sampled", extrapolate = F,  quant = 0.5)
SchinusShape_1970_basic@data[,"psi_50"] = psiShape_1970_basic@data[,"psi_Sampled"]
psiShape_1970_basic = lopodShape(mLopodShape_1970_basic, "psi_Sampled", extrapolate = F,  quant = 0.95)
SchinusShape_1970_basic@data[,"psi_95"] = psiShape_1970_basic@data[,"psi_Sampled"]
spplot( SchinusShape_1970_basic,
        zcol = c("psi_05", "psi_50", "psi_95"),
        names.attr = c("Psi (5% quantile)", "Psi (median)", "Psi (95% quantile)"), main = "Occupancy (Psi): 1970")
