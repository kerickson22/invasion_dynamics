1870
2020
florida_df[[1]] #1880
florida_df[[2]] #1890
florida_df[[3]] #1900
florida_df[[4]] #1910
florida_df[[5]] #1920
florida_df[[6]] #1930
florida_df[[7]] #1940
florida_df[[8]] #1950
florida_df[[9]] #1960
florida_df[[10]] #1970
florida_df[[11]] #1980
florida_df[[12]] #1990
florida_df[[13]] #2000
florida_df[[14]] #2010
florida_df[[15]] #2020

new_recs_1880 <-florida_df[[1]]
new_recs_1890 <-florida_df[[2]]
new_recs_1890$numAnacard <- new_recs_1890$numAnacard - new_recs_1880$numAnacard
new_recs_1890$numSchinus <- new_recs_1890$numSchinus - new_recs_1890$numSchinus
#.... 

new_recs_1970 <- florida_df[[10]]
new_recs_1970$numAnacard <- florida_df[[10]]$numAnacard - florida_df[[9]]$numAnacard
new_recs_1970$numSchinus <- florida_df[[10]]$numSchinus - florida_df[[10]]$numSchinus

new_recs_1980 <- florida_df[[11]]
new_recs_1980$numAnacard <- florida_df[[11]]$numAnacard - florida_df[[10]]$numAnacard
new_recs_1980$numSchinus <- florida_df[[11]]$numSchinus - florida_df[[10]]$numSchinus

new_recs_1990 <- florida_df[[12]]
new_recs_1990$numAnacard <- florida_df[[12]]$numAnacard - florida_df[[11]]$numAnacard
new_recs_1990$numSChinus <- florida_df[[12]]$numSchinus - florida_df[[11]]$numSchinus

new_recs_2000 <- florida_df[[13]]
new_recs_2000$numAnacard <- florida_df[[13]]$numAnacard - florida_df[[12]]$numAnacard
new_recs_2000$numSchinus <- florida_df[[13]]$numSchinus - florida_df[[12]]$numSchinus

new_recs_2010 <- florida_df[[14]]
new_recs_2010$numAnacard <- florida_df[[14]]$numAnacard - florida_df[[13]]$numAnacard
new_recs_2010$numSchinus <- florida_df[[14]]$numSchinus - florida_df[[13]]$numSchinus

new_recs_2020 <- florida_df[[15]]
new_recs_2020$numAnacard <- florida_df[[15]]$numAnacard - florida_df[[14]]$numAnacard
new_recs_2020$numSchinus <- florida_df[[15]]$numSchinus - florida_df[[14]]$numSchinus 


ld_Shape_new1970 = shapeLopodData( Shapefile = new_recs_1970,
                                fieldN = "numAnacard",
                                fieldY = "numSchinus",  
                                Adjacency = T,
                                keepFields = F)

ld_Shape_new1980 = shapeLopodData( Shapefile = new_recs_1980,
                                   fieldN = "numAnacard",
                                   fieldY = "numSchinus",  
                                   Adjacency = T,
                                   keepFields = F)

ld_Shape_new1990 = shapeLopodData( Shapefile = new_recs_1990, #something weird more detections than effort
                                   fieldN = "numAnacard",
                                   fieldY = "numSchinus",  
                                   Adjacency = T,
                                   keepFields = F)
ld_Shape_new2000 = shapeLopodData( Shapefile = new_recs_2000,
                                   fieldN = "numAnacard",
                                   fieldY = "numSchinus",  
                                   Adjacency = T,
                                   keepFields = F)
ld_Shape_new2010 = shapeLopodData( Shapefile = new_recs_2010,
                                   fieldN = "numAnacard",
                                   fieldY = "numSchinus",  
                                   Adjacency = T,
                                   keepFields = F)

ld_Shape_new2020 = shapeLopodData( Shapefile = new_recs_2020,
                                   fieldN = "numAnacard",
                                   fieldY = "numSchinus",  
                                   Adjacency = T,
                                   keepFields = F)




  schinusSEff_1970 = spplot(new_recs_1970, zcol = c("numAnacard"), main = "Sampling Effort-1970", border=NA)
  schinusDetect_1970 = spplot(new_recs_1970, zcol = c("numSchinus"), main = "Detections-1970", border=NA)
  plot(schinusSEff_1970, split = c(1,1,2,1), xlab="1970", more = T, border=NA)
  plot(schinusDetect_1970,  split = c(2,1,2,1), more = F, border=NA)
  
  schinusSEff_1980 = spplot(new_recs_1980, zcol = c("numAnacard"), main = "Sampling Effort-1980", border=NA)
  schinusDetect_1980 = spplot(new_recs_1980, zcol = c("numSchinus"), main = "Detections-1980", border=NA)
  plot(schinusSEff_1980, split = c(1,1,2,1), xlab="1980", more = T, border=NA)
  plot(schinusDetect_1980,  split = c(2,1,2,1), more = F, border=NA)

  schinusSEff_2000 = spplot(new_recs_2000, zcol = c("numAnacard"), main = "Sampling Effort-2000", border=NA)
  schinusDetect_2000 = spplot(new_recs_2000, zcol = c("numSchinus"), main = "Detections-2000", border=NA)
  plot(schinusSEff_2000, split = c(1,1,2,1), xlab="2000", more = T, border=NA)
  plot(schinusDetect_2000,  split = c(2,1,2,1), more = F, border=NA)
  
  schinusSEff_2010 = spplot(new_recs_2010, zcol = c("numAnacard"), main = "Sampling Effort-2010", border=NA)
  schinusDetect_2010 = spplot(new_recs_2010, zcol = c("numSchinus"), main = "Detections-2010", border=NA)
  plot(schinusSEff_2010, split = c(1,1,2,1), xlab="2010", more = T, border=NA)
  plot(schinusDetect_2010,  split = c(2,1,2,1), more = F, border=NA)
  
 # schinusSEff_2020 = spplot(new_recs_2020, zcol = c("numAnacard"), main = "Sampling Effort-2020", border=NA)
  #schinusDetect_2020 = spplot(new_recs_2020, zcol = c("numSchinus"), main = "Detections-2020", border=NA)
  #plot(schinusSEff_2020, split = c(1,1,2,1), xlab="2020", more = T, border=NA)
  #plot(schinusDetect_2020,  split = c(2,1,2,1), more = F, border=NA)

  
  mLopodShape_new1980 = modelLopod(LopodData = ld_Shape_new1980,
                                varP = T,
                                q = NULL,
                                pmin = 0,
                                CAR = F,
                                nChains = 3,
                                warmup = 500,
                                sampling = 1000,
                                nCores =3)
  lopodTrace( mLopodShape_new1980, inc_warmup = T)
  
  ###Look at detections
  SchinusShape_new1980 = mLopodShape_new1980@LopodData@geoDataObject
  psiShape_new1980 = lopodShape(mLopodShape_new1980, "psi_Sampled", extrapolate = F,  quant = 0.05)
  SchinusShape_new1980@data[,"psi_Sampled_05"] = psiShape_new1980@data[,"psi_Sampled"]
  
  psiShape_new1980 = lopodShape(mLopodShape_new1980, "psi_Sampled", extrapolate = F,  quant = 0.5)
  SchinusShape_new1980@data[,"psi_Sampled_50"] = psiShape_new1980@data[,"psi_Sampled"]
  
  psiShape_new1980 = lopodShape(mLopodShape_new1980, "psi_Sampled", extrapolate = F,  quant = 0.95)
  SchinusShape_new1980@data[,"psi_Sampled_95"] = psiShape_new1980@data[,"psi_Sampled"]
  spplot( SchinusShape_new1980,
          zcol = c("psi_Sampled_05", "psi_Sampled_50", "psi_Sampled_95"),
          names.attr = c("Psi (5% quantile)", "Psi (median)", "Psi (95% quantile)"), main = "Psi:1980")
  
  
  mLopodShape_new2000 = modelLopod(LopodData = ld_Shape_new2000,
                                   varP = T,
                                   q = NULL,
                                   pmin = 0,
                                   CAR = F,
                                   nChains = 3,
                                   warmup = 2000,
                                   sampling = 20000,
                                   nCores =3)
  lopodTrace( mLopodShape_new2000, inc_warmup = T)
  save(mLopodShape_new2000, file="mLopodShape_new2000")
  
  ###Look at detections
  SchinusShape_new1980 = mLopodShape_new1980@LopodData@geoDataObject
  psiShape_new1980 = lopodShape(mLopodShape_new1980, "psi_Sampled", extrapolate = F,  quant = 0.05)
  SchinusShape_new1980@data[,"psi_Sampled_05"] = psiShape_new1980@data[,"psi_Sampled"]
  
  psiShape_new1980 = lopodShape(mLopodShape_new1980, "psi_Sampled", extrapolate = F,  quant = 0.5)
  SchinusShape_new1980@data[,"psi_Sampled_50"] = psiShape_new1980@data[,"psi_Sampled"]
  
  psiShape_new1980 = lopodShape(mLopodShape_new1980, "psi_Sampled", extrapolate = F,  quant = 0.95)
  SchinusShape_new1980@data[,"psi_Sampled_95"] = psiShape_new1980@data[,"psi_Sampled"]
  spplot( SchinusShape_new1980,
          zcol = c("psi_Sampled_05", "psi_Sampled_50", "psi_Sampled_95"),
          names.attr = c("Psi (5% quantile)", "Psi (median)", "Psi (95% quantile)"), main = "Psi:1980")
  
  
  
  mLopodShape_new2010 = modelLopod(LopodData = ld_Shape_new2010,
                                   varP = T,
                                   q = NULL,
                                   pmin = 0,
                                   CAR = F,
                                   nChains = 3,
                                   warmup = 1000,
                                   sampling = 2000,
                                   nCores =3)
  lopodTrace( mLopodShape_new2010, inc_warmup = T)
  
  SchinusShape_new2010 = mLopodShape_new2010@LopodData@geoDataObject
  psiShape_new2010 = lopodShape(mLopodShape_new2010, "psi_Sampled", extrapolate = F,  quant = 0.05)
  SchinusShape_new2010@data[,"psi_Sampled_05"] = psiShape_new2010@data[,"psi_Sampled"]
  
  psiShape_new2010 = lopodShape(mLopodShape_new2010, "psi_Sampled", extrapolate = F,  quant = 0.5)
  SchinusShape_new2010@data[,"psi_Sampled_50"] = psiShape_new2010@data[,"psi_Sampled"]
  
  psiShape_new2010 = lopodShape(mLopodShape_new2010, "psi_Sampled", extrapolate = F,  quant = 0.95)
  SchinusShape_new2010@data[,"psi_Sampled_95"] = psiShape_new2010@data[,"psi_Sampled"]
  spplot( SchinusShape_new2010,
          zcol = c("psi_Sampled_05", "psi_Sampled_50", "psi_Sampled_95"),
          names.attr = c("Psi (5% quantile)", "Psi (median)", "Psi (95% quantile)"), main = "Psi:2010")
  
  
  pShape_new2010 = lopodShape(mLopodShape_new2010, "pCorr", extrapolate = F,  quant = 0.05)
  SchinusShape_new2010p@data[,"pCorr_05"] = pShape_new2010@data[,"pCorr"]
  
  pShape_new2010 = lopodShape(mLopodShape_new2010, "pCorr", extrapolate = F,  quant = 0.5)
  SchinusShape_new2010p@data[,"pCorr_50"] = pShape_new2010@data[,"pCorr"]
  
  pShape_new2010 = lopodShape(mLopodShape_new2010, "pCorr", extrapolate = F,  quant = 0.95)
  SchinusShape_new2010p@data[,"pCorr_95"] = pShape_new2010@data[,"pCorr"]
  spplot( SchinusShape_new2010p,
          zcol = c("pCorr_05", "pCorr_50", "pCorr_95"),
          names.attr = c("p (5% quantile)", "p (median)", "p (95% quantile)"), main = "p:2010")
  
  