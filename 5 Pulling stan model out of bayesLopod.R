# psiipq: Global p and q, psi for each sampling unit
#   if(LopodData@geoType == "Shapefile"){
#     stanData = list( nSampledCells = length(LopodData@geoInfo$sampledId$featureShape),
#                      sampledId = LopodData@geoInfo$sampledId$cellStan,
#                      N = N,
#                      y = y,
#                      minP = pmin
#     )
#
#   }
#
# psiipiq: Global q, psi and p for each sampling unit
# if(LopodData@geoType == "Shapefile"){
#
#   stanData = list( nSampledCells = length(LopodData@geoInfo$sampledId$featureShape),
#                    sampledId = LopodData@geoInfo$sampledId$cellStan,
#                    N = N,
#                    y = y,
#                    minP = pmin
#   )
#
# }

pmin=0
getParamsN <- function(LopodObject) {
  N = LopodObject@geoDataObject$sampEffort[LopodObject@geoInfo$sampledId$featureShape]
  return(N)
}

getParamsy <- function(LopodObject) {
  y = LopodObject@geoDataObject$detections[LopodObject@geoInfo$sampledId$featureShape]
  return(y)
}


N <- getParamsN(ld_Shape_1970)
y <- getParamsy(ld_Shape_1970)


stanData = list( nSampledCells = length(ld_Shape_1970@geoInfo$sampledId$featureShape),
                                       sampledId = ld_Shape_1970@geoInfo$sampledId$cellStan,
                                       N = N,
                                       y = y,
                                       minP = pmin
                      )
fit <- stan(file = 'psiipq.stan', data = stanData,
            iter = 1000, chains = 3)

traceplot(fit, inc_warmup=T)

stan_trace(fit, pars="lp__", inc_warmup=T)


finalShapefile = ld_Shape_1970@geoDataObject
FeatureName="psi"
finalShapefile@data[,FeatureName] = NA
if((param %in% c(modelPar$allCellsPars,modelPar$sampledPars))==F) stop(paste("For this model only the following  parameters can be mapped:",toString(c(modelPar$allCellsPars,modelPar$sampledPars))))

CellsID = ld_Shape_1970@geoInfo$sampledId
ParObjects=paste(fit@model_pars,"[",CellsID$cellStan,"]",sep="")


if ((is.null(LopodModel@LopodData@geoInfo$W_sparse)==F)&(LopodModel@modelInfo$CAR==F)){

  ParObjects=paste(fit@model_pars,"[",1:dim(CellsID)[1],"]",sep="")

}


ParValues = rstan::summary(fit)$summary[,columnName]
#grab the ones related to psi_Sampled
myPars2 <- rstan::summary(fit)$summary[1:48, ]
myPars_mean <- myPars2[, 1]

finalShapefile@data[CellsID$featureShape,FeatureName] = myPars_mean


}

return(finalShapefile)

