library(raster)
library(data.table)
library(caret)
library(corrplot)
library(caretSDM)
library(usdm)

# setup

spps <- "Scep1" # "Loph1"
gdrive <- "G:/Meu Drive"

folder_mx <- "D:/Users/paulo/maxent_temp" # "~/Documents/maxent_temp"
folder_mx_final <- "maxent_model"
maxentPath <- "D:/Users/paulo/maxent.jar" # "~/maxent.jar

##---- load data -------------------
source("load_data.R")
source("functions.R")

# load("workspace.RData")
# save.image("workspace.RData")

# set parallel workers
# library(doParallel)
# cl <- makePSOCKcluster(16)
# registerDoParallel(cl)

library(doSNOW)
cl <- makeSOCKcluster(16)
registerDoSNOW(cl)

# SET species AND enviroment DATA FOR MODELS
envs <- names(r.hybis)
spp <- "Scep1"



# check for corralations and VIFs
M <- vifcor(hybis[, ..envs])
M


# create model control
form.spp <- as.formula(paste(spp, "~", paste(envs, collapse = "+")))

set.seed(12345)
index <- dive_cv(hybis$Dive, number = 25, p = 0.8)
# check if all cvs have both presence and ausence
# check <- lapply(index, function(x, y) table(y[x]), y = hybis[[spp]])
# all(sapply(check, function(x) x[1] > 15))

control <- trainControl(method = "spatialcv",
                          index = index,
                          classProbs = TRUE,
                          savePredictions = "final",
                          seeds = NULL,
                          summaryFunction = twoClassSDM)
                          #verboseIter = TRUE) # won't work in parallel

metrictune <- "ROC"
metricthr <- "Sens=Spec"

# shinkai <- shinkai[Dive == "6K1338"] # filter only first Dive
r.hybis.clamp <- clamp_data(as.data.frame(hybis), r.hybis)

retrain_formula <- function(profile) {
  tmpenv <- profile$optVariables
  tmpenv <- gsub("north_south1", "north_south", tmpenv)
  as.formula(paste(spp, "~", paste(tmpenv, collapse = "+")))
}

subsets <- 2:9
caretFuncs$summary <- twoClassSummary
rfeCtrl <- rfeControl(functions = caretFuncs,
                      index = control$index)


# Random Forest
set.seed(12345)
model.rf <- train(form.spp,
                data = hybis,
                method = "rf",
                metric = metrictune,
                trControl = control,
                tuneGrid = expand.grid(mtry = seq(1, 9, 1)),
                ntree = 1000)
model.rf


rfProfile <- rfe(form.spp,
                  data = hybis,
                  sizes = subsets,
                  rfeControl = rfeCtrl,
                  metric = metrictune,
                  method = "rf",
                  trControl = trainControl(method = "none", classProb = TRUE),
                  tuneGrid = model.rf$bestTune
)
rfProfile


thr.rf <- thresholder2(model.rf, final = TRUE, thr.interval = 0.001)
model.rf <- setThreshold(model.rf, summary(thr.rf, which.method = metricthr))

evaluate(model.rf)
# 
# r.rf <- response(model.rf)
# plot(r.rf, plot_thr = FALSE)
# 
p.rf <- predict2(r.hybis.clamp, model.rf, type = "both1")
plot(p.rf)






# Boosted Regression Tree
set.seed(12345)
model.brt <- train(form.spp,
                  data = hybis,
                  method = "gbm",
                  metric = metrictune,
                  trControl = control,
                  tuneGrid = expand.grid(n.trees = seq(100, 3000, 50), # seq(100, 3000, 50),
                                         interaction.depth = seq(1, 16, 3), # seq(1, 19, 3),
                                         shrinkage = 10^seq(-1, -4), # 10^seq(-1, -4),
                                         n.minobsinnode = c(5, 15)), # c(5, 10, 15)),
                  verbose = FALSE)

model.brt


gbmModel <- getModelInfo("gbm")$gbm
gbmModel$varImp <- function(object, numTrees = NULL, ...) {
  if (is.null(numTrees)) 
    numTrees <- object$tuneValue$n.trees
  varImp <- gbm::relative.influence(object, n.trees = numTrees)
  out <- data.frame(varImp)
  colnames(out) <- "Overall"
  rownames(out) <- object$var.names
  out
}
brtProfile <- rfe(form.spp,
                  data = hybis,
                  sizes = subsets,
                  rfeControl = rfeCtrl,
                  metric = metrictune,
                  method = gbmModel,
                  trControl = trainControl(method = "none", classProb = TRUE),
                  tuneGrid = model.brt$bestTune,
                  verbose = FALSE
)
brtProfile



thr.brt <- thresholder2(model.brt, final = TRUE, thr.interval = 0.001)
model.brt <- setThreshold(model.brt, summary(thr.brt, which.method = metricthr))

evaluate(model.brt)
# 
# r.brt <- response(model.brt)
# plot(r.brt, plot_thr = FALSE)
# 
p.brt <- predict2(r.hybis.clamp, model.brt, type = "both1")
plot(p.brt)


# # Maxent setup - for background data as absence data - not used in the final paper 
# nbg <- 10000
# 
# # select background points
# dens <- MASS::kde2d(hybis$Lon, hybis$Lat, n = c(ncol(r.hybis), nrow(r.hybis)),
#                     lims = as.vector(extent(r.hybis)))
# dens <- raster(apply(dens$z, 1, rev), template = r.hybis)
# 
# dens[sum(is.na(r.hybis)) > 0] <- 0
# 
# # a implementacao do R e muito lenta
# set.seed(12345)
# bg <- wrswoR::sample_int_crank(n = ncell(dens), size = nbg, prob = as.vector(dens))
# 
# bg.xy <- xyFromCell(r.hybis, bg)
# bg.val <- extract(r.hybis[[envs]], bg)
# bg.y <- rep("ausence", nbg)
# 
# # update index
# tmp <- which(hybis[[spp]] == "presence")
# index.bg <- lapply(index, intersect, y = tmp)
# 
# n_hybis <- nrow(hybis)
# indexOut.bg <- lapply(index, setdiff, x = seq_len(n_hybis))
# indexOut.bg <- lapply(indexOut.bg, intersect, y = tmp)
# 
# index.bg2 <- createMultiFolds(rep("ausence", nbg), k = 5, times = 2) # 10
# indexOut.bg2 <- lapply(index.bg2, function(x, y) setdiff(y, x), y = seq_len(nbg))
# 
# index.bg2 <- lapply(index.bg2, `+`, n_hybis)
# indexOut.bg2 <- lapply(indexOut.bg2, `+`, n_hybis)
# 
# index.bg <- mapply(c, index.bg, index.bg2)
# indexOut.bg <- mapply(c, indexOut.bg, indexOut.bg2)
# 
# hybis.bg <- rbind(hybis[, c(spp, "Lon", "Lat", envs), with = FALSE],
#                      data.table(bg.y, bg.xy, bg.val), use.names = FALSE)
# 
# control.bg <- control
# control.bg$index <- index.bg
# control.bg$indexOut <- indexOut.bg
# rm(dens, bg, bg.xy, bg.val, bg.y, tmp, index.bg, indexOut.bg, index.bg2, indexOut.bg2)


# Maxent
set.seed(12345)

dir.create(folder_mx_final, showWarnings = FALSE)
dir.create(folder_mx, showWarnings = FALSE)
unlink(list.files(folder_mx_final, full.names = TRUE), recursive = TRUE)

model.maxent <- train(form.spp,
                      data = hybis,
                      method = maxentCaret,
                      metric = metrictune,
                      trControl = control,
                      filesPath = folder_mx,
                      maxentPath = maxentPath,
                      thrtype = metricthr,
                      # flags = 'addsamplestobackground=false',
                      tuneGrid = expand.grid(reg = c("lqpth", "lqt", "lqpt", "lqh", "lqph", "lqp", "lq", "l"),
                                             beta = c(seq(0.2, 2, 0.2), 2.5, 3, 3.5, 4, 4.5, 5))
                    )

# save final Maxent model
file.copy(model.maxent$finalModel$path, folder_mx_final, recursive = TRUE)
unlink(list.files(folder_mx, full.names = TRUE), recursive = TRUE)
model.maxent$finalModel$path <- file.path(folder_mx_final, basename(model.maxent$finalModel$path))
model.maxent$finalModel$lambda <- file.path(model.maxent$finalModel$path, basename(model.maxent$finalModel$lambda))

model.maxent



maxentProfile <- rfe(form.spp,
                  data = hybis,
                  sizes = subsets,
                  rfeControl = rfeCtrl,
                  metric = metrictune,
                  method = maxentCaret,
                  filesPath = folder_mx,
                  maxentPath = maxentPath,
                  trControl = trainControl(method = "none", classProb = TRUE),
                  tuneGrid = model.maxent$bestTune
)
maxentProfile


thr.maxent <- thresholder2(model.maxent, final = TRUE, thr.interval = 0.001)
model.maxent <- setThreshold(model.maxent, summary(thr.maxent, which.method = metricthr))

evaluate(model.maxent)
# 
# r.maxent <- response(model.maxent)
# plot(r.maxent, plot_thr = FALSE)
# 
p.maxent <- predict2(r.hybis.clamp, model.maxent, type = "both1")
plot(p.maxent)


# Artificial Neural Networks
set.seed(12345)
model.ann <- train(form.spp,
                    data = hybis,
                    method = "nnet",
                    metric = metrictune,
                    trControl = control,
                    preProcess = c("center", "scale"),
                    tuneGrid = expand.grid(size = c(1, 3, 5, 7, 9, 11, 16, 21), # c(1, 3, 5, 7, 9, 11, 13, 21, 31),
                                           decay = c(0, 10^seq(-1, -5))), # c(0, 10^seq(-1, -5))),
                    maxit = 500
)
model.ann



annProfile <- rfe(form.spp,
                  data = hybis,
                  sizes = subsets,
                  rfeControl = rfeCtrl,
                  metric = metrictune,
                  method = "nnet",
                  preProcess = c("center", "scale"),
                  trControl = trainControl(method = "none", classProb = TRUE),
                  tuneGrid = model.ann$bestTune,
                  maxit = 500,
                  trace = FALSE
)

annProfile


thr.ann <- thresholder2(model.ann, final = TRUE, thr.interval = 0.001)
model.ann <- setThreshold(model.ann, summary(thr.ann, which.method = metricthr))


evaluate(model.ann)
# 
# r.ann <- response(model.ann)
# plot(r.ann, plot_thr = FALSE)
# 
p.ann <- predict2(r.hybis.clamp, model.ann, type = "both1")
plot(p.ann)


# Generalized Addictive Modes - to test a custom formula
# gamcustom <- getModelInfo("gam")$gam
# 
# gamcustom$fit <- function(x, y, wts, param, lev, last, classProbs, ...) {
#                     require(mgcv)
# 
#                     # create formula
#                     modForm <- y ~ s(cos_aspect, k = 10) + s(depth, k = 10) + s(finebpi_3_30, k = 10) +
#                     s(rift_dist, m = 1, k = 6) + s(sin_aspect, k = 10) + s(slope, k = 10) + s(vrm_011, k = 10)
# 
#                     # adjust data set
#                     dat <- if (is.data.frame(x)) x else as.data.frame(x, stringsAsFactors = TRUE)
#                     dat <- cbind(dat, y = y)
# 
#                     out <- mgcv::gam(modForm, data = dat, family = binomial(),
#                                      select = param$select, method = as.character(param$method), ...)
#                     out
# 
#                     }

set.seed(12345)
model.gam <- train(form.spp,
                    data = hybis,
                    method = "gam",
                    metric = metrictune,
                    trControl = control,
                    tuneGrid = expand.grid(method = c("REML", "ML"),
                                           select = c(TRUE, FALSE))
)
# model.gam$method <- "gam"
model.gam


# gamProfile is not needed, because GAM has a "select" parameter


# check select variables by gam
summary(model.gam$finalModel)

# gam.check(model.gam$finalModel)

thr.gam <- thresholder2(model.gam, final = TRUE, thr.interval = 0.001)
model.gam <- setThreshold(model.gam, summary(thr.gam, which.method = metricthr))

evaluate(model.gam)
# 
# r.gam <- response(model.gam)
# plot(r.gam, plot_thr = FALSE)
# 
p.gam <- predict2(r.hybis.clamp, model.gam, type = "both1")
plot(p.gam)




# Ensemble model
models <- list(model.rf, model.brt, model.maxent, model.ann, model.gam)

# # evaluate models first for ensemble
e.models <- evaluate(models, errorFunction = se)
dot_plot(e.models, metric = "ROC")
# pairs_plot(e.models, metric = "ROC")
# dot_plot(e.models, metric = "TSS")
# pairs_plot(e.models, metric = "TSS")
# 
r.models <- response(models, errorFunction = se)
# plot(r.models)
# ?plot.response.train


# Do the ensemble
ens <- createEnsemble(models, calc.pred = TRUE, ensemble_method = "mean")


thr.ens <- thresholder2(ens, final = TRUE, thr.interval = 0.001)
ens <- setThreshold(ens, summary(thr.ens, which.method = metricthr))


evaluate(ens)
# 
# r.ens <- response(ens)
# plot(r.ens, plot_thr = FALSE)
# 
p.ens <- predict2(r.hybis.clamp, ens, type = "both1")
plot(p.ens)


# Variables importance
set.seed(12345)
v.models <- varImp2(models, nperm = 100, errorFunction = se)
# plot(v.models)

# save final  models
names(ens$model.list) <- c('rf', 'brt', 'maxent', 'ann', 'gam')
save(M, envs, spp, ens, hybis, shinkai, r.hybis, r.shinkai,
     metrictune, control, metricthr, form.spp, e.models, r.models, v.models,
     annProfile, brtProfile, maxentProfile, rfProfile,
     thr.ann, thr.brt, thr.gam, thr.maxent, thr.rf, thr.ens, file = "modelos.RData")


# confidence maps
set.seed(12345)
conf.hybis <- confidence_map(ens, r.hybis, nrep = 200, doclamp = TRUE, return.all = FALSE)
# plot(conf.hybis)
writeRaster(conf.hybis, "rasters/hybis_ensemble_confidence.tif")

set.seed(12345)
conf.shinkai <- confidence_map(ens, r.shinkai, nrep = 200, doclamp = TRUE, return.all = FALSE)
# plot(conf.shinkai)
writeRaster(conf.shinkai, "rasters/shinkai_ensemble_confidence.tif")


# model predictions and save
p.models.hybis <- predict2(r.hybis, c(models, list(ens)), doclamp = TRUE, type = "both1")
names(p.models.hybis)[c(5:6, 11:12)] <- c("maxent_pred", "maxent_presence", "ensemble_pred", "ensemble_presence")
writeRaster(p.models.hybis, "rasters/hybis.tif", bylayer = TRUE, suffix = "names")

p.models.shinkai <- predict2(r.shinkai, c(models, list(ens)), doclamp = TRUE, type = "both1")
names(p.models.shinkai)[c(5:6, 11:12)] <- c("maxent_pred", "maxent_presence", "ensemble_pred", "ensemble_presence")
writeRaster(p.models.shinkai, "rasters/shinkai.tif", bylayer = TRUE, suffix = "names")


stopCluster(cl)
registerDoSEQ()
