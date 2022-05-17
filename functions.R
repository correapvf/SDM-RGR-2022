# Function to split ROV dives in sequential folds for cross validation
# dive_cv <- function(dives, number = 5, repeats = 10, numLimit = 15) {
#     message("This function assumes that data is ordered by Dive and Time")
# 
#     tmp <- split(seq_along(dives), dives)
#     tmp2 <- lapply(tmp, function(x) split(x, ceiling(seq_along(x) / (length(x) / number))))
# 
#     out <- vector("list", number * repeats)
#     i <- 1
#     for (r in 1:repeats) {
#         tmp3 <- lapply(tmp2, function(x) x[sample(1:number)])
# 
#         for (n in 1:number) {
#             x <- lapply(tmp3, `[[`, n)
#             out[[i]] <- unname(unlist(mapply(function(x, y) x[!(x %in% y)], tmp, x)))
#             i <- i + 1
#         }
#     }
# 
#     names(out) <- paste0("Rep", rep(1:repeats, each = number), "_CV", rep(1:number, repeats))
# 
#     return(out)
# }




dive_cv <- function(dives, number = 25, p = 0.8) {
  message("This function assumes that data is ordered by Dive and Time")
  
  tmp <- split(seq_along(dives), dives)
  ndives <- length(tmp)
  
  tails <- lapply(tmp, function(x) split(x, ceiling(seq_along(x) / (length(x) * p))))
  heads <- lapply(tmp, function(x) split(x, rev(ceiling(seq_along(x) / (length(x) * p)))))
  ht <- list(lapply(heads, `[[`, 1), lapply(tails, `[[`, 1))
  
  out <- vector("list", number)
  for (n in 1:number) {
    ht_sel <- sample(1:2, ndives, replace = TRUE)
    out[[n]] <- mapply(function(x, y) ht[[y]][[x]], 1:ndives, ht_sel)
  }

  names(out) <- paste0("Rep", 1:number)
  out <- lapply(out, unlist)
  return(out)
}



# for ploting
guide_externalticks <- function(...) {
  guide <- guide_colorbar(...)
  class(guide) <- c("guide", "guide_externalticks", "colorbar")
  guide
}

guide_gengrob.guide_externalticks <- function(guide, theme) {
  # browser()
  dir <- guide$direction
  guide <- NextMethod()
  is_ticks <- grep("^ticks$", guide$layout$name)
  ticks <- guide$grobs[is_ticks][[1]]
  if (dir == "vertical") {
    
    ticks$x0 <- rep(tail(ticks$x1, 1), length(ticks$x0)/2)
    ticks$x1 <- ticks$x0 + head(ticks$x1, 1)
  } else {
    ticks$y1 <- rep(tail(ticks$y1, 1), length(ticks$y1))
  }
  
  guide$grobs[[is_ticks]] <- ticks
  guide
}


define_breaks <- function(df, n) {
  deg2dm <- function(x, lat) {
    d1 <- abs(x)
    d2 <- floor(d1)
    l <- if (lat) ifelse(x > 0, "N", "S") else ifelse(x > 0, "E", "W")
    sprintf("%d°%04.1f'%s", d2, (d1-d2)*60, l)
  }
  
  minmax <- range(df$x)
  xutm <- seq(from = minmax[1], to = minmax[2], length.out = n*2+1)[(1:n)*2]
  
  minmax <- range(df$y)
  yutm <- seq(from = minmax[1], to = minmax[2], length.out = n*2+1)[(1:n)*2]
  
  sputm <- SpatialPoints(data.frame(x=xutm, y=yutm),
                         proj4string = CRS("+proj=utm +zone=25 +south +datum=WGS84 +units=m +no_defs"))
  spgeo <- spTransform(sputm, CRS("+proj=longlat +datum=WGS84"))
  
  xdeg <- deg2dm(spgeo@coords[, 1], FALSE)
  ydeg <- deg2dm(spgeo@coords[, 2], TRUE)
  
  return(list(breakx = xutm, breaky = yutm, labelx = xdeg, labely = ydeg))
}


plot_clamp <- function(x, probs = c(0.025, 0.975)) {
  q <- quantile(x, probs = probs)
  x[x < q[1]] <- q[1]
  x[x > q[2]] <- q[2]
  return(x)
}


theme_mapa <- function() {
  theme_gray(base_size = 8, base_line_size = 0.25) %+replace%    #replace elements we want to change
    theme(
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_rect(fill = NA, colour = "black"),
      axis.title = element_blank(),
      axis.text.y = element_text(angle = 90, size = 5),
      axis.text.x = element_text(size = 5),
      strip.text = element_blank(),
      strip.background = element_blank()
    )
}


theme_artigo <- function() {
  theme_classic(base_size = 8, base_line_size = 0.25) %+replace%    #replace elements we want to change
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = rel(0.9)),
      plot.tag = element_text(size = 9.5, face = "bold", hjust = -0.1),
      plot.tag.position = c(0, 1)
    )
}



library(splines)
calibration.plot <- function(model, newdata = NULL, bins = 10, smoothingdf = 6) {
  pred <- predict2(model, newdata, type = "prob1")
  obs <- if (is.null(newdata)) model$trainingData$.outcome else caretSDM:::get.response(model, newdata, NULL)
  obs2 <- abs(as.numeric(obs) - 2)
  
  
  # based on functions provided by XXXXX
  # smoothdist functions
  require(splines)
  gam1 <- glm(obs2 ~ ns(pred, df=smoothingdf), weights=rep(1, length(pred)), family=binomial)
  x <- seq(min(pred), max(pred), length = 512)
  y <- predict(gam1, newdata = data.frame(pred = x), se.fit = TRUE,
               type = "response")
  predd <- data.frame(x=x, y=y$fit, se=y$se.fit)
  
  return(list(predd=predd, obs=data.frame(obs = obs, pred = pred)))
  
  # # calibplot
  # ylow <- predd$y - 2 * predd$se
  # ylow[ylow<0] <- 0
  # yhigh <- predd$y + 2 * predd$se
  # yhigh[yhigh>1] <- 1
  # plot(predd$x, ylow, type="l", col="orange", ylim=c(0,1), xlim=c(0,1), lwd=2)
  # lines(predd$x, yhigh, col="orange")
  # abline(0, 1, lty="dashed")
  # lines(predd$x, predd$y, col="deepskyblue")
  # rug(pred[obs==0])
  # rug(pred[obs==1], col = "orange")
  
  # ecalp
  # g <- floor(pred*bins)
  # b <- 0:(bins-1)
  # p <- sapply(b, function(x) if (length(obs[g==x])==0) -1 else sum(obs[g==x]) / length(obs[g==x]))
  # mx <- sapply(b, function(x,g) mean(pred[g==x]), g)
  # points(mx, p, xlim=c(0,1), ylim=c(0,1))
  
}

twoClassSDMs <- function(data, lev = NULL, model = NULL) {
  # modifield from caret to include TSS metric
  
  if (length(lev) > 2) {
    stop(paste("Your outcome has", length(lev), "levels. The twoClassSummary() function isn't appropriate."))
  }
  if (!all(levels(data[, "pred"]) == lev)) {
    stop("levels of observed and predicted data do not match")
  }
  
  dataComplete <- data[stats::complete.cases(data), ]
  kapp <- unlist(e1071::classAgreement(table(dataComplete[, 
                                                          "obs"], dataComplete[, "pred"])))["kappa"]
  
  rocObject <- try(pROC::roc(data$obs, data[, lev[1]], direction = ">",
                             quiet = TRUE), silent = TRUE)
  rocAUC <- if (inherits(rocObject, "try-error")) NA else rocObject$auc
  
  Se <- caret::sensitivity(data[, "pred"], data[, "obs"], lev[1])
  Sp <- caret::specificity(data[, "pred"], data[, "obs"], lev[2])
  
  y <- ifelse(data$obs == lev[1], 1, 0)
  p <- data[, lev[1]]
  
  # install_github('meeliskull/prg/R_package/prg')
  prg <- prg::calc_auprg(
    prg::create_prg_curve(labels = ifelse(data$obs == lev[1], 1, 0), 
                          pos_scores = data[, lev[1]])
  )
  
  out <- c(rocAUC, prg, Se, Sp, Se + Sp - 1)
  names(out) <- c("ROC", "PRG", "Sens", "Spec", "TSS")
  return(out)
}
