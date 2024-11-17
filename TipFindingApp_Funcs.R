#-------------------------------------------------------------------------------
# FUNCTION DEFINITIONS----------------------------------------------------------
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
### MAIN FUNCTIONS--------------------------------------------------------------
#-------------------------------------------------------------------------------
TreatKymo <- function(kymo){
  kymo <- apply(kymo, 2, as.numeric)
  kymo.out <- kymo
  kymo.fault.indx <- c()
  for(i in 1:dim(kymo)[1]){
    t.slice <- kymo[i, ]
    bad.entries <- which(is.na(t.slice) | is.nan(t.slice) | is.infinite(t.slice))
    if((length(bad.entries) > 0) & (length(bad.entries) < length(t.slice))){
      t.slice <- SubstituteWithNeighbor(t.slice, subs.indx = bad.entries)
      kymo.out[i, ] <- t.slice 
    }
    if(length(bad.entries) == length(t.slice)){
      kymo.fault.indx <- c(kymo.fault.indx, i)
    } else{
      if(all(diff(diff(t.slice)) == 0)){
        kymo.fault.indx <- c(kymo.fault.indx, i)  
      }
    }
  }
  if(length(kymo.fault.indx) > 0){
    kymo.out[kymo.fault.indx, ] <- apply(kymo.out, 2, SubstituteWithNeighbor, kymo.fault.indx)
  }
  
  return(kymo.out)  
}

SubstituteWithNeighbor <- function(vec, subs.indx){
  return(approx(x = (1:length(vec))[-subs.indx], y = vec[-subs.indx], xout = subs.indx)$y)
}

#------------------------------------------------------------------------------- 
TreatKymos <- function(k1, k2){
  k1.out <- k1 <- apply(k1, 2, as.numeric)
  k2.out <- k2  <- apply(k2, 2, as.numeric)
  
  kymo.fault.indx <- c()
  for(i in 1:dim(k1)[1]){
    t.slice.1 <- k1[i, ]
    t.slice.2 <- k2[i, ]
    
    bad.entries.1 <- which(is.na(t.slice.1) | is.nan(t.slice.1) | is.infinite(t.slice.1))
    bad.entries.2 <- which(is.na(t.slice.2) | is.nan(t.slice.2) | is.infinite(t.slice.2))
    bad.entries <- unique(bad.entries.1, bad.entries.2)
    if((length(bad.entries) > 0) & (length(bad.entries) < length(t.slice))){
      t.slice.1 <- SubstituteWithNeighbor(t.slice.1, subs.indx = bad.entries)
      t.slice.2 <- SubstituteWithNeighbor(t.slice.2, subs.indx = bad.entries)
      k1.out[i, ] <- t.slice.1
      k2.out[i, ] <- t.slice.2
    }
    if(length(bad.entries) == length(t.slice)){
      kymo.fault.indx <- c(kymo.fault.indx, i)
    } else{
      if(all(diff(diff(t.slice.1)) == 0) | all(diff(diff(t.slice.2)) == 0)){
        kymo.fault.indx <- c(kymo.fault.indx, i)  
      }
    }
  }
  if(length(kymo.fault.indx) > 0){
    k1.out[kymo.fault.indx, ] <- apply(k1.out, 2, SubstituteWithNeighbor, kymo.fault.indx)
    k2.out[kymo.fault.indx, ] <- apply(k2.out, 2, SubstituteWithNeighbor, kymo.fault.indx)
  }
  
  return(list(k1.out, k2.out))  
}

#--------------------------------------------------------------------------------

# FindTip <- function(kymo, use.smooth = TRUE, kymo.span.n = 7, kymo.loess.dg = 1, 
#                     n.pts = 5, mad.tol = 3.5, qntl = 0.95, fit.right = FALSE, fix.slope = FALSE){
#   kymo <- as.matrix(apply(kymo, 2, as.numeric))
#   
#   # Smooth image
#   if(use.smooth){
#     imaj <- t(apply(kymo, 1, SmoothWithLoess, #INTEGRATE VERSION?
#                     spn = kymo.span.n / (dim(kymo)[2]),
#                     degree = kymo.loess.dg))
#   } else{imaj <- kymo}
#   
#   # Subtract minimium value
#   imaj <- imaj - min(imaj, na.rm = TRUE)
#   
#   tip.tbl <- t(as.data.frame(apply(imaj, 1, FindTipLn, n.pts, mad.tol, qntl, 
#                                    fit.right)))
#   #colnames(tip.tbl) <- c("intercept", "slope", "sharpest.slope", "min.local", 
#   #                       "max.local", "max.global", "max.peak", "fit.ind.ini", "fit.ind.end")
#   colnames(tip.tbl) <- c("subpixel", "slope", "pixel", "min.local", 
#                          "max.local", "max.global", "max.peak", "fit.ind.ini", "fit.ind.end")
#   
#   if(fix.slope){
#     slope <- median(tip.tbl[, 2], na.rm = TRUE)
#     intercept <- sapply(1:dim(kymo)[1], function(i) mean(tip.tbl[i, 8]:tip.tbl[i, 9] - slope * imaj[i, tip.tbl[i, 8]:tip.tbl[i, 9]], na.rm = TRUE))
#     tip.tbl[, 1] <- intercept
#     tip.tbl[, 2] <- slope
#     # WHY mean(fit.indxs - slope * t.slice[fit.indxs])??
#   }
#   out.lst <- list("tip.tbl" = tip.tbl, "kymo.smth" = imaj)
#   return(out.lst)
# }
FindTip <- function(kymo, use.smooth = TRUE, kymo.span.n = 7, kymo.loess.dg = 2, 
  n.pts = 5, mad.tol = 3.5, qntl = 0.95, fit.right = FALSE, fix.slope = FALSE,
  fluo.frac = 0.25, min.chunk.size = 7, coarse.kymo.span.n = 50, coarse.kymo.loess.dg = 2,
  use.coarse = FALSE, rmv.min = TRUE, algorithm = "simple"){
  
  
  if(algorithm == "simple"){
    if(use.coarse){
      imaj <- t(apply(kymo, 1, SmoothWithLoess,
        spn = coarse.kymo.span.n / (dim(kymo)[2]),
        degree = coarse.kymo.loess.dg))
      
      if(use.smooth){
        imaj.raw <- t(apply(kymo, 1, SmoothWithLoess,
          spn = kymo.span.n / (dim(kymo)[2]),
          degree = kymo.loess.dg))
      } else{
        imaj.raw <- kymo
      }
      
      if(rmv.min){
        imaj <- t(apply(imaj, 1, function(vec) vec - min(vec)))
        imaj.raw <- t(apply(imaj.raw, 1, function(vec) vec - min(vec)))
      }
      
      tip.tbl <- t(sapply(1:dim(imaj)[1], function(i) SimpleCoarse(imaj[i, ], imaj.raw[i, ], n.pts, min.chunk.size, fluo.frac)))
      imaj <- imaj.raw
    } else{
      if(use.smooth){
        imaj <- t(apply(kymo, 1, SmoothWithLoess,
          spn = kymo.span.n / (dim(kymo)[2]),
          degree = kymo.loess.dg))
      } else{
        imaj <- kymo
      }
      
      if(rmv.min){
        imaj <- t(apply(imaj, 1, function(vec) vec - min(vec)))
      }
      
      tip.tbl <- t(as.data.frame(apply(imaj, 1, Simple, n.pts, min.chunk.size, fluo.frac)))
      
    }
  } else if(algorithm == "mad"){
    
  } else if(algorithm == "oldonline"){
    if(use.smooth){
      imaj <- t(apply(kymo, 1, SmoothWithLoess,
        spn = kymo.span.n / (dim(kymo)[2]),
        degree = kymo.loess.dg))
    } else{
      imaj <- kymo
    }
    
    if(rmv.min){
      imaj <- t(apply(imaj, 1, function(vec) vec - min(vec)))
    }
    
    tip.tbl <- t(as.data.frame(apply(imaj, 1, OldOnline, n.pts, mad.tol, qntl, 
      fit.right)))
    
  } else if(algorithm == "oldgit"){
    
  } else{
    # Ideas
    # Best fit
    # Constant region fit
    # Closest intercept to local average
  }
  
  colnames(tip.tbl) <- c("subpixel", "slope", "pixel", "min.local", 
    "max.local", "max.global", "max.peak", "fit.ind.ini", "fit.ind.end")
  
  # 
  out.lst <- list("tip.tbl" = tip.tbl, "kymo.smth" = imaj)
  return(out.lst)
}

Simple <- function(t.slice, n.pts = 5, min.chunk.size = 7, fluo.frac = 0.25){
  
  # Fluorescence difference series
  fluo.diff <- diff(t.slice)
  fluo.diff.global.min <- which.min(fluo.diff)
  # Find inis and ends of decrease chunks
  inis <- which(diff(fluo.diff < 0) == 1) + 1
  ends <- which(diff(fluo.diff < 0) == -1) + 1
  chunk.ind <- 1:length(inis)
  
  if(length(ends) == 0){
    ends <- length(t.slice)
  }
  
  # Monotonic signal?
  if(length(inis) == 0){
    # Ever-decreasing?
    if(all(fluo.diff < 0)){
      inis <- fluo.diff.global.min #1                    
      ends <- fluo.diff.global.min #length(fluo.diff)    
    } else{
      # Ever-increasing
      warning("No fluorescence decrease in the time slice!")
      return(cbind("intercept" = NA, 
        "slope" = NA,
        "sharpest.slope" = NA,
        "chunk.min" = NA,
        "chunk.max" = NA,
        "max.global" = NA,
        "max.peak" = NA,
        "fit.ind.ini" = NA,
        "fit.ind.end" = NA))
    }
  }
  
  # Last initial fall greater than end? 
  if(inis[length(inis)] > ends[length(ends)]){
    ends <- c(ends, length(t.slice))
  }
  
  # First initial fall larger than end?
  if(inis[1] > ends[1]){
    inis <- c(1, inis)
  }
  
  ###
  if(length(inis) == 1){
    big.drop <- 1
  } else{
    # Remove Zero size chunks
    # SIZE
    chunk.sizes <- ends - inis
    # Remove small chunks 
    size.min <- which(chunk.sizes >= min.chunk.size)
    #abline(v=inis[size.min])
    #abline(v=ends[size.min])
    
    # FLUO
    chunk.max.fluo <- sapply(1:length(inis), function(i) max(t.slice[inis[i]:ends[i]], na.rm = T))
    #fluo.min <- which(chunk.max.fluo > fluo.frac * max(t.slice))
    fluo.min <- which(chunk.max.fluo > fluo.frac * max(chunk.max.fluo))
    
    # SLOPE
    # Calculate fluo decrease chunk slope
    chunk.slope <- sapply(1:length(inis), function(i) min(fluo.diff[inis[i]:ends[i]], na.rm = T)) #median
    # Remove small slope  chunks 
    slope.min <- which(chunk.slope < as.numeric(lm(t.slice ~ c(1:length(t.slice)))$coefficients[2]))
    
    if((length(size.min) > 0) & (length(fluo.min) > 0)){
      filter.1 <- fluo.min[fluo.min%in%size.min]
    } else{
      filter.1 <- union(size.min, fluo.min)
    }
    
    if(length(filter.1) > 1){
      if(length(slope.min) > 0){
        filter.2 <- slope.min[slope.min%in%filter.1]
      } else{
        filter.2 <- filter.1
      }
      if(length(filter.2) > 0){
        #big.drop <- filter.2[which.min(chunk.min.diff[filter.2])]
        big.drop <- filter.2[length(filter.2)]
      } else{
        big.drop <- filter.1[which.max(chunk.sizes[filter.1])]
        #big.drop <- filter.1[length(filter.1)]
      }
      
    } else if(length(filter.1) == 1){
      big.drop <- filter.1
    } else{
      filter.1 <- which.max(chunk.sizes)
    }
    
  }
  
  if((length(big.drop) < 1) | is.na(big.drop) | is.null(big.drop)){
    message("Giving up, using highest difference signal")
    tip.ini <- which.min(fluo.diff) - n.pts
    tip.end <- which.min(fluo.diff) + n.pts
    
    #fit.center <- which.min(fluo.diff)
    fit.vec.fine <- tip.ini:tip.end
    
  } else{
    # Get indexes from big drop
    tip.ini <- inis[big.drop]
    tip.end <- ends[big.drop]
    fluo.indxs <- tip.ini:tip.end
    
    # CHECK INDEX LIMITS
    coarse.fit.vec <- tip.ini:(tip.end - (min.chunk.size - 1))
    coarse.fit.diff <- sapply(coarse.fit.vec, function(i) as.numeric(median(fluo.diff[i:(i + (min.chunk.size - 1))])))
    #coarse.diff.thrsh <- median(coarse.fit.diff)
    
    first.min.slope <- max(which(coarse.fit.diff == coarse.fit.diff[which.min(coarse.fit.diff)]))
    #}
    fit.vec.fine <- max(tip.ini, (coarse.fit.vec[first.min.slope] - n.pts)):min(tip.end, coarse.fit.vec[first.min.slope] + n.pts)#+ n.pts) #remove n.pts from ends? 
  }
  
  fine.fit.diff <- sapply(fit.vec.fine, function(i) as.numeric(mean(fluo.diff[i:(i + (n.pts - 1))])))
  
  fit.indxs.ini <- fit.vec.fine[max(which(fine.fit.diff == fine.fit.diff[which.min(fine.fit.diff)]))] + 1#+ 1 #+1??
  
  # FIT A LINEAR MODEL
  fit.indxs <- fit.indxs.ini:(fit.indxs.ini + (n.pts - 1))
  vals <- t.slice[fit.indxs]
  mdl <- lm(fit.indxs ~ vals)
  
  peaks <- which(diff(sign(fluo.diff)) == -2) + 1
  tip.tbl <- cbind("intercept" = mdl$coefficients[1], 
    "slope" = mdl$coefficients[2],
    "sharpest.slope" =  fit.indxs[which.min(fluo.diff[fit.indxs])] + 1,
    "chunk.min" = max(fluo.indxs),
    "chunk.max" = min(fluo.indxs),
    "max.global" = which.max(t.slice),
    "max.peak" = peaks[which.max(t.slice[peaks])],
    "fit.ind.ini" = min(fit.indxs),
    "fit.ind.end" = max(fit.indxs)
  )
  return(tip.tbl)
}

SimpleCoarse <- function(t.slice, t.slice.raw, n.pts = 5, min.chunk.size = 7, fluo.frac = 0.25){
  
  # Fluorescence difference series
  fluo.diff <- diff(t.slice)
  fluo.diff.global.min <- which.min(fluo.diff)
  # Find inis and ends of decrease chunks
  inis <- which(diff(fluo.diff < 0) == 1) + 1
  ends <- which(diff(fluo.diff < 0) == -1) + 1
  chunk.ind <- 1:length(inis)
  
  if(length(ends) == 0){
    ends <- length(t.slice)
  }
  
  # Monotonic signal?
  if(length(inis) == 0){
    # Ever-decreasing?
    if(all(fluo.diff < 0)){
      inis <- fluo.diff.global.min #1                    
      ends <- fluo.diff.global.min #length(fluo.diff)    
    } else{
      # Ever-increasing
      warning("No fluorescence decrease in the time slice!")
      return(cbind("intercept" = NA, 
        "slope" = NA,
        "sharpest.slope" = NA,
        "chunk.min" = NA,
        "chunk.max" = NA,
        "max.global" = NA,
        "max.peak" = NA,
        "fit.ind.ini" = NA,
        "fit.ind.end" = NA))
    }
  }
  
  # Last initial fall greater than end? 
  if(inis[length(inis)] > ends[length(ends)]){
    ends <- c(ends, length(t.slice))
  }
  
  # First initial fall larger than end?
  if(inis[1] > ends[1]){
    inis <- c(1, inis)
  }
  
  ###
  if(length(inis) == 1){
    big.drop <- 1
  } else{
    # Remove Zero size chunks
    # SIZE
    chunk.sizes <- ends - inis
    # Remove small chunks 
    size.min <- which(chunk.sizes >= min.chunk.size)
    #abline(v=inis[size.min])
    #abline(v=ends[size.min])
    
    # FLUO
    chunk.max.fluo <- sapply(1:length(inis), function(i) max(t.slice[inis[i]:ends[i]], na.rm = T))
    #fluo.min <- which(chunk.max.fluo > fluo.frac * max(t.slice))
    fluo.min <- which(chunk.max.fluo > fluo.frac * max(chunk.max.fluo))
    
    # SLOPE
    # Calculate fluo decrease chunk slope
    chunk.slope <- sapply(1:length(inis), function(i) min(fluo.diff[inis[i]:ends[i]], na.rm = T)) #median
    # Remove small slope  chunks 
    slope.min <- which(chunk.slope < as.numeric(lm(t.slice ~ c(1:length(t.slice)))$coefficients[2]))
    
    if((length(size.min) > 0) & (length(fluo.min) > 0)){
      filter.1 <- fluo.min[fluo.min%in%size.min]
    } else{
      filter.1 <- union(size.min, fluo.min)
    }
    
    if(length(filter.1) > 1){
      if(length(slope.min) > 0){
        filter.2 <- slope.min[slope.min%in%filter.1]
      } else{
        filter.2 <- filter.1
      }
      if(length(filter.2) > 0){
        #big.drop <- filter.2[which.min(chunk.min.diff[filter.2])]
        big.drop <- filter.2[length(filter.2)]
      } else{
        big.drop <- filter.1[which.max(chunk.sizes[filter.1])]
        #big.drop <- filter.1[length(filter.1)]
      }
      
    } else if(length(filter.1) == 1){
      big.drop <- filter.1
    } else{
      filter.1 <- which.max(chunk.sizes)
    }
    
  }
  
  if((length(big.drop) < 1) | is.na(big.drop) | is.null(big.drop)){
    message("Giving up, using highest difference signal")
    tip.ini <- which.min(fluo.diff) - n.pts
    tip.end <- which.min(fluo.diff) + n.pts
    
    #fit.center <- which.min(fluo.diff)
    fit.vec.fine <- tip.ini:tip.end
    
  } else{
    # Get indexes from big drop
    tip.ini <- inis[big.drop]
    tip.end <- ends[big.drop]
    fluo.indxs <- tip.ini:tip.end
    
    t.slice <- t.slice.raw
    fluo.diff <- diff(t.slice)
    # CHECK INDEX LIMITS
    coarse.fit.vec <- tip.ini:(tip.end - (min.chunk.size - 1))
    coarse.fit.diff <- sapply(coarse.fit.vec, function(i) as.numeric(median(fluo.diff[i:(i + (min.chunk.size - 1))])))
    #coarse.diff.thrsh <- median(coarse.fit.diff)
    
    first.min.slope <- max(which(coarse.fit.diff == coarse.fit.diff[which.min(coarse.fit.diff)]))
    #}
    fit.vec.fine <- max(tip.ini, (coarse.fit.vec[first.min.slope] - n.pts)):min(tip.end, coarse.fit.vec[first.min.slope] + n.pts)#+ n.pts) #remove n.pts from ends? 
  }
  
  fine.fit.diff <- sapply(fit.vec.fine, function(i) as.numeric(mean(fluo.diff[i:(i + (n.pts - 1))])))
  
  fit.indxs.ini <- fit.vec.fine[max(which(fine.fit.diff == fine.fit.diff[which.min(fine.fit.diff)]))] + 1#+ 1 #+1??
  
  # FIT A LINEAR MODEL
  fit.indxs <- fit.indxs.ini:(fit.indxs.ini + (n.pts - 1))
  vals <- t.slice[fit.indxs]
  mdl <- lm(fit.indxs ~ vals)
  #}
  
  peaks <- which(diff(sign(fluo.diff)) == -2) + 1
  tip.tbl <- cbind("intercept" = mdl$coefficients[1], 
    "slope" = mdl$coefficients[2],
    "sharpest.slope" =  fit.indxs[which.min(fluo.diff[fit.indxs])] + 1,
    "chunk.min" = max(fluo.indxs),
    "chunk.max" = min(fluo.indxs),
    "max.global" = which.max(t.slice),
    "max.peak" = peaks[which.max(t.slice[peaks])],
    "fit.ind.ini" = min(fit.indxs),
    "fit.ind.end" = max(fit.indxs)
  )
  #print(dim(tip.tbl))
  return(tip.tbl)
}

OldOnline <- function(kymo, n.pts = 5, mad.tol = 3.5, qntl = 0.95, fit.right = FALSE){
  
  kymo <- as.matrix(apply(kymo, 2, as.numeric))
  use.smooth <- TRUE
  kymo.span.n <- 7
  kymo.loess.dg <- 1
  fix.slope <- FALSE
  
  # Smooth image
  if(use.smooth){
    imaj <- t(apply(kymo, 1, SmoothWithLoess, #INTEGRATE VERSION?
      spn = kymo.span.n / (dim(kymo)[2]),
      degree = kymo.loess.dg))
  } else{imaj <- kymo}
  
  # Subtract minimium value
  imaj <- imaj - min(imaj, na.rm = TRUE)
  
  tip.tbl <- t(as.data.frame(apply(imaj, 1, FindTipLn, n.pts, mad.tol, qntl, 
    fit.right)))
  #colnames(tip.tbl) <- c("intercept", "slope", "sharpest.slope", "min.local", 
  #                       "max.local", "max.global", "max.peak", "fit.ind.ini", "fit.ind.end")
  colnames(tip.tbl) <- c("subpixel", "slope", "pixel", "min.local", 
    "max.local", "max.global", "max.peak", "fit.ind.ini", "fit.ind.end")
  
  if(fix.slope){
    slope <- median(tip.tbl[, 2], na.rm = TRUE)
    intercept <- sapply(1:dim(kymo)[1], function(i) mean(tip.tbl[i, 8]:tip.tbl[i, 9] - slope * imaj[i, tip.tbl[i, 8]:tip.tbl[i, 9]], na.rm = TRUE))
    tip.tbl[, 1] <- intercept
    tip.tbl[, 2] <- slope
    # WHY mean(fit.indxs - slope * t.slice[fit.indxs])??
  }
  #out.lst <- list("tip.tbl" = tip.tbl, "kymo.smth" = imaj)
  return(tip.tbl)
}

#-------------------------------------------------------------------------------
FindTipLn <- function(t.slice, n.pts = 5, mad.tol = 3.5, qntl = 0.95, fit.right = FALSE){
  
  # Fluorescence difference series
  fluo.diff <- diff(t.slice)
  
  # Find inis and ends of decrease chunks
  inis <- which(diff(fluo.diff < 0) == 1) + 1
  ends <-  which(diff(fluo.diff < 0) == -1) + 1
  
  # Monotonic signal?
  if(length(inis) == 0){
    # Ever-decreasing?
    if(all(fluo.diff < 0)){
      inis <- which.min(fluo.diff) #1                    
      ends <- which.min(fluo.diff) #length(fluo.diff)    
    } else{
      # Ever-increasing
      warning("No fluorescence decrease in the time slice!")
    }
  }
  
  # Last initial fall greater than end? 
  if(inis[length(inis)] > ends[length(ends)]){
    ends <- c(ends, length(fluo.diff))
  }
  
  # First initial fall larger than end?
  if(inis[1] > ends[1]){
    inis <- c(1, inis)
  }
  
  # Calculate fluo decrease chunk size
  chunk.sizes <- ends - inis
  # Calculate fluo decrease chunk max fluo
  chunk.max.fluo <- sapply(1:length(inis), function(i) max(t.slice[inis[i]:ends[i]]))
  # Calculate fluo decrease chunk slope
  chunk.slope <- sapply(1:length(inis), function(i) mean(diff(t.slice[inis[i]:ends[i]])))
  
  # Find chunk indexes of large drops
  sz.indx <- FilterObs(chunk.sizes, mad.tol, qntl)
  fluo.indx <- FilterObs(chunk.max.fluo, mad.tol, qntl)
  slope.indx <- FilterObs(-chunk.slope, mad.tol, qntl)
  
  # Establish minimum values
  sz.min <- which(chunk.sizes > n.pts)
  fluo.min <- which(chunk.max.fluo > max(t.slice[ends[length(ends)]:length(t.slice)]))
  slope.min <- which(chunk.slope < as.numeric(lm(t.slice ~ c(1:length(t.slice)))$coefficients[2]))
  
  # Count reocurring chunk index
  union.count <- table(c(sz.indx, fluo.indx, slope.indx, sz.min, fluo.min, slope.min))
  # Labels of reocurring chunk index
  unique.indxs <- as.numeric(names(union.count))
  # Get the maximum of the most frequently reocurring index
  big.drop <- unique.indxs[max(which(union.count == max(union.count)))]
  
  if(is.na(big.drop) | is.null(big.drop)){
    message("Giving up, using highest difference signal")
    big.drop <- which.min(fluo.diff)
  }
  
  # Get indexes from big drop
  tip.ini <- inis[big.drop]
  tip.end <- ends[big.drop]
  fluo.indxs <- tip.ini:tip.end
  
  if(length(fluo.indxs) > 1){
    #max.local <- fluo.indxs[1]
    inflct.pts.down <- which(diff(sign(diff(t.slice[fluo.indxs], differences = 2))) == 2) + 1
    sharpest.slope <- fluo.indxs[inflct.pts.down[which.min(diff(t.slice[fluo.indxs])[inflct.pts.down])]]
    if(fit.right){
      sharpest.slope <- fluo.indxs[max(inflct.pts.down)]
    }
  } else{sharpest.slope <- which.min(fluo.diff)}
  fit.vec.min <- max(1, (sharpest.slope - n.pts))
  fit.vec.max <- min((sharpest.slope + 1), length(t.slice) - n.pts)
  fit.vec <- fit.vec.min:fit.vec.max
  fit.slope <- sapply(fit.vec, function(i) sum(diff(t.slice[i:(i + (n.pts - 1))])))
  fit.indxs <- fit.vec[which.min(fit.slope)]:(fit.vec[which.min(fit.slope)] + (n.pts - 1))
  
  # FIT A LINEAR MODEL
  vals <- t.slice[fit.indxs]
  mdl <- lm(fit.indxs ~ vals)
  
  peaks <- which(diff(sign(fluo.diff)) == -2) + 1
  tip.tbl <- cbind("intercept" = mdl$coefficients[1], 
    "slope" = mdl$coefficients[2],
    "sharpest.slope" = sharpest.slope,
    "chunk.min" = max(fluo.indxs),
    "chunk.max" = min(fluo.indxs),
    "max.global" = which.max(t.slice),
    "max.peak" = peaks[which.max(t.slice[peaks])],
    "fit.ind.ini" = min(fit.indxs),
    "fit.ind.end" = max(fit.indxs)
  )
  
  return(tip.tbl)
}
#------------------------------------------------------------------------------- 
# RefineTipLoc <- function(tip.loc, tip.span.n = 7, tip.span.dg = 1, 
#                          rm.tip.out = FALSE, rm.tip.out.dg = 1, 
#                          rm.tip.out.spn = 0.4, tol = 3.5, px.tol = 7){
#   
#   if(rm.tip.out){
#     tip.loc <- RemoveTipOutliers(tip.loc, dg = rm.tip.out.dg, spn = rm.tip.out.spn, tol, px.tol)
#   }
#   tip.loc.smth <- SmoothWithLoess(tip.loc,  spn = tip.span.n / length(tip.loc), degree = tip.span.dg)
#   
#   #tip.loc.tbl <- cbind.data.frame("tip.loc.raw" = tip.loc, 
#   #"tip.loc.smth" = tip.loc.smth)
#   return(tip.loc.smth)
# }
RefineTipLoc <- function(tip.loc, tip.span.n = 7, tip.span.dg = 1, 
  rm.tip.out = FALSE, rm.tip.out.dg = 1, 
  rm.tip.out.spn = 0.4, tol = 3.5, px.tol = 7, algorithm = "maddiff"){
  
  #tip.loc <- RemoveGrossTipOutliers(tip.loc, dg = rm.tip.out.dg, spn = rm.tip.out.spn, tol, px.tol)
  out.indx <- c()
  if(rm.tip.out){
    # Perform 1st loess fit
    tip.fit <- tip.loc
    x.fit <- c(1:length(tip.loc))
    loess.mdl <-  loess(tip.fit ~ x.fit, span = rm.tip.out.spn, control = loess.control(surface = "direct"),
      family = "gaussian", degree = rm.tip.out.dg)
    
    #Algorithms: "maddiff", "pixeltol", tso", 
    if(algorithm == "maddiff"){
      mdl.diff <- loess.mdl$fitted - tip.loc
      if(is.character(px.tol) | is.null(px.tol)){
        out.indx <- which(abs(mdl.diff) > tol * mad(mdl.diff))
      } else{
        
        px.indx <- which(abs(diff(tip.loc)) > px.tol)
        if(length(px.indx) >= 2){
          reps <- which(diff(px.indx) == 1)
          if(length(reps) > 0){
            px.indx <- px.indx[-reps]
          }
        }
        
        out.indx <- unique(c(which(abs(mdl.diff) > tol * mad(mdl.diff)), 
          px.indx))
      }
      
    } else if(algorithm == "pixeltol"){
      px.indx <- which(abs(diff(tip.loc)) > px.tol)
      if(length(px.indx) >= 2){
        reps <- which(diff(px.indx) == 1)
        if(length(reps) > 0){
          px.indx <- px.indx[-reps]
        }
      }
      
      out.indx <- px.indx
    }
    
    if(length(out.indx) > 0){
      tip.fit <- tip.loc[-out.indx]
      x.fit <- c(1:length(tip.loc))[-out.indx]
      # Perform 2ND loess fit
      loess.mdl <-  loess(tip.fit ~ x.fit, span = tip.span.n / length(tip.loc), 
        control = loess.control(surface = "direct"),
        family = "gaussian", degree = tip.span.dg)
      
      # Approximate to interpolation times
      mdl.interp <- predict(loess.mdl, out.indx, se = F)
      tip.loc[out.indx] <- mdl.interp
      message(paste("Outlier removed from", out.indx,"; "))
    }
    
  }
  if(tip.span.n > 0){
    tip.loc.smth <- SmoothWithLoess(tip.loc,  spn = tip.span.n / length(tip.loc), degree = tip.span.dg)
  } else{
    tip.loc.smth <- tip.loc 
  }
  
  #tip.loc.tbl <- cbind.data.frame("tip.loc.raw" = tip.loc, 
  #"tip.loc.smth" = tip.loc.smth)
  return(list("tip.loc" = tip.loc, "tip.loc.smth" = tip.loc.smth, "out.indx" = out.indx))
} 
#-------------------------------------------------------------------------------
ExtractFluoTimeSeries <-function(kymo.align, max.tip, roi.px, avg.width = 5, tip.mrgn = 0, use.median = TRUE){
  roi.mid.px <- round(roi.px, digits = 0)
  roi.up.px <- roi.mid.px + floor((avg.width - 1) / 2)
  roi.down.px <- roi.mid.px - ceiling((avg.width - 1) / 2)
  edge.roi.up <- max(roi.up.px, floor(max.tip) - 1)
  edge.roi.down <- edge.roi.up - (avg.width - 1)
  center.roi.up <- roi.up.px + round((edge.roi.up - roi.up.px) / 2, digits = 0) - 1
  center.roi.down <- center.roi.up - (avg.width - 1)
  
  roi.ts <- GetFluoTimeSeries(kymo.align, ini.ind = roi.down.px, avg.width = avg.width, use.median = use.median)
  center.ts <- GetFluoTimeSeries(kymo.align, ini.ind = center.roi.down, avg.width = avg.width, use.median = use.median)
  edge.ts <- GetFluoTimeSeries(kymo.align, ini.ind = edge.roi.down, avg.width = avg.width, use.median = use.median)
  
  ini.px <- cbind("roi" = roi.down.px, "center" = center.roi.down, "edge" = edge.roi.down)
  fluo.ts <- cbind("roi.ts" = roi.ts, "center.ts" = center.ts, "edge.ts" = edge.ts)
  
  return(list("ini.px" = ini.px,
    "fluo.ts" = fluo.ts))
  
  
}
#-------------------------------------------------------------------------------
FilterKymo <- function(kymo, time.step, low.per, high.per){
  sampling.freq <- 1 / time.step
  series.length <- dim(kymo)[1] * time.step
  smth.vec <- GetSmoothVec(sampling.freq, series.length, low.per, high.per)
  lvls <- length(smth.vec)
  imaj.mra <- waveslim::mra.2d(kymo, wf = "la8", J = lvls, method = "modwt",
    boundary = "periodic")
  
  selected.level <- grepl("HL", names(imaj.mra))
  
  out.imaj <- matrix(0, nrow = dim(kymo)[1], ncol = dim(kymo)[2])
  
  i <- 1
  j <- 1
  for(detail in imaj.mra){
    if(selected.level[i]){
      if(smth.vec[j]){out.imaj <- out.imaj + detail}
      j <- j + 1
    }
    i <- i + 1
  }
  return(out.imaj)
}
#-------------------------------------------------------------------------------
### PLOTTING FUNCTIONS----------------------------------------------------------
#-------------------------------------------------------------------------------
PlotKymoInput <- function(kymo, px.sz = 1, px.unit = "pixel", 
  time.step = 1, time.unit = "frame", 
  brks = NULL, red.bias = 1){
  if(!is.null(brks)){
    if(length(brks) == 2){
      brks <- seq(from = min(brks), to = max(brks), length.out = 257)
    } else{
      brks <- seq(from = min(c(kymo)), to = max(c(kymo)), length.out = 257)
    }
  } else{
    brks <- seq(from = min(c(kymo)), to = max(c(kymo)), length.out = 257)
  }
  
  if(px.unit == "pixel"){
    xs <- 1:dim(kymo)[2]
  } else{
    xs <- (1:dim(kymo)[2] - 1) * px.sz
  }
  
  if(time.unit == "frame"){
    ys <- 1:dim(kymo)[1]
  } else{
    ys <- (1:dim(kymo)[1] - 1) * time.step
  }
  
  pal <- colorRampPalette(tim.colors(256), bias = red.bias)
  
  plt.out <-  fields::image.plot(x = xs, 
    y = ys, 
    z = t(kymo), 
    col = pal(256), #fields::tim.colors
    ylim = rev(range(ys)),
    xlab = paste("Length (", px.unit, ")", sep = ""),
    ylab = paste("Time (", time.unit, ")", sep = ""), 
    breaks = brks)
  
  return(plt.out)
}
#-------------------------------------------------------------------------------
PlotKymoWithAllTip <- function(kymo, tip.tbl, px.sz = 1, px.unit = "pixel", 
  time.step = 1, time.unit = "frame", 
  brks = NULL, red.bias = 1, restrict.x = TRUE, cx = 0.35){
  if(!is.null(brks)){
    if(length(brks) == 2){
      brks <- seq(from = min(brks), to = max(brks), length.out = 257)
    } else{
      brks <- seq(from = min(c(kymo)), to = max(c(kymo)), length.out = 257)
    }
  } else{
    brks <- seq(from = min(c(kymo)), to = max(c(kymo)), length.out = 257)
  }
  
  if(px.unit == "pixel"){
    xs <- 1:dim(kymo)[2]
  } else{
    xs <- (1:dim(kymo)[2] - 1) * px.sz
  }
  
  if(time.unit == "frame"){
    ys <- 1:dim(kymo)[1]
  } else{
    ys <- (1:dim(kymo)[1] - 1) * time.step
  }
  
  if(restrict.x){
    xlm <- range(c(tip.tbl[,-2])) + c(-1,1)
  } else{
    xlm <- c(1,dim(kymo)[2])
  }
  
  
  pal <- colorRampPalette(tim.colors(256), bias = red.bias)
  
  fields::image.plot(x = xs, y = ys, 
    z = t(kymo), 
    col = pal(256), #fields::tim.colors
    ylim = rev(range(ys)),
    xlab = paste("Length (", px.unit, ")", sep = ""),
    ylab = paste("Time (", time.unit, ")", sep = ""), 
    breaks = brks, 
    xlim = xlm * px.sz)
  
  indxs <- c(1, 3:7)
  lbls <- colnames(tip.tbl)[indxs]
  clrs <- c("magenta", "black", "grey40", adjustcolor("white", alpha.f = 0.7), "cyan", "purple")
  cxs <- c(cx*1, rep(cx, length(lbls) - 1))
  
  invisible(sapply(c(1:3,5:length(indxs),4),
    function(indx, tip.tbl, px.sz, ys, cxs, clrs) points(tip.tbl[, which(colnames(tip.tbl) == lbls[indx])] * px.sz,
      ys, col =clrs[indx], cex=cxs[indx], pch = 19),
    tip.tbl, px.sz, ys, cxs, clrs))
  
  legend("topright", legend = lbls, col = "black", pch = 21, pt.cex = cxs * 5, pt.bg = clrs)
  return(NULL)
}
#-------------------------------------------------------------------------------
# PlotKymoWithTip <- function(kymo, tip.loc.raw, tip.loc.smth, px.sz = 1, px.unit = "pixel", 
#                             time.step = 1, time.unit = "frame", 
#                             brks = NULL, red.bias = 1, cx = 0.5, cex.ax = 1, ln = 0.75, rmv.lgnd = FALSE){
#   if(!is.null(brks)){
#     if(length(brks) == 2){
#       brks <- seq(from = min(brks), to = max(brks), length.out = 257)
#     } else{
#       brks <- seq(from = min(c(kymo)), to = max(c(kymo)), length.out = 257)
#     }
#   } else{
#     brks <- seq(from = min(c(kymo)), to = max(c(kymo)), length.out = 257)
#   }
#   
#   if(px.unit == "pixel"){
#     xs <- 1:dim(kymo)[2]
#   } else{
#     xs <- (1:dim(kymo)[2] - 1) * px.sz
#   }
#   
#   if(time.unit == "frame"){
#     ys <- 1:dim(kymo)[1]
#   } else{
#     ys <- (1:dim(kymo)[1] - 1) * time.step
#   }
#   
#   pal <- colorRampPalette(tim.colors(256), bias = red.bias)
#   
#   margn.top <- c(0,6,5,8)
#   #margn.middle <- c(0,6,0,8)
#   margn.bottom <- c(6,6,0,8)
#   par(mfrow = c(2, 1))
#   
#   par(mar = margn.top)
#   
#   image(x = ys, y = xs, 
#                      z = kymo, 
#                      col = pal(256), #fields::tim.colors
#                      xlim = range(ys),
#                      #ylab = paste("Length (", px.unit, ")", sep = ""),
#                      #xlab = paste("Time (", time.unit, ")", sep = ""), 
#                      breaks = brks, 
#                      ylim = (range(c(tip.loc.raw, tip.loc.smth), na.rm = TRUE) + c(-2, 2)) * px.sz,
#                      xlab = "",
#                      ylab ="",
#                      xaxt = "n", yaxt = "n")
#   
#   axis(3, cex.axis = cex.ax)
#   axis(2, cex.axis = cex.ax)
#   mtext(paste("Time (", time.unit, ")", sep = ""), side = 3, line = 3, cex = cex.ax)
#   mtext(paste("Length (", px.unit, ")", sep = ""), side =2, line = 3, cex = cex.ax)
#   points(ys, tip.loc.raw * px.sz, col = "white", cex = cx, pch = 19)
#   lines(ys, tip.loc.smth * px.sz, col = "black", lwd = ln)
#   
#   if(!rmv.lgnd){
#     legend("topleft", legend = c("Tip estimate raw", "Tip estimate smooth"), col = "black", pch = 21, pt.cex = c(cx*1.5, NA), pt.bg = c("white", NA), lwd = c(NA, ln))
#   }
#   
#   par(mar = margn.bottom)
#   plot(ys[-1], diff(tip.loc.smth * px.sz) / time.step, type = "n", xlab = "", ylab ="", xaxt = "n", yaxt = "n", lwd = 2.5)
#   points(ys[-1], diff(tip.loc.raw * px.sz) / time.step, type = "o", col = "grey", lwd = ln, pch = 19, cex = cx)
#   lines(ys[-1], diff(tip.loc.smth * px.sz) / time.step, type = "l", col = "black", lwd = ln * 2.5)
#   axis(1, cex.axis = cex.ax)
#   axis(2, cex.axis = cex.ax)
#   mtext(paste("Time (", time.unit, ")", sep = ""), side = 1, line = 3, cex = cex.ax)
#   mtext(paste("Growth rate (", px.unit, "/", time.unit, ")", sep = ""), side =2, line = 3, cex = cex.ax)
#   if(!rmv.lgnd){
#   legend("bottomleft", c("Raw", "Smooth"), pch = c(19, NA), pt.cex = c(cx*1.5, NA), lwd = c(ln, ln * 2.5), col = c("grey", "black"))
#   }
#   return(NULL)
# }
PlotKymoWithTip <- function(kymo, tip.loc.raw, tip.loc.smth, tip.loc.out,
  out.indx, px.sz = 1, px.unit = "pixel", 
  time.step = 1, time.unit = "frame", 
  brks = NULL, red.bias = 1, cx = 0.5, cex.ax = 1, ln = 0.75, rmv.lgnd = FALSE){
  if(!is.null(brks)){
    if(length(brks) == 2){
      brks <- seq(from = min(brks), to = max(brks), length.out = 257)
    } else{
      brks <- seq(from = min(c(kymo)), to = max(c(kymo)), length.out = 257)
    }
  } else{
    brks <- seq(from = min(c(kymo)), to = max(c(kymo)), length.out = 257)
  }
  
  if(px.unit == "pixel"){
    xs <- 1:dim(kymo)[2]
  } else{
    xs <- (1:dim(kymo)[2] - 1) * px.sz
  }
  
  if(time.unit == "frame"){
    ys <- 1:dim(kymo)[1]
  } else{
    ys <- (1:dim(kymo)[1] - 1) * time.step
  }
  
  pal <- colorRampPalette(tim.colors(256), bias = red.bias)
  
  margn.top <- c(0,6,5,8)
  #margn.middle <- c(0,6,0,8)
  margn.bottom <- c(6,6,0,8)
  par(mfrow = c(2, 1))
  
  par(mar = margn.top)
  
  image(x = ys, y = xs, 
    z = kymo, 
    col = pal(256), #fields::tim.colors
    xlim = range(ys),
    #ylab = paste("Length (", px.unit, ")", sep = ""),
    #xlab = paste("Time (", time.unit, ")", sep = ""), 
    breaks = brks, 
    ylim = (range(c(tip.loc.raw, tip.loc.smth), na.rm = TRUE) + c(-2, 2)) * px.sz,
    xlab = "",
    ylab ="",
    xaxt = "n", yaxt = "n")
  
  axis(3, cex.axis = cex.ax)
  axis(2, cex.axis = cex.ax)
  mtext(paste("Time (", time.unit, ")", sep = ""), side = 3, line = 3, cex = cex.ax)
  mtext(paste("Length (", px.unit, ")", sep = ""), side =2, line = 3, cex = cex.ax)
  
  points(ys, tip.loc.raw * px.sz, col = "white", cex = cx, pch = 19, type = "o")
  lines(ys, tip.loc.smth * px.sz, col = "darkgrey", lwd = ln + 1)
  
  if(length(out.indx) > 0){
    points(ys[out.indx], tip.loc.out[out.indx] * px.sz, col = "magenta", cex = cx, pch = 19)
    points(ys[out.indx], tip.loc.raw[out.indx] * px.sz, col = "black", cex = cx*2, pch = 4)
  }
  
  if(!rmv.lgnd){
    if(length(out.indx) > 0){
      legend("topleft", legend = c("Tip estimate raw", "Tip estimate smooth", "Outlier", "New value"), col = c("white","darkgrey", "black", "magenta"), pch = c(19,19, 4, 19), pt.cex = c(cx*1.5,0,cx*2.5,cx*1.5), lwd = c(ln,ln + 1,1,0), bty = "n")
      
    } else{
      legend("topleft", legend = c("Tip estimate raw", "Tip estimate smooth"), col = "darkgrey", pch = 21, pt.cex = c(cx*1.5, NA), pt.bg = c("white", NA), lwd = c(ln, ln + 1))
    }
  }
  
  par(mar = margn.bottom)
  par(xaxs = "i")
  plot(ys[-1], diff(tip.loc.smth * px.sz) / time.step, type = "n", xlab = "", ylab ="", xaxt = "n", yaxt = "n", lwd = 2.5)
  points(ys[-1], diff(tip.loc.raw * px.sz) / time.step, type = "o", col = "grey", lwd = ln, pch = 19, cex = cx)
  lines(ys[-1], diff(tip.loc.smth * px.sz) / time.step, type = "l", col = "black", lwd = ln * 2.5)
  
  if(length(out.indx) > 0){
    points(ys[unique(c(out.indx - 1, out.indx))], (diff(tip.loc.out * px.sz) / time.step)[unique(c(out.indx - 1, out.indx))], col = "magenta", cex = cx, pch = 19)
    points(ys[unique(c(out.indx - 1, out.indx))], (diff(tip.loc.raw * px.sz) / time.step)[unique(c(out.indx - 1, out.indx))], col = "black", cex = cx *2, pch = 4)
    abline( v = ys[unique(c(out.indx - 1, out.indx))], col = "magenta", cex = cx)
  }
  
  
  axis(1, cex.axis = cex.ax)
  axis(2, cex.axis = cex.ax)
  mtext(paste("Time (", time.unit, ")", sep = ""), side = 1, line = 3, cex = cex.ax)
  mtext(paste("Growth rate (", px.unit, "/", time.unit, ")", sep = ""), side =2, line = 3, cex = cex.ax)
  if(!rmv.lgnd){
    if(length(out.indx) > 0){
      legend("bottomleft", c("Raw", "Smooth", "Outlier", "New value"), pch = c(19,19, 4, 19), pt.cex = c(cx*1.5,0,cx*2.5,cx*1.5), lwd = c(ln, ln * 2.5), col = c("grey", "black", "black", "magenta"))
      #col = c("white","darkgrey", "black", "magenta"), pch = c(19,19, 4, 19), pt.cex = c(cx*1.5,0,cx*2.5,cx*1.5), lwd = c(ln,ln + 1,1,0)
    } else{
      legend("bottomleft", c("Raw", "Smooth"), pch = c(19, NA), pt.cex = c(cx*1.5, NA), lwd = c(ln, ln * 2.5), col = c("grey", "black")) 
    }
    
  }
  return(NULL)
}

PlotKymoWithDygraph <- function(kymo, tip.loc.raw, tip.loc.smth, px.sz = 1, px.unit = "pixel", 
  time.step = 1, time.unit = "frame", 
  brks = NULL, red.bias = 1, cx = 0.5, cex.ax = 1, ln = 0.75, rmv.lgnd = FALSE) {
  if(time.unit == "frame"){
    ys <- 1:dim(kymo)[1]
  } else{
    ys <- (1:dim(kymo)[1] - 1) * time.step
  }
  
  # First, creating a data.frame with first column being the x-axis
  dygraph.df <- as.data.frame(cbind(x.ax = ys[-1],
    y.smth = diff(tip.loc.smth * px.sz) / time.step,
    y.raw = diff(tip.loc.raw * px.sz) / time.step))
  # Second, creating the plots
  x.label <- paste("Time (", time.unit, ")", sep = "")
  y.label <- paste("Growth rate (", px.unit, "/", time.unit, ")", sep = "")
  
  y.smth_min <- min(dygraph.df$y.smth)
  if(y.smth_min < 0) {
    y.min <- y.smth_min + y.smth_min /3
  }else{
    if(y.smth_min > 0) {
      y.min <- y.smth_min - y.smth_min /3
    }else{
      y.min <- -0.4
    }
  } 
  y.max <- max(dygraph.df$y.smth) + max(dygraph.df$y.smth) / 3
  
  dygraph(dygraph.df, xlab = x.label, ylab = y.label) %>%
    dySeries("y.raw", label = "Raw", color = "lightgrey", drawPoints = T, pointSize = 3) %>%
    dySeries("y.smth", label = "Smoothed", color = "black", strokeWidth = 2) %>%
    dyAxis("x", valueRange = range(dygraph.df$x.ax), rangePad = 10) %>%
    dyAxis("y", valueRange = c(y.min, y.max), rangePad = 10) %>% 
    dyRangeSelector() %>%
    dyOptions(drawGrid = FALSE)
}

#-------------------------------------------------------------------------------
PlotTimeSlice <- function(t.slice.raw, t.slice.smth, tip.ests,restrict.x = FALSE, cx = 0.5){
  if(restrict.x){
    xlm <- range(c(tip.ests[-2])) + c(-1,1)
  } else{
    xlm <- c(1,length(t.slice.raw))
  }
  plot(t.slice.raw, xlim = xlm, xlab = "Pixel (number)", ylab = "Fluorescence (AU)")
  points(t.slice.smth, xlim = xlm, col = "black", pch = 19, cex = cx)
  indxs <- c(1, 3:7)
  lbls <- names(tip.ests)[indxs]
  clrs <- c("magenta", "black", "grey40", "lightgrey", "cyan", "purple")
  #cxs <- c(cx*1, rep(cx, length(lbls) - 1))
  ln <- 2.5 
  
  invisible(sapply(1:length(lbls), 
    function(indx, tip.ests, clrs, ln) abline(v=tip.ests[which(names(tip.ests) == lbls[indx])], col =clrs[indx], lwd = ln), 
    tip.ests,clrs, ln))
  abline(lm(t.slice.smth[tip.ests[8]:tip.ests[9]]~c(tip.ests[8]:tip.ests[9])),col = "red")
  points(tip.ests[8]:tip.ests[9],t.slice.smth[tip.ests[8]:tip.ests[9]], col = "red", lwd = 1.5)
  abline(h = 0, lty = "dotdash")
  
  legend("topright", legend = c("raw signal", "smoothed", "linear fit",lbls), col = c("black", "black", "red",clrs), lwd = c(NA,NA,1,rep(ln, length(lbls))), pch =  c(1,19,1,rep(NA, length(lbls))))
  
  return(NULL)
}
#-------------------------------------------------------------------------------
PlotKymoAndFluo <- function(kymo.align, tip.fluo.ts, roi.ini.pxs, avg.width = 5, 
  max.tip = NULL, px.sz = 1, px.unit = "pixel", 
  time.step = 1, time.unit = "pixel", brks = NULL, 
  red.bias = 1, cex.ax = 1.15, ln = 2, y.lab = "Fluorescence (AU)", rmv.lgnd = FALSE){ 
  
  roi.ini.pos <- roi.ini.pxs * px.sz
  avg.width.pos <- (avg.width - 1) * px.sz
  
  xs <- (1:dim(kymo.align)[1] - 1) * time.step
  ys <- (1:dim(kymo.align)[2] - 1) * px.sz
  
  if(is.null(max.tip)){
    max.tip.pos <- max(ys)
  } else{
    max.tip.pos <- max.tip * px.sz
  }
  
  if(!is.null(brks)){
    if(length(brks) == 2){
      brks <- seq(from = min(brks), to = max(brks), length.out = 257)
    } else{
      brks <- seq(from = min(c(kymo.align)), to = max(c(kymo.align)), length.out = 257)
    }
  } else{
    brks <- seq(from = min(c(kymo.align)), to = max(c(kymo.align)), length.out = 257)
  }
  
  pal <- colorRampPalette(tim.colors(256), bias = red.bias)
  
  margn.top <- c(0,6,5,8)
  margn.middle <- c(0,6,0,8)
  margn.bottom <- c(6,6,0,8)
  par(mfrow = c(3, 1))
  
  par(mar = margn.top)
  
  image(x = xs, y = ys, 
    z = kymo.align, 
    col = pal(256), #fields::tim.colors
    breaks = brks, 
    xlab = "",
    ylab ="",
    xaxt = "n", yaxt = "n",
    ylim = c(0, max.tip.pos))
  
  axis(3, cex.axis = cex.ax, cex.axis = cex.ax)
  axis(2, cex.axis = cex.ax, cex.axis = cex.ax)
  mtext(paste("Time (", time.unit, ")", sep = ""), side = 3, line = 3, cex = cex.ax)
  mtext(paste("Length from tip (", px.unit, ")", sep = ""), side =2, line = 3, cex = cex.ax)
  
  rect(xleft = min(xs) - time.step * 2, 
    ybottom = roi.ini.pos[1],
    xright = max(xs) + time.step * 2,
    ytop = roi.ini.pos[1] + avg.width.pos,
    density = 10,
    col = adjustcolor("black", alpha.f = 0.5)
  )
  
  text(x = min(xs) + (max(xs) - min(xs)) / 2, 
    y = (roi.ini.pos[1] + avg.width.pos) + 4 * avg.width.pos,
    labels = "Tip",
    cex = cex.ax)
  
  rect(xleft = min(xs) - time.step * 2, 
    ybottom = roi.ini.pos[2],
    xright = max(xs) + time.step * 2,
    ytop = roi.ini.pos[2] + avg.width.pos,
    density = 10,
    col = adjustcolor("black", alpha.f = 0.5)
  )
  
  text(x = min(xs) + (max(xs) - min(xs)) / 2, 
    y = (roi.ini.pos[2] + avg.width.pos) + 4 * avg.width.pos,
    labels = "Center",
    cex = cex.ax)
  
  rect(xleft = min(xs) - time.step * 2, 
    ybottom = roi.ini.pos[3],
    xright = max(xs) + time.step * 2,
    ytop = roi.ini.pos[3] + avg.width.pos,
    density = 10,
    col = adjustcolor("black", alpha.f = 0.5)
  )
  
  text(x = min(xs) + (max(xs) - min(xs)) / 2, 
    y = (roi.ini.pos[3] + avg.width.pos) + 4 * avg.width.pos,
    labels = "Edge",
    cex = cex.ax)
  
  if(!rmv.lgnd){
    legend("topleft", "Tip aligned kymograph", text.col = "white", bty = "n")
  }
  
  par(mar = margn.middle)
  par(xaxs = "i")
  
  plot(xs, tip.fluo.ts[, 1], type = "o", pch = 19, cex = 0., col = "firebrick3", lwd = ln, xlab = "", ylab ="", xaxt = "n", yaxt = "n")
  axis(2, cex = cex.ax, cex.axis = cex.ax)
  mtext(y.lab, side =2, line = 3, cex = cex.ax)
  if(!rmv.lgnd){
    legend("bottomleft", "Tip signal", lwd = ln, col = "firebrick3")
  }
  
  par(mar = margn.bottom)
  
  plot(xs, tip.fluo.ts[, 2], type = "o", pch = 19, cex = 0., col = "royalblue3", lwd = ln, xlab = "", ylab ="", xaxt = "n", yaxt = "n", ylim = range(c(tip.fluo.ts[, 2:3]), na.rm = TRUE))
  lines(xs, tip.fluo.ts[, 3], type = "o", pch = 19, cex = 0., col = "darkorchid2", lwd = ln)
  axis(2, cex = cex.ax, cex.axis = cex.ax)
  mtext(y.lab, side =2, line = 3, cex = cex.ax, cex.axis = cex.ax)
  axis(1, cex = cex.ax, cex.axis = cex.ax)
  mtext(paste("Time (", time.unit, ")", sep = ""), side = 1, line = 3, cex = cex.ax)
  
  if(!rmv.lgnd){
    legend("bottomleft", c("Center signal", "Edge signal"), lwd = ln, col = c("royalblue3", "darkorchid2"))
  }
  
  
  return(NULL)
  
}

PlotKymoAndFluoWithDygraph1 <- function(kymo.align, tip.fluo.ts, roi.ini.pxs, avg.width = 5, 
  max.tip = NULL, px.sz = 1, px.unit = "pixel", 
  time.step = 1, time.unit = "pixel", brks = NULL, 
  red.bias = 1, cex.ax = 1.15, ln = 2, y.lab = "Fluorescence (AU)", rmv.lgnd = FALSE) {
  
  xs <- (1:dim(kymo.align)[1] - 1) * time.step
  
  # First, creating a data.frame with first column being the x-axis
  dygraph.df1 <- as.data.frame(cbind(x.ax = xs, 
    y.ax = tip.fluo.ts[, 1]))
  
  # Second, creating the plots
  dygraph(dygraph.df1, ylab = y.lab, group = "fluorescence") %>%
    dySeries("y.ax", label = "Tip Signal", color = "red") %>%
    dyAxis("x", rangePad = 10) %>%
    dyAxis("y", rangePad = 10) %>%
    dyLegend(show = "always") %>%
    dyRangeSelector() %>%
    dyOptions(drawGrid = FALSE)
}

PlotKymoAndFluoWithDygraph2 <- function(kymo.align, tip.fluo.ts, roi.ini.pxs, avg.width = 5, 
  max.tip = NULL, px.sz = 1, px.unit = "pixel", 
  time.step = 1, time.unit = "pixel", brks = NULL, 
  red.bias = 1, cex.ax = 1.15, ln = 2, y.lab = "Fluorescence (AU)", rmv.lgnd = FALSE) {
  
  xs <- (1:dim(kymo.align)[1] - 1) * time.step
  
  # First, creating a data.frame with first column being the x-axis
  dygraph.df2 <- as.data.frame(cbind(x.ax = xs, 
    y.ax1 = tip.fluo.ts[, 2],
    y.ax2 = tip.fluo.ts[, 3]))
  
  # Second, creating the plots
  dygraph(dygraph.df2, ylab = y.lab, group = "fluorescence") %>%
    dySeries("y.ax1", label = "Center Signal", color = "blue") %>%
    dySeries("y.ax2", label = "Edge Signal", color = "magenta") %>%
    dyAxis("x", rangePad = 10) %>%
    dyAxis("y", rangePad = 10) %>%
    dyOptions(drawGrid = FALSE)
}

#-------------------------------------------------------------------------------
PlotFilteredKymo <- function(kymo, kymo.filt, max.tip, px.sz, px.unit, 
  time.step, time.unit, brks = NULL, red.bias.raw = 1, 
  red.bias.filt = 3, rmv.lgnd = FALSE, cex.ax = 1.15){ 
  
  xs <- (1:dim(kymo)[1] - 1) * time.step
  ys <- (1:dim(kymo)[2] - 1) * px.sz
  
  if(is.null(max.tip)){
    max.tip.pos <- max(ys)
  } else{
    max.tip.pos <- max.tip * px.sz
  }
  
  if(!is.null(brks)){
    if(length(brks) == 2){
      brks <- seq(from = min(brks), to = max(brks), length.out = 257)
    } else{
      brks <- seq(from = min(c(kymo)), to = max(c(kymo)), length.out = 257)
    }
  } else{
    brks <- seq(from = min(c(kymo)), to = max(c(kymo)), length.out = 257)
  }
  
  pal.raw <- colorRampPalette(tim.colors(256), bias = red.bias.raw)
  
  margn.top <- c(0,6,5,8)
  margn.bottom <- c(6,6,0,8)
  par(mfrow = c(2, 1))
  
  par(mar = margn.top)
  
  image(x = xs, y = ys, 
    z = kymo, 
    col = pal.raw(256), #fields::tim.colors
    breaks = brks, 
    xlab = "",
    ylab ="",
    xaxt = "n", yaxt = "n",
    ylim = c(0, max.tip.pos))
  
  axis(3, cex.axis = cex.ax, cex.axis = cex.ax)
  axis(2, cex.axis = cex.ax, cex.axis = cex.ax)
  mtext(paste("Time (", time.unit, ")", sep = ""), side = 3, line = 3, cex = cex.ax)
  mtext(paste("Length from base (", px.unit, ")", sep = ""), side =2, line = 3, cex = cex.ax)
  
  
  if(!rmv.lgnd){
    legend("topleft", "Raw kymograph", text.col = "white", bty = "n")
  }
  
  pal.filt <- colorRampPalette(tim.colors(256), bias = red.bias.filt)
  
  par(mar = margn.bottom)
  
  brks <- seq(from = min(c(kymo.filt)), to = max(c(kymo.filt)), length.out = 257)
  
  image(x = xs, y = ys, 
    z = kymo.filt, 
    col = pal.filt(256), #fields::tim.colors
    breaks = brks, 
    xlab = "",
    ylab ="",
    xaxt = "n", yaxt = "n",
    ylim = c(0, max.tip.pos))
  
  axis(1, cex.axis = cex.ax, cex.axis = cex.ax)
  axis(2, cex.axis = cex.ax, cex.axis = cex.ax)
  mtext(paste("Time (", time.unit, ")", sep = ""), side = 1, line = 3, cex = cex.ax)
  mtext(paste("Length from base (", px.unit, ")", sep = ""), side =2, line = 3, cex = cex.ax)
  
  if(!rmv.lgnd){
    legend("topleft", "Filtered kymograph", text.col = "white", bty = "n")
  }
  
  return(NULL)
  
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
PlotTrimKymoAndFluo <- function(kymo.trim, trim.length, trim.time, tip.filt.fluo.ts, 
  roi.ini.pxs, avg.width, max.tip, px.sz, px.unit, time.step, time.unit, rmv.lgnd, brks = NULL, 
  red.bias = 1, y.lab = "Fluorescence (AU)", cex.ax = 1.15, ln = 1.5){
  
  roi.ini.pos <- roi.ini.pxs * px.sz
  avg.width.pos <- (avg.width - 1) * px.sz
  
  xs <- (trim.time[1]:trim.time[2] - 1) * time.step
  ys <- (trim.length[1]:trim.length[2] - 1) * px.sz
  
  if(is.null(max.tip)){
    max.tip.pos <- max(ys)
  } else{
    max.tip.pos <- max.tip * px.sz
  }
  max.tip.pos <- min(max.tip.pos, max(xs))
  
  if(!is.null(brks)){
    if(length(brks) == 2){
      brks <- seq(from = min(brks), to = max(brks), length.out = 257)
    } else{
      brks <- seq(from = min(c(kymo.trim)), to = max(c(kymo.trim)), length.out = 257)
    }
  } else{
    brks <- seq(from = min(c(kymo.trim)), to = max(c(kymo.trim)), length.out = 257)
  }
  
  pal <- colorRampPalette(tim.colors(256), bias = red.bias)
  
  margn.top <- c(0,6,5,8)
  margn.middle <- c(0,6,0,8)
  margn.bottom <- c(6,6,0,8)
  par(mfrow = c(3, 1))
  
  par(mar = margn.top)
  
  image(x = xs, y = ys, 
    z = kymo.trim, 
    col = pal(256), #fields::tim.colors
    breaks = brks, 
    xlab = "",
    ylab ="",
    xaxt = "n", yaxt = "n",
    ylim = c(0, max.tip.pos))
  
  axis(3, cex.axis = cex.ax, cex.axis = cex.ax)
  axis(2, cex.axis = cex.ax, cex.axis = cex.ax)
  mtext(paste("Time (", time.unit, ")", sep = ""), side = 3, line = 3, cex = cex.ax)
  mtext(paste("Length from tip (", px.unit, ")", sep = ""), side =2, line = 3, cex = cex.ax)
  
  rect(xleft = min(xs) - time.step * 2, 
    ybottom = roi.ini.pos[1],
    xright = max(xs) + time.step * 2,
    ytop = roi.ini.pos[1] + avg.width.pos,
    density = 10,
    col = adjustcolor("black", alpha.f = 0.5)
  )
  
  text(x = min(xs) + (max(xs) - min(xs)) / 2, 
    y = (roi.ini.pos[1] + avg.width.pos) + 4 * avg.width.pos,
    labels = "Tip",
    cex = cex.ax)
  
  rect(xleft = min(xs) - time.step * 2, 
    ybottom = roi.ini.pos[2],
    xright = max(xs) + time.step * 2,
    ytop = roi.ini.pos[2] + avg.width.pos,
    density = 10,
    col = adjustcolor("black", alpha.f = 0.5)
  )
  
  text(x = min(xs) + (max(xs) - min(xs)) / 2, 
    y = (roi.ini.pos[2] + avg.width.pos) + 4 * avg.width.pos,
    labels = "Center",
    cex = cex.ax)
  
  rect(xleft = min(xs) - time.step * 2, 
    ybottom = roi.ini.pos[3],
    xright = max(xs) + time.step * 2,
    ytop = roi.ini.pos[3] + avg.width.pos,
    density = 10,
    col = adjustcolor("black", alpha.f = 0.5)
  )
  
  text(x = min(xs) + (max(xs) - min(xs)) / 2, 
    y = (roi.ini.pos[3] + avg.width.pos) + 4 * avg.width.pos,
    labels = "Edge",
    cex = cex.ax)
  
  if(!rmv.lgnd){
    legend("topleft", "Tip aligned filtered kymograph", text.col = "white", bty = "n")
  }
  
  par(mar = margn.middle)
  
  plot(xs, tip.filt.fluo.ts[, 1], type = "o", pch = 19, cex = 0., col = "firebrick3", lwd = ln, xlab = "", ylab ="", xaxt = "n", yaxt = "n")
  axis(2, cex = cex.ax, cex.axis = cex.ax)
  mtext(y.lab, side =2, line = 3, cex = cex.ax)
  if(!rmv.lgnd){
    legend("bottomleft", "Tip filtered signal", lwd = ln, col = "firebrick3")
  }
  
  par(mar = margn.bottom)
  
  plot(xs, tip.filt.fluo.ts[, 2], type = "o", pch = 19, cex = 0., col = "royalblue3", lwd = ln, xlab = "", ylab ="", xaxt = "n", yaxt = "n", ylim = range(c(tip.filt.fluo.ts[, 2:3]), na.rm = TRUE))
  lines(xs, tip.filt.fluo.ts[, 3], type = "o", pch = 19, cex = 0., col = "darkorchid2", lwd = ln)
  axis(2, cex = cex.ax, cex.axis = cex.ax)
  mtext(y.lab, side =2, line = 3, cex = cex.ax, cex.axis = cex.ax)
  axis(1, cex = cex.ax, cex.axis = cex.ax)
  mtext(paste("Time (", time.unit, ")", sep = ""), side = 1, line = 3, cex = cex.ax)
  
  if(!rmv.lgnd){
    legend("bottomleft", c("Center signal", "Edge signal"), lwd = ln, col = c("royalblue3", "darkorchid2"))
  }
  
  
  return(NULL)
}

#-------------------------------------------------------------------------------
### AUXILIARY FUNCTIONS---------------------------------------------------------
#-------------------------------------------------------------------------------
# FUNCTION TO READ KYMOGRAPH DEPENDING ON THE EXTENSION
SaveTable <- function(tbl, fl.nm = "output", extension = "txt", sep = " "){
  if(extension == "txt"){
    if(!(sep %in% c("\t", " "))){
      message("Something wrong with the separator!")
    }
    dat <- write.table(tbl, paste(fl.nm, extension, sep = "."), sep = sep)
  } else if(extension == "csv"){
    if(!(sep %in% c(",", ";"))){
      message("Something wrong with the separator!")
    }
    dat <- write.csv(tbl, paste(fl.nm, extension, sep = "."), sep = sep)
  }else{
    message("Something wrong with the file extension!")
  }
}
#-------
# FUNCTION TO READ KYMOGRAPH DEPENDING ON THE EXTENSION
ReadKymo <- function(fl.nm){
  # Find out extension
  fl.split <- unlist(strsplit(fl.nm, split = "[.]"))
  extension <- fl.split[length(fl.split)]
  
  
  # Find out separator
  fl.con <- file(fl.nm, "r")
  fl.ln <- readLines(fl.con, n = 1)
  close(fl.con)
  seps <- c("\t", " ", ",", ";")
  seps.res <- unlist(lapply(seps, function(sep, fl.ln) length(unlist(strsplit(fl.ln, split = sep))), fl.ln))
  #unlist(lapply(seps, function(sep) length(unlist(strsplit(fl.ln, split = sep))))) # MORE ECONOMICAL
  
  if(any(seps.res != 1)){
    sep <- seps[which.max(seps.res)]
  } else{
    message("Something wrong with the data separator!")
    sep <- NULL
  }
  
  if(extension == "txt"){
    dat <- read.table(fl.nm, sep = sep, header = FALSE)
  } else if(extension == "csv"){
    dat <- read.csv(fl.nm, sep = sep, header = FALSE)
  } else{
    message("Something wrong with the file extension!")
    dat <- NULL
  }
  
  
  dat <- as.matrix(apply(dat, 2, as.numeric))
  return(dat)
}
#-------------------------------------------------------------------------------
RemoveTipOutliers <- function(tip.loc, dg = 2, spn = 0.2, tol = 3.5, px.tol = 5){
  tip.loc.diff <- diff(tip.loc)
  px_score <- which(tip.loc.diff >= px.tol)
  if(length(px_score) > 0){
    input.diff <- which(tip.loc.diff < px.tol)[findInterval(px_score,  which(tip.loc.diff < px.tol))]
    tip.loc[px_score + 1] <- tip.loc[input.diff + 1]
  }
  
  tip.loc.coarse <- SmoothWithLoess(tip.loc, spn, dg)
  MZ_score <- abs((0.6745 * ((tip.loc.coarse - tip.loc) - median(tip.loc.coarse - tip.loc))) / mad(tip.loc.coarse - tip.loc))
  if(any(MZ_score > tol)){
    tip.loc[which(MZ_score > tol)] <- tip.loc.coarse[which(MZ_score > tol)]
  }
  return(tip.loc)
}
#-------------------------------------------------------------------------------
FilterObs <- function(obs, mad.tol = 3.5, qntl = 0.95){
  # ALWAYS calculated the threshold ABOVE (average must be positive!)
  mad.cutoff <- (abs(median(obs, na.rm = T)) + mad.tol * mad(obs, na.rm = T))
  indxs <- which(obs > mad.cutoff)
  
  if(length(indxs) != 0){
    return(indxs)
  } else{
    qntl.cutoff <- median(obs, na.rm = T) + 3 * IQR(obs, na.rm = T)
    indxs <- which(obs > qntl.cutoff)
    #return(indxs) # There can't be a 'return(indxs)' here (it'll never pass through the next ifs)
  }
  if(length(indxs) == 0){
    qntl.cutoff <- median(obs, na.rm = T) + 1.5 * IQR(obs, na.rm = T)
    indxs <- which(obs > qntl.cutoff)
    return(indxs)
  }
  if(length(indxs) == 0){
    indxs <- which(obs > quantile(obs, probs = qntl, na.rm = T))
    return(indxs)
  }
  return(indxs)
}
#-------------------------------------------------------------------------------
SmoothWithLoess <- function(y, spn = 0.03, degree = 2){
  x <- 1:length(y)
  smth <- loess(y ~ x, span = spn, control = loess.control(surface = "direct"), 
    family = "gaussian", degree = degree) # span controls smoothing
  #to increase accuracy, pass further arguments to loess.control,
  #statistics = "exact", trace.hat = "exact", iterations = 10
  return(smth$fitted)
}
#-------------------------------------------------------------------------------
AlignByTip <- function(tip.loc, imaj){
  exceed <- which(tip.loc > dim(imaj)[2])
  if(length(exceed) > 0){
    row.limit <- min(which(tip.loc > dim(imaj)[2])) - 1
    message(paste("Kymograph frames", (row.limit + 1), "to", dim(imaj)[2], 
      "excluded"))
  } else{row.limit <- length(tip.loc)}
  
  imaj.nw <- t(sapply(1:row.limit, GetNewFracRow, imaj = imaj,
    tip.loc = tip.loc))
  
  return(imaj.nw[, dim(imaj.nw)[2]:1])
}
#-------------------------------------------------------------------------------
GetNewFracRow <- function(row.n, imaj, tip.loc){
  tip.round <- floor(tip.loc[row.n])
  tip.frac <- (tip.loc%%1)[row.n]
  #if(length(2:tip.round) != 0){}
  ind.vec <- 2:tip.round
  zvals.vec <- imaj[row.n, ]
  new.row <- sapply(ind.vec, CalcNewFracVal, frac = tip.frac, vec = zvals.vec)
  val <- c(rep(0, (dim(imaj)[2]) - length(new.row)), new.row)
  
  return(val)
}
#-------------------------------------------------------------------------------
CalcNewFracVal <- function(i, frac, vec){
  return(weighted.mean(c(vec[i], vec[i - 1]), c(frac, 1 - frac)))
}
#-------------------------------------------------------------------------------
GetFluoTimeSeries <- function(kymo, ini.ind, avg.width, use.median = TRUE){
  if(use.median){
    return(apply(kymo[, ini.ind:(ini.ind + avg.width - 1)], 1, median))
  } else{
    return(apply(kymo[, ini.ind:(ini.ind + avg.width - 1)], 1, mean))
  }
}
#-----------------------------------------------------------------------------
GetFreqOrPerBands <- function(sampling.freq, band.n = 11, per = FALSE){
  
  freq.band.min <- sampling.freq / (2^(1:band.n + 1)) #2^(1:band.n))
  freq.band.max <- sampling.freq / (2^(1:(band.n))) #(2^(0:(band.n - 1)))
  approx.freq.band.max <- sampling.freq / (2^(1:band.n + 1)) ##2^(1:band.n))
  
  out.bands <- cbind.data.frame(freq.band.max, 
    freq.band.min, 
    approx.freq.band.max)
  
  if(per){
    out.bands <- 1 / out.bands
    colnames(out.bands) <- c("per.band.min", 
      "per.band.max", 
      "approx.per.band.min")
  }
  
  return(out.bands)
}
#-----------------------------------------------------------------------------
GetSmoothVec <- function(sampling.freq, series.length, low.per = 16, 
  high.per = 256){
  #in seconds!
  per.bands <- GetFreqOrPerBands(sampling.freq, 
    band.n = 11, 
    per = TRUE)
  
  per_band_max <- per.bands$per.band.max
  per_band_min <- per.bands$per.band.min
  
  vec.length <- sum(per_band_max <= series.length)
  smth.vec <- rep(TRUE, vec.length)
  
  smth.vec[(per_band_max[1:vec.length] <= low.per)] <- FALSE
  smth.vec[(per_band_min[1:vec.length] >= high.per)] <- FALSE
  
  return(smth.vec)
}

# 'Batch Mode' function: -----------------------------------------------------
RunBatch <- function(df, df.name, parameters = all.par.tbl(), extension = extension.out(), sep = sep.out()) { #, decision.Tab6 = decision.Tab6()
  
  kymo <- df
  # Getting the value of all parameters: 
  fl.nm <- as.numeric(parameters[which(parameters == "fl.nm"), "value"])
  px.sz <- as.numeric(parameters[which(parameters == "px.sz"), "value"]) 
  px.unit <- as.character(parameters[which(parameters == "px.unit"), "value"]) 
  
  time.step <- as.numeric(parameters[which(parameters == "time.step"), "value"]) 
  time.unit <- as.character(parameters[which(parameters == "time.unit"), "value"]) 
  red.bias.tab2 <- as.numeric(parameters[which(parameters == "red.bias.tab2"), "value"])
  
  use.smooth <- as.character(parameters[which(parameters == "use.smooth"), "value"]) 
  if(use.smooth == "TRUE"){ use.smooth <- TRUE }else{ use.smooth <- FALSE}
  kymo.span.n <- as.numeric(parameters[which(parameters == "kymo.span.n"), "value"])
  kymo.loess.dg <- as.numeric(parameters[which(parameters == "kymo.loess.dg"), "value"]) 
  
  n.pts <- as.numeric(parameters[which(parameters == "n.pts"), "value"]) 
  mad.tol <- as.numeric(parameters[which(parameters == "mad.tol"), "value"])
  fix.slope <- as.character(parameters[which(parameters == "fix.slope"), "value"]) 
  if(fix.slope == "TRUE"){ fix.slope <- TRUE }else{ fix.slope <- FALSE}
  
  qntl <- as.numeric(parameters[which(parameters == "qntl"), "value"]) 
  fit.right <- as.character(parameters[which(parameters == "fit.right"), "value"])
  if(fit.right == "TRUE"){ fit.right <- TRUE }else{ fit.right <- FALSE}
  
  
  
  fluo.frac <- as.numeric(parameters[which(parameters == "fluo.frac"), "value"])
  min.chunk.size <- as.numeric(parameters[which(parameters == "min.chunk.size"), "value"])
  coarse.kymo.loess.dg <- as.numeric(parameters[which(parameters == "coarse.kymo.loess.dg"), "value"])
  coarse.kymo.span.n <- as.numeric(parameters[which(parameters == "coarse.kymo.span.n"), "value"])
  
  use.coarse <- as.character(parameters[which(parameters == "use.coarse"), "value"])
  if(use.coarse == "TRUE"){ use.coarse <- TRUE }else{ use.coarse <- FALSE}
  rmv.min <- as.character(parameters[which(parameters == "rmv.min"), "value"])
  if(rmv.min == "TRUE"){ rmv.min <- TRUE }else{ rmv.min <- FALSE}
  algorithm <- as.character(parameters[which(parameters == "algorithm"), "value"])
  
  
  
  tip.estimate <- as.character(parameters[which(parameters == "tip.est"), "value"])
  
  t.slice.ind <- as.numeric(parameters[which(parameters == "t.slice.ind"), "value"])
  red.bias.tab3 <- as.numeric(parameters[which(parameters == "red.bias.tab3"), "value"]) 
  restrict.x <- as.character(parameters[which(parameters == "restrict.x"), "value"])
  if(restrict.x == "TRUE"){ restrict.x <- TRUE }else{ restrict.x <- FALSE}
  
  tip.span.n <- as.numeric(parameters[which(parameters == "tip.span.n"), "value"])
  tip.span.dg <- as.numeric(parameters[which(parameters == "tip.span.dg"), "value"])
  
  
  out.rm.algorithm <- as.character(parameters[which(parameters == "out.rm.algorithm"), "value"])
  
  
  rm.tip.out <- as.character(parameters[which(parameters == "rm.tip.out"), "value"])
  if(rm.tip.out == "TRUE"){ rm.tip.out <- TRUE }else{ rm.tip.out <- FALSE}
  
  rm.tip.out.dg <- as.numeric(parameters[which(parameters == "rm.tip.out.dg"), "value"]) 
  rm.tip.out.spn <- as.numeric(parameters[which(parameters == "rm.tip.out.spn"), "value"])
  tol <- as.numeric(parameters[which(parameters == "tol"), "value"])
  
  px.tol <- as.numeric(parameters[which(parameters == "px.tol"), "value"])
  rmv.lgnd <- as.character(parameters[which(parameters == "rmv.lgnd"), "value"])
  if(rmv.lgnd == "TRUE"){ rmv.lgnd <- TRUE }else{ rmv.lgnd <- FALSE}
  red.bias.tab4 <- as.numeric(parameters[which(parameters == "red.bias.tab4"), "value"])
  
  tip.mrgn <- as.numeric(parameters[which(parameters == "tip.mrgn"), "value"])
  manual.ROI <- as.character(parameters[which(parameters == "manual.ROI"), "value"])
  if(manual.ROI == "TRUE"){ manual.ROI <- TRUE }else{ manual.ROI <- FALSE}
  roi.px <- as.numeric(parameters[which(parameters == "roi.px"), "value"]) 
  
  avg.width <- as.numeric(parameters[which(parameters == "avg.width"), "value"])
  use.median.for.avg <- as.character(parameters[which(parameters == "use.median.for.avg"), "value"])
  if(use.median.for.avg == "TRUE"){ use.median.for.avg <- TRUE }else{ use.median.for.avg <- FALSE}
  red.bias.tab5 <- as.numeric(parameters[which(parameters == "red.bias.tab5"), "value"]) 
  
  y.lab.tab5 <- as.character(parameters[which(parameters == "y.lab"), "value"])
  
  
              # per.min <- as.numeric(parameters[which(parameters == "per.min"), "value"]) 
              # per.max <- as.numeric(parameters[which(parameters == "per.max"), "value"])
              # 
              # trim.time <- c()
              # trim.length <- c()
              # 
              # trim.time[1] <- as.numeric(parameters[which(parameters == "trim.time.1"), "value"])
              # trim.time[2] <- as.numeric(parameters[which(parameters == "trim.time.2"), "value"])
              # trim.length[1] <- as.numeric(parameters[which(parameters == "trim.length.1"), "value"])
              # 
              # trim.length[2] <- as.numeric(parameters[which(parameters == "trim.length.2"), "value"])
              # tip.mrgn.filt <- as.numeric(parameters[which(parameters == "tip.mrgn.filt"), "value"])
              # manual.ROI.tab6 <- as.character(parameters[which(parameters == "manual.ROI.tab6"), "value"]) 
              # if(manual.ROI.tab6 == "TRUE"){ manual.ROI.tab6 <- TRUE }else{ manual.ROI.tab6 <- FALSE}
              # 
              # roi.px.tab6 <- as.numeric(parameters[which(parameters == "roi.px.tab6"), "value"]) 
              # avg.width.tab6 <- as.numeric(parameters[which(parameters == "avg.width.tab6"), "value"])
              # use.median.for.avg.tab6 <- as.character(parameters[which(parameters == "use.median.for.avg.tab6"), "value"])
              # if(use.median.for.avg.tab6 == "TRUE"){ use.median.for.avg.tab6 <- TRUE }else{ use.median.for.avg.tab6 <- FALSE}
              # 
              # red.bias.tab6 <- as.numeric(parameters[which(parameters == "red.bias.tab6"), "value"])
              # red.bias.tab6.trim <- as.numeric(parameters[which(parameters == "red.bias.tab6.trim"), "value"])
              # y.lab.tab6 <- as.character(parameters[which(parameters == "y.lab.tab6"), "value"])
              # 
  
  # Running calls from Third Tab:
  kymo <- TreatKymo(kymo)
  
  tip.lst <- FindTip( kymo, #NEW!!!
    use.smooth, kymo.span.n, kymo.loess.dg, 
    n.pts , mad.tol , qntl, fit.right, fix.slope,     
    fluo.frac , min.chunk.size , 
    coarse.kymo.span.n , coarse.kymo.loess.dg ,
    use.coarse , rmv.min , algorithm) #TODO Change FindTip function
  
  kymo.smth <- tip.lst$kymo.smth 
  tip.tbl <-  tip.lst$tip.tbl
  tip.loc <-  tip.tbl[, which(colnames(tip.tbl) == tip.estimate )] 
  
  
  # Running calls from Fourth Tab:
  tip.loc.smth.lst <-  RefineTipLoc(tip.loc , tip.span.n , tip.span.dg , rm.tip.out , rm.tip.out.dg , 
    rm.tip.out.spn , tol , px.tol , out.rm.algorithm ) 
  
  tip.loc.out <-  tip.loc.smth.lst$tip.loc 
  tip.loc.smth <-  tip.loc.smth.lst$tip.loc.smth 
  out.indx <- tip.loc.smth.lst$out.indx 
  
  # Running calls from Fifth Tab:
  kymo.align <-  AlignByTip(tip.loc = tip.loc.smth  + tip.mrgn , imaj = kymo ) 
  
  tip.fluo.lst <- ExtractFluoTimeSeries(kymo.align , roi.px , avg.width , use.median = use.median.for.avg , max.tip = min(tip.loc.smth , na.rm = TRUE)) # Wait for Dani to finish this and put '()' in all of them
  
  
  # Running calls from Sixth Tab:
              # kymo.filt <- FilterKymo(kymo, time.step, low.per = per.min, high.per = per.max)
              # kymo.filt.align <- AlignByTip(tip.loc = tip.loc.smth + tip.mrgn.filt, imaj = kymo.filt)
              # kymo.trim <- kymo.filt.align[trim.time[1]:trim.time[2], trim.length[1]:trim.length[2]] 
              # tip.filt.fluo.lst <- ExtractFluoTimeSeries(kymo.trim, roi.px.tab6, avg.width.tab6, 
              #   use.median = use.median.for.avg.tab6, 
              #   max.tip = min(tip.loc.smth, na.rm = TRUE)) 
              # 
  
  # Running calls from 'Save' Tab:
  all.ts <- cbind("time" = 0:(dim(kymo)[1] - 1) * time.step,
    "tip.loc.raw" = tip.loc * px.sz,
    "tip.loc.smth" = tip.loc * px.sz,
    "growth" = c(NA, diff(tip.loc.smth * px.sz)/time.step),
    tip.fluo.lst[[2]]
    )
  
  # all.trim.fluo <- cbind("time" = ((trim.time[1]:trim.time[2]) - 1) * time.step,
  #   tip.filt.fluo.lst[[2]])
  # 
  
  # Creating the documents for saving (it should be in the wd where RunBatch() was called ) and getting their names (for the 'fs' argument in downloadHandler)
  
  pdf1_name <- file.path("TipDetection_", df.name, ".pdf", fsep = "")
  pdf(pdf1_name)
  PlotKymoWithAllTip(kymo.smth , tip.tbl , px.sz , px.unit , time.step , 
    time.unit , brks = NULL, red.bias = red.bias.tab3 , restrict.x = restrict.x )
  dev.off()
  
  pdf2_name <- file.path("TipLocation_", df.name, ".pdf", fsep = "")
  pdf(pdf2_name)
  PlotKymoWithTip(kymo.smth , tip.loc.raw = tip.loc , tip.loc.smth = tip.loc.smth , tip.loc.out = tip.loc.out , out.indx = out.indx ,
    px.sz , px.unit , time.step , time.unit , brks = NULL, red.bias = red.bias.tab4 , rmv.lgnd = rmv.lgnd )
  dev.off()
  
  pdf3_name <- file.path("FluorescenceSeries_", df.name, ".pdf", fsep = "")
  pdf(pdf3_name)
  PlotKymoAndFluo(kymo.align , tip.fluo.ts = tip.fluo.lst[[2]], roi.ini.pxs = tip.fluo.lst[[1]], avg.width , max.tip = max(tip.loc.smth , na.rm = TRUE), 
    px.sz , px.unit , time.step , time.unit , brks = NULL, red.bias = red.bias.tab5 , y.lab = y.lab.tab5 )
  dev.off()
  
                # pdf4_name <- file.path("FilterKymograph_", df.name, ".pdf", fsep = "")
                # pdf(pdf4_name)
                # if (decision.Tab6 == 0){  # & input$tab.save == 0
                #   PlotFilteredKymo(kymo, kymo.filt, max.tip = max(tip.loc.smth, na.rm = TRUE), px.sz, px.unit, 
                #     time.step, time.unit, brks = NULL, red.bias.raw = red.bias.tab4, 
                #     red.bias.filt = red.bias.tab6, rmv.lgnd.tab6, cex.ax = 1.15)
                #   PlotTrimKymoAndFluo(kymo.trim, trim.length, trim.time, tip.filt.fluo.ts = tip.filt.fluo.lst[[2]],
                #     roi.ini.pxs = tip.filt.fluo.lst[[1]], avg.width = avg.width.tab6, max.tip = max(tip.loc.smth, na.rm = TRUE), 
                #     px.sz, px.unit, time.step, time.unit, rmv.lgnd.tab6, brks = NULL, red.bias = red.bias.tab6.trim, y.lab = y.lab.tab6)
                # }else{
                #   wrng.txt <- paste("The Parameters per.min, per.max, trim.time.1, trim.time.2, 
                #     trim.length.1, trim.length.2, tip.mrgn.filt, manual.ROI.tab6, 
                #     roi.px.tab6, avg.width.tab6, use.median.for.avg.tab6, 
                #     red.bias.tab6, red.bias.tab6.trim and y.lab.tab6 present in 
                #     'Parameters.txt(.csv)' all pertain to the sixth tab ('Filter Kymograph') 
                #     and since it wasn't run, such paramenters do not have any reliable value")
                #   plot(NA, xlim=c(0,50), ylim=c(0,50), bty='n', xaxt='n', yaxt='n', xlab='', ylab='')
                #   text(20, 40, labels = wrng.txt, cex = 0.75)
                # }
                # dev.off()
  
  table1_name <- file.path("Tip_aligned_kymograph_", df.name, fsep = "")
  SaveTable(kymo.align, fl.nm = table1_name, extension, sep)
  
      # table2_name <- file.path("Tip_aligned_filtered_kymograph_", df.name, fsep = "")
      # SaveTable(kymo.filt.align, fl.nm = table2_name, extension, sep)
      # 
      # table3_name <- file.path("Tip_aligned_filtered_trimmed_kymograph_", df.name, fsep = "")
      # SaveTable(kymo.trim, fl.nm = table3_name, extension, sep)
  
  table4_name <- file.path("All_time_series_", df.name, fsep = "")
  SaveTable(all.ts, fl.nm = table4_name, extension, sep)
  
      # table5_name <- file.path("All_filtered_time_series_", df.name, fsep = "")
      # SaveTable(all.trim.fluo, fl.nm = table5_name, extension, sep)
  
  table6_name <- file.path("Parameters_", df.name, fsep = "")
  SaveTable(parameters, fl.nm = table6_name, extension, sep)
  
  
  all_names <- c(pdf1_name, pdf2_name, pdf3_name, file.path(table1_name, extension, fsep = "."), #pdf4_name, 
        #file.path(table2_name, extension, fsep = "."), file.path(table3_name, extension, fsep = "."),
    file.path(table4_name, extension, fsep = "."), #file.path(table5_name, extension, fsep = "."),
    file.path(table6_name, extension, fsep = "."))
  
  return(all_names)
}

# 'Rotations' functions: ----

clock.rotation <- function(x) {
  t(apply(x, 2, rev))
}

counter.rotation <- function(x) {
  apply(t(x), 2, rev)
}

invert.column <- function(m) {
  m <- m[, dim(m)[2]:1]
}

invert.row <- function(m) {
  m <- m[dim(m)[1]:1, ]
}