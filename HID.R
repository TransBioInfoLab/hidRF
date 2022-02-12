library("randomForestSRC")
#library("stats")
#library("spatstat")
#library("RColorBrewer")
#############################################################################################
### MDII based weight method to update subRF                                              ###
#############################################################################################
wt.itr <- function(formula, data, ntime = 100,
                   pmd.initial,
                   mtry = round(ncol(data)/3), nodesize = 3,
                   wt = function(btpmd,digpmd){log(1/btpmd/digpmd)}, iteration = 3){
  pmd <- pmd.initial
  rm(pmd.initial)
  gc()
  bd.kink <- function(pmd){
    o.max <- nrow(pmd)
    sapply(seq(o.max),function(x)
      which.max(sort(pmd[x,-x],decreasing = FALSE)-
                  smooth.spline(seq(o.max-1),sort(pmd[x,-x],decreasing = FALSE))$y))
  }
  lapply(1:iteration,function(i){
    sim.intcopy <- pmd
    k <- bd.kink(sim.intcopy)
    btw <- mapply(function(x,k)
      mean(sort(sim.intcopy[x,-x], decreasing = FALSE)[1:k], na.rm=TRUE)
      , seq(nrow(sim.intcopy)), k)
    digw <- diag(pmd)
    
    if (is.null(names(diag))){id <- 1:ncol(pmd)} else{
      yvar_name=colnames(data)[which(sapply(1:ncol(data),function(i) {
        grepl(colnames(data)[i],as.character(formula)[[2]])})==TRUE)]
      id <- match(colnames(data[,which(colnames(data)!=yvar_name)]),names(digw))}
    btpmd <- btw[id]; digpmd <- digw[id]
    w <- wt(btpmd,digpmd)
    
    if (iteration!=1){
      sim.out <- rfsrc(formula = formula, data, 
                       ntime = ntime,
                       predictorWt = w, statistics = F,
                       mtry = mtry, nodesize = nodesize, seed = -1)
      gc()
      pmd <<- find.interaction(sim.out, method="maxsubtree", sorted = F, verbose = F) 
      rm(sim.out)
      gc()
      }
    out <- cbind(pmd, digw, btw, w)
    colnames(out) <- c(paste("pmd", rownames(pmd)),"digpmd", "btpmd", "wt")
    data.frame(out)
  })
}

##################
### Grow subRF ###
##################
subRF <- function(formula, data, 
                  ntime = 100,
                  w0, 
                  subvars = function(data){ceiling((ncol(data)-1)/5)}, 
                  n.RF = 5, wtRF = T, 
                  subRF = subRF, 
                  forest.turn = FALSE){
  var.columns <- (1:ncol(data))[which(colnames(data) %in% names(w0))]
  x <- which(names(w0) %in% colnames(data) )
  w0 <- as.vector(w0)
  w0[which(w0 <= 0)] <- 0
  size <- subvars(data)
  lapply(1:n.RF, function(i){
    if (size == length(var.columns)) {var.pt <- var.columns} else {
      set.seed(i)
      var.pt <- sample(x = var.columns,size = size,replace = FALSE, prob = w0[x])}
    subdata <- data[,c(var.pt,setdiff(1:ncol(data),var.columns))]
    wts1 <- w0[var.pt]
    if (wtRF == T){
      sim.out <- rfsrc(formula, data = subdata, ntime = ntime,
                       predictorWt = wts1, seed = -1)
      gc()} else {
        if (forest.turn) {
        turn.o <- tune(formula, data = subdata)
        sim.out <- rfsrc(formula,data = subdata,
                         ntime = ntime, 
                         mtry = turn.o$rf$mtry, nodesize = turn.o$rf$nodesize, 
                         seed = -1)
        rm(trun.o)
      } else {
        sim.out <- rfsrc(formula,data = subdata,
                         ntime = ntime, seed = -1)
      }
        gc()
        }
    pmd <- find.interaction(sim.out, method = "maxsubtree", sorted = F, verbose = F)
    err.sub <- sim.out$err.rate 
    mtry <- sim.out$mtry
    nodesize <- sim.out$nodesize
    list(pmd = pmd, mtry = mtry, nodesize = nodesize, err.sub = err.sub)
  })
}

##########################
### Variable Selection ###
##########################
var.sel <- function(formula, data,
                    ntime = 100,
                    subvars = function(data){ceiling((ncol(data)-1)/5)},
                    w.initial = "vimp", max.var = 50,
                    wt = function(btpmd,digpmd){log(1/btpmd/digpmd)},
                    itrSub = 5, wtSub = T,
                    itrWt = 1, 
                    subQtl = function(btpmd,digpmd){which(btpmd < quantile(btpmd, probs = 0.1))}, 
                    verbose = TRUE,
                    obj = NULL,
                    forest.turn = FALSE){
  if (w.initial == "vimp"){
    w.initial <- function(formula, data, obj = NULL){
      if (is.null(obj)){
        if (forest.turn) {
       turn.o <- tune(formula, data)
        sim.out <- rfsrc(formula, data, ntime = ntime,
                         importance = "random",
                         mtry = turn.o$rf$mtry, nodesize = turn.o$rf$nodesize, 
                         seed = -1)
        rm(turn.o)
        } else {
          sim.out <- rfsrc(formula, data, ntime = ntime,
                           importance = "random", seed = -1) 
        }
        gc()
}
      else {sim.out <- obj}
      if (sim.out$family == "class"){
        return(list(sim.out$err.rate, pmax(sim.out$importance[,1],0)))
      } 
      if (sim.out$family == "regr+"){
        return(list(sum(get.mv.error(sim.out)), 
                    rowSums(pmax(get.mv.vimp(sim.out, standardize=TRUE),0))))
      } 
      if (sim.out$family == "surv"){
        return(list(na.omit(sim.out$err.rate)[[1]], pmax(sim.out$importance,0)))
      } 
      if (sim.out$family == "regr") {
        return(list(round((1-na.omit(sim.out$err.rate)/var(sim.out$yvar,na.rm = TRUE)),2),
             pmax(sim.out$importance,0)))
      }}
    w0 <- w.initial(formula, data, obj = obj)[[2]]
  } else if (w.initial == "md"){
    w.initial <- function(formula, data, obj = NULL){
      if (is.null(obj)){
        if (forest.turn) {
        turn.o <- tune(formula, data)
        sim.out <- rfsrc(formula,data,
                         ntime = ntime,
                         importance = FALSE,
                         mtry = turn.o$rf$mtry, nodesize = turn.o$rf$nodesize, 
                         seed = -1)
        rm(turn.o)
        } else {
          sim.out <- rfsrc(formula,data,
                           ntime = ntime,
                           importance = FALSE, seed = -1) 
        }
        gc()
        }
      else {sim.out <- obj}
      1/(max.subtree(sim.out, sub.order = FALSE, max.order = 1)$order+0.001)}
    w0 <- w.initial(formula, data, obj = obj)
  } else { w0 <- w.initial}
  err.init <- w.initial(formula,data,obj = obj)[[1]]
  rm(obj)
  gc()
  
  subchain0 <- subRF(formula, data, 
                     ntime = ntime,
                     w0 = w0, 
                     subvars = subvars, n.RF = itrSub, 
                     wtRF = wtSub, 
                     forest.turn = forest.turn)
  gc()
  err.sub <- unlist(lapply(1:itrSub, function(x){subchain0[[x]]$err.sub}))
  chainwts <- lapply(1:itrSub, function(i){
    if (verbose == TRUE) {cat("------", i,"/",itrSub, "------ \n")}
    if (itrWt == 1) {
      wt.itr(formula = formula, ntime = ntime,
             pmd.initial = subchain0[[i]]$pmd,
             wt = wt, iteration = itrWt)
    } else{
    wt.itr(formula = formula,data = subchain0[[i]]$subdata, 
           ntime = ntime,
           pmd.initial = subchain0[[i]]$pmd,
           mtry = subchain0[[i]]$mtry,
           nodesize = subchain0[[i]]$nodesize,
           wt = wt, iteration = itrWt)}
  })
  rm(subchain0)
  gc()
  var.sl.list <- lapply(1:itrSub, function(i){
    names <- rownames(chainwts[[i]][[itrWt]])
    id <- subQtl(chainwts[[i]][[itrWt]]$btpmd,chainwts[[i]][[itrWt]]$digpmd)
    list(sl.w = chainwts[[i]][[itrWt]]$wt[id], sl.var = names[id])
  })
  rm(chainwts)
  gc()
  #var.sl <- unique(unlist(lapply(1:itrSub, function(i){var.sl.list[[i]]$sl.var})))
  varl<-lapply(lapply(var.sl.list,`[`,2),`[[`,1)
  if (length(unique(unlist(varl)))<max.var){
    var.sl<-unique(unlist(varl))
  }
  else {
    freql<-table(unlist(varl))[order(table(unlist(varl)),decreasing=TRUE)]
    freql.stepsum<-sapply(unique(freql),function(i){
      sum(freql>=i)
    })
    freql.thresh<-unique(freql)[which(abs(freql.stepsum-max.var)==min(abs(freql.stepsum-max.var)))]
    var.sl<-names(freql)[which(freql>=freql.thresh)]
  }
  
  list(var.sl = var.sl,var.sl.list = var.sl.list, # itrwt = chainwts,
       err.init = err.init, 
       err.sub = err.sub)
}

##########################################################################
### Calculate PBII for all interactions, i.e, select 10 variables,    ###
### minimum order is defined by function "nterms". Take 3 for         ###
### example, all interactions will be (10 choose 3)+(10 choose 2)     ###
##########################################################################
pmdhvp <- function(obj, inter, verbose = T){
  pmd <- xx <- find.interaction(obj, sorted = F , method = "maxsubtree", verbose = F)
  rm(obj)
  gc()
  name <- rownames(pmd)
  pmd[lower.tri(xx)] <- (xx[lower.tri(xx)]+xx[upper.tri(xx)])/2
  pmd[upper.tri(xx,diag=T)] <- 0
  pmd = pmd + t(pmd) + diag(diag(xx))
  
  k <- length(inter)
  pmdvp <- matrix(unlist(lapply(1:k,function(i){
    if (((i%%100 == 0)|(i == k))&(verbose == T)){cat("\r","----",i,"/",k,100*(i/k),"%----")}
    
    (mean(apply(as.matrix(combn(inter[[i]],2)),2,function(x){
      pmd[x[1],x[2]] }))/max(1,length(inter[[i]])/5))  }) ),k ,1)
  rownames(pmdvp) <- unlist(lapply(1:k,function(i){ paste(name[inter[[i]]], collapse = "_") }))

  list(pmdvp = pmdvp, pmd = pmd)
}

##########################################################################
### This function is particularly used after selecting top 150 pmdvimp ###
### ranked interactions. It is to get the HIVIMP and normalized HIVIMP ###
### for the 150 interactions.                                          ###
##########################################################################

hvp <- function(joint, obj, importance = "anti",
                block = 1)
{ allobj <- obj
rm(obj)
gc()
x.var <- allobj$xvar.names

o <- do.call(rbind,lapply(1:nrow(joint),function(i){
    if (allobj$family == "class"){
      joinvimp <- vimp(allobj,x.var[joint[i,]],importance = importance,joint = T, seed = -1,block=block)$importance[1]
      hivimp <- (joinvimp - sum(allobj$importance[joint[i,],1]))
    }  else{
      if (allobj$family == "regr+") {
        joinvimp <- vimp(allobj,x.var[joint[i,]],importance = importance,joint = T, seed = -1,block=block)
        joinvimp <- rowSums(pmax(get.mv.vimp(joinvimp, standardize=TRUE),0))
      } else {
        joinvimp <- vimp(allobj,x.var[joint[i,]],importance = importance,joint = T, seed = -1,block=block)$importance
      }
        hivimp <- (joinvimp - sum(allobj$importance[joint[i,]]))
      }
    joinvimp <- joinvimp
    nmHIvimp <- hivimp/ncol(joint)
 # rm(obj,allobj)
  gc()
  out <- data.frame(PBII = abs(nmHIvimp),
                    jointVIMP = joinvimp)
  out <- out[,!is.na(out)]
  out <- as.data.frame(out)
  rownames(out) <- paste(x.var[joint[i,]],collapse ="_")
  out
}))
o
}

##############################################################
### Function to determine up to which order we will search ###
### for the interactions.                                  ###
##############################################################

nterms <- function(obj){
  if (nrow(obj$xvar) >= (100*ncol(obj$xvar))){
    stat.obj <- stat.split(obj)
    list <- rapply(stat.obj,classes="matrix",how="list",
                   f = function(x) x[,6,drop=FALSE])
    listidx <- seq(length(list))
    tree.depths <- sapply(listidx,function(x) max((unlist(list[[x]]))))
    Dt <- mean(tree.depths)
    All.mdp <- mean(max.subtree(obj, max.order = 0)$order)
    terms <- min(c(round(Dt-All.mdp),ncol(obj$xvar)))
    rm(obj,stat.obj)
    gc()
    list(tree.depth = Dt, min.depth = All.mdp, terms = terms)
  } else {
    x <- list(terms = min(c(floor(log2((nrow(obj$xvar))/obj$nodesize)),ncol(obj$xvar)))
    )
    rm(obj)
    gc() 
    x
    }
}

################################################
### Main function to get everything together ###
################################################

hidRF <- function(formula, data, 
                   ntime = 100,
                   terms = NULL,
                   importance = "anti", block = 1,
                   subRF = TRUE, subvars = NULL,
                   w.initial = "vimp",
                   nSubRF = 10, wtSub = TRUE,
                   itrWt = 1, max.var = 50,
                   nPBII = 50,
                   verbose = TRUE,
                   var.sel.only = FALSE,
                   MDIIonly = FALSE,
                   MDIIorder = TRUE,
                   forest.tune = FALSE){
  wt <- function(btpmd,digpmd){log(1/btpmd/digpmd)}
  nPBii = function(x){min(length(x), nPBII)}
  if (forest.tune){
    turn.o <- tune(formula, data)
    Allobj <-  rfsrc(formula, data, ntime = ntime,
                     importance = "random", statistics = TRUE,
                     mtry = turn.o$rf$mtry, nodesize = turn.o$rf$nodesize, 
                     seed = -1)
    rm(turn.o)
  } else {
    Allobj <-  rfsrc(formula, data, ntime = ntime,
                     importance = "random", statistics = TRUE, seed = -1)   
  }
    gc()
  if (is.null(terms)){
    terms <- nterms(Allobj)$terms
    if (verbose == TRUE) cat(paste("Automatically calculating upto ", terms, "-way interactions", sep =""))
  }
  subQtl <- NULL
  if (is.null(subQtl)){
    subQtl <- function(btpmd,digpmd){
      if (length(btpmd)<=140) {
        which(btpmd < sort(btpmd,decreasing=FALSE)
              [floor(0.5*which.max(sort(btpmd,decreasing = FALSE)-
                                     (((max(btpmd)-min(btpmd))/(length(btpmd)-1))*
                                        (seq(length(btpmd))-1)+min(btpmd)))+
                       0.5*which.max(sort(btpmd,decreasing=FALSE)-
                                       predict(lm(sort(btpmd,decreasing=FALSE)~seq(length(btpmd))))))])}
      else {
        which(btpmd < sort(btpmd,decreasing=FALSE)
              [which.max(sort(btpmd,decreasing = FALSE)-
                           (((max(btpmd)-min(btpmd))/(length(btpmd)-1))*
                              (seq(length(btpmd))-1)+min(btpmd)))])}
    }
  }
  
  if (is.null(subvars)){
    subvars <- function(data){
      yvar.num<-length(which(sapply(1:ncol(data),function(i) {
        grepl(colnames(data)[i],as.character(formula)[[2]])})==TRUE)
      )
      if((ncol(data)-yvar.num)<=50) {ncol(data)-yvar.num}
      else {ceiling((ncol(data)-yvar.num)/log(ncol(data)-yvar.num))}}
  }
  
  if (subRF == T){
    itrSub <- nSubRF
    if (is.null(itrSub)){itr=terms*ceiling(log10(ncol(data)-length(
      which(sapply(1:ncol(data),function(i) {
        grepl(colnames(data)[i],as.character(formula)[[2]])})==TRUE)
    )))} else{itr=itrSub}
    if (verbose == TRUE) {cat("Fitting", itr, "subRFs \n")}
    sub.x.obj <- var.sel(formula, data,
                         ntime = ntime,
                         subvars = subvars,
                         w.initial = w.initial, 
                         max.var=max.var, wt = wt,
                         itrSub = itr, wtSub = wtSub,
                         itrWt = itrWt,
                         subQtl = subQtl, verbose = verbose,
                         obj = Allobj,
                         forest.turn = forest.turn)
    if (var.sel.only) {return(list(subRF.obj = sub.x.obj, rf.obj = Allobj))} else {
    gc()
    rm(Allobj)
    gc()
    yvar_name=colnames(data)[which(sapply(1:ncol(data),function(i) {
      grepl(colnames(data)[i],as.character(formula)[[2]])})==TRUE)]
    data <- data[,c(sub.x.obj$var.sl, yvar_name)]
    
    var.sl <- sub.x.obj$var.sl
    var.sl.list <- sub.x.obj$var.sl.list
    itrwt <- sub.x.obj$itrwt
    err.init <- sub.x.obj$err.init
    err.sub <- sub.x.obj$err.sub
    rm(sub.x.obj)
    gc()
    if (forest.tune) {
    turn.o <- tune(formula, data, mtryStart = (ncol(data) - 1))
    allobj <- rfsrc(formula, data = data, ntime = ntime,
                    importance = importance, seed = -1,
                    mtry = (ncol(data)-1), nodesize = turn.o$rf$nodesize)
    rm(turn.o)
    } else {
      allobj <- rfsrc(formula, data = data, ntime = ntime,
                      importance = importance, seed = -1)  
    }
    gc() 
    
 
    cmbn <- lapply(2:terms,function(i){t(combn(1:length(allobj$xvar.names),i))})
    inter <- unlist(lapply(1:length(cmbn), function(i){
      lapply(1:nrow(cmbn[[i]]),function(j){ matrix(cmbn[[i]][j,],1,length(cmbn[[i]][j,])) })
    }), recursive = F )
    
      
      if (MDIIorder){
        n.vimp <- nPBii(inter)
        sub.itr.obj <- pmdhvp(allobj, inter, verbose = verbose)
        pmd <- sub.itr.obj$pmd
        if (length(inter) > n.vimp) {
          if (verbose == TRUE) {cat("\n Selected", n.vimp, "from", length(inter), "interactions based on Minimal Depth Interaction Importance (MDII)\n")}
          it.sel <- which(rank(sub.itr.obj$pmdvp) <= n.vimp)
          inter <- inter[it.sel]
        } else {it.sel <- 1:length(inter)}
        pmdvimp.long <- sub.itr.obj$pmdvp
      }else{
        pmd <- pmdvimp.long <- NULL
        n.vimp <- length(inter)
      }
      
    if (MDIIonly == TRUE) {
      return(list(terms = terms,
                  rf.obj = allobj, pmd = pmd,
                  var.sl = var.sl,
                  var.sl.list = var.sl.list,
                  itrwt = itrwt,
                  MDII = pmdvimp.long,
                  err.init = err.init,
                  err.sub = err.sub,
                  err.all = round(100*(1-allobj$err.rate[1000]),2)
                  ))
      
    } else{
    
  #  if (verbose == TRUE) {cat("\n Calculating", n.vimp, "High Order Interaction VIMPs \n")}
      
      temp <- vimp(allobj,allobj$xvar.names,importance = importance,joint = F, seed = -1, block = block)
      if (allobj$family == "regr+") {
        allobj$importance <- rowSums(pmax(get.mv.vimp(temp, standardize=TRUE),0))
      } else {
        allobj$importance <- temp$importance
      }
      rm(temp)
      gc()
      
    itrcts <- as.data.frame(do.call(rbind,lapply(1:length(inter),function(i){
      
      if (verbose == TRUE) {cat("\r","Calculating",i,"/",n.vimp, "Permutation-Based Interaction Importance (PBII) indices")}
      
      as.data.frame(hvp(joint = inter[[i]], allobj, importance = importance, block = block
      ))
    })))

      
      if (MDIIorder){
        itrcts$MDII <- as.vector(sub.itr.obj$pmdvp[it.sel,1])
        rm(sub.itr.obj)
        gc() 
        }
    
    return(list(terms = terms,
                interaction = itrcts[order(itrcts$PBII, decreasing = TRUE),],
                rf.obj = allobj, 
                pmd = pmd,
                var.sl = var.sl,
                var.sl.list = var.sl.list,
                itrwt = itrwt,
                MDII = sort(pmdvimp.long),
                err.init = err.init,
                err.sub = err.sub,
                err.all = na.omit(allobj$err.rate) ))
    
    
    
    }}
  } else {
  #  sub.x.obj <- list()
  #  sub.x.obj$var.sl <- sub.x.obj$var.sl.list <- sub.x.obj$itrwt <- NULL
    var.sl <- 
    var.sl.list <- 
    itrwt <- 
    err.init <- 
    err.sub <- NULL
    if (verbose == TRUE) {cat("Fitting RF \n")}
    if (forest.tune) {
    turn.o <- tune(formula, data, mtryStart = (ncol(data) - 1))
    allobj <- rfsrc(formula, data = data, ntime = ntime,
                    importance = importance, seed = -1,
                    mtry = (ncol(data)-1), nodesize = turn.o$rf$nodesize)
    rm(turn.o)
    } else {
      allobj <- rfsrc(formula, data = data, ntime = ntime,
                      importance = importance, seed = -1) 
    }
    gc()
  
  if (!var.sel.only) {
  cmbn <- lapply(2:terms,function(i){t(combn(1:length(allobj$xvar.names),i))})
  inter <- unlist(lapply(1:length(cmbn), function(i){
    lapply(1:nrow(cmbn[[i]]),function(j){ matrix(cmbn[[i]][j,],1,length(cmbn[[i]][j,])) })
  }), recursive = F )
    
    if (MDIIorder){
      n.vimp <- nPBii(inter)
      sub.itr.obj <- pmdhvp(allobj, inter, verbose = verbose)
      pmd <- sub.itr.obj$pmd
      if (length(inter) > n.vimp) {
        if (verbose == TRUE) {cat("\n Selected", n.vimp, "from", length(inter), "interactions based on Minimal Depth Interaction Importance (MDII)\n")}
        it.sel <- which(rank(sub.itr.obj$pmdvp) <= n.vimp)
        inter <- inter[it.sel]
      } else {it.sel <- 1:length(inter)}
      pmdvimp.long <- sub.itr.obj$pmdvp
    }else{
      pmd <- pmdvimp.long <- NULL
      n.vimp <- length(inter)
    }
    
    
  
  if (MDIIonly == TRUE) {
    return(list(terms = terms,
                rf.obj = allobj, pmd = pmd,
                var.sl = var.sl,
                var.sl.list = var.sl.list,
                itrwt = itrwt,
                MDII = pmdvimp.long,
                err.init = err.init,
                err.sub = err.sub,
                err.all = na.omit(allobj$err.rate)
                ))
    
  } else{

 # if (verbose == TRUE) {cat("\n Calculating", n.vimp, "High Order Interaction VIMPs \n")}
    
    temp <- vimp(allobj,allobj$xvar.names,importance = importance,joint = F, seed = -1, block = block)
    if (allobj$family == "regr+") {
      allobj$importance <- rowSums(pmax(get.mv.vimp(temp, standardize=TRUE),0))
    } else {
      allobj$importance <- temp$importance
    }
    rm(temp)
    gc()
    
  itrcts <- as.data.frame(do.call(rbind,lapply(1:length(inter),function(i){
    
    if (verbose == TRUE) {cat("\r","Calculating",i,"/",n.vimp, "Permutation-Based Interaction Importance (PBII) indices")}
    as.data.frame(hvp(joint = inter[[i]], allobj, importance = importance, block = block
    ))
  })))

  
  if (MDIIorder){
    itrcts$MDII <- as.vector(sub.itr.obj$pmdvp[it.sel,1])
    rm(sub.itr.obj)
    gc() 
    }
  


  return(list(terms = terms,
       interaction = itrcts[order(itrcts$PBII, decreasing = TRUE),],
       rf.obj = allobj, pmd = pmd,
       var.sl = var.sl,
       var.sl.list = var.sl.list,
       itrwt = itrwt,
       MDII = sort(pmdvimp.long),
       err.init = err.init,
       err.sub = err.sub,
       err.all = na.omit(allobj$err.rate) ))
  }
  }}
}


Pvalue <- function(o, n.null = 100, importance = "anti", block = 1, verbose = TRUE){
  terms <- o$terms
  allobj <- o$rf.obj
  interactions <- o$interaction
  interactions$PBII <- abs(o$interaction$PBII)
  PBII <- unlist(abs(o$interaction$PBII)) #= o$interaction$nmHIvimp
  rm(o)
  gc()
  
  if (is.null(allobj$importance)){
    temp <- vimp(allobj,allobj$xvar.names,importance = importance,joint = F, seed = -1, block = block)
    if (allobj$family == "regr+") {
      allobj$importance <- rowSums(pmax(get.mv.vimp(temp, standardize=TRUE),0))
    } else {
      allobj$importance <- temp$importance
    }
    rm(temp)
    gc()
  }
  
  cmbn <- lapply(2:terms,function(i){t(combn(1:length(allobj$xvar.names),i))})
  inter <- unlist(lapply(1:length(cmbn), function(i){
    lapply(1:nrow(cmbn[[i]]),function(j){ matrix(cmbn[[i]][j,],1,length(cmbn[[i]][j,])) })
  }), recursive = F )
  
  set.seed(-1)
  inter <- inter[sample(1:length(inter),min(length(inter),n.null))]
  
  null.distr <- as.data.frame(do.call(rbind,lapply(1:length(inter),function(i){
    
    if (verbose == TRUE) {cat("\r","Calculating ",round(100*i/n.null), "% of the NULL distribution", sep = "")}
    as.data.frame(hvp(joint = inter[[i]],allobj, importance = importance, block = block
    ))
  })))
  
  interactions$Pvalue <-  unlist(lapply(1:length(PBII), function(i){
    length(which(null.distr$PBII>=PBII[i]))/n.null
  }))

  list(interaction = interactions, NULLdistr = null.distr$PBII)
}
