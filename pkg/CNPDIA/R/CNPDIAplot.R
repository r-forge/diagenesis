# ==============================================================================
# ==============================================================================

# ------------------------------------------------------------------------------
# Images of output
# ------------------------------------------------------------------------------


image2DDIA <- function (x, select = NULL, which = select, subset = NULL, ...) {

  if (!"deSolve"%in% class(x)) stop ("cannot make an image from output that is not dynamic")
  ldots <- list(...)
  nmdots <- names(ldots)

# subsetting
      if (!missing(subset)) {
        e <- substitute(subset)
        r <- eval(e, as.data.frame(x), parent.frame())
        if (is.numeric(r)) {
            isub <- r
        }
        else {
            if (!is.logical(r)) 
                stop("'subset' must evaluate to logical or be a vector with integers")
            isub <- r & !is.na(r)
        }
    }
    else isub <- TRUE

    att    <- attributes(x)
    nspec  <- att$nspec
    dimens <- att$dimens
    proddim <- prod(dimens)
    if ((ncol(x) - nspec * proddim) < 1) 
        stop("ncol of 'x' should be > 'nspec' * dimens if x is a vector")

  if ("MPBDIAdyn" %in% class(x))
      Dec <- rbind(MPBDIAsvar(), MPBDIA1D())
  else
      Dec <- rbind(CNPDIAsvar(), CNPDIA1D())

# variables to plot
  if ("grid" %in% nmdots){
    grid <- ldots$grid
    if (is.list(grid))
      grid <- grid$x.mid
    ldots$grid <- NULL
  } else if (att$model == "CNPDIA_model_2"){
     grid <- c(-mean(x[, "Hwater"]), cumsum(att$dx))
   } else { 
    grid <- c(0, cumsum(att$dx))
  }
  Which <- which
  if (is.null(Which) | missing(which))
    Which <- 1:nspec

  if (is.numeric(Which))
    Which <- Dec$names[Which]

  nv <- length(Which)

# the ellipsis
  
    if (!any(match(nmdots, c("mfrow", "mfcol"), nomatch = 0))) {
        nc <- min(ceiling(sqrt(nv)), 3)
        nr <- min(ceiling(nv/nc), 3)
        mfrow <- c(nr, nc)
    }
    else if ("mfcol" %in% nmdots) 
        mfrow <- rev(ldots$mfcol)
    else mfrow <- ldots$mfrow
    if (!is.null(mfrow)) 
        mf <- par(mfrow = mfrow)
    
    ask <- ldots$ask
    ldots$ask <- NULL
    if (is.null(ask)) 
        ask <- prod(par("mfrow")) < nv && dev.interactive()
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    
    repdots <- function (dots, n) 
       if (is.function(dots)) dots else rep(dots, length.out = n)
    expanddots <- function (dots, default, n) {
      dots <- if (is.null(dots)) default else dots
      rep(dots, length.out = n)
    }
    expanddotslist <- function (dots, n) {
     if (is.null(dots)) 
        return(dots)
     dd <- if (!is.list(dots)) 
        list(dots)
     else dots
     rep(dd, length.out = n)
    }
    extractdots<- function (dots, index) {
     ret <- lapply(dots, "[", index)
     ret <- lapply(ret, unlist)
     return(ret)
    }

   if ("legend" %in% nmdots){
    if ("colkey" %in% nmdots) stop ("'legend' and 'colkey' cannot both be specified")
       if (ldots$legend == FALSE) ldots$colkey <- FALSE
    ldots$legend <- NULL
   } 
    if ("colkey" %in% nmdots)
     colkey <- ldots$colkey
    else
     colkey <- NULL
    ldots$colkey <- NULL

        
    Dotmain      <- lapply(ldots, repdots, nv)
    Dotmain$main <- expanddots(ldots$main, Which, nv)
    Dotmain$xlab <- expanddots(ldots$xlab, "times", nv)
    Dotmain$ylab <- expanddots(ldots$ylab, "cm", nv)
    if ("clab" %in% nmdots)
      Clab <- rep(ldots$clab, length.out = nv)
    else{
      Un <- Dec[match(Which, Dec$names),"units"]
      Clab <- .CNPDIA$labels[match(Un, .CNPDIA$labels$Units),"Labels"]
    }    
    Dotmain$clab <- Clab

    xxlim <- expanddotslist(ldots$xlim, nv)
    yylim <- expanddotslist(ldots$ylim, nv)
    zzlim <- expanddotslist(ldots$zlim, nv)
    times <- x[isub, 1]

    for (ip in 1:nv) {
        z <- subset(x, which = Which[ip])[isub,]

        dotmain <- extractdots(Dotmain, ip)
        if (!is.null(xxlim)) 
            dotmain$xlim <- xxlim[[ip]]
        if (!is.null(yylim)) 
            dotmain$ylim <- yylim[[ip]]
        if (!is.null(zzlim)) 
            dotmain$zlim <- zzlim[[ip]]
        else dotmain$zlim <- range(z, na.rm = TRUE)
        List <- alist(z = z, x = times, y = grid, colkey = colkey)

        do.call("image2D", c(List, dotmain))
    }
}



image2D.CNPDIAdyn <- function (z, which, ylim = c(20, 0), 
    colkey = list(cex.clab = 0.8, line.clab = 0.5, cex.axis = 0.8), ...) {
   if (missing(which)) which <- NULL
   image2DDIA(x = z, which = which, ylim = ylim, colkey = colkey, ...)
  }


image2D.MPBDIAdyn <- function (z, which, ylim = c(20, 0), 
    colkey = list(cex.clab = 0.8, line.clab = 0.5, cex.axis = 0.8), ...) {
   if (missing(which)) which <- NULL
   image2DDIA(x = z, which = which, ylim = ylim, colkey = colkey, ...)
  }
  


# ------------------------------------------------------------------------------
# S4 methods
# ------------------------------------------------------------------------------

setGeneric("image2D", function(z, ...) plot3D::image2D(z, ...))
setOldClass("CNPDIAdyn")
setMethod("image2D", signature("CNPDIAdyn"), image2D.CNPDIAdyn)
setOldClass("MPBDIAdyn")
setMethod("image2D", signature("MPBDIAdyn"), image2D.MPBDIAdyn)


matplot.1D.dia <- function (z, which, ylim = c(20, 0), ...) {

  if (!"deSolve"%in% class(z)) stop ("cannot make an image from output that is not dynamic")
  ldots <- list(...)
  nmdots <- names(ldots)
  
  if ("MPBDIAdyn" %in% class(z))
      Dec <- rbind(MPBDIAsvar(), MPBDIA1D())
  else
      Dec <- rbind(CNPDIAsvar(), CNPDIA1D())

  if ("grid" %in% nmdots){
    grid <- ldots$grid
    if (is.list(grid))
      grid <- grid$x.mid
    ldots$grid <- NULL
  } else 
    grid <- attributes(z)$Depth
  
  
  Which <- which
  if (is.null(Which) | missing(Which))
    Which <- attributes(z)$ynames

  if (is.numeric(Which))
    Which <- Dec$names[Which]
  
  nv <- length(Which)
  if ("main" %in% nmdots)
    Main <- rep(ldots$main, length.out = nv)
  else
    Main <- NULL
  ldots$main <- NULL

  if ("xlab" %in% nmdots)
    Xlab <- rep(ldots$xlab, length.out = nv)
  else{
    Un <- Dec[match(Which, Dec$names),"units"]
    Xlab <- .CNPDIA$labels[match(Un, .CNPDIA$labels$Units),"Labels"]
  }    
  ldots$xlab <- NULL

  if ("ylab" %in% nmdots)
    Ylab <- rep(ldots$ylab, length.out = nv)
  else{
    Ylab <- "cm"
  }    
  ldots$ylab <- NULL  
  
  if (!any(match(nmdots, c("mfrow", "mfcol"), nomatch = 0))) {
     nc <- min(ceiling(sqrt(nv)), 3)
     nr <- min(ceiling(nv/nc), 3)
     mfrow <- c(nr, nc)
  } else if ("mfcol" %in% nmdots)
     mfrow <- rev(ldots$mfcol)
  else mfrow <- ldots$mfrow

  if (! is.null(mfrow))  mf <- par(mfrow = mfrow)
  ## interactively wait if there are remaining figures

  if (is.null(ldots$ask))
    ldots$ask <- prod(par("mfrow")) < nv && dev.interactive()
  par(ask = ldots$ask)
  ldots$ask <- NULL

  if ("xyswap" %in% nmdots) {
    xyswap <- ldots$xyswap
    ldots$xyswap <- NULL
  } else xyswap <- TRUE
    
  class(z) <- class(z)[-1]  
  do.call("matplot.1D",
            c(alist(x = z, which = Which, grid = grid, ylim = ylim, 
                    xlab = Xlab, ylab = Ylab, xyswap = xyswap), ldots))
}

matplot1D.CNPDIAdyn <- function (z, which, ylim = c(20, 0), 
   type = "l", col = "grey", lty = 1, ...) {
   if (missing(which)) which <- NULL
   matplot.1D.dia(z = z, which = which, ylim = ylim, 
     type = type, col = col, lty = lty, ...)
  }


matplot1D.MPBDIAdyn <- function (z, which, ylim = c(20,0), 
   type = "l", col = "grey",  lty = 1, ...) {
  if (missing(which)) which <- NULL
   matplot.1D.dia(z = z, which = which, ylim = ylim, 
     type = type, col = col,  lty = lty, ...)
  }
  

# ------------------------------------------------------------------------------
# S4 methods
# ------------------------------------------------------------------------------
matplot1D <- function (z, ...) UseMethod("matplot1D")

matplot1D.default <- function (z, ...) {
if (inherits (z, "CNPDIAdyn"))
  matplot1D.CNPDIAdyn(z,...)
else
  deSolve::matplot.1D(x = z,...)
#  graphics::matplot(x,...)
  NextMethod()
}

setGeneric("matplot1D", function(z, ...) matplot1D(z, ...))
setOldClass("CNPDIAdyn")
setMethod("matplot1D", signature("CNPDIAdyn"), matplot1D.CNPDIAdyn)
setOldClass("MPBDIAdyn")
setMethod("matplot1D", signature("MPBDIAdyn"), matplot1D.MPBDIAdyn)









plot.dia <- function (x, ..., which, ylim = c(20, 0)) {

  if (!"steady1D"%in% class(x)) stop ("cannot make a 1D plot from output that is not steady-state")
  ldots <- list(...)
  nmdots <- names(ldots)
  
  if ("MPBDIAstd" %in% class(x))
      Dec <- rbind(MPBDIAsvar(), MPBDIA1D())
  else
      Dec <- rbind(CNPDIAsvar(), CNPDIA1D())

  if ("grid" %in% nmdots){
    grid <- ldots$grid
    if (is.list(grid))
      grid <- grid$x.mid
    ldots$grid <- NULL
  } else 
    grid <- x$Depth
  
  
  Which <- which
  if (is.null(Which) | missing(Which))
    Which <- attributes(x)$ynames

  if (is.numeric(Which))
    Which <- Dec$names[Which]
  
  nv <- length(Which)
  if ("main" %in% nmdots)
    Main <- rep(ldots$main, length.out = nv)
  else
    Main <- NULL
  ldots$main <- NULL

  if ("xlab" %in% nmdots)
    Xlab <- rep(ldots$xlab, length.out = nv)
  else{
    Un <- Dec[match(Which, Dec$names),"units"]
    Xlab <- .CNPDIA$labels[match(Un, .CNPDIA$labels$Units),"Labels"]
  }    
  ldots$xlab <- NULL

  if ("ylab" %in% nmdots)
    Ylab <- rep(ldots$ylab, length.out = nv)
  else{
    Ylab <- "cm"
  }    
  ldots$ylab <- NULL  
  
  if (!any(match(nmdots, c("mfrow", "mfcol"), nomatch = 0))) {
     nc <- min(ceiling(sqrt(nv)), 3)
     nr <- min(ceiling(nv/nc), 3)
     mfrow <- c(nr, nc)
  } else if ("mfcol" %in% nmdots)
     mfrow <- rev(ldots$mfcol)
  else mfrow <- ldots$mfrow

  if (! is.null(mfrow))  mf <- par(mfrow = mfrow)
  ## interactively wait if there are remaining figures

  if (is.null(ldots$ask))
    ldots$ask <- prod(par("mfrow")) < nv && dev.interactive()
  par(ask = ldots$ask)
  ldots$ask <- NULL

  if ("xyswap" %in% nmdots) {
    xyswap <- ldots$xyswap
    ldots$xyswap <- NULL
  } else xyswap <- TRUE
    
  
  X <- x
  class(X) <- class(X)[-1]
  do.call("plot",
            c(alist(x = X, xyswap = xyswap,
                    which = Which, grid = grid, ylim = ylim, 
                    xlab = Xlab, ylab = Ylab), ldots))
}

plotstd.generic <- function(x, ...)
  plot(x, ...)

plot.CNPDIAstd <- function (x, ..., which, ylim = c(20, 0)) {
   if (missing(which)) which <- NULL
   plot.dia(x = x, ..., which = which, ylim = ylim)
  }


plot.MPBDIAstd <- function (x, ..., which, ylim = c(20,0)) {
  if (missing(which)) which <- NULL
   plot.dia(x = x, ..., which = which, ylim = ylim)
  }
  

# ------------------------------------------------------------------------------
# S4 methods
# ------------------------------------------------------------------------------

#setGeneric("plot", function(x, ...) plot(x, ...))
setOldClass("CNPDIAstd")
setMethod("plot", signature("CNPDIAstd"), plot.CNPDIAstd)
setOldClass("MPBDIAstd")
setMethod("plot", signature("MPBDIAstd"), plot.MPBDIAstd)



