
#### From http://ase.tufts.edu/gsc/gradresources/guidetomixedmodelsinr/mixed%20model%20guide.html

overdisp_fun <- function(model) {
  ## number of variance parameters in an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m) * (nrow(m) + 1)/2
  }
  # The next two lines calculate the residual degrees of freedom
  model.df <- sum(sapply(VarCorr(model), vpars)) + length(fixef(model))
  rdf <- nrow(model.frame(model)) - model.df
  # extracts the Pearson residuals
  rp <- residuals(model, type = "pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  # Generates a p-value. If less than 0.05, the data are overdispersed.
  pval <- pchisq(Pearson.chisq, df = rdf, lower.tail = FALSE)
  c(chisq = Pearson.chisq, ratio = prat, rdf = rdf, p = pval)
}






CSQPlot<-function(vars,label="Chi-Square Quantile Plot"){
  #usually, vars is xxx$residuals or data from one group and label is for plot
  x<-cov(scale(vars),use="pairwise.complete.obs")
  squares<-sort(diag(as.matrix(scale(vars))%*%solve(x)%*%as.matrix(t(scale(vars)))))
  quantiles<-quantile(squares)
  hspr<-quantiles[4]-quantiles[2]
  cumprob<-c(1:length(vars[,1]))/length(vars[,1])-1/(2*length(vars[,1]))
  degf<-dim(x)[1]
  quants<-qchisq(cumprob,df=degf)
  gval<-(quants**(-1+degf/2))/(exp(quants/2)*gamma(degf/2)*(sqrt(2)**degf))
  scale<-hspr / (qchisq(.75,degf)-qchisq(.25,degf))
  se<-(scale/gval)*sqrt(cumprob*(1-cumprob)/length(squares))
  lower<-quants-2*se
  upper<-quants+2*se
  
  plot(quants,squares,col='red',pch=19,cex=1.2,xlab="Chi-Square Quantiles",
       ylab=label,main=paste("Chi-Square Quantiles for",label),ylim=range(upper,lower, squares) , xlim=range(c(0,quants)))
  lines(c(0,100),c(0,100),col=1)
  lines(quants,upper,col="blue",lty=2,lwd=2)
  lines(quants,lower,col="blue",lty=2,lwd=2)
  legend(0,range(upper,lower)[2]*.9,c("Data","95% Conf Limits"),lty=c(0,2),col=c("red","blue"),lwd=c(2,2),
         pch=c(19,NA))
}


modelassump <- function(model){
  d <- par(mfrow=c(2,2), mar=c(4,4,1,1))
  qqPlot(residuals(model))
  hist(residuals(model, type = "response"), breaks = 100, main = "")
  plot(residuals(model) ~ fitted(model))
  abline(h=0, col="red")
  #rootogram(model, style = "hanging")
  #plot(0,type='n',axes=FALSE,ann=FALSE)
  acf(residuals(model))
  par(d)
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


#Zuur et al. 2009
##################################################################
#Below are some functions that we took from the pairs help file and
#modified, or wrote ourselves. You need to copy and paste all these commands
#into R. It is perhaps better not to try and understand what it all does.



## put correlations on the panels,
## with size proportional to the correlations.
panel.cor <- function(x, y, digits=1, prefix="", cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r1=cor(x,y,use="pairwise.complete.obs")
  r <- abs(cor(x, y,use="pairwise.complete.obs"))
  
  txt <- format(c(r1, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex <- 0.9/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex * r)
}

panel.smooth2=function (x, y, col = par("col"), bg = NA, pch = par("pch"),
                        cex = 1, col.smooth = "red", span = 2/3, iter = 3, ...)
{
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok))
    lines(stats::lowess(x[ok], y[ok], f = span, iter = iter),
          col = 1, ...)
}


panel.lines2=function (x, y, col = par("col"), bg = NA, pch = par("pch"),
                       cex = 1, ...)
{
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)){
    tmp=lm(y[ok]~x[ok])
    abline(tmp)}
  
}




panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="white", ...)
}



#VIF
myvif <- function(mod) {
  v <- vcov(mod)
  assign <- attributes(model.matrix(mod))$assign
  if (names(coefficients(mod)[1]) == "(Intercept)") {
    v <- v[-1, -1]
    assign <- assign[-1]
  } else warning("No intercept: vifs may not be sensible.")
  terms <- labels(terms(mod))
  n.terms <- length(terms)
  if (n.terms < 2) stop("The model contains fewer than 2 terms")
  if (length(assign) > dim(v)[1] ) {
    diag(tmp_cor)<-0
    if (any(tmp_cor==1.0)){
      return("Sample size is too small, 100% collinearity is present")
    } else {
      return("Sample size is too small")
    }
  }
  R <- cov2cor(v)
  detR <- det(R)
  result <- matrix(0, n.terms, 3)
  rownames(result) <- terms
  colnames(result) <- c("GVIF", "Df", "GVIF^(1/2Df)")
  for (term in 1:n.terms) {
    subs <- which(assign == term)
    result[term, 1] <- det(as.matrix(R[subs, subs])) * det(as.matrix(R[-subs, -subs])) / detR
    result[term, 2] <- length(subs)
  }
  if (all(result[, 2] == 1)) {
    result <- data.frame(GVIF=result[, 1])
  } else {
    result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
  }
  invisible(result)
}

corvif <- function(dataz) {
  dataz <- as.data.frame(dataz)
  #correlation part
  cat("Correlations of the variables\n\n")
  tmp_cor <- cor(dataz,use="complete.obs")
  print(tmp_cor)
  
  #vif part
  form    <- formula(paste("fooy ~ ",paste(strsplit(names(dataz)," "),collapse=" + ")))
  dataz   <- data.frame(fooy=1,dataz)
  lm_mod  <- lm(form,dataz)
  
  cat("\n\nVariance inflation factors\n\n")
  print(myvif(lm_mod))
}

myvif <- function(mod) {
  v <- vcov(mod)
  assign <- attributes(model.matrix(mod))$assign
  if (names(coefficients(mod)[1]) == "(Intercept)") {
    v <- v[-1, -1]
    assign <- assign[-1]
  } else warning("No intercept: vifs may not be sensible.")
  terms <- labels(terms(mod))
  n.terms <- length(terms)
  if (n.terms < 2) stop("The model contains fewer than 2 terms")
  if (length(assign) > dim(v)[1] ) {
    diag(tmp_cor)<-0
    if (any(tmp_cor==1.0)){
      return("Sample size is too small, 100% collinearity is present")
    } else {
      return("Sample size is too small")
    }
  }
  R <- cov2cor(v)
  detR <- det(R)
  result <- matrix(0, n.terms, 3)
  rownames(result) <- terms
  colnames(result) <- c("GVIF", "Df", "GVIF^(1/2Df)")
  for (term in 1:n.terms) {
    subs <- which(assign == term)
    result[term, 1] <- det(as.matrix(R[subs, subs])) * det(as.matrix(R[-subs, -subs])) / detR
    result[term, 2] <- length(subs)
  }
  if (all(result[, 2] == 1)) {
    result <- data.frame(GVIF=result[, 1])
  } else {
    result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
  }
  invisible(result)
}



corvif <- function(dataz) {
  dataz <- as.data.frame(dataz)
  #correlation part
  cat("Correlations of the variables\n\n")
  tmp_cor <- cor(dataz,use="complete.obs")
  print(tmp_cor)
  
  #vif part
  form    <- formula(paste("fooy ~ ",paste(strsplit(names(dataz)," "),collapse=" + ")))
  dataz   <- data.frame(fooy=1,dataz)
  lm_mod  <- lm(form,dataz)
  
  cat("\n\nVariance inflation factors\n\n")
  print(myvif(lm_mod))
}



scatterhist = function(x, y, xlab="", ylab=""){
  zones=matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
  layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
  xhist = hist(x, plot=FALSE)
  yhist = hist(y, plot=FALSE)
  top = max(c(xhist$counts, yhist$counts))
  par(mar=c(3,3,1,1))
  plot(jitter(x),y)
  par(mar=c(0,3,1,1))
  barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0)
  par(mar=c(3,0,1,1))
  barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE)
  par(oma=c(3,3,0,0))
  mtext(xlab, side=1, line=1, outer=TRUE, adj=0, 
        at=.8 * (mean(x) - min(x))/(max(x)-min(x)))
  mtext(ylab, side=2, line=1, outer=TRUE, adj=0, 
        at=(.8 * (mean(y) - min(y))/(max(y) - min(y))))
}

lm_eqn = function(m) {
  
  l <- list(a = format(coef(m)[1], digits = 2),
            b = format(abs(coef(m)[2]), digits = 2),
            r2 = format(summary(m)$r.squared, digits = 3));
  
  if (coef(m)[2] >= 0)  {
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
  } else {
    eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(r)^2~"="~r2,l)    
  }
  
  as.character(as.expression(eq));                 
}


####################################################
## LEGEND
####################################################

LEGEND <- function (x, y = NULL, legend, fill = NULL, 
                    col = par("col"), pt.col=col, line.col=col,
                    border = "black", lty, lwd, pch, angle = 45, density = NULL,
                    bty = "o", bg = par("bg"), box.lwd = par("lwd"), box.lty = par("lty"),
                    box.col = par("fg"), pt.bg = NA, cex = 1, pt.cex = cex, pt.lwd = lwd,
                    xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1, adj = c(0,
                                                                                0.5), text.width = NULL, text.col = par("col"), text.font = NULL,
                    merge = do.lines && has.pch, trace = FALSE, plot = TRUE,
                    ncol = 1, horiz = FALSE, title = NULL, inset = 0, xpd, title.col = text.col,
                    title.adj = 0.5, seg.len = 2)
{
  if (missing(legend) && !missing(y) && (is.character(y) ||
                                         is.expression(y))) {
    legend <- y
    y <- NULL
  }
  mfill <- !missing(fill) || !missing(density)
  if (!missing(xpd)) {
    op <- par("xpd")
    on.exit(par(xpd = op))
    par(xpd = xpd)
  }
  title <- as.graphicsAnnot(title)
  if (length(title) > 1)
    stop("invalid 'title'")
  legend <- as.graphicsAnnot(legend)
  n.leg <- if (is.call(legend))
    1
  else length(legend)
  if (n.leg == 0)
    stop("'legend' is of length 0")
  auto <- if (is.character(x))
    match.arg(x, c("bottomright", "bottom", "bottomleft",
                   "left", "topleft", "top", "topright", "right", "center"))
  else NA
  if (is.na(auto)) {
    xy <- xy.coords(x, y)
    x <- xy$x
    y <- xy$y
    nx <- length(x)
    if (nx < 1 || nx > 2)
      stop("invalid coordinate lengths")
  }
  else nx <- 0
  xlog <- par("xlog")
  ylog <- par("ylog")
  rect2 <- function(left, top, dx, dy, density = NULL, angle,
                    ...) {
    r <- left + dx
    if (xlog) {
      left <- 10^left
      r <- 10^r
    }
    b <- top - dy
    if (ylog) {
      top <- 10^top
      b <- 10^b
    }
    rect(left, top, r, b, angle = angle, density = density,
         ...)
  }
  segments2 <- function(x1, y1, dx, dy, ...) {
    x2 <- x1 + dx
    if (xlog) {
      x1 <- 10^x1
      x2 <- 10^x2
    }
    y2 <- y1 + dy
    if (ylog) {
      y1 <- 10^y1
      y2 <- 10^y2
    }
    segments(x1, y1, x2, y2, ...)
  }
  points2 <- function(x, y, ...) {
    if (xlog)
      x <- 10^x
    if (ylog)
      y <- 10^y
    points(x, y, ...)
  }
  text2 <- function(x, y, ...) {
    if (xlog)
      x <- 10^x
    if (ylog)
      y <- 10^y
    text(x, y, ...)
  }
  if (trace)
    catn <- function(...) do.call("cat", c(lapply(list(...),
                                                  formatC), list("\n")))
  cin <- par("cin")
  Cex <- cex * par("cex")
  if (is.null(text.width))
    text.width <- max(abs(strwidth(legend, units = "user",
                                   cex = cex, font = text.font)))
  else if (!is.numeric(text.width) || text.width < 0)
    stop("'text.width' must be numeric, >= 0")
  xc <- Cex * xinch(cin[1L], warn.log = FALSE)
  yc <- Cex * yinch(cin[2L], warn.log = FALSE)
  if (xc < 0)
    text.width <- -text.width
  xchar <- xc
  xextra <- 0
  yextra <- yc * (y.intersp - 1)
  ymax <- yc * max(1, strheight(legend, units = "user", cex = cex)/yc)
  ychar <- yextra + ymax
  if (trace)
    catn("  xchar=", xchar, "; (yextra,ychar)=", c(yextra,
                                                   ychar))
  if (mfill) {
    xbox <- xc * 0.8
    ybox <- yc * 0.5
    dx.fill <- xbox
  }
  do.lines <- (!missing(lty) && (is.character(lty) || any(lty >
                                                            0))) || !missing(lwd)
  n.legpercol <- if (horiz) {
    if (ncol != 1)
      warning(gettextf("horizontal specification overrides: Number of columns := %d",
                       n.leg), domain = NA)
    ncol <- n.leg
    1
  }
  else ceiling(n.leg/ncol)
  has.pch <- !missing(pch) && length(pch) > 0
  if (do.lines) {
    x.off <- if (merge)
      -0.7
    else 0
  }
  else if (merge)
    warning("'merge = TRUE' has no effect when no line segments are drawn")
  if (has.pch) {
    if (is.character(pch) && !is.na(pch[1L]) && nchar(pch[1L],
                                                      type = "c") > 1) {
      if (length(pch) > 1)
        warning("not using pch[2..] since pch[1L] has multiple chars")
      np <- nchar(pch[1L], type = "c")
      pch <- substr(rep.int(pch[1L], np), 1L:np, 1L:np)
    }
    if (!is.character(pch))
      pch <- as.integer(pch)
  }
  if (is.na(auto)) {
    if (xlog)
      x <- log10(x)
    if (ylog)
      y <- log10(y)
  }
  if (nx == 2) {
    x <- sort(x)
    y <- sort(y)
    left <- x[1L]
    top <- y[2L]
    w <- diff(x)
    h <- diff(y)
    w0 <- w/ncol
    x <- mean(x)
    y <- mean(y)
    if (missing(xjust))
      xjust <- 0.5
    if (missing(yjust))
      yjust <- 0.5
  }
  else {
    h <- (n.legpercol + (!is.null(title))) * ychar + yc
    w0 <- text.width + (x.intersp + 1) * xchar
    if (mfill)
      w0 <- w0 + dx.fill
    if (do.lines)
      w0 <- w0 + (seg.len + x.off) * xchar
    w <- ncol * w0 + 0.5 * xchar
    if (!is.null(title) && (abs(tw <- strwidth(title, units = "user",
                                               cex = cex) + 0.5 * xchar)) > abs(w)) {
      xextra <- (tw - w)/2
      w <- tw
    }
    if (is.na(auto)) {
      left <- x - xjust * w
      top <- y + (1 - yjust) * h
    }
    else {
      usr <- par("usr")
      inset <- rep_len(inset, 2)
      insetx <- inset[1L] * (usr[2L] - usr[1L])
      left <- switch(auto, bottomright = , topright = ,
                     right = usr[2L] - w - insetx, bottomleft = ,
                     left = , topleft = usr[1L] + insetx, bottom = ,
                     top = , center = (usr[1L] + usr[2L] - w)/2)
      insety <- inset[2L] * (usr[4L] - usr[3L])
      top <- switch(auto, bottomright = , bottom = , bottomleft = usr[3L] +
                      h + insety, topleft = , top = , topright = usr[4L] -
                      insety, left = , right = , center = (usr[3L] +
                                                             usr[4L] + h)/2)
    }
  }
  if (plot && bty != "n") {
    if (trace)
      catn("  rect2(", left, ",", top, ", w=", w, ", h=",
           h, ", ...)", sep = "")
    rect2(left, top, dx = w, dy = h, col = bg, density = NULL,
          lwd = box.lwd, lty = box.lty, border = box.col)
  }
  xt <- left + xchar + xextra + (w0 * rep.int(0:(ncol - 1),
                                              rep.int(n.legpercol, ncol)))[1L:n.leg]
  yt <- top - 0.5 * yextra - ymax - (rep.int(1L:n.legpercol,
                                             ncol)[1L:n.leg] - 1 + (!is.null(title))) * ychar
  if (mfill) {
    if (plot) {
      if (!is.null(fill))
        fill <- rep_len(fill, n.leg)
      rect2(left = xt, top = yt + ybox/2, dx = xbox, dy = ybox,
            col = fill, density = density, angle = angle,
            border = border)
    }
    xt <- xt + dx.fill
  }
  if (plot && (has.pch || do.lines)) {
    pt.COL <- rep_len(pt.col, n.leg)
    line.COL <- rep_len(line.col, n.leg)
  }
  if (missing(lwd) || is.null(lwd))
    lwd <- par("lwd")
  if (do.lines) {
    if (missing(lty) || is.null(lty))
      lty <- 1
    lty <- rep_len(lty, n.leg)
    lwd <- rep_len(lwd, n.leg)
    ok.l <- !is.na(lty) & (is.character(lty) | lty > 0) &
      !is.na(lwd)
    if (trace)
      catn("  segments2(", xt[ok.l] + x.off * xchar, ",",
           yt[ok.l], ", dx=", seg.len * xchar, ", dy=0, ...)")
    if (plot)
      segments2(xt[ok.l] + x.off * xchar, yt[ok.l], dx = seg.len *
                  xchar, dy = 0, lty = lty[ok.l], lwd = lwd[ok.l],
                col = line.COL[ok.l])
    xt <- xt + (seg.len + x.off) * xchar
  }
  if (has.pch) {
    pch <- rep_len(pch, n.leg)
    pt.bg <- rep_len(pt.bg, n.leg)
    pt.cex <- rep_len(pt.cex, n.leg)
    pt.lwd <- rep_len(pt.lwd, n.leg)
    ok <- !is.na(pch)
    if (!is.character(pch)) {
      ok <- ok & (pch >= 0 | pch <= -32)
    }
    else {
      ok <- ok & nzchar(pch)
    }
    x1 <- (if (merge && do.lines)
      xt - (seg.len/2) * xchar
      else xt)[ok]
    y1 <- yt[ok]
    if (trace)
      catn("  points2(", x1, ",", y1, ", pch=", pch[ok],
           ", ...)")
    if (plot)
      points2(x1, y1, pch = pch[ok], col = pt.COL[ok], cex = pt.cex[ok],
              bg = pt.bg[ok], lwd = pt.lwd[ok])
  }
  xt <- xt + x.intersp * xchar
  if (plot) {
    if (!is.null(title))
      text2(left + w * title.adj, top - ymax, labels = title,
            adj = c(title.adj, 0), cex = cex, col = title.col)
    text2(xt, yt, labels = legend, adj = adj, cex = cex,
          col = text.col, font = text.font)
  }
  invisible(list(rect = list(w = w, h = h, left = left, top = top),
                 text = list(x = xt, y = yt)))
}




###########################################################
## Scatter plot matrix
###########################################################

# Correlation panel
panel.cor <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y), digits=2)
  txt <- paste0("R = ", r)
  cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = ifelse((cex.cor * r) > 0.6 | (cex.cor*r) < -0.6, (cex.cor * r), 0.6), col = ifelse(r > 0.5 | r < -0.5, "red", "black"))
}

