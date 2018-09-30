###=============================================================================
## SIGDA:  Scale Invariant Geometric Data Analysis
#
#  Author:  Max Robinson
#            Institute for Systems Biology
#            max.robinson@systemsbiology.org
#
#  first version:  10 October 2015
#  latest version: 30 September 2018
###=============================================================================
library(Matrix)
library(irlba)

###-----------------------------------------------------------------------------
### Interface functions
###
### projective.decomposition()  -- a matrix scaling function that equalizes RMS
### sigda()                     -- scale invariant geometric data analysis
###
### High-level plotting functions:
###     espalier.plot()         -- plot a scale invariant axis against scale
###     espalier.row.plot()     -- plot row points only on a scale invariant axis against scale
###     espalier.col.plot()     -- plot column points only on a scale invariant axis against scale
###
###     sigda.biplot()          -- plot on two scale invariant axes
###     sigda.row.biplot()      -- plot row points only on two scale invariant axes
###     sigda.column.biplot()   -- plot column points only on two scale invariant axes
###
### Low-level plotting functions:
###     double.biplot()         -- sigda.biplot() for other data
###     single.biplot()         -- sigda.{row,column}.biplot() for other data
###-----------------------------------------------------------------------------

#==============================================================================
# This algorithm is iterative, and produces an approximation to a
# given error tolerance.  The code appears to be numerically stable,
# and using 64-bit floating point representations can achieve a
# coefficient of variation tolerance of 1e-14, sometimes 1e-15,
# but in practice 1e-6 has appeared to be sufficient.
###-------------------------------------------------------------------
# M =    rms(A)  *  Diag{row.factor}    %*%  W  %*%  Diag{col.factor}
# W = (1/rms(A)) *  Diag{1/row.factor}  %*%  A  %*%  Diag{1/col.factor}
# scalar = rms(A) = Frobenius(A) / sqrt(m*n)
###-------------------------------------------------------------------
# Methods: Only the first letter matters, capitalization doesn't
#   Sinkhorn-Knopp  row, column, repeat, full steps
#   Ruiz            row and column, half steps (sqrt)
#   Hybrid          Sinkhorn (row), Sinkhorn (column), Ruiz; repeat
# Sinkhorn-Knopp appears to be the fastest without parallelization
# Ruiz is more parallelizable
# Hybrid _might_ sometimes be faster than Sinkhorn-Knopp, but I doubt it.
###-------------------------------------------------------------------
projective.decomposition <- function(
        A,
        D.start  = NULL,
        method   = 'Sinkhorn',
        tol      = 1e-6,
        max.iter = 1000,
        verbose  = FALSE,
        digits   = 3) {
    m <- dim(A)[1]
    n <- dim(A)[2]
    sqrt.m  <- sqrt(m)
    sqrt.n  <- sqrt(n)
    sqrt.mn <- sqrt.m * sqrt.n
    D <- list('scalar'     = Frobenius(A) / sqrt.mn,
              'row.factor' = vector(mode='numeric', length=m),
              'col.factor' = vector(mode='numeric', length=n),
              'tolerance'  = tol
    )

    # Initial scaling factors
    W2   <- as(A,'dgTMatrix')
    W2@x <- W2@x / D$scalar
    W2@x <- W2@x * W2@x
    a2   <- as.vector(W2@x)
    if (  is.null(D.start)
       || (  (m != length(D.start$row.factor))
          && (n != length(D.start$col.factor)))) {
        D$row.factor <- sqrt(rowSums(W2)/n)
        D$col.factor <- sqrt(colSums(W2)/m)
    } else {
        if (m == length(D.start$row.factor)) {
            if (n == length(D.start$col.factor)) {
                if (verbose) { print("Reusing row and column factors") }
                D$row.factor <- D.start$row.factor * D.start$row.factor
                D$col.factor <- D.start$col.factor * D.start$col.factor
            } else {
                if (verbose) { print("Reusing row factors") }
                D$row.factor <- D.start$row.factor * D.start$row.factor
                W2@x <- a2/D$row.factor[1+W2@i]
                D$col.factor <- colSums(W2)/m
                W2@x <- a2/(D$row.factor[1+W2@i]*D$col.factor[1+W2@j])
                D$row.factor <- D$row.factor*(rowSums(W2)/n)
            }
        } else { # (n == length(D.start$col.factor))
            if (verbose) { print("Reusing column factors") }
            D$col.factor <- D.start$col.factor * D.start$col.factor
            W2@x <- a2/D$col.factor[1+W2@j]
            D$row.factor <- rowSums(W2)/n
            W2@x <- a2/(D$row.factor[1+W2@i]*D$col.factor[1+W2@j])
            D$col.factor <- D$col.factor*(colSums(W2)/m)
        }
    }
    W2@x <- a2/(D$row.factor[1+W2@i]*D$col.factor[1+W2@j])

    # Balancing the RMS via mean L1 and Hadamard-squaring
    style <- toupper(substr(method,1,1))
    for (iter in c(1:max.iter)) {
        rL <- rowSums(W2)/n
        cL <- colSums(W2)/m
        coeff.v <- max(cv(cL[cL > 0]), cv(rL[rL > 0]))
        if (verbose) {
            print(paste(iter,signif(coeff.v,digits),date()))
        }
        if (coeff.v < 2*tol) { break }
        # Stopping condition; 2*tol comes from scaling rms() via the L1 norm

        if ('S' == style) { # Sinkhorn-style algorithm
            D$row.factor <- D$row.factor * rL
            W2@x <- a2/(D$row.factor[1+W2@i]*D$col.factor[1+W2@j])
            D$col.factor <- D$col.factor * cL
        } else if ('R' == style) { # Ruiz-style algorithm
            D$col.factor <- D$col.factor * sqrt(cL)
            D$row.factor <- D$row.factor * sqrt(rL)
        } else { # Hybrid Ruiz/Sinkhorn-style algorithm
            if (1 == (iter %% 2)) {
                D$col.factor <- D$col.factor * sqrt(cL)
                D$row.factor <- D$row.factor * sqrt(rL)
            } else {
                D$row.factor <- D$row.factor * rL
                W2@x <- a2/(D$row.factor[1+W2@i]*D$col.factor[1+W2@j])
                D$col.factor <- D$col.factor * cL
            }
        }
        W2@x <- a2/(D$row.factor[1+W2@i]*D$col.factor[1+W2@j])
    }
    D$row.factor <- sqrt(D$row.factor)
    D$col.factor <- sqrt(D$col.factor)
    D$iterations <- iter
    D
}

###
## sigda():  Perform a scale invariant geometric analysis of a data matrix <A>,
#  returning <dimensions> scale invariant dimensions.
#
sigda <- function(
    A,                # m x n real-valued data matrix, preferably with row and column labels
    dimensions=NULL,  # How many dimensions to compute; defaults to min(dim(A))
    D.start = NULL,   # Starting projective decomposition, if available
    tol = 1e-6,       # goal for coefficient of variance
    ...)              # additional parameters will be passed to projective.decomposition()
{
    m <- dim(A)[1]
    n <- dim(A)[2]
    max.dimensions <- min(m,n) - 1
    if (is.null(dimensions)) {
        dimensions <- max.dimensions
    }
    K <- 1+dimensions

    # Step 1:  Project the data into an embedded projective space
    W <- as(A, "dgTMatrix")
    if ((is.null(D.start)) || (tol < D.start$tol)) {
        D <- projective.decomposition(W, tol=tol, D.start = D.start, ...)
    } else {
        D <- D.start
    }
    W@x <- W@x / (D$scalar * D$col.factor[1+W@j] * D$row.factor[1+W@i])

    # Step 2:  Singular value decomposition in the embedding space
    if (dimensions + 3 < max.dimensions) { # Use a fast sparse algorithm
        svd <- irlba(as(W,"dgCMatrix"), nu=K, nv=K)
    } else { # use R's standard algorithm
        svd <- svd(W, nu=K, nv=K)
    }
    rm(W)
    # heuristic preferring a right-handed coordinate system
    UV    <- reorient(svd$u,svd$v,K)
    svd$u <- UV$U
    svd$v <- UV$V

    # Prepare the structure of the object to return
    if (! is.null(dimnames(A)[[1]])) {
        rowNames <- dimnames(A)[[1]]
    } else {
        rowNames <- paste("r",c(1:m),sep='')
    }
    if (! is.null(dimnames(A)[[2]])) {
        colNames <- dimnames(A)[[2]]
    } else {
        colNames <- paste("c",c(1:n),sep='')
    }
    dims      <- c(2:K)
    dim.names <- paste('EC',c(1:dimensions),sep='')
    est.mean  <- ( (sum(svd$d[dims]) + (max.dimensions - length(dims))*(min(svd$d)/2))
                 / max.dimensions)
    sig <- list(
        'D' = D,
    'singular' = list(
        'scale' = est.mean,
        'pole'  = svd$d[1] / est.mean,
        'value' = svd$d[dims]) / est.mean,

    'row.scale' = rowSums(A)/n,
    'row.pole'  = svd$u[,1],
    'row'       = svd$u[,-1:0],

    'col.scale' = colSums(A)/m,
    'col.pole'  = svd$v[,1],
    'col'       = svd$v[,-1:0]
    )
    dimnames(sig$row)[[1]] <- rowNames
    dimnames(sig$col)[[1]] <- colNames
    names(sig$singular$value) <- dim.names
    dimnames(sig$row)[[2]]    <- dim.names
    dimnames(sig$col)[[2]]    <- dim.names

    ###
    # Step 3: Projection to Espalier Coordinates (Spherical --> Euclidean)
    ###
    # This step includes accounting for the relative scaling of the axes in the
    # row vector and column vector spaces, but now along the EC coordinate axes.
    ###
    if (m > n) { # columns are higher-dimensional than rows
      for (d in dim.names) {
        sig$row[,d] <- ( (sig$singular$value[d] * sig$row[,d]
                       / (sig$singular$pole     * sig$row.pole)))
        sig$col[,d] <- sig$col[,d] / sig$col.pole
      }
    } else if (n > m) { # rows are higher-dimensional than columns
      for (d in dim.names) {
        sig$row[,d] <- sig$row[,d] / sig$row.pole
        sig$col[,d] <- ( (sig$singular$value[d] * sig$col[,d])
                       / (sig$singular$pole     * sig$col.pole))
      }
    } else { # m == n; rows, columns have the same dimensionality
      for (d in dim.names) {
        sig$row[,d] <- sig$row[,d] / sig$row.pole
        sig$col[,d] <- sig$col[,d] / sig$col.pole
      }
    }
    sig
}

###
# Basic math:
###
# M = { m_ij } = { a r_i w_ij c_j } = a Diag{r} W Diag{c} (Projective Decomposition)
# W = { w_ij } = sum(k: u_ik d_k v_jk) = U Diag{d} t(V)   (Singular Value Decomposition)
#
# U and V are orthogonal (i.e. rotation) matrices, for which
#     inverse(U) = t(U) = transpose(U),
#     inverse(V) = t(V) = transpose(V).
# We can therefore solve for S and V:
#     S  = W Diag{1/d} V,    s_ik = sum(j: w_ij v_jk) / d_k
#   t(V) = Diag{1/d} t(S) W, v_jk = sum(i: s_ik w_ij) / d_k
#
# singular.scale    = E[d_*]
# singular.pole     = d_1     / singular.scale
# singular.value[k] = d_{k+1} / singular.scale
#
# row.pole[i]  = u_i1
# row.scale[i] = mean(m_i*)
# row[i,k]     = u_i{k+1} / u_i1 { * something[k] if n > m }
#
# col.pole[j]  = v_j1
# col.scale[j] = mean(m_*j)
# col[j,k]     = v_j{k+1} / v_j1 { * something[k] if m > n }
#
###

###-----------------------------------------------------------------------------
### High-level plotting functions
###-----------------------------------------------------------------------------

# Plot two scale invariant dimensions against each other.
sigda.biplot <- function(sig, x, y,
                         row.pch=1, row.size=1, row.color = 'gray',
                         row.labels=NULL, row.pos = NULL, row.cex=1,
                         col.pch=2, col.size=1, col.color = NULL,
                         col.labels=NULL, col.pos = NULL, col.cex=1,
                         xlim=NULL, ylim=NULL, zoom=1, ...)
{
    double.biplot(sig$row, sig$col, x, y, type='EC',
                  row.pch=row.pch, row.size=row.size, row.color=row.color,
                  row.labels=row.labels, row.pos=row.pos, row.cex=row.cex,
                  col.pch=col.pch, col.size=col.size, col.color=col.color,
                  col.labels=col.labels, col.pos=col.pos, col.cex=col.cex,
                  xlim=xlim, ylim=ylim, zoom=zoom, ...)
}

# Plot two scale invariant dimensions against each other; row points only.
sigda.row.biplot <- function(sig, x, y,
                             row.pch=1, row.color = 'gray', row.size=1,
                             row.labels=NULL, row.cex=1, row.pos = NULL,
                             pch=NULL, color=NULL, size=NULL, labels=NULL, pos=NULL,
                             zoom=1, ...)
{
    if (is.null(pch))    { pch    <- row.pch }
    if (is.null(color))  { color  <- row.color }
    if (is.null(size))   { size   <- row.size }
    if (is.null(labels)) { labels <- row.labels }
    if (is.null(pos))    { pos    <- row.pos }
    single.biplot(v=sig$row, x=x, y=y,
                  pch=pch, size=size, color=color,
                  labels=labels, pos=pos, single.cex=row.cex,
                  zoom=zoom, ...)
}

# Plot two scale invariant dimensions against each other; column points only.
sigda.col.biplot <- function(sig, x, y,
                             col.pch=2, col.color = 'gray', col.size=1,
                             col.labels=NULL, col.cex=1, col.pos = NULL,
                             pch=NULL, color=NULL, size=NULL, labels=NULL, pos=NULL,
                             zoom=1, ...)
{
    if (is.null(pch))    { pch    <- col.pch }
    if (is.null(color))  { color  <- col.color }
    if (is.null(size))   { size   <- col.size }
    if (is.null(labels)) { labels <- col.labels }
    if (is.null(pos))    { pos    <- col.pos }
    single.biplot(v=sig$col, x=x, y=y,
                  pch=pch, size=size, color=color,
                  labels=labels, pos=pos, single.cex=col.cex,
                  zoom=zoom, ...)
}

# Plot a scale invariant dimension (horizontal axis) against scale (vertical axis).
# The scale invariant dimensions are called an "Espalier Coordinate" due to
# this kind of plot.
espalier.plot <- function(sig, ec,
    row.pch = 1, row.size = 1, row.color = 'gray', row.dither=FALSE,
    row.labels = NULL, row.pos = NULL, row.cex=1,
    row.height=TRUE,
    col.pch = 2, col.size = 1, col.color = 'black', col.dither=FALSE,
    col.labels = NULL, col.pos = NULL, col.cex=1,
    col.height=TRUE,
    y.scale = 1, x.zoom=1, y.zoom=1,
    xlim = NULL, ylim = NULL,
    type='EC', ylab='scale', add=FALSE,
    ...)
{
    # Horizontal axis orientation
    i0 <- abs(ec)
    si <- ec/i0

    # Horizontal axis labels
    xtype <- paste(type,i0,sep='')
    i0 <- match(xtype,dimnames(sig$col)[[2]])
    rows <- dim(sig$row)[1]
    columns  <- dim(sig$col)[1]

    # Vertical dithering
    if (row.dither) {
        dh.row <- runif(rows) - 0.5
        dha.max <- 1
    } else {
        dh.row <- rep(0,rows)
        dha.max <- 0
    }
    if (col.dither) {
        dh.col <- runif(columns) - 0.5
        dho.max <- 1
    } else {
        dh.col <- rep(0,columns)
        dho.max <- 0
    }
    # Plot dimensions
    if (is.null(xlim)) {
        xr <- range(si*c(sig$row[,i0],sig$col[,i0]))/x.zoom
    } else {
        xr <- xlim
    }
    if (row.height) {
        y.row <- y.scale*(columns*sig$row.scale + dh.row)/(columns+dha.max)
    } else {
        min.scale <- min(c(sig$row.scale,sig$col.scale))
        y.row <- visibly.lower(y.scale*(columns*min.scale + dh.row)/(columns+dha.max))
    }
    if (col.height) {
        y.col <- y.scale*(rows*sig$col.scale + dh.col)/(rows+dho.max)
    } else {
        min.scale <- min(c(sig$row.scale,sig$col.scale))
        y.col <- visibly.lower(y.scale*(rows*min.scale + dh.col)/(rows+dho.max))
    }
    if (is.null(ylim)) {
        yr <- range(c(y.row,y.col))*c(y.zoom,1/y.zoom)
    } else {
        yr <- ylim
    }
    if (! add) {
        plot(x = si * sig$row[,i0], y = y.row, type = 'n',
             xlab=ifelse(si < 0, paste('-',xtype, sep=''), xtype), ylab=ylab,
             xlim=xr,  ylim=yr, ...)
        abline(v=0,lty=3)
    }

    # Rows
    points(si*sig$row[,i0], y.row, pch=row.pch, cex=row.size, col=row.color)
    if ((! is.null(row.labels)) && (FALSE != row.labels)) {
        if (rows != length(row.labels)) {
          row.labels <- dimnames(sig$row)[[1]]
        }
        if (is.null(row.pos)) {
          row.pos <- rep(1,rows)
        }
        text(x = si*sig$row[,i0], y = y.row, labels = row.labels,
             pos = row.pos, cex=row.cex, col = row.color)
    }

    # Columns
    points(si*sig$col[,i0], y.col, pch=col.pch, cex=col.size, col=col.color)
    if ((! is.null(col.labels)) && (FALSE != col.labels)) {
        if (columns != length(col.labels)) {
          col.labels <- dimnames(sig$col)[[1]]
        }
        if (is.null(col.pos)) {
          col.pos <- rep(2,columns)
        }
        text(x = si*sig$col[,i0], y = y.col, labels = col.labels,
             pos = col.pos, cex = col.cex, col = col.color)
    }
}

visibly.lower <- function(x) {
  exp(14*log(x)/12)
}

# Plot a scale invariant dimension (horizontal axis) against scale (vertical axis).
# The scale invariant dimensions are called an "Espalier Coordinate" due to
# this kind of plot.
espalier.row.plot <- function(sig, ec,
                          row.pch = 1, row.size = 1, row.color = 'gray', row.dither=FALSE,
                          row.labels = NULL, row.pos = NULL, row.cex=1,
                          row.height=TRUE,
                          y.scale = 1, x.zoom=1, y.zoom=1,
                          xlim = NULL, ylim = NULL,
                          type='EC', ylab='scale', add=FALSE,
                          ...)
{
  # Horizontal axis orientation
  i0 <- abs(ec)
  si <- ec/i0
  
  # Horizontal axis labels
  xtype <- paste(type,i0,sep='')
  i0 <- match(xtype,dimnames(sig$col)[[2]])
  rows    <- dim(sig$row)[[1]]
  columns <- dim(sig$col)[[1]]
  
  # Vertical dithering
  if (row.dither) {
    dh.row <- runif(rows) - 0.5
    dha.max <- 1
  } else {
    dh.row <- rep(0,rows)
    dha.max <- 0
  }

  # Plot dimensions
  if (is.null(xlim)) {
    xr <- range(si*sig$row[,i0])/x.zoom
  } else {
    xr <- xlim
  }
  if (row.height) {
    y.row <- y.scale*(columns*sig$row.scale + dh.row)/(columns+dha.max)
  } else {
    min.scale <- min(c(sig$row.scale,sig$col.scale))
    y.row <- y.scale*visibly.lower((columns*min.scale + dh.row)/(columns+dha.max))
  }
  if (is.null(ylim)) {
    yr <- range(y.row)*c(y.zoom,1/y.zoom)
  } else {
    yr <- ylim
  }
  if (! add) {
    plot(x = si * sig$row[,i0], y = y.row, type = 'n',
         xlab=ifelse(si < 0, paste('-',xtype, sep=''), xtype), ylab=ylab,
         xlim=xr,  ylim=yr, ...)
    abline(v=0,lty=3)
  }
  
  # Rows
  points(si*sig$row[,i0], y.row, pch=row.pch, cex=row.size, col=row.color)
  if ((! is.null(row.labels)) && (FALSE != row.labels)) {
    if (rows != length(row.labels)) {
      row.labels <- dimnames(sig$row)[[1]]
    }
    if (is.null(row.pos)) {
      row.pos <- rep(1,rows)
    }
    text(x = si*sig$row[,i0], y = y.row, labels = row.labels,
         pos = row.pos, cex=row.cex, col = row.color)
  }
}

# Plot a scale invariant dimension (horizontal axis) against scale (vertical axis).
# The scale invariant dimensions are called an "Espalier Coordinate" due to
# this kind of plot.
espalier.col.plot <- function(sig, ec,
                          col.pch = 2, col.size = 1, col.color = 'black', col.dither=FALSE,
                          col.labels = NULL, col.pos = NULL, col.cex=1,
                          col.height=TRUE,
                          y.scale = 1, x.zoom=1, y.zoom=1,
                          xlim = NULL, ylim = NULL,
                          type='EC', ylab='scale', add=FALSE,
                          ...)
{
  # Horizontal axis orientation
  i0 <- abs(ec)
  si <- ec/i0
  
  # Horizontal axis labels
  xtype <- paste(type,i0,sep='')
  i0 <- match(xtype,dimnames(sig$col)[[2]])
  rows <- dim(sig$row)[1]
  columns  <- dim(sig$col)[1]
  
  # Vertical dithering
  if (row.dither) {
    dh.row <- runif(rows) - 0.5
    dha.max <- 1
  } else {
    dh.row <- rep(0,rows)
    dha.max <- 0
  }
  if (col.dither) {
    dh.col <- runif(columns) - 0.5
    dho.max <- 1
  } else {
    dh.col <- rep(0,columns)
    dho.max <- 0
  }
  # Plot dimensions
  if (is.null(xlim)) {
    xr <- range(si*c(sig$row[,i0],sig$col[,i0]))/x.zoom
  } else {
    xr <- xlim
  }
  if (row.height) {
    y.row <- y.scale*(columns*sig$row.scale + dh.row)/(columns+dha.max)
  } else {
    min.scale <- min(c(sig$row.scale,sig$col.scale))
    y.row <- y.scale*(columns*min.scale + dh.row)/(columns+dha.max)
  }
  if (col.height) {
    y.col <- y.scale*(rows*sig$col.scale + dh.col)/(rows+dho.max)
  } else {
    min.scale <- min(c(sig$row.scale,sig$col.scale))
    y.col <- y.scale*(rows*min.scale + dh.col)/(rows+dho.max)
  }
  if (is.null(ylim)) {
    yr <- range(c(y.row,y.col))*c(y.zoom,1/y.zoom)
  } else {
    yr <- ylim
  }
  if (! add) {
    plot(x = si * sig$row[,i0], y = y.row, type = 'n',
         xlab=ifelse(si < 0, paste('-',xtype, sep=''), xtype), ylab=ylab,
         xlim=xr,  ylim=yr, ...)
    abline(v=0,lty=3)
  }
  
  # Rows
  points(si*sig$row[,i0], y.row, pch=row.pch, cex=row.size, col=row.color)
  if ((! is.null(row.labels)) && (FALSE != row.labels)) {
    if (rows != length(row.labels)) {
      row.labels <- dimnames(sig$row)[[1]]
    }
    if (is.null(row.pos)) {
      row.pos <- rep(1,rows)
    }
    text(x = si*sig$row[,i0], y = y.row, labels = row.labels,
         pos = row.pos, cex=row.cex, col = row.color)
  }
  
  # Columns
  points(si*sig$col[,i0], y.col, pch=col.pch, cex=col.size, col=col.color)
  if ((! is.null(col.labels)) && (FALSE != col.labels)) {
    if (columns != length(col.labels)) {
      col.labels <- dimnames(sig$col)[[1]]
    }
    if (is.null(col.pos)) {
      col.pos <- rep(2,columns)
    }
    text(x = si*sig$col[,i0], y = y.col, labels = col.labels,
         pos = col.pos, cex = col.cex, col = col.color)
  }
}

###-----------------------------------------------------------------------------
### Low-level plotting functions
###-----------------------------------------------------------------------------
### These lower-level interfaces are provided for special purposes,
### e.g. plotting PCA results in the same manner as SIGDA results.

# Plot one set of points (row or column) on two scale invariant dimensions
single.biplot <- function(v, x, y, type='EC',
    pch=3, size=1, color = 'black',
    labels=NULL, pos = NULL, single.cex=1, arrows=FALSE,
    xlim=NULL, ylim=NULL, zoom=1, add=FALSE, ...)
{
    # Axis orientations
    i0 <- abs(x)
    j0 <- abs(y)
    si <- x/i0
    sj <- y/j0

    # Axis labels
    xtype <- paste(type,i0,sep='')
    ytype <- paste(type,j0,sep='')
    i0 <- match(xtype,dimnames(v)[[2]])
    j0 <- match(ytype,dimnames(v)[[2]])
    xtype <- ifelse(si >= 0,xtype,paste('-',xtype,sep=''))
    ytype <- ifelse(sj >= 0,ytype,paste('-',ytype,sep=''))
    n     <- dim(v)[1]

    # Plot dimensions
    if (is.null(xlim)) {
        xr <- range(si*v[,i0])/zoom
    } else {
        xr <- xlim
    }
    if (is.null(ylim)) {
        yr <- range(sj*v[,j0])/zoom
    } else {
        yr <- ylim
    }

    # Default colors
    if (is.null(color)) {
        cr <- order(order(atan2(si*v[,i0],sj*v[,j0])))-1
        color <- hsv((2/3)*cr/n,1,2/3)
    }

    # Drawing
    if (! add) {
        plot(x = si*v[,i0], y = sj*v[,j0], type='n',
             xlim=xr, ylim=yr, xlab = xtype, ylab = ytype, ...)
        abline(h=0,lty=3)
        abline(v=0,lty=3)
    }

    # Arrows or points
    if (arrows) {
        arrows(rep(0,n),rep(0,n),
               si*v[,i0], sj*v[,j0], lwd=size,
               length=1/10, angle=30, col = color)
    } else {
        points(x = si*v[,i0], y = sj*v[,j0],
               pch=pch, cex = size, col = color)
    }

    # Labels
    if ((! is.null(labels)) && (labels != FALSE)) {
        if (n != length(labels)) {
            labels <- dimnames(v)[[1]]
        }
        if (is.null(pos)) {
            cr <- order(order(atan2(si*v[,i0],sj*v[,j0])))-1
            pos <- 1 + (trunc(4*cr/n) %% 4)
        }
        text(x = si*v[,i0], y = sj*v[,j0], labels = labels,
             pos = pos, cex=single.cex, col = color)
    }
}

# Plot both row and column points on two scale invariant dimensions
double.biplot <- function(row.v, col.v, x, y, type='EC',
                      row.pch=1, row.size=1, row.color = 'gray',
                      row.labels=NULL, row.pos = NULL, row.cex=1,
                      col.pch=2, col.size=1, col.color = 'black',
                      col.labels=NULL, col.pos = NULL, col.cex=1,
                      row.arrows=FALSE, col.arrows=FALSE,
                      xlim=NULL, ylim=NULL, zoom=1, add=FALSE, ...)
{
    # Axis orientations
    i0 <- abs(x)
    j0 <- abs(y)
    si <- x/i0
    sj <- y/j0

    # Axis labels
    xtype <- paste(type,i0,sep='')
    ytype <- paste(type,j0,sep='')
    i0 <- match(xtype,dimnames(col.v)[[2]])
    j0 <- match(ytype,dimnames(col.v)[[2]])
    xtype <- ifelse(si >= 0,xtype,paste('-',xtype,sep=''))
    ytype <- ifelse(sj >= 0,ytype,paste('-',ytype,sep=''))
    rows    <- dim(row.v)[1]
    columns <- dim(col.v)[1]

    # Plot dimensions
    if (is.null(xlim)) {
        xr <- range(si*c(row.v[,i0],col.v[,i0]))/zoom
    } else {
        xr <- xlim
    }
    if (is.null(ylim)) {
        yr <- range(sj*c(row.v[,j0],col.v[,j0]))/zoom
    } else {
        yr <- ylim
    }
    if (! add) {
        plot(x = si*row.v[,i0], y = sj*row.v[,j0], type='n',
             xlim=xr, ylim=yr, xlab = xtype, ylab = ytype, ...)
        abline(h=0,lty=3)
        abline(v=0,lty=3)
    }

    # Rows
    if (row.arrows) {
        arrows(rep(0,rows),rep(0,rows),
               si*row.v[,i0], sj*row.v[,j0], lwd=row.size,
               length=1/10, angle=30, col = row.color)
    } else {
        points(x = si*row.v[,i0], y = sj*row.v[,j0],
               pch=row.pch, col=row.color, cex=row.size)
    }
    if ((! is.null(row.labels)) && (row.labels != FALSE)) {
        if (rows != length(row.labels)) {
            row.labels <- dimnames(row.v)[[1]]
        }
        if (is.null(row.pos)) {
            row.round <- order(order(atan2(si*row.v[,i0],sj*row.v[,j0])))-1
            row.pos   <- 1 + (round(4*row.round/rows) %% 4)
        }
        text(x = si*row.v[,i0], y = sj*row.v[,j0], labels = row.labels,
             pos = row.pos, cex=row.cex, col = row.color)
    }

  # Columns
    if (is.null(col.color)) {
        col.round <- order(order(atan2(si*col.v[,i0],sj*col.v[,j0])))-1
        col.color <- hsv((2/3)*col.round/columns,1,2/3)
    }
    if (col.arrows) {
        arrows(rep(0,columns),rep(0,columns),
               si*col.v[,i0], sj*col.v[,j0], lwd=col.size,
               length=1/10, angle=30, col = col.color)
    } else {
        points(x = si*col.v[,i0], y = sj*col.v[,j0],
               pch=col.pch, cex = col.size, col = col.color)
    }
    if ((! is.null(col.labels)) && (col.labels != FALSE)) {
        if (columns != length(col.labels)) {
            col.labels <- dimnames(col.v)[[1]]
        }
        if (is.null(col.pos)) {
            col.round <- order(order(atan2(si*col.v[,i0],sj*col.v[,j0])))-1
            col.pos <- 1 + (trunc(4*col.round/columns) %% 4)
        }
        text(x = si*col.v[,i0], y = sj*col.v[,j0], labels = col.labels,
             pos = col.pos, cex=col.cex, col = col.color)
    }
}

###-----------------------------------------------------------------------------
### Utility functions
###-----------------------------------------------------------------------------

# Coefficient of variation of a vector
cv <- function(x) { sqrt(var(x))/mean(x) }

rms <- function(x) { sqrt(mean(x*x)) }

# Matrix norm
Frobenius <- function(A) {
    if (   ('dgTMatrix' %in% class(A))
        || ('dgCMatrix' %in% class(A))) {
        sqrt(sum(A@x*A@x))
    } else {
        sqrt(sum(as.vector(A)*as.vector(A)))
    }
}

###
## reorient()
#
# The relative orientation of vectors in N-dimensional space is
# difficult to work with for N > 3.  Here we want the basis
# vectors produced by SVD to have the same chirality (in the new
# order produced by SVD) across different runs and even different
# data.  We approach this problem by
# (1) ensuring that the first basis vector is oriented more positive
#     than negative (i.e. has a nonnegative projection onto 1_k,
#     where k = min{m,n}), and
# (2) considering 1_k to be the "zeroth basis vector, we attempt
#     to orient each consecutive triplet of basis vectors into
#     the same (right-handed) chirality.
# Note that the singular vector pairs are fully corresponding
# directions, so the sign of one determines the sign of other.
#
# If the basis vectors were of length 3 (if k = 3), we could
# simply take the determinant of the 3x3 matrix of three basis
# vectors and switch the sign of the third basis vector whenever
# the determinant is negative; however for k > 3, the
# matrix is not square and therefore does not have a
# determinant.  As a heuristic, then, we combine entries
# in the matrix into a 3 x 3 matrix, decimating the
# longer side modulo 3, and then compute the determinant.
# This appears to have the desired effect.
###
reorient <- function(U,V,K)
{
    m <- dim(U)[[1]]
    n <- dim(V)[[1]]
    if (m <= n) {
        M <- matrix(c(1,1,1,
	              0,1,1,
		      0,0,1),
		    nrow=3,byrow=TRUE)
	Is <- (c(1:m)-1) %% 3
	UV <- list('U' = U, 'V' = V)
	for (k in c(1:K)) {
	    j <- 1 + ((k+1) %% 3)
	    for (i in c(0:2)) {
	        M[1+i,j] <- sum(U[Is == i,k])
	    }
	    if (determinant(M)$sign < 0) {
	        UV$U[,k] <- -UV$U[,k]
	        UV$V[,k] <- -UV$V[,k]
		M[,j]    <- -M[,j]
	    }
	}
    } else {
        VU <- reorient(V,U,K)
	UV <- list('U' = VU$V, 'V' = VU$U, 'L' = VU$L)
    }
    UV
}

# end of sigda.R
