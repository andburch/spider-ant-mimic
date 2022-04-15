#==================
# "Private"

.JRunExpr <- function(expr) {
  if (is.function(expr))
    expr()
  else
    eval(expr)
}

.JplotToDevice <- function(filename, plotFn, onlyIfDoesntExist, openDeviceFn) {
  if (!onlyIfDoesntExist || !file.exists(filename)) {
    openDeviceFn()
    tryCatch({
      if (is.function(plotFn))
        plotFn()
      else
        invisible(eval(plotFn))
    }, finally = {
      dev.off()
    })
  }
}
# Default resolution seems to be 72, so to increase image pixel size without decreasing text size,
# line width etc, increase resolution accordingly.
# Doesn't work with ggplot since it must be evaluated to work.
# Try using ggsave or print(plotFn())
PlotToPng <- function(filename, plotFn, width=800, height=600, res = NA, onlyIfDoesntExist = F) {
  .JplotToDevice(filename, plotFn, onlyIfDoesntExist, function () {
    # type = 'cairo' seems to produce _much_ nicer graphics with antialiasing
    png(filename, width = width, height = height, type = 'cairo', res = res)
  })
}

# Beware: width and height are in inches, not pixels!
PlotToPDF <- function(filename, plotFn, width=8, height=6, onlyIfDoesntExist = F) {
  .JplotToDevice(filename, plotFn, onlyIfDoesntExist, function () {
    pdf(filename, width = width, height = height)
  })
}

WriteToFile <- function(filename, expr) {
  sink(filename)
  on.exit(sink())
  
  .JRunExpr(expr)
}

# Plots a list of densities as lines
PlotDensities <- function(densities, cols, lty = 1, lwd = 2, add = FALSE, xlim = NULL, ylim = NULL, includeInX = numeric(0), ...) {
  
  # Create empty plot
  if (!add) {
    if (is.null(xlim))
      xlim <- range(lapply(densities, function(d) d$x), na.rm = TRUE)
    xlim <- range(c(xlim, includeInX))
    if (is.null(ylim))
      ylim <- range(lapply(densities, function(d) d$y), na.rm = TRUE)
    plot(NULL, xlim = xlim, ylim = ylim, ...)
  }
  
  # Recycle lty if it's a single number
  if (length(lty) == 1)
    lty <- rep(lty, length.out = length(densities))
  
  # Plot densities as lines
  i <- 1
  for(d in densities) {
    lines(d, col = cols[i], lty = lty[i], lwd = lwd)
    i <- i + 1
  }
}

# Sets some options, runs an expression, then restores the original option values
TemporarilySetOptions <- function(expr, ...) {
  oldOptions <- options(...)
  # Restore original options when the function exits
  on.exit(options(oldOptions))
  
  .JRunExpr(expr)
}

# Captialises the first word in a single string
# Copied (and changed) from the doc for toupper
capsentence <- function(s, strict = FALSE) {
  s <- as.character(s)
  paste(toupper(substring(s, 1, 1)),
        {s <- substring(s, 2); if(strict) tolower(s) else s},
        sep = "", collapse = " " )
}

capwords <- function(s, strict = FALSE) {
  s <- as.character(s)
  sapply(strsplit(s, split = " "), capsentence, USE.NAMES = !is.null(names(s)))
}

# Italicise parts of strings for display, e.g. in a plot legend.
# Only 1 section of each string will be italicised.
italicisePattern <- function(str, pattern) {
  
  # Italicise 1 part of a single string
  .italiciseBit <- function(str, start, length) {
    if (start == -1) return(str)
    
    pre <- substr(str, 1, start - 1)
    it <- substr(str, start, start + length)
    suf <- substr(str, start + length + 1, length(str))
    # That old black magic has me in its spell
    as.expression(bquote(.(pre) * italic(.(it)) * .(suf)))
  }
  
  # Evaluate the regexp pattern against all the strings
  matches <- regexec(pattern, str)
  # For each string...
  sapply(1:length(matches), function(i) {
    # Replace the string with an expression with the matching substring italicised
    m <- matches[[i]]
    .italiciseBit(str[i], m[1], attr(m, "match.length"))
  })
}

# Adds transparency to a colour, as a value between 0 (completely transparent) and 255 (not transparent)
addTransparency <- function(colour, alpha) {
  c <- col2rgb(colour)
  rgb(c[1,], c[2,], c[3,], alpha, maxColorValue = 255)
}

# Runs and reports dunns test, but outputs the results differently than the dunns.test function - 
# full row/column names are reported in the order specified by the order argument.
# Note that the results are printed rather than returned - this is due to the quirky way
# that dunn.test is implemented. It only prints out some of the kruskal wallis results 
# rather than returning them.
customDunnsTest <- function(params, stat, order, testGroups = params$species, alpha = 0.05, adjustment = NULL) {
  
  # Bizarrely, dunn.test prints but doesn't return the results of the
  # kruskal-wallis test (other than chi-squared), so we still need some output
  # from the function
  dn <- dunn.test::dunn.test(params[, stat], g = testGroups, method = adjustment, wrap = TRUE, table = FALSE, alpha = alpha)
  cat("\n")
  
  nsp <- length(order)
  cols <- order[1:(nsp - 1)]
  rows <- order[2:nsp]
  
  # Build a results table containing rows and columns in the order they should be output
  tb <- do.call(cbind, lapply(cols, function(col) {
    # Get all pairs of rows and columns in the format used by dn$comparisons, "a - b"
    cellsFwd <- paste(col, rows, sep = " - ")
    # Find indices of matches
    idxs <- match(cellsFwd, dn$comparisons)
    # Find missing indices from forward pairs in the reverse pairs
    missingIdxs <- which(is.na(idxs))
    if (length(missingIdxs) > 0) {
      # Try to get missing comparisons using "b - a"
      cellsRev <- paste(rows, col, sep = " - ")
      matchRev <- match(cellsRev, dn$comparisons)
      idxs[missingIdxs] <- matchRev[missingIdxs]
    }
    
    # Now we have the indices of the comparisons for "col" and all rows
    # Create a column with (Z, P, " ") for each species row.
    # Round Z to 6 places, P to 4.
    Z <- round(dn$Z[idxs], digits = 6)
    P <- round(dn$P.adjusted[idxs], digits = 4)
    P <- ifelse(is.na(P), P, paste0(P, ifelse(P < alpha / 2, "*", " ")))
    # The tricky c(rbind()) thing interleaves vectors
    c(rbind(Z, P, rep('', length(idxs))))
  }))
  tb[is.na(tb)] <- ''
  # Clear out the upper triangle
  tb[row(tb) <= (col(tb) - 1) * 3] <- ''
  tb <- as.data.frame(tb)
  tb <- cbind(c(rbind(paste(rows, '|'), '|', '|')), tb)
  colnames(tb) <- c('', cols)
  print(tb, row.names = FALSE)
}
