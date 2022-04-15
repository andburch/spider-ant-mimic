# Produces a box plot of a single parameter for the various bee/mimic and
# wasp/mimc species
pretty_box_plot <- function(values, species, all_species, param_label, title) {
  species_factor <- factor(species, all_species)
  names <- as.character(levels(species_factor))
  names <- italicisePattern(names, ".*")
  par(mar = c(4, 5, 4, 1))
  yRange <- range(values)
  yRange[2] <- yRange[1] + diff(yRange) * 1.1
  breaks <- pretty(yRange)
  ylim <- range(breaks)
  formula <- reformulate("species_factor", "values")
  plot(formula, names = names, range = 0, 
       xlab = "", ylab = param_label, yaxt = "n",
       xlim = c(0.8, 8.2), ylim = ylim, 
       main = title, cex.axis = 1.1, cex.main = 2, cex.lab = 2)
  # Draw background colours
  y0 <- extendrange(ylim)[1]
  y1 <- extendrange(ylim)[2]
  rect(0.1, y0, 3.5, y1, col = addTransparency(.plot_colours("red stingless bee"), 40))
  rect(3.5, y0, 6.5, y1, col = addTransparency(.plot_colours("bee-mimic P. cruentata"), 40))
  rect(6.5, y0, 7.5, y1, col = addTransparency(.plot_colours("wasp-mimic Pyrophleps"), 40))
  rect(7.5, y0, 9, y1, col = addTransparency(.plot_colours("wasp"), 40))
  # Replot to get white boxes
  plot(formula, names = rep('', length(names)), range = 0, 
       xlab = "", ylab = "", add = TRUE, axes = FALSE, ylim = c(0, 0.8))
  means <- tapply(values, species_factor, mean)
  points(1:length(means), means, pch = 16)
  #axis(2, seq(0, 0.8, by = 0.1), cex.axis = 1.1)
  axis(2, breaks, cex.axis = 1.1)
  padj <- .75
  mtext("Bee mimics", line = -1, at = 2, cex = 1.4, padj = padj)
  mtext("Bees", line = -1, at = 5, cex = 1.4, padj = padj)
  mtext("Wasp\nmimics", line = -1, at = 7, cex = 1.4, padj = padj)
  mtext("Wasps", line = -1, at = 8, cex = 1.4, padj = padj)
}

# Given a category, return the colour to be used to plot it in PCA (types are groupd by coloour)
.colFromCategory <- function(category) {
  ifelse(category == "bee-mimic H. pahangensis", "red",
         ifelse(category == "bee-mimic A. argentifasciata", "red",
                ifelse(category == "bee-mimic P. cruentata", "red",
                       ifelse(category == "dwarf honey bee", "red",
                              ifelse(category == "black stingless bee", "red",
                                     ifelse(category == "red stingless bee", "red",
                                            ifelse(category == "wasp-mimic Pyrophleps", "black",
                                                   "black") # "wasp" 
                                     ))))))
}

# Given a category, return the pch (point symbol) to be used to plot it
.pchFromCategory <- function(category) {
  ifelse(category == "bee-mimic H. pahangensis", 16,
         ifelse(category == "bee-mimic A. argentifasciata", 15,
                ifelse(category == "bee-mimic P. cruentata", 17,
                       ifelse(category == "dwarf honey bee", 0,
                              ifelse(category == "black stingless bee", 1,
                                     ifelse(category == "red stingless bee", 2,
                                            ifelse(category == "wasp-mimic Pyrophleps", 18,
                                                   5) # "wasp" 
                                     ))))))
}

# Given a category, return the colour to be used to plot it in bar plot, 
# each category has a different colour
.plot_colours <- function(category) {
  ifelse(category == "bee-mimic H. pahangensis", "green",
         ifelse(category == "bee-mimic A. argentifasciata", "violetred",
                ifelse(category == "bee-mimic P. cruentata", "orange",
                       ifelse(category == "dwarf honey bee", "brown",
                              ifelse(category == "black stingless bee", "black",
                                     ifelse(category == "red stingless bee", "red1",
                                            ifelse(category == "wasp-mimic Pyrophleps", "blue",
                                                   "turquoise") # "wasp" 
                                     ))))))
}

###

# Reports results to text files and plots to png files
report_result <- function(params, name){

  make_file_name <- function(base) sprintf(base, name)
  
  draw_legend <- function(cats, col, lwd = NA, pch = NULL, position = "topright", cex = 1.5){
    labels <- sapply(cats, capsentence)
    species_pattern <- " [A-Z](\\.|[a-z]+)( [a-z]+)?$"
    labels <- italicisePattern(labels, species_pattern)
    legend(position, labels, col = col, lwd = lwd, pch = pch, inset = .01, cex = cex)
  }
  
  # Get names of parameters to include in PCA
  in_pca <- sapply(PARAMS_INFO, function(pi) pi$in_pca)
  param_names <- names(PARAMS_INFO)[in_pca]
  x_labels <- sapply(PARAMS_INFO[in_pca], function(pi) pi$label)

  # Impute NA values to means for parameters, but don't create a flag column
  pca_able <- TrajsStatsReplaceNAs(params[,param_names], "first_min_deltaS")
  pca_able <- TrajsStatsReplaceNAs(pca_able, "first_min_C")
  number_components <- 5
  PCA <- FactoMineR::PCA(pca_able, ncp = number_components, graph = FALSE)
  
  ## Output various PCA results
  WriteToFile(make_file_name("data_output/%s_PCA.txt"), {
    cat(sprintf("Information about principal component analysis (first %d components)\n\n", number_components))
    print(t(PCA$eig[1:number_components, ]))
    cat("\nCorrelation loadings:\n")
    print(PCA$var$coord)
  })
  
  PlotToPng(make_file_name("plots/%s_PCA_plot.png"), 
            function(){
              # Add a top margin to the plot so that the legend doesn't obscure any points
              ylim <- range(PCA$ind$coord[,2])
              ylim[2] <- extendrange(PCA$ind$coord[,2], f = 0.12)[2]
              plot(PCA$ind$coord[,1:2], 
                   col = .colFromCategory(params$category), pch = .pchFromCategory(params$category), 
                   ylim = ylim,
                   xaxt='n', yaxt='n', xlab = "", ylab = "",
                   main = "Principal Component Analysis", cex.main = 2, cex = 2)
              title(xlab = "PC1", ylab = "PC2", line = 1, cex.lab = 2)
              cats <- allcategories
              draw_legend(cats, col = .colFromCategory(cats), pch = .pchFromCategory(cats), cex = 1.4) 
            },
            width = FIGURE_WIDTH, height = FIGURE_HEIGHT, res = FIGURE_RESOLUTION)
  
  univariate_analysis(PCA$ind$coord[, 1], "PC1", "", params$species, cat_order = all_species, make_file_name)
  
  PlotToPng(make_file_name("plots/%s_PC1_bar_plot.png"), 
            pretty_box_plot(PCA$ind$coord[,1], params$species, all_species, "PC1", "PCA Dimension 1"),
            width = FIGURE_WIDTH, height = FIGURE_HEIGHT, res = FIGURE_RESOLUTION)
  
  # Direction autocorrelation plot
  PlotToPng(make_file_name("plots/%s_scatter_plot_directionauto.png"), 
            function(){ 
              plot(first_min_C ~ first_min_deltaS, data = params, 
                   col = .colFromCategory(params$category), 
                   pch = .pchFromCategory(params$category),
                   ylab = expression("C(" * Delta * s * ")"),
                   xlab = expression(Delta * s))
              cats <- allcategories
              draw_legend(cats, col = .colFromCategory(cats), pch = .pchFromCategory(cats), cex = 1.2)          
              })
  
  # Probability density plots for each (possibly transformed) parameter
  .plotSpeciesDensities <- function(stat, data = params, xlim = NULL, log = "") {
    #A4 page width in inches minus 1 inch margins x 300 dpi
    PlotToPng(make_file_name(sprintf("plots/%%s_plot_%s.png", stat)),
              function(){
                cats <- allcategories
                densities <- lapply(cats, function(sp) {
                  density(data[data$category == sp, stat], na.rm = TRUE)
                })
                col <- .plot_colours(cats)
                title <- capsentence(stat, TRUE)
                title <- gsub("_", " ", title)
                xlab <- x_labels[which(param_names == stat)]
                tr <- PARAMS_INFO[[stat]]$transform
                if (tr != "") {
                  if (tr == "log")
                    tr <- "ln"
                  xlab <- sprintf("%s %s", tr, xlab)
                }
                par(mar = c(5, 5, 4, 1))
                PlotDensities(densities, col, main = title, lwd = 3, xlab = xlab, ylab = "Probability density", 
                              log = log, cex.main = 2, cex.axis = 1.1, cex.lab = 2)
                position <- ifelse(stat %in% c("straightness", "Emax", "mean_speed", "min_speed", "sd_speed"), "topleft", "topright")
                draw_legend(cats, position = position, col = col, lwd = 3)
              },
              width = FIGURE_WIDTH, height = FIGURE_HEIGHT, res = FIGURE_RESOLUTION)
    
  }
  
  for(stat in param_names){
    if(!any(is.na(params[,stat]))){
      .plotSpeciesDensities(stat)
    }
  }
  
  manova_analysis(pca_able, params$species, make_file_name, cat_order = all_species)
}

# Writes a CSV file with all parameters for each trajectory, and with min, mean
# and max for each species and parameter.
report_params <- function(params, name, species_order) {

  make_file_name <- function(base) sprintf(base, name)

  # Write out a complete list of all parameters for each trajectory
  write.csv(params, make_file_name("data_output/%s_parameters_complete.csv"), row.names = FALSE)
  
  # Write out a list of parameter mean, min and max values for each species
  agg_fn <- function(x) { 
    if (any(!is.na(x)))
      c(mean = mean(x, na.rm = TRUE), min = min(x, na.rm = TRUE), max = max(x, na.rm = TRUE))
    else
      c(mean = NA, min = NA, max = NA)
  }
  agg <- aggregate(params[, names(PARAMS_INFO)], by = list(params$species), agg_fn)
  names(agg) <- c("Species", names(agg)[2:ncol(agg)])
  # Sort rows
  agg <- agg[match(species_order, agg$Species),]
  write.csv(agg, make_file_name("data_output/%s_parameters_species_summary.csv"), row.names = FALSE)
  
  # Shows maximums and minimums for each parameter from parameters_complete
  # WriteToFile(make_file_name("data_output/%s_max_min_parameters.txt"),
  #             for(stat in param_names){
  #               cat(sprintf("%s:\tmaximum,\t%g,\t%s\n", stat, params[which.max(params[,stat]), stat], params[which.max(params[,stat]), ]$csv.file.name))
  #               cat(sprintf("%s:\tminimum,\t%g,\t%s\n", stat, params[which.min(params[,stat]), stat], params[which.min(params[,stat]), ]$csv.file.name))
  #             }
  # )
}

# Reports kruskall-wallace and dunn tests for the named parameters
report_dunn <- function(params, param_names, param_info, file_prefix, dunnTestGroups = params$species, dunnGroupName = 'species', dunnTestOrder = all_species, p.adjustment = "holm", alpha = 0.05) {
  kw_statistics <- function(stat, label){
    WriteToFile(sprintf("data_output/%s_%s.txt", file_prefix, stat), 
                TemporarilySetOptions({
                  cat(sprintf("Comparison of %s by %s (%s adjustment)\n\n", label, dunnGroupName, p.adjustment))
                  customDunnsTest(params, stat, order = dunnTestOrder, testGroups = dunnTestGroups, adjustment = p.adjustment, alpha = alpha)
                },
                width = 10000)
    )
  }
  
  for(name in param_names) {
    kw_statistics(name, param_info[[name]]$label)
  }
}
