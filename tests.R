#install.packages("tibble")
library(tibble)

manova_analysis <- function(params, categories, make_file_name, cat_order, alpha=0.05){
  manova_results <- manova(as.matrix(params) ~ as.factor(categories))
  WriteToFile(make_file_name("data_output/%s_manova.txt"), 
              print(summary(manova_results)))
  
  for(param_name in names(params)){
    transform <-  PARAMS_INFO[[param_name]]$transform
    if(!is.null(transform))
      univariate_analysis(params[, param_name], param_name, transform, categories, cat_order, make_file_name, alpha=alpha)
  }
}

# Apply reverse transform 
reverse_transform <- function(val, transform) {
  if(transform == "log")
    exp(val)
  else if (transform == "sqrt")
    val ^ 2
  else if (transform == "")
    val
  else
    stop(sprintf("Unknown data transform '%s'", transform))
}

# Apply ANOVA followed by Tukey-Kramer and report the result
univariate_analysis <- function(param, param_name, transform, categories, cat_order, make_file_name, alpha = 0.05){
  cat_fac <- as.factor(categories)
  anova <- aov(param ~ cat_fac)
  means <- aggregate(reverse_transform(param, transform), list(cat_fac), mean)
  tukey <- TukeyHSD(anova, conf.level = 1 - alpha)
  WriteToFile(make_file_name(sprintf("data_output/%%s_%s_anova.txt", param_name)),
              TemporarilySetOptions({
                cat(sprintf("ANOVA followed by Tukey-Kramer tests for %s%s\n\n", 
                            ifelse(transform == "", "(un-transformed) ", paste0(transform, "-transformed ")), param_name))
                print(summary(anova))
                cat("\n\n")
                cat(sprintf("Tukey multiple comparisons of means \n  %g%% family-wise confidence level of (column - row)%s\n", 
                            100 * (1 - alpha),
                            ifelse(transform == "", "", sprintf(" (%s-transformed data)", transform))))
                cat("  Adjusted p-value\n")
                cat("  Relative sample means (untransformed) as percent (column / row)\n\n")
                print(report_tukey(tukey, transform, means, categories, cat_order, alpha), row.names = FALSE)
                
                # cat("\n\nSample means:\n")
                # print(means, row.names = FALSE)
              }, width = 1000)
  )
}

# Formats the results of a Tukey test in a more meaningful way than the default
report_tukey <- function(tukey, transform, means, categories, cat_order, alpha = 0.05) {
  # Separate out the category names from the paired names
  tukey <- tukey[[1]]
  # This expression is unique to this data set - assume that "-" is only used as
  # a separator or else within "-mimic"
  pairs <- sub("XXX", "-mimic", sub("-", ":", sub("-mimic", "XXX", rownames(tukey))))
  pairs <- strsplit(pairs, ":")
  # Convert ratio to a percent
  valToPC <- function(val, relativeTo) {
    round(100 * val / relativeTo)
  }

  ncats <- length(cat_order)
  cpr <- 4 # Cells per species-row
  mtbl <- matrix("", nrow = cpr * ncats, ncol = ncats)
  
  .addCell <- function(first_name, second_name, diff, lower, upper, first_mean, second_mean, p, ss) {
    c <- match(first_name, cat_order)
    r <- match(second_name, cat_order)
    # Always fill in the lower-left diagonal
    if (r < c) {
      # Swap indices
      tmp <- r
      r <- c
      c <- tmp
      # Negate CI
      tmp <- lower
      lower <- -upper
      upper <- -tmp
      
      diff <- -diff
      tmp <- first_mean
      first_mean <- second_mean
      second_mean <- tmp
    }
    ri <- (r - 1) * cpr + 1
    mtbl[[ri, c]] <<- sprintf("%g,%g", lower, upper)
    mtbl[[ri + 1, c]] <<- paste0(p, ss)
    mtbl[[ri + 2, c]] <<- sprintf("%g%%", valToPC(first_mean, second_mean))
  }
  
  for (i in 1:nrow(tukey)) {
    first <- pairs[[i]][1]
    second <- pairs[[i]][2]
    p <- tukey[i, 4]
    pr <- round(p, 3)
    ss <- ifelse(p < alpha, "*", "")
    first_mean <- means[which(means[, 1] == first), 2]
    second_mean <- means[which(means[, 1] == second), 2]
    .addCell(first, second, tukey[i,1], tukey[i,2], tukey[i,3], first_mean, second_mean, pr, ss)
    #cat(sprintf("%-19s%-19s: diff %g (%g%%)\t 95%% CI %s:%s\tp-value %g%s\n", first, second, tukey[i, 1], valToPC(tukey[i, 1], first_mean), lower, upper, pr, ss))
  }
  
  tbl <- as.tibble(mtbl) # Prevents row names being printed
  tbl <- cbind(c(rbind(cat_order, rep("", ncats), rep("", ncats), rep("", ncats))), tbl)
  #tbl <- as.tibble(tbl)
  colnames(tbl) <- c("", cat_order)
  # Chop off first row (actually 4 rows) and last column
  tbl <- tbl[(cpr+1):nrow(tbl), 1:(ncol(tbl)-1)]
  tbl
}
