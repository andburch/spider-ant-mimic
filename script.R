#' ---
#' title: "Spider Ant Mimics"
#' author: "Andrew Burchill"
#' date: "2019-06-19"
#' output:
#'    html_document:
#'       theme: lumen
#'       includes:
#'         after_body: footer.html
#'         before_body: style.html
#'       code_folding: "hide"
#'       toc: true
#'       toc_float: 
#'         collapsed: false
#'         smooth_scroll: true
#'       
#' always_allow_html: yes  
#' ---
#' 


#+ setup, include=FALSE

#___Packages and Functions################
library(trajr)
library(GeneCycle)
library(readr)
library(tidyverse)
library(magrittr)
library(gsubfn)
library(ggfortify)
library(segclust2d)
library(factoextra)

#setwd("C:/Users/Andrew Burchill/Dropbox (ASU)/Ant mimic strobing/Finished")
setwd("D:/Weaver ants/Dropbox (ASU)/Ant mimic strobing/Finished")

#Functions below
{ 
  
  
remove_outliers <- function(x, na.rm = TRUE, ...) {
    qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
    H <- 1.5 * IQR(x, na.rm = na.rm)
    y <- x
    y[x < (qnt[1] - H)] <- NA
    y[x > (qnt[2] + H)] <- NA
    y
  }  
  
cvs_reader_func <- function(filename, stringsAsFactors = FALSE) {
  read_csv(filename,
           col_types = cols_only(frame = col_guess(), 
                                 t = col_guess(), 
                                 x = col_guess(),  
                                 y = col_guess()
           ),
           skip = 1) %>% return()
}

nobu <- function(speed.vector, fraction = 5){
  ### look for ultradian rhythms using this
  
  ################# Modified from code by Nobuaki Mizumoto
  
  da <- speed.vector
  
  N <- length(da)
  xmean <- mean(da)
  how.long <- round(N / fraction) #don't want it to have less than 5 cycles
  
  #qp1per is Qp 1%. for statistical significant. can also use 5%
  P = Qp1per = Qp = rep(0, how.long - 1)    
  for (j in 2:how.long) {
    P[j - 1] <- j  #P is the period
    
    #divides the data into a matrix with j-1 rows
    #I don't want it to keep saying 'data length is not a sub-multiple or multiple of the number of rows'
    suppressWarnings(   
      rhythm.tab <- matrix(da, ncol = P[j - 1], byrow = T)
    )   
    y <- rep(0, P[j - 1])
    for (i in 1:P[j - 1]) {
      y[i] <- mean(rhythm.tab[, i])   #this, like, averages each col?
    }
    K <- length(rhythm.tab[, 1])
    
    Qp[j - 1] <- K * N * sum((y - xmean) ^ 2) / sum((da - xmean) ^ 2)
    Qp1per[j - 1] <- qchisq(0.99, P[j - 1] - 1)
  }
  res <- data.frame(P = P, Qp1per = Qp1per, Qp = Qp)
  
  diff <- Qp - Qp1per
  
  biggest <- P[diff == max(diff)] / 240  #divide by 240 SECONDS fps
  
  #par(pin = c(4, 3))
  # plot(P / 240, Qp, type = "l", xlab = "time (second)")
  # points(P / 240, Qp1per, type = "l", col = 2)
  # abline(v = biggest, col = "blue")
  
  
  return(max(diff))
  
  
  #find top 5% of points
  #perc <- 5
  #which(diff > quantile(diff, prob = 1 - perc/100))
  
  
}

intervaler <- function(traj, greatthan = 0, ps = 1, ns = 11,
                       cutoffspeed = 0.1, display = FALSE, getdbs=FALSE) {
  filename <- attr(traj, "filename")
  
  #calculate the intervals of stopping
  traj %>%
    TrajSmoothSG(p = ps, n = ns) %>% 
    TrajSpeedIntervals(slowerThan = cutoffspeed) -> intervals 
  #remove PAUSES that are not really part of the gait
  intervals<-intervals[intervals$duration<.25,]
  plot(intervals, main = filename)
  
  #calculate the intervals of MOVING (less used)
  traj %>%
    TrajSmoothSG(p = ps, n = ns) %>% 
    TrajSpeedIntervals(fasterThan = cutoffspeed) -> REVERSEintervals 
  
  #print the number of intervals of each type
  if (display) print(
    paste0(
      "stopped intervals: ", length(intervals$duration),
      ". MOVING intervals: ", length(REVERSEintervals$duration)
    ))
  
  if (getdbs==T) {
    list(
      intervals =  intervals,
      REVERSEintervals = REVERSEintervals,
      freqs = intervals$startTime %>% diff() %>% magrittr::extract(. > greatthan)
    ) %>% 
      return()
    
  } else {
  #stop if there aren't enough intervals
  if (nrow(intervals) > 1 & nrow(REVERSEintervals) > 1) {
    
    #get the "frequencies" of stop-start "waves"
    freqs <- intervals$startTime %>% diff() %>% magrittr::extract(. > greatthan)
    
    #then get the durations of stops and starts
    stopdurs <- intervals$duration %>% magrittr::extract(. > greatthan)
    movedurs <- REVERSEintervals$duration %>% magrittr::extract(. > greatthan)
    
    #plot the densities of these durations
    if (display) {
      print(
        ggplot() + labs(fill = "type", title = filename) +
          geom_density(data = data.frame(x = freqs), 
                       aes(x = x, fill = "freq"), alpha = 0.4) + 
          geom_density(data = data.frame(y = stopdurs),
                       aes(x = y, fill = "stop"), alpha = 0.4) +
          geom_density(data = data.frame(z = movedurs),
                       aes(x = z, fill = "move"), alpha = 0.4) +
          xlim(-0.01, 1.5) + 
          scale_fill_manual(
            name = "vars", guide = "legend",
            values = c("freq" = "red", "stop" = "yellow", "move" = "blue"),
            labels = c("freq" = "stop-start duration", "stop" = "stop durations", "move" =  "movement durations")) 
        
      )
    }
    
    combined <- NULL
    #include a thing that makes a DB of durations of stops and starts
    combined <- gdata::combine(intervals, REVERSEintervals)
    combined <- combined[order(combined$startFrame),]
    #remove first and last breaks (they could be cut off)
    combined <- combined %>% slice(-n()) %>% slice(-1)
    combined <- combined %>% select(duration, source) %>% mutate(filename = filename)
    
    # INTERVAL.DB <<- rbind(INTERVAL.DB, combined)
    
    
    #return all the values
    list(
      med_freqs = median(freqs),
      sd_freqs = sd(freqs),
      med_stopdurs = median(stopdurs),
      sd_stopdurs = sd(stopdurs),
      med_movedurs = median(movedurs),
      sd_movedurs = sd(movedurs)
    ) %>% return()
  } else {
    
    list(
      med_freqs = NA,
      sd_freqs = NA,
      med_stopdurs = NA,
      sd_stopdurs = NA,
      med_movedurs = NA,
      sd_movedurs = NA
    ) %>% return()
  }}
  
}

stat_generator <- function(traj){
  cutoffspeed = 0.1
  filename <- attr(traj, "filename")
  
  #create a list of speeds over time
  ############ Hey! Tracker speeds and TrajR-calculated speeds differ!#######
  traj %>% 
    TrajSmoothSG(p = 1, n = 11) %>%
    TrajDerivatives() %>% pluck("speed") -> smoothspeed
  traj %>% TrajDerivatives() %>% pluck("speed") -> rawspeed
  
  #get the speed ONLY during movement
  movespeed <- rawspeed %>% magrittr::extract(. > cutoffspeed)
  
  #get the stop-start interval variables
  int_vars <- intervaler(traj, display = FALSE)
  
  #Measures of straightness
  straightness <- TrajStraightness(traj)
  resampled <- TrajRediscretize(traj, .001)  #should I pick a different value?
  Emax <- TrajEmax(resampled)
  
  #attempts to get a spectral value...
  
  #detspeed <- pracma::detrend(rawspeed)
  #detspeed <- ts(detspeed, frequency = 31536000, start = 0)
  #I'm only getting the index? not useful really
  #spect <- GeneCycle::robust.spectrum(detspeed) 
  dfreqs <- dominant.freqs(rawspeed, 3)
  #acf(detspeed, main = filename)
  nobud <- nobu(rawspeed)
  title(filename)
  
  
  list(
    int_vars,
    med_movespeed = median(movespeed),
    sd_movespeed = sd(movespeed),
    straightness = straightness,
    Emax = Emax,
    dfreq1 = dfreqs[1],
    dfreq2 = dfreqs[2],
    dfreq3 = dfreqs[3],
    nobu = nobud,
    avg.speed = mean(rawspeed)
  ) %>% flatten() %>% return()
  
}

# Custom PCA plotting function
customPcaPlot <- function(x, xlabs, xcols, choices = 1L:2L, ycol = "#ff2222aa", ...) {
  # Draw points
  pts <- t(t(x$x[, choices]))
  plot(pts, type = "p", 
       xlim = extendrange(pts[, 1L]), ylim = extendrange(pts[, 2L]), 
       asp = 1,
       xlab = "PC1", ylab = "PC2", pch = 16, col = xcols, ...)
  text(pts, labels = xlabs, pos = 1, ...)
  
  # Draw arrows
  axs <- t(t(x$rotation[, choices])) * 3.5
  text(axs, labels = dimnames(axs)[[1L]], col = ycol, ...)
  arrows(0, 0, axs[, 1L] * .8, axs[, 2L] * .8, length = .1, col = ycol)
}

pca_maker <- function(stats, cutdownsource = F, keep.nas = FALSE, excluded.species = NULL, ...){
  col.list <- quos(...)
  
  #remove excluded species
  `%not_in%` <- purrr::negate(`%in%`)
  if (!is.null(excluded.species)) {stats <- stats %>% filter(species %not_in% excluded.species)}
  
  #remove or keep NAs
  if (keep.nas == FALSE) stats <- na.omit(stats);
  
  #use this for plotting
  all.pca.stats <- stats
  
  #remove columns that shouldn't be included
  excluded.cols <- c("individual", "species", "run")
  pca.stats <- stats %>% select(-one_of(excluded.cols))
  
  
  #select (or not) subset of columns to use
  if (length(col.list) >= 1) pca.stats <- pca.stats %>% select(!!!col.list);
  
  # Perform the PCA
  PCA <- prcomp(pca.stats, scale. = TRUE, center=TRUE)
  # Plot it using custom plotting function. Could just call biplot instead
  customPcaPlot(PCA, stats$individual, stats$species, cex = .8)
  #legend("bottomleft", c("Spider", "Mimic", "Ant"), pch = 16, 
  #       col = c('red', 'blue', 'black'), inset = c(0.01, .02))
  
  if (cutdownsource==F)  return(PCA) else return(all.pca.stats)
  
}


}

#___Begin################

path = "D:/Weaver ants/Dropbox (ASU)/Ant mimic strobing/Finished/"
files <- list.files(
  path= path,
  pattern = "[_][ABCDEF][.]csv", recursive = T
                    ) %>% 
  paste0("D:/Weaver ants/Dropbox (ASU)/Ant mimic strobing/Finished/",.)

trajs <- NULL
for (i in files) {
  print(i)
  
  #create a temp trajectory from the filename
  i %>% cvs_reader_func() %>%
    TrajFromCoords(xCol = "x", yCol = "y", fps = 240, timeUnits = "s") -> traj
  
  #give them a filename attribute (for plotting labels)
  
  attr(traj, "full.filename") <- i
  attr(traj, "filename") <- i %>% 
    strapplyc(paste0(
                    str_sub(path,-5,-2),   #this is the "shed" from "Finished"
                    "[/]([[:print:]]+[_][ABCDEF])[.]csv"
      ), simplify = T)
  
  trajs[[attr(traj, "filename")]] <- traj
  
}



#' # Intro for Lochlan
#' This is just a document to kinda display some of our data and
#' preliminary results. It'll give you an idea of what we have,
#' what we need to tinker with, and what we still should do.
#' 
#' Also, I realize that much of this might need my explanation, or we should
#' go over it together. Writing EVERYTHING out would just take too long.
#' This is like a visual aid, really.
#'
#' # Runs that definitely need some polishing up:
#' * M.bicolor4: clusters weird in PCA
#' * Omajor1_B: not really stop-start
#' * Ohaddoni5_B: a biiig pause
#' * Ohaddoni4_B: last half is a run
#' * Ohaddoni2_B: literally one burst
#' * Ohaddoni2_A: Hmmmm, kinda irregular but kinda good
#' * Mluctuosa6_A: big pause
#' * Mluctuosa5_A: needs smoothing/increase in threshold
#' * Mluctuosa4_D: sucks
#'   + 4_C: big pauses
#'   + 4_B: big pause
#'   + 4_A: uggggh
#' * Mluctuosa1_A & B: a bit of noise
#' * Mbicolor4_D: hmm. pause and needs smoothing/increase in threshold
#' * Mbicolor1_A: HUGE pause
#' 


big.db <-
  tibble(name = names(trajs),
         traj = trajs
        ) %>% 
  mutate(species = name %>% strapplyc("(.+)[0-9]_[A-Z]$", s=T) %>% as.factor(),
         individual = name %>% strapplyc("(.+[0-9])_[A-Z]$", s=T) %>% as.factor()
         )

big.db <- 
  TrajsMergeStats(big.db$traj, stat_generator) %>%
  bind_cols(big.db, .)

big.db$intervals <-
  lapply(big.db$traj, function(x) intervaler(x, getdbs = T) %>% .$intervals)
big.db$REVERSEintervals <-
  lapply(big.db$traj, function(x) intervaler(x, getdbs = T) %>% .$REVERSEintervals)
big.db$freqs <-
  lapply(big.db$traj, function(x) intervaler(x, getdbs = T) %>% .$freqs)


per.ind <-
big.db %>%
  group_by(individual, species) %>%
  summarize(
    med_stopdurs = bind_rows(intervals) %>%
      filter(duration<.25) %>% 
      .$duration %>% median(),
    sd_stopdurs = bind_rows(intervals) %>%
      filter(duration<.25) %>%
      .$duration %>% sd(),
    mad_stopdurs = bind_rows(intervals) %>%
      filter(duration<.25) %>% 
      .$duration %>% mad(),
    med_movedurs = bind_rows(REVERSEintervals) %>%
      .$duration %>% median(),
    sd_movedurs = bind_rows(REVERSEintervals) %>%
      .$duration %>% sd(),
    mad_movedurs = bind_rows(REVERSEintervals) %>%
      .$duration %>% mad(),
    med_freqs = unlist(freqs) %>%
      median(),
    sd_freqs = unlist(freqs) %>%
      sd(),
      mad_freqs = unlist(freqs) %>%
      mad(),
    ) %>% ungroup()


pca.stated <- per.ind %>% filter(species!="Goldcampo", species!="Smallgoldpoly")
prcomped <- 
  per.ind %>% filter(species!="Goldcampo", species!="Smallgoldpoly") %>%
  select(-individual,-species) %>% 
  prcomp(scale. = TRUE, center=TRUE) 

autoplot(prcomped,
         data=pca.stated %>% set_rownames(pca.stated$individual),
         colour="species", label=T,
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3, size=5, x=1,y=2)

fviz_pca_ind(prcomped,
             label="none",
             col.ind = pca.stated$species, # color by groups
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "convex",
             ellipse.level=0.99,
             legend.title = "Groups",
             repel = TRUE
)


cbind(spec=pca.stated$species,val=prcomped$x[,1]) %>% boxplot(val ~ spec, data=.)



tapply(big.db$intervals, big.db$individual, bind_rows)

#+ stats, out.width=c('50%', '50%'), fig.show='hold'
INTERVAL.DB <- NULL
stats <- TrajsMergeStats(trajs, stat_generator)
rownames(stats) <- names(trajs) %>% unlist()
stats$species <- row.names(stats) %>%
  strapplyc("(.+)[0-9]_[A-Z]$") %>% unlist() %>% as.factor()
stats$individual <- row.names(stats) %>%
  strapplyc("(.+[0-9])_[A-Z]$") %>% unlist() %>% as.factor()
stats$run <- row.names(stats)


#' # Here's a PCA of the runs
#' 
#' This sort of thing is the real clincher; it'll probably be the
#' most important figure in our paper. It will BE the paper, haha.

#+ pca-dont-print1, include=FALSE
pca_maker(stats[c(1,2,3,4,5,6,7,8,15:17)], cutdownsource = F, FALSE, c("Goldcampo","Smallgoldpoly")) -> prcomped
pca_maker(stats[c(1,2,3,4,5,6,7,8,15:17)] %>% mutate(run=rownames(.)), cutdownsource = T, FALSE, c("Goldcampo","Smallgoldpoly")) -> pca.stated

#+ pcas1, fig.cap="This is made only using variables about the stopping and starting"
autoplot(prcomped,
         data=pca.stated %>% set_rownames(pca.stated$run),
         colour="species", label=T,
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3, size=5, x=1,y=2)


#+ pca-dont-print2, include=FALSE
pca_maker(stats[c(1,3,5,7,15:17)], cutdownsource = F, FALSE, c("Goldcampo","Smallgoldpoly")) -> prcomped
pca_maker(stats[c(1,3,5,7,15:17)] %>% mutate(run=rownames(.)), cutdownsource = T, FALSE, c("Goldcampo","Smallgoldpoly")) -> pca.stated

#+ pcas2, fig.cap="This was made excluding the standard deviation variables from above."
autoplot(prcomped,
         data=pca.stated %>% set_rownames(pca.stated$run),
         colour="species", label=T,
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3, size=5, x=1,y=2)



#+ pca-dont-print3, include=FALSE
pca_maker(stats[c(1:9,11:17)], cutdownsource = F, FALSE, c("Goldcampo","Smallgoldpoly")) -> prcomped
pca_maker(stats[c(1:9,11:17)] %>% mutate(run=rownames(.)), cutdownsource = T, FALSE, c("Goldcampo","Smallgoldpoly")) -> pca.stated

#+ pcas3, fig.cap="This last one includes EVERYTHING I have, like straightness of the trajectory and the estimated 'dominant frequencies' of their speed over time, etc."
autoplot(prcomped,
         data=pca.stated %>% set_rownames(pca.stated$run),
         colour="species", label=T,
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3, size=5, x=1,y=2)



#+ plot
ggplot(stats %>% filter(species!="Goldcampo") %>% filter(species!="Smallgoldpoly"),
       aes(y=med_movespeed, x=med_movedurs, color=species)) +
  geom_point(size=6)

#' # Which variables matter?
#' Okay, now let's look at which variables they actually differ in.
#' 
#' * "med_" means "median value of _" and "durs" means duration. "freqs" refers
#' to the combined stop-start period.
#' * "dfreq[123]" refers to the "dominant frequencies" of the time series. It's
#' some confusing spectral analysis stuff, but it should theoretically
#' describe how the periods of stopping and starting act. I've included
#' the first three dominant frequencies.
#' * "nobu" (from my labmate Nobuaki) represents the time length of highest
#' rhythmicity, so it's kinda similar.

#+ var-compare, out.width=c('50%', '50%'), fig.show='hold', warning=FALSE
stats2 <- stats %>%
  filter(species!="Goldcampo") %>% filter(species!="Smallgoldpoly") %>% 
  mutate(species=factor(species))
for (i in colnames(stats2)) {
    try(plot(stats2[[i]] ~ stats2[['species']], 
             main=paste0("Comparing '", i, "' across strobing species"),
             ylab=i, xlab="species" ))
}

stats %>% 
  filter(species!="Goldcampo") %>% filter(species!="Smallgoldpoly") %>% 
  plot(data=., med_movedurs ~ species)



# 
# 
INTERVAL.DB$species <-
  INTERVAL.DB$filename %>% strapplyc("(.+)[0-9]_[A-Z]$") %>% unlist() %>% as.factor()

stops <- subset(INTERVAL.DB, source == "intervals" )
moves <- subset(INTERVAL.DB, source == "REVERSEintervals" )
stops %>% group_by(species) %>% mutate(duration =  remove_outliers(duration)) %>% na.omit() %>% ungroup() -> stopsgood


# INTERVAL.DB %>% 
#   #filter(duration>0.03) %>% 
#   group_by(species) %>% mutate(duration =  remove_outliers(duration)) %>% na.omit() %>% ungroup() %>% filter(species!="Goldcampo") %>% filter(species!="Smallgoldpoly")  %>% mutate(species = species %>% as.numeric() %>% as.factor())->hey
# 
# 
# 
# ggplot(subset(hey, source=="REVERSEintervals"), aes(x=species, y=duration, fill=filename)) + 
#   geom_boxplot() 



#trajs[[20]] %>% TrajDerivatives() %>% {cbind(.$speed,c(0,.$acceleration))} %>% as.data.frame() %>%  segclust(seg.var=c("V1","V2"), lmin=5, Kmax=50, scale.variable=F, ncluster = c(2)) %>% plot()

#' # Other notes
#' Also, I think we found (as other people have) that when M.luctuosa pauses,
#' they wave their ("extra") front legs around like antennae to mimic ants better.
#' However, I *think* M.bicolor *doesn't* do that: they keep their front legs
#' still. Excitingly, when the Opisthopsis stop, unlike other ants, they
#' don't move their antennae at all. Mimicry!
#' 

#___Session info###############

#' # Session info
#' 
#+ session.info
sessionInfo()

