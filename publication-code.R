#' ---
#' title: "Spider Ant Mimics"
#' author: "Andrew Burchill"
#' date: "2019-06-19"
#' output:
#'    html_document:
#'       theme: lumen
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
library(multcomp)
library(sandwich)





dendro_analysis <- function(per.ind) {
  
  mimicry.relatedness <- per.ind %>% 
    mutate(
      group = case_when(
        .$species == "Goldcampo" ~ "1",
        .$species == "Smallgoldpoly" ~ "1",
        .$species == "O.haddoni" ~ "2",
        .$species == "O.major" ~ "2",
        .$species == "M.luctuosa" ~ "1",
        .$species == "M.bicolor" ~ "2"
      )) %>% .$group
  
  
  
  clustd <-per.ind %>% as.data.frame() %>% 
    set_rownames(per.ind %>%
                   mutate(species = species %>% 
                            recode_factor("M.bicolor" = "M. bicolor",
                                          "O.haddoni" = "O. haddoni",
                                          "O.major" = "O. major",
                                          "M.luctuosa" = "M. luctuosa",
                                          "Smallgoldpoly" = "P. aurea",
                                          "Goldcampo" = "C. bigenus"
                            ),
                          individual = paste0(species, str_sub(individual, -1,-1))
                   ) %>% .$individual
    ) %>% .[-1:-2] %>% scale %>% dist() %>% hclust(method="ward")
  
  
  fviz_dend(rev(clustd), palette = "jco", rect = T,  rect_lwd=.5,
            show_labels = T, labels_track_height = 10, k=2, lwd=1, 
            lower_rect = -12, main= "B) Cluster Dendogram",
            label_cols = mimicry.relatedness[clustd$order] %>% map(.,~ifelse(.=="1","#0073C2FF","#EFC000FF")) %>% unlist() %>% rev())
  
  #return(clustd)
  
  
}
  

########

MATLAB_results_path <- "sync/ant_time_series.txt"


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
      ". MOVING intervals: ", length(REVERSEintervals$duration)))
  
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
      freqs <- 
        intervals$startTime %>% 
        diff() %>% 
        magrittr::extract(. > greatthan)
      
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
      
      
      #return all the values
      list(
        med_freqs = median(freqs),
        sd_freqs = sd(freqs),
        med_stopdurs = median(stopdurs),
        sd_stopdurs = sd(stopdurs),
        med_movedurs = median(movedurs),
        sd_movedurs = sd(movedurs),
        prop_time_stopped = sum(stopdurs)/max(traj$time)+0
      ) %>% return()
    } else {
      
      list(
        med_freqs = NA,
        sd_freqs = NA,
        med_stopdurs = NA,
        sd_stopdurs = NA,
        med_movedurs = NA,
        sd_movedurs = NA,
        prop_time_stopped = 0
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
  
  title(filename)
  
  
  list(
    int_vars,
    med_movespeed = median(movespeed),
    sd_movespeed = sd(movespeed),
    straightness = straightness,
    avg.speed = mean(rawspeed),
    sd.speed = sd(rawspeed),
    total_time = traj$time %>% tail(1)
  ) %>% flatten() %>% return()
  
}
PCA_it <- function(pca.prepped, title="", ...){
  
  #exclude certain columns before PCA
  col.list <- quos(...)
  prcomped <-
    pca.prepped %>% 
    select(-individual,-species) %>% 
    prcomp(scale. = TRUE, center=TRUE) 
  
  
  file<-paste0(title,"PCA.png", collapse = "-")
  ragg::agg_png(file, width= 1900, height = 1080, res=300)
  
  fviz_pca_ind(prcomped,
               label="none",
               title=title,
               col.ind = pca.prepped$species, # color by groups
               addEllipses = TRUE, # Concentration ellipses
               ellipse.type = "convex",
               ellipse.level=0.99,
               legend.title = "Groups",
               repel = TRUE,
               axes.linetype=NA,
               palette = palette(c( "#76480d",
                                    "#dc3c07",
                                    "#f8ba7c",
                                    "#256676",
                                    "#4bd6fd",
                                    "#6a9fee"
               )),
  ) %>% 
    plot()
  dev.off()
  invisible(readline(prompt="Press [enter] to continue"))
  
  
  PC1_results <-
    bind_cols(Species = pca.prepped$species, PC1 = prcomped$x[, 1])
  
  
  file<-paste0(title,"boxplot.png", collapse = "-")
  ragg::agg_png(file, width= 1900, height = 1080, res=300)
  
  
  print(
    ggplot(PC1_results, aes(x = Species, y = PC1)) +
      geom_boxplot(lwd = 1.5) + scale_fill_manual(values = palette(
        c(
          "#76480d",
          "#f8ba7c",
          "#dc3c07",
          "#4bd6fd",
          "#6a9fee",
          "#256676"
        )
      )) + theme_classic(base_size = 15) + theme(legend.position = "none")
  )
  dev.off()
  # boxplot(PC1 ~ Species, data=PC1_results, main="PC1 values across species")
  
  
  #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2847912/pdf/pone.0009788.pdf
  #this method is good
  model <- aov(PC1 ~ Species, data=PC1_results)
  
  mc = glht(model,
            mcp(Species = "Tukey"),vcov=vcovHC)
  
  mcs = summary(mc, test=adjusted("single-step"))
  cld(mcs,
      level=0.05,
      decreasing=TRUE) %>% print()
  return(PC1_results)
}
start_stop_PCA <- function(temp.db) {
  
  
  temp.db$intervals <-
    lapply(temp.db$traj, function(x) intervaler(x, getdbs = T) %>% .$intervals)
  temp.db$REVERSEintervals <-
    lapply(temp.db$traj, function(x) intervaler(x, getdbs = T) %>% .$REVERSEintervals)
  temp.db$freqs <-
    lapply(temp.db$traj, function(x) intervaler(x, getdbs = T) %>% .$freqs)
  
  per.individual <-
    temp.db %>% 
    group_by(individual, species) %>%
    summarize(
      med_stopdurs = bind_rows(intervals) %>%
        filter(duration<.25) %>%  #pauses longer than this aren't part of strobing
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
      # med_freqs = unlist(freqs) %>%
      #   median(),
      sd_freqs = unlist(freqs) %>%
        sd(),
      mad_freqs = unlist(freqs) %>%
        mad()
    ) %>% ungroup() %>% select(-starts_with("sd"))
  
  pca.prep <- per.individual %>% na.omit()
  PCA_it(pca.prep,"Movement and Pause Durations",-individual,-species)
  
}
import_MATLAB <- function(temp.db, path) {
  
  #get the rhythmicity stuff
  ant_time_series <- read_csv(path) 
  ant_time_series$Colony %<>% substr(.,1,nchar(.)-4)
  ant_time_series%<>% rename(name=Colony)
  return(full_join(ant_time_series, temp.db))
  
}
MATLAB_PCA <- function(temp.db) {
  
  per.individual <-
    temp.db %>%
    group_by(individual, species) %>%
    summarize(
      sync = log(weighted.mean(Synchrony, total_time)),
      rhy = log(weighted.mean(Rhythmicity, total_time)),
      per = (weighted.mean(Period, total_time)),
      str = car::logit(weighted.mean(straightness, total_time)),
      stopped = weighted.mean(prop_time_stopped, total_time),
      speed = log(weighted.mean(avg.speed, total_time))
      
    ) %>% ungroup() %>%
    mutate(per=bestNormalize::boxcox(per)$x.t)
  
  #make sure they're all normal
  per.individual %>%
    summarize_at(vars(-individual,-species), ~shapiro.test(.)$p.value)
  
  PCA_it(per.individual, "Full Wavelet Analysis", -individual, -species) -> test
  #agricolae::kruskal(test$PC1, test$species,group=TRUE,p.adj="bonferroni")$groups
  
  
  return(per.individual)
}
temp.db <- NULL

raw_data <- readRDS("raw_data.RDS") 

big.db <- 
  TrajsMergeStats(raw_data$traj, stat_generator) %>%
  bind_cols(raw_data, .)

big.db %>% start_stop_PCA()

big.db %>% import_MATLAB(MATLAB_results_path) %>% MATLAB_PCA()



