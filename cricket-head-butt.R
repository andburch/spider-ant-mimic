#' ---
#' title: "Weaver Ants"
#' author: "Andrew Burchill"
#' date: "2019-06-07"
#' output:
#'    html_document:
#'       theme: lumen
#'       includes:
#'         after_body: footer.html
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
#___Load Packages################

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4806549/
#   https://rviews.rstudio.com/2017/09/25/survival-analysis-with-r/


list.of.packages <- c("magrittr", 
                      "trajr", 
                      "gganimate", 
                      "data.table",
                      "gsubfn",
                      "grid",
                      "zoo",
                      "circular",
                      "survminer",
                      "survival",
                      "ggfortify",
                      "coxphw",
                      "knitr",
                      "kableExtra",
                      "tidyverse" 
                      )
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
for (i in list.of.packages) library(i, character.only = T)
select <- dplyr::select    #because FUCK the MASS packages

# 
# Todo.
# 1) Figure out a cut-off for slow-starting horizontal runs
#   - either min distance based off the "furthest away" starting point, or
#   - cutoff based on the max starting speed
# 1.5) Also, don't forget about the ending cutoff too!
# 2) Characterize the speed-over-time graphs for the functions?
# 3) Compare m1s and m2s to see if there's consistency


#___Functions#################

#gets data from the mdf files
extractor <- function(file.name, type=c("cricket","sequence")) {
  
  if(type == "cricket") {
    
    #reads the file and extracts the right parts
    f <- file(paste0(getwd(),"/",file.name))
    x <- readLines(f)
    close(f)
    
    #this gets the px/cm conversion for each video.
    #they're annoyingly going to be stored in the display values for opacity and point size
    # px.cm.ratio <- x[2] %>%
    #   strapplyc("^Displaying\ \\w+\ \\w+\ \\w+\ \\w+\ \\w+\ \\w+\ \\w+\ (\\w+)\ (\\w+)") %>% 
    #   unlist() %>% paste(sep=".", collapse = ".") %>% as.numeric()
    # 
    #x[2] %>%
    #strapplyc("^Displaying\ \\w+\ \\w+\ \\w+\ \\w+\ \\w+\ \\w+\ \\w+\ \\w+\ \\w+\ \\w+\ \\w+\ \\w+\ (\\w+)\ \\w+\ (\\w+)\ \\w+")
    
    #jesus christ what a bad way to select the right ratio
    px.ratio <- px.cm.ratios[which(px.cm.ratios$file.name == file.name),]$px.cm.ratio
    if (length(px.ratio)==0) px.ratio <- 1
    
    track.lines <- union(grep("^Track", x),grep("^End", x))
    
    head <- data.table::fread(text=x, skip = track.lines[1], nrows=track.lines[2]-track.lines[1]-1) %>%
      as_tibble() %>% select(ID = V2, x = V3, y = V4, frame = V6) %>% mutate(x = x/px.ratio, y = y/px.ratio)
    
    butt <- data.table::fread(text=x, skip = track.lines[2], nrows=track.lines[3]-track.lines[2]-1) %>%
      as_tibble() %>% select(ID = V2, x = V3, y = V4, frame = V6) %>% mutate(x = x/px.ratio, y = y/px.ratio)
    
    
    #This is for videos that shift....
    if (length(track.lines) > 3) {
      #this should be the THIRD track, and it should have only two points
      jump.ref <- data.table::fread(text=x, skip = track.lines[3], nrows=track.lines[4]-track.lines[3]-1) %>%
        as_tibble() %>% select(ID = V2, x = V3, y = V4, frame = V6)
      try(if(nrow(jump.ref) != 2) stop("The third track doesn't have only 2 rows"))
      
      shift.x <- jump.ref$x[2] - jump.ref$x[1]
      shift.y <- jump.ref$y[2] - jump.ref$y[1]
      shift.frame <- jump.ref$frame[1]
      
      #use magrittr's crazy two-way pipe
      head[head$frame > shift.frame,] %<>% mutate(x = x - shift.x,
                                                  y = y - shift.y)
      butt[butt$frame > shift.frame,] %<>% mutate(x = x - shift.x,
                                                  y = y - shift.y)
    }
    
    
    #turns em into a trajectory. why? idk
    th <- head %>% TrajFromCoords(xCol = "x", yCol = "y", fps = 240, timeUnits = "s")
    tb <- butt %>% TrajFromCoords(xCol = "x", yCol = "y", fps = 240, timeUnits = "s")
    
    #h is for head, b is for butt, m is for midpoint
    combined <- cbind(h = th, b = tb) %>% mutate(length = d(h.x,h.y,b.x,b.y), m.x = (h.x+b.x)/2, m.y = (h.y+b.y)/2 )
    combined <- rowid_to_column(combined)
    
    #get the reference point info
    ref <- data.table::fread(text=x, skip = grep("^Ref", x)-1, nrows=1) %>% 
      mutate(V2=V2/px.ratio, V3 = V3/px.ratio) %>% #converts based on image size
      as_tibble() 
    
    #calculate the angle between ref, head, and butt
    combined$angle.from.ref <- 
      combined[c("rowid", "h.x","h.y","b.x","b.y")] %>%
      apply(1, find_angles, ref=ref)
    combined$angle.from.ref <- 
      suppressWarnings(circular::as.circular(combined$angle.from.ref, units="degrees"))
    
    
    #get distances between each step and the reference + the start point
    combined <- combined %>% mutate(dist.from.ref = d(m.x,m.y, ref[[2]],ref[[3]]), dist.from.start = d(m.x,m.y, combined$m.x[1],combined$m.y[1]))
    
    #pass along the name (mostly for plotting)
    attr(combined, "file.name") <- file.name
    #also pass along the reference point
    attr(combined, "ref.x") <- ref[[2]]
    attr(combined, "ref.y") <- ref[[3]]
    
    
    return(combined)
    
  } else if(type == "sequence") {
    #reads the file and extracts the right parts
    f <- file(paste0(getwd(),"/",file.name))
    x <- readLines(f)
    close(f)
    
    # #jesus christ what a bad way to select the right ratio
    # px.ratio <- px.cm.ratios[which(px.cm.ratios$file.name == file.name),]$px.cm.ratio
    # if (length(px.ratio)==0) px.ratio <- 1
    
    track.lines <- union(grep("^Track", x),grep("^End", x))
    
    tracks<-NULL
    #go through all the tracks, but don't get the last one (it's "End file")
    for (i in seq_along(track.lines) %>% head(-1)) {
      tracks[[as.character(i)]] <- data.table::fread(text=x, skip = track.lines[i], nrows=track.lines[i+1]-track.lines[i]-1) %>%
        as_tibble() %>% select(ID = V2, x = V3, y = V4, frame = V6) 
    }
    
    #get the start and end frames of the tracks
    times <- tibble(min= numeric(), max=numeric())
    for (i in seq_along(tracks)) {
      times[i,]$max <- tracks[[i]]$frame %>% max()
      times[i,]$min <- tracks[[i]]$frame %>% min()
    }
    #create duration attached, noting if they're holding on at the end
    ENDING <- file.name %>% strapplyc(., "[-]([0-9]+).mdf", simplify = TRUE) %>% as.numeric()
    times %<>% mutate(length = max-min, right.censored = ifelse(max >= ENDING, 0, 1))
    return(times)
  } 
  
}


#calculate the angle between the reference point, the head, and the butt
find_angles <- function(x,ref) {
  
  temp.ref <- ref
  #transform points so head is at 0,0
  #transform the x's
  temp.ref[2] <- temp.ref[[2]] - x[2]
  x[4] <- x[4] - x[2]
  
  #transform the y's
  temp.ref[3] <- temp.ref[[3]] - x[3]
  x[5] <- x[5] - x[3]
  
  #calculates angle from x axis and butt and then subtracts the angle from x axis to REF
  #this gets the angle (in the clockwise direction) from REF to head to butt.
  # rad <- atan2(x[4], x[5]) - atan2(temp.ref[[2]], temp.ref[[3]])
  
  #This gets the angle between the butt (the head) and the vertical
  rad <- atan2(x[4], x[5])
  #turn into degrees
  degree <- rad*180/pi
  #make everything positive values 
  #(also add 90deg cuz our plot are horizontal)
  degree <- (degree+360 +90) %% 360
  
  #if the ref is on the right side of the ~4000px video
  #then you'll need to invert the values 180deg
  if (ref[[2]] > 2000) {degree <- (degree+180) %% 360}
  
  return(degree)
  
  
}

#gets distance
d <- function(x1,y1,x2,y2) {
  sqrt( (x1-x2)^2 + (y1-y2)^2 ) %>% return()
}

#plots the crickets journey
plot_cricket <- function(combined) {
  
  ref <- c(attr(combined, "ref.x"), attr(combined, "ref.y"))
  
  combined[c("rowid", "h.x","h.y","b.x","b.y")] %>%
    gather(variable, value, -rowid) %>%
    mutate(spot = substr(variable,1,1), type=substr(variable,3,3) ) %>%
    select(-variable) %>% 
    spread(type ,value = value) -> df
  
  df[dim(df)[1]+1,] <- c(max(df$rowid)+1, "R", ref[1], ref[2])
  df$x %>% as.numeric() -> df$x
  df$y %>% as.numeric() -> df$y
  plot.cricket <- 
    ggplot(df, aes(x=x, y=y)) +
    geom_line(aes(group=rowid), color = "tan", cex=2, alpha =.2) +
    geom_line(aes(group= spot), alpha=.1) + 
    labs(captionn=attr(combined,"file.name")) +
    geom_point(aes(colour = factor(spot)), cex=3) + #transition_reveal(rowid) + 
    xlim(0,4096) + ylim(0,2160) #+    geom_text(aes(label=rowid),hjust=0, vjust=0)
  
  return(plot.cricket)
 # gganimate(beauty, "OUTPUT.gif")
}

#only works with rowid... 
#plots angle change over timesteps
plot_angle <- function(combined, ...) {
  
  sequence.var <- quo(...)
  
  #used for polar angle plotting
  #gets how many times it crosses 0/360 line
  which_revolution <- function(x) {
    yo<-NULL
    yo[1] <-0
    angle.diffs <- x %>% diff() %>% add(180) %>% mod(360) %>% subtract(180)
    angle.diffs[length(angle.diffs)+1] <- 0
    
    for (i in seq_along(angle.diffs)) {
      
      if ((angle.diffs[i] + x[i] > 360) | (angle.diffs[i] + x[i] < 0)) {
        yo[i+1] <- yo[i] +1
      } else {
        yo[i+1] <- yo[i]
      }
      
    }
    return(head(yo,-1))
  }
  
  #labels which points were made up
  fake_point_labeller <- function(df) {
    df$fake.points <- "no"
    df$fake.points[
      df$angle.from.ref %>% sapply(function(x) identical(x,0)) %>% which()
      ] <- "yes"
    df$fake.points[
      df$angle.from.ref %>% sapply(function(x) identical(x,360)) %>% which()
      ] <- "yes"
    return(df)
  }
  
  threshold <- 360
  df <- combined %>% select(angle.from.ref, !!sequence.var)
  df$revolution <- df$angle.from.ref %>% which_revolution()
  
  #complicated, but makes fake points, one on the 0+ and one on the 360- side of the line each time the path crosses it
  if(n_distinct(df$revolution) > 1){        
    split.df <- split(df, df$revolution)
    df <- split.df[[1]]
    for(i in seq_along(split.df)[-1]){
      interp.df <- rbind(df[df[[2]] == max(df[[2]]), ],
                         split.df[[i]][split.df[[i]][[2]] == min(split.df[[i]][[2]]), ])
      interp.df[[2]] <- interp.df[[2]][[1]] + 
        (threshold - interp.df[[1]][1]) / 
        (threshold - interp.df[[1]][1] + interp.df[[1]][2]) *
        diff(interp.df[[2]])
      
      if (interp.df[[1]][1] < interp.df[[1]][2]){
        interp.df[[1]] <- c( 0,threshold) 
      } else {
        interp.df[[1]] <- c(threshold,  0) 
      }
      
      df <- rbind(df, interp.df, split.df[[i]])
    }
  }
  
  df <- fake_point_labeller(df)
  
  plot.angs <- ggplot(df, aes(x=angle.from.ref, y=!!sequence.var, group=revolution)) +
    coord_polar(theta="x", direction=1) +
    geom_path(aes(color = factor(revolution)), size = 1) + 
    scale_x_continuous(limits = c(0,360), breaks = seq(0, 360, 45)) +
    labs(caption = attr(combined,"file.name")) +
    geom_point(data=subset(df, fake.points == "no")) 
  return(plot.angs)
  
}


stat_generator <- function(traj){
  #cutoffspeed <- 0.001
  resamp.rate <- .3
  filename <- attr(traj, "file.name")
  
  total.distance <- traj %>% tail(n=1) %>% pluck("dist.from.start")
  total.time <- traj %>% tail(n=1) %>% pluck("displacementTime")
  
  
  total.avg.speed <- traj %>%
    tail(n = 1) %>% 
    transmute(speed = dist.from.start / displacementTime) %>%
    as.numeric() %>%
    multiply_by(24) #frames per second
  
  ##This works for data without reference points
  #get the signage right for forward/backward speed
  total.x.disp <- sum(traj$displacement) %>% Re()
  signed.x.disp <- (total.x.disp + Re(traj$displacement)) %>%
    abs() %>%
    subtract(abs(total.x.disp))
  
  
  ####assuming that the biggest speeds DON'T happen during the switches between movies
  ####e.g. they only happen between slices that are 100 frames apart!
  #top.forward.speed <- max(signed.x.disp) / (100/24)    #100 frames at 24fps
  #top.backward.speed <- min(signed.x.disp) / (100/24)   #100 frames at 24fps
  
  #some of the small trajectories don't have 12 data points!
  if (dim(traj)[1] < 14) w = 3 else w = 12
  ###Fuck it, assume all slices are 100 frames, okay? we gotta HURRY!
  #SLIDING WINDOW OF 1/20TH THE LENGTH OF TRAJECTORY
  top.forward.speed <- 
    signed.x.disp %>% 
    rollapply(w = w, mean, partial = F) %>% 
    max() %>% divide_by(100/24) #100 frames at 24fps
  
  #SLIDING WINDOW OF 5 FRAMES
  #WRITING OVER THIS
  top.backward.speed <- 
    signed.x.disp %>% 
    rollapply(w = w, mean, partial = F) %>% 
    min() %>% divide_by(100/24) #100 frames at 24fps
  
  #actual proportion of going backward
  top.backward.speed <- 
    signed.x.disp %>% 
    rollapply(w = w, mean, partial = F) %>%
    is_less_than(0) %>% 
    sum(na.rm = TRUE) %>% #counts how many frames are going backwards
    divide_by(length(signed.x.disp)/100)
  
  
  signed.x.disp %>% 
    rollapply(w = w, mean, partial = F) %>%
    plot(main=filename)
  
  avg.speed <- traj %>% TrajDerivatives() %>% pluck('speed') %>% mean() %>% multiply_by(24)
  
  #Measures of straightness
  straightness <- TrajStraightness(traj)
  sinuosity <- TrajSinuosity2(traj)
  Emax <- TrajEmax(traj)
  resampled <- TrajRediscretize(traj, resamp.rate)  #should I pick a different value?
  rEmax <- TrajEmax(resampled)
  rsinuosity <- TrajSinuosity2(resampled)
  rstraightness <- TrajStraightness(resampled)
  
  #plot(traj)
  #plot(resampled)
  
  #circ <- circuclar::circular(traj$angle.from.ref, units="degrees") 
  #Variance/sd/mean of cricket angle
  angle.var <- var.circular(traj$angle.from.ref)
  angle.sd <- sd.circular(traj$angle.from.ref)
  angle.mean <- mean.circular(traj$angle.from.ref)
  plot(traj$angle.from.ref, points.plot = TRUE); lines(density(traj$angle.from.ref, bw = 25, adjust = 2), shrink=1.5, offset=.6, col="red")
  
  
  list(
    total.avg.speed = total.avg.speed,
    top.forward.speed = top.forward.speed,
    top.backward.speed = top.backward.speed,
    straightness = straightness,
    sinuosity = sinuosity,
    #Emax = Emax,
    #rEmax = rEmax,
    #rsinuosity = rsinuosity,
    #rstraightness = rstraightness,
    #angle.var =  angle.var,
    #why use var if we have sd?
    angle.sd = angle.sd,
    angle.mean = angle.mean,
    avg.speed = avg.speed,
    total.time = total.time,
    total.distance = total.distance
  ) %>% flatten() %>% return()
  
}

get_metadata <- function(files, type = c("cricket","sequence","ant.number")) {
  
  if(type == "cricket") {
    angle = files %>% 
      strapplyc("^[[:print:]]+/([0-9]+)[msl][0-9]", simplify = T) %>% 
      ordered(., levels=c("0","45","90"))
    date = files %>% 
      strapplyc("^[[:print:]]+/([[:print:]]+)/[0-9]+[msl][0-9]", simplify = T) %>%
      substr(1,6) %>% as.factor()
    size = files %>% 
      strapplyc("^[[:print:]]+/[0-9]+([msl])[0-9]", simplify = T) %>% 
      ordered(., levels=c("s","m","l"))
    run.num = files %>% 
      strapplyc("^[[:print:]]+/[[:print:]]+/[0-9]+[msl]([0-9])", simplify = T) %>%
      ordered(., levels=c("0","1","2","3"))
    colony.id = date %>% tibble(date =.) %>% 
      transmute(colony = case_when(
        .$date %in% c("oct29", "oct30", "oct31") ~ "narrowparkend",
        .$date %in% c("oct23", "oct24", "oct25") ~ "treeline",
        .$date %in% c("oct20", "oct22", "nov03") ~ "firstcolony")) %>% 
      unname() %>% unlist() %>% as.factor()
    
  } else if(type == "sequence") {
    angle = files %>% 
      strapplyc("([0-9]+)[a-z]+[-][sm][0-9][-][0-9]+.mdf", simplify = T) %>% 
      ordered(., levels=c("0","45","90"))
    date = files %>% 
      strapplyc("^[[:print:]]+/([[:print:]]+)/", simplify = T) %>%
      substr(1,6) %>% as.factor()
    size = files %>% 
      strapplyc("([sml])[0-9][-][0-9]+.mdf", simplify = T) %>% 
      ordered(., levels=c("s","m","l"))
    run.num = files %>% 
      strapplyc("[sml]([0-9])[-][0-9]+.mdf", simplify = T) %>%
      ordered(., levels=c("0","1","2","3"))
    colony.id = date %>% tibble(date =.) %>% 
      transmute(colony = case_when(
        .$date %in% c("oct29", "oct30", "oct31") ~ "narrowparkend",
        .$date %in% c("oct23", "oct24", "oct25") ~ "treeline",
        .$date %in% c("oct20", "oct22", "nov03") ~ "firstcolony")) %>% 
      unname() %>% unlist() %>% as.factor()
#this is unique for the attachment data
    seq.type = files %>% 
      str_extract_all("(seq|timepoint)") %>% unlist() %>% as.factor()
    run.num = bind_cols(run.num=run.num, seq.type = seq.type)
  } else if(type=="ant.number") {
    
    angle = files %>% 
      strapplyc("([0-9]+)num[-]", simplify = T) %>% 
      ordered(., levels=c("0","45","90"))
    date = files %>% 
      strapplyc("^[[:print:]]+/([[:print:]]+)/", simplify = T) %>%
      substr(1,6) %>% as.factor()
    size = files %>% 
      strapplyc("[0-9]+num[-]([sml])", simplify = T) %>% 
      ordered(., levels=c("s","m","l"))
    run.num = files %>% 
      strapplyc("[0-9]+num[-][sml]([0-9])", simplify = T) %>%
      ordered(., levels=c("0","1","2","3"))
    colony.id = date %>% tibble(date =.) %>% 
      transmute(colony = case_when(
        .$date %in% c("oct29", "oct30", "oct31") ~ "narrowparkend",
        .$date %in% c("oct23", "oct24", "oct25") ~ "treeline",
        .$date %in% c("oct20", "oct22", "nov03") ~ "firstcolony")) %>% 
      unname() %>% unlist() %>% as.factor()
  }
  
  
  return(bind_cols(angle=angle,
                   date=date,
                   size=size,
                   run.num=run.num,
                   colony.id=colony.id)
  )
}


help_special_snowflakes <- function(all.of.em, type = c("cricket","sequence")) {
  
  #because of a hand occluding things, #5 in this run couldn't be followed
  #so it needs to be right-censored
  all.of.em[which(
    all.of.em$file.name %>% grepl("45seq-s1-1907",.) &
      all.of.em$min == 502)
    ,] %<>% mutate(right.censored = 0)
  
  return(all.of.em)
  
}

#tells you what data you need to collect
whats_needed <- function() {
  cat("\n===Runs with less than two trajectories===\n")
  big.db$file.name %>%
    get_metadata("cricket") %>%
    bind_cols(file.names=big.db$file.name, .) %>%
    group_by(date,angle, size) %>% 
    summarize(n=n()) %>% 
    subset(n<2) %>% print.data.frame()
  
  cat("\n===Runs without count data===\n") 
  big.db %>% filter(!map(ant.counts, length %>% unlist()) > 1) %>% group_by(date,angle, size) %>% 
    summarize(n=n()) %>% arrange(desc(n)) %>% print.data.frame()
  
  cat("\n===Runs without attachment data===\n") 
  big.db %>% filter(!map(indvd.seq, length %>% unlist()) > 1) %>% group_by(date,angle, size) %>% 
    summarize(n=n()) %>% arrange(desc(n)) %>% print.data.frame()
  
  cat("\n**********************************\n") 
  
  cat("\n===Which groups have COUNT data?===\n") 
  big.db %>% filter(map(ant.counts, length %>% unlist()) > 1) %>% group_by(angle, size) %>% 
    summarize(n=n()) %>% arrange(desc(n)) %>% print.data.frame()
  
  cat("\n===Which groups have ATTACHMENT data?===\n") 
  big.db %>% filter(map(indvd.seq, length %>% unlist()) > 1) %>% group_by(angle, size) %>% 
    summarize(n=n()) %>% arrange(desc(n)) %>% print.data.frame()
  
  cat("\n===================================\n") 
  
  cat("\n===Runs with pixel ratios but no trajectory===\n") 
  px.cm.ratios$file.name[which(!(px.cm.ratios$file.name %in% big.db$file.name))] %>% 
    print()
  
  cat("\n===Runs with trajectories but no pixel ratio information===\n") 
  big.db$file.name[which(!(big.db$file.name %in% px.cm.ratios$file.name))] %>% 
    print()
  
}

#assign_speeds also can assign count data
assign_speeds <-function(x, assigned.var = c("counts","speeds")) {
  
  #ugh, because they're stored in lists of 1?
  cricket.traj <-  x$cricket.traj
  indvd.seq <- x$indvd.seq
  ant.counts <- x$ant.counts
  
  #currently using just movement speed
  speeds <-
    cricket.traj %>% TrajFromCoords() %>% TrajSmoothSG(p=1,n=5) %>% TrajDerivatives() %>% `[`(1:2) %>% as_tibble()
  plot(speeds[[2]], speeds[[1]])

    
  
  assigner <- switch(assigned.var,
                    counts = ant.counts %>% select(Ant.Number, Frame),
                    speeds = speeds)
  
  #note, I might want to apply a ROLLAPPLY() to get more averaged speeds
  indvd.seq %>% mutate(
    start.int = cut(min, c(0, assigner[[2]]), labels=F), #this finds which cricket speed intervals each ant's starting time was in
    end.int =
      cut(max, c(0, assigner[[2]]), labels=F) %>%
      ifelse(is.na(.), length(assigner[[2]]), .),  #this is in case it goes PAST the end
    start.speed = assigner[[1]][start.int],
    end.speed = assigner[[1]][end.int],
    avg.speed = map2_dbl(start.int, end.int,
                         ~seq(.x,.y) %>% assigner[[1]][.] %>% mean()
    )# this calculates the speed points between stop and start, calculates average
  ) %>% as_tibble() -> db
  
  if(assigned.var == "counts") names(db) <- gsub("speed", "count", names(db))
  
  db %>% return()
}


file_mover <- function(type = c("cricket","sequence","ant.number")) {
  
  setwd("Z:/Users/Andrew_Burchill/Weaver ant vids")
  replete <- NULL
  replete <- list.files(recursive = T,
                        pattern = switch(type,
                                         cricket = "[msl][0-9]*.mdf",
                                         sequence = "(seq|timepoint)",
                                         ant.number = "num"))
  
  replete <- as.data.frame(replete, stringsAsFactors = F)
  colnames(replete) <- "full.name"
  
  db <- switch(type,
               cricket = replete %>% 
                 mutate(
                   angle = strapplyc(full.name, "^[[:print:]]+/([0-9]{1,2})[[:print:]]*/[[:print:]]+/[[:print:]]+\\.mdf") %>% unlist(),
                   date = strapplyc(full.name, "^([[:print:]]{3}[[:blank:]][0-9][0-9])[[:print:]]+") %>% unlist(),
                   name = strapplyc(full.name, "^[[:print:]]+/[[:print:]]+/[[:print:]]+/([[:print:]]+)\\.mdf") %>% unlist(),
                   full.name =  paste0("Z:/Users/Andrew_Burchill/Weaver ant vids/", full.name)
                 ),
               sequence = replete %>% 
                 mutate(
                   angle = strapplyc(full.name, "^[[:print:]]+/([0-9]{1,2})[[:print:]]*/[[:print:]]+/[[:print:]]+\\.mdf") %>% unlist(),
                   date = strapplyc(full.name, "^([[:print:]]{3}[[:blank:]][0-9][0-9])[[:print:]]+") %>% unlist(),
                   name = strapplyc(full.name, "^[[:print:]]+/[[:print:]]+/[[:print:]]+/([[:print:]]+)\\.mdf") %>% unlist(),
                   full.name =  paste0("Z:/Users/Andrew_Burchill/Weaver ant vids/", full.name)
                 ),
               ant.number = replete %>% 
                 mutate(
                   angle = strapplyc(full.name, "^[[:print:]]+/([0-9]{1,2})[[:print:]]*/[[:print:]]+/[[:print:]]+\\.csv") %>% unlist(),
                   date = strapplyc(full.name, "^([[:print:]]{3}[[:blank:]][0-9][0-9])[[:print:]]+") %>% unlist(),
                   name = strapplyc(full.name, "^[[:print:]]+/[[:print:]]+/[[:print:]]+/([[:print:]]+)\\.csv") %>% unlist(),
                   full.name =  paste0("Z:/Users/Andrew_Burchill/Weaver ant vids/", full.name)
                 )
  )
  
  db$date %<>%
    gsub(pattern = '([[:upper:]])', perl = TRUE, replacement = '\\L\\1', .) %>%
    gsub("\\s", "", .) 
  db %<>% mutate(destination = paste0("data/", date, "/", angle, name,
                                      switch(type,
                                             cricket = ".mdf",
                                             sequence = ".mdf",
                                             ant.number = ".csv")
  )
  )
  
  setwd("C:/Users/aburchil/Dropbox (ASU)/Weaver ant cloud storage/weaver-ant-project")
  #okay, now move the files
  db %>% apply(1, function(x) file.copy(x[1],x[5], overwrite = FALSE)) %>% db$full.name[.] %>% print()
  
  
  
}

#+ begin, include=FALSE
#___Begin########

#setwd("C:/Users/aburchil/Dropbox (ASU)/Weaver ant cloud storage/weaver-ant-project")
#setwd("D:/Weaver ants/Dropbox (ASU)/weaver ant cloud storage/weaver-ant-project")


#make a big database to hold everything
px.cm.ratios <- read.csv("data/px.cm.ratios.csv", stringsAsFactors=FALSE)
{
cricket.files <- 
  list.files(path= "data/", pattern = "[0-9]{1,2}[msl][0-9]", recursive = T) %>% 
  paste0("data/",.)
  
# Remove "troublesome" runs
  remove<- c("data/oct29/90l1.mdf",
             "data/oct29/90l1.mdf.mdf",
             "data/oct29/90m1.mdf",
             "data/oct29/90s1.mdf",
             "data/oct31/90s2.mdf",
             "data/oct23/45m2-raw-needs-adjusting.mdf")
  
  cricket.files[!(cricket.files %in% remove)] ->cricket.files
  
trajs <- 
  map(cricket.files, ~extractor(.,"cricket")) %>%  #go through all the names of the files and extract the goodies
  map(~select(.,
              X = m.x, Y = m.y, frame = h.frame,
              length, angle.from.ref, dist.from.ref, dist.from.start)) %>%   #now select only what's needed
  map(~TrajFromCoords(.,
                      xCol = "X", yCol = "Y", timeCol = "frame", fps = 5, 
                      timeUnits = "s/24", spatialUnits = "cm") %>% as_tibble()) %>%  #turn it into a TrajR format
  setNames(cricket.files) #and give em names

#make a giant database from trajectories to hold everything
big.db <- 
  tibble(file.name = "", cricket.traj = trajs) %>%
  mutate(file.name = map_chr(cricket.traj, ~attr(., "file.name")))  #give em fiiiiiile names
  
#now get the metadata (size, angle, etc etc)
big.db <- big.db$file.name %>% get_metadata("cricket") %>% bind_cols(big.db, .)

#get the stats from each run
big.db <- TrajsMergeStats(big.db$cricket.traj, stat_generator) %>%
  bind_cols(big.db, .) %>%  #below are just a few tweaks
  mutate(angle.mean =  as.circular(angle.mean, type="angles",units="degrees", modulo = "2pi"),
         size.angle = paste0(size,"-",angle) %>% as.factor())
}



## now load the attachment duration sequence data
{
seq.files <- 
  list.files(path= "data/", pattern = "[a-z]+[-][sml][0-9][-][0-9]+[.]mdf", recursive = T ) %>%
  paste0("data/",.)

indvd.sequences <-
  map(seq.files, ~extractor(.,"sequence")) %>% #extract goodies from files
  setNames(seq.files) %>% 
  bind_rows(.id="file.name") #make sure they're all identified by filename

big.db <-
  indvd.sequences %>% 
  pluck("file.name") %>%
  get_metadata("sequence") %>%  #get metadata from files
  bind_cols(indvd.sequences, .) %>% 
  help_special_snowflakes("sequence") %>%  #tidy it up, some files are weird
  group_by(angle,date,size,run.num,colony.id) %>%  #begin joining seqs and the big.db
  nest() %>%  #first, put each file's sequences into a nested DB
  full_join(big.db, ., by = c("angle","date","size","run.num","colony.id"), keep=T) %>%
  rename(indvd.seq = data) %>%  #"data" isn't informative
  mutate_at(c("date", "colony.id", "run.num"), ~as.factor(.)) %>% #make sure these are still factors
  replace(.=="NULL", NA) #for some reason, non sequenced files get NULLS, not NAs



}

## now load the attached ant COUNT data
{

ant.number.files <-
  list.files(path= "data/", pattern = "[0-9].csv", recursive = T ) %>%
  paste0("data/",.)

ant.counts <-
  map(ant.number.files, ~read.csv(., stringsAsFactors=FALSE, skip=1)) %>% 
  setNames(ant.number.files) %>% bind_rows(.id="file.name")

big.db$ant.counts<-NULL
#add the ant counts as a tibble in the database
big.db <-
  ant.counts %>% 
  pluck("file.name") %>%
  get_metadata("ant.number") %>%  #get metadata from files
  bind_cols(ant.counts, .) %>% 
  group_by(angle,date,size,run.num,colony.id) %>%  #begin joining seqs and the big.db
  nest() %>%  #first, put each file's sequences into a nested DB
  full_join(big.db, ., by = c("angle","date","size","run.num","colony.id"), keep=T) %>%
  rename(ant.counts = data) %>%  #"data" isn't informative
  mutate_at(c("date", "colony.id", "run.num"), ~as.factor(.)) %>% #make sure these are still factors
  replace(.=="NULL", NA) #for some reason, non sequenced files get NULLS, not NAs

#big.db %<>% mutate(avg.ant.num=)
}


#mark which runs start out of frame
big.db %<>% full_join(., px.cm.ratios %>% select(file.name, start.in.frame))
#only keep in-frame beginnings
big.db %<>% filter(start.in.frame == "good") 
#remove extra "NULL" entries at the bottom
big.db %<>% filter(map(cricket.traj, length %>% unlist()) > 1) 
  

#' # Weaver Ant Project {.tabset .tabset-fade}
#' 

#___Speed Stuff?#####################
#' ## Looking at Speed and Ant Number
#' 
#' This stuff is still developing
#' 
#' ### Visualizing number of ants vs. speed
#' 
#' Both speed and the number of ants are in time series for each run. I haven't counted the ants
#' in many of the runs yet, however. Location (and thus speed) was sampled
#' every 60 frames and the ant count was sampled every 240 frames (its much more difficult).
#' 
#' Here we have how the number of ants changes over time.

#+ dev-stuff
speed_smusher <- function(x) {
  
  cricket.traj <-  x$cricket.traj
  indvd.seq <- x$indvd.seq
  ant.counts <- x$ant.counts
  
  #calculate pure MOVEMENT speed
  cricket.traj %>% TrajDerivatives() %>% `[`(1:2) %>% as_tibble() %>%
    mutate(speedTimes=speedTimes+1) -> trajspeed
  
  #calculate "speed" TOWARD THE REF
  cricket.traj %>% select(dist.from.ref, time) %>%
    apply(2, diff) %>% as_tibble() %>% transmute(ref.speed=-1*dist.from.ref/time) %>% #need the negative one, because.
    #bind them together now
    bind_cols(.,trajspeed) %>%
    #add the ant counts
    right_join(.,ant.counts, by=c("speedTimes" = "Frame")) %>% 
    return()
}

big.db %>% filter(map(ant.counts, length %>% unlist()) > 1) %>%
  apply(1, function(x) speed_smusher(x)) %>% bind_rows() %>%
  { bind_cols(.,pluck(.,"file.name") %>% get_metadata("ant.number"))} %>%
  mutate(
    size.angle=paste(size,angle,sep=".") %>% as.factor(),
    size=factor(size, ordered = F)) %>% 
  na.omit() -> time.speed.count 
time.speed.count$weight <- 0
time.speed.count$weight[time.speed.count$size=="s"]<-0.25
time.speed.count$weight[time.speed.count$size=="m"]<-1.25
time.speed.count$weight[time.speed.count$size=="l"]<-2.25


grouped <- 
  time.speed.count %>% group_by(file.name, angle, size, run.num, date) %>% 
  summarize(mean=mean(ref.speed), med=median(ref.speed),
            meana=mean(Ant.Number), meda=median(Ant.Number),
            size.angle=unique(size.angle), weight=unique(weight),
            max=max(ref.speed), maxa= Ant.Number[ref.speed == max(ref.speed)])

#+ ants-over-time
big.db %>%
  filter(map(ant.counts, length %>% unlist()) > 1) %>% 
  select(ant.counts,size,angle) %>% unnest() %>%
  ggplot(aes(x=(Frame), y=Ant.Number, color=file.name)) + 
  geom_line() + facet_wrap(~ as.factor(angle):size) +
  ggtitle("Number of ants over time") +theme(legend.position="none") 

#' Here is just a crazy plot that give you a very general pattern of how speed
#' and ant count are related. (Big negative speeds are when slipping occurred)
#+ all-points, fig.cap="Figure: the points here are each time I recorded the number of ants on a trajectory (every 240 frames). Each run recorded has many points on here."
ggplot(time.speed.count, aes(x=weight/Ant.Number, y=ref.speed)) +
  geom_jitter(size=3, aes(color = angle, shape=size), width = .1, height = 0) + 
  scale_color_brewer(palette = "Dark2") + 
  ggtitle("Number of ants vs speed towards goal")  

# 
# ggplot(time.speed.count, aes(x=Ant.Number, y=ref.speed, color=size.angle)) + geom_point(size=5, aes(color = size.angle, shape=size))  + stat_ellipse()
# 
# ggplot(time.speed.count, aes(x=Ant.Number, y=ref.speed)) + 
#   geom_jitter(size=3, aes(color = angle, shape=size),  width=.2)  + 
#   stat_ellipse(mapping=aes(color=size), linetype=3, show.legend = F) +stat_ellipse(mapping=aes(color=angle), size=1) +  scale_color_brewer(palette = "Dark2")

#+ stuff, out.width=c('50%', '50%'), fig.show='hold', fig.cap="Figure: each individual run only has one point here"

ggplot(grouped, aes(x=angle:size, y=max, label=file.name)) + 
  geom_boxplot() +geom_jitter(position=position_jitter(0.2)) +
  ylim(c(0,.2)) + ggtitle("Maximum 'forward' speed") #+  geom_text(fontface = "bold")

ggplot(grouped, aes(x=angle:size, y=meda, label=file.name)) + 
  geom_boxplot() + geom_jitter(position=position_jitter(0.2)) +
  ggtitle("Median number of ants attached") #+  geom_text(fontface = "bold")

#' Hmmmm, it looks like angle doesn't affect median number of ants, at least from this figure.
#' Let's visualize it differently.

#+
ggplot(grouped, aes( y=meda, label=file.name)) + 
  geom_boxplot() + 
  facet_wrap(~ size:angle) + #geom_jitter(position=position_jitter(0.2)) +
  ggtitle("Median number of ants attached")

#' I feel like angle matters in some size classes more than others... Let's test all the options
#' in a linear model! I "dredge" all the submodels of "median.ant.count ~ weight*angle" and see
#' which are best.

#+ lms, include=TRUE, echo=TRUE
model <- 
  grouped %>% ungroup %>%
  mutate(angle=as.numeric(paste(angle))) #turn the angle factor into a number

lm(meda ~ weight*angle, data=model, na.action="na.fail") %>% MuMIn::dredge()
lm(meda ~ weight*angle, data=model, na.action="na.fail") %>%
  MuMIn::dredge() %>%
  MuMIn::get.models(subset = T) %>%
  MuMIn::importance()

#' Huzzah! Looks like angle DOES matter for the number of ants involved. It's
#' certainly not quite as important as the size of the cricket, but still.


#+ more-stuff, out.width=c('50%', '50%'), fig.show='hold', fig.cap="Figure: each individual run only has one point here"


ggplot(grouped, aes(x=weight/meana, y=mean, color=angle)) + 
  geom_point(size=5, aes(color = angle, shape=size)) + xlab("Weight/Ant.Number")+
  ggtitle("The average weight carried per ant vs average speed")  

ggplot(grouped, aes(x=(maxa), y=max, color=angle)) + 
  geom_point(size=5, aes(color = angle, shape=size))  + xlab("Ant.Number")+
  ggtitle("The number of ants attached at max speed")  

#' Pretty obviously, it looks like small ants carry little weight per ant and
#' move much much faster. Also, the more ants, the slower they go. And the more
#' weight per ants, the slower they go. 

#+ more-more-stuff



# 
# ggplot(big.db, aes( y=straightness)) + 
#   geom_boxplot() + 
#   facet_wrap(~ size:angle) + #geom_jitter(position=position_jitter(0.2)) +
#   ggtitle("Median number of ants attached") + #+  geom_text(fontface = "bold")
#   theme_pubclean()

big.db %<>% 
  mutate(
    angle = as.numeric(paste(angle)),
    weight = case_when(
  .$size == "s" ~ 0.25,
  .$size == "m" ~ 1.25,
  .$size == "l" ~ 2.25)
  )

#___Survival stuff####################
#' 

#' 
#' ## Survival Analysis 
#' 
#' So the plan was to investigate the duration that individual ants remain attached
#' to the cricket. **Do ants "hold on" harder on vertical surfaces??** However, collecting
#' this data is a huuuuuge pain: you have to manually track individual ants
#' and see when the let go of (or grab on to) the cricket. We mus t
#'  subsample within the runs and get a random subsample of the data.  But how?!
#'  

#+ attachment-size-angle, include=FALSE

#prep data for survival analysis

all.seqs <- 
  big.db %>% 
  filter(map(indvd.seq, length %>% unlist()) > 1) %>% #this filters out the list(NA)s
  filter(map(ant.counts, length %>% unlist()) > 1) %>% 
  apply(1, function(x) assign_speeds(x, "speeds")) %>%  #puts speeds on the survival data
  bind_rows() %>% 
  {bind_cols(get_metadata(.$file.name, "sequence"),.)} %>% #then re-generate the metadata and add it
  select(-seq.type1) %>% 
  rowid_to_column() #yeeeeeeeeeeeeeech. it gets confusing because multiple ants during
#the "timepoint" sequences are identical!


#now lets repeat it then add the new columns?
big.db %>% 
  filter(map(indvd.seq, length %>% unlist()) > 1) %>% #this filters out the list(NA)s
  filter(map(ant.counts, length %>% unlist()) > 1) %>% 
  apply(1, function(x) assign_speeds(x, "counts")) %>%  #puts speeds on the survival data
  bind_rows() %>% 
  {bind_cols(get_metadata(.$file.name, "sequence"),.)} %>% #then re-generate the metadata and add it
  select(-seq.type1, -start.int, -end.int) %>% 
  rowid_to_column() %>% 
  full_join(all.seqs) -> all.seqs

#' ### Sampling techniques
#' 
#' We've considered two methods of subsampling:
#' 
#' 1. **Attachment events:** Sample across attachment/detachment events 
#' over time within a run. From some early starting point, find the next ant that attaches to the cricket. 
#' track it until it detaches, wait until you see another ant attach, then follow this new ant too. 
#' Repeat! These attachment events never overlap in time.
#' 2. **The entire cohort from one random timepoint:** Sample the entire transporter cohort at a
#' given time point. Randomly choose a frame in the run, mark each ant currently attached, then 
#' follow those ants until they let go. Durations overlap in time.
#' 
#' Regardless of their merits or demerits, I collected data using both types. Let's look
#' at the difference between sampling techniques that we talked about.
#'  Do they yield the same results? 

#+ sampling, fig.align="center", fig.width=6, fig.height=6, fig.cap="Figure: Does it matter how I sample the ants?"

seqs <- all.seqs

attachments <-
  Surv(time = seqs$length,  #this is the ending AGE that I have.
       event = seqs$right.censored)  #this is whether they were still grabbing on at the end

(survfit(attachments ~seq.type, data = seqs) -> fit5) %>% 
  ggsurvplot(., data=seqs, conf.int =TRUE, 
             # risk.table = "abs_pct",  risk.table.y.text = TRUE, risk.table.col = "strata", 
             surv.median.line = "hv", break.time.by = 250, conf.int.alpha=0.15, pval=T )

#' Turns out that they aren't the same, as odd as that is. Ants chosen by being already attached at a
#' random time point hold on for longer. 
#' 
#' This corroborates what I've seen anecdotally. There seem to be a limited number of places
#' on the cricket that ants can get a good grip. The ants at those spots hold on for a long
#' time and the ants trying to attach to other spots are pretty short-term. These 
#' "didn't really grab on" ants have a high turn-over. 
#' 
#' 
#' 

#' ### Angles and Sizes on Attachment Duration
#' 
#' So, I have tons of data from different angles of inclination and different cricket sizes.
#' Let's see how attachment durations differ!



#+ inclines-sizes-vers-events, out.width=c('50%', '50%'), fig.show='hold', fig.cap="Figure: Using attachment event subsampling"

seqs <- all.seqs %>% filter(seq.type=="seq")

seqs <- seqs %>% na.omit()
attachments <-
  Surv(time = seqs$length,  #this is the ending AGE that I have.
       event = seqs$right.censored)  #this is whether they were still grabbing on at the end

(survfit(attachments ~ size, data = seqs) -> fit1) %>% 
  ggsurvplot(., data=seqs, conf.int =TRUE, 
             surv.median.line = "hv", break.time.by = 250, conf.int.alpha=0.15, pval=T )

(survfit(attachments ~ angle, data = seqs) -> fit2) %>% 
  ggsurvplot(., data=seqs, conf.int =TRUE, 
             surv.median.line = "hv", break.time.by = 250, conf.int.alpha=0.15, pval=T )

#' Whelp. As you can see, neither cricket size nor angle (or the interaction) affect survival here.
#' What about with time-point subsampling?


#+ inclines-sizes-vers-slice, out.width=c('50%', '50%'), fig.show='hold', fig.cap="Figure: Using time-slice event subsampling"

seqs <- all.seqs %>% filter(seq.type=="timepoint")

seqs <- seqs %>% na.omit()
attachments <-
  Surv(time = seqs$length,  #this is the ending AGE that I have.
       event = seqs$right.censored)  #this is whether they were still grabbing on at the end

(survfit(attachments ~ size, data = seqs) -> fit1) %>% 
  ggsurvplot(., data=seqs, conf.int =TRUE, 
             surv.median.line = "hv", break.time.by = 250, conf.int.alpha=0.15, pval=T )

(survfit(attachments ~ angle, data = seqs) -> fit2) %>% 
  ggsurvplot(., data=seqs, conf.int =TRUE, 
             surv.median.line = "hv", break.time.by = 250, conf.int.alpha=0.15, pval=T )

#' Yeah, neither seems to show much...
#' 
#' <hr />
#' 
#' **However**, I also can calculate the (approximate) speed that the cricket was moving when these ants
#' attach, release, and the average speed during the attachment time (I have the whole time series).
#'  I *also* have the count data (number of ants attached) over time as well, so I can smush that 
#'  trajectory data with this attachment data and do some "regression" to see
#'   if speed, number of ants, or anything else affects attachment duration.
#'   
#'   I'll use speed as an example, because I already have it type up, and because there are soooo
#'   many posibilities here.
#' 
#' But before I do that, I need to have my speed data looking pretty.
#' 
#' ### Speed Data and which transformation to use? {.tabset}
#'  
#'  Should I just use the normal speeds? Should I log them? Not sure...
#'  The covariates are right-skewed.
all.seqs %>% filter(seq.type=="seq") %>% .$start.speed %>%
  density %>% plot(main="starting speed distribution")
#'  But I'm not sure it needs to be normal for these survival tests.
#'  

#' #### It's probably boring, don't read
#' 
#' I don't want to get bogged down in these details. I just picked log, and it seems to work.
#'   
#' #### Log, sqrt, nothing
#' 
all.seqs %>% filter(seq.type=="seq") %>% .$start.speed %>% log() %>% 
  density %>% plot(main="logged starting speed distribution")
all.seqs %>% filter(seq.type=="seq") %>% .$start.speed %>% sqrt() %>% 
  density %>% plot(main="logged starting speed distribution")
#' Square root seems to make it the most normal... BUT!
#' 
#' #### Cox assumptions for proportional hazards
#' 
#' Uh, the "martingale residuals" should apparently be roughly linear to satisfy the proportional hazards
#' assumption, so we need to pick which transformation
#' of the data makes the loess line straightest?

ggcoxfunctional(attachments ~ start.speed + log(start.speed) + sqrt(start.speed), data = seqs)
ggcoxfunctional(attachments ~ avg.speed + log(avg.speed) + sqrt(avg.speed), data = seqs)
ggcoxfunctional(attachments ~ end.speed + log(end.speed) + sqrt(end.speed), data = seqs)





#' Yeech, they all look bad, but log makes things more significant? Also it's more commonly used and
#' interpretable?
#' 
#' <hr />
#' 
#' ### "Regressing" with the covariates
#' 
#' Unlike normal regression, where it seems relatively straightforward to find out which factors
#' are significant, the interactions between variables are less "linear" here. None of those three factors are significant by
#' themselves, but in combination, two or more often both become significant. I also don't know how to
#' interpret the "causality" of *ending speed* on the length of time the ant stays attached, so I just 
#' use the average speed and the starting speed.

(coxph(attachments ~ log(avg.speed) + log(start.speed) , data=seqs) -> fit3) %>% print()
fit3 %>% ggforest(data=seqs)

#' Yay it looks significant and cool! However..
#' 

cox.zph(fit3) 

#' The *big* important assumption of proportional hazards (HRs) hasn't been met...
#' 
#' Turns out the data def-o isn't going to meet the requirements, hazards aren't constant over
#' attachment time.


#+ time-stuff

(coxph(attachments ~ log(start.speed):length , data=seqs) -> fit3) %>% print()

#' Notice how just the interaction term of start speed and total time is hellllla significant.
#' 
#' #### We have a few options in this situation though. 
#' 
#' 1. Ignore the violated assumption -> Ehhhh, apparently there are some issues in interpretation here.
#' 2. Estimate piece-wise HRs -> It does look like around t=110, something different may be happening.
#' 3. Include a time-by-covariate interaction -> Boooo! Too many options, too hard to interpret.
#' 4. Calculate a population-averaged HR -> Easy to understand, robust to stuff, huzzah!
#' 

#ggcoxdiagnostics(fit, type = "dfbetas", linear.predictions = TRUE)
#https://stat.ethz.ch/pipermail/r-help/2011-June/282065.html
#https://www.jstatsoft.org/article/view/v084i02/v84i02.pdf

#+ coxphw, background="#03A678"
coxphw(Surv(time = seqs$length,  #this is the ending AGE that I have.
            event = seqs$right.censored) ~  log(avg.speed)+log(start.speed), data=seqs, template = "AHR") -> fit4
fit4 %>% summary()

#' Okay, so interpretting the generalized concordance probability... If the speeds were categorical
#' variables, you could say that the an ant with a "higher" average speed would let go earlier than
#' an ant with a "lower" average speed with a probability of 
{{concord(fit4)[1]*100}}
#', and if the ant has a "higher" starting speed, its probability of letting go earlier than a 
#' lower-start-speed ant is
{{concord(fit4)[2]*100}}
#' . 

#' 
#' <hr />
#' 

#' 

#___Making Cut-offs-------------

#' ##Fixing out-of-frame beginnings
#' 
#' Since some of the trajectories start with the cricket OUT of frame, I miss
#' their "ramp-up" period, which is probably different than
#' the rest of the run. Damn. And this happens/happened more with 90 and 45 degree
#' inclines, because the cricket would fall/roll out of the frame.
#' 
#+ plots, echo= FALSE, warning= FALSE

initial_speed <- function(x) {
  #gets the first one [OR TWO] "speed" points and averages them, then adds it to a tibble
  x$cricket.traj %>% TrajFromCoords() %>% 
    TrajDerivatives() %>% pluck("speed") %>% .[1] %>% #:2] %>% mean() %>%
    as_tibble() %>% 
    mutate(file.name = x$cricket.traj %>% attr("file.name")) %>% 
    return()
}

big.db %>% apply(1, function(x) initial_speed(x)) %>% bind_rows %>% arrange(desc(value)) -> init.speed
init.speed %>% { bind_cols(.,pluck(.,"file.name") %>% get_metadata("cricket"))} -> init.speed

init.speed %>% 
  group_by(angle,size) %>%
  summarize(number=n(), percent.bad= sum(value>0.04)/n()) %>% 
  knitr::kable(caption = "Proportion of runs starting faster than 0.04", digits=2) %>% 
  kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),full_width = F)
#' Below I just plot the initial starting speeds of the trajectories I have. As you can see,
#' it looks like there are some outliers; they are mostly likely ones where the 
#' cricket starts out of the frame.

ggplot(init.speed, aes(x=value, fill=size)) + geom_histogram(bins=40) +
  facet_grid(size~.) + 
  geom_vline(data = filter(init.speed, size=="s"), aes(xintercept=0.04)) +
  geom_vline(data = filter(init.speed, size=="m"), aes(xintercept=0.02)) +
  geom_vline(data = filter(init.speed, size=="l"), aes(xintercept=0.01)) +
  ggtitle("Frequencies of starting speeds")

#' Okay, so I have two options, remove the videos that start out of frame...

#+ lots-plots, out.width=c('50%', '50%'), fig.show='hold', fig.height=3
thresh<-0.04
for (i in subset(big.db, size=="s")$cricket.traj) {
  
  
  i %>% TrajFromCoords() %>% TrajDerivatives() %>% pluck("speed") -> speeds
  plot(speeds, type="l", main= attr(i, 'file.name'))
  points(speeds, main= attr(i, 'file.name'))
  abline(h=thresh, col="blue")
  abline(v=min( which(speeds>thresh))-1, col="blue")
  #THIS IS SUCH A DUMB WAY TO DO IT
  lines(
    c(rep(-10,min(which(speeds>thresh))-1),
      speeds[(min(which(speeds>thresh))-1):length(speeds)+1]
      ), col="red")
  
}


# i[which(i$dist.from.start > di),]
# 
# 
inner_join(big.db,
           grouped %>% select(-file.name), by=c("angle","size","run.num", 'date')) -> wow
wow %>% colnames()


df.wow <- wow[c(8,9,11,10,13,15,24,29)]
autoplot(prcomp(df.wow, scale. = T), data=wow, colour="angle", shape="size",
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3, size=5, x=3,y=2)

plot(wow$size.angle.y, prcomp(df.wow, scale. = T)$x[,1])

lm( prcomp(df.wow, scale. = T)$x[,1] ~ wow$angle:wow$size+wow$size , REML = F) %>% MuMIn::AICc()
lmerTest::lmer( prcomp(df.wow, scale. = T)$x[,1]
                ~ wow$angle*wow$size + (1|wow$colony.id), na="na.fail") %>%
  MuMIn::dredge()  %>% attributes() %>% .$model.calls %>% unname() %>% unlist()

#___Session info###############

#' # Session info
#' 
#+ session.info
sessionInfo()


#___Style and CSS##############

#' <style>
#' body {
#' color: #594F57;
#' font-family: 'Vollkorn', "Source Sans Pro","Helvetica Neue",Helvetica,Arial,sans-serif;
#' }
#' 
#' h1, .h1, h2, .h2, h3, .h3 {
#'   margin-top: 50px;
#' }
#' 
#' h1 {
#'   padding-top: 60px;
#' }
#' .hljs {
#'   color: #03A678;
#' }
#' 
#' h1.title {
#'   text-align: center;
#' }
#' 
#' p.caption {
#'   font-size: 0.9em;
#'   font-style: italic;
#'   color: grey;
#'   margin-right: 10%;
#'   margin-left: 10%;  
#'   text-align: center;
#' }
#' 
#' pre {
#'   border: 0px solid #cccccc;
#'   border-radius: 0px;
#' }
#' 
#' 
#' a {
#'     color: #03A678;
#'   }
#' a:hover {
#'   color: #03A678;
#' }
#' 
#'   #TOC {
#'   top: 20%;
#' position: fixed;
#' width: 210px;
#' }
#' .tocify {
#'   border: none;
#'   color: grey;
#'   border-radius: 0px;
#' }
#' .list-group-item.active {
#'   color: #03A678;
#'     background-color: #FBF9FA;
#'   border-left: solid;
#'   border-color: #03A678;
#' }
#' .list-group-item:hover {
#'   color: #03A678;
#'     background-color: #FBF9FA;
#'   border-left: solid;
#'   border-color: #03A678;
#' }
#' .tocify-extend-page {
#'   display: none;
#' }
#' 
#' 
#'   .toc-content {
#'     margin-left: 210px
#'   }
#'   
#' @media (max-width: 1200px) {
#' #TOC {
#' width: auto;
#' display: none;
#' }
#' .toc-content {
#' margin-left: auto
#' }
#' }
#' 
#' 
#' 
#' .btn {
#' border-width: 0 0px 0px 0px;
#' font-weight: normal;
#' text-transform: ;
#' color: #03A678;
#' }
#' .btn-default {
#' color: #03A678;
#' background-color: #FBF9FA;
#' border-color: #FBF9FA;
#' }
#' 
#' .nav-pills > li.active > a, .nav-pills > li.active > a:hover, .nav-pills > li.active > a:focus {
#' background-color: #03A678;
#' }
#'  
#' #TOC {
#' background: url("https://andburch.github.io///images/ant.png");
#' background-size: contain;
#' padding-top: 80px !important;
#' background-repeat: no-repeat;
#' }
#' </style>

#_______END_______####
#_______END_______####

##########################################
##########################################

 
# #-------------------------------------------------
# for (i in seq_len(dim(hog)[1])){
#   
#   
#   
#   cricket.traj <-  hog$cricket.traj[[i]]
#   indvd.seq <- hog$indvd.seq[[i]]
#   
#   #currently using just movement speed
#   speeds <- 
#     cricket.traj %>% TrajFromCoords() %>% TrajDerivatives() %>% `[`(1:2) %>% as_tibble()
#   
#   #note, I might want to apply a ROLLAPPLY() to get more averaged speeds
#   indvd.seq %>% mutate(
#     start.int = cut(min, c(0, speeds[[2]]), labels=F), #this finds which cricket speed intervals each ant's starting time was in
#     end.int = 
#       cut(max, c(0, speeds[[2]]), labels=F) %>%
#       ifelse(is.na(.),length(speeds[[2]]),.),
#     start.speed = speeds[[1]][start.int],
#     end.speed = speeds[[1]][end.int],
#     avg.speed = map2_dbl(start.int, end.int,
#                          ~seq(.x,.y) %>% speeds[[1]][.] %>% mean()
#     ) # this calculates the speed points between stop and start, calculates average
#   ) %>% as_tibble() %>% print()
# }
# 
# 
# #--------------------
# 
# stats2<- stats_subsetter(stats, ord=F)
# 
# 
# 
# stats_subsetter <- function(stats, ord = T, ...) {
#   stats2 <- stats
#   
#   stats2 <- subset(stats2, size != "l") %>% mutate(size = ordered(size, levels=c("s","m")))
#   stats2$size.angle <- paste0(stats2$size,"-",stats2$angle) %>% as.factor()
#   
#   stats2 <- stats2[which(stats2$date %in% c("oct29", "oct30", "oct31", "oct25", "oct24", "oct23")),]
#   
#   stats2<-stats2 %>%
#     mutate(angle = factor(angle, ordered = ord),
#            size = factor(size, ordered =ord))
#   return(stats2)
# }
#   
# lmer(total.avg.speed ~ angle*size + (1|colony.id/date), data=stats2[-6,]) %>% AIC
# 
# 
# 
# 
# 
# 
# 
# 
# a <- 0
# for (i in big.db$cricket.traj) {
#   
#   if (120 %in% diff(i$time)) {
#     print("DO THIS ONE!")
#     print(attr(i,"file.name"))
#     print(dim(i)[1])
#     print((i$time))
#     i$time %>% diff() %>% print()
#     a<-a+1}
#   
#   #print("---------------------")
# 
# }
# 
# 
# 
# 
# real_speed <- function(traj) {
#   
#   total.x.disp <- sum(traj$displacement) %>% Re()
#   signed.x.disp <- (total.x.disp + Re(traj$displacement)) %>%
#     abs() %>%
#     subtract(abs(total.x.disp))
#   return(signed.x.disp)
#   
#   
#   
# }
# 
# #this is a made-up cutoff for how many PIXELS they have to move before... something
# di<-100
# a <- NULL
# for (i in trujs) {
#   i %>% real_speed() %>% rollapply(w=3, mean, partial = T) %>%
#     plot(type="l", main= attr(i, 'file.name'), ylim=c(-400,400))
#   
#   abline(h=0, col="red")
#   abline(v=which(i$dist.from.start > di)[1], col="blue")
#   i %>% TrajFromCoords() %>% TrajDerivatives() %>% pluck("speed") %>% plot(type="l")
#   i %>% TrajFromCoords() %>% TrajDerivatives() %>% pluck("speed") %>%
#     acf(main=attr(i, "file.name"), demean = T, lag.max = 100, ci.type="white")
#   i %>% TrajFromCoords() %>% plot(col="red", main= attr(i, 'file.name'))
#   #i[which(i$dist.from.start > di),] %>% TrajFromCoords() %>% lines()
#   #turning <- i$angle.from.ref*pi/180
#   #CircStats::rose.diag(turning, bins=12)
#   
#   i %>% TrajFromCoords() %>% TrajDerivatives() %>% pluck("speed") %>%
#     acf(main=attr(i, "file.name"), demean = T, lag.max = 100, plot=FALSE) -> ts.acf
#   alpha <- 0.95
#   conf.lims <- c(-1,1)*qnorm((1 + alpha)/2)/sqrt(ts.acf$n.used)
#   (conf.lims[1] < ts.acf$acf & ts.acf$ acf< conf.lims[2]) %>% 
#     not() %>% which() %>% max() -> max.time
#   paste0(attr(i, "file.name"),".  ",
#          "max lag: ", max.time,
#          ",  percent series: ", (max.time/ts.acf$n.used*100) %>% round()) %>% print()
#   
#   a[attr(i, "file.name")]<-(max.time/ts.acf$n.used*100)
#   }
# 
# 
# 
# for (i in seq_along(trujs)) {
#   trujs[[i]] <- trujs[[i]][which(trujs[[i]]$dist.from.start > di),]
#   
#   
# }
# 
# 
# alpha <- 0.95
# conf.lims <- c(-1,1)*qnorm((1 + alpha)/2)/sqrt(ts.acf$n.used)
# 
# 
# 
# 
# 
# for (i in trujs) {
#   dim(i)[1] %>% print()
# }
# 
# 
# #-------------------------------------------------------------
# 
# ##############################################################
# 
# set.seed(43)
# x <- rnorm(50)
# y <- rnorm(50)
# angles <- runif(50, min=-pi, max=pi)
# Plot the locations:
#   
#   plot(x, y, pch=19, cex=0.8, col="Blue")
# Add arrows to show the orientations at these points:
#   
#   length <- 0.2
# arrows(x, y, x1=x+length*cos(angles), y1=y+length*sin(angles), 
#        length=0.05, col="Gray")
# 
# 
