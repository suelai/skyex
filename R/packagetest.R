
#' @title SkyEx-F labeling
#' @description Pair Labeling using SkyEx-F
#'

#' @param data A dataframe of pairs
#' @param p A preference function (e.g. p<-high("sim1")*high("sim2"))
#' @param label Name of the label column (e.g. "Class")
#' @param posclass How is the positive class expressed (e.g. 1)
#' @param negclass How is the negative class expressed (e.g. 0)
#'
#' @return of A skyexf object
#'
#' @export


skyexf <- function (data, p, label, posclass, negclass){
  if(is.null(p)) {
    stop('the preference is NULL')
  }

  if(is.null(label)) {
    stop('select the column that contains the label')
  }

  if(is.na(match (label, names(data)))){
    stop('enter a valid column name for the variable label')
  }

  if(!is.element(posclass, data[[label]])){
    stop(paste0("There is no ", posclass, " in ", label))
  }

  if(!is.element(negclass, data[[label]])){
    stop(paste0("There is no ", negclass, " in ", label))
  }

  data$rrow <- seq.int(nrow(data))

  res<- psel(data, p, top  = nrow(data))

  df <- data.frame(levels = numeric(), positives = numeric(),
                   precision = numeric(), recall= numeric(),
                   stringsAsFactors = FALSE)

  k = 1
  while (k<=max(res$.level)) {
    tempres<-subset(res, res$.level<=k)
    df <- rbind(df, data.frame(levels = k,
                               positives = as.numeric(nrow(tempres)),
                               precision = as.numeric(nrow(subset(tempres, tempres[[label]]==posclass)))/as.numeric(nrow(tempres)),
                               recall = as.numeric(nrow(subset(tempres, tempres[[label]]==posclass)))/as.numeric(nrow(subset(res, res[[label]]==posclass)))


    ))
    k=k+1
  }

  df$fmeasure <- (2 * df$precision * df$recall / (df$precision + df$recall))

  cut<-df[which.max(df$fmeasure),1]
  res$pred_class <- ifelse(res$.level<=cut,yes = posclass, no=negclass)
  res <- res[order(match(res$rrow, data$rrow)), ]
  of <- list(classes = res$pred_class, analysis = df,
            k = df[which.max(df$fmeasure),1], fmeasure = max(df$fmeasure))
  class(of) <- "skyexf"
  return (of)
}


#' @title Plot SkyEx-F cut-offs
#' @description Plot different SkyEx-F cut-offs and the corresponding F-measure

#' @param skyexf.object A skyexf ovject

#'
#' @export


plot.skyexf.cutoffs <- function (skyexf.object, metric, xlab, ylab) {
  if(missing(metric)){
    metric="fmeasure"
  }
  if(missing(xlab)){
    xlab="k"
  }
  if(missing(ylab)){
    ylab=metric
  }
  if(identical(metric,"all")){

    plot(skyexf.object$analysis$levels, skyexf.object$analysis$fmeasure,
         xlab = 'k', ylab = 'metric',
         ylim=c(0,1), type='l', col=6,
         cex.axis=1.3, cex.lab=1.3, lwd=2)
    abline(v = skyexf.object$k, lty = 2, col=2)
    lines(skyexf.object$analysis$levels, skyexf.object$analysis$precision, col=3, lty=2, lwd=2)
    lines(skyexf.object$analysis$levels, skyexf.object$analysis$recall, col=4, lty=3, lwd=2)
    axis(3, at=skyexf.object$k,labels=skyexf.object$k, col.axis="red", las=1, cex.axis=0.9)
    legend(0.75*max(skyexf.object$analysis$levels), 0.92, legend=c("F-measure", "Precision",  "Recall"),
           col=c(6, 3, 4), lty=1:3, cex=0.9)

  }
  else{

  plot(skyexf.object$analysis$levels, skyexf.object$analysis[[metric]],
       xlab = xlab, ylab = ylab, type='l', cex.axis=1.3, cex.lab=1.3)
  abline(v = skyexf.object$k, lty = 2, col=2)
  p1<-c(skyexf.object$k,skyexf.object$analysis[skyexf.object$k,][[metric]])
  points(t(p1), pch=16, col=2)
  text(t(p1), paste0("(",skyexf.object$k,",",
                     round(skyexf.object$analysis[skyexf.object$k,][[metric]],2),")"), col=2, pos = c(4,2))
  }
}




#' @title SkyEx-D labeling
#' @description Pair Labeling using SkyEx-D

#' @param data A dataframe of pairs
#' @param p A preference function (e.g. p<-high("sim1")*high("sim2"))
#' @param simlist Similarities to be taken into consideration for the density function (e.g. simlist= c("sim1", "sim2"))
#' @param smooth.coefficient A smoothing coefficient
#' @param posclass How is the positive class expressed (e.g. 1)
#' @param negclass How is the negative class expressed (e.g. 0)
#'
#' @return od A skyexd object
#'
#' @export
#'

skyexd <- function (data, p, simlist, smooth.coefficient, posclass, negclass){
  if(is.null(p)) {
    stop('the preference is NULL')
  }


  data$rrow <- seq.int(nrow(data))

  options(rPref.parallel = TRUE)
  res<- psel(data, p, top  = nrow(data))

  df <- data.frame(levels = numeric(), positives = numeric(),
                   meandist = numeric(),
                   stringsAsFactors = FALSE)



  k = 1
  while (k<=max(res$.level)) {
    dist<-rdist(subset(res[simlist], res$.level==k),
                subset(res[simlist], res$.level>k))

    tempres<-subset(res, res$.level<=k)
    df <- rbind(df, data.frame(levels = k,
                               positives = as.numeric(nrow(tempres)),
                               meandist=mean(dist)

    ))
    k=k+1
  }

  df$sumdist<-df$meandist*df$positives
  df$cumsumdist<-as.numeric(cumsum(as.numeric(df[order(df$levels), "sumdist"])))
  df$CumDist<-df$cumsumdist/df$positives

  n <- length(df$levels)
  fdx <- vector(length = n)

  for (i in 2:n) {
    fdx[i-1] <- (df$CumDist[i-1] - df$CumDist[i]) / (df$levels[i-1] - df$levels[i])
  }

  fdx[n] <- (df$CumDist[n] - df$CumDist[n - 1]) / (df$levels[n] - df$levels[n - 1])

  df$deriv <- fdx

  smoothderiv<-smth(x=as.numeric(df$deriv), window = smooth.coefficient, method="gaussian")

  kd<-which(smoothderiv<0.0)[1]

  res$pred_class <- ifelse(res$.level<=kd,yes = posclass, no=negclass)
  res <- res[order(match(res$rrow, data$rrow)), ]
  df$smooth.deriv <-smoothderiv

  od <- list(classes = res$pred_class, analysis = df[c("levels", "positives", "deriv", "smooth.deriv")],
            k = kd)
  class(od) <- "skyexd"
  return (od)

}

#' @title Plot SkyEx-D cut-offs
#' @description Plot the first derivative of the density function for different SkyEx-D cut-offs

#' @param skyexd.object A skyexd object

#'
#' @export


plot.skyexd.cutoffs <- function(skyexd.obj, xlab, ylab){
  if(missing(xlab)){
    xlab='k'
  }
  if(missing(ylab)){
    ylab=expression(paste(mu, "'"[k]))
  }
  plot(skyexd.obj$analysis$levels,skyexd.obj$analysis$smooth.deriv, type='l',  lwd=1,
       xlab = xlab, ylab=ylab, cex.axis=1.3, cex.lab=1.3)
  abline(h=0.0, lty = 2, lwd=1, col=1)
  abline(v = skyexd.obj$k, lty = 1, lwd=1, col=2)

}


#' @title Plot SkyEx-D for different smoothing coefficients
#' @description Plot the first derivative of the density function for different SkyEx-D cut-offs

#' @param skyexd.object A skyexd object
#' @param smooth.coefficient A value for the smoothing coefficient
#'
#' @export

plot.skyexd.smooth <- function (skyexd.obj, smooth.coefficient)
{

  smoothderiv<-smth(x=as.numeric(skyexd.obj$analysis$deriv), window = smooth.coefficient, method="gaussian")

  kd<-which(smoothderiv<0.0)[1]

  plot(skyexd.obj$analysis$levels, smoothderiv, type='l',  lwd=1,
       xlab = 'k', ylab=expression(paste(mu, "'"[k])), cex.axis=1.3, cex.lab=1.3)
  abline(h=0.0, lty = 2, lwd=1, col=1)
  abline(v = kd, lty = 1, lwd=1, col=2)
  axis(3, at=kd,labels=kd, col.axis="red", las=1, cex.axis=0.9)
}

#' @title Evaluate SkyEx algorithms
#' @description Calculates precision, recall and F-measure

#' @param prediction A vector with the predicted classes
#' @param labels A vector of the actual classes
#' @param posclass How is the positive class expressed (e.g. 1)

#'
#' @export
evaluate.skyex <- function(prediction, labels, posclass){
  ev <- data.frame(prediction, labels)
  prec = nrow(subset(ev, ev$prediction==posclass
                     & ev$labels==posclass))/nrow(subset(ev, ev$prediction==posclass))
  rec =  nrow(subset(ev, ev$prediction==posclass
                     & ev$labels==posclass))/nrow(subset(ev, ev$labels==posclass))
  f = (2*prec*rec)/(prec+rec)
  e <- list(precision = prec, recall = rec,
            fmeasure = f)
  class(e) <- "eval"
  return (e)
}


#' @title Plot SkyEx pairs 2D
#' @description Two dimentional plotting of the pairs and marking the true positives (TP), true negative (TN), false positives (FP), and false negatives (FN).

#' @param data A dataframe of pairs
#' @param sim1 Name of the column in the x axis
#' @param sim2 Name of the column in the y axis
#' @param prediction A vector of the predicted classes
#' @param labels A vector of the actual labels
#' @param posclass How the positive class is expressed (e.g. posclass=1)
#' @param colTP Color of the true positives
#' @param colTN Color of the true negatives
#' @param colFP Color of the false positives
#' @param colFN Color of the false negatives
#' @param legend A boolean argument if the legend should be present or not (e.g. legend=TRUE)
#' @param leg.x The x coordinate of the legend
#' @param leg.y The y coordinate of the legend
#' @param leg.font The font size of the legend
#' @param title The title of the plot
#' @param xlab Naming the x-axis
#' @param ylab Naming the y-axis

#'
#' @export


plot.pairs2D <- function(data, sim1, sim2, prediction, labels, posclass, colTP, colTN,
                         colFP, colFN, legend, leg.x, leg.y, leg.font,
                         title, xlab, ylab){

  if(missing(colTP)) {
    colTP<-"lightpink"
  }
  if(missing(colTN)) {
    colTN<-"lightskyblue1"
  }
  if(missing(colFP)) {
    colFP<-"red4"
  }
  if(missing(colFN)) {
    colFN<-"blue3"
  }
  if(missing(title)) {
    title<-"SkyEx Labeling"
  }

  if(missing(xlab)) {
    xlab<-sim1
  }

  if(missing(ylab)) {
    ylab<-sim2
  }

  if(missing(legend)) {
    legend<-TRUE
  }
  if(missing(leg.x)){
    leg.x<-min(data[[sim1]])
  }
  if(missing(leg.y)){
    leg.y<-max(data[[sim2]])
  }
  if(missing(leg.font)){
    leg.font<-0.6
  }
  if(leg.x<0 | leg.y<0){
    par(xpd=TRUE)
  }
  else{
    if(leg.x>max(data[[sim1]]) | leg.y > max(data[[sim2]]) )
    {
      par(xpd=TRUE)
    }
    else
    {
      par(xpd=FALSE)
    }
  }
  if(legend){

    plot(data[[sim1]], data[[sim2]], col=ifelse(prediction==posclass,
                                                ifelse(labels==posclass, colTP, colFP),
                                                ifelse(labels==posclass, colFN, colTN)),
         pch=ifelse(prediction==labels, 1, 19),
         xlab=xlab, ylab=ylab, main = title)
    legend(leg.x,leg.y,c("TP","TN","FP", "FN"),
           col = c(colTP, colTN, colFP, colFN), pch=c(1,1,19,19),
           cex=leg.font)
  }
  else
  {

    plot(data[[sim1]], data[[sim2]], col=ifelse(prediction==posclass,
                                                ifelse(labels==posclass, colTP, colFP),
                                                ifelse(labels==posclass, colFN, colTN)),
         pch=ifelse(prediction==labels, 1, 19),
         xlab=xlab, ylab=ylab, main = title)
  }
}




#' @title Plot SkyEx pairs 3D
#' @description Three dimentional plotting of the pairs and marking the true positives (TP), true negative (TN), false positives (FP), and false negatives (FN).

#' @param data A dataframe of pairs
#' @param sim1 Name of the column in the x axis
#' @param sim2 Name of the column in the y axis
#' @param sim3 Name of the column in the z axis
#' @param prediction A vector of the predicted classes
#' @param labels A vector of the actual labels
#' @param posclass How the positive class is expressed (e.g. posclass=1)
#' @param colTP Color of the true positives
#' @param colTN Color of the true negatives
#' @param colFP Color of the false positives
#' @param colFN Color of the false negatives
#' @param legend A boolean argument if the legend should be present or not (e.g. legend=TRUE)
#' @param leg.x The x coordinate of the legend
#' @param leg.y The y coordinate of the legend
#' @param leg.font The font size of the legend
#' @param title The title of the plot
#' @param xlab Naming the x-axis
#' @param ylab Naming the y-axis
#' @param zlab Naming the z-axis

#'
#' @export


plot.pairs3D <- function(data, sim1, sim2, sim3, prediction, labels, posclass,
                         colTP, colTN,
                         colFP, colFN, legend, leg.x, leg.y, leg.font,
                         title, xlab, ylab, zlab){
  if(missing(colTP)) {
    colTP<-"lightpink"
  }
  if(missing(colTN)) {
    colTN<-"lightskyblue1"
  }
  if(missing(colFP)) {
    colFP<-"red4"
  }
  if(missing(colFN)) {
    colFN<-"blue3"
  }
  if(missing(title)) {
    title<-"SkyEx Labeling"
  }

  if(missing(xlab)) {
    xlab<-sim1
  }

  if(missing(ylab)) {
    ylab<-sim2
  }

  if(missing(zlab)) {
    zlab<-sim3
  }

  if(missing(legend)) {
    legend<-TRUE
  }
  if(missing(leg.x)){
    leg.x<-mean(data[[sim1]])
  }
  if(missing(leg.y)){
    leg.y<-mean(data[[sim2]])
  }

  if(missing(leg.font)){
    leg.font<-0.6
  }
  if(leg.x<0 | leg.y<0){
    par(xpd=TRUE)
  }
  else{
    if(leg.x>max(data[[sim1]]) | leg.y > max(data[[sim2]]) )
    {
      par(xpd=TRUE)
    }
    else
    {
      par(xpd=FALSE)
    }
  }

  if(legend){
    scatter3D(data[[sim1]], data[[sim2]], data[[sim3]], colvar = NULL,
              col=ifelse(prediction==posclass,
                         ifelse(labels==posclass, colTP, colFP),
                         ifelse(labels==posclass, colFN, colTN)),
              pch=ifelse(prediction==labels, 1, 19), bty='f',
              phi = 0, xlab=xlab, ylab=ylab, zlab=zlab,
              main = title)

    legend(leg.x,leg.y, col = c(colTP, colTN, colFP, colFN),
           pch=c(1,1,19,19),
           c("TP","TN","FP", "FN"),
           cex=leg.font)
  }
  else
  {
    scatter3D(data[[sim1]], data[[sim2]], data[[sim3]], colvar = NULL,
              col=ifelse(prediction==posclass,
                         ifelse(labels==posclass, colTP, colFP),
                         ifelse(labels==posclass, colFN, colTN)),
              pch=ifelse(prediction==labels, 1, 19), bty='f',
              phi = 0, xlab=sim1, ylab=sim2, zlab=sim3,
              main = title)
  }
}



#' @title Interactive Plot SkyEx pairs 3D
#' @description Three dimentional interactive plotting of the pairs and marking the true positives (TP), true negative (TN), false positives (FP), and false negatives (FN).

#' @param data A dataframe of pairs
#' @param sim1 Name of the column in the x axis
#' @param sim2 Name of the column in the y axis
#' @param sim3 Name of the column in the z axis
#' @param prediction A vector of the predicted classes
#' @param labels A vector of the actual labels
#' @param posclass How the positive class is expressed (e.g. posclass=1)
#' @param colTP Color of the true positives
#' @param colTN Color of the true negatives
#' @param colFP Color of the false positives
#' @param colFN Color of the false negatives
#' @param legend A boolean argument if the legend should be present or not (e.g. legend=TRUE)
#' @param leg.x The x coordinate of the legend
#' @param leg.y The y coordinate of the legend
#' @param leg.font The font size of the legend
#' @param title The title of the plot
#' @param xlab Naming the x-axis
#' @param ylab Naming the y-axis
#' @param zlab Naming the z-axis

#'
#' @export


plot.pairs.interactive.3D <- function(data, sim1, sim2, sim3,
                                      prediction, labels, posclass,
                                      colTP, colTN,
                                      colFP, colFN,
                                      legend, leg.x, leg.y,
                                      leg.font, title,
                                      xlab, ylab, zlab)
{
  if(missing(colTP)) {
    colTP<-"lightpink"
  }
  if(missing(colTN)) {
    colTN<-"lightskyblue1"
  }
  if(missing(colFP)) {
    colFP<-"red4"
  }
  if(missing(colFN)) {
    colFN<-"blue3"
  }
  if(missing(title)) {
    title<-""
  }

  if(missing(xlab)) {
    xlab<-sim1
  }

  if(missing(ylab)) {
    ylab<-sim2
  }

  if(missing(zlab)) {
    zlab<-sim3
  }

  if(missing(legend)) {
    legend<-TRUE
  }
  if(missing(leg.x)){
    leg.x<-min(data[[sim1]])
  }
  if(missing(leg.y)){
    leg.y<-max(data[[sim2]])
  }

  if(missing(leg.font)){
    leg.font<-1
  }
  if(leg.x<0 | leg.y<0){
    par(xpd=TRUE)
  }
  else{
    if(leg.x>max(data[[sim1]]) | leg.y > max(data[[sim2]]) )
    {
      par(xpd=TRUE)
    }
    else
    {
      par(xpd=FALSE)
    }
  }

  if(legend){
    plot3d(data[[sim1]], data[[sim2]], data[[sim2]],
           col=ifelse(prediction==posclass,
                      ifelse(labels==posclass, colTP, colFP),
                      ifelse(labels==posclass, colFN, colTN)),
           pch=ifelse(prediction==labels, 1, 19),
           xlab = xlab, ylab=ylab, zlab=zlab, main=title)

    legend3d(leg.x,leg.y, col = c(colTP, colTN, colFP, colFN),
             pch=c(1,1,19,19),
             c("TP","TN","FP", "FN"),
             cex=leg.font)
  }
  else
  {
    plot3d(data[[sim1]], data[[sim2]], data[[sim2]],
           col=ifelse(prediction==posclass,
                      ifelse(labels==posclass, colTP, colFP),
                      ifelse(labels==posclass, colFN, colTN)),
           pch=ifelse(prediction==labels, 1, 19),
           xlab = sim1, ylab=sim2, zlab=sim3, main=title)
  }

}



#' @title Textual blocking
#' @description Creates blocks of entities that have textual similarity. Returns the pairs
#'

#' @param data A dataframe of entities
#' @param column The column name of the attribute that will be considered for blocking
#' @param label Method for textual blocking; choose among levenshtein, cosine, jaccard, jarowinker, qgram
#' @param max_distance The maximal distance allowed in a block
#'
#' @return blocks A dataframe of pairs
#'
#' @export


textual.blocking <- function(data, column, method, max_distance) {
  blocks<-NULL
  data$row <- rownames(data)

  if(is.null(column)) {
    stop('select the column for textual blocking')
  }

  if(is.na(match (column, names(data)))){
    stop('enter a valid column name')
  }

  if(is.null(method)) {
    stop('select a method: levenshtein, cosine, jaccard, jarowinker, qgram')
  }
  if(is.null(max_distance)) {
    stop('set max_distance to the maximal distance allowed for the selected column')
  }

  switch(method,
         levenshtein={
           blocks<-stringdist_inner_join(data, data,
                                         by=column, method = 'lv',
                                         max_dist = max_distance,
                                         ignore_case = TRUE)
         },
         cosine={
           blocks<-stringdist_inner_join(data, data,
                                         by=column, method = 'cosine',
                                         max_dist = max_distance,
                                         ignore_case = TRUE)
         },
         jaccard={
           blocks<-stringdist_inner_join(data, data,
                                         by=column, method = 'jaccard',
                                         max_dist = max_distance,
                                         ignore_case = TRUE)
         },
         jarowinker={
           blocks<-stringdist_inner_join(data, data,
                                         by=column, method = 'jw',
                                         max_dist = max_distance,
                                         ignore_case = TRUE)
         },
         qgram={
           blocks<-stringdist_inner_join(data, data,
                                         by=column, method = 'qgram',
                                         max_dist = max_distance,
                                         ignore_case = TRUE)
         },
         {
           stop('choose a valid method from levenshtein, cosine, jaccard, jaro-winker, qgram')
         }
  )

  blocks<- subset(blocks, blocks$row.x!=blocks$row.y)
  blocks<- subset(blocks, select=-c(row.x,row.y))

  return(blocks)
}



#' @title Spatial blocking
#' @description Creates blocks of entities that near spatially. Returns the pairs
#'

#' @param data A dataframe of entities
#' @param longitude The column name that contains the longitudes
#' @param latitude The column name that contains the latitudes
#' @param max_distance The maximal distance in meters allowed in a block
#'
#' @return blocks A dataframe of pairs
#'
#' @export


spatial.blocking <- function(data, longitude, latitude, max_distance) {
  blocks<-NULL
  data$row <- rownames(data)

  if(is.null(longitude)  | is.null(longitude)) {
    stop('specify the columns of the longitude and the latitude')
  }

  if(is.na(match (longitude, names(data)))){
    stop('enter a valid column name for the longitude')
  }

  if(is.na(match (latitude, names(data)))){
    stop('enter a valid column name for the latitude')
  }

  if(is.null(max_distance)) {
    stop('set max_distance to the maximal distance in meters to create spatial
         blocks')
  }
  blocks<-geo_inner_join(x=data, y=data,
                         by=c(longitude, latitude), method="haversine",
                         max_dist = max_distance/1000,
                         unit = "km")
  blocks<- subset(blocks, blocks$row.x!=blocks$row.y)
  blocks<- subset(blocks, select=-c(row.x,row.y))

  return(blocks)
  }



#' @title Prefix blocking
#' @description Creates blocks of entities that have the same prefix. Returns the pairs
#'

#' @param data A dataframe of entities
#' @param column The column name of the attribute on which the prefix should be calculated
#' @param prefix_size The maximal number of characters for prefix blocking
#'
#' @return blocks A dataframe of pairs
#'
#' @export

prefix.blocking <- function(data, column, prefix_size) {
  if(is.null(column)) {
    stop('select the column for prefix blocking')
  }

  if(is.na(match (column, names(data)))){
    stop('enter a valid column name')
  }


  if(is.null(prefix_size)) {
    stop('choose the prefix size (number of characters)')
  }

  blocks<-NULL
  data$row <- rownames(data)
  data$prefix<-substring(data[[column]], 1, prefix_size)
  blocks<-merge(x=data, y=data, by="prefix")
  blocks<- subset(blocks, blocks$row.x!=blocks$row.y)
  blocks<- subset(blocks, select=-c(row.x,row.y))
  blocks <- subset(blocks, select=-c(prefix))
  return(blocks)
}


#' @title Suffix blocking
#' @description Creates blocks of entities that have the same suffix. Returns the pairs
#'

#' @param data A dataframe of entities
#' @param column The column name of the attribute on which the suffix should be calculated
#' @param suffix_size The maximal number of characters for suffix blocking
#'
#' @return blocks A dataframe of pairs
#'
#' @export


suffix_blocking <- function(data, column, suffix_size) {
  if(is.null(column)) {
    stop('select the column for suffix blocking')
  }
  if(is.na(match (column, names(data)))){
    stop('enter a valid column name')
  }


  if(is.null(suffix_size)) {
    stop('choose the suffix size (number of characters)')
  }

  blocks<-NULL
  data$row <- rownames(data)
  data$suffix<-str_sub(data[[column]], - suffix_size,
                       -1)

  blocks<-merge(x=data, y=data, by="suffix")
  blocks<- subset(blocks, blocks$row.x!=blocks$row.y)
  blocks<- subset(blocks, select=-c(row.x,row.y))
  blocks <- subset(blocks, select=-c(suffix))

  return(blocks)
}



#' @title Pairwise textual similarity
#' @description Compares the pairs pairwise regarding a textual attribute. Returns a vector of text similarity
#'

#' @param data A dataframe of pairs
#' @param method A method for the text similarity, choose among levenshtein, cosine, jaccard, jaro-winker
#' @param column1 The first column name of the attribute that will be compared
#' @param column2 The second column name of the attribute that will be compared
#'
#' @return sim A vector of text similarities
#'
#' @export

text.similarity <- function (data,method, column1, column2){
  if(is.null(column1) | is.null(column2)) {
    stop('select the columns to compare')
  }

  if(is.na(match (column1, names(data)))){
    stop('enter a valid column name for column1')
  }

  if(is.na(match (column2, names(data)))){
    stop('enter a valid column name for column2')
  }

  switch(method,
         levenshtein={

           return(1 - (stringdist(data[[column1]], data[[column2]],
                                  method = 'lv')/max(nchar(data[[column1]]),
                                                     nchar( data[[column2]]))))
         },
         cosine={
           return(1 - stringdist(data[[column1]], data[[column2]],
                                 method = 'cosine'))
         },
         jaccard={
           return(1 - stringdist(data[[column1]], data[[column2]],
                                 method = 'jaccard'))
         },
         jarowinker={
           return(1 - stringdist(data[[column1]], data[[column2]],
                                 method = 'jw'))
         },

         {
           stop('choose a valid method from levenshtein, cosine, jaccard, jaro-winker')
         }
  )
}



#' @title Pairwise spatial similarity
#' @description Compares the pairs pairwise regarding a spatial attribute. Returns a vector of spatial similarity
#'

#' @param data A dataframe of pairs
#' @param lat1 The column name of the latitude of the first entity
#' @param long1 The column name of the longitude of the first entity
#' @param lat2 The column name of the latitude of the second entity
#' @param long2 The column name of the longitude of the second entity
#' @param max_distance The maximal distance allowed
#'
#' @return sim A vector of spatial similarities
#'
#' @export

spatial.similarity <- function (data,lat1, long1, lat2, long2, max_distance){
  if(is.null(lat1) | is.null(lat2) | is.null(long1) | is.null(long2)) {
    stop('select the columns to compare')
  }

  if(is.na(match (lat1, names(data)))){
    stop('enter a valid column name for lat1')
  }

  if(is.na(match (long1, names(data)))){
    stop('enter a valid column name for long1')
  }

  if(is.na(match (lat2, names(data)))){
    stop('enter a valid column name for lat2')
  }

  if(is.na(match (long2, names(data)))){
    stop('enter a valid column name for long2')
  }

  distance<-distVincentyEllipsoid(data[,c(long1,lat1)], data[,c(long2,lat2)])

  distance <- 1-(distance/max_distance)

  distance <- ifelse(distance <0, 0, distance)

}

#' @title Pairwise semantic similarity
#' @description Compares the pairs pairwise regarding their semantic attributes. Returns a vector of semantic similarity
#'

#' @param data A dataframe of pairs
#' @param column1 The column name containing the semantics for the first entity
#' @param column2 The column name containing the semantics for the second entity
#' @param pythonpath The path where the python file is saved
#' @param method The method for the semantic similarity; choose between path (Path similarity), and wup (Wu&Palmer)

#'
#' @return sim A vector of semantic similarities
#'
#' @export

semantic.similarity <- function (data, column1, column2, pythonpath, method){
  if(is.null(column1) | is.null(column2)) {
    stop('select the columns to compare')
  }

  if(is.na(match (column1, names(data)))){
    stop('enter a valid column name for column1')
  }

  if(is.na(match (column2, names(data)))){
    stop('enter a valid column name for column2')
  }

  use_python(pythonpath, required = T)
  nltk <- import("nltk")

  switch(method,
         wup={
            wup_python<- function(x,y){
             source_python('pythonscript/wup.py')
             return (add(x, y))
           }
         },
         path={
           wup_python<- function(x,y){
             source_python('pythonscript/path.py')
             return (add(x, y))
           }
         },

         {
           stop('choose a valid method from wup and path')
         }


)

  wup_python_list <- function (text1, text2)

  {
    if (is.na(text1) || text1 == ''){
      return (0.0)
    }

    if (is.na(text2) || text2 == ''){
      return (0.0)
    }

    text1<- gsub('[[:punct:] ]+',' ',text1)
    text2<-gsub('[[:punct:] ]+',' ',text2)
    if(strcmpi(text1,text2)){
      return (1.0)
    }
    p1<-as.list(strsplit(text1, '\\s+')[[1]])
    p2<-as.list(strsplit(text2, '\\s+')[[1]])

    if(length(intersect(p1, p2))!=0){
      return (1.0)
    }
    else {
      S<-c(0.0)
      for (i in p1){
        for (j in p2){
          S<- c(S, wup_python(i, j))
        }

        return (max(S))
      }
    }
  }

  return(mapply(wup_python_list, stringi::stri_enc_toutf8(data[[column1]]),
                stringi::stri_enc_toutf8(data[[column2]])))

}




#' Dataset of automatically labeled pairs of spatial entities that are at most 50 m far from each other
#'
#' @docType data
#'
#' @usage data(pairs50)
#'
#'
#' @keywords datasets
#'
#' @references Isaj, Suela, Esteban Zimányi, and Torben Bach Pedersen. "Multi-Source Spatial Entity Linkage." Proceedings of the 16th International Symposium on Spatial and Temporal Databases. ACM, 2019.

#'
#' @examples
#' data(pairs50)
#'
#'
#' #' Dataset of manually labeled pairs of spatial entities
#'
#' @docType data
#'
#' @usage data(pairsManual)
#'
#'
#' @keywords datasets
#'
#' @references Isaj, Suela, Esteban Zimányi, and Torben Bach Pedersen. "Multi-Source Spatial Entity Linkage." Proceedings of the 16th International Symposium on Spatial and Temporal Databases. ACM, 2019.

#'
#' @examples
#' data(pairsManual)

