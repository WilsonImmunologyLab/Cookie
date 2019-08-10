#' plotReduction
#'
#' Plot all samples on a 2D reduction, samples can be colored by factors.
#'
#' @param object Cookie object
#' @param reduction Reduction results used.
#' @param group.factor Color samples by which factor.
#' @param sampling Hightlight the selected sample in a speficic sampling.
#'
#' @importFrom ggplot2 ggplot geom_point aes
#'
#' @export
#'

plotReduction <- function(
  object,
  reduction = NULL,
  group.factor = NULL,
  sampling = NULL
) {
  if(is.null(reduction)){
    stop("Please provide a reduction name! can be umap, tsne, cmds")
  }
  if(!is.null(object)){
    if(!is.null(object@reduction[[reduction]])){
      data <- data.frame(dim1=object@reduction[[reduction]]@embedding[,1],dim2=object@reduction[[reduction]]@embedding[,2],id=rownames(object@raw.data))
      if(!is.null(group.factor)){
        vector <- object@raw.data[,group.factor]
        data <- cbind(data,vector)
        colnames(data)[length(colnames(data))] <- "group"
      }

      if(!is.null(sampling)){
        vector <- object@samplings[[sampling]]@sampling
        size <- rep(2, times=dim(vector)[1])
        levels(vector[,1]) <- c(levels(vector[,1]),"Unselected")
        size[is.na(vector[,1])] <- 1
        vector[is.na(vector[,1]),1] <- "Unselected"
        data <- cbind(data,vector)
        colnames(data)[length(colnames(data))] <- "selected"
        data <- cbind(data,size)
        colnames(data)[length(colnames(data))] <- "size"
      }

      if(!is.null(group.factor) && !is.null(sampling)) {
        p <- ggplot(data = data, aes(x=dim1,y=dim2)) + geom_point(aes(shape=selected, color=group, size=size))
      } else if (!is.null(group.factor)) {
        p <- ggplot(data = data, aes(x=dim1,y=dim2)) + geom_point(aes(color=group))
      } else if (!is.null(sampling)) {
        p <- ggplot(data = data, aes(x=dim1,y=dim2)) + geom_point(aes(shape=selected, size=size))
      } else {
        p <- ggplot(data = data, aes(x=dim1,y=dim2)) + geom_point()
      }
      return(p)
    } else {
      stop("The reduction name you provided does not exist!")
    }
  } else {
    stop("Please provide a Cookie object!")
  }
}




#' plotSizeTest
#'
#' Plot coverage on each factor for a size test to determine an appropriate sample size
#'
#' @param object Cookie object
#' @param type Figure type, could be "radar" or "line"
#' @param test.name The name for a size test
#'
#' @importFrom ggplot2 ggplot geom_point aes geom_line theme labs ylab xlab scale_x_continuous theme_bw
#' @importFrom reshape2 melt
#' @importFrom ggradar ggradar
#'
#' @export
#'

plotSizeTest <- function(
  object,
  type = "line",
  test.name = NULL
) {
  if(is.null(test.name)){
    stop("Please provide a size test name!")
  }
  if(!is.null(object)){
    if(!is.null(object@sample.size.test[[test.name]])){
      data <- object@sample.size.test[[test.name]]@coverage
      if(type == "radar") {
        p <- ggradar(data)
      } else {
        melt.data <- melt(data, id.vars = "Size")
        p <- ggplot(melt.data, aes(Size, value, group=variable, color=factor(variable))) + geom_point(color="black") +
          geom_line(linetype="solid", size=1) + theme(legend.position = "right") + labs(color = "Factors") + ylab("Coverage") + xlab("Sampling sizes") +
          scale_x_continuous(breaks=data$Size, labels=data$Size) + theme_bw()
      }
      return(p)
    } else {
      stop("The sample size test name you provided does not exist!")
    }
  } else {
    stop("Please provide a Cookie object!")
  }
}


#' plotSampling
#'
#' Plot coverage on each factor for a sampling
#'
#' @param object Cookie object
#' @param name The name for a size test
#'
#' @importFrom ggradar ggradar
#'
#' @export
#'

plotSampling <- function(
  object,
  name = NULL
) {
  if(is.null(name)){
    stop("Please provide a size test name!")
  }
  if(!is.null(object)){
    if(!is.null(object@samplings[[name]])){
      data <- object@samplings[[name]]@coverage
      p <- ggradar(data)
      return(p)
    } else {
      stop("The sample size test name you provided does not exist!")
    }
  } else {
    stop("Please provide a Cookie object!")
  }
}



#' plotDistribution
#'
#' Plot coverage on each factor for a sampling
#'
#' @param object Cookie object
#' @param sampling choose original data or a specific sampling subset for analysis
#' @param factor The factor for plot
#'
#' @importFrom ggplot2 ggplot aes geom_rect coord_polar theme_bw labs theme element_blank geom_text
#'
#' @export
#'

plotDistribution <- function(
  object,
  sampling = "original",
  factor = NULL
) {
  if(is.null(factor)){
    stop("Please provide a factor name!")
  }
  if(!is.null(object)){
    if(sampling == "original") {
      data <- object@raw.data
    } else if(!is.null(object@samplings[[sampling]])){
      data <- object@raw.data
      sel <- object@samplings[[sampling]]@sampling
      sel.index <- !is.na(sel[,1])
      data <- data[sel.index,]
    } else {
      stop("The sampling name you provided does not exist!")
    }
    if(factor %in% colnames(data)){
      data <- data[,factor]

      plot.data <- as.data.frame(table(data))
      plot.data$fraction = plot.data$Freq / sum(plot.data$Freq)
      plot.data$ymax = cumsum(plot.data$fraction)
      plot.data$ymin = c(0, head(plot.data$ymax, n = -1))
      plot.data$label <- plot.data$data

      for (n in 1:dim(plot.data)[1]) {
        cur <- plot.data[n,]
        if(cur$fraction < 0.02) {
          plot.data$label[n] <- " "
        }
      }

      p <- ggplot(data = plot.data, aes(fill = label, ymax = ymax, ymin = ymin,xmax = 4, xmin = 3)) +
        geom_rect(colour = "grey30") +
        coord_polar(theta = "y") +
        theme_bw() +
        labs(x = "", y = "", title = factor) +
        theme(panel.grid=element_blank()) +
        theme(axis.text=element_blank()) +
        theme(axis.ticks=element_blank()) +
        theme(panel.border=element_blank()) +
        geom_text(aes(x = 3.9, y = ((ymin+ymax)/2), label = label)) + theme(legend.position = "none")

      return(p)
    } else {
      stop("The factor name you provided does not exist!")
    }
  } else {
    stop("Please provide a Cookie object!")
  }
}
