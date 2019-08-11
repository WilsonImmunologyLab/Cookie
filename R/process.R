#' normalization
#'
#' normalize data for numerical factors
#'
#' @param object Cookie object
#'
#' @export
#'
#'

normalization <- function (
  object = NULL
) {
  if(!is.null(object)) {
    types <- object@factor.type
    num.index <- which(types == "num")
    raw.data <- object@raw.data
    if(length(num.index) > 0) {
      for (index in num.index) {
        col <- raw.data[,index]
        na.index <- which(col == "NA")
        col[na.index] <- min(col)
        raw.data[,index] <- (col - min(col))/(max(col) - min(col))
        raw.data[na.index,index] <- -1
      }
    } else {
      cat("Didn't find any numerical factor, will directly copy data from raw.data slot!\n")
    }
    object@normalize.data <- raw.data
    return(object)
  } else {
    stop("Please provide Cookie object")
  }
}


#' distCalculation
#'
#' calculate distance from normalized data
#'
#' @param object Cookie object
#'
#' @export
#'
#'

distCalculation <- function (
  object = NULL
) {
  if(!is.null(object)) {
    types <- object@factor.type
    num.index <- which(types == "num")
    data <- object@normalize.data
    if(length(num.index) > 0) {
      # numerical matrix
      data1 <- data[,num.index]
      # bool and char matrix
      data2 <- data[,-num.index]

      dist.matrix1 <- hammingCodingCpp(data1)
      dist.matrix2 <- binaryCodingCpp(data2)

      dist.matrix <- dist.matrix1 + dist.matrix2
    } else {
      dist.matrix <- binaryCodingCpp(data)
    }

    object@dist.matrix <- dist.matrix

    return(object)
  } else {
    stop("Please provide Cookie object")
  }
}


#' binaryCoding
#'
#' calculate binary distance for bool and char distance
#'
#' @param data data matrix
#'
#' @export
#'
#'

binaryCoding <- function(
  data = NULL
) {
  n <- dim(data)[1]
  m <- dim(data)[2]
  res <- matrix(0,n,n)
  for (i in 1:n) {
    cat("i = ", i, " \n")
    for (j in i:n) {
      #cat("j = ", j, " \n")
      vec.i <- data[i,]
      vec.j <- data[j,]

      vec.res <- rep(0,m)
      for (k in 1:m) {
        if(vec.i[k] == vec.j[k]) {
          vec.res[k] <- 0
        } else {
          vec.res[k] <- 1
        }
      }
      a <- sum(vec.res)
      res[i,j] <- a
      res[j,i] <- a
    }
  }
  return(res)
}

#' hammingCoding
#'
#' calculate hamming distance for bool and char distance
#'
#' @param data data matrix
#'
#' @export
#'
#'

hammingCoding <- function(
  data = NULL
){
  n <- dim(data)[1]
  m <- dim(data)[2]
  res <- matrix(0,n,n)
  for (i in 1:n) {
    cat("i = ", i, " \n")
    for (j in i:n) {
      #cat("j = ", j, " \n")
      vec.i <- data[i,]
      vec.j <- data[j,]

      vec.res <- abs(vec.i - vec.j)
      a <- sum(vec.res)
      res[i,j] <- a
      res[j,i] <- a
    }
  }
  return(res)
}


#' reductionTSNE
#'
#' Run t-SNE to project samples into 2D map using pairwise distances
#'
#' @param object (For Seurat) Seurat object
#' @param assay (For Seurat) run t-SNE for which assay, choose from RNA, ADT, Joint or All
#' @param perplexity numeric; Perplexity parameter (should not be bigger than 3 * perplexity < nrow(X) - 1, see details for interpretation)
#' @param dim integer; Output dimensionality (default: 2)
#' @param seed integer; seed for reproducible results.
#' @param theta numeric; Speed/accuracy trade-off (increase for less accuracy), set to 0.0 for exact TSNE (default: 0.5)
#'
#' @importFrom Rtsne Rtsne
#'
#' @export
#'

reductionTSNE <- function(
  object,
  perplexity = 30,
  dim = 2,
  seed = 42,
  theta = 0.5
) {
  if(!is.null(object)){
    set.seed(seed = seed)

    cat("Start run t-SNE from distances...\n")

    data <- object@dist.matrix

    # run tsne with pairwise distances
    my.tsne <- Rtsne(data, perplexity = perplexity, is_distance = TRUE, dims = dim, theta = theta)

    object@reduction[['tsne']] <- createDimReductionObject(data = my.tsne$Y, method = "t-SNE")
    return(object)
  } else {
    stop("Please provide a Cookie object!")
  }
}




#' reductionUMAP
#'
#' Run UMAP to project samples into 2D space using pairwise distances
#'
#' @param object Cookie object
#' @param seed see number. default is 42
#' @param method could be "naive" or "umap-learn"
#' @param n.neighbors integer; number of nearest neighbors
#' @param n.components  integer; dimension of target (output) space
#' @param metric character or function; determines how distances between data points are computed. When using a string, available metrics are: euclidean, manhattan. Other available generalized metrics are: cosine, pearson, pearson2. Note the triangle inequality may not be satisfied by some generalized metrics, hence knn search may not be optimal. When using metric.function as a function, the signature must be function(matrix, origin, target) and should compute a distance between the origin column and the target columns
#' @param verbose logical or integer; determines whether to show progress messages
#' @param n.epochs  integer; number of iterations performed during layout optimization
#' @param min.dist  numeric; determines how close points appear in the final layout
#' @param spread numeric; used during automatic estimation of a/b parameters.
#' @param set.op.mix.ratio numeric in range [0,1]; determines who the knn-graph is used to create a fuzzy simplicial graph
#' @param local.connectivity  numeric; used during construction of fuzzy simplicial set
#' @param negative.sample.rate  integer; determines how many non-neighbor points are used per point and per iteration during layout optimization
#'
#' @importFrom umap umap umap.defaults
#'
#' @export
#'
reductionUMAP <- function(
  object,
  seed = 42,
  method = "umap-learn",
  n.neighbors = 15,
  n.components = 2,
  metric = "euclidean",
  verbose = TRUE,
  n.epochs = 200,
  min.dist = 0.1,
  spread = 1,
  set.op.mix.ratio = 1,
  local.connectivity = 1L,
  negative.sample.rate = 5L
) {
  if(!is.null(object)){
    set.seed(seed = seed)

    my.umap.conf <- umap.defaults
    my.umap.conf$input <- "dist"
    my.umap.conf$n_neighbors <- n.neighbors
    my.umap.conf$n_components <- n.components
    my.umap.conf$metric <- metric
    my.umap.conf$verbose <- verbose
    my.umap.conf$n_epochs <- n.epochs
    my.umap.conf$min_dist <- min.dist
    my.umap.conf$spread <- spread
    my.umap.conf$set_op_mix_ratio <- set.op.mix.ratio
    my.umap.conf$local_connectivity <- local.connectivity
    my.umap.conf$negative_sample_rate <- negative.sample.rate

    cat("Start run UMAP from distances...\n")

    dist <- object@dist.matrix
    my.umap <- umap(dist,my.umap.conf,method = method)

    object@reduction[['umap']] <- createDimReductionObject(data = my.umap$layout, method = "UMAP")
    return(object)
  } else {
    stop("Please provide a Cookie object!")
  }
}


#' reductionMDS
#'
#' Run Classical MDS to project samples into 2D space using pairwise distances
#'
#' @param object Cookie object
#' @param eig Indicates whether eigenvalues should be returned.
#' @param k The maximum dimension of the space which the data are to be represented in; must be in {1, 2, …, n-1}.
#' @param add Logical indicating if an additive constant c* should be computed, and added to the non-diagonal dissimilarities such that the modified dissimilarities are Euclidean.
#' @param x.ret Indicates whether the doubly centred symmetric distance matrix should be returned.
#'
#' @export
#'

reductionMDS <- function(
  object,
  eig = FALSE,
  k = 2,
  add = FALSE,
  x.ret = FALSE
) {
  if(!is.null(object)){
    cat("Start run Classical MDS from distances...\n")

    data <- object@dist.matrix

    # run tsne with pairwise distances
    my.mds <- cmdscale(data,eig=eig, k=k,add=add, x.ret = x.ret)
    my.mds <- as.matrix(my.mds)
    object@reduction[['cmds']] <- createDimReductionObject(data = my.mds, method = "Classical MDS")
    return(object)
  } else {
    stop("Please provide a Cookie object!")
  }
}



#' sampleSizeTest
#'
#' Run sample size test for current dataset
#'
#' @param object Cookie object
#' @param prime.factor The unique prime factor.
#' @param size.range Sample size range
#' @param name A name for this run. e.g. test1
#'
#' @importFrom cluster pam
#'
#' @export
#'

sampleSizeTest <- function(
  object = NULL,
  prime.factor = NULL,
  size.range = NULL,
  name = NULL
) {
  if(!is.null(object)){
    if(!is.null(name)) {
      n.sample <- dim(object@normalize.data)[1]
      n.factor <- dim(object@normalize.data)[2]
      n.size <- length(size.range)
      dist.matrix <- object@dist.matrix
      data <- object@normalize.data
      type <- object@factor.type

      coverage <- matrix(data = NA, nrow = n.size, ncol = (n.factor + 1))
      selection <- matrix(data = NA, nrow = n.sample, ncol = n.size)

      if(!is.null(prime.factor)) {
        # sampling from each level of prime factor
        factors <- colnames(object@normalize.data)
        if(prime.factor %in% factors) {
          subject.list <- unique(data[,prime.factor])
          i = 1
          for (n in size.range) {
            cat("test sample size = ",n,"for each subject in prime factor... \n")
            for (subject in subject.list) {
              index <- which(data[,prime.factor] == subject)
              if(length(index) > n) {
                sub.dist.matrix <- dist.matrix[index,index]
                res <- pam(x = sub.dist.matrix, k = n, diss = TRUE)

                orig.index <- index[res[["id.med"]]]
                selection[orig.index,i] <- "Selected"
              } else {
                cat("Number of samples in current subject = ",subject," is less than n = ",n,", will select all samples in this subject!\n")
                selection[index,i] <- "Selected"
              }
            }

            a <- selection[,i]
            index.a <- which(!is.na(a))
            subset <- data[index.a,]
            coverage[i,1] = n
            for (j in 1:n.factor) {
              if(type[j] != "num") {
                coverage[i,(j+1)] <- length(unique(subset[,j]))/length(unique(data[,j]))
              } else {
                coverage[i,(j+1)] <- length(unique(floor(subset[,j]*10)))/length(unique(floor(data[,j]*10)))
              }
            }
            i <- i + 1
          }
        } else {
          stop("The prime factor you provided is not exist! Please check your input!")
        }
      } else {
        # sampling from the entire population
        i = 1
        for (n in size.range) {
          cat("test sample size = ",n," for the entire population... \n")
          res <- pam(x = dist.matrix, k = n, diss = TRUE)
          selection[res[["id.med"]],i] <- "Selected"

          subset <- data[res[["id.med"]],]
          for (j in 1:n.factor) {
            if(type[j] != "num") {
              coverage[i,j] <- length(unique(subset[,j]))/length(unique(data[,j]))
            } else {
              coverage[i,j] <- length(unique(floor(subset[,j]*10)))/length(unique(floor(data[,j]*10)))
            }
          }
          i <- i + 1
        }
      }
      coverage <- as.data.frame(coverage)
      rownames(coverage) <- size.range
      colnames(coverage) <- c("Size",colnames(data))

      selection <- as.data.frame(selection)
      rownames(selection) <- rownames(data)
      colnames(selection) <- size.range

      object@sample.size.test[[name]] <- createSampleSizeTestObject(prime.factor = prime.factor,coverage = coverage,selection = selection)
      return(object)
    } else {
      stop("Please provide a name for this test (e.g. test1)!")
    }
  } else {
    stop("Please provide a Cookie object!")
  }
}



#' sampling
#'
#' Sampling from current dataset
#'
#' @param object Cookie object
#' @param prime.factor The unique prime factor.
#' @param important.factor the important factors
#' @param sample.size Sample size
#' @param name A name for this run. e.g. test1
#'
#' @importFrom cluster pam
#'
#' @export
#'

sampling <- function(
  object = NULL,
  prime.factor = NULL,
  important.factor = NULL,
  sample.size = NULL,
  name = NULL
) {
  if(!is.null(object)){
    if(!is.null(name)) {
      n.sample <- dim(object@normalize.data)[1]
      n.factor <- dim(object@normalize.data)[2]
      dist.matrix <- object@dist.matrix
      data <- object@normalize.data
      type <- object@factor.type

      coverage <- matrix(data = NA, nrow = 1, ncol = (n.factor + 1))
      selection <- matrix(data = NA, nrow = n.sample, ncol = 1)
      # step 1
      if(!is.null(prime.factor)) {
        # sampling from each level of prime factor
        factors <- colnames(object@normalize.data)
        if(prime.factor %in% factors) {
          subject.list <- unique(data[,prime.factor])

          for (subject in subject.list) {
            index <- which(data[,prime.factor] == subject)
            if(length(index) > sample.size) {
              sub.dist.matrix <- dist.matrix[index,index]
              res <- pam(x = sub.dist.matrix, k = sample.size, diss = TRUE)

              orig.index <- index[res[["id.med"]]]
              selection[orig.index,1] <- "Selected"
            } else {
              cat("Number of samples in current subject = ",subject," is less than n = ",sample.size,", will select all samples in this subject!\n")
              selection[index,1] <- "Selected"
            }
          }
        } else {
          stop("The prime factor you provided is not exist! Please check your input!")
        }
      } else {
        # sampling from the entire population
        res <- pam(x = dist.matrix, k = sample.size, diss = TRUE)
        selection[res[["id.med"]],1] <- "Selected"
      }

      # step 2
      index.a <- which(!is.na(selection))
      subset <- data[index.a,]
      if(!is.na(important.factor)) {
        for (important in important.factor) {
          original.important <- unique(data[,important])
          cur.important <- unique(subset[,important])
          diff <- setdiff(original.important, cur.important)

          if(length(diff) > 0) {
            for (var in diff) {
              candidate.index <- which(data[,important] == var)
              a <- dist.matrix[candidate.index, index.a]
              a <- rowSums(a)
              sel.index <- which(a == min(a))
              sel.index <- candidate.index[sel.index]

              # set selected marker
              selection[sel.index] <- "Selected"
            }
          }
        }
      }

      # summary
      a <- selection[,1]
      index.a <- which(!is.na(a))
      subset <- data[index.a,]
      coverage[1,1] <- sample.size
      for (j in 1:n.factor) {
        if(type[j] != "num") {
          coverage[1,(j+1)] <- length(unique(subset[,j]))/length(unique(data[,j]))
        } else {
          coverage[1,(j+1)] <- length(unique(floor(subset[,j]*10)))/length(unique(floor(data[,j]*10)))
        }
      }

      coverage <- as.data.frame(coverage)
      rownames(coverage) <- sample.size
      colnames(coverage) <- c("Size",colnames(data))

      selection <- as.data.frame(selection)
      rownames(selection) <- rownames(data)
      colnames(selection) <- sample.size

      object@samplings[[name]] <- createSamplingObject(prime.factor = prime.factor,important.factor = important.factor, coverage = coverage, sampling = selection, size = sample.size)
      return(object)
    } else {
      stop("Please provide a name for this test (e.g. test1)!")
    }
  } else {
    stop("Please provide a Cookie object!")
  }
}
