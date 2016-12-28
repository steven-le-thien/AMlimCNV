#A function that is used to split the coverages into their corresponding PCR pool
pool_splitting <- function(coverage) {
  # pool_cluster <- apply(coverage[, 1:2],1, function(x) return(as.integer(readline(prompt = paste("Enter the pool of", x[1], "position", x[2], ": ")))))
  pool_cluster <- apply(coverage[, 1:2],1, function(x) return(1))
  coverage <- cbind(coverage, pool_cluster)
  return(lapply(unique(pool_cluster), function(x) return(coverage[coverage$pool_cluster == x,])))
}

#A function that returns ratio matrix with rij = amplicon_i/amplicon_j
ratio_matrix <- function(single_coverage) {
  output <- lapply(unlist(single_coverage), function(x) return(unlist(single_coverage)/x))
  output <- do.call(cbind, output)
  return (as.data.frame(output))
}

#Performing am-lim correction to estimate true ratio matrix
am_lim <- function(reference_coverage, test_coverage) {
  #Preparing the matrices
  reference_matrix <- apply(reference_coverage, 2, ratio_matrix)
  test_matrix <- ratio_matrix(test_coverage)
  target_matrix <- as.data.frame(test_matrix/(Reduce("+", reference_matrix)/length(reference_matrix)))

  #Perform AM-lim correction
  while (abs(target_matrix[1,2] - 1/target_matrix[2,1]) > 0.000005)
    target_matrix <- (target_matrix + 1/t(target_matrix))/2
  return(as.data.frame(target_matrix))
}

#Main function
ratio_analysis <- function(cov) {
  #Splitting into pools
  pools <- pool_splitting(cov)

  #Perform am_lim on each pool
  am_lim_matrices <- lapply(pools, function(x) return(am_lim(as.data.frame(x)[,3:(dim(as.data.frame(x))[2] - 2)], as.data.frame(x)[,dim(as.data.frame(x))[2]-1])))

  #Merge pools
  current_pool <- cbind(as.data.frame(pools[[1]])[, 1:2], as.data.frame(am_lim_matrices[[1]]))

  #Take only the first column of the matrix
  current_pool <- current_pool[,c(1,3)]

  #Take combined result in each bin and determine whether each bin has enough data
  unbinned <- tapply(as.numeric(current_pool[,2]), current_pool[,1], function(x) return (c(mean(x), length(x) > bin.threshold.1.max_amplicons_per_bin)))
  binned <- do.call(rbind, unbinned)
  rownames(binned) <- names(unbinned)
  rownames(binned) <- as.numeric(lapply(rownames(binned), chrom_position))
  colnames(binned) <- c('Ratio', 'Confident')
  binned <- as.data.frame(binned[order(as.integer(rownames(binned))),])
  binned <- binned[binned$Confident == 1,]
  output <- binned[,1]
  names(output) <- rownames(binned)

  return(output)
}




