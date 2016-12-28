pool_splitting <- function(coverage) {
  # pool_cluster <- apply(coverage[, 1:2],1, function(x) return(as.integer(readline(prompt = paste("Enter the pool of", x[1], "position", x[2], ": ")))))
  pool_cluster <- apply(coverage[, 1:2],1, function(x) return(1))
  coverage <- cbind(coverage, pool_cluster)
  return(lapply(unique(pool_cluster), function(x) return(coverage[coverage$pool_cluster == x,])))
}

ratio_matrix <- function(single_coverage) {
  output <- lapply(unlist(single_coverage), function(x) return(unlist(single_coverage)/x))
  output <- do.call(cbind, output)
  return (as.data.frame(output))
}

am_lim <- function(reference_coverage, test_coverage) {
  reference_matrix <- apply(reference_coverage, 2, ratio_matrix)
  test_matrix <- ratio_matrix(test_coverage)
  target_matrix <- as.data.frame(test_matrix/(Reduce("+", reference_matrix)/length(reference_matrix)))
  while (abs(target_matrix[1,2] - 1/target_matrix[2,1]) > 0.000005)
    target_matrix <- (target_matrix + 1/t(target_matrix))/2
  return(as.data.frame(target_matrix))
}

merging_pools <- function(pools_results) {
  
}

