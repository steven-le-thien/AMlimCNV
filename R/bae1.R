#BAE1 or Bad Amplicon Elimination 1 is a process in which SNP amplicons with low average coverage 
#read counts are being identified and disposed off. The method used during SRP2014-2015 was to
#identify manually those amplicons whose average Original Coverage or Coverage is clearly far 
#lower than the rest. 

#This script offers an automatic way of doing so through 2-means clustering with one cluster with
#means 0 and another with means = the average of the top half of the sorted data-set
#(these means are heuristics and does not need to be so accurate). This is done for both Coverage 
#and Original Average index

bae1_kmeans <- function(coverage, mean_array) {
  #Calculate heuristics for true coverage average (higher cluster mean)
  true_mean_heuristics <- mean(sort(mean_array, method = 'radix', decreasing = TRUE)
                               [length(mean_array)/2:length(mean_array)])
  
  #Assigning cluster through 2-means clustering 
  kmeans_result <- kmeans(mean_array, centers =  c(BAE1.threshold.kmeans.lower.cluster.center.0.to.400
                                                   , true_mean_heuristics))
  
  #Remove amplicons that are in cluster 1 (cluster around 0)
  clustered_coverage <- cbind(coverage, kmeans_result["cluster"])
  filtered_coverage <- clustered_coverage[clustered_coverage$cluster == 2,]
  filtered_coverage <- filtered_coverage[ ,1:dim(filtered_coverage)[2] - 1]
  
  return (filtered_coverage)
}

bae1_absolute <- function (coverage, mean_array) {
  clustered_coverage <- cbind(coverage, mean_array)
  filtered_coverage <- clustered_coverage[clustered_coverage$mean_array > BAE1.threshold.absolute.0.to.400,]
  filtered_coverage <- filtered_coverage[ ,1:dim(filtered_coverage)[2] - 1]
  
  return (filtered_coverage)
}

bae1 <- function (coverage) {
  #Calculate the average coverage across all samples
  numeric_array <- coverage[, 3:dim(coverage)[2]]
  mean_array <- apply(numeric_array, 1, mean)
  
  if (is_test < 0) {
    plot(mean_array, xlab = "Amplicon", ylab = "Mean Coverage", main = "Bad Amplicon Elimination 1 combined Report", ylim = c(0,400))
  }
  
  if (BAE1.threshold.mode.0.or.1 == 0)
    return (bae1_kmeans(coverage, mean_array))
  else if (BAE1.threshold.mode.0.or.1 == 1)
    return (bae1_absolute(coverage, mean_array))
}



