#Setting file.
#Only change what you need

#The method used to determine BAE1 threshold, enter either 0 for “kmeans” or 1 for “absolute”
#kmeans clustering will automatically determine data that is close to a certain 
#cluster center to be eliminated
#absolute clustering takes a absolute coverage threshold; any data whose coverage is
#lower than that value will be eliminated
BAE1.threshold.mode.0.or.1 <- 1
BAE1.threshold.absolute.0.to.400 <- 100
BAE1.threshold.kmeans.lower.cluster.center.0.to.400 <- 50

#The method used to determine BAE2 threshold, enter either 0 for “kmeans” or 1 for “absolute”
BAE2.threshold.mode.0.or.1 <- 1
BAE2.threshold.absolute.0.to.1 <- 0.95
BAE2.threshold.kmeans.lower.cluster.center.0.to.1 <- 0.95

#How many snps per bin in order for it to be count as reliable ?
bin.threshold.1.max_amplicons_per_bin <- 4

#P-value threshold (2-tailed z-test)
p.value.threshold <- 0.005

