#BAE2 or Bad Amplicon Elimination 2 is the process in which inconsistent SNP amplicons are removed.
#BAE2. The method used is polynomial regression over a Coverage data that has been normalized with
#Quality

#A function used to generate relative quality, which is the ratio of any quality to the greatest
#quality for that amplicon
relative_quality <- function(quality) {
  max_quality_vector <- apply(quality, 1, max)
  return(quality/max_quality_vector)
}

#A function that returns the formula used to normalize coverage against quality (see paper)
normalize_definition <- function(raw_reference, normalized_target, relative_quality) {
  if (normalized_target > raw_reference)
    return (normalized_target - abs(normalized_target - raw_reference) * relative_quality)
  else
    return (normalized_target + abs(normalized_target - raw_reference) * relative_quality)
}

#A function that returns does the normalizing of the coverages
normalizing_reference <- function(reference_coverage, normalized_target, relative_quality) {
  num_of_reference <- dim(reference_coverage)[2]
  num_of_amplicons <- dim(reference_coverage)[1]
  output <- data.frame(matrix(nrow = num_of_amplicons, ncol = num_of_reference))

  #Loop through all the amplicons and references data
  for (amplicon in 1:num_of_amplicons)
    for (reference in 1:num_of_reference)
      output[amplicon, reference] <- normalize_definition(reference_coverage[amplicon, reference],
                                                          normalized_target[amplicon]/num_of_reference,
                                                          relative_quality[amplicon, reference])
  colnames(output) <- colnames(reference_coverage)
  return(output)
}

#A function that generate the formula for multivariate polynoial regression
multivariate_poly_reg_formula_generator <- function(reference) {
  num_of_reference <- dim(reference)[2]
  stub <- paste("bae2_target_coverage ~ 0 + normalized_reference$", colnames(reference)[1], sep = "")
  for (counter in 2:num_of_reference)
    stub <- paste(stub, " + normalized_reference$", colnames(reference)[counter], sep = "")
  return(stub)
}

#A function that generates expected coverage using coefficients from regression
expecting_coverage <- function(normalized_reference, coefficient) {
  return(apply(t(normalized_reference) * unlist(coefficient), 2, sum))
}

#A function that generates the harmonic means cost from actual and expected coverage of each amplicon
harmonic_mean_cost <- function(expected, actual) {
  forward_ratio <- expected/actual
  backward_ratio <- actual/expected
  return (2/(forward_ratio + backward_ratio))
}

#A function that returns the cost for a single BAE2 run
bae2_report <- function(i, num_cov, num_qual) {
  print(paste(round(i/dim(num_cov)[2]*100, 2), "%", " ...", sep = ""))

  #Housekeeping
  bae2_target_coverage <- num_cov[,i]
  bae2_target_quality <- num_qual[,i]
  bae2_reference_coverage <- num_cov[, -i]
  bae2_reference_quality <- num_qual[, -i]

  normalized_quality <- relative_quality(bae2_reference_quality)
  normalized_reference <- normalizing_reference(bae2_reference_coverage, bae2_target_coverage, normalized_quality)

  #Generate regression formula
  reg_formula <- as.formula(multivariate_poly_reg_formula_generator(normalized_reference))

  #Does the regression
  reg_result <- lm(formula = reg_formula)


  #Obtain coefficient
  reg_coefficient <- reg_result['coefficients']

  #Generate cost
  expected <- expecting_coverage(normalized_reference, reg_coefficient)
  cost <- harmonic_mean_cost(expected, bae2_target_coverage)
  names(cost) <- names(bae2_target_coverage)
  if (is_test < 0) {
    plot(cost, xlab = "Amplicon", ylab = "Ind. Cost", main = "Bad Amplicon Elimination 2 individual Report", sub = colnames(num_cov)[i],ylim = c(0,1) )
  }
  return (cost)
}

#Removing amplicons based on kmeans clustering
bae2_kmeans <- function(coverage, mean_vector) {
  kmeans_result <- tryCatch({
    kmeans_report <- kmeans(mean_vector, centers = c(0.98,1))
    clustered_coverage <- cbind(coverage, kmeans_report["cluster"])
    filtered_coverage <- clustered_coverage[clustered_coverage$cluster == 2,]
    filtered_coverage <- filtered_coverage[ ,1:dim(filtered_coverage)[2] - 1]
    return (filtered_coverage)
  }, warning = function(war) {

  }, error = function(err) {
    return (-1)
  })

  if (kmeans_result == -1) {
    print("No amplicons eliminated through BAE2")
    return (coverage)
  } else {
    return (kmeans_result)
  }
}

#Removing amplicons based on absolute cut-off
bae2_absolute <- function(coverage, mean_vector) {
  clustered_coverage <- cbind(coverage, mean_vector)
  filtered_coverage <- clustered_coverage[clustered_coverage$mean_vector > BAE2.threshold.absolute.0.to.1,]
  filtered_coverage <- filtered_coverage[ ,1:dim(filtered_coverage)[2] - 1]

  return (filtered_coverage)
}

#Main function
bae2 <- function(coverage, quality) {
  #Housekeeping
  numeric_coverage <- coverage[, 3:dim(coverage)[2]]
  numeric_quality <- quality[, 3:dim(coverage)[2]]

  #Taking 1 reference as target at a time, rotate to ensure no bias in target coverage
  output <- data.frame(matrix(nrow = dim(numeric_coverage)[2], ncol = dim(numeric_coverage)[2]))
  report <- lapply(1:dim(numeric_coverage)[2], bae2_report, numeric_coverage, numeric_quality)


  #Uncomment this line to generate a report file
  # write.csv2(as.data.frame(report), "report_report.csv", row.names = FALSE)
  mean_vector <- apply(as.data.frame(report), 1, mean)

  #Plot result
  if (is_test < 0) {
    plot(mean_vector, xlab = "Amplicon", ylab = "Mean Cost", main = "Bad Amplicon Elimination 2 combined Report", ylim = c(0,1))
  }

  #Eliminate amplicon
  if (BAE2.threshold.mode.0.or.1 == 0)
    return (bae2_kmeans(coverage, mean_vector))
  else if (BAE2.threshold.mode.0.or.1 == 1)
    return (bae2_absolute(coverage, mean_vector))
}
