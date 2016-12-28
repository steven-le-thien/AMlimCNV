enough_data_in_bin <- function(x) {
  if (x[2] == 1)
    return (x[1])
  else 
    return (NA)
}
function2 <- function(x) {
  bin_no <- substr(x, 4, nchar(x))
  if (bin_no == 'X')
    bin_no <- 24
  if (bin_no == 'Y')
    bin_no <- 25
  return (bin_no)
}
single_test <- function (c, reference, test, ref_quality, test_quality) {
  is_test <<- c
  if (c > 0) {
    input_test <- cbind(reference[, 1:2], reference[,2 + c])
    colnames(input_test) <- c("Chrom", "Position", "numeric_array")
    input_test_quality <- cbind(ref_quality[, 1:2], ref_quality[,2 + c])
    input_reference <- reference[,-(c+2)]
    input_ref_quality <- ref_quality[, -(c+2)]
  } else {
    input_test <- test
    input_test_quality <- test_quality
    input_reference <- reference
    input_ref_quality <- ref_quality
  }

  #Bad Amplicon Elimination. Noting that in principle, SNP amplicons that is identified as 'bad amplicons'
  #in reference samples would also be 'bad' in test samples. Thus, at this stage, we are pooling both
  #Run BAE1
  print("Performing BAE1")
  setwd(script.dir)
  source("bae1.R")
  
  pooled_coverage <- merge(input_reference, input_test, by = c("Chrom", "Position"))
  pooled_quality <- merge(input_ref_quality, input_test, input_test_quality = c("Chrom", "Position"))
  
  coverage_bae1 <- bae1(pooled_coverage)
  quality_bae1 <-merge(pooled_quality, coverage_bae1[1:2], by = c("Chrom", "Position"))
  
  print("...done")
  
  #Run BAE2
  print("Performing BAE2")
  setwd(script.dir)
  source("bae2.R")
  
  coverage_bae1_bae2 <- bae2(coverage_bae1, quality_bae1)
  print("...done")
  
  #Run ratio analysis
  print("Performing ratio analysis")
  setwd(script.dir)
  source("ratio_analysis.R")
  
  pools <<- pool_splitting(coverage_bae1_bae2)
  am_lim_matrices <- lapply(pools, function(x) return(am_lim(as.data.frame(x)[,3:(dim(as.data.frame(x))[2] - 2)], as.data.frame(x)[,dim(as.data.frame(x))[2]-1])))
    
  current_pool <- cbind(as.data.frame(pools[[1]])[, 1:2], as.data.frame(am_lim_matrices[[1]]))
  current_pool <- current_pool[,c(1,3)]
  unbinned <- tapply(as.numeric(current_pool[,2]), current_pool[,1], function(x) return (c(mean(x), length(x) > bin.threshold.1.max_amplicons_per_bin)))
  binned <- do.call(rbind, unbinned)
  rownames(binned) <- names(unbinned)
  rownames(binned) <- as.numeric(lapply(rownames(binned), function2 ))
  colnames(binned) <- c('Ratio', 'Confident')
  # binned <- as.numeric(lapply(binned, enough_data_in_bin))
  binned <- as.data.frame(binned[order(as.integer(rownames(binned))),])
  binned <- binned[binned$Confident == 1,]
  output <- binned[,1]
  names(output) <- rownames(binned)
  print("...done")
  
  return(output)
}
report_cnv <- function (order, p_value, z_score, test_report, pop_mean, pop_sd, names) {
  
  if (!is.na(p_value[order])) {
    if (p_value[order] < p.value.threshold) {
      print(paste("Reject disomic hypothesis, bin", names[order], "does not have copy number 2"))
      if (z_score[order] > 0) {
        print ("Moving on to trisomic hypothesis testing")
        tester <- test_report[order] / 1.5
        if (z_test(tester, pop_mean, pop_sd) < p.value.threshold) {
          print(paste("Reject trisomic hypothesis, bin", names[order], "does not have copy number 3"))
          print("Test is inconclusive")
        } else {
          print(paste("Bin", names[order], "has copy number 3"))
        }
      }
      else {
        print ("Moving on to monosomic hypothesis testing")
        tester <- test_report[order] / 0.5
        if (z_test(tester, pop_mean, pop_sd) < p.value.threshold) {
          print(paste("Reject monosomic hypothesis, bin", names[order], "does not have copy number 1"))
          print("Test is inconclusive")
        } else {
          print(paste("Bin", names[order], "has copy number 1"))
        }
      }
    } else {
      print(paste("Bin", names[order], "has copy number 2"))
    }
  }
}
z_test <- function(testand, pop_mean, pop_sd) {
  return (2*pnorm(-abs((testand-pop_mean)/pop_sd)))
}

run_full_test <- function(reference, test, ref_quality, test_quality) {
  print("Investigating population...")
  normal_sampling_report <- lapply(1:(dim(reference)[2]-2), single_test, reference, test, ref_quality, test_quality)
  
  #Sampling analysis (inferential)
  sampling <- unlist(normal_sampling_report)[!is.na(unlist(normal_sampling_report))]
  sampling_size <- length(normal_sampling_report)
  pop_mean <- mean(sampling)
  pop_sd <- sd(sampling)/sampling_size #Correction factor ~1 since population size is enormous
  
  print("Investigating test data")
  test_report <- single_test(-1, reference, test, ref_quality, test_quality)
  test_report <- test_report/test_report[!is.na(test_report)][1]
  barplot(test_report, xlab = "Bin", ylab = "Copy Number Ratio", main = "Copy Number of Test Sample", ylim = c(0,2))
  
  z_score <- (test_report-pop_mean)/pop_sd
  
  plot(y = z_score, x = names(z_score), xlab = "Bin", ylab = "Z-score", main = "Z-score for Test Sample Bins", xlim = c(0,23))
  
  pvalue2sided <- 2*pnorm(-abs(z_score))
  
  print("Testing disomic hypothesis")
  lapply(1:length(pvalue2sided), report_cnv, p_value = pvalue2sided, z_score, test_report, pop_mean, pop_sd, names(z_score))
  
  print("...done")
  print("...done")
}