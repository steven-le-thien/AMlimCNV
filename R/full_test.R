#Stand-in function that map the chromosome to its numeric position, X is 24 and Y is 25
chrom_position <- function(x) {
  bin_no <- substr(x, 4, nchar(x))
  if (bin_no == 'X')
    bin_no <- 24
  if (bin_no == 'Y')
    bin_no <- 25
  return (bin_no)
}

#A single iteration of the test, including BAE1, BAE2, ratio analysis
single_test <- function (c, reference, test, ref_quality, test_quality) {

  #Determine whether it is in the process of sampling or it is the final test case
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

  output <- ratio_analysis(coverage_bae1_bae2)
  print("...done")

  return(output)
}

#A logic gate that determines which hypothesis to test for and print corresponding call
report_cnv <- function (order, p_value, z_score, test_report, pop_mean, pop_sd, names) {
  if (!is.na(p_value[order])) {#skip NA values
    if (p_value[order] < p.value.threshold) {#if disomic hypo is not accepted
      print(paste("Reject disomic hypothesis, bin", names[order], "does not have copy number 2"))
      if (z_score[order] > 0) {#if there seems to be more than 2 copy
        print ("Moving on to trisomic hypothesis testing")
        tester <- test_report[order] / 1.5
        if (2*pnorm(-abs((tester-pop_mean)/pop_sd)) < p.value.threshold) {#if trisomic hypo is not accepted
          print(paste("Reject trisomic hypothesis, bin", names[order], "does not have copy number 3"))
          print("Test is inconclusive")
        } else {#do not reject trisomic hypothesis
          print(paste("Bin", names[order], "has copy number 3"))
        }
      }
      else {#if there seems to be less than 2 copy
        print ("Moving on to monosomic hypothesis testing")
        tester <- test_report[order] / 0.5
        if (2*pnorm(-abs((tester-pop_mean)/pop_sd)) < p.value.threshold) {#do not accept monosomic hypothesis
          print(paste("Reject monosomic hypothesis, bin", names[order], "does not have copy number 1"))
          print("Test is inconclusive")
        } else {#do not reject monosomic hypothesis
          print(paste("Bin", names[order], "has copy number 1"))
        }
      }
    } else {
      print(paste("Bin", names[order], "has copy number 2"))
    }
  }
}
