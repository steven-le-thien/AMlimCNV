#Basic logistics
source("setting.R", chdir = TRUE)
script.dir <- getwd()
reference.data.dir <- paste(script.dir, "/data/reference_raw_data", sep = "")
test.data.dir <- paste(script.dir, "/data/test_raw_data", sep = "")

setwd(script.dir)
source("full_test.R")

#Generate data frame for each attribute (this step can be optimized further with a streaming-like process)
source("parsing_from_excel.R")
reference_quality <- read_data(reference.data.dir, "Quality")
reference_coverage <- read_data(reference.data.dir, "Coverage")
reference_original_coverage <- read_data(reference.data.dir, "Original.Coverage")
# reference_pool_cluster <- read_data(reference_data_dir, "Pool.Cluster")
test_quality <- read_data(test.data.dir, "Quality")
test_coverage <- read_data(test.data.dir, "Coverage")
test_original_coverage <- read_data(test.data.dir, "Original.Coverage")

#Sampling
print("Investigating population...")
normal_sampling_report <- lapply(1:(dim(reference_coverage)[2]-2), single_test, reference_coverage, test_coverage, reference_quality, test_quality)
print("...done")

#Sampling analysis (inferential)
sampling <- unlist(normal_sampling_report)[!is.na(unlist(normal_sampling_report))]
sampling_size <- length(normal_sampling_report)
pop_mean <- mean(sampling)
pop_sd <- sd(sampling)/sampling_size #Correction factor ~1 since population size is enormous

#Test disomic analysis
print("Investigating test data")
test_report <- single_test(-1, reference_coverage, test_coverage, reference_quality, test_quality)
test_report <- test_report/test_report[!is.na(test_report)][1]

#Plot copy number ratio
barplot(test_report, xlab = "Bin", ylab = "Copy Number Ratio", main = "Copy Number of Test Sample", ylim = c(0,2))

#Plot z_score
z_score <- (test_report-pop_mean)/pop_sd
plot(y = z_score, x = names(z_score), xlab = "Bin", ylab = "Z-score", main = "Z-score for Test Sample Bins", xlim = c(0,23))

#Generate p value
pvalue2sided <- 2*pnorm(-abs(z_score))

#Begin test series
print("Testing disomic hypothesis")
result <- lapply(1:length(pvalue2sided), report_cnv, p_value = pvalue2sided, z_score, test_report, pop_mean, pop_sd, names(z_score))

print("...done")

gc()

