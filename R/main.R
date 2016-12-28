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

result <- run_full_test(reference_coverage, test_coverage, reference_quality, test_quality)
