# load libraries
suppressPackageStartupMessages({
  library(URD)
})

# Find the most recent test URD object
test_files <- list.files("test_data", pattern = "test_urd_object_.*\\.rds", full.names = TRUE)
if (length(test_files) == 0) {
  stop("No test URD object found.")
}
latest_test_file <- test_files[which.max(file.info(test_files)$mtime)]

# Load the test URD object
message(sprintf("Loading test URD object from: %s", latest_test_file))
urd_object <- readRDS(latest_test_file)

# Check object structure
message("\nURD Object Structure:")
message("Count data dimensions: ", paste(dim(urd_object@count.data), collapse=" x "))
message("Logupx data dimensions: ", paste(dim(urd_object@logupx.data), collapse=" x "))
message("Number of cells: ", ncol(urd_object@count.data))
message("Number of genes: ", nrow(urd_object@count.data))

# Check if data slots are empty
message("\nData slot status:")
message("Count data empty: ", is.null(urd_object@count.data))
message("Logupx data empty: ", is.null(urd_object@logupx.data))
message("Group IDs empty: ", is.null(urd_object@group.ids))

# Print first few gene names if available
if (!is.null(rownames(urd_object@count.data))) {
    message("\nFirst few genes:")
    message(paste(head(rownames(urd_object@count.data)), collapse=", "))
}

# Print first few cell names if available
if (!is.null(colnames(urd_object@count.data))) {
    message("\nFirst few cells:")
    message(paste(head(colnames(urd_object@count.data)), collapse=", "))
} 