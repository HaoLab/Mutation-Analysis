library(getopt)
opt_spec <- matrix(c('help', 'h', 0, 'logical', 'help manual',
                     'folder', 'f', 1, 'character', 'sample path file, 2 column, name and folder',
                     'output', 'o', 1, 'character', 'output path'),
                   byrow=TRUE, ncol=5)
opt = getopt(opt_spec, commandArgs(TRUE))
if (!is.null(opt$help)){
  cat(getopt(opt_spec, usage=TRUE))
  q(save='no', status=1)
}
stopifnot(file.exists(opt$folder, opt$output))

folder <- read.delim(opt$folder, header = F, stringsAsFactors = F)
table <- data.frame(PID = folder$V1, depth = NA, depth_dedup = NA, dup_ratio = NA, uniformity = NA,
                    coverage_100x = NA, coverage_500x = NA, coverage_1000x = NA, 
                    coverage_5000x = NA, coverage_10000x = NA,
                    stringsAsFactors = F)
cat("qc_table\n")
for(j in 1:nrow(table)){
  filepath <- paste(folder$V2[j], "/qc_sample", sep = "")
  filepath <- list.files(filepath, pattern = "qc_sample[.]txt$", full.names = T)
  qc <- read.delim(filepath, header = T, stringsAsFactors = F)
  table$depth[j] <- qc$average_depth
  if("dedup_average_depth" %in% colnames(qc)){
  table$depth_dedup[j] <- qc$dedup_average_depth
}
  table$coverage_100x[j] <- paste(round(qc$coverage_100x*100, 2), "%", sep = "")
  table$coverage_500x[j] <- paste(round(qc$coverage_500x*100, 2), "%", sep = "")
  table$coverage_1000x[j] <- paste(round(qc$coverage_1000x*100, 2), "%", sep = "")
  table$coverage_5000x[j] <- paste(round(qc$coverage_5000x*100, 2), "%", sep = "")
  table$coverage_10000x[j] <- paste(round(qc$coverage_10000x*100, 2), "%", sep = "")
}
cat("dedup_depth\n")
for(j in 1:nrow(table)){
  filepath <- paste(folder$V2[j], "/qc_sample", sep = "")
  filepath <- list.files(filepath, pattern = "dedup_depth[.]txt$", full.names = T)
  if(length(filepath) != 0){
    qc <- read.delim(filepath, header = T, stringsAsFactors = F)
    table$depth_dedup[j] <- qc$average_depth
  }
}


cat("uniformity\n")
for(j in 1:nrow(table)){
  filepath <- paste(folder$V2[j], "/qc_sample", sep = "")
  filepath <- list.files(filepath, pattern = "uniformity[.]txt$", full.names = T)
  if(length(filepath) != 0){
  qc <- read.delim(filepath, header = T, stringsAsFactors = F)
  table$uniformity[j] <- paste(round(qc$uniformity*100, 2), "%", sep = "")
  }
}
cat("metrics\n")
for(j in 1:nrow(table)){
  filepath <- paste(folder$V2[j], "/basic_analysis/align", sep = "")
  filepath <- list.files(filepath, pattern = "sorted[.]mkdup[.]metrics$", full.names = T)
  if(length(filepath) != 0){
    qc <- read.delim(filepath, header = T, stringsAsFactors = F, nrows = 1, fill = T, comment.char="#")
    table$dup_ratio[j] <- paste(round(qc$PERCENT_DUPLICATION*100, 2), "%", sep = "")
  }
}

write.table(table, paste(opt$output, "/QC_pipeline.txt", sep = ""), 
            quote = F, row.names = F, sep = "\t")


