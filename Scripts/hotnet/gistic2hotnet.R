library(getopt)
opt_spec <- matrix(c('help', 'h', 0, 'logical', 'help manual',
                     'samplesize', 'n', 1, 'numeric', 'Sample size',
                     'input', 'i', 1, 'character', 'scores.gistic',
                     'output', 'o', 1, 'character', 'output path',
                     'prefix', 'p', 1, 'character', 'prefix of your output files'),
                   byrow=TRUE, ncol=5)
opt = getopt(opt_spec, commandArgs(TRUE))
if (!is.null(opt$help)){
  cat(getopt(opt_spec, usage=TRUE))
  q(save='no', status=1)
}
stopifnot(file.exists(opt$input, opt$output))


# opt <- data.frame(input = "../Results/hotnet2/tumor/gistic_tumor_001.scores", 
#                   output = "../Results/hotnet2/tumor",
#                   prefix = "tumor_001",
#                   stringsAsFactors = F)
gistic <- read.delim(opt$input, header = F, stringsAsFactors = F)
colnames(gistic) <- c("Chromosome", "Start", "End", "Type", "logq_negative", "gscore", "averageAMP", "frequency", "Gene")

gistic_amp <- subset(gistic, Type == "Amp" & frequency >= 3/opt$samplesize)
gistic_del <- subset(gistic, Type == "Del" & frequency >= 3/opt$samplesize)

genelist <- unique(gistic_amp$Gene)[-which(unique(gistic_amp$Gene) == "-")]
data <- data.frame(gene = genelist, freq = NA, logq = NA, stringsAsFactors = F)
for(i in 1:nrow(data)){
  tmp <- subset(gistic_amp, Gene == genelist[i])
  if(length(unique(tmp$logq_negative)) == 1){
    data$logq[i] <- unique(tmp$logq_negative)
  }else{
    logq_tmp <- data.frame(logq = unique(tmp$logq_negative), len = NA, stringsAsFactors = F)
    for(k in 1:nrow(logq_tmp)){
      logq_tmp$len[k] <- max(as.numeric(tmp$End[which(tmp$logq_negative == logq_tmp$logq[k])])) - min(as.numeric(tmp$Start[which(tmp$logq_negative == logq_tmp$logq[k])]))
    }
    data$logq[i] <- logq_tmp$logq[which.max(logq_tmp$len)]
  }
  if(length(unique(tmp$frequency)) == 1){
    data$freq[i] <- unique(tmp$frequency)
  }else{
    freq_tmp <- data.frame(freq = unique(tmp$freq), len = NA, stringsAsFactors = F)
    for(k in 1:nrow(freq_tmp)){
      freq_tmp$len[k] <- max(as.numeric(tmp$End[which(tmp$freq == freq_tmp$freq[k])])) - min(as.numeric(tmp$Start[which(tmp$freq == freq_tmp$freq[k])]))
    }
    data$freq[i] <- freq_tmp$freq[which.max(freq_tmp$len)]
  }
}

write.table(data[,c(1,3)], paste(opt$output, "/", opt$prefix, "_gistic_amp.txt", sep = ""),
            quote = F, col.names = F, row.names = F, sep = "\t")

data <- subset(data, freq >= 3/opt$samplesize)
write.table(data[,c(1,2)], paste(opt$output, "/", opt$prefix, "_freq_amp.txt", sep = ""),
            quote = F, col.names = F, row.names = F, sep = "\t")
			
			data <- data.frame(gene = genelist, freq = NA, logq = NA, stringsAsFactors = F)
			
genelist <- unique(gistic_del$Gene)[-which(unique(gistic_del$Gene) == "-")]	
data <- data.frame(gene = genelist, freq = NA, logq = NA, stringsAsFactors = F)		
for(i in 1:nrow(data)){
  tmp <- subset(gistic_del, Gene == genelist[i])
  if(length(unique(tmp$logq_negative)) == 1){
    data$logq[i] <- unique(tmp$logq_negative)
  }else{
    logq_tmp <- data.frame(logq = unique(tmp$logq_negative), len = NA, stringsAsFactors = F)
    for(k in 1:nrow(logq_tmp)){
      logq_tmp$len[k] <- max(as.numeric(tmp$End[which(tmp$logq_negative == logq_tmp$logq[k])])) - min(as.numeric(tmp$Start[which(tmp$logq_negative == logq_tmp$logq[k])]))
    }
    data$logq[i] <- logq_tmp$logq[which.max(logq_tmp$len)]
  }
  if(length(unique(tmp$frequency)) == 1){
    data$freq[i] <- unique(tmp$frequency)
  }else{
    freq_tmp <- data.frame(freq = unique(tmp$freq), len = NA, stringsAsFactors = F)
    for(k in 1:nrow(freq_tmp)){
      freq_tmp$len[k] <- max(as.numeric(tmp$End[which(tmp$freq == freq_tmp$freq[k])])) - min(as.numeric(tmp$Start[which(tmp$freq == freq_tmp$freq[k])]))
    }
    data$freq[i] <- freq_tmp$freq[which.max(freq_tmp$len)]
  }
}

write.table(data[,c(1,3)], paste(opt$output, "/", opt$prefix, "_gistic_del.txt", sep = ""),
            quote = F, col.names = F, row.names = F, sep = "\t")

data <- subset(data, freq >= 3/opt$samplesize)
write.table(data[,c(1,2)], paste(opt$output, "/", opt$prefix, "_freq_del.txt", sep = ""),
            quote = F, col.names = F, row.names = F, sep = "\t")

