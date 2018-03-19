library(getopt)
opt_spec <- matrix(c('help', 'h', 0, 'logical', 'help manual',
                     'samplesize', 'n', 1, 'numeric', 'Sample size',
                     'maf', 'm', 1, 'character', 'maf file',
                     'mutsig', 's', 1, 'character', 'out put of MutSigCV',
                     'output', 'o', 1, 'character', 'output path',
                     'prefix', 'p', 1, 'character', 'prefix of your output files'),
                   byrow=TRUE, ncol=5)
opt = getopt(opt_spec, commandArgs(TRUE))
if (!is.null(opt$help)){
  cat(getopt(opt_spec, usage=TRUE))
  q(save='no', status=1)
}
stopifnot(file.exists(opt$maf, opt$mutsig, opt$output))


# opt <- data.frame(maf = "../Analysis/lymp_001.maf", 
#                   mutsig = "../Analysis/lymp/lymp_001.sig_genes.txt", 
#                   output = "../Analysis/hotnet2/",
#                   prefix = "Lymp_001",
#                   stringsAsFactors = F)
maf <- read.delim(opt$maf, header = T, stringsAsFactors = F)
mutsig <- read.delim(opt$mutsig, header = T, stringsAsFactors = F)
if (is.null(opt$samplesize)){
  opt$samplesize <- length(unique(maf$Tumor_Sample_Barcode))
}

freq <- table(maf$Hugo_Symbol)
freq <- data.frame(gene = names(freq), freq = c(freq)/opt$samplesize, stringsAsFactors = F)
freq <- subset(freq, freq >= 3/opt$samplesize)
write.table(freq, paste(opt$output, "/", opt$prefix, "_freq.txt", sep = ""),
            quote = F, col.names = F, row.names = F, sep = "\t")

mutsig <- data.frame(gene = mutsig$gene, sig = -log10(mutsig$q))
write.table(mutsig, paste(opt$output, "/", opt$prefix, "_mutsig.txt", sep = ""),
            quote = F, col.names = F, row.names = F, sep = "\t")
