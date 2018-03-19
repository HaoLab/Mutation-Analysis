dir.create("../Results/mutect/tumor_bedonly/")
dir.create("../Results/mutect/lymp_bedonly/")

mutect_list <- read.table("../Summary/mutect.list", header = F, stringsAsFactors = F)
chr <- paste("chr", c(1:20), sep = "")
tumor_list <- subset(mutect_list, V1 == "tumor")
lymp_list <- subset(mutect_list, V1 == "lymp")

for(j in 1:nrow(tumor_list)){
  file_list <- list.files(path = paste("../Results/mutect/tumor_chr_bedonly/", tumor_list$V3[j], sep = ""),
                          pattern = "vcf$", full.names = T)
  chrfile_list <- paste("../Results/mutect/tumor_chr_bedonly/", tumor_list$V3[j], "/", 
                        tumor_list$V3[j], "_", chr, ".dedup.vcf", sep = "")
  chrfile_list <- chrfile_list[chrfile_list %in% file_list]
  for(n in 1:length(chrfile_list)){
    vcf <- read.table(chrfile_list[n], header = T, stringsAsFactors = F)
    if(n == 1){
      vcf_all <- vcf
    }else{
      vcf_all <- rbind(vcf_all, vcf)
    }
  }
  write.table(vcf_all, paste("../Results/mutect/tumor_bedonly/", tumor_list$V3[j], ".dedup.vcf", sep = ""), 
              sep = "\t", quote = F, row.names = F)
}


for(j in 1:nrow(lymp_list)){
  file_list <- list.files(path = paste("../Results/mutect/lymp_chr_bedonly/", lymp_list$V3[j], sep = ""),
                          pattern = "vcf$", full.names = T)
  chrfile_list <- paste("../Results/mutect/lymp_chr_bedonly/", lymp_list$V3[j], "/", 
                        lymp_list$V3[j], "_", chr, ".dedup.vcf", sep = "")
  chrfile_list <- chrfile_list[chrfile_list %in% file_list]
  for(n in 1:length(chrfile_list)){
    vcf <- read.table(chrfile_list[n], header = T, stringsAsFactors = F)
    if(n == 1){
      vcf_all <- vcf
    }else{
      vcf_all <- rbind(vcf_all, vcf)
    }
  }
  write.table(vcf_all, paste("../Results/mutect/lymp_bedonly/", lymp_list$V3[j], ".dedup.vcf", sep = ""), 
              sep = "\t", quote = F, row.names = F)
}


