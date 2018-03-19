tumor_list_snp <- list.files(path = "../Results/mutect/tumor_bedonly", pattern = "snp[.]oncotator[.]vcf$")
tumor_list_indel <- gsub("snp", "indel", tumor_list_snp)
tumor_patient <- substr(tumor_list_snp, 1, 5)
dir.create("../Results/mutect/tumor_foroncotator/")

for(i in 1:length(tumor_list_snp)){
  snp <- read.delim(paste("../Results/mutect/tumor_bedonly", tumor_list_snp[i], sep = "/"), header = T, stringsAsFactors = F)
  snp <- snp[, 1:10]
  colnames(snp)[10] <- tumor_patient[i]
  if(file.exists(paste("../Results/mutect/tumor_bedonly", tumor_list_indel[i], sep = "/"))){
    indel <- read.delim(paste("../Results/mutect/tumor_bedonly", tumor_list_indel[i], sep = "/"), header = T, stringsAsFactors = F)
    indel <- indel[, 1:10]
    colnames(indel)[10] <- tumor_patient[i]
    snp <- rbind(snp, indel)
  }
  write.table(snp, paste("../Results/mutect/tumor_foroncotator/", tumor_list_snp[i], sep = ""), 
              quote = F, row.names = F, sep = "\t")
}


lymp_list_snp <- list.files(path = "../Results/mutect/lymp_bedonly", pattern = "snp[.]oncotator[.]vcf$")
lymp_list_indel <- gsub("snp", "indel", lymp_list_snp)
lymp_patient <- substr(lymp_list_snp, 1, 5)
dir.create("../Results/mutect/lymp_foroncotator/")

for(i in 1:length(lymp_list_snp)){
  snp <- read.delim(paste("../Results/mutect/lymp_bedonly", lymp_list_snp[i], sep = "/"), header = T, stringsAsFactors = F)
  snp <- snp[, 1:10]
  colnames(snp)[10] <- lymp_patient[i]
  if(file.exists(paste("../Results/mutect/lymp_bedonly", lymp_list_indel[i], sep = "/"))){
    indel <- read.delim(paste("../Results/mutect/lymp_bedonly", lymp_list_indel[i], sep = "/"), header = T, stringsAsFactors = F)
    indel <- indel[, 1:10]
    colnames(indel)[10] <- lymp_patient[i]
    snp <- rbind(snp, indel)
  }
  write.table(snp, paste("../Results/mutect/lymp_foroncotator/", lymp_list_snp[i], sep = ""), 
              quote = F, row.names = F, sep = "\t")
}



tumor_list <- list.files(path = "../Results/mutect/tumor_bedonly", pattern = "dedup[.]vcf$")
tumor_list <- tumor_list[-c(grep(pattern = "snp", tumor_list), grep(pattern = "indel", tumor_list))]
tumor_patient <- substr(tumor_list, 1, 5)
dir.create("../Results/mutect/tumor_foroncotator_all/")

for(i in 1:length(tumor_list)){
  snp <- read.delim(paste("../Results/mutect/tumor_bedonly", tumor_list[i], sep = "/"), header = T, stringsAsFactors = F)
  snp <- subset(snp, FILTER == "PASS")
  snp <- snp[, 1:10]
  colnames(snp)[10] <- tumor_patient[i]
  write.table(snp, paste("../Results/mutect/tumor_foroncotator_all/", tumor_list[i], sep = ""), 
              quote = F, row.names = F, sep = "\t")
}




lymp_list <- list.files(path = "../Results/mutect/lymp_bedonly", pattern = "dedup[.]vcf$")
lymp_list <- lymp_list[-c(grep(pattern = "snp", lymp_list), grep(pattern = "indel", lymp_list))]
lymp_patient <- substr(lymp_list, 1, 5)
dir.create("../Results/mutect/lymp_foroncotator_all/")

for(i in 1:length(lymp_list)){
  snp <- read.delim(paste("../Results/mutect/lymp_bedonly", lymp_list[i], sep = "/"), header = T, stringsAsFactors = F)
  snp <- subset(snp, FILTER == "PASS")
  snp <- snp[, 1:10]
  colnames(snp)[10] <- lymp_patient[i]
  write.table(snp, paste("../Results/mutect/lymp_foroncotator_all/", lymp_list[i], sep = ""), 
              quote = F, row.names = F, sep = "\t")
}




tumor_list_snp <- list.files(path = "../Results/mutect/tumor_bedonly", pattern = "snp[.]oncotator[.]vcf$")
tumor_list_indel <- gsub("snp", "indel", tumor_list_snp)
tumor_patient <- substr(tumor_list_snp, 1, 5)
dir.create("../Results/mutect/tumor_foroncotator_001/")

for(i in 1:length(tumor_list_snp)){
  snp <- read.delim(paste("../Results/mutect/tumor_bedonly", tumor_list_snp[i], sep = "/"), header = T, stringsAsFactors = F)
  snp <- subset(snp, T_AF >= 0.01)
  snp <- snp[, 1:10]
  colnames(snp)[10] <- tumor_patient[i]
  if(file.exists(paste("../Results/mutect/tumor_bedonly", tumor_list_indel[i], sep = "/"))){
    indel <- read.delim(paste("../Results/mutect/tumor_bedonly", tumor_list_indel[i], sep = "/"), header = T, stringsAsFactors = F)
    indel <- subset(indel, T_AF >= 0.01)
    indel <- indel[, 1:10]
    colnames(indel)[10] <- tumor_patient[i]
    snp <- rbind(snp, indel)
  }
  write.table(snp, paste("../Results/mutect/tumor_foroncotator_001/", tumor_list_snp[i], sep = ""), 
              quote = F, row.names = F, sep = "\t")
}


lymp_list_snp <- list.files(path = "../Results/mutect/lymp_bedonly", pattern = "snp[.]oncotator[.]vcf$")
lymp_list_indel <- gsub("snp", "indel", lymp_list_snp)
lymp_patient <- substr(lymp_list_snp, 1, 5)
dir.create("../Results/mutect/lymp_foroncotator_001/")

for(i in 1:length(lymp_list_snp)){
  snp <- read.delim(paste("../Results/mutect/lymp_bedonly", lymp_list_snp[i], sep = "/"), header = T, stringsAsFactors = F)
  snp <- subset(snp, T_AF >= 0.01)
  snp <- snp[, 1:10]
  colnames(snp)[10] <- lymp_patient[i]
  if(file.exists(paste("../Results/mutect/lymp_bedonly", lymp_list_indel[i], sep = "/"))){
    indel <- read.delim(paste("../Results/mutect/lymp_bedonly", lymp_list_indel[i], sep = "/"), header = T, stringsAsFactors = F)
    indel <- subset(indel, T_AF >= 0.01)
    indel <- indel[, 1:10]
    colnames(indel)[10] <- lymp_patient[i]
    snp <- rbind(snp, indel)
  }
  write.table(snp, paste("../Results/mutect/lymp_foroncotator_001/", lymp_list_snp[i], sep = ""), 
              quote = F, row.names = F, sep = "\t")
}