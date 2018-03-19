tumor_list_snp <- list.files(path = "../Results/mutect/tumor_bedonly", pattern = "snp[.]oncotator[.]vcf$")
tumor_list_indel <- gsub("snp", "indel", tumor_list_snp)
tumor_patient <- substr(tumor_list_snp, 1, 5)
dir.create("../Results/mutect/tumor_clone/")

for(i in 1:length(tumor_list_snp)){
	  snp <- read.delim(paste("../Results/mutect/tumor_bedonly", tumor_list_snp[i], sep = "/"), header = T, stringsAsFactors = F)
  snp <- snp[, 1:11]
    colnames(snp)[10] <- tumor_patient[i]
    if(file.exists(paste("../Results/mutect/tumor_bedonly", tumor_list_indel[i], sep = "/"))){
	        indel <- read.delim(paste("../Results/mutect/tumor_bedonly", tumor_list_indel[i], sep = "/"), header = T, stringsAsFactors = F)
        indel <- indel[, 1:11]
	    colnames(indel)[10] <- tumor_patient[i]
	    snp <- rbind(snp, indel)
	      }
      write.table(snp, paste("../Results/mutect/tumor_clone/", tumor_list_snp[i], sep = ""), 
		                quote = F, row.names = F, sep = "\t")
}


lymp_list_snp <- list.files(path = "../Results/mutect/lymp_bedonly", pattern = "snp[.]oncotator[.]vcf$")
lymp_list_indel <- gsub("snp", "indel", lymp_list_snp)
lymp_patient <- substr(lymp_list_snp, 1, 5)
dir.create("../Results/mutect/lymp_clone/")

for(i in 1:length(lymp_list_snp)){
	  snp <- read.delim(paste("../Results/mutect/lymp_bedonly", lymp_list_snp[i], sep = "/"), header = T, stringsAsFactors = F)
  snp <- snp[, 1:11]
    colnames(snp)[10] <- lymp_patient[i]
    if(file.exists(paste("../Results/mutect/lymp_bedonly", lymp_list_indel[i], sep = "/"))){
	        indel <- read.delim(paste("../Results/mutect/lymp_bedonly", lymp_list_indel[i], sep = "/"), header = T, stringsAsFactors = F)
        indel <- indel[, 1:11]
	    colnames(indel)[10] <- lymp_patient[i]
	    snp <- rbind(snp, indel)
	      }
      write.table(snp, paste("../Results/mutect/lymp_clone/", lymp_list_snp[i], sep = ""), 
		                quote = F, row.names = F, sep = "\t")
}

