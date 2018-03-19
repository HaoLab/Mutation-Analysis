sample <- read.table("../Summary/paired.info", stringsAsFactors = F, header = T)

if(file.exists("/lustre/rdi/user/yujn/data/Normal_30_CF/Normal_CF_Dedup.txt")){
  Normal <- read.delim("/lustre/rdi/user/yujn/data/Normal_30_CF/Normal_CF_Dedup.txt", header = T, stringsAsFactors = F)
}else{
  Normal <- read.delim("/work-d/user/hantc/lustre/rdi/user/yujn/data/Normal_30_CF/Normal_CF_Dedup.txt", header = T, stringsAsFactors = F)
}
#Normal <- read.delim("../Summary/Normal_CF_Dedup.txt", header = T, stringsAsFactors = F)
Normal <- subset(Normal, Freq >= 0.2)
Normal$SNV <- paste(Normal$Chrom,Normal$Pos,Normal$Ref,Normal$Alt, sep = ":")

black <- readRDS("../../../tools/blacklist.rds")

table_snp <- data.frame(tumor = sample$tumor, tumor_all = 0, tumor_f1 = 0, tumor_f2 = 0, 
                        lymp = sample$lymp, lymp_all = 0, lymp_f1 = 0, lymp_f2 = 0,
                        overlap_all = 0, overlap_filtered1 = 0, overlap_filtered2 = 0)
table_indel <- data.frame(tumor = sample$tumor, tumor_all = 0, tumor_f1 = 0, tumor_f2 = 0, 
                          lymp = sample$lymp, lymp_all = 0, lymp_f1 = 0, lymp_f2 = 0,
                          overlap_all = 0, overlap_filtered1 = 0, overlap_filtered2 = 0)

for (i in 1:nrow(sample)){
  if(file.exists(paste("../Results/mutect/tumor_bedonly/", sample$tumor[i], ".dedup.vcf", sep = ""))){
    tumor <- read.delim(paste("../Results/mutect/tumor_bedonly/", sample$tumor[i], ".dedup.vcf", sep = ""),
                        header = T, stringsAsFactors = F)
    colnames(tumor)[1] <- "CHROM"
    tumor <- subset(tumor, FILTER == "PASS")
    tumor$SNV <- paste(tumor$CHROM,tumor$POS,tumor$REF,tumor$ALT, sep = ":")
    if (nrow(tumor) > 0){
      if(length(grep(tumor$FORMAT, pattern = "PID")) != 0){
        tumor_set1 <- tumor[grep(tumor$FORMAT, pattern = "PID"),]
        tumor_set2 <- tumor[-grep(tumor$FORMAT, pattern = "PID"),]
      }else{
        tumor_set1 <- tumor[0,]
        tumor_set2 <- tumor
      }
      n <- ncol(tumor)
      if(nrow(tumor_set1) > 0){
        colnamelist <- unlist(strsplit(tumor_set1$FORMAT[1], split = ":"))
        colnamelist_N <- paste("N", colnamelist, sep = "_")
        colnamelist_T <- paste("T", colnamelist, sep = "_")
        for (j in 1:11){
          tumor_set1[,(j+n)] <- unlist(strsplit(tumor_set1$NORMAL, split = ":"))[seq(j, (11*(nrow(tumor_set1)-1)+j), 11)]
        }
        for (j in 1:11){
          tumor_set1[,(j+n+11)] <- unlist(strsplit(tumor_set1$TUMOR, split = ":"))[seq(j, (11*(nrow(tumor_set1)-1)+j), 11)]
        }
        colnames(tumor_set1)[(n+1):(n+22)] <- c(colnamelist_N, colnamelist_T)
      }
      if(nrow(tumor_set2) > 0){
        colnamelist <- unlist(strsplit(tumor_set2$FORMAT[1], split = ":"))
        colnamelist_N <- paste("N", colnamelist, sep = "_")
        colnamelist_T <- paste("T", colnamelist, sep = "_")  
        for (j in 1:6){
          tumor_set2[,(j+n)] <- unlist(strsplit(tumor_set2$NORMAL, split = ":"))[seq(j, (9*(nrow(tumor_set2)-1)+j), 9)]
        }
        for (j in 7:8){
          tumor_set2[,(j+n)] <- NA
        }
        for (j in 9:11){
          tumor_set2[,(j+n)] <- unlist(strsplit(tumor_set2$NORMAL, split = ":"))[seq((j-2), (9*(nrow(tumor_set2)-1)+j-2), 9)]
        }
        for (j in 1:6){
          tumor_set2[,(j+n+11)] <- unlist(strsplit(tumor_set2$TUMOR, split = ":"))[seq(j, (9*(nrow(tumor_set2)-1)+j), 9)]
        }
        for (j in 7:8){
          tumor_set2[,(j+n+11)] <- NA
        }
        for (j in 9:11){
          tumor_set2[,(j+n+11)] <- unlist(strsplit(tumor_set2$TUMOR, split = ":"))[seq((j-2), (9*(nrow(tumor_set2)-1)+j-2), 9)]
        }
        colnames(tumor_set2)[c((n+1):(n+6), (n+9):(n+11), (n+12):(n+17), (n+20):(n+22))] <- c(colnamelist_N, colnamelist_T)
        colnames(tumor_set2)[c((n+7):(n+8), (n+18):(n+19))] <- c("N_PGT", "N_PID", "T_PGT", "T_PID")
      }
      tumor <- rbind(tumor_set1, tumor_set2)
      tumor <- tumor[order(as.integer(rownames(tumor))),]
      for(j in c(15:18,22,23,25:29,33,34)){
        tumor[,j] <- as.numeric(tumor[,j])
      }
    }
    if (nrow(tumor) > 0){
      tumor2 <- read.delim(paste("../Results/mutect/tumor_bedonly/anno/dedup_", sample$tumor[i], ".hg19_multianno.txt",
                                 sep = ""),header = F, stringsAsFactors = F)
      colnames(tumor2)[1:58] <- tumor2[1,1:58]
      tumor2 <- tumor2[-1,6:58]
      rownames(tumor2) <- c(1:nrow(tumor2))
      tumor2 <- tumor2[rownames(tumor),]
      tumor <- cbind(tumor, tumor2)
      
      tumor_snp <- subset(tumor, nchar(REF) == 1 & nchar(ALT) == 1)
      tumor_indel <- subset(tumor, nchar(REF) != 1 | nchar(ALT) != 1)
      write.table(tumor_snp, paste("../Results/mutect/tumor_bedonly/", sample$tumor[i], ".snp.oncotator.vcf", sep = ""),
                  quote = F, row.names = F, sep = "\t")
      tumor_snp$ExAC_ALL <- as.numeric(tumor_snp$ExAC_ALL)
      tumor_snp$thousand <- as.numeric(tumor_snp$`1000g2015aug_all`)
      tumor_snp <- subset(tumor_snp, N_AF < 0.4)
      tumor_snp <- subset(tumor_snp, SNV %in% black$MUT == FALSE)
      tumor_snp <- subset(tumor_snp, SNV %in% Normal$SNV == FALSE & (ExAC_ALL < 0.005 | is.na(ExAC_ALL) == TRUE) &
                            (thousand < 0.005 | is.na(thousand) == TRUE))
      
      tumor_snp <- subset(tumor_snp, ExonicFunc.refGene != "synonymous SNV" & Func.refGene == "exonic" & 
                            ExonicFunc.refGene != "." & ExonicFunc.refGene != "unknown")
      tumor_snp_filtered1 <- subset(tumor_snp, T_AF >= 0.03)
      tumor_snp_filtered2 <- subset(tumor_snp_filtered1, T_AF >= 0.05)
      table_snp$tumor_all[i] = nrow(tumor_snp)
      table_snp$tumor_f1[i] = nrow(tumor_snp_filtered1)
      table_snp$tumor_f2[i] = nrow(tumor_snp_filtered2)
      
      tumor_indel$ExAC_ALL <- as.numeric(tumor_indel$ExAC_ALL)
      tumor_indel$thousand <- as.numeric(tumor_indel$`1000g2015aug_all`)
      tumor_indel <- subset(tumor_indel, N_AF < 0.4)
      tumor_indel <- subset(tumor_indel, SNV %in% black$MUT == FALSE)
      tumor_indel <- subset(tumor_indel, SNV %in% Normal$SNV == FALSE & (ExAC_ALL < 0.005 | is.na(ExAC_ALL) == TRUE) &
                              (thousand < 0.005 | is.na(thousand) == TRUE))
      write.table(tumor_indel, paste("../Results/mutect/tumor_bedonly/", sample$tumor[i], ".indel.oncotator.vcf", sep = ""),
                  quote = F, row.names = F, sep = "\t")
      tumor_indel <- subset(tumor_indel, Func.refGene == "exonic" &
                              ExonicFunc.refGene != "." & ExonicFunc.refGene != "unknown")
      tumor_indel_filtered1 <- subset(tumor_indel, T_AF >= 0.03)
      tumor_indel_filtered2 <- subset(tumor_indel_filtered1, T_AF >= 0.05)
      table_indel$tumor_all[i] = nrow(tumor_indel)
      table_indel$tumor_f1[i] = nrow(tumor_indel_filtered1)
      table_indel$tumor_f2[i] = nrow(tumor_indel_filtered2)
      write.table(tumor_snp, paste("../Results/mutect/tumor_bedonly/", sample$tumor[i], ".snp.dedup.vcf", sep = ""),
                  quote = F, row.names = F, sep = "\t")
      write.table(tumor_indel, paste("../Results/mutect/tumor_bedonly/", sample$tumor[i], ".indel.dedup.vcf", sep = ""),
                  quote = F, row.names = F, sep = "\t")
      write.table(tumor_snp_filtered2, paste("../Results/mutect/tumor_bedonly/", sample$tumor[i], ".snp.dedup.filtered.vcf", sep = ""),
                  quote = F, row.names = F, sep = "\t")
      write.table(tumor_indel_filtered2, paste("../Results/mutect/tumor_bedonly/", sample$tumor[i], ".indel.dedup.filtered.vcf", sep = ""),
                  quote = F, row.names = F, sep = "\t")
    }
  }
  
  if(file.exists(paste("../Results/mutect/lymp_bedonly/", sample$lymp[i], ".dedup.vcf", sep = ""))){
    lymp <- read.delim(paste("../Results/mutect/lymp_bedonly/", sample$lymp[i], ".dedup.vcf", sep = ""),
                       header = T, stringsAsFactors = F)
    colnames(lymp)[1] <- "CHROM"
    lymp <- subset(lymp, FILTER == "PASS")
    lymp$SNV <- paste(lymp$CHROM,lymp$POS,lymp$REF,lymp$ALT, sep = ":")
    if (nrow(lymp) > 0){
      if(length(grep(lymp$FORMAT, pattern = "PID")) != 0){
        lymp_set1 <- lymp[grep(lymp$FORMAT, pattern = "PID"),]
        lymp_set2 <- lymp[-grep(lymp$FORMAT, pattern = "PID"),]
      }else{
        lymp_set1 <- lymp[0,]
        lymp_set2 <- lymp
      }
      n <- ncol(lymp)
      if(nrow(lymp_set1) > 0){
        colnamelist <- unlist(strsplit(lymp_set1$FORMAT[1], split = ":"))
        colnamelist_N <- paste("N", colnamelist, sep = "_")
        colnamelist_T <- paste("T", colnamelist, sep = "_")
        for (j in 1:11){
          lymp_set1[,(j+n)] <- unlist(strsplit(lymp_set1$NORMAL, split = ":"))[seq(j, (11*(nrow(lymp_set1)-1)+j), 11)]
        }
        for (j in 1:11){
          lymp_set1[,(j+n+11)] <- unlist(strsplit(lymp_set1$lymp, split = ":"))[seq(j, (11*(nrow(lymp_set1)-1)+j), 11)]
        }
        colnames(lymp_set1)[(n+1):(n+22)] <- c(colnamelist_N, colnamelist_T)
      }
      if(nrow(lymp_set2) > 0){
        colnamelist <- unlist(strsplit(lymp_set2$FORMAT[1], split = ":"))
        colnamelist_N <- paste("N", colnamelist, sep = "_")
        colnamelist_T <- paste("T", colnamelist, sep = "_")  
        for (j in 1:6){
          lymp_set2[,(j+n)] <- unlist(strsplit(lymp_set2$NORMAL, split = ":"))[seq(j, (9*(nrow(lymp_set2)-1)+j), 9)]
        }
        for (j in 7:8){
          lymp_set2[,(j+n)] <- NA
        }
        for (j in 9:11){
          lymp_set2[,(j+n)] <- unlist(strsplit(lymp_set2$NORMAL, split = ":"))[seq((j-2), (9*(nrow(lymp_set2)-1)+j-2), 9)]
        }
        for (j in 1:6){
          lymp_set2[,(j+n+11)] <- unlist(strsplit(lymp_set2$lymp, split = ":"))[seq(j, (9*(nrow(lymp_set2)-1)+j), 9)]
        }
        for (j in 7:8){
          lymp_set2[,(j+n+11)] <- NA
        }
        for (j in 9:11){
          lymp_set2[,(j+n+11)] <- unlist(strsplit(lymp_set2$lymp, split = ":"))[seq((j-2), (9*(nrow(lymp_set2)-1)+j-2), 9)]
        }
        colnames(lymp_set2)[c((n+1):(n+6), (n+9):(n+11), (n+12):(n+17), (n+20):(n+22))] <- c(colnamelist_N, colnamelist_T)
        colnames(lymp_set2)[c((n+7):(n+8), (n+18):(n+19))] <- c("N_PGT", "N_PID", "T_PGT", "T_PID")
      }
      lymp <- rbind(lymp_set1, lymp_set2)
      lymp <- lymp[order(as.integer(rownames(lymp))),]
      for(j in c(15:18,22,23,25:29,33,34)){
        lymp[,j] <- as.numeric(lymp[,j])
      }
    }
    if (nrow(lymp) > 0){
      lymp2 <- read.delim(paste("../Results/mutect/lymp_bedonly/anno/dedup_", sample$lymp[i], ".hg19_multianno.txt",
                                sep = ""),header = F, stringsAsFactors = F)
      colnames(lymp2)[1:58] <- lymp2[1,1:58]
      lymp2 <- lymp2[-1,6:58]
      rownames(lymp2) <- c(1:nrow(lymp2))
      lymp2 <- lymp2[rownames(lymp),]
      lymp <- cbind(lymp, lymp2)
      
      lymp_snp <- subset(lymp, nchar(REF) == 1 & nchar(ALT) == 1)
      lymp_indel <- subset(lymp, nchar(REF) != 1 | nchar(ALT) != 1)
      
      lymp_snp$ExAC_ALL <- as.numeric(lymp_snp$ExAC_ALL)
      lymp_snp$thousand <- as.numeric(lymp_snp$`1000g2015aug_all`)
      lymp_snp <- subset(lymp_snp, N_AF < 0.4)
      lymp_snp <- subset(lymp_snp, SNV %in% black$MUT == FALSE)
      lymp_snp <- subset(lymp_snp, SNV %in% Normal$SNV == FALSE & (ExAC_ALL < 0.005 | is.na(ExAC_ALL) == TRUE) &
                           (thousand < 0.005 | is.na(thousand) == TRUE))
      write.table(lymp_snp, paste("../Results/mutect/lymp_bedonly/", sample$lymp[i], ".snp.oncotator.vcf", sep = ""),
                  quote = F, row.names = F, sep = "\t")
      
      lymp_snp <- subset(lymp_snp, ExonicFunc.refGene != "synonymous SNV" & Func.refGene == "exonic" & 
                           ExonicFunc.refGene != "." & ExonicFunc.refGene != "unknown" )
      lymp_snp_filtered1 <- subset(lymp_snp, T_AF >= 0.03)
      lymp_snp_filtered2 <- subset(lymp_snp_filtered1, T_AF >= 0.05)
      table_snp$lymp_all[i] = nrow(lymp_snp)
      table_snp$lymp_f1[i] = nrow(lymp_snp_filtered1)
      table_snp$lymp_f2[i] = nrow(lymp_snp_filtered2)
      
      lymp_indel$ExAC_ALL <- as.numeric(lymp_indel$ExAC_ALL)
      lymp_indel$thousand <- as.numeric(lymp_indel$`1000g2015aug_all`)
      lymp_indel <- subset(lymp_indel, N_AF < 0.4)
      lymp_indel <- subset(lymp_indel, SNV %in% black$MUT == FALSE)
      lymp_indel <- subset(lymp_indel, SNV %in% Normal$SNV == FALSE & (ExAC_ALL < 0.005 | is.na(ExAC_ALL) == TRUE) &
                             (thousand < 0.005 | is.na(thousand) == TRUE))
      write.table(lymp_indel, paste("../Results/mutect/lymp_bedonly/", sample$lymp[i], ".indel.oncotator.vcf", sep = ""),
                  quote = F, row.names = F, sep = "\t")
      lymp_indel <- subset(lymp_indel, Func.refGene == "exonic" &
                             ExonicFunc.refGene != "." & ExonicFunc.refGene != "unknown")
      lymp_indel_filtered1 <- subset(lymp_indel, T_AF >= 0.03)
      lymp_indel_filtered2 <- subset(lymp_indel_filtered1, T_AF >= 0.05)
      table_indel$lymp_all[i] = nrow(lymp_indel)
      table_indel$lymp_f1[i] = nrow(lymp_indel_filtered1)
      table_indel$lymp_f2[i] = nrow(lymp_indel_filtered2)
      write.table(lymp_snp, paste("../Results/mutect/lymp_bedonly/", sample$lymp[i], ".snp.dedup.vcf", sep = ""),
                  quote = F, row.names = F, sep = "\t")
      write.table(lymp_indel, paste("../Results/mutect/lymp_bedonly/", sample$lymp[i], ".indel.dedup.vcf", sep = ""),
                  quote = F, row.names = F, sep = "\t")
      write.table(lymp_snp_filtered2, paste("../Results/mutect/lymp_bedonly/", sample$lymp[i], ".snp.dedup.filtered.vcf", sep = ""),
                  quote = F, row.names = F, sep = "\t")
      write.table(lymp_indel_filtered2, paste("../Results/mutect/lymp_bedonly/", sample$lymp[i], ".indel.dedup.filtered.vcf", sep = ""),
                  quote = F, row.names = F, sep = "\t")
    }
  }
  
  if(file.exists(paste("../Results/mutect/lymp_bedonly/", sample$lymp[i], ".dedup.vcf", sep = "")) &
     file.exists(paste("../Results/mutect/tumor_bedonly/", sample$tumor[i], ".dedup.vcf", sep = ""))){
    table_snp$overlap_all[i] = length(intersect(tumor_snp$SNV, lymp_snp$SNV))
    table_snp$overlap_filtered1[i] = length(intersect(tumor_snp_filtered1$SNV, lymp_snp_filtered1$SNV))
    table_snp$overlap_filtered2[i] = length(intersect(tumor_snp_filtered2$SNV, lymp_snp_filtered2$SNV))
    table_indel$overlap_all[i] = length(intersect(tumor_indel$SNV, lymp_indel$SNV))
    table_indel$overlap_filtered1[i] = length(intersect(tumor_indel_filtered1$SNV, lymp_indel_filtered1$SNV))
    table_indel$overlap_filtered2[i] = length(intersect(tumor_indel_filtered2$SNV, lymp_indel_filtered2$SNV))
  }
}

write.table(table_snp, "../Analysis/Overlap_snp_mutect.txt", quote = F, row.names = F, sep = "\t")
write.table(table_indel, "../Analysis/Overlap_indel_mutect.txt", quote = F, row.names = F, sep = "\t")




