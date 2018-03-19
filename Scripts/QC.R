Args <- commandArgs()

#print(Args)

if(length(Args) != 7){
  print("Arguments not correct!")}

if(Args[6] == "-h"){
  print("Usage:")
  print("the following arguments should be provided")
  print("Path to the txt containing patientId")
  print("Name of your project.")
  #print("column number of patientId, tissue, sequencingId, sampleId, indexId, tissueType.")
  print("Example: Rscript QC.R /lustre/rdi/user/hantc/tools/QC/patientId.txt Test")
}else{
  
  path <- dirname(Args[6])
  setwd(path)
  
  sample <- read.delim(Args[6], header = F, stringsAsFactors = F)
  list <- read.delim("/work-a/shared/GeneticTest/list.txt", header = T, stringsAsFactors = F)
  
  lst <- list[0,]
  
  for (i in 1:nrow(sample)){
    lst <- rbind(lst, subset(list, patient == sample$V1[i]))
  }
  QC <- lst[,c(1,4,6:8,13)]
  
  QC$PATH <- NA
  for (i in 1:nrow(QC)){
    test_mappedFileName1 <- paste("/work/shared/Analysis/", QC$sequencingId[i], "/mpileup/summary.txt.gz", sep = "")
    test_mappedFileName2 <- paste("/work-a/shared/Mapping/", QC$sequencingId[i], "/mpileup/summary.txt.gz", sep = "")
    if(file.exists(test_mappedFileName1) == TRUE){
      QC$PATH[i] <- paste("/work/shared/Analysis/", QC$sequencingId[i], sep = "")
    }else{
      if(file.exists(test_mappedFileName2) == TRUE){
        QC$PATH[i] <- paste("/work-a/shared/Mapping/", QC$sequencingId[i], sep = "")
      }
    }
  }
  for (i in which(is.na(QC$PATH) == T)){
    test_mappedFileName1 <- paste("/work/shared/Analysis/", QC$sequencingId[i], "/fqstat/summary.txt.gz", sep = "")
    test_mappedFileName2 <- paste("/work-a/shared/Mapping/", QC$sequencingId[i], "/fqstat/summary.txt.gz", sep = "")
    if(file.exists(test_mappedFileName1) == TRUE){
      QC$PATH[i] <- paste("/work/shared/Analysis/", QC$sequencingId[i], sep = "")
    }else{
      if(file.exists(test_mappedFileName2) == TRUE){
        QC$PATH[i] <- paste("/work-a/shared/Mapping/", QC$sequencingId[i], sep = "")
      }
    }
  }
  
  
  QC$Mreads <- NA
  QC$Mbases <- NA
  QC$GC <- NA
  QC$N <- NA
  QC$Q20 <- NA
  QC$Q30 <- NA
  QC$Q40 <- NA
  QC$totalReads <- NA
  QC$mapped <- NA
  QC$mapped_p <- NA
  QC$paired <- NA
  QC$paired_p <- NA
  QC$singleton <- NA
  QC$singleton_p <- NA
  QC$interchrom <- NA
  QC$interchrom_p <- NA
  QC$dupratio <- NA
  QC$mappedBases <- NA
  QC$targetedBases <- NA
  QC$targetedRatio <- NA
  QC$depth <- NA
  QC$coverage <- NA
  QC$coverage10 <- NA
  QC$coverage100 <- NA
  QC$coverage500 <- NA
  QC$coverage1000 <- NA
  QC$coverage5000 <- NA
  QC$mappedBases_mkdup <- NA
  QC$targetedBases_mkdup <- NA
  QC$targetedRatio_mkdup <- NA
  QC$depth_mkdup <- NA
  QC$coverage_mkdup <- NA
  QC$coverage10_mkdup <- NA
  QC$coverage100_mkdup <- NA
  QC$coverage500_mkdup <- NA
  QC$coverage1000_mkdup <- NA
  QC$coverage5000_mkdup <- NA
  QC$panel <- lst$panel
  QC$note <- lst$desc
  
  for (i in 1:nrow(QC)){
    test_Filefolder <- QC$PATH[i]
    if (is.na(test_Filefolder) == F){
      test_fqstatFileName<- paste(test_Filefolder, "/fqstat/summary.txt.gz", sep = "")
      if (file.exists(test_fqstatFileName) == TRUE){
        fqstat <- read.delim(test_fqstatFileName, header = T, stringsAsFactors = F)
        if(length(which(fqstat$id == paste(QC$sampleId[i], "_", QC$indexId[i], sep = ""))) != 0){QC[i, 8:14] <- fqstat[which(fqstat$id == paste(QC$sampleId[i], "_", QC$indexId[i], sep = "")), 2:8]}
      }
      
      test_mappedFileName<- paste(test_Filefolder, "/mapped/summary.txt.gz", sep = "")
      if (file.exists(test_mappedFileName) == TRUE){
        map <- read.delim(test_mappedFileName, header = T, stringsAsFactors = F)
        if(length(which(map$id == paste(QC$sampleId[i], "_", QC$indexId[i], sep = ""))) != 0){QC[i, 15:23] <- map[which(map$id == paste(QC$sampleId[i], "_", QC$indexId[i], sep = "")), 2:10]}
      }
      
      test_metricsFileName <- paste(test_Filefolder, "/mapped/", QC$sampleId[i], "_", QC$indexId[i], ".sorted.filtered.mkdup.metrics", sep = "")
      if (file.exists(test_metricsFileName) == TRUE){
        QC[i,24] = read.delim(test_metricsFileName, header = F, stringsAsFactors = F)[14,2]
      }
      
      test_coverageFileName <- paste(test_Filefolder, "/mpileup/summary.txt.gz", sep = "")
      if (file.exists(test_coverageFileName) == TRUE){
        coverage <- read.delim(test_coverageFileName, header = T, stringsAsFactors = F)
        if(length(which(coverage$id == paste(QC$sampleId[i], "_", QC$indexId[i], sep = ""))) != 0){QC[i, 25:34] <- coverage[which(coverage$id == paste(QC$sampleId[i], "_", QC$indexId[i], sep = "")), c(2,5:7,9,11,13,15,17,19)]}
      }
      
      test_coveragemkdupFileName <- paste(test_Filefolder, "/mpileup.mkdup/summary.txt.gz", sep = "")
      if (file.exists(test_coveragemkdupFileName) == TRUE){
        coverage_mkdup <- read.delim(test_coveragemkdupFileName, header = T, stringsAsFactors = F)
        if (length(which(coverage_mkdup$id == paste(QC$sampleId[i], "_", QC$indexId[i], sep = ""))) != 0){QC[i, 35:44] <- coverage_mkdup[which(coverage_mkdup$id == paste(QC$sampleId[i], "_", QC$indexId[i], sep = "")), c(2,5:7,9,11,13,15,17,19)]}
      }
    }
    
    
    if (i %% 20 == 0){print(paste(i, "patients have been processed."))}
    
  }
  
  write.table(QC, paste(Args[7], "_QC.txt", sep = ""), row.names = F, quote = F, sep = "\t")
  write.csv(QC, paste(Args[7], "_QC.csv", sep = ""), row.names = F, quote = T, fileEncoding = "GBK")
}
