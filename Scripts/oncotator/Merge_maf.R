lymp_list <- list.files(path = "../Results/mutect/lymp_oncotator", pattern = "maf$")

for(i in 1:length(lymp_list)){
  maf <- read.delim(paste("../Results/mutect/lymp_oncotator", 
                          lymp_list[i], sep = "/"), header = T, comment.char = "#", 
                    stringsAsFactors = F, skipNul = TRUE)
  if(i == 1){
    all <- maf
  }else{
    all <- rbind(all, maf)
  }
}
all$Tumor_Sample_Barcode <- all$Matched_Norm_Sample_Barcode

write.table(all, "../Analysis/lymp.maf", row.names = F, quote = F, sep = "\t")

tumor_list <- list.files(path = "../Results/mutect/tumor_oncotator", pattern = "maf$")

for(i in 1:length(tumor_list)){
  maf <- read.delim(paste("../Results/mutect/tumor_oncotator", 
                          tumor_list[i], sep = "/"), header = T, comment.char = "#", 
                    stringsAsFactors = F, skipNul = TRUE)
  if(i == 1){
    all <- maf
  }else{
    all <- rbind(all, maf)
  }
}
all$Tumor_Sample_Barcode <- all$Matched_Norm_Sample_Barcode
write.table(all, "../Analysis/tumor.maf", row.names = F, quote = F, sep = "\t")


lymp_list <- list.files(path = "../Results/mutect/lymp_oncotator_all", pattern = "maf$")

for(i in 1:length(lymp_list)){
  maf <- read.delim(paste("../Results/mutect/lymp_oncotator_all", 
                          lymp_list[i], sep = "/"), header = T, comment.char = "#", 
                    stringsAsFactors = F, skipNul = TRUE)
  if(i == 1){
    all <- maf
  }else{
    all <- rbind(all, maf)
  }
}
all$Tumor_Sample_Barcode <- all$Matched_Norm_Sample_Barcode

write.table(all, "../Analysis/lymp_all.maf", row.names = F, quote = F, sep = "\t")

tumor_list <- list.files(path = "../Results/mutect/tumor_oncotator_all", pattern = "maf$")

for(i in 1:length(tumor_list)){
  maf <- read.delim(paste("../Results/mutect/tumor_oncotator_all", 
                          tumor_list[i], sep = "/"), header = T, comment.char = "#", 
                    stringsAsFactors = F, skipNul = TRUE)
  if(i == 1){
    all <- maf
  }else{
    all <- rbind(all, maf)
  }
}
all$Tumor_Sample_Barcode <- all$Matched_Norm_Sample_Barcode
write.table(all, "../Analysis/tumor_all.maf", row.names = F, quote = F, sep = "\t")


tumor_list <- list.files(path = "../Results/mutect/tumor_oncotator", pattern = "maf$")

for(i in 1:length(tumor_list)){
  maf <- read.delim(paste("../Results/mutect/tumor_oncotator", 
                          tumor_list[i], sep = "/"), header = T, comment.char = "#", 
                    stringsAsFactors = F, skipNul = TRUE)
  if(i == 1){
    all <- maf
  }else{
    all <- rbind(all, maf)
  }
}
all$Tumor_Sample_Barcode <- all$Matched_Norm_Sample_Barcode
write.table(all, "../Analysis/tumor.maf", row.names = F, quote = F, sep = "\t")


lymp_list <- list.files(path = "../Results/mutect/lymp_oncotator_001", pattern = "maf$")

for(i in 1:length(lymp_list)){
  maf <- read.delim(paste("../Results/mutect/lymp_oncotator_001", 
                          lymp_list[i], sep = "/"), header = T, comment.char = "#", 
                    stringsAsFactors = F, skipNul = TRUE)
  if(i == 1){
    all <- maf
  }else{
    all <- rbind(all, maf)
  }
}
all$Tumor_Sample_Barcode <- all$Matched_Norm_Sample_Barcode

write.table(all, "../Analysis/lymp_001.maf", row.names = F, quote = F, sep = "\t")

tumor_list <- list.files(path = "../Results/mutect/tumor_oncotator_001", pattern = "maf$")

for(i in 1:length(tumor_list)){
  maf <- read.delim(paste("../Results/mutect/tumor_oncotator_001", 
                          tumor_list[i], sep = "/"), header = T, comment.char = "#", 
                    stringsAsFactors = F, skipNul = TRUE)
  if(i == 1){
    all <- maf
  }else{
    all <- rbind(all, maf)
  }
}
all$Tumor_Sample_Barcode <- all$Matched_Norm_Sample_Barcode
write.table(all, "../Analysis/tumor_001.maf", row.names = F, quote = F, sep = "\t")
