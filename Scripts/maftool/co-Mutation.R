createOncoMatrix = function(m, g = NULL, chatty = TRUE){
  
  if(is.null(g)){
    stop("Please provde atleast two genes!")
  }
  
  subMaf = subsetMaf(maf = m, genes = g, includeSyn = FALSE)
  
  if(nrow(subMaf) == 0){
    return(NULL)
  }
  
  oncomat = data.table::dcast(data = subMaf[,.(Hugo_Symbol, Variant_Classification, Tumor_Sample_Barcode)], formula = Hugo_Symbol ~ Tumor_Sample_Barcode,
                              fun.aggregate = function(x){
                                x = unique(as.character(x))
                                xad = x[x %in% c('Amp', 'Del')]
                                xvc = x[!x %in% c('Amp', 'Del')]
                                
                                if(length(xvc)>0){
                                  xvc = ifelse(test = length(xvc) > 1, yes = 'Multi_Hit', no = xvc)
                                }
                                
                                x = ifelse(test = length(xad) > 0, yes = paste(xad, xvc, sep = ';'), no = xvc)
                                x = gsub(pattern = ';$', replacement = '', x = x)
                                x = gsub(pattern = '^;', replacement = '', x = x)
                                return(x)
                              } , value.var = 'Variant_Classification', fill = '', drop = FALSE)
  
  #convert to matrix
  data.table::setDF(oncomat)
  rownames(oncomat) = oncomat$Hugo_Symbol
  oncomat = as.matrix(oncomat[,-1, drop = FALSE])
  
  variant.classes = as.character(unique(subMaf[,Variant_Classification]))
  variant.classes = c('',variant.classes, 'Multi_Hit')
  names(variant.classes) = 0:(length(variant.classes)-1)
  
  #Complex variant classes will be assigned a single integer.
  vc.onc = unique(unlist(apply(oncomat, 2, unique)))
  vc.onc = vc.onc[!vc.onc %in% names(variant.classes)]
  names(vc.onc) = rep(as.character(as.numeric(names(variant.classes)[length(variant.classes)])+1), length(vc.onc))
  variant.classes2 = c(variant.classes, vc.onc)
  
  oncomat.copy <- oncomat
  #Make a numeric coded matrix
  for(i in 1:length(variant.classes2)){
    oncomat[oncomat == variant.classes2[i]] = names(variant.classes2)[i]
  }
  
  #If maf has only one gene
  if(nrow(oncomat) == 1){
    mdf  = t(matrix(as.numeric(oncomat)))
    rownames(mdf) = rownames(oncomat)
    colnames(mdf) = colnames(oncomat)
    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes))
  }
  
  #convert from character to numeric
  mdf = as.matrix(apply(oncomat, 2, function(x) as.numeric(as.character(x))))
  rownames(mdf) = rownames(oncomat.copy)
  
  
  #If MAF file contains a single sample, simple sorting is enuf.
  if(ncol(mdf) == 1){
    mdf = as.matrix(mdf[order(mdf, decreasing = TRUE),])
    colnames(mdf) = sampleId
    
    oncomat.copy = as.matrix(oncomat.copy[rownames(mdf),])
    colnames(oncomat.copy) = sampleId
    
    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes))
  } else{
    #Sort by rows as well columns if >1 samples present in MAF
    #Add total variants per gene
    mdf = cbind(mdf, variants = apply(mdf, 1, function(x) {
      length(x[x != "0"])
    }))
    #Sort by total variants
    mdf = mdf[order(mdf[, ncol(mdf)], decreasing = TRUE), ]
    #colnames(mdf) = gsub(pattern = "^X", replacement = "", colnames(mdf))
    nMut = mdf[, ncol(mdf)]
    
    mdf = mdf[, -ncol(mdf)]
    
    mdf.temp.copy = mdf #temp copy of original unsorted numeric coded matrix
    
    mdf[mdf != 0] = 1 #replacing all non-zero integers with 1 improves sorting (& grouping)
    tmdf = t(mdf) #transposematrix
    mdf = t(tmdf[do.call(order, c(as.list(as.data.frame(tmdf)), decreasing = TRUE)), ]) #sort
    
    mdf.temp.copy = mdf.temp.copy[rownames(mdf),] #organise original matrix into sorted matrix
    mdf.temp.copy = mdf.temp.copy[,colnames(mdf)]
    mdf = mdf.temp.copy
    
    #organise original character matrix into sorted matrix
    oncomat.copy <- oncomat.copy[,colnames(mdf)]
    oncomat.copy <- oncomat.copy[rownames(mdf),]
    
    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes))
  }
}


createmutMat = function(maf, top = 25, genes = NULL, samples = NULL){
  # 找到出现人群突变频率最高的基因
  if(is.null(genes)){
    genes = getGeneSummary(x = maf)[1:top, Hugo_Symbol]
  }
  if(is.null(samples)){
    samples = getGeneSummary(x = maf)[1:top, Tumor_Sample_Barcode]
  }
  # 至少需要2个基因
  if(length(genes) < 2){
    stop("Minimum two genes required!")
  }
  
  om = createOncoMatrix(m = maf, g = genes)
  # patient列表
  all.tsbs = as.character(getSampleSummary(x = maf)[,Tumor_Sample_Barcode])
  # patient为行，基因为列，统计出现次数  
  mutMat = t(om$numericMatrix)
  # 无top突变基因的病人列表
  missing.tsbs = samples[!samples %in% rownames(mutMat)]
  # 至少需要两个病人  
  if(nrow(mutMat) < 2){
    stop("Minimum two samples required!")
  }
  # 将出现次数转化为0-1变量，出现为1，为出现为0
  mutMat[mutMat > 0 ] = 1
  # 将未包括在mutMat中的病人（无top突变的）添加至mutMat中
  if(length(missing.tsbs) > 0){
    missing.tsbs = as.data.frame(matrix(data = 0, nrow = length(missing.tsbs), ncol = ncol(mutMat)),
                                 row.names = missing.tsbs)
    colnames(missing.tsbs) = colnames(mutMat)
    mutMat = rbind(mutMat, missing.tsbs)
  }
  # 以上输入文件仅要求提供人群突变频率最高的基因上，每个病人是否有该基因的突变（0-1）
  # maftools默认的方法有一个缺陷，如果人群中，某一样本无突变，则该样本会被忽略，为避免该情况使用samples变量作为样本列表
  return(mutMat)
}

somaticInteractions_mutMat = function(mutMat, top = 25, genes = NULL, pvalue = c(0.05, 0.01), 
                                      returnAll = FALSE, findPathways = TRUE, kMax = 3, 
                                      fontSize = 0.8, verbose = TRUE){
  #pairwise fisher test source code borrowed from: https://www.nature.com/articles/ncomms6901
  
  if(is.null(genes)){
    genes = colnames(mutMat)   
  }
  mutMat <- mutMat[, genes]
  if(top <= ncol(mutMat)){
    mutMat <- mutMat[,names(sort(apply(mutMat,2,sum), decreasing = T))[1:top]]
  }else{
    stop(paste("less than", top, "genes were selected from input matrix!"))
  }
  
  if(ncol(mutMat) < 2){
    stop("Minimum two genes required!")
  }
  if(nrow(mutMat) < 2){
    stop("Minimum two samples required!")
  }
  interactions = sapply(1:ncol(mutMat), function(i) sapply(1:ncol(mutMat), function(j) {f<- try(fisher.test(mutMat[,i], mutMat[,j]), silent=TRUE); if(class(f)=="try-error") NA else ifelse(f$estimate>1, -log10(f$p.val),log10(f$p.val))} ))
  oddsRatio <- oddsGenes <- sapply(1:ncol(mutMat), function(i) sapply(1:ncol(mutMat), function(j) {f<- try(fisher.test(mutMat[,i], mutMat[,j]), silent=TRUE); if(class(f)=="try-error") f=NA else f$estimate} ))
  rownames(interactions) = colnames(interactions) = rownames(oddsRatio) = colnames(oddsRatio) = colnames(mutMat)
  
  if(returnAll){
    sigPairs = which(x = 10^-abs(interactions) < 1, arr.ind = TRUE)
  }else{
    sigPairs = which(x = 10^-abs(interactions) < max(pvalue), arr.ind = TRUE)
  }
  
  sigPairsTbl = data.table::rbindlist(
    lapply(X = seq_along(1:nrow(sigPairs)), function(i) {
      x = sigPairs[i,]
      g1 = rownames(interactions[x[1], x[2], drop = FALSE])
      g2 = colnames(interactions[x[1], x[2], drop = FALSE])
      tbl = as.data.frame(table(apply(X = mutMat[,c(g1, g2), drop = FALSE], 1, paste, collapse = "")))
      combn = data.frame(t(tbl$Freq))
      colnames(combn) = tbl$Var1
      pval = 10^-abs(interactions[x[1], x[2]])
      fest = oddsRatio[x[1], x[2]]
      d = data.table::data.table(gene1 = g1,
                                 gene2 = g2,
                                 pValue = pval, oddsRatio = fest)
      d = cbind(d, combn)
      d
    }), fill = TRUE)
  
  sigPairsTbl = sigPairsTbl[!gene1 == gene2] #Remove doagonal elements
  sigPairsTbl$Event = ifelse(test = sigPairsTbl$oddsRatio > 1, yes = "Co_Occurance", no = "Mutually_Exclusive")
  sigPairsTbl$pair = apply(X = sigPairsTbl[,.(gene1, gene2)], MARGIN = 1, FUN = function(x) paste(sort(unique(x)), collapse = ", "))
  sigPairsTblSig = sigPairsTbl[order(as.numeric(pValue))][!duplicated(pair)]
  
  #Source code borrowed from: https://www.nature.com/articles/ncomms6901
  if(nrow(interactions) >= 5){
    interactions[10^-abs(interactions) > max(pvalue)] = 0
    diag(interactions) <- 0
    m <- nrow(interactions)
    n <- ncol(interactions)
    
    interactions[interactions < -4] = -4
    interactions[interactions > 4] = 4
    r = interactions
    rd = hclust(dist(r))$order
    cd = hclust(dist(t(r)))$order
    interactions = interactions[rd, , drop = FALSE]
    interactions = interactions[,rd, drop = FALSE]
    
    interactions[lower.tri(x = interactions)] = NA
    
    par(bty="n", mgp = c(2,.5,0), mar = c(2, 4, 3, 5)+.1, las=2, tcl=-.33)
    image(x=1:n, y=1:m, interactions, col=RColorBrewer::brewer.pal(9,"PiYG"),
          breaks = c(-4:0-.Machine$double.eps,0:4), xaxt="n", yaxt="n",
          xlab="",ylab="", xlim=c(0, n+4), ylim=c(0, n+4))
    abline(h=0:n+.5, col="white", lwd=.5)
    abline(v=0:n+.5, col="white", lwd=.5)
    
    mtext(side = 2, at = 1:m, text = colnames(interactions), cex = fontSize, font = 2)
    mtext(side = 3, at = 1:n, text = colnames(interactions), las = 2, line = -2, cex = fontSize, font = 2)
    
    #q <- p.adjust(10^-abs(interactions), method="BH")
    #p <- p.adjust(10^-abs(interactions), method="holm")
    #w = arrayInd(which(interactions < .05), rep(m,2))
    #points(w, pch=".", col="white", cex=1.5)
    w = arrayInd(which(10^-abs(interactions) < min(pvalue)), rep(m,2))
    points(w, pch="*", col="black")
    w = arrayInd(which(10^-abs(interactions) < max(pvalue)), rep(m,2))
    points(w, pch=".", col="black")
    #image(y = 1:8 +6, x=rep(n,2)+c(2,2.5)+1, z=matrix(c(1:8), nrow=1), col=brewer.pal(8,"PiYG"), add=TRUE)
    image(y = seq(0.5*nrow(interactions), 0.9*nrow(interactions), length.out = 8), x=rep(n,2)+c(2,2.5)+1, z=matrix(c(1:8), nrow=1), col = RColorBrewer::brewer.pal(8,"PiYG"), add=TRUE)
    #axis(side = 4, at = seq(1,7) + 6.5,  tcl=-.15, label=seq(-3, 3), las=1, lwd=.5)
    atLims = seq(0.5*nrow(interactions), 0.9*nrow(interactions), length.out = 7)
    axis(side = 4, at = atLims,  tcl=-.15, label=seq(-3, 3), las=1, lwd=.5)
    mtext(side=4, at = median(atLims), "log10 (p-value)", las=3, cex = 0.9, line = 3, font = 2)
    
    par(xpd=NA)
    text(x=n+2.2, y= max(atLims)+1.2, "Co-occurance", pos=4, cex = 0.9, font = 2)
    text(x=n+2.2, y = min(atLims)-1.2, "Exclusive", pos=4, cex = 0.9, font = 2)
    
    points(x = n+1, y = 0.2*n, pch = "*", cex = 2)
    text(x = n+1, y = 0.2*n, paste0(" p < ", min(pvalue)), pos=4, cex = 0.9, font = 2)
    points(x = n+1, y = 0.1*n, pch = ".", cex = 2)
    text(x = n+1, y = 0.1*n, paste0("p < ", max(pvalue)), pos=4, cex = 0.9)
  }
  
  
  sig.genes.pvals = NULL
  
  if(findPathways){
    if(nrow(sigPairsTblSig[Event %in% 'Mutually_Exclusive']) > 2){
      if(verbose){
        message("Checking for Gene sets.. ")
      }
      sig.genes = unique(c(sigPairsTblSig[Event %in% 'Mutually_Exclusive', gene1], sigPairsTblSig[Event %in% 'Mutually_Exclusive', gene2]))
      sig.genes.pvals = c()
      
      for(k in 3:kMax){
        sig.genes.combn = combn(x = sig.genes, m = k)
        
        if(verbose){
          message(paste0("k = ", k, ": ", ncol(sig.genes.combn), " combinations.."))
        }
        
        sps = lapply(seq_along(1:ncol(sig.genes.combn)), function(i){
          x = sig.genes.combn[,i]
          mm = mutMat[,x, drop = FALSE]
          grid.mat = t(expand.grid(rep(list(0:1), k)))
          #colllapse grid and get all the levels (all posiible combinations)
          lvls = names(table(apply(grid.mat, 2, paste, collapse = '')))
          mm.lvls = data.frame(table(apply(mm, 1, paste, collapse = '')))
          
          #check if for any missing combinations
          lvls.missing = lvls[!lvls %in% mm.lvls[,1]]
          
          if(length(lvls.missing) > 0){
            mm.lvls = rbind(mm.lvls, data.frame(Var1 = lvls.missing, Freq = 0)) #add missing combinations with zero values
          }
          
          #reorder
          mm.lvls = mm.lvls[order(mm.lvls$Var1),]
          if(verbose){
            message("Geneset: ", paste(x, collapse = ", "))
          }
          xp = cometExactTest::comet_exact_test(tbl = as.integer(as.character(mm.lvls$Freq)), mutmatplot = FALSE, pvalthresh = 0.1)
          data.table::data.table(gene_set = paste(x, collapse = ", "), pvalue = xp)
        })
        
        sig.genes.pvals = rbind(sig.genes.pvals, data.table::rbindlist(sps))
      }
      
      sig.genes.pvals = sig.genes.pvals[pvalue > 0][order(pvalue, decreasing = FALSE)]
      if(nrow(sig.genes.pvals[pvalue < 0.05]) > 0){
        if(verbose){
          message("Signifcantly altered gene-sets:")
          print(sig.genes.pvals[pvalue < 0.05])
        }
      }
      
    }
  }
  
  sigPairsTblSig = sigPairsTblSig[,pair := NULL]
  
  return(list(pairs = sigPairsTblSig, gene_sets = sig.genes.pvals))
}

require(maftools)


mutMat_CNV <- read.delim("mutMat.txt", header = T, stringsAsFactors = F)
somaticInteractions_mutMat(mutMat_CNV, pvalue = 0.1, top = 25)

# same as maf@data$Tumor_Sample_Barcode
samples <- read.delim("../../Summary/paired.info", header = T, stringsAsFactors = F)
samples <- substr(samples$blood, 1, 5)

laml = read.maf(maf = "../../Analysis/tumor_001.maf")
mutMat_SNV <- createmutMat(laml, top = 25, genes = NULL, samples = samples)
somaticInteractions_mutMat(mutMat_SNV, pvalue = 0.1, top = 25)

mutMat <- cbind(mutMat_SNV, mutMat_CNV)
somaticInteractions_mutMat(mutMat, pvalue = 0.1)