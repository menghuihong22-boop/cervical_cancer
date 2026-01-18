xmmc <- 'cervical_cancer'
#--------------------------------------------------------------
'****************************************************************'
library(AnnoProbe)
library(GEOquery) 
library(tidyverse)
GSE_file <- 'GSE9750_series_matrix.txt'
GSE_accession <- "GSE9750"
GPL <- 'GPL96'

exp.df=read.table(GSE_file,
                  header=T,
                  sep="\t", #Error in scan(file = file, what = what, sep = sep, quote = quote, dec = dec
                  comment.char="!")
row.names(exp.df) <- exp.df[,1]
head(exp.df)[,1:6]
idprob = AnnoProbe::idmap(GPL,type = 'soft')#
head(idprob)
expr_1 <- filterEM(exp.df,idprob)
expr_1$ID_REF <- row.names(expr_1)
expr_1[1:4,1:4]
exp1 <- expr_1 %>% as_tibble() %>%
  separate_rows("ID_REF", sep = " /// ")#///
colnames(exp1)[1] <- "Symbol"
gsm <- colnames(exp1[,-1])
for( i in 1:length(gsm) ){
  
  exp1[,colnames(exp1) == as.name(gsm[i])] <- 
    as.numeric(unlist(exp1[,colnames(exp1) == as.name(gsm[i])])   )
  
}#
exp1[is.na(exp1)]<-0#na0
exp1[1:4,1:4]
exp2 <- aggregate(.~Symbol,exp1,mean)#
rownames(exp2) <- exp2[,1]
exp2 <- exp2[,-1]
range(exp2)


##file
##nn
##nm
##header
readfile<-function(file,m = 42, n=1000, header=T){
  pt <- file(file, "r")
  name <- NULL
  if(header){
    #name <- strsplit(readLines(pt, 1), split="\t")[[1]];  #
    f1 <- readLines(pt, n)
    data <- read.table(text=f1, sep="\t", 
                       #col.names=name,
                       skip = m)
  }else{
    data <- read.table(text=f1, sep="\t")
  }
  close(pt)
  data 
}
exp.gp <- readfile(file=GSE_file,m = 32,n=35, header=T) %>% t() %>% as.data.frame
exp.gp <- exp.gp[,-3] %>% select(2,1) 
exp.gp <- exp.gp[-1,]
colnames(exp.gp) <- c('sample','group')

write.csv(exp.gp,paste0(GSE_accession,'_group.csv'),row.names = FALSE)


#
'***********************************************************************'
exp.gp <- read.csv(paste0(GSE_accession,'_group.csv'),header = T,row.names = 1)
exp3 <- exp2 %>% select(row.names(exp.gp)) %>% t() %>% as.data.frame()
exp4 <- cbind(exp.gp,exp3) %>% arrange(desc(group))  #arrange(desc(assists))
exp.gp1 <- exp4 %>% select(1)
exp4 <- exp4[,-1] %>% t() %>% as.data.frame
#log2
ex <- exp4
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC) { #ex[which(ex <= 0)] <- NaN
  exp4 <- log2(ex+1)
  print("log2 transform finished")}else{print("log2 transform not needed")}
range(exp4)
write.csv(exp4,paste0(GSE_accession,'_exp.csv'))
write.csv(exp.gp1,paste0(GSE_accession,'_group.csv'))




#DEG -------------------------------------------------------------------
setwd('../01.DEGs/DEGs')
GSE_accession <- "GSE9750"
exp <- read.csv(paste0(GSE_accession,'_exp.csv'),header = T,row.names = 1)
group <- read.csv(paste0(GSE_accession,'_group.csv'),header = T,row.names = 1)


pacman::p_load(limma, DEqMS, DESeq2, edgeR, ggplot2, ggrepel, 
               pheatmap, circlize, ComplexHeatmap, heatmap3,tidyverse) # 


source("E:/A /A //DEG_limma/DEGs_limma_right.R")

logFoldChange <- 0.5 # 
adj.P.Val <- 0.05 # P(FDR)
p_choose <- "P.Value"  ## adj.P.Val  P.Value

DEG_limma <- function(exp, group){
  
  design <- model.matrix(~ 0 + group$group)
  colnames(design) <- c("Case", "Control")
  fit <- lmFit(exp, design)
  contrast.matrix <- makeContrasts(Case - Control, levels = design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  options(digits = 4) # 输出结果的小数点后保留4位
  result_DEG <- topTable(fit2, coef = 1, adjust = "BH", number = dim(exp)[1])
  
  # 所有基因的差异情况
  write.table(data.frame(symbol = rownames(result_DEG), result_DEG), 
              file = paste0(GSE_accession,"_DEG_all_information.csv"), sep = ",", row.names = F)
  
  
  if(! file.exists(paste0("0_", GSE_accession, '_DEGs.log'))){
    
    diffSig_DEG <- result_DEG[(result_DEG$adj.P.Val < 0.05 & abs(result_DEG$logFC) >= 1), ] 
    DEG <- rownames(diffSig_DEG) 
    DEG_up <- rownames(diffSig_DEG[diffSig_DEG$logFC > 0,]) 
    DEG_down <- rownames(diffSig_DEG[diffSig_DEG$logFC < 0,]) 
    cat(paste0(GSE_accession, '( adj.p < 0.05, log|FC| >= 1 ): DEGs: ', length(DEG), "; DEGs_up: ", length(DEG_up), "; DEGs_down: ", length(DEG_down)),
        file = paste0("0_", GSE_accession, '_DEGs.log'), append = T, sep = '\n')
    
    diffSig_DEG <- result_DEG[(result_DEG$P.Value < 0.05 & abs(result_DEG$logFC) >= 1), ] 
    DEG <- rownames(diffSig_DEG)
    DEG_up <- rownames(diffSig_DEG[diffSig_DEG$logFC > 0,])
    DEG_down <- rownames(diffSig_DEG[diffSig_DEG$logFC < 0,])
    cat(paste0(GSE_accession, '( p.val < 0.05, log|FC| >= 1 ): DEGs: ', length(DEG), "; DEGs_up: ", length(DEG_up), "; DEGs_down: ", length(DEG_down)),
        file = paste0("0_", GSE_accession, '_DEGs.log'), append = T, sep = '\n')
    
    diffSig_DEG <- result_DEG[(result_DEG$P.Value < 0.05 & abs(result_DEG$logFC) >= 0.5), ] 
    DEG <- rownames(diffSig_DEG)
    DEG_up <- rownames(diffSig_DEG[diffSig_DEG$logFC > 0,])
    DEG_down <- rownames(diffSig_DEG[diffSig_DEG$logFC < 0,])
    cat(paste0(GSE_accession, '( p.val < 0.05, log|FC| >= 0.5 ): DEGs: ', length(DEG), "; DEGs_up: ", length(DEG_up), "; DEGs_down: ", length(DEG_down)),
        file = paste0("0_", GSE_accession, '_DEGs.log'), append = T, sep = '\n')
  }
  
  # 显著差异的基因
  
  diffSig_DEG <- result_DEG[(result_DEG[,p_choose] < adj.P.Val & abs(result_DEG$logFC) >= logFoldChange), ] 
  diffSig_DEG <- diffSig_DEG[order(diffSig_DEG$logFC),]
  DEG <- rownames(diffSig_DEG)
  DEG_up <- rownames(diffSig_DEG[diffSig_DEG$logFC > 0,])
  DEG_down <- rownames(diffSig_DEG[diffSig_DEG$logFC < 0,])
  write.table(DEG, file = paste0(GSE_accession, "_DEGs_all.txt"),row.names = F, col.names = F, quote=F)
  write.table(DEG_up, file = paste0(GSE_accession,"_DEGs_up.txt"),row.names = F, col.names = F, quote=F)
  write.table(DEG_down, file = paste0(GSE_accession, "_DEGs_down.txt"),row.names = F, col.names = F, quote=F)
  write.table(data.frame(symbol = rownames(diffSig_DEG), diffSig_DEG), 
              file = paste0(GSE_accession, "_DEG_sig_information.csv"), sep = ",", row.names = F)
  cat(paste0("Filter criteria: ", p_choose, " < 0.05, log|FC| >= ", logFoldChange),
      file = paste0("0_", GSE_accession, '_DEGs.log'), append = T, sep = '\n')
  message("Table completed")
  
  # 火山图
  
  DEGs_file <- result_DEG %>% as.data.frame() # 差异基因统计信息
  DEGs_file$color <- ifelse(DEGs_file[,p_choose] < adj.P.Val & abs(DEGs_file$logFC) >= logFoldChange,
                            ifelse(DEGs_file$logFC > logFoldChange, "Up expression", "Down expression"), 
                            "Non significant")
  DEGs_file$color <- factor(DEGs_file$color, levels = c("Up expression", "Down expression", "Non significant"))
  colnames(DEGs_file)[colnames(DEGs_file) == p_choose] <- "p_choose"
  # table(DEGs_file$color)
  label <- c(rownames(diffSig_DEG %>% slice_max(logFC, n = 5)), rownames(diffSig_DEG %>% slice_min(logFC, n = 5)))
  color <- c("Up expression" = "red", "Non significant" = "gray", "Down expression" = "blue")
  shape <- c("Up expression" = 24, "Non significant" = 16, "Down expression" = 25)
  p <- ggplot(DEGs_file, aes(logFC, -log10(p_choose), col = color)) +
    geom_point(aes(shape = color, fill = color), size=2) +
    theme_bw() +
    scale_color_manual(values = color) +
    scale_shape_manual(values = shape) +
    labs(x = paste0("log2(Fold-change Case VS Control)"), y = paste0("-log10(", p_choose, ")")) +
    geom_label_repel(
      data = DEGs_file[label,],
      aes(label = label),
      size = 5,
      # fill = "darkred", color = "white",
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines"),
      show.legend = F
    ) +
    geom_hline(yintercept = -log10(adj.P.Val), lty = 4, col = "#E39A35", lwd = 1) + # lty：设置水平线的类型（直线，虚线），lwd“设置水平线宽度
    geom_vline(xintercept = c(-logFoldChange, logFoldChange), lty = 4, col = "#E39A35", lwd = 1) +
    # annotate("text", x = 2, y = 0, label = "|logFC|>1", size = 5, col = "black") +
    # annotate("text", x = 4, y = 2, label = "adj.p<0.05", size = 5, col = "black") +
    theme(
      # plot.margin = unit(rep(3, 4), "lines"), # plot.margin:图片边缘添加一点空间。边距的顺序是上（top）、右（right）、下（bottom）、左（left）
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      text = element_text(size = 12),
      axis.title.x = element_text(face = "bold", size = 16),
      axis.title.y = element_text(face = "bold", size = 16),
      axis.text.x = element_text(size = 14, angle = 0),
      axis.text.y = element_text(size = 14, angle = 0),
      legend.text = element_text(face = "bold", size = 13),
      legend.title = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(linewidth = 1, colour = "black"),
      panel.border = element_blank(),
      legend.position = "bottom"
    )
  
  pdf(file = paste0(GSE_accession, "_DEGs_volcano.pdf"), height = 8, width = 8, onefile = FALSE)
  print(p)
  dev.off()
  message("Volcano completed")
  
  
  ## 热图
  
  red_de_expr <- exp[DEG, ] 
  ha1 = HeatmapAnnotation(Sample = c(rep("Control", table(group$group)[2]), rep("Case", table(group$group)[1])),
                          col = list(Sample = c("Control" = "#68A180", "Case" = "#E63863")))
  ha2 = rowAnnotation(Gene = c(rep("Down", sum(diffSig_DEG$logFC < 0)), rep("Up", sum(diffSig_DEG$logFC > 0))),
                      col = list(Gene = c("Down" = "#295599", "Up" = "#E95C59")))
  # draw(ha1)
  madt <- as.matrix(red_de_expr)
  madt2 <- t(scale(t(madt)))
  p <- 
    densityHeatmap(madt2,
                   # ylim = c(-2, 2),
                   title = "Distribution as heatmap", 
                   ylab = "Expression",
                   # col = topo.colors(10)
                   col = c('#295599', '#3e94c0', '#78c6d0', '#b4d9e4', '#fffef0', '#f9cdac', '#ec7d92', '#bc448c')
    ) %v%
    Heatmap(madt2, 
            col = colorRampPalette(c("blue", "white", "red"))(100),
            top_annotation = ha1,
            left_annotation = ha2,
            show_row_names = F ,
            show_column_names = F ,
            name = "Expr", 
            cluster_rows = F,
            height = unit(16, "cm"))
  
  pdf(file=paste0(GSE_accession, "_DEGs_pheatmap.pdf"),width=9,height=9,onefile=FALSE)
  print(p)
  dev.off()
  
  top_20_up <- DEGs_file %>% slice_max(n = 20, order_by = logFC)
  top_20_down <- DEGs_file %>% slice_min(n = 20, order_by = logFC)
  top_20 <- rbind(top_20_up, top_20_down)
  top_20 <- rownames(top_20)
  
  red_de_expr <- exp[top_20, ] 
  ha1 = HeatmapAnnotation(Sample = c(rep("Control", table(group$group)[2]), rep("Case", table(group$group)[1])),
                          col = list(Sample = c("Control" = "#68A180", "Case" = "#E63863")))
  ha2 = rowAnnotation(Gene = c(rep("Up", 20), rep("Down", 20)),
                      col = list(Gene = c("Up" = "#295599", "Down" = "#E95C59")))
  # draw(ha1)
  madt <- as.matrix(red_de_expr)
  madt2 <- t(scale(t(madt)))
  p <- 
    densityHeatmap(madt2,
                   # ylim = c(-2, 2),
                   title = "Distribution as heatmap", 
                   ylab = "Expression",
                   # col = topo.colors(10)
                   col = c('#295599', '#3e94c0', '#78c6d0', '#b4d9e4', '#fffef0', '#f9cdac', '#ec7d92', '#bc448c')
    ) %v%
    Heatmap(madt2, 
            col = colorRampPalette(c("blue", "white", "red"))(100),
            top_annotation = ha1,
            left_annotation = ha2,
            show_row_names = T ,
            show_column_names = F ,
            name = "Expr", 
            cluster_rows = F,
            height = unit(16, "cm"))
  
  pdf(file=paste0(GSE_accession, "_DEGs_pheatmap_top20.pdf"),width=9,height=9,onefile=FALSE)
  print(p)
  dev.off()
  
  message("Pheatmap completed")
  
}


DEG_limma(exp, group)




'********************************************************'




'********************************************************'
#SVM-RFE---------------------
library(e1071)
#library(Rmpi)
library(snow)
library(parallel)

#
# source("E:/A /A ///svmRFE.R")

svmRFE <- function(X, k=1, halve.above=5000) {
  # Feature selection with Multiple SVM Recursive Feature Elimination (RFE) algorithm
  n = ncol(X) - 1
  
  # Scale data up front so it doesn't have to be redone each pass
  cat('Scaling data...')
  X[, -1] = scale(X[, -1])
  cat('Done!\n')
  flush.console()
  
  pb = txtProgressBar(1, n, 1, style=3)
  
  i.surviving = 1:n
  i.ranked    = n
  ranked.list = vector(length=n)
  
  # Recurse through all the features
  while(length(i.surviving) > 0) {
    if(k > 1) {
      # Subsample to obtain multiple weights vectors (i.e. mSVM-RFE)            
      folds = rep(1:k, len=nrow(X))[sample(nrow(X))]
      folds = lapply(1:k, function(x) which(folds == x))
      
      # Obtain weights for each training set
      w = lapply(folds, getWeights, X[, c(1, 1+i.surviving)])
      w = do.call(rbind, w)
      
      # Normalize each weights vector
      w = t(apply(w, 1, function(x) x / sqrt(sum(x^2))))
      
      # Compute ranking criteria
      v    = w * w
      vbar = apply(v, 2, mean)
      vsd  = apply(v, 2, sd)
      c    = vbar / vsd
    } else {
      # Only do 1 pass (i.e. regular SVM-RFE)
      w = getWeights(NULL, X[, c(1, 1+i.surviving)])
      c = w * w
    }
    
    # Rank the features
    ranking = sort(c, index.return=T)$ix
    if(length(i.surviving) == 1) {
      ranking = 1
    }
    
    if(length(i.surviving) > halve.above) {
      # Cut features in half until less than halve.above
      nfeat = length(i.surviving)
      ncut  = round(nfeat / 2)
      n     = nfeat - ncut
      
      cat('Features halved from', nfeat, 'to', n, '\n')
      flush.console()
      
      pb = txtProgressBar(1, n, 1, style=3)
      
    } else ncut = 1
    
    # Update feature list
    ranked.list[i.ranked:(i.ranked-ncut+1)] = i.surviving[ranking[1:ncut]]
    i.ranked    = i.ranked - ncut
    i.surviving = i.surviving[-ranking[1:ncut]]
    
    setTxtProgressBar(pb, n-length(i.surviving))
    flush.console()
  }
  
  close(pb)
  
  return (ranked.list)
}

svmRFE.wrap <- function(test.fold, X, ...) {
  # Wrapper to run svmRFE function while omitting a given test fold
  train.data = X[-test.fold, ]
  test.data  = X[test.fold, ]
  
  # Rank the features
  features.ranked = svmRFE(train.data, ...)
  
  return(list(feature.ids=features.ranked, train.data.ids=row.names(train.data), test.data.ids=row.names(test.data)))
}


#
expt <- t(exp) %>% as.data.frame %>% mutate(group = group$group) 
''
improtant_gene <- read.table("clipboard",header=F)
improtant_gene <- improtant_gene$V1
'na0'
improtant_gene_exp <- expt[,improtant_gene] %>% drop_na
'0'
improtant_gene_exp=improtant_gene_exp[,which(colSums(improtant_gene_exp) > 0)]
head(improtant_gene_exp)
write.csv(improtant_gene_exp,'improtant_gene_exp.csv')
group
improtant_gene_exp <-  improtant_gene_exp %>% mutate(group = group$group) %>% 
  select(group,everything())
improtant_gene_exp$group <- as.factor(improtant_gene_exp$group)
#group

input <- improtant_gene_exp
# R4.0factor
#devtools::install_github("Tong-Chen/ImageGP")
library(ImageGP)

# class
group = "group"
# group - 
# group - 
if(numCheck(input[[group]])){
  if (is.numeric(input[[group]])) {
    input[[group]] <- mixedToFloat(input[[group]])
  }
} else{
  input[[group]] <- as.factor(input[[group]])
}
class(input[[group]])
set.seed(2058)

nfold = 10 #10  
nrows = nrow(input)
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds = lapply(1:nfold, function(x) which(folds == x))



results = lapply(folds, svmRFE.wrap, input, k=10, halve.above=100)
#
top.features = WriteFeatures(results, input, save=F)
write.csv(top.features,'SVM_top.features.csv')
#10min15
featsweep = lapply(1:20, FeatSweep.wrap, results, input)

no.info = min(prop.table(table(input[,1])))
errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))

#
pdf("svm_rfe.pdf", height = 8, width = 10)
PlotErrors(errors,no.info=no.info) #
dev.off()
plot(top.features)#



#LASSO------------------------------------------------
library(glmnet)
set.seed(2000)
#glmnetfamilycox
mod <- glmnet(input[,(2:length(colnames(input)))],input[,1],family = "binomial") #
#familygaussian, binomial, poisson, multinomial, cox, mgaussian
#binomiallogisticscox
mod


##Lassocv.glmnet
cvmod <- cv.glmnet(as.matrix(input[,(2:length(colnames(input)))]),
                   as.matrix(input[,1]),family = "binomial")  # matrix
cvmod

plot(mod,label = T,lwd=2)
plot(mod,xvar = "lambda",label = T,lwd=2)
plot(cvmod)
pdf(file = "lasso.Binomial.Deviance.pdf",height = 5,width = 7)
par(mgp = c(2.5,1,0),mar=c(4,5,3,3))
plot(cvmod,xlab='Log Lambda',cex.lab = 1)+
  text(x = log(cvmod$lambda.min)-0.3 ,y = 0.9,
       paste('Lambda.min\n',round(cvmod$lambda.min,3)),cex=1,adj=0.5)+
  text(x = log(cvmod$lambda.1se),y = 1,
       paste('Lambda.lse\n',round(cvmod$lambda.1se,3)),cex=1)
dev.off()

pdf(file = "lasso.voefficients.venalty.pdf",height = 5,width = 7)
#par(mgp = c(4,1,0),mai=c(2,2,1,1))
plot(mod, xvar="lambda",cex.lab = 1)+
  abline(v = c(log(cvmod$lambda.min), log(cvmod$lambda.1se)),lty=2)+
  text(x = log(cvmod$lambda.min),y =-0.5,
       paste('Lambda.min\n',round(cvmod$lambda.min,3)),cex=1,adj=0.9)+
  text(x = log(cvmod$lambda.1se),y = -0.5,
       paste('Lambda.lse\n',round(cvmod$lambda.1se,3)),cex=1,adj=0.9)
dev.off()
#cvmod$lambda.min cvmod$lambda.1se
cvmod$lambda.min
##coefgenelambda.min
##
coef.min <- coef(cvmod,s="lambda.min")
coef.min
coef.min_data <- as.matrix(coef.min) %>% as.data.frame()
write.csv(coef.min_data,'LASSO_result.csv')
#

#



#Boruta-------------------------------------------
#install.packages("Boruta")
library(Boruta)
library(dplyr)
library(ImageGP)
set.seed(500)

boruta <- Boruta(input[,(2:length(colnames(input)))],input[,1], pValue=0.05, mcAdj=T, 
                 maxRuns=300)

boruta
table(boruta$finalDecision)
boruta.imp <- function(x){
  imp <- reshape2::melt(x$ImpHistory, na.rm=T)[,-1]
  colnames(imp) <- c("Variable","Importance")
  imp <- imp[is.finite(imp$Importance),]
  
  variableGrp <- data.frame(Variable=names(x$finalDecision), 
                            finalDecision=x$finalDecision)
  
  showGrp <- data.frame(Variable=c("shadowMax", "shadowMean", "shadowMin"),
                        finalDecision=c("shadowMax", "shadowMean", "shadowMin"))
  
  variableGrp <- rbind(variableGrp, showGrp)
  
  boruta.variable.imp <- merge(imp, variableGrp, all.x=T)
  
  sortedVariable <- boruta.variable.imp %>% group_by(Variable) %>% 
    summarise(median=median(Importance)) %>% arrange(median)
  sortedVariable <- as.vector(sortedVariable$Variable)
  
  
  boruta.variable.imp$Variable <- factor(boruta.variable.imp$Variable, levels=sortedVariable)
  
  invisible(boruta.variable.imp)
}

boruta.variable.imp <- boruta.imp(boruta)
#
head(boruta.variable.imp)
write.csv(boruta.variable.imp,'Boruta_result.csv')
#
pdf('Boruta_importance_score.pdf',width = 8 ,height = 6)
sp_boxplot(boruta.variable.imp, melted=T, xvariable = "Variable", yvariable = "Importance",
           legend_variable = "finalDecision", 
           legend_variable_order = c("shadowMax", "shadowMean", "shadowMin", "Confirmed","Tentative"),
           xtics_angle = 90)
dev.off()
boruta.finalVarsWithTentative <- 
  data.frame(Item=getSelectedAttributes(boruta, withTentative = T), 
             Type="Boruta_with_tentative")
pdf('Boruta_importance_box.pdf',width = 8 ,height = 6)
caret::featurePlot(input[,boruta.finalVarsWithTentative$Item], 
                   input[,1], plot="box")
dev.off()




#KM-------------------------
'***********************************************************'
dir <- '06.TCGA'
if (!dir.exists(paste0("../",dir,"/"))) {dir.create(paste0("../",dir,"/"))}
setwd(paste0("../",dir,"/"))
#install.packages('survminer')
library(survminer)
library(survival)
library(tidyverse)
library(clusterProfiler) # R
library(stringr) # 
library(AnnotationDbi)
library(org.Hs.eg.db) # hg19
library(R.utils)
setwd('E:/A /A /cervical_cancer/06.TCGA')

##B
TCGA_exp=read.table(file="TCGA-CESC.htseq_fpkm.tsv",sep = "\t",header = T)
head(TCGA_exp)[1:6,1:5]
TCGA_exp1 <- separate(data = TCGA_exp , 
                      col = Ensembl_ID, into = c("Ensembl_ID", "xxx"), 
                      sep = "[.]") %>% dplyr::select(everything(),-("xxx"))
colnames(TCGA_exp1)
head(TCGA_exp1)[1:6,1:5]
gene_tcga <- bitr(TCGA_exp1$Ensembl_ID, 
                fromType = "ENSEMBL", #fromTypeID
                toType = c("SYMBOL"), #toTypeID
                OrgDb = org.Hs.eg.db)#Orgdb
head(gene_tcga)
TCGA_exp2 = data.frame(gene_tcga,
                       TCGA_exp1[match(gene_tcga$ENSEMBL,TCGA_exp1$Ensembl_ID),]) %>%#
  dplyr::select(SYMBOL,everything() ) %>% 
  #distinct(ENTREZID, .keep_all =TRUE) %>% 
  drop_na()  
head(TCGA_exp2)[1:6,1:5]
TCGA_exp2 <- TCGA_exp2[,-c(2:3)]
head(TCGA_exp2)[1:6,1:5]
dim(TCGA_exp2)
colnames(TCGA_exp2)
class(TCGA_exp2$TCGA.DS.A7WI.01A)
table(TCGA_exp2$SYMBOL)
TCGA_exp2 <- aggregate(.~SYMBOL,TCGA_exp2,mean)#
rownames(TCGA_exp2)=TCGA_exp2[,1]
TCGA_exp2=TCGA_exp2[,-1]
head(TCGA_exp2)[1:6,1:5]
#TCGA
TCGA_exp3 <- TCGA_exp2 %>% t() %>% as.data.frame() %>% 
  mutate(sample = row.names(.)) %>% 
  dplyr::select(sample,everything()) %>% 
  separate(sample,sep = '[.]',
           into = c("o", "oo",'ooo','oooo') ) %>% 
 # mutate(group1 =  str_remove(TCGA_exp3[,4],pattern = c("A"))) %>% 
  mutate(group1 =  as.numeric(str_sub(.[,4],1,nchar(.[,4])-1))) %>% 
  dplyr::select(group1,everything(),-o,-oo,-ooo,-oooo) %>% 
  filter(group1 < 10) %>% 
  dplyr::select(everything(),-group1) %>% 
  t() %>% as.data.frame
##
clinical=read.table(file="TCGA-CESC.survival.tsv",sep = "\t",header = T)
#Clinical=Clinical[,c(1:8,10:13)]
##-.
t=clinical

T1=data.frame(gsub("-",".",t$sample),t)

colnames(T1)[1]="SampleID"
T1
T1=T1[,-c(2,4)]
Clinical=T1


key_gene_00
nian <- 15
for (i in 1:length(key_gene_00)){
  scjy = key_gene_00[i]
  GENE=t(TCGA_exp3[which(row.names(TCGA_exp3)==scjy),]) %>% as.data.frame()
  #B
  
  med=median(GENE[,scjy])
  B = data.frame()
  GENE$type <-  ifelse(GENE[,scjy ] > med ,'High',"Low")
  table(GENE$type)
  
  ##B
  ##,
  gene_clinical=data.frame(Clinical[match(rownames(GENE),Clinical$SampleID),],GENE) %>% 
    mutate(OS.time = .[,3]/365) %>% 
    filter(OS.time < nian)
  
  ##
  ##B12
  gene_clinical$type<- ifelse((gene_clinical$type =="High"), 1,2)
  
  GC=gene_clinical
  
  ##
  fit.surv <-Surv(GC$OS.time,GC$OS)
  
  km<-survfit(fit.surv~1,data = GC)
  
  km_2<- survfit(fit.surv~type,data=GC)
  print(km_2$logse) 
  if(surv_pvalue( km_2, method = "Log-rank")$pval < 0.05){# Gehan-Breslow  Log-rank TW PP
    p <- ggsurvplot(km_2, 
               title= paste0(scjy,"_Kaplan-Meier"),#main
               #palette = 'jco',
               pval=TRUE,  #P
               pval.method = TRUE,
               conf.int = F, # 
               #pval = paste0("p: ",surv_pvalue( km_2, method = "Gehan-Breslow")$pval) , # P
               surv.median.line = "none",  #  hv  h  v
               risk.table = TRUE, # 
               risk.table.col = "strata", # 
               xlab = "Follow up time(y)", # x
               legend = c(0.9,0.2), # 
               #legend.title = , # 
               #legend = "top", #  "bottom"
               xlim = c(0,nian),
               legend.title =paste0(scjy)  ,#
               legend.labs = c("High", "Low"),
               linetype = "strata",# 
               palette = c("#E7B800", "#2E9FDF"),##
               #xlab="Days",
               #ylab="OS"
               ggtheme = theme_bw()
    ) 
   
    pdf(paste0(scjy,"_Survival.pdf"),width = 6,height = 5 )
    print(p)
    #withTimeout(p, timeout = 5)#stop execution after one second
    dev.off()
    message('PDF_done')
    png(paste0(scjy,"_Survival.png"),width = 530,height = 450 )
    print(p)
    #withTimeout(p, timeout = 5)#stop execution after one second
    dev.off()
    message('PNG_done')
  }else{
    message(paste0(scjy,"_Survival not significant") )
  }
 
}





#-Nomogram----------------------------------------------------------
'****************************************************************'
dir <- '07.Nomogram'
setwd(paste0("E:/A /A /",xmmc))
if (!dir.exists(dir)) {dir.create(dir)}
setwd(dir)
range(exp)
#exp <- log2(exp+1)
key_gene_01 <- c('VEGFA','CALML3')
RF.v <- exp[key_gene_01,] %>% t() %>% as.data.frame()
table(group$group)
group.v <- c(rep("Control",table(group$group)[[2]]),rep('Case',table(group$group)[[1]])) %>% 
  factor(.,levels = c("Control","Case"),ordered = F) %>% as.data.frame()
design.v <- model.matrix(~ 0 + group.v$.)
colnames(design.v) <- c( "Control","Case")
#design.v <- model.matrix(~factor(group.v))    #designDEG
design.v
RF.v<-cbind(design.v[,2],RF.v)
colnames(RF.v)[1]<- 'group'
rownames(RF.v) <- rownames(group)
write.table(RF.v,"Nomogram_input.txt",sep="\t",row.names = T,col.names = T)

pbc<-read.table("Nomogram_input.txt",sep="\t",header=TRUE,row.names=1)
#pbc$group <- factor(pbc$group,levels=c(0,1),labels=c('No','Yes'))
library(rms)
dd <- datadist(pbc)
options(datadist = "dd")
LR <- lrm(group~.,pbc,x=T,y=T,maxit=1000)
#LR <- lrm(group~.,pbc,x=F,y=F)
print(LR)
#
x<-nomogram(LR,fun = plogis,lp=F,
            fun.at = c(0.001,0.5,0.999),
            funlabel="Risk of PE")

pdf(file = "fig4-1A.Gene_Nomogram.pdf",width=10,height=6,onefile=FALSE)
plot(x)
dev.off()

#
set.seed(44)
cal1<-calibrate(LR, method="boot",B=40)
pdf(file = "fig4-1B.Nomogram_cal.pdf",width=6,height=6,onefile=FALSE)
plot(cal1,lwd = 2,lty = 1,xlim = c(0,1),ylim= c(0,1),xlab = "Nomogram-prediced OS (%)",
     ylab = "Observed OS (%)",cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
dev.off()

#DCA
#install.packages("rmda")
library(rmda)
# install_github("mdbrown/DecisionCurve")
paste(key_gene_01,collapse = "+")

Nomogram <- decision_curve(group ~ VEGFA+CALML3, data = pbc
                           #,policy = "opt-in"
                           ,study.design = 'cohort')

VEGFA <- decision_curve(group ~ VEGFA, data = pbc
                       #,policy = "opt-in"
                       ,study.design = 'cohort')
CALML3 <- decision_curve(group ~ CALML3, data = pbc
                        #,policy = "opt-in"
                        ,study.design = 'cohort')
paste(key_gene_01,collapse = ",")
list <- list(Nomogram, VEGFA,CALML3)

pdf("DCA_.pdf",width=5,height=5)
par(mgp = c(2.5, 1, 0), mar = c(4, 4.5, 1, 1))
plot_decision_curve(list,
                    curve.names= c("Nomogram",'VEGFA', "CALML3"),
                    cost.benefit.axis =FALSE,
                    # col= c('red', 'blue', 'green'),
                    confidence.intervals=FALSE,
                    standardize = FALSE,
                    legend.position="bottomleft") #
dev.off()

##ROC
pred_f_training<-predict(LR,pbc)
#Death
modelroc <- roc(pbc$group,pred_f_training)
#ROC
pdf(file=paste(GSE_accession, "_ROC_Nomogram.pdf", sep = ""),width=6,height=6) 
plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     print.thres=TRUE)
dev.off()





#ssGSEA####--------------------
'**********************ssGSEA**************************'
# rm(list = ls()); gc()
dir <- '08.ssGSEA'
setwd(paste0("E:/A /A /",xmmc))
if (!dir.exists(dir)) {dir.create(dir)}
setwd(dir)
library(magrittr)
library(stringr)
library(dplyr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(msigdbr)
library(ggrepel)
library(enrichplot)
library(aplot)
# df.deg <- read.table("../01_DEG/01.DEG_all_res.xls", header = T)

#

df.expr <- exp
all.gene <- AnnotationDbi::select(org.Hs.eg.db, rownames(df.expr), "ENTREZID", "SYMBOL")
exprSet <- merge(df.expr, all.gene, by.x = "row.names", by.y  = "SYMBOL") %>% 
  dplyr::distinct(ENTREZID, .keep_all = T) %>% na.omit()
rownames(exprSet) <- NULL
exprSet <- exprSet %>% tibble::column_to_rownames(var = "ENTREZID") %>% dplyr::select(-Row.names)


#
#fea_gene <- c('ACTA1','KYNU')
hub_gene <- key_gene
# hub_gene <- hub_gene$Name
hub_gene = AnnotationDbi::select(org.Hs.eg.db, hub_gene, "ENTREZID", "SYMBOL")
entrez.gene <- hub_gene$ENTREZID
df.gene <- exprSet[entrez.gene, ] %>% as.data.frame()
df.gene2 <- cbind(GeneID = hub_gene$SYMBOL, df.gene)
write.table(df.gene2, file = "01.hub_gene_expr.xls", sep = "\t", quote = F, col.names = T, row.names = F)



gsea.plot = function(res.kegg, top.hall, gene){
  gsdata <- do.call(rbind, lapply(top.hall, enrichplot:::gsInfo, object = res.kegg))
  gsdata$Description = factor(gsdata$Description, levels = top.hall)
  p1 = ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) + theme_classic(14) + 
    theme(panel.grid.major = element_line(colour = "grey92"), 
          panel.grid.minor = element_line(colour = "grey92"), 
          panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) + 
    scale_x_continuous(expand = c(0, 0)) +
    scale_color_brewer(palette = "Set1") +
    ggtitle(gene) +
    geom_hline(yintercept = 0, color = "black", size = 0.8) +
    geom_line(aes_(y = ~runningScore, color = ~Description), size = 1) +
    theme(legend.position = "right", legend.title = element_blank(), legend.background = element_rect(fill = "transparent")) +
    ylab("Running Enrichment Score") + 
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.line.x = element_blank(),
          text = element_text(face = "bold", family = "Times"))
  i = 0
  for (term in unique(gsdata$Description)) {
    idx <- which(gsdata$ymin != 0 & gsdata$Description == term)
    gsdata[idx, "ymin"] <- i
    gsdata[idx, "ymax"] <- i + 1
    i <- i + 1 }
  p2 = ggplot(gsdata, aes_(x = ~x)) + geom_linerange(aes_(ymin = ~ymin, ymax = ~ymax, color = ~Description)) + xlab(NULL) + ylab(NULL) + 
    theme_classic(14) + theme(legend.position = "none", 
                              axis.ticks = element_blank(), 
                              axis.text = element_blank(), 
                              axis.line.x = element_blank(),
                              text = element_text(face = "bold", family = "Times")) + 
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) + 
    scale_color_brewer(palette = "Set1")
  p = aplot::insert_bottom(p1, p2, height = 0.15)
  return(p)
}

df.m = msigdbr()
df.kegg = subset(df.m, gs_subcat == "CP:KEGG")[c(3,5)]
df.exp <- exprSet %>% t %>% as.data.frame()


if (!dir.exists("02.GSEA")) {dir.create("02.GSEA")}
if (!dir.exists("03.GSEA_Res")) {dir.create("03.GSEA_Res")}
lapply(1:nrow(hub_gene), function(i){
  hub = rownames(df.gene2)[i] %>% as.character()
  hub.exp = df.exp[[hub]]
  hub.cor = cor(df.exp, hub.exp, method = "spearman") %>% as.data.frame %>% na.omit
  hub.coreff = hub.cor[[1]]
  names(hub.coreff) = rownames(hub.cor)
  hub.coreff = hub.coreff[order(hub.coreff, decreasing = T)]
  res.kegg = GSEA(hub.coreff, TERM2GENE = df.kegg, pvalueCutoff = 0.2, seed = 1, pAdjustMethod = "BH", eps = 0)
  write.table(res.kegg, file = paste0("03.GSEA_Res/0", i, ".", hub_gene$SYMBOL[which(hub_gene$ENTREZID==hub)], ".res.xls"),
              sep = "\t", row.names = T, col.names = T, quote = F)
  top.kegg = res.kegg@result
  top.kegg = top.kegg[order(top.kegg$p.adjust, decreasing = F),]$Description[1:5]
  p = gsea.plot(res.kegg, top.kegg, hub_gene$SYMBOL[which(hub_gene$ENTREZID==hub)])
  fn1 = paste0("02.GSEA/", sprintf("%02d",i), ".", hub_gene$SYMBOL[which(hub_gene$ENTREZID==hub)], ".png")
  fn2 = paste0("02.GSEA/", sprintf("%02d",i), ".", hub_gene$SYMBOL[which(hub_gene$ENTREZID==hub)], ".pdf")
  ggsave(fn1, p, width = 12, height = 6, units = "in", limitsize = 300)
  ggsave(fn2, p, width = 12, height = 6, units = "in", limitsize = 300)
  return(0)
})




#Immunity--------------------------------
'******************************************************'
dir <- '09.Immunity'
setwd(paste0("E:/A /A /",xmmc))
if (!dir.exists(dir)) {dir.create(dir)}
setwd(dir)
#expgroup
#
require(cowplot)
require(tidyverse)
require(ggplot2)
require(ggsci)
require(ggpubr)

#R
''
key_gene_00 <- read.table("clipboard",header=F)
key_gene_00 <- key_gene_00$V1
key_gene <- key_gene_00
'na0'
key_gene_00_exp <- expt[,key_gene_00] %>% drop_na
'0'
key_gene_00_exp=key_gene_00_exp[,which(colSums(key_gene_00_exp) > 0)]
head(key_gene_00_exp)
write.csv(key_gene_00_exp,'key_gene_00_exp.csv')
group
key_gene_00_exp <-  key_gene_00_exp %>% mutate(group = group$group) %>% 
  select(group,everything())
key_gene_00_exp$group <- as.factor(key_gene_00_exp$group)
rm(c)
c = as.character()
for (i in 1:length(key_gene)) {
  a <- key_gene[i]
  if (i<2){c = paste0(c,a)}else{c = paste0(c,',',a)}
}
c#



mydata<-expt %>% 
  dplyr::select('group',key_gene) %>% 
  ## gather,gather
  gather(key="gene",value="Expression",VEGFA,PLCB3,CALML3,ITPR3,ITPR2,ERBB2,BDKRB2,TBXA2R) %>% #c
  ##
  dplyr::select(group,gene,Expression,everything()) 
head(mydata)  ## 


p <- ggboxplot(mydata, x = "gene", y = "Expression",
               color = "group", palette = "jama",
               add = "jitter",
               title = paste0(GSE_accession,'_exp_boxplot')
)
#  Add p-value
pdf(paste0(GSE_accession,'_exp_boxplot.pdf'),width = length(key_gene),height = 5 )
p + stat_compare_means( label = "p.signif",
                        aes(group = group),
                        #label = "p.format"
                        method = 't.test'#"t.test"
)

dev.off()

#Immunity--
'********************************************************'
dir <- '09.Immunity'
setwd(paste0("E:/A /A /",xmmc))
if (!dir.exists(dir)) {dir.create(dir)}
setwd(dir)
write.csv(exp,paste0(GSE_accession,'_exp.csv'))
source('E:/A /A ///Cibersort.R')
imu_results <- CIBERSORT('E:/A /A ///LM22.txt',
                         paste0(GSE_accession,'_exp.csv'),#CSV
                         perm = 1000, 
                         QN = T)  #perm=1000QN=TRUE,TF

my_theme <- function(){
  theme(panel.grid = element_blank(),       # 
        panel.border = element_blank(),     # 
        legend.position="right",            # legend
        legend.text = element_text(size=8), # legend
        legend.title = element_text(size=8),# legend
        axis.line = element_line(size=1),   # 
        text = element_text(family="Times"),# 
        axis.text.y = element_text(size = 8,face='bold',color='black'),# y
        axis.text.x = element_text(size = 8,face='bold',color='black',angle=90,hjust=1),        # xangle=45  45 
        axis.title = element_text(size=10,face="bold"),  # 
        plot.title = element_text(hjust=0.5,size=10))    # 
}  


#
group 
res1 <- read.table('CIBERSORT-Results.txt',sep = '\t' ,header = T,row.names = 1)
res1 <- res1[,-(23:25)]
head(res1)[,1:4]
res2 <- t(res1) %>% as.data.frame()

res2 <- dplyr::mutate(res2,cell_type = row.names(res2)) %>% 
  dplyr::select(cell_type,everything()) 
head(res2)[,1:4]

# 
pdf(paste0(GSE_accession,'_immu.pdf'),width = 12,height = 6)
# a = c(sample(rainbow(22),22))
# b=a
# good_rainbow <- as.data.frame(b)
# write.csv(good_rainbow,'good_rainbow.csv')
good_rainbow <- read.csv('E:/A /A ///good_rainbow.csv',
                         header = T,row.names = 1)
b <- good_rainbow$b
res2 %>%
  gather(sample, fraction,-cell_type) %>%
  # 
  ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
  geom_bar(stat='identity',position="fill",alpha = 0.7) +
  ggtitle("title")+
  coord_flip() +
  scale_fill_manual(values= b) +
  #scale_x_discrete(limits = rev(levels(res2)))+
  my_theme()
dev.off()
#R
dim(res1)
box_data_imu <- res1  %>% 
  mutate(group1 = group$group) 
head(box_data_imu)[,1:4]

mydata <- box_data_imu %>% 
  dplyr::select(group1,everything()) %>% 
  ## gather,gather
  gather(key="cell.type",value="value",-group1) %>% #c
  ##
  dplyr::select(group1,cell.type,value,everything()) 
head(mydata)  ## 


mydata$group1 <- factor(mydata$group1,levels = c("Control","Case"))
levels(mydata$group1)

p <- ggboxplot(mydata, x = "cell.type", y = "value",
               color = "group1", palette = "jama",
               add = "jitter",
               title = paste0(GSE_accession,'_immu_boxplot')
)+
  theme (axis.text.x = element_text (angle = 45, vjust = 1, hjust=1))
#  Add p-value
pdf(paste0(GSE_accession,'_immu_boxplot.pdf'),
    width = length(res2$cell_type)/2,
    height = 5.5 )
p + stat_compare_means( label = "p.signif",
                        aes(group = group1),
                        #label = "p.format"
                        method = 't.test'#"t.test"
)

dev.off()


#
library(corrplot)
library(WGCNA)

key_gene_exp <- exp[key_gene,] %>% t() %>% as.data.frame
head(key_gene_exp)
dim(key_gene_exp)
imu_key <- res1
dim(res1)
imu_cor <- cor(key_gene_exp,imu_key,method = 'spearman')
nsample = length(row.names(res1))
imu_cor_P = corPvalueStudent(imu_cor ,nsample)
imu <- cbind(imu_cor,imu_cor_P)
write.csv(imu,'Gene-immune associations.csv')
# signif
textMatrix = paste(signif(imu_cor, 2), "\n(", signif(imu_cor_P, 1), ")", sep = "")
dim(textMatrix) = dim(imu_cor)

# 
pdf("Gene-immune associations.pdf",width = 8, height=length(key_gene)+0.5)
par(mar =c(5, 6.5, 3, 3));
# 
labeledHeatmap(Matrix = imu_cor, xLabels = colnames(imu_key), 
               yLabels = colnames(key_gene_exp), 
               cex.lab = 0.5, 
               ySymbols = colnames(key_gene_exp), colorLabels = FALSE, 
               colors = blueWhiteRed(30), 
               textMatrix = textMatrix, setStdMargins = FALSE, 
               cex.text = 0.5, zlim = c(-1,1),
               main = paste("Gene-immune relationships"))

dev.off()





#--kegg、GO enrich-------------------------------------------

'**********************************************'
dir <- '0.enrich'
setwd(paste0("E:/A /A /",xmmc))
if (!dir.exists(dir)) {dir.create(dir)}
setwd(dir)
library(clusterProfiler) # R
library(stringr) # 
library(AnnotationDbi)
library(org.Hs.eg.db) # hg19
library(DOSE)
library(ggplot2) # 
library(ggrepel) # 
library(tidyr)
library(tidyverse)
# gene symbolEntrez ID,
rich_gene = read.table("clipboard",header=T)
gene_df <- bitr(rich_gene$elements, 
                fromType = "SYMBOL", #fromTypeID
                toType = c("ENSEMBL", "ENTREZID"), #toTypeID
                OrgDb = org.Hs.eg.db)#Orgdb
head(gene_df)


# GO
go <- enrichGO(gene = gene_df$ENTREZID, # Entrez ID
               OrgDb = org.Hs.eg.db, # 
               keyType = "ENTREZID", # 
               ont = "ALL", # ,BP()/CC()/MF()/ALL()
               pAdjustMethod = "BH", # P,fdrnone
               pvalueCutoff = 1,
               qvalueCutoff = 0.2, # p/q qp0-10.2
               readable = T # IDsymbol
)
go.res <- data.frame(go) # GO
write.csv(go.res,"Table_GO_result.csv",quote = F) # GO
df <- gene_df$ENTREZID

kegg <- enrichKEGG(gene = df, 
                   organism = "hsa",
                   keyType = "kegg", 
                   pAdjustMethod = "none",
                   pvalueCutoff = 10,
                   qvalueCutoff = 0.2,#p0-1
                   minGSSize = 10,
                   maxGSSize = 500,
                   use_internal_data = F)
# IDsymbol
kk = setReadable(kegg,OrgDb = "org.Hs.eg.db",keyType = "ENTREZID")
head(kegg)
names(kegg)

kegg_1 = data.frame(kegg) %>% 
  as_tibble() %>% 
  separate(GeneRatio,into = c('gr1','gr2')) %>% 
  separate(BgRatio,into = c('br1','br2')) %>% 
  mutate_at(vars(gr1,gr2),as.numeric) %>% 
  mutate_at(vars(br1,br2),as.numeric) %>% 
  mutate(gr = gr1/gr2) %>% 
  mutate(br = br1/br2) %>% 
  mutate(ratio = gr/br)


head(kegg_1)
write.csv(kegg,'Table_KEGG_result.csv')
write.csv(kegg_1,'Table_KEGG_result_Calculated.csv')
 

###GO_tree.
source('E:/A /A ///.R')
library(tidygraph)
library(ggraph)
library(treemap)
library(ggsci)
#library(ccgraph)
pAdjustMethod <- "BH"  ## BH none
ego <- go
dim(ego) # 801
write.csv(as.data.frame(ego),"GO-enrich.csv",row.names = F)
sum(ego$ONTOLOGY == "CC") # Cellular component  106
sum(ego$ONTOLOGY == "MF") # Molecular function  160
sum(ego$ONTOLOGY == "BP") # Biological process  1567

cat(paste0("pAdjustMethod: ", pAdjustMethod, "\nGO: ", dim(ego)[1], "\nCC: ", sum(ego$ONTOLOGY == "CC"), "\nMF: ", sum(ego$ONTOLOGY == "MF"), "\nBP: ", sum(ego$ONTOLOGY == "BP")),
    file = paste0("0_Enrichment.log"), append = T, sep = '\n')

go_pic <- as.data.frame(ego)
go_pic <- rbind(go_pic[go_pic$ONTOLOGY == "CC",][1:10,], go_pic[go_pic$ONTOLOGY == "MF",][1:10,], go_pic[go_pic$ONTOLOGY == "BP",][1:10,])
go_pic$Description <- paste0(go_pic$ONTOLOGY, ": ", go_pic$Description)
# GO_bar

pdf(file = "GO_tree.pdf", width = 18, height = 10, onefile = FALSE)
treemap(go_pic,
        index="Description", #
        vSize="Count", #
        vColor="p.adjust", #
        type="value", #
        palette='Set3', #
        fontsize.labels=c(15, 15), #
        align.labels=list(c("center", "center"), c("left", "top")), #
        border.col="black", #  
        border.lwds=c(2,2),#
        title = "GO"
)
dev.off()


top5_kegg <- as.data.frame(kk)
top5_kegg <- top5_kegg %>% slice_min(n = 5, order_by = p.adjust)
top5_kegg <- top5_kegg %>% dplyr::select(Description, geneID, Count)
top5_kegg <- top5_kegg %>% as_tibble() %>% separate_rows(geneID, sep = "/") %>% as.data.frame()
top5_kegg$Count <- 1

index<- c("Description", "geneID")
nodes<- gather_graph_node(top5_kegg, index=index, value = "Count", root = "KEGG")
edges<- gather_graph_edge(top5_kegg, index=index, root = "KEGG")
graph<- tbl_graph(nodes = nodes, edges = edges)

p <- ggraph(graph,layout = 'dendrogram', circular = TRUE) +
  geom_edge_diagonal(aes(color=node1.node.branch),alpha=1/3) +
  geom_node_point(aes(size=node.size,color=node.branch),alpha=1/3) +
  coord_fixed()+
  theme_void()+
  theme(legend.position = "none")+
  scale_size(range = c(0.5,80)) +
  geom_node_text(aes(x = 1.15 * x,y = 1.15 * y,
                     label = node.short_name,
                     angle = -((-node_angle(x, y) + 90) %% 180) + 90,
                     filter = leaf,color = node.branch), size = 10,
                 # hjust = 'outward'
  ) +
  # scale_colour_manual(values= rep( brewer.pal(9,"Paired") , 30))+
  scale_colour_manual(values= pal_d3("category20")(20))+
  geom_node_text(aes(label=node.short_name,filter = !leaf,color = node.branch),
                 fontface="bold",size=7)

pdf(file = "KEGG_clus.pdf", width = 15, height = 15)
print(p)
dev.off()






#ceRNA------------------------------------------------------------
'*******************************ceRNA*****************************'
dir <- '10.ceRNA'
setwd(paste0("E:/A /A /",xmmc))
if (!dir.exists(dir)) {dir.create(dir)}
setwd(dir)


miRNAlncRNA 
#miRNA-lncRNA(Starbase)
library(data.table)
library(tidyverse)
starbase_lnc <- fread("E:/A /A //ceRNA/ENCORI_miRNA_lncRNA_lite.csv") %>% 
  as.data.frame()
head(starbase_lnc)[,1:6]
starbase_lnc <- starbase_lnc[,c(1,2,4,9)]
colnames(starbase_lnc)[1:4] <- c("miID","miRNA","lncRNA","pancancerNum")
#miRNA3miRNArR
miRNA <- read.table("ALL_gene_miRNA.csv",sep=",",header=T,check.names = FALSE) 
colnames(miRNA)[1:2] <- c('Gene_id',"miRNA")
#miRNA <- miRNA[,-1]
#miRNA <- miRNA %>% dplyr::select(miRNA=2) #
miRNA
result_starbase <- merge(starbase_lnc,miRNA,by="miRNA",all=F)
write.csv(result_starbase,"ALL_gene_miRNA_lncRNA_all.csv")

result_starbase1 <- result_starbase[result_starbase$pancancerNum>10,] #clipExpNum
length(table(result_starbase1$miRNA)) #miRNA
length(table(result_starbase1$lncRNA)) #mRNA
nrow(result_starbase1) #
write.csv(result_starbase1,file="imortant_miRNA_lnc.csv",row.names=F)

aa = as.name(paste0('top5_lnc'))
aa  = tibble()
for (i in 1:length(miRNA$miRNA)) {
  a = result_starbase[result_starbase$miRNA == miRNA$miRNA[i],] %>% 
    arrange(desc(pancancerNum))  
  #filter(pancancerNum > result_starbase$pancancerNum[6])
  a = a[1:5,]
  aa = rbind(aa,a)
}
write.csv(aa,paste0('allgene_top5_miRNA_lncRNA','.csv'))






