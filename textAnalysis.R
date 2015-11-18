##installation of packages and loading all required packages.
#source("http://bioconductor.org/biocLite.R")
#biocLite("Rgraphviz")
library("Rgraphviz")

for (package in c('tm', 'SnowballC', 'wordcloud', 'topicmodels', 'ggplot2', 'wordcloud', 'colorspace', 'fpc', 'cluster', 'FactoMineR', 'plot3D', 'plot3Drgl', "MASS")) {
  if (!require(package, character.only=T)) {
    install.packages(package)
    library(package, character.only=T)
  }
}

rm(package)

##Tools Functions, run these functions first
#frequent terms
GetFreqTerms <- function(tdm){
  low.freq <- as.numeric(winDialogString("Enter low frequency amount", "100"))
  
  term.freq <- rowSums(as.matrix(tdm))
  term.freq <- subset(term.freq, term.freq >= low.freq)
  df <- data.frame(term = names(term.freq), freq = term.freq)
  df$term <- factor(df$term, levels = df$term)
  df <- df[order(df$freq), ]
  
  ggplot(df, aes(reorder(term, freq), freq)) + geom_bar(stat = "identity") + 
    xlab("Terms") + ylab("Count") + coord_flip()  
}

#word associations
GetWordAssoc <- function(tdm){
  term.interest <- winDialogString("Enter the term of interest to find association", "")
  cor.limit <- as.numeric(winDialogString("Enter lower correlation bound limit to find association", "0.1"))
  term.assocs <- data.frame(corr = findAssocs(tdm, term.interest, cor.limit)[,1],
                            terms = row.names(findAssocs(tdm, term.interest, cor.limit)))
  term.assocs$terms <- factor(term.assocs$terms ,levels = term.assocs$terms)
  
  x <- ggplot(term.assocs, aes( y = terms  ) ) +
    geom_point(aes(x = corr), data = term.assocs) +
    xlab(paste0("Correlation with the term ", "\"", term.interest, "\""))
  
  print(x)
  
  assoc.terms <- as.character(term.assocs$terms)
  assoc.terms[length(assoc.terms)+1] <- term.interest
  return(assoc.terms)
}

#plot term document matrix by frequency
PlotTDMbyFreq <- function(tdm){
  low.freq <- as.numeric(winDialogString("Enter low frequency amount", "100"))
  freq.terms <- findFreqTerms(tdm, lowfreq = low.freq)
  
  cor.threshold <- as.numeric(winDialogString("Enter correlation threshold to plot term document matrix", "0.1"))

  fontsize <- length(freq.terms)
  if (fontsize < 25) fontsize <- 20
  
  for (i in 1:10){
    cor.threshold <- i/10
    plot(tdm, "dot", term = freq.terms, main = paste("Correlation Threshold", cor.threshold, sep=" "), corThreshold = cor.threshold, weighting = F, attrs=list(node=list(shape="plaintext", fontsize=fontsize), edge=list(color="grey", minlen=2)))
  }
}

#plot term document matrix
PlotTDMbyAssocTerms <- function(tdm, assoc.terms, magnify){
  cor.threshold <- as.numeric(winDialogString("Enter correlation threshold to plot term document matrix", "0.1"))
  for (i in 1:10){
    cor.threshold <- i/10
    #plot(tdm, "dot", term = assoc.terms, main = paste("Correlation Threshold", cor.threshold), corThreshold = cor.threshold, weighting = F, 
    #     attrs=list(node=list(shape="plaintext", fontsize=length(assoc.terms)), edge=list(color="grey", minlen=2)))    
    plot(tdm, "dot", term = assoc.terms, main = paste("Correlation Threshold", cor.threshold), corThreshold = cor.threshold, weighting = F, 
         attrs=list(node=list(shape="plaintext", fontsize=length(assoc.terms)*magnify), edge=list(color="grey", minlen=2)))    
  }
}

#wordcloud
PlotWordcloud <- function(tdm){
  m <- as.matrix(tdm)
  word.freq <- sort(rowSums(m), decreasing = T)
  wordcloud(words = names(word.freq), freq = word.freq, random.order = F, max.words = 100, colors=brewer.pal(8, "Dark2"))
}

#Clustering
#hierarchical clustering
GetHClust <- function(tdm){
  sparse.num <- as.numeric(winDialogString("Enter maximal allowed sparsity (0-99.99) of terms to be included in clustering (the smaller, the lesser the terms)", "90"))/100
  
  tdm2 <- removeSparseTerms(tdm, sparse = sparse.num)
  m2 <- as.matrix(tdm2)
  
  distMatrix <- dist(scale(m2))
  fit <- hclust(distMatrix, method = "ward.D")
  
  plot(fit, cex = 0.9, hang = -1)
  Sys.sleep(5)
  
  k <- as.numeric(winDialogString("Enter number of clusters", "5"))
  rect.hclust(fit, k = k)
}

#kmeans clustering
GetKMeansClust <- function(tdm){
  sparse.num <- as.numeric(winDialogString("Enter maximal allowed sparsity (0-99.99) of terms to be included in clustering (the smaller, the lesser the terms)", "90"))/100
  
  tdm2 <- removeSparseTerms(tdm, sparse = sparse.num)
  m2 <- as.matrix(tdm2)
  
  m3 <- t(m2)
  k <- as.numeric(winDialogString("Enter number of clusters", "3"))
  kmeansResults <- kmeans(m3, k)
  round(kmeansResults$centers, digits = 3)
  
  for (i in 1:k){
    cat(paste("cluster", i, ": ", sep = ""))
    s <- sort(kmeansResults$centers[i, ], decreasing = T)
    cat(names(s)[1:10], "\n")
  }  
}

# partitioning around metroids
GetPamClust <- function(tdm){
  sparse.num <- as.numeric(winDialogString("Enter maximal allowed sparsity (0-99.99) of terms to be included in clustering (the smaller, the lesser the terms)", "90"))/100
  
  tdm2 <- removeSparseTerms(tdm, sparse = sparse.num)
  m2 <- as.matrix(tdm2)
  
  m3 <- t(m2)
  
  pamResult <- pamk(m3, metric = "manhattan")
  
  k <- pamResult$nc
  pamResult <- pamResult$pamobject
  
  for (i in 1:k){
    cat("cluster", i, ": ", 
        colnames(pamResult$medoids)[which(pamResult$medoids[i, ] == 1)], "\n")
  }
  
  layout(matrix(c(1, 2), 1, 2))
  plot(pamResult, col.p = pamResult$clustering)
  
  layout(matrix(1))
}

# correspondance analysis
PlotCorrespondenceAnalysis <- function(tdm, k){
  sparse.num <- as.numeric(winDialogString("Enter maximal allowed sparsity (0-99.99) of terms to be included in clustering (the smaller, the lesser the terms)", "90"))/100
  
  tdm2 <- removeSparseTerms(tdm, sparse = sparse.num)
  m1 <- as.matrix(tdm2)
  
  # remove columns with all zeros
  m1 = m1[,colSums(m1)!=0]
  
  # for convenience, every matrix entry must be binary (0 or 1)
  m1[m1 > 1] = 1
  
  tdm.ca = CA(m1, ncp = 3, graph=FALSE)
  
  # default plot of words
  #plot(tdm.ca$row$coord, type="n", xaxt="n", yaxt="n", xlab="", ylab="")
  #text(tdm.ca$row$coord[,1], tdm.ca$row$coord[,2], labels=rownames(m1),
  #     col=hsv(0,0,0.6,0.5))
  #title(main="Correspondence Analysis of words", cex.main=1)
  
  # partitioning around medoids 
  if (k == 0){
    pamResult <- pamk(m1, metric = "manhattan")
    k <- pamResult$nc    
  }
  
  # pam clustering
  tdm.pam = pam(tdm.ca$row$coord[,1:2], k)
  
  # get clusters
  clusters = tdm.pam$clustering
  
  
  # create data frame
  tdm.df = data.frame(
    words = rownames(m1),
    dim1 = tdm.ca$row$coord[,1],
    dim2 = tdm.ca$row$coord[,2],
    dim3 = tdm.ca$row$coord[,3],
    freq = rowSums(m1),
    cluster = as.factor(clusters))
  
  # plot
  x <- ggplot(tdm.df, aes(x=dim1, y=dim2, label=words)) +
    geom_text(aes(size=freq, colour=cluster), alpha=0.7) +
    scale_size_continuous(breaks=seq(20,80,by=10), range=c(3,8)) +
    scale_colour_manual(values=brewer.pal(12, "Paired")) +
    theme_bw() + 
    labs(x="Dim 1", y="Dim 2")  
  print(x)

  x <- ggplot(tdm.df, aes(x=dim2, y=dim3, label=words)) +
    geom_text(aes(size=freq, colour=cluster), alpha=0.7) +
    scale_size_continuous(breaks=seq(20,80,by=10), range=c(3,8)) +
    scale_colour_manual(values=brewer.pal(12, "Paired")) +
    theme_bw() + 
    labs(x="Dim 2", y="Dim 3")  
  print(x)
  
  x <- ggplot(tdm.df, aes(x=dim1, y=dim3, label=words)) +
    geom_text(aes(size=freq, colour=cluster), alpha=0.7) +
    scale_size_continuous(breaks=seq(20,80,by=10), range=c(3,8)) +
    scale_colour_manual(values=brewer.pal(12, "Paired")) +
    theme_bw() + 
    labs(x="Dim 1", y="Dim 3")
  print(x)
}

# correspondance analysis, 3Dplot
PlotCorrespondenceAnalysis3D <- function(tdm, k){
  sparse.num <- as.numeric(winDialogString("Enter maximal allowed sparsity (0-99.99) of terms to be included in clustering (the smaller, the lesser the terms)", "90"))/100
  
  tdm2 <- removeSparseTerms(tdm, sparse = sparse.num)
  m1 <- as.matrix(tdm2)
  
  # remove columns with all zeros
  m1 = m1[,colSums(m1)!=0]
  
  # for convenience, every matrix entry must be binary (0 or 1)
  m1[m1 > 1] = 1
  
  tdm.ca = CA(m1, ncp = 3, graph=FALSE)
  
  # partitioning around medoids 
  if(k == 0){
    pamResult <- pamk(m1, metric = "manhattan")
    k <- pamResult$nc
  }
  
  # pam clustering
  tdm.pam = pam(tdm.ca$row$coord[,1:3], k)
  
  # get clusters
  clusters = tdm.pam$clustering
  
  
  # create data frame
  tdm.df = data.frame(
    words = rownames(m1),
    dim1 = tdm.ca$row$coord[,1],
    dim2 = tdm.ca$row$coord[,2],
    dim3 = tdm.ca$row$coord[,3],
    freq = rowSums(m1),
    cluster = as.factor(clusters))
  
  
  with(tdm.ca, text3D(tdm.ca$row$coord[, 1], tdm.ca$row$coord[, 2], tdm.ca$row$coord[, 3], 
                      labels = rownames(tdm.ca$row$coord), #colvar = cluster, 
                      col = gg.col(100), theta = 60, phi = 20,
                      xlab = "Dim 1", ylab = "Dim 2", zlab = "Dim 3", 
                      main = "Text", cex = 0.6, 
                      bty = "g", ticktype = "detailed", d = 2)) #,
  #clab = c("Urban","Pop"), adj = 0.5, font = 2))
  
  with(tdm.df, text3D(dim2, dim1, dim3, 
                      labels = words, colvar = as.numeric(cluster), 
                      col = gg.col(100), theta = 60, phi = 20,
                      xlab = "Dim 1", ylab = "Dim 2", zlab = "Dim 3", 
                      main = "Text", cex = 0.6, 
                      bty = "g", ticktype = "detailed", d = 2)) #,
  
  # Make the rgl version
  
  plotrgl()
  
  fix(tdm.df)
  write.csv(tdm.df, "temp.csv")
}


####Scripts here####
#read csv file into data frame
text.file <- file.choose()
text.all <- read.csv(text.file, header = TRUE, sep = ",")
influencers <- read.csv("influencers.csv", header = TRUE, sep = ",")

text.all$CONTENT <- gsub('http\\S+\\s*', '', text.all$CONTENT) #remove links


#read corpus from column in data frame. change column name when needed
winDialog(type = "ok", "Opening a view of the first rows of the text dataset\n\nPlease note the column name of the column containing the text to be analyzed.\n\nCopy the exact name (case sensitive) of the column and paste it into the next dialog box")
View(head(text.all))

Sys.sleep(5)
col.num <- winDialogString("Enter the name of the column containing text to be analyzed", "")
#text <- text.all[, which( colnames(text.all)==col.num)]
text <- text.all[text.all$RT == 0, which( colnames(text.all)==col.num)]
text <- tolower(text)
text <- unique(text)

#text <- gsub('juice', "", text)
#for (i in influencers$Name){
#  text <- text[grep(i, text, invert = T)]
#}
#text <- text[grep("rockyandmayur", text, invert = T)]
#text <- text[grep("techotweet", text, invert = T)]
#text <- text[grep("murray", text, invert = T)]
#text <- text[grep("vickychef", text, invert = T)]
#text <- text[grep("vantaskigoli", text, invert = T)]


text.corpus <- Corpus(VectorSource(text))


####sentiments analysis WIP####
#library(devtools)
#install_github("mannau/tm.plugin.sentiment")
#library(tm.plugin.sentiment)

#text.corpus2 <- score(text.corpus)

#text.all.sent <- cbind(text.all, data.frame(meta(text.corpus2)))
#text.all.sent$sentiment <- 0

#text.all.sent[text.all.sent$polarity == 0 | is.na(text.all.sent$polarity), "sentiment"] <- "NEUTRAL"
#text.all.sent[text.all.sent$polarity < 0 & !is.na(text.all.sent$polarity), "sentiment"] <- "NEGATIVE"
#text.all.sent[text.all.sent$polarity > 0 & !is.na(text.all.sent$polarity), "sentiment"] <- "POSITIVE"

#plot distribution of polarity
#ggplot(text.all.sent, aes(x = sentiment)) +
#  geom_bar(aes(y = ..count.., fill = polarity)) +
#  scale_fill_brewer(palette = "RdGy") + 
#  labs(x = "polarity categories", y = "number of comments")

#seperate words by polarity
#pol = levels(factor(text.all.sent$sentiment))
#npol = length(pol)
#pol.docs = rep("", npol)

#for (i in 1:npol){
#  tmp = text.all.sent[text.all.sent$sentiment == pol[i], which( colnames(text.all)==col.num)]
#  pol.docs[i] = paste(tmp, collapse = " ")
#}




#create corpus
#corpus.sent <- Corpus(VectorSource(pol.docs))

#dtm.control <- list(tolower = TRUE, 
#                    removePunctuation = TRUE, 
#                    removeNumbers = TRUE, 
#                    stopwords = c(stopwords("english")), 
#                    stemming = FALSE, #stemming reduces complexity and duplicated words from grammar etc 
#                    wordLengths = c(3, Inf), 
#                    weighting = weightTf)

#tdm.sent <- TermDocumentMatrix(corpus.sent, control = dtm.control)
#tdm.sent <- as.matrix(tdm.sent)
#colnames(tdm.sent) <- pol

#comparison word cloud
#comparison.cloud(tdm.sent, colors = brewer.pal(npol, "Dark2"), 
#                 random.order = FALSE, title.size = 1.5)


#extract term document matrix
dtm.control <- list(tolower = TRUE, 
                    removePunctuation = TRUE, 
                    removeNumbers = TRUE, 
                    stopwords = c(stopwords("english")), 
                    stemming = TRUE, #stemming reduces complexity and duplicated words from grammar etc 
                    wordLengths = c(3, Inf), #usually 3 to inf 
                    weighting = weightTf)


tdm <- TermDocumentMatrix(text.corpus, control = dtm.control)
inspect(tdm[1:10,1:5])


##Run desired functions below
set.seed(12345) #set seed for reproducible results

#wordcloud
PlotWordcloud(tdm)

#exploration by frequency
GetFreqTerms(tdm)
PlotTDMbyFreq(tdm)

#exploration by word association
assoc.terms <- GetWordAssoc(tdm)
PlotTDMbyAssocTerms(tdm, assoc.terms, 2)

#clustering 
GetHClust(tdm)
GetKMeansClust(tdm)
GetPamClust(tdm)
PlotCorrespondenceAnalysis(tdm, 3)
PlotCorrespondenceAnalysis3D(tdm, 7)

