setwd("~/Downloads/Gene_expression_with_Distance_and_Inensity/TNBC_1142243F")
source("HID.R")

################################################
### merge data
################################################
x <- as.data.frame(t(read.csv("1142243F_SCTnormalized_genes_cytokine_or_chemokine.csv", row.names=1)))
y <- read.csv("1142243F_tumor_tissues_dist_and_intensity.csv", row.names=1)
rownames(y) <- sub("-", ".", rownames(y))
length(rownames(y)) ## check match
length(which((rownames(y) ==rownames(x))==TRUE)) ## check match

################################################
### create outcome
################################################
x$y <- y$dist_to_cluster.1
x$y <- x$y/(log(25)+1)

################################################
### Step 1: variable selection
################################################
formula <- as.formula(y~.)
intero <- hidRF(formula, data = x, subRF = TRUE,
                itrSub = 10,
                var.sel.only = TRUE)
intero$rf.obj

sel <- intero[[1]]$var.sl

#write.csv(sel, file = "1142243F_majorSingnalTNBC_cytokine_or_chemokine.csv", row.names = FALSE)

dat2 = x[,c(paste("y", 1:12, sep =""), sel)]
rm(intero)
gc()

################################################
### Step 2-3: interaction selection
################################################
intero <- hidRF(formula, data = dat2, 
                terms = 6, block = 1,
                subRF = FALSE)
intero$rf.obj
head(intero$interaction)

result <- Pvalue(intero)
head(result$interaction)
# write.csv(result$interaction[,c("PBII","Pvalue")], file = "1142243F_topInteractionsTNBC_cytokine_or_chemokine.csv")

################################################
### Model Comparison: iRF
################################################
library(iRF)
o2 <- iRF(x = as.matrix(dat2[,-1]), y = dat2[, 1],
          interactions.return=TRUE)
o2$interaction

# write.csv(o2$interaction[[1]], file = "iRF1142243F_topInteractionsTNBC_cytokine_or_chemokine.csv")

