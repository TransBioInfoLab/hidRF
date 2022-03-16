setwd("~/Downloads/Gene_expression_with_Distance_and_Inensity/TNBC_1160920F")
source("HID.R")
################################################
### merge data
################################################
x <- as.data.frame(t(read.csv("1160920F_SCTnormalized_genes_cytokine_or_chemokine.csv", row.names=1)))
y <- read.csv("1160920F_tumor_tissues_dist_and_intensity.csv", row.names=1)
rownames(y) <- sub("-", ".", rownames(y))

length(rownames(y)) ## check match
length(which((rownames(y) ==rownames(x))==TRUE)) ## check match

################################################
### create outcome
################################################
cluster <- paste("dist_to_cluster",1:12, sep =".")
wt <- paste("cluster_intensity",1:12, sep =".")
yname <- paste("y",1:12,sep ="")
for (i in 1:12){
  x[,yname[i]] <- y[,cluster[i]]/log(1+y[,wt[i]])
}

################################################
### Step 1: variable selection
################################################
formula <- as.formula(paste("cbind(",paste(yname,collapse = ","),")~.", collapse = ""))
intero <- hidRF(formula, data = x, subRF = TRUE,
                nSubRF = 10,
                 var.sel.only = TRUE)
intero$rf.obj

sel <- intero[[1]]$var.sl
# write.csv(sel, file = "1160920F_majorSingnalTNBC_cytokine_or_chemokine.csv", row.names = FALSE)

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

# write.csv(result$interaction[,c("PBII","Pvalue")], file = "1160920F_topInteractionsTNBC_cytokine_or_chemokine.csv")
