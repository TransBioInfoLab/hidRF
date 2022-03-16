################################################
### Step 1: variable selection
################################################
load("~/Downloads/TF/BRCA.RData")
source("HID.R")
formula <- as.formula("Surv(PFI.time, PFI) ~ .")
intero <- hidRF(formula, data = dat, subRF = TRUE,
                nSubRF = 8,
                var.sel.only = TRUE,
                forest.tune = TRUE)
intero$rf.obj

sel <- intero[[1]]$var.sl
save(intero,file = "BRCAstep1.RData")

dat2 = dat[,c("PFI.time", "PFI", sel)]
rm(intero)
gc()
write.csv(sel, file = "BRCAmajorSingnal.csv")

################################################
### Step 2-3: interaction selection
################################################
intero <- hidRF(Surv(PFI.time, PFI) ~ ., data = dat2, 
                terms = 4, block = 1,
                subRF = FALSE,
                forest.tune = TRUE)
intero$rf.obj
intero$interaction

result <- Pvalue(intero)

save(intero,result, file = "BRCAstep2.RData")
write.csv(result$interaction[,c("PBII","Pvalue")], file = "BRCAtopInteractions.csv")
