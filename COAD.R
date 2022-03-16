################################################
### Step 1: variable selection
################################################
load("~/Downloads/TF/COAD.RData")
source("HID.R")
formula <- as.formula("Surv(PFI.time, PFI) ~ .")
intero <- hidRF(formula, data = dat, subRF = TRUE,
                nSubRF = 8,
                var.sel.only = TRUE,
                forest.tune = TRUE)
intero$rf.obj

sel <- intero[[1]]$var.sl
save(intero,file = "COADstep1.RData")

dat2 = dat[,c("PFI.time", "PFI", sel)]
rm(intero)
gc()
write.csv(sel, file = "COADmajorSingnal.csv")

################################################
### Step 2-3: interaction selection
################################################
intero <- hidRF(Surv(PFI.time, PFI) ~ ., data = dat2, 
                terms = 6, block = 1,
                subRF = FALSE,
                forest.tune = TRUE)
intero$rf.obj
intero$interaction

result <- Pvalue(intero)

save(intero,result, file = "COADstep2.RData")
write.csv(result$interaction[,c("PBII","Pvalue")], file = "COADtopInteractions.csv")
