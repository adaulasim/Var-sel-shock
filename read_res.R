require(VarSelLCM)
load("res_clust_AllVars_m2.rda")

## Avec le critère BIC
length(res_bic)
bic_val = rep(NA, length(res_bic))
for(k in 1:length(res_bic))
  bic_val[k] = res_bic[[k]]@criteria@BIC 

## Choix du nombre de classes
plot(2:18, bic_val)
abline(v= which.max(bic_val) + 1, col="red")

## Les effectifs des classes
table(res_bic[[which.max(bic_val)]]@partitions@zMAP)

## variables pertinentes 
best_bic = res_bic[[which.max(bic_val)]]
best_bic@model@names.relevant

VarSelShiny(res_bic[[which.max(bic_val)]])


## Avec le critère MICL
micl_val = rep(NA, length(res_micl))
for(k in 1:length(res_micl))
  micl_val[k] = res_micl[[k]]@criteria@MICL

plot(2:18, micl_val)
abline(v= which.max(micl_val) + 1, col="red")
## Les effectifs des classes
table(res_micl[[which.max(micl_val)]]@partitions@zOPT)

## variables pertinentes 
best_micl = res_micl[[which.max(micl_val)]]
best_micl@model@names.relevant


VarSelShiny(res_micl[[which.max(micl_val)]])
