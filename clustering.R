require(VarSelLCM)
xbrut = read.csv2("hemopred.csv", na.strings = c("NA", "A", "D", "K"), dec=".", stringsAsFactors=F)

# remove ID column
#xbrut = xbrut[,-1]

colnames(xbrut)
categ_vars = c("SEXE", "ATCD1", "ATCD2", "ATCD3", "ATCD4", "ATCD5", "ATCD6", "ADMIS_MOTIF", "CHIR_ABDO", "CHIR_THOR", "SDRA", "SEPTUM")
num_vars = setdiff(colnames(xbrut), categ_vars)
colstypes = rep(NA, length(ncol(xbrut)))
colstypes[which(colnames(xbrut) %in% categ_vars)] = "factor"
colstypes[which(colnames(xbrut) %in% num_vars)] = "numeric"
x = read.csv2("hemopred_sedki.csv", na.strings = c("NA", "A", "D", "K"), dec=".", stringsAsFactors=F, colClasses=colstypes)
x <- x[,-1]
summary(x)

## clustering par MICL 
# Cluster analysis with variable selection (with parallelisation)
res_micl  = list()
for(k in 2:10)
   res_micl[[k-1]] <- VarSelCluster(x, k, nbcores = 40, initModel=100, crit.varsel = "MICL")


res_bic  = list()
for(k in 2:10)
  res_bic[[k-1]] <- VarSelCluster(x, k, nbcores = 40, nbSmall=500, nbKeep = 100, crit.varsel = "MICL")

save(res_micl, res_bic, file="res_clust_hemopred.rda")