#Initialisation
library(dplyr)
library(gtsummary)
library(car)
library(ggplot2)
library(scales)
library(tibble)
library(VarSelLCM)

load("hemopred_data.rda")
library(readxl)
x_complete <- read_excel("hemopred_complete_db.xlsx")
etio_outcome <- x_complete[,c(1,138,140,142,144,146,147,149,151)]
x_total <- x %>% left_join(etio_outcome,by="ID")
x_total[,c(45:53)] <- lapply(x_total[,c(45:53)],factor)
x_total <- x_total[,c(-1,-3)]

#Table 1 
tbl_summary(x_total)
summary(x_total)

#Figure 1
load("res_clust_AllVars_m2.rda")

micl_val = rep(NA, length(res_micl))
for(k in 1:length(res_micl))
  micl_val[k] = res_micl[[k]]@criteria@MICL

plot(2:15, micl_val)
abline(v= which.max(micl_val) + 1, col="red")

#Table 3
database_clustered_hd <- cbind(x_total,res_micl[[4]]@partitions@zMAP)
colnames(database_clustered_hd)[52] <- "Clusters_MICL"
tbl_summary(database_clustered_hd[,c(1:43,52)],by=Clusters_MICL)

#Figure 2
load("resultats_M2_M4.rda")
bestm = res_micl_M2[[4]]
#group = paste("class-",1:5, sep="") 
group = c("Défaillance gauche","Défaillance droite","Réanimation optimale","Hypovolémie sévère","Hypovolémie modérée")
res_radar = data.frame(cbind(group, 
                             t(bestm@param@paramContinuous@mu),
                             bestm@param@paramCategorical@alpha[[1]],
                             bestm@param@paramCategorical@alpha[[2]]))
colnames(res_radar) = c("group",colnames(t(bestm@param@paramContinuous@mu)), 
                        paste(names(bestm@param@paramCategorical@alpha)[1], "-",c("0","1"),sep=""),
                        paste(names(bestm@param@paramCategorical@alpha)[2], "-",c("0","1"), sep=""))
rownames(res_radar)=NULL
res_radar = res_radar[, which(colnames(res_radar) %in% c("group",
                                                         bestm@model@names.relevant,
                                                         paste(names(bestm@param@paramCategorical@alpha)[1], "-",c("0","1"),sep=""),
                                                         paste(names(bestm@param@paramCategorical@alpha)[2], "-",c("0","1"), sep="")))]
res_radar[,-1] = apply(res_radar[,-1], 2, function(x){rescale(as.numeric(x))})
res_radar <- res_radar[,c(-12,-14)]
res_radar <- res_radar[,c(1,3,4,8,9,11,10,12,2,5,6,7,13)]
colnames(res_radar) <- c("group","Base Excess","Lactate","LVEF","LVRF","Vmax E","RV/LV Ratio","Paradox. Septum", "PP Variations",
                         "SVC Index","IVC Index","VC Variations","PLR Response")

# Library
library(fmsb)

par(mfrow = c(3, 2))

#RADAR LV FAILURE
LV_fail <- res_radar[1,-1]
LV_fail <- rbind(rep(1,12),rep(0,12),LV_fail)
radarchart(LV_fail, axistype=1 , 
           pcol=rgb(0.84,0.1,0.11,0.9), pfcol=rgb(0.84,0.353,0.353,0.5), plwd=3, 
           cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,20,5), cglwd=0.8,
           vlcex=1,title=paste("Cluster 1 : LV Failure"))

#RADAR RV FAILURE
RV_fail <- res_radar[2,-1]
RV_fail <- rbind(rep(1,12),rep(0,12),RV_fail)
radarchart(RV_fail, axistype=1 , 
           pcol=rgb(0.169,0.514,0.729,0.9), pfcol=rgb(0.314,0.588,0.729,0.5), plwd=3, 
           cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,20,5), cglwd=0.8,
           vlcex=1,title=paste("Cluster 2 : RV Failure"))

#RADAR WELL.RESC
Well_resc <- res_radar[3,-1]
Well_resc <- rbind(rep(1,12),rep(0,12),Well_resc)
radarchart(Well_resc, axistype=1 , 
           pcol=rgb(0.2,0.655,0.357,0.9), pfcol=rgb(0.671,0.867,0.643,0.5), plwd=3, 
           cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,20,5), cglwd=0.8,
           vlcex=1,title=paste("Cluster 3 : Well-resuscitated"))

#RADAR SEVERE HYPOVOL
Severe_hypo <- res_radar[4,-1]
Severe_hypo <- rbind(rep(1,12),rep(0,12),Severe_hypo)
radarchart(Severe_hypo, axistype=1 , 
           pcol=rgb(0.867,0.867,0.392,0.9), pfcol=rgb(1,1,0.749,0.5), plwd=3, 
           cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,20,5), cglwd=0.8,
           vlcex=1,title=paste("Cluster 4 : Severe Hypovolemia"))

#RADAR MODERATE HYPOVOL
Mod_hypo <- res_radar[5,-1]
Mod_hypo <- rbind(rep(1,12),rep(0,12),Mod_hypo)
radarchart(Mod_hypo, axistype=1 , 
           pcol=rgb(0.992,0.5,0.259,0.9), pfcol=rgb(0.992,0.682,0.38,0.5), plwd=3, 
           cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,20,5), cglwd=0.8,
           vlcex=1,title=paste("Cluster 5 : Moderate Hypovolemia"))

#Table 4
tbl_summary(database_clustered_hd[,c(44:52)],by=Clusters_MICL)

#Figure 3
require(circlize)

test_mat <- matrix(c(32,1,3,54,0,1,1,
                     49,2,1,14,8,3,14,
                     79,1,3,29,0,0,24,
                     29,6,28,11,2,0,7,
                     106,13,54,18,0,1,18),nrow=5,byrow=T)

colnames(test_mat) <- c("SEPSIS","HEMORRHAGE","HYPOVOLEMIA","CARDIOGENIC","PULM.EMBOL.","ANAPHYLAXIS","OTHER")
rownames(test_mat) <- c("LV.FAILURE","RV.FAILURE","WELL.RESUSC","SEV.HYPOVOL","MODER.HYPOVOL")

par(mfrow = c(2, 4))

# GENERAL CHORD
grid.col = c(SEPSIS = "#e6f598",HEMORRHAGE = "#d53e4f",HYPOVOLEMIA = "#fc8d59", CARDIOGENIC="#3288bd",
             PULM.EMBOL. = "#99d594", ANAPHYLAXIS = "#fee08b",OTHER = "#ffffbf",
             LV.FAILURE = "grey",RV.FAILURE = "grey",WELL.RESUSC = "grey", SEV.HYPOVOL = "grey", MODER.HYPOVOL = "grey")

chordDiagram(t(test_mat), annotationTrack = "grid",grid.col=grid.col)

circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, 1.5, CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA)

# CHORD SEPSIS
grid.col.sepsis = c(SEPSIS = "#e6f598",HEMORRHAGE = "grey",HYPOVOLEMIA = "grey", CARDIOGENIC="grey",
                    PULM.EMBOL. = "grey", ANAPHYLAXIS = "grey",OTHER = "grey",
                    LV.FAILURE = "grey",RV.FAILURE = "grey",WELL.RESUSC = "grey", SEV.HYPOVOL = "grey", MODER.HYPOVOL = "grey")

chordDiagram(t(test_mat), annotationTrack = "grid",grid.col=grid.col.sepsis)

circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, 1.5, CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA)

# CHORD HEMORRHAGE
grid.col.hemor = c(SEPSIS = "grey",HEMORRHAGE = "#d53e4f",HYPOVOLEMIA = "grey", CARDIOGENIC="grey",
                   PULM.EMBOL. = "grey", ANAPHYLAXIS = "grey",OTHER = "grey",
                   LV.FAILURE = "grey",RV.FAILURE = "grey",WELL.RESUSC = "grey", SEV.HYPOVOL = "grey", MODER.HYPOVOL = "grey")

chordDiagram(t(test_mat), annotationTrack = "grid",grid.col=grid.col.hemor)

circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, 1.5, CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA)

# CHORD HYPOVOLEMIA
grid.col.hypovol = c(SEPSIS = "grey",HEMORRHAGE = "grey",HYPOVOLEMIA = "#fc8d59", CARDIOGENIC="grey",
                     PULM.EMBOL. = "grey", ANAPHYLAXIS = "grey",OTHER = "grey",
                     LV.FAILURE = "grey",RV.FAILURE = "grey",WELL.RESUSC = "grey", SEV.HYPOVOL = "grey", MODER.HYPOVOL = "grey")

chordDiagram(t(test_mat), annotationTrack = "grid",grid.col=grid.col.hypovol)

circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, 1.5, CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA)

# CHORD CARDIOGENIC
grid.col.cardio = c(SEPSIS = "grey",HEMORRHAGE = "grey",HYPOVOLEMIA = "grey", CARDIOGENIC="#3288bd",
                    PULM.EMBOL. = "grey", ANAPHYLAXIS = "grey",OTHER = "grey",
                    LV.FAILURE = "grey",RV.FAILURE = "grey",WELL.RESUSC = "grey", SEV.HYPOVOL = "grey", MODER.HYPOVOL = "grey")

chordDiagram(t(test_mat), annotationTrack = "grid",grid.col=grid.col.cardio)

circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, 1.5, CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA)

# CHORD PE
grid.col.pe = c(SEPSIS = "grey",HEMORRHAGE = "grey",HYPOVOLEMIA = "grey", CARDIOGENIC="grey",
                PULM.EMBOL. = "#99d594", ANAPHYLAXIS = "grey",OTHER = "grey",
                LV.FAILURE = "grey",RV.FAILURE = "grey",WELL.RESUSC = "grey", SEV.HYPOVOL = "grey", MODER.HYPOVOL = "grey")

chordDiagram(t(test_mat), annotationTrack = "grid",grid.col=grid.col.pe)

circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, 1.5, CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA)

# CHORD ANAPHYLAXIS
grid.col.ana = c(SEPSIS = "grey",HEMORRHAGE = "grey",HYPOVOLEMIA = "grey", CARDIOGENIC="grey",
                 PULM.EMBOL. = "grey", ANAPHYLAXIS = "#fee08b",OTHER = "grey",
                 LV.FAILURE = "grey",RV.FAILURE = "grey",WELL.RESUSC = "grey", SEV.HYPOVOL = "grey", MODER.HYPOVOL = "grey")

chordDiagram(t(test_mat), annotationTrack = "grid",grid.col=grid.col.ana)

circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, 1.5, CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA)

# CHORD OTHER
grid.col.other = c(SEPSIS = "grey",HEMORRHAGE = "grey",HYPOVOLEMIA = "grey", CARDIOGENIC="grey",
                   PULM.EMBOL. = "grey", ANAPHYLAXIS = "grey",OTHER = "#ffffbf",
                   LV.FAILURE = "grey",RV.FAILURE = "grey",WELL.RESUSC = "grey", SEV.HYPOVOL = "grey", MODER.HYPOVOL = "grey")

chordDiagram(t(test_mat), annotationTrack = "grid",grid.col=grid.col.other)

circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, 1.5, CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA)

#Figure 4
load("res_clust_AllVars_m4.rda")

micl_val = rep(NA, length(res_micl))
for(k in 1:length(res_micl))
  micl_val[k] = res_micl[[k]]@criteria@MICL

plot(2:15, micl_val)
abline(v= which.max(micl_val) + 1, col="red")

#Table 5
database_clustered_hd_echo <- cbind(x_total,res_micl[[4]]@partitions@zMAP)
colnames(database_clustered_hd_echo)[52] <- "Clusters_MICL"
tbl_summary(database_clustered_hd_echo[,c(1:43,52)],by=Clusters_MICL)

#Table 6
tbl_summary(database_clustered_hd_echo[,c(44:52)],by=Clusters_MICL) 

#Figure 5
library(networkD3)
library(tidyr)
library(tibble)
library(dplyr)

links <- data.frame(
  source=c("Sev. Hypovol. (4)","Sev. Hypovol. (4)","Sev. Hypovol. (4)",
           "RV Fail. (4)","RV Fail. (4)","RV Fail. (4)",
           "LV Fail. (4)","LV Fail. (4)","LV Fail. (4)","LV Fail. (4)","LV Fail. (4)",
           "Well-resc (4)","Well-resc (4)","Well-resc (4)","Well-resc (4)",
           "Mod. Hypovol. (4)","Mod. Hypovol. (4)","Mod. Hypovol. (4)","Mod. Hypovol. (4)","Mod. Hypovol. (4)"),
  target=c("Sev. Hypovol. (2)","RV Fail. (2)","Mod. Hypovol. (2)",
           "RV Fail. (2)","Well-resc (2)","Mod. Hypovol. (2)",
           "Sev. Hypovol. (2)","RV Fail. (2)","LV Fail. (2)","Well-resc (2)","Mod. Hypovol. (2)",
           "RV Fail. (2)","LV Fail. (2)","Well-resc (2)","Mod. Hypovol. (2)",
           "Sev. Hypovol. (2)","RV Fail. (2)","LV Fail. (2)","Well-resc (2)","Mod. Hypovol. (2)"),
  value=c(28,2,4,54,7,2,1,2,73,9,4,20,5,100,10,42,4,3,9,161)
)

nodes <- data.frame(name=c(as.character(links$source),as.character(links$target)) %>% unique())

links$IDsource <- match(links$source,nodes$name)-1
links$IDtarget <- match(links$target,nodes$name)-1

sankey_net <- sankeyNetwork(Links = links, Nodes = nodes, Source = "IDsource", Target = "IDtarget", Value = "value",
                            NodeID = "name",sinksRight = FALSE)
sankey_net

#Table 7
load("resultats_M2_M4.rda")
clustered_x <- cbind(x_total,Clustering_MICL_M2,Clustering_MICL_M4)
non_migrating_cluster1 <- clustered_x %>% filter(Clustering_MICL_M2 == 1 & Clustering_MICL_M4 == 1)
tbl_summary(non_migrating_cluster1)
non_migrating_cluster5 <- clustered_x %>% filter(Clustering_MICL_M2 == 5 & Clustering_MICL_M4 == 5)
tbl_summary(non_migrating_cluster5)
migrating_cluster <- clustered_x %>% filter(Clustering_MICL_M2 == 5 & Clustering_MICL_M4 == 1)
tbl_summary(migrating_cluster)

#Figure S1
load("res_clust_AllVars_m1.rda")

micl_val = rep(NA, length(res_micl))
for(k in 1:length(res_micl))
  micl_val[k] = res_micl[[k]]@criteria@MICL

plot(2:15, micl_val)
abline(v= which.max(micl_val) + 1, col="red")

#Table S1
database_clustered_fullbase <- cbind(x_total,res_micl[[which.max(micl_val)]]@partitions@zMAP)
colnames(database_clustered_fullbase)[52] <- "Clusters_MICL"
tbl_summary(database_clustered_fullbase[],by="Clusters_MICL")

#Figure S2
load("res_clust_AllVars_m3.rda")

micl_val = rep(NA, length(res_micl))
for(k in 1:length(res_micl))
  micl_val[k] = res_micl[[k]]@criteria@MICL

plot(2:15, micl_val)
abline(v= which.max(micl_val) + 1, col="red")

#Table S2
database_clustered_hd_clinical <- cbind(x_total,res_micl[[5]]@partitions@zMAP)
colnames(database_clustered_hd_clinical)[52] <- "Clusters_MICL"
tbl_summary(database_clustered_hd_clinical,by=Clusters_MICL)