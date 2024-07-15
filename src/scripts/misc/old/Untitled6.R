tbl=read.table("~/projects/mito_asd/admixture/postimpute.merged.ldpruned.2.Q")
indTable = read.table("~/Data/MyProject/admixture/MyProject.HO.merged.ind",
                      col.names = c("Sample", "Sex", "Pop"))
popGroups = read.table("~/Google Drive/GA_Workshop Jena/HO_popGroups.txt", col.names=c("Pop", "PopGroup"))

mergedAdmixtureTable = cbind(tbl, indTable)
mergedAdmWithPopGroups = merge(mergedAdmixtureTable, popGroups, by="Pop")
ordered = mergedAdmWithPopGroups[order(mergedAdmWithPopGroups$PopGroup),]
barplot(t(as.matrix(subset(ordered, select=V1:V6))), col=rainbow(6), border=NA)