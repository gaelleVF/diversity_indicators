setwd("/home/deap/Documents/Gaelle/Selection_Participative/Ecoagri/BILAN_PROJET_E442_DEAP")
library(rebus)
library(reshape2)


# 0. Functions -----------
get_allelic_distrib <- function(M, mkr){
  
  for (i in mkr){M[,i]=as.factor(M[,i])}
  
  tab = lapply(grep(paste(mkr,collapse="|"),colnames(M)), function(x){return(cbind(names(table(M[,x])),table(M[,x])))})
  n <- max(sapply(tab, nrow)) 
  t = do.call(cbind, lapply(tab, function (x)
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
  
  colnames(t) = rep(mkr,each=2)
  write.table(t,"summary_marqueurs.csv")
}

format_data <- function(ssr, list_mkr, save="microsatellite_markers_reformated.csv"){
  
  ssr[ssr == 0] <- NA
  ssr[ssr == "."] <- NA
  
  if(length(grep("[_]",ssr$Variete)) == nrow(ssr)){
    ssr$code <- unlist(lapply(as.character(ssr$Variete),function(x) strsplit(x,"_")[[1]][1]))
    ssr$Variete <- unlist(lapply(as.character(ssr$Variete),function(x) strsplit(x,"_")[[1]][2]))
  }
  
  if(!("individu")%in% colnames(ssr)){
    ssr <- ssr[order(ssr$Variete),]
    ssr$individu <- unlist(lapply(unique(ssr$Variete), function(x){
      return(seq(1,table(ssr$Variete)[x],1)) 
    }))
  }
  ssr[is.na(ssr$individu),"individu"] <- 1
  
  # mkr <- gsub("a|b","",list_mkr)
  num_mkr <- table(list_mkr)
  
  # check corrections
  get_allelic_distrib(ssr,list_mkr)
  
  
  M <- melt(ssr, measure.vars = list_mkr)
  M <- M[,grep("Variete|individu|TYPE|Type|variable|value",colnames(M))]
  colnames(M)[grep("TYPE|Type",colnames(M))] = "type"
  colnames(M)[grep("Variete",colnames(M),fixed=TRUE)] = "variety"
  colnames(M)[grep("individu",colnames(M))] = "individual"
  colnames(M)[grep("variable",colnames(M))] = "locus"
  colnames(M)[grep("value",colnames(M))] = "allele_value"
  
  M <- M[order(M$variety,M$locus),]
  
  m <- rbind(M[rep(which(M$locus %in% names(num_mkr[num_mkr==1])), each = 2),], 
             M[which(!(M$locus %in% names(num_mkr[num_mkr==1]))),])
  
  m <- m[order(m$variety,m$individual,m$locus),]
  m$position <- rep(c(1,2),times = nrow(m)/2)
  m$individual <- as.numeric(as.character(m$individual))
  
  num_ind <- as.matrix(by(m$individual, m$variety, max))
  num_ind <- as.data.frame(cbind(rownames(num_ind),num_ind))
  
  #m$allelic_weighing_coefficient <- unlist(lapply(m$variety, function(x){1/(2*as.numeric(num_ind[x,2]))}))
  #m$allelic_weighing_coefficient <- unlist(lapply(m$variety, function(x){1/(2*as.numeric(as.character(num_ind[x,2])))}))
  
  # prendre en compte les NA dans le calcul du coefficient de ponderation par marqueur
  coef <- NULL
  for (i in 1:nrow(m)){
    print(i)
    #b <- m[grep(m[i,"variety"],m$variety),]
    b <- m[grep(paste("^",m[i,"variety"],"$",sep=""),m$variety),]
    b <- b[b$locus %in% m[i,"locus"],]
    b <- b[!is.na(b$allele_value),]
    if(nrow(b)>0){a <- 1/nrow(b)}else{a=0}
    coef <- c(coef,a)
  }
  m$allelic_weighing_coefficient <- coef
  
  
  # m$locus = gsub("a|b","",m$locus)
  
  write.table(m,file=save, row.names=FALSE)
  
}

# 1.1. Format data ---------------- 
ssr <- read.table("06_Resultats_marqueurs_edites_2.csv",header=T,sep=";")
#ssr <- read.table("./Rémi/v%c3%a9rification_ssr_to_matrix/06_Resultats_marqueurs_edites_RP.csv",header=T,sep=";")

ssr <- ssr[!is.na(ssr$Numero),]  # delete empty rows
ssr <- ssr[,grep("Numero|Variete|TYPE_VAR|Type|individu|ok",colnames(ssr))]  # keep only edited markers
colnames(ssr) <- unlist(lapply(colnames(ssr), function(x) strsplit(x,"_")[[1]][1]))
list_mkr <- unique(colnames(ssr)[grep("gwm|Cfd",colnames(ssr))])

pb_allelic_weighing_coefficient <- format_data(ssr,list_mkr)
pb_allelic_weighing_coefficient$V2 <- as.numeric(as.character(pb_allelic_weighing_coefficient$V2))

# RP: les deux lignes precedentes ne fonctionnent pas toute seule
# mais la fonction formar_data fonctionne, donc possible de retirer ces lignes


format_data(ssr,list_mkr,save="microsatellite_markers_reformated.csv")

# comme indique dans le mail, il y a un probleme avec les noms de varietes inclus dans
# d'autres noms de varietes (ex. LONA, AMI), j'ai donc ajoute varietyNOMvariety, ce qui evite
# tout probleme et ne prend pas beaucoup de temps, car malgre differents essais sur les fonctions
# je n'arrive pas e savoir pour cela ne fonctionne pas dans la version actuelle
# il y a encore un probleme de doublon FLECHEDOR

# GVF : il faut coller "^" nom de la variété et "$" pour ne conserver que les chaines de caractères qui correspondent exactement.
# C'est modifié dans la fonction
# Pourquoi il y a deux variétés qui s'appelent FLECHEDOR ?


ssr_old <- read.table("./Rémi/v%c3%a9rification_ssr_to_matrix/base_donnees_mol_08_07_08_RP.csv",header=T,sep=";")
#ssr_old <- read.table("base_donnees_mol_08_07_08_RP_exempleAMI.csv",header=T,sep=";")
ssr_old$variete <- gsub("variety","",ssr_old$variete) 
colnames(ssr_old)[1:4] <- c("Numero","Variete","TYPE","individu")
list_mkr <- unique(colnames(ssr_old)[grep("gwm|cfd",colnames(ssr_old))])
format_data(ssr_old,list_mkr,"old_microsatellite_markers_reformated.csv")


# RP: j'ai rebricole la sortie (je ne sais pourquoi il n'y a pas de correspondance entre les deux)
# lorsque je lance le script (peut-etre as-tu fais des modifs dans base_donnees_mol_08_07_08)e

# 1.2. Merge data
old <- read.table("old_microsatellite_markers_reformated.csv",header=T,sep=" ")
new <- read.table("microsatellite_markers_reformated.csv",header=T,sep=" ")
new$locus = gsub("C","c",new$locus)
new$type = gsub("L","l",new$type)
new$type = gsub("V","v",new$type)

references = c("SOISSONS","APACHE","SHANGO","VIVANT","ETOILEDECHOISY","APEXAL","PROGRESS","MAGDALENA","LONA","ARCANE","PRINQUAL","FURIO","FLOREAL","VICTO","MAVERICK","NEWTON")
# PROGRESS ET maverick
old <- old[!(old$variety %in% references),]
old <- old[old$locus %in% unique(new$locus),]

#write.table(old, file="old_verification_retrait_reference.csv", sep=";")
# RP: ok avec ce choix

M <- rbind(old,new)
M <- M[order(M$variety, M$locus),]
M$locus <- droplevels(M$locus)

save(M, file="microsatellite_reformatted_all.RData")
write.table(M, file="microsatellite_reformatted_all.csv", sep=";")


