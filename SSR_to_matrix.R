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
  if(!("individu")%in% colnames(ssr)){
    ssr <- ssr[order(ssr$Variete),]
    ssr$individu <- unlist(lapply(unique(ssr$Variete), function(x){
      return(seq(1,table(ssr$Variete)[x],1)) 
    }))
  }
  ssr[is.na(ssr$individu),"individu"] <- 1
  
  if(length(grep("[_]",ssr$Variete)) == nrow(ssr)){
    ssr$code <- unlist(lapply(as.character(ssr$Variete),function(x) strsplit(x,"_")[[1]][1]))
    ssr$Variete <- unlist(lapply(as.character(ssr$Variete),function(x) strsplit(x,"_")[[1]][2]))
  }

  # mkr <- gsub("a|b","",list_mkr)
  num_mkr <- table(list_mkr)
  
  # check corrections
  get_allelic_distrib(ssr,list_mkr)
  
  
  M <- melt(ssr, measure.vars = list_mkr)
  M <- M[,grep("Variete|individu|TYPE|variable|value",colnames(M))]
  colnames(M)[grep("TYPE",colnames(M))] = "type"
  colnames(M)[grep("Variete",colnames(M))] = "variety"
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
  
  m$allelic_weighing_coefficient <- unlist(lapply(m$variety, function(x){1/(2*as.numeric(num_ind[x,2]))}))
  
  # m$locus = gsub("a|b","",m$locus)
  
  write.table(m,file=save, row.names=FALSE)
  
}

# 1.1. Format data ---------------- 
ssr <- read.table("06_Resultats_marqueurs_edites.csv",header=T,sep=";")
ssr <- ssr[!is.na(ssr$Numero),]  # delete empty rows
ssr <- ssr[,grep("Numero|Variete|TYPE_VAR|individu|ok|ec",colnames(ssr))]  # keep only edited markers
colnames(ssr) <- unlist(lapply(colnames(ssr), function(x) strsplit(x,"_")[[1]][1]))
list_mkr <- unique(colnames(ssr)[5:(ncol(ssr)-1)])
format_data(ssr,list_mkr,save="microsatellite_markers_reformated.csv")

ssr_old <- read.table("base_donnees_mol_08_07_08.csv",header=T,sep=";")
colnames(ssr_old)[1:4] <- c("code","annee_inscription","TYPE_VAR","Variete")
list_mkr <- unique(colnames(ssr_old)[5:ncol(ssr_old)])
format_data(ssr_old,list_mkr,"old_microsatellite_markers_reformated.csv")


# 1.2. Merge data
old <- read.table("old_microsatellite_markers_reformated.csv",header=T,sep=" ")
new <- read.table("microsatellite_markers_reformated.csv",header=T,sep=" ")
new$locus = gsub("C","c",new$locus)
new$type = gsub("L","l",new$type)
new$type = gsub("V","v",new$type)

references = c("SOISSONS","APACHE","SHANGO","VIVANT","ETOILEDECHOISY","APEXAL","PROGRESS","MAGDALENA","LONA","ARCANE","PRINQUAL","FURIO","FLOREAL","VICTO","MAVERICK","NEWTON")
old <- old[!(old$variety %in% references),]
old <- old[old$locus %in% unique(new$locus),]


M <- rbind(old,new)
M <- M[order(M$variety, M$locus),]
M$locus <- droplevels(M$locus)

save(M, file="microsatellite_reformatted_all.RData")
write.table(M, file="microsatellite_reformatted_all.csv", sep=";")


