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
  
  #m$allelic_weighing_coefficient <- unlist(lapply(m$variety, function(x){1/(2*as.numeric(num_ind[x,2]))}))
  #m$allelic_weighing_coefficient <- unlist(lapply(m$variety, function(x){1/(2*as.numeric(as.character(num_ind[x,2])))}))
  
  # Prendre en compte les NA dans le calcul du coef de pondération, du coup le faire marqueur par marqueur
  coef <- NULL
  for (i in 1:nrow(m)){
    b <- m[grep(m[i,"variety"],m$variety),]
    b <- b[b$locus %in% m[i,"locus"],]
    b <- b[!is.na(b$allele_value),]
    if(nrow(b)>0){a <- 1/(2*nrow(b))}else{a=0}
    coef <- c(coef,a)
  }
  m$allelic_weighing_coefficient <- coef

  
  
  # m$locus = gsub("a|b","",m$locus)
  
  write.table(m,file=save, row.names=FALSE)
  
}

# 1.1. Format data ---------------- 
ssr <- read.table("06_Resultats_marqueurs_edites.csv",header=T,sep=";")
ssr <- ssr[!is.na(ssr$Numero),]  # delete empty rows
ssr <- ssr[,grep("Numero|Variete|TYPE_VAR|individu|ok|ec",colnames(ssr))]  # keep only edited markers
colnames(ssr) <- unlist(lapply(colnames(ssr), function(x) strsplit(x,"_")[[1]][1]))
list_mkr <- unique(colnames(ssr)[grep("gwm|Cfd",colnames(ssr))])

# RP: pourquoi ne pas inclure gwm99 dans la liste en mettant -1 sur la derniere lignee
# ce que tu ne fais pas pour le old dataset

# GVF : en effet c'est une erreur, il y avait avant une colonne "code" en dernière colonne que j'ai du retirer par la suite.
# je corrige ça

pb_allelic_weighing_coefficient <- format_data(ssr,list_mkr)
pb_allelic_weighing_coefficient$V2 <- as.numeric(as.character(pb_allelic_weighing_coefficient$V2))

format_data(ssr,list_mkr,save="microsatellite_markers_reformated.csv")

# RP: pourquoi par variete, la somme des donnees disponibles pour chaque microsatellite n'est-elle pas de 1
# ex. pour JAPHABELLE, tu as 10 individus et potentiellement heterozygotie pour chaque individu
# donc si tu les individus ont ete genotypees pour le marqueur, tu devrais avoir 0.05 de ponderation
# dans la colonne allelic weighing coefficient
# pour les varietes lignees pures cela fonctionne
# Rem: j'ai trouve le probleme que j'ai modifie (ligne de calcul m$allelic_weighing_coefficient dans la fct format_data)

# GVF : ok, j'ai passé le strplit au début de la fonction format.data pour n'avoir que le nom de la variété et pas
# le code du génotypage afin d'avoir bien le même nom de variété pour chaque individu


# RP: par contre, le script suppose, pour les varietes avec plusieurs individus, qu'il n'y a
# pas de NA ou qu'il y a des NA pour l'ensemble des individus, cependant pour certains marqueurs
# certains individus ont ete genotypes et d'autres pas. Dans ce cas, il me semble tout de meme
# important que le allelic_weighing_coefficient somme e 1, donc il faudrait calculer
# le coefficient par marqueur
# Rem: le plus rapide est peut-etre de le faire "e la main" sur excel, car nous n'avons pas
# beaucoup de varietes dans cette situation

# verification donnees brutes
# RP: j'ai verifie les correspondances entre mon fichier BDD et le tien.
# Le tien est plus complet car:
# - je n'ai pas le marqueur gwm415 (pas utilise dans notre etude oe nous avons utilise les memes marqueurs que l'etude de Bonnin)
# - je n'ai pas un certain nombre de varietes anciennes utilisees dans l'etude de Bonneuil, car je me suis focalise sur 1980-2006
# Par ailleurs, dans certains cas, il existe des differences (d'un pas) entre les deux fichiers sur gwm120, mais cela est rare
# et je ne vois pas e quelle etape cette "correction par regroupement a ete realisee
# enfin certaines varietes n'ont pas la meme nomenclature (parenthese ou tiret lorsque les varietes sont en plusieurs noms, des espaces derriere certains noms)
# mais elles sont presentes dans les deux fichiers =>
# Je pense donc qu'il est preferable de se baser sur ton fichier
# et pour gwm120 de trouver un moyen d'homogeneiser avec les nouvelles donnees, mais le cas se
# presente aussi pour les autres marqueurs, Harry ayant ete plus precis que les donnees
# des analyses precedentes

# mise en forme base de donnees
#ssr_old <- read.table("base_donnees_mol_08_07_08.csv",header=T,sep=";")
#colnames(ssr_old)[1:4] <- c("code","annee_inscription","TYPE","Variete")
#list_mkr <- unique(colnames(ssr_old)[5:ncol(ssr_old)])
#format_data(ssr_old,list_mkr,"old_microsatellite_markers_reformated.csv")

ssr_old <- read.table("base_donnees_mol_08_07_08_RP.csv",header=T,sep=";")
colnames(ssr_old)[1:4] <- c("Numero","Variete","TYPE","individu")
list_mkr <- unique(colnames(ssr_old)[5:ncol(ssr_old)])
format_data(ssr_old,list_mkr,"old_microsatellite_markers_reformated.csv")

# RP: j'avais des erreurs dans la sortie, j'ai donc modifie le fichier d'entree de la fonction format_data
# pour qu'il soit similaire au precedent fichier 06_Resultats_marqueurs_edites
# cela marche bien ensuite

# 1.2. Merge data
old <- read.table("old_microsatellite_markers_reformated.csv",header=T,sep=" ")
new <- read.table("microsatellite_markers_reformated.csv",header=T,sep=" ")
new$locus = gsub("C","c",new$locus)
new$type = gsub("L","l",new$type)
new$type = gsub("V","v",new$type)

references = c("SOISSONS","APACHE","SHANGO","VIVANT","ETOILEDECHOISY","APEXAL","PROGRESS","MAGDALENA","LONA","ARCANE","PRINQUAL","FURIO","FLOREAL","VICTO","MAVERICK","NEWTON")
old <- old[!(old$variety %in% references),]
old <- old[old$locus %in% unique(new$locus),]

#write.table(old, file="old_verification_retrait_reference.csv", sep=";")
# RP: ok avec ce choix

M <- rbind(old,new)
M <- M[order(M$variety, M$locus),]
M$locus <- droplevels(M$locus)

save(M, file="microsatellite_reformatted_all.RData")
write.table(M, file="microsatellite_reformatted_all.csv", sep=";")


