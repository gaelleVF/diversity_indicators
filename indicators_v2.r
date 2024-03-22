### packages & sources

### data
# only varieties with acreage and molecular data should be included in the analysis, i.e. other varieties must be removed before computing indicators
# in the case where all the varieties of a tested scenario do not represent 100% of the surfaces of the area studied (eg. due to non-exhaustive survey)
# or in the case where alleles of the microsatellites used are not informed for each variety (eg. due to technical problems arising during the genotyping procedure)
# the calculation of the indicators are based on relative proportions and allelic frequencies (see the main text for details)
# However, it should be noted that the interpretation of the results obtained should be considered with caution in the case of large proportions of missing data 

if(FALSE){
  # create files for proportions of varieties on R script (because it will be necessary to create a large number of scenarios with given structures)
  nrow_variety_proportion = 20
  ncol_variety_proportion = 3
  variety_data <- data.frame(matrix(NA,nrow=nrow_variety_proportion,ncol=ncol_variety_proportion))
  colnames(variety_data)  <- c("scenario","variety","proportion")
  variety_data$scenario   <- c(rep("scenario_1",nrow_variety_proportion/4),rep("scenario_2",nrow_variety_proportion/4),rep("scenario_3",nrow_variety_proportion/4),rep("scenario_4",nrow_variety_proportion/4))
  variety_data$variety    <- c("variety_1","variety_2","variety_3","variety_4","variety_5","variety_1","variety_2","variety_3","variety_4","variety_5","variety_1","variety_2","variety_3","variety_4","variety_5","variety_1","variety_2","variety_3","variety_4","variety_5")
  variety_data$proportion <- c(50,5,0,20,25,
                               10,10,5,25,50,
                               5,30,45,10,10,
                               0,50,47,2,1)
  variety_data$scenario <- as.factor(variety_data$scenario)
  variety_data$variety <- as.factor(variety_data$variety)
  str(variety_data)
  
  
  
  
  
  
  number_varieties <- length(which_varieties)
  scenarios <- list(
    c(50,5,0,20,25),
    c(10,10,5,25,50),
    c(5,30,45,10,10),
    c(0,50,47,2,1)
  )
  number_scenarios <- length(scenarios)
  
  variety_data <- data.frame(matrix(NA,nrow=number_varieties*number_scenarios,ncol=3))
  colnames(variety_data)  <- c("scenario","variety","proportion")
  
  variety_data$scenario   <- as.factor(rep(seq(1,number_scenarios,1), each=number_varieties))
  variety_data$variety    <- as.factor(rep(which_varieties,number_scenarios))
  variety_data$proportion <- unlist(scenarios)
  
  str(variety_data)
}

# 0.1. functions --------
get_indicators <- function(variety_data, ssr_data){
  
  sum_allelic_weighing_coefficient <- tapply(ssr_data$allelic_weighing_coefficient,ssr_data$variety, FUN=sum)/length(levels(ssr_data$locus))
  sum_allelic_weighing_coefficient # to verify that the sum of allelic weighing coefficients = 1, i.e. that the genetic information characterizing each variety has the same weight in the calculation of the different indicators
  
  ### indicator values
  # information on the data used to calculate indicators
#  scenario_levels <- levels(variety_data$scenario) # number of scenarios
  microsatellite_levels <- levels(ssr_data$locus) # number of microsatellites
  
  # storage dataframe of indicator values
  values_indicators <- data.frame(matrix(NA,nrow=length(grep("proportion",colnames(variety_data))),ncol=9))
  rownames(values_indicators) <- colnames(variety_data)[grep("proportion",colnames(variety_data))]
  colnames(values_indicators) <- c("varietal_richness","proportion_top5","N2","H_Bonneuil","Kimura_and_Crow",
                                   "Hstar_Bonneuil","Hstar_effective_number_of_alleles",
                                   "HTstar_Bonneuil","HTstar_effective_number_of_alleles")
  
  # loop to calculate indicator values for each scenario
  for (i in 1:nrow(values_indicators))
  {
    temporary_scenario <- variety_data[,grep(paste("variety",rownames(values_indicators)[i],sep="|"),colnames(variety_data))] # to choose proportion data for a specific scenario
    colnames(temporary_scenario)[grep(rownames(values_indicators)[i],colnames(temporary_scenario))] = "proportion"
    temporary_scenario$proportion = as.numeric(as.character(temporary_scenario$proportion))
    
    reorder_proportions   <- temporary_scenario[order(temporary_scenario$proportion,decreasing=TRUE),] # to reorder data (descending order of proportions)
    relative_proportion <- reorder_proportions$proportion/sum(reorder_proportions$proportion) # to calculate the relative proportions of the varieties of the scenario (i.e. sum of relative proportions = 1)
    varietal_richness   <- length(which(reorder_proportions$proportion != 0)) # Varietal richness
    proportion_top5     <- round(sum(reorder_proportions$proportion[1:5]),digits=2) # ProportionTop5
    N2                  <- 1/sum((relative_proportion)^2) # N2
    
    list_varieties <- reorder_proportions$variety # to list the names of the varieties of the specific scenario
    temporary_ssr <- subset(ssr_data, ssr_data$variety %in% list_varieties) # to choose varieties for the same specific scenario
    temporary_ssr_and_proportion <- merge(reorder_proportions,temporary_ssr,by = 'variety') # to combine proportion and microsatellite data in a single dataframe
    
    H_per_locus = NULL
    Kimura_and_Crow_per_locus = NULL
    Hstar_per_locus = NULL
    Hstar_effective_number_of_alleles_per_locus = NULL
    for (j in 1:length(microsatellite_levels))
    {
      temporary_per_locus <- temporary_ssr_and_proportion[temporary_ssr_and_proportion$locus == microsatellite_levels[j],] # to choose one locus
      temporary_per_locus$allele_value <- as.factor(temporary_per_locus$allele_value) # to transform quantitative allele values as a factor with different levels
      temporary_per_locus <- temporary_per_locus[!is.na(temporary_per_locus$allele_value) & temporary_per_locus$proportion != 0,] # get rid of missing data and varieties which proportion = 0
      
      temporary_per_locus <- droplevels(temporary_per_locus)
      
      if(nrow(temporary_per_locus)>0){
        # indicators unweighted by the relative proportions of the varieties
        pij_no_relative_proportion <- tapply(temporary_per_locus$allelic_weighing_coefficient, temporary_per_locus$allele_value, FUN=sum) # allelic frequencies to calculate the proportions of the different alleles of the locus (no relative proportions of alleles)
        pij <- pij_no_relative_proportion/sum(pij_no_relative_proportion) # allelic frequencies: to calculate the relative proportions of the different alleles of the locus (i.e. sum of relative proportions of alleles = 1)
        H_per_locus[j] <- 1 - (sum((pij)^2)) # to calculate a part of Nei's gene diversity (= H)
        Kimura_and_Crow_per_locus[j] <- 1 / (sum((pij)^2)) # to calculate a part of Kimura & Crow's effective number of alleles
        
        
        # indicators weighted by the relative proportions of the varieties
        temporary_per_locus$relative_proportion_for_ssr <- (temporary_per_locus$proportion*temporary_per_locus$allelic_weighing_coefficient)/(sum(temporary_per_locus$proportion*temporary_per_locus$allelic_weighing_coefficient)) # to calculate the relative proportions of the varieties of the scenario per locus (i.e. sum of relative proportions = 1)
        pstarij_no_relative_proportion <- tapply(temporary_per_locus$relative_proportion_for_ssr, temporary_per_locus$allele_value, FUN=sum) # weighted allelic frequencies to calculate the proportions of the different alleles of the locus weighted by the relative proportions of the varieties (no relative proportions of alleles)
        pstarij <- pstarij_no_relative_proportion/sum(pstarij_no_relative_proportion) # weighted allelic frequencies to calculate the relative proportions of the different alleles of the locus weighted by the relative proportions of the varieties (i.e. sum of relative proportions of alleles = 1)
        Hstar_per_locus[j] <- 1 - (sum((pstarij)^2)) # to calculate a part of Hstar index (H*)
        Hstar_effective_number_of_alleles_per_locus[j] <- 1 / (sum((pstarij)^2)) # to calculate a part of Hstar effective number of alleles index
      }
    }
    
    H_Bonneuil <- sum(na.omit(H_per_locus))/length(na.omit(H_per_locus)) # Nei's gene diversity (= H)
    Kimura_and_Crow <- sum(na.omit(Kimura_and_Crow_per_locus))/length(na.omit(Kimura_and_Crow_per_locus)) # Kimura & Crow's effective number of alleles
    Hstar_Bonneuil <- sum(na.omit(Hstar_per_locus))/length(na.omit(Hstar_per_locus)) # Hstar index (H*)
    Hstar_effective_number_of_alleles <- sum(na.omit(Hstar_effective_number_of_alleles_per_locus))/length(na.omit(Hstar_effective_number_of_alleles_per_locus)) # Hstar effective number of alleles index
    
    # indicators weighted by the relative proportions of the varieties and coefficient of intravarietal diversity
    GstL  <- 0.4  # coefficient for landraces (L)
    GstOL <- 0.94 # coefficient for old commercial lines (OL)
    GstML <- 1    # coefficient for modern pure lines (ML)
    GstPV <- 0.4    # coefficient for peasant varieties (PV)
    HTstar_denominator_no_relative_proportion <- tapply(temporary_per_locus$relative_proportion_for_ssr,temporary_per_locus$type, FUN=sum) # proportions of the different types of varieties weighted by the relative proportions of the varieties (no relative proportions of varieties)
    HTstar_denominator_no_relative_proportion[is.na(HTstar_denominator_no_relative_proportion)] <- 0 # to replace NA by 0 in the contexte of an absence of a type of variety
    HTstar_denominator <- as.data.frame(t(HTstar_denominator_no_relative_proportion/sum(HTstar_denominator_no_relative_proportion))) # proportions of the different types of varieties weighted by the relative proportions of the varieties (i.e. sum of relative proportions of varieties = 1)
    proportionL  <- ifelse(!is.null(HTstar_denominator$"variete de pays"),HTstar_denominator$"variete de pays",0)
    proportionOL <- ifelse(!is.null(HTstar_denominator$"lignee ancienne"),HTstar_denominator$"lignee ancienne",0)
    proportionML <- ifelse(!is.null(HTstar_denominator$"lignee pure moderne"),HTstar_denominator$"lignee pure moderne",0)
    proportionPV <- ifelse(!is.null(HTstar_denominator$"variete population"),HTstar_denominator$"variete population",0)
    HTstar_Bonneuil <- Hstar_Bonneuil / ((GstL * proportionL) + (GstOL * proportionOL) + (GstML * proportionML) + (GstPV * proportionPV))
    HTstar_effective_number_of_alleles <- Hstar_effective_number_of_alleles / ((GstL * proportionL) + (GstOL * proportionOL) + (GstML * proportionML) + (GstPV * proportionPV))
    
    
    
    values_indicators[i,1] <- varietal_richness
    values_indicators[i,2] <- proportion_top5
    values_indicators[i,3] <- N2
    values_indicators[i,4] <- H_Bonneuil
    values_indicators[i,5] <- Kimura_and_Crow
    values_indicators[i,6] <- Hstar_Bonneuil
    values_indicators[i,7] <- Hstar_effective_number_of_alleles
    values_indicators[i,8] <- HTstar_Bonneuil
    values_indicators[i,9] <- HTstar_effective_number_of_alleles
  }
  
  return(values_indicators)
}


# Créer les mélanges en amont, comme ça si on veut créer des mélanges spécifiques c'est possible
create_scenar <- function(rep_data, ssr_data, departement, annee, type_scenar, type_var, prop){
  ## Fonction creant les scenarios
  # rep_data = donnees de repartition des varietes dans les departements
  # departement = departement depuis lequel on veut creer le scenario
  # annee = vecteur avec les années sur lesquelles on veut créer le scénario
  # type_scenar = 
      # "proportion" = remplacer les variétés proportionnellement à leur fréquence dans le paysage
      # "random" = remplace aleatoirement les varietes du paysage sans prendre en compte les frequences
      # "domin" = remplacer les varietes dominantes (replacement en fonction de la fréquence de la variété)
      # "rare" = remplacer les varietes minoritaires
  # type_var = "MIX","OL","L","PV"
  # prop = proportion de remplacement. Vector dont la longueur est egal au nombre d'annees d'evolution du scenario
  
  
  
  # A verifier pour les scénarios : 
  # - qu'on a pas trop de données manquantes --> fixer un seuil + sortir les surfaces pour lesquelles on n'a pas d'infos et les ssr manquants
  # - qu'on a de l'info pour assez de ssr --> fixer un seuil (15 ?)
  
  d = rep_data[rep_data$departement %in% departement & rep_data$annee %in% annee,]
  variete = unique(d$variete)
  
  
  
  
  # verif proportion de varietes manquantes au ssr
  
  
  # si type_var == "m?lange", cr?er des m?langes ? partir des vari?t?s dominantes ? Des m?langes ? combien ? 
  # Faire par exemple 2-3 m?langes diff?rents et les tirer al?atoirement
  # donner une ann?e de d?part et puis remplacer les vari?t?s d'ann?e en ann?e dans les vari?t?s existantes
  
  # choix des d?partements = ceux contrast?s pour lesquels on a des donn?es pour beaucoup de vari?t?s
  
}

lcm <- function(x, y) {
  # find the least common multiple of 2 numbers
  # choose the greater number
  if(x > y) {
    greater = x
  } else {
    greater = y
  }
  while(TRUE) {
    if((greater %% x == 0) && (greater %% y == 0)) {
      lcm = greater
      break
    }
    greater = greater + 1
  }
  return(lcm)
}

create_mixture <- function(varieties, num_var, ssr_data, name_mix){
  ### function that creates ssr data for mixtures
        # varieties : vector of varieties to mix
        # num_var : number of individuals for each variety
        # name_mix : name of the mixture
        
  if(name_mix %in% ssr_data$variety){stop("This mixture name is already in ssr_data")}
  
  temp_ssr <- ssr_data[ssr_data$variety %in% varieties,]
  temp_ssr <- droplevels(temp_ssr)
  
  # if one of the component of the mixture already has several individuals in ssr_data, multiply the number for the other components
  nb_ind <- unique(temp_ssr[,c("variety","individual")])
  nb_ind <- as.matrix(by(nb_ind$individual,list(nb_ind$variety),max))
  nb_ind[,1] <- as.numeric(as.character(nb_ind[,1]))
  
  if(length(unique(nb_ind[,1])) > 1){
    LCM=1
    for (i in 1:nrow(nb_ind)){
      LCM = lcm(LCM,nb_ind[i])
    }
    nb_ind <- LCM/nb_ind[,1]
  }else{
    nb_ind <- nb_ind[,1]
  }

  
  nb_microsat <- length(unique(temp_ssr$locus))
    
  ssr <- NULL
  comp=1
  for (i in 1:length(varieties)){
    v = varieties[i]
    ssr_v <- temp_ssr[temp_ssr$variety %in% v,]
    ssr_v <- ssr_v[order(ssr_v$individual,ssr_v$locus),]
    if(nrow(ssr_v)>0){
      ssr <- rbind(ssr, cbind(ssr_v[rep(1:nrow(ssr_v),times=num_var[i]*nb_ind[v]),],rep(seq(comp,comp+num_var[i]*nb_ind[v]*max(ssr_v$individual)-1,1),each=nb_microsat)))
      comp <- comp + num_var[i]*nb_ind[v]*max(ssr_v$individual)
    }else{
      warning(paste("No ssr data for ", v,sep=""))
    }
  }
  
  ssr$variety <- name_mix
  ssr$individual <- ssr[,ncol(ssr)]
  ssr <- ssr[,-ncol(ssr)]
  ssr$allelic_weighing_coefficient <- 1/(2*max(ssr$individual))
  
  return(ssr)
}



# 0.2. get data ----------
repartition_data <- read.csv2("./script_indicateurs/surfaces_varietes/surface_variete_1981_2017_vf.csv",header=T,sep=";")
str(repartition_data)

# use a predefined file for microsatellite data and import on R script (because the characterization of the varieties must not be modified according to the scenario, i.e. imputation of missing data must be done before)
# ssr_data <- read.csv2("microsatellite_markers_reformatted_example.csv",header=TRUE,dec=".",row.names=1)
ssr_data <- get(load("microsatellite_reformatted_all.RData"))
str(ssr_data)

# 1.0. create scenarios ---------

mix1 = create_mixture(c("SOISSONS","ALTIGO","OVALO"),c(1,1,1),ssr_data,name_mix="MIX1")
ssr_data <- rbind(ssr_data,mix1)

# Create automaticaly the matrix with the proportion of varieties
which_varieties <- list(
                      "1" = as.character(repartition_data[repartition_data$departement %in% "ain" & repartition_data$annee %in% 2006,"variete"]),
                      "2" = c("MIX1","ROCALOEX","ROUGEDUROC"),
                      "3" = c("SOISSONS","ALTIGO","OVALO","ROCALOEX","ROUGEDUROC"),
                      "4" = c("SOISSONS","ALTIGO","OVALO","ROCALOEX","DAUPHIBOIS")
                    )

proportions <- list(
                  "1" = matrix(as.numeric(as.character(repartition_data[repartition_data$departement %in% "ain" & repartition_data$annee %in% 2006,"repartition_pct"],ncol=1))),
                  "2" = cbind(c(15,35,50),c(30,40,30),c(90,5,5)),
                  "3" = cbind(c(5,5,5,35,50),c(10,10,10,40,30),c(30,30,30,5,5)),
                  "4" = matrix(c(0,50,47,2,1),ncol=1)
                )

# Check if same number of varieties and proportions for each scenario
    a <- sapply(which_varieties,length)
    b <- sapply(proportions,nrow)
    d <- which(a != b)
    if(length(d)>0){stop("Attention : le nombre de variété et de proportions pour le(s) scénario(s) ",d," ne concordent pas.",sep="")}
    
# Check if varieties' names are in ssr_data, give proportion of missing data
    noms <- as.character(unique(ssr_data$variety))
    a <- sapply(which_varieties,function(x){x %in% noms})
    b <- grep(FALSE,a)
    prop <- sapply(b,function(x){
      d <- a[[x]]
      prop <- round(sum(proportions[[x]][d == FALSE]),2)
      return(prop)
    })
    if(length(prop)>0){warning(paste("Attention : la proportion de l'assolement pour laquelle il manque des données ssr pour le(s) scenario(s) ",paste(b,collapse=" et "),
                                     " est de ", paste(prop,collapse=" et "),sep=""))}
    
    
# create scenarios    
scenarios <- lapply(1:length(which_varieties), function(i){
  a <- as.data.frame(cbind(which_varieties[[i]], proportions[[i]]))
  colnames(a) <- c("variety",paste("proportion_",seq(1,ncol(a)-1,1)))
  return(a)
} )
    

# 2.0. Calculate indicators ----
indicators <- lapply(scenarios,function(x) get_indicators(x, ssr_data))


# 3.0. save data ------
OUT <- lapply(1:length(indicators),function(i){
  a <- list(scenarios[[i]],indicators[[i]])
  names(a) <- c("scenarios","indicators")
  return(a)
  })


save(OUT, file="./resultats/OUT_indicators.RData")

### export indicator values
#write.table(values_indicators,"values_indicators.csv", sep=";", row.names=TRUE, col.names=TRUE)



### literature references
# Kimura and Crow 1964. Genetics 49:725-738.
# Nei 1973. Proc. Natl. Acad. Sci. USA 70:3321-3323.
# Jost 2008. Molecular Ecology 17:4015-4026.
# Bonneuil et al. 2012. Ecological Indicators 23:280-289.
# Perronne et al. 2017. Agriculture, Ecosystems and Environment 236:12-20.

