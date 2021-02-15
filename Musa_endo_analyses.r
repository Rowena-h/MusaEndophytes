#####################################################
#####################################################
#####################################################
######                                         ######
######    Fungal Endophytes in Musa CWR Seeds  ######
######                                         ######
#####################################################
#####################################################
#####################################################


library(ape)
library(colorspace)
library(dplyr)
library(egg)
library(eulerr)
library(ggplot2)
library(ggplotify)
library(ggpubr)
library(ggrepel)
library(ggstance)
library(ggtree)
library(ggwordcloud)
library(grid)
library(metR)
library(multcompView)
library(plyr)
library(reshape2)
library(RVAideMemoire)
library(stringr)
library(vegan)


########################
## OTU CLASSIFICATION ##
########################

#Read in all data
df.controls <- read.csv("data/All_data.csv")
#Read in metadata
metadata <- read.csv("data/metadata.csv")
#Remove controls
df.all <- subset(df.controls, df.controls$Direct.Sequencing.or.Culture != "Control")
#Add Musa species
df.all$MusaSpecies <- metadata$MusaSpecies[match(substring(df.all$Serial.No, 1, 6), metadata$Serial.No)]

#Read in UNITE blastn results
unite <- read.csv("data/unite_blast_otus.tsv", sep="\t", header=FALSE, col.names=c("otu", "id", "title", "identity", "evalue", "bitscore"))
#Make empty dataframe for top UNITE hits and associated taxonomy
unite.class <- data.frame(OTU=unique(unite$otu), top.hit=NA, identity=NA, sp.=NA, gen.=NA, fam.=NA, ord.=NA, class.=NA, phy.=NA, king.=NA)
#Make vector of taxonomy levels
tax.levels <- c("sp.", "gen.", "fam.", "ord.", "class.", "phy.", "king.")

#For each OTU...
for (i in 1:length(unite.class$OTU)) {
  
  #Pull top hit and its identity and add to dataframe
  unite.class$identity[i] <- unite$identity[unite$otu == unite.class$OTU[i]][which.max(unite$identity[unite$otu == unite.class$OTU[i]])]
  unite.class$top.hit[i] <- str_split(unite$title[unite$otu == unite.class$OTU[i]][which.max(unite$identity[unite$otu == unite.class$OTU[i]])], pattern="\\|")[[1]][1]
  
  #For each taxonomy level...
  for (j in 1:length(tax.levels)) {
    
    #Add taxonomy data from UNITE
    unite.class[i, tax.levels[j]] <- gsub(paste0('^.*', substring(tax.levels, 1, 1)[j], '__\\s*|\\s*;.*$'), "", str_split(unite$title[unite$otu == unite.class$OTU[i]][which.max(unite$identity[unite$otu == unite.class$OTU[i]])], pattern="\\|")[[1]][5])
  }
  
  #Replace taxonomy level classification based on following thresholds:
  ## >99% identity = species
  ## >= 98 < 99 = genus
  ## >= 96 < 98 = family
  ##  >= 94 < 96 = order
  ## >= 92 < 94 = class
  ## <92 = phylum
  if (unite.class$identity[i] >= 98 & unite.class$identity[i] < 99) {
    unite.class$sp.[i] <- paste(gsub('^.*g__\\s*|\\s*;.*$', '', str_split(unite$title[unite$otu == unite.class$OTU[i]][which.max(unite$identity[unite$otu == unite.class$OTU[i]])], pattern="\\|")[[1]][5]), "sp.")
  }
  if (unite.class$identity[i] >= 96 & unite.class$identity[i] < 98) {
    for (j in 1:2) {
      unite.class[i, tax.levels[j]] <- paste(gsub('^.*f__\\s*|\\s*;.*$', "", str_split(unite$title[unite$otu == unite.class$OTU[i]][which.max(unite$identity[unite$otu == unite.class$OTU[i]])], pattern="\\|")[[1]][5]), tax.levels[j])
    }
  }
  if (unite.class$identity[i] >= 94 & unite.class$identity[i] < 96) {
    for (j in 1:3) {
      unite.class[i, tax.levels[j]] <- paste(gsub('^.*o__\\s*|\\s*;.*$', "", str_split(unite$title[unite$otu == unite.class$OTU[i]][which.max(unite$identity[unite$otu == unite.class$OTU[i]])], pattern="\\|")[[1]][5]), tax.levels[j])
    }
  } 
  if (unite.class$identity[i] >= 92 & unite.class$identity[i] < 94) {
    for (j in 1:4) {
      unite.class[i, tax.levels[j]] <- paste(gsub('^.*c__\\s*|\\s*;.*$', "", str_split(unite$title[unite$otu == unite.class$OTU[i]][which.max(unite$identity[unite$otu == unite.class$OTU[i]])], pattern="\\|")[[1]][5]), tax.levels[j])
    }
  } 
  if (unite.class$identity[i] < 92) {
    for (j in 1:5) {
      unite.class[i, tax.levels[j]] <- paste(gsub('^.*p__\\s*|\\s*;.*$', "", str_split(unite$title[unite$otu == unite.class$OTU[i]][which.max(unite$identity[unite$otu == unite.class$OTU[i]])], pattern="\\|")[[1]][5]), tax.levels[j])
    }
  }
  
  #If UNITE taxonomy contains 'unidentified'...
  if (length(grep("unidentified", unite.class[i,])) > 0) {
    #Correct the classification from higher levels
    for (j in grep("unidentified", unite.class[i,])) {
      unite.class[i, j] <- paste(unite.class[i, max(grep("unidentified", unite.class[i,])) + 1],
                                 colnames(unite.class)[j])
    }
    
  }
  
}

#Remove underscores from species names
unite.class$sp. <- gsub("_", " ", unite.class$sp.)
#Make sp. endings uniform
unite.class$sp. <- sub(" sp$", " sp.", unite.class$sp.)


##Ascomycota

#Read in OTU assignments from T-BAS
tbas.class.asco <- read.csv("data/TBAS_assignments_asco_081220.csv")
#Adapt taxon level naming to match UNITE dataframe
colnames(tbas.class.asco)[c(3, 5, 7, 9, 11)] <- c("phy.", "class.", "ord.", "fam.", "gen.")
#Filter UNITE dataframe for only Ascomycota
unite.class.asco <- unite.class[unite.class$phy. == "Ascomycota",]
#Combine results for Supplementary Table
asco.table <- cbind(unite.class.asco, tbas.class.asco[match(unite.class.asco$OTU, tbas.class.asco$Query.sequence), c(13:11, 9, 7, 5, 3)])
#Add column for classification clash between UNITE and T-BAS
unite.class.asco$clash <- NA

#For each ascomycete OTU...
for (i in 1:length(unite.class.asco$OTU)) {
  
  classified.level <- min(which(sapply(str_split(unite.class.asco[i,4:9], " "), length) == 1))
  
  #If the assignments don't match...
  if (!unite.class.asco[i, tax.levels[classified.level]] %in% tbas.class.asco[tbas.class.asco$Query.sequence == unite.class.asco$OTU[i], tax.levels[classified.level]]) {
    #Add to the clash column
    unite.class.asco$clash[i] <- "Y"
    
    min.level <- min(which(!is.na(match(unite.class.asco[i,], tbas.class.asco[tbas.class.asco$Query.sequence == unite.class.asco$OTU[i],])))[-1] - 3)
    
    for (j in 1:min.level - 1) {
      unite.class.asco[i, tax.levels[j]] <- paste(gsub(paste0('^.*', substring(tax.levels[min.level], 1, 1), '__\\s*|\\s*;.*$'), "", str_split(unite$title[unite$otu == unite.class.asco$OTU[i]][which.max(unite$identity[unite$otu == unite.class.asco$OTU[i]])], pattern="\\|")[[1]][5]), tax.levels[j])
    }
    
  }
}

#Update supplementary table with clashes and consensus
asco.table <- cbind(asco.table, unite.class.asco[match(asco.table$OTU, unite.class.asco$OTU), c(11, 4:10)])
colnames(asco.table) <- c(colnames(asco.table)[1:12], "gen.tbas", "fam.tbas", "ord.tbas", "class.tbas", "phy.tbas", "clash", "sp.consensus", "gen.consensus", "fam.consensus", "ord.consensus", "class.consensus", "phy.consensus", "king.consensus")

##Non-Ascomycota

#Read in OTU assignments from T-BAS
tbas.class.nonasco <- read.csv("data/TBAS_assignments_allfungi_081220.csv")
#Filter for non-ascomycetes
tbas.class.nonasco <- tbas.class.nonasco[tbas.class.nonasco$Phylum.level.assignment != "Ascomycota",]
#Adapt taxon level naming to match UNITE dataframe
colnames(tbas.class.nonasco)[c(3, 5)] <- c("phy.", "ord.")
#Filter UNITE dataframe for non-Ascomycota
unite.class.nonasco <- unite.class[unite.class$phy. != "Ascomycota",]
#Combine results for Supplementary Table
nonasco.table <- cbind(unite.class.nonasco, tbas.class.nonasco[match(unite.class.nonasco$OTU, tbas.class.nonasco$Query.sequence), c(7:5, 3)])
#Add column for classification clash between UNITE and T-BAS
unite.class.nonasco$clash <- NA

#For each non-ascomycete OTU...
for (i in 1:length(unite.class.nonasco$OTU)) {
  
  if (1 %in% sapply(str_split(unite.class.nonasco[i,4:9], " "), length)) {
    classified.level <- min(which(sapply(str_split(unite.class.nonasco[i,4:9], " "), length) == 1))
    
    #If the assignments don't match...
    if (!unite.class.nonasco[i, tax.levels[classified.level]] %in% tbas.class.nonasco[tbas.class.nonasco$Query.sequence == unite.class.nonasco$OTU[i], tax.levels[classified.level]]) {
      #Add to the clash column
      unite.class.nonasco$clash[i] <- "Y"
      
      if (length(which(!is.na(match(unite.class.nonasco[i,], tbas.class.nonasco[tbas.class.nonasco$Query.sequence == unite.class.nonasco$OTU[i],])))) > 0) {
        min.level <- min(which(!is.na(match(unite.class.nonasco[i,], tbas.class.nonasco[tbas.class.nonasco$Query.sequence == unite.class.nonasco$OTU[i],])))[-1] - 3)
      } else {
        min.level <- 7
      }
      
      for (j in 1:min.level - 1) {
        unite.class.nonasco[i, tax.levels[j]] <- paste(gsub(paste0('^.*', substring(tax.levels[min.level], 1, 1), '__\\s*|\\s*;.*$'), "", str_split(unite$title[unite$otu == unite.class.nonasco$OTU[i]][which.max(unite$identity[unite$otu == unite.class.nonasco$OTU[i]])], pattern="\\|")[[1]][5]), tax.levels[j])
      }
      
    }
  }
}

classification.df <- rbind(unite.class.asco, unite.class.nonasco)

#Update supplementary table with clashes and consensus
nonasco.table <- cbind(nonasco.table, unite.class.nonasco[match(nonasco.table$OTU, unite.class.nonasco$OTU), c(11, 4:10)])
colnames(nonasco.table) <- c(colnames(nonasco.table)[1:12], "ord.tbas", "phy.tbas", "clash", "sp.consensus", "gen.consensus", "fam.consensus", "ord.consensus", "class.consensus", "phy.consensus", "king.consensus")
write.csv(rbind.fill(asco.table, nonasco.table), "Supplementary_Table_3.csv", row.names=FALSE)


#Read in ITS OTU mapping from USEARCH
otus <- read.table("data/USEARCH_OTUS.txt", sep='\t', col.names=c("Serial.No","OTU"))

#Make main OTU dataframe
df <- df.all
#Add OTU dataframe
df["OTU"] <- otus$OTU[match(df$Serial.No, otus$Serial.No)]
#Add classification to dataframe
df[,c("Phylum","Class","Order","Genus", "Species")] <- classification.df[,c("phy.","class.","ord.","gen.", "sp.")][match(df$OTU, classification.df$OTU),]
#Remove accessions which couldn't cluster into OTUs
df <- subset(df, OTU != "#N/A")
#Remove empty rows
df <- df[!apply(is.na(df) | df == "", 1, all),]
#Remove duplicate clones
df <- subset(df, !duplicated(data.frame(Serial.No=sub("\\--.*", "", df$Serial.No), OTU=df$OTU)))
#Remove controls
df <- df[df$Direct.Sequencing.or.Culture != "Control",]
#Read in sequences flagged as mixed by GenBank
flagged <- read.csv("data/GenBank_flagged.tsv", header=FALSE)$V1
#Remove flagged sequences
df <- df[-which(!is.na(match(df$Serial.No, flagged))),]


#############################################################
## FIGURE 1 - OTUs PER SEED AND SPECIES ACCUMULATION CURVE ##
#############################################################

## FIGURE 1A - UNIQUE OTUs PER SEED ##

for (method in c("Culture", "Direct")) {
  
  #Extract vector of seeds with OTUs
  seed <- substring(df$Serial.No[df$Direct.Sequencing.or.Culture == method],1,9)[!duplicated(substring(df$Serial.No[df$Direct.Sequencing.or.Culture == method], 1, 9))]
  
  seed.df <- data.frame(Serial.No=seed, count=NA)
  
  #Count number of fungi per seed (all Musa accessions)
  for (i in 1:length(seed)) {
    seed.df$count[i] <- length(grep(seed[i], df$Serial.No))
  }
  
  seed.count <- data.frame(x=sort(unique(seed.df$count)), y=NA)
  
  for (i in 1:length(seed.count$x)) {
    seed.count$y[i] <- length(seed.df$count[seed.df$count == seed.count$x[i]])
  }
  
  #Count number of fungi per seed (by Musa species)
  
  #Add Musa species column to seed dataframe
  seed.df["musa"] <- metadata$MusaSpecies[match(gsub("-.*", "", seed.df$Serial.No), metadata$Serial.No)]
  #To species level
  seed.df$musa <- word(seed.df$musa, 1, 2)
  
  #Make vector of Musa species
  musa <- sort(unique(word(unique(seed.df$musa), 1, 2)))
  
  #Make dataframe to add counts to
  seed.count <- data.frame(x=sort(unique(seed.df$count)), Macu=NA, Mbal=NA, Miti=NA, Mvel=NA, Mvio=NA, method=method)
  
  for (i in 1:length(musa)) {
    for (j in 1:length(seed.count$x)) {
      seed.count[j, i+1] <- length(seed.df$count[seed.df$count == seed.count$x[j] & seed.df$musa == musa[i]])
    }
  }
  
  assign(paste0("seed.count.", method), seed.count)
  assign(paste0("seed.df.", method), seed.df)
  
}

#Combine culture and direct seed counts
seed.count <- rbind(seed.count.Culture, seed.count.Direct)

#Create dataframe for total number of seeds per species
df.tot <- data.frame(musa=rep(musa, each=2), x=0, y=NA, method=rep(c("Culture", "Direct"), times=length(musa)))

#Add total seed count
for (i in df.tot$musa) {
  df.tot[df.tot$method == "Culture" & df.tot$musa == i, "y"] <- length(unique(substring(df.all$Serial.No[grep(i, df.all$MusaSpecies)][df.all$Direct.Sequencing.or.Culture[grep(i, df.all$MusaSpecies)] == "Culture"], 1, 9)))
  df.tot[df.tot$method == "Direct" & df.tot$musa == i, "y"] <- length(unique(substring(df.all$Serial.No[grep(i, df.all$MusaSpecies)][df.all$Direct.Sequencing.or.Culture[grep(i, df.all$MusaSpecies)] == "Direct"], 1, 9)))
}

#Remove numbers of seeds with isolates from the total seed count
df.tot$y[df.tot$method == "Culture"] <- df.tot$y[df.tot$method == "Culture"] - colSums(seed.count.Culture[2:6])
df.tot$y[df.tot$method == "Direct"] <- df.tot$y[df.tot$method == "Direct"] - colSums(seed.count.Direct[2:6])

#Combine counts for all accessions and for separate Musa species into one dataframe for plotting
seed.df <- rbind(df.tot, 
                 data.frame(x=seed.count$x, y=seed.count$Macu, musa=musa[1], method=seed.count$method),
                 data.frame(x=seed.count$x, y=seed.count$Mbal, musa=musa[2], method=seed.count$method),
                 data.frame(x=seed.count$x, y=seed.count$Miti, musa=musa[3], method=seed.count$method),
                 data.frame(x=seed.count$x, y=seed.count$Mvel, musa=musa[4], method=seed.count$method),
                 data.frame(x=seed.count$x, y=seed.count$Mvio, musa=musa[5], method=seed.count$method),
                 data.frame(x=0, y=rep(20, 2), musa="Musa gracilis", method=c("Culture", "Direct")))

#Custom order of Musa species in plot
seed.df$musa <- factor(seed.df$musa, levels=c("Musa acuminata", "Musa balbisiana", "Musa itinerans", "Musa velutina", "Musa violascens", "Musa gracilis"))
#Replace labels
seed.df$method <- sub("Culture", "Culture-dependent", seed.df$method)
seed.df$method <- sub("Direct", "Culture-independent", seed.df$method)

#Plot seed count bargraphs
gg.seed <- ggplot(seed.df, aes(x=x, y=y)) +
  geom_bar(stat="identity") +
  geom_text(data=subset(seed.df, y > 0), aes(label=y), nudge_y=15, size=2) +
  facet_grid(method ~ musa) +
  scale_x_discrete(name="Number of unique OTUs",
                   limits=c(seq(0,7,1))) +
  scale_y_continuous(name="Number of seeds",
                     limits=c(0, 256),
                     breaks=seq(0, 255, 50),
                     expand=c(0, 5)) +
  theme_minimal() +
  theme(axis.title.x=element_text(margin=unit(c(3, 0, 0, 0), "mm")),
        axis.title.y=element_text(margin=unit(c(0, 3, 0, 0), "mm")),
        panel.spacing.y=unit(1, "lines"),
        strip.text.x=element_text(face="italic", size=7),
        strip.text.y=element_text(size=9),
        plot.title.position="plot") +
  labs(subtitle=expression(bold("A")))


## FIGURE 1B - SPECIES ACCUMULATION CURVE ##

#Make dataframe of species
spec.df <- metadata[c("Serial.No", "MusaSpecies")]

#Create OTU abundance matrix
otu.count <- matrix(nrow=length(unique(spec.df$Serial.No)),ncol=length(unique(df$OTU)))
colnames(otu.count) <- unique(df$OTU)
rownames(otu.count) <- unique(spec.df$Serial.No)

for (i in 1:length(rownames(otu.count))) {
  for (j in 1:length(colnames(otu.count))) {
    otu.count[i,j] <- length(df$OTU[grep(rownames(otu.count)[i], df$Serial.No)][df$OTU[grep(rownames(otu.count)[i], df$Serial.No)] == colnames(otu.count)[j]])
  }
}

#Accumulation curves including singletons

#Create dataframe to collect rarefaction accumulation curve data for acuminata, balbisians and itinerans
spec.accum.df <- data.frame(individual=NA, richness=NA, error=NA, musa=sort(spec.df$MusaSpecies[spec.df$MusaSpecies == musa[1] | spec.df$MusaSpecies == musa[2] | spec.df$MusaSpecies == musa[3]]), rare="Including singletons")

#Rarefaction for acuminata, balbisians and itinerans
for (i in 1:length(musa[1:3])) {
  #Subset OTU abundance matrix for top three host species
  tmp.df <- subset(otu.count, rownames(otu.count) %in% rownames(otu.count)[spec.df$MusaSpecies[match(rownames(otu.count), spec.df$Serial.No)] == musa[i]])
  #Species accumulation
  tmp.accum <- specaccum(tmp.df, method="rarefaction")
  #Add to results dataframe
  spec.accum.df[spec.accum.df$musa == musa[i],] <- data.frame(individual=tmp.accum$sites, richness=tmp.accum$richness, error=tmp.accum$sd, musa=musa[i], rare="Including singletons")
}

#Rarefaction for all Musa accessions
tmp.accum <- specaccum(otu.count, method="rarefaction")
tmp.accum.df <- data.frame(individual=tmp.accum$sites, richness=tmp.accum$richness, error=tmp.accum$sd, musa="All", rare="Including singletons")

#Combine species accumulation data including singletons
spec.accum.df <- rbind(spec.accum.df, tmp.accum.df)


#Accumulation curves excluding singletons

#Create dataframe to collect rarefaction accumulation curve data for acuminata, balbisians and itinerans
spec.accum.excl.df <- data.frame(individual=NA, richness=NA, error=NA, musa=sort(spec.df$MusaSpecies[spec.df$MusaSpecies == musa[1] | spec.df$MusaSpecies == musa[2] | spec.df$MusaSpecies == musa[3]]), rare="Excluding singletons")

#Remove singleton OTUs
otu.count.excl <- subset(otu.count, select=-which(colSums(otu.count) == 1))

#Rarefaction for acuminata, balbisians and itinerans
for (i in 1:length(musa[1:3])) {
  #Subset for top three host species
  tmp.df <- subset(otu.count.excl, rownames(otu.count.excl) %in% rownames(otu.count.excl)[spec.df$MusaSpecies[match(rownames(otu.count.excl), spec.df$Serial.No)] == musa[i]])
  #Species accumulation
  tmp.accum <- specaccum(tmp.df, method="rarefaction")
  #Add to results dataframe
  spec.accum.excl.df[spec.accum.excl.df$musa == musa[i],] <- data.frame(individual=tmp.accum$sites, richness=tmp.accum$richness, error=tmp.accum$sd, musa=musa[i], rare="Excluding singletons")
}

#Rarefaction for all Musa accessions
tmp.accum <- specaccum(otu.count.excl, method="rarefaction")
tmp.accum.df <- data.frame(individual=tmp.accum$sites, richness=tmp.accum$richness, error=tmp.accum$sd, musa="All", rare="Excluding singletons")

#Combine species accumulation data including and excluding singletons 
spec.accum.df <- rbind(spec.accum.excl.df, tmp.accum.df, spec.accum.df)

#Plot accumulation curves
gg.accum <- ggplot(spec.accum.df, aes(x=individual, y=richness, color=musa)) +
  facet_wrap(. ~ rare, scales="free", nrow=1) +
  geom_ribbon(aes(ymax=richness+error, ymin=richness-error), linetype=0, fill="grey", alpha=0.4, show.legend=FALSE) +
  geom_line(size=0.5) +
  scale_color_manual(values=c("black", "#66c2a5","#fc8d62","#8da0cb"),
                     name=expression(paste(italic("Musa")," species")),
                     labels=c("All", expression(italic("Musa acuminata")), expression(italic("Musa balbisiana")), expression(italic("Musa itinerans")))) +
  labs(x=expression(paste("Number of ",italic("Musa")," accessions sampled")), y="Number of unique OTUs") +
  theme_minimal() +
  theme(axis.title.x=element_text(margin=unit(c(3, 0, 0, 0), "mm")),
        axis.title.y=element_text(margin=unit(c(0, 3, 0, 0), "mm")),
        legend.text.align=0,
        plot.title.position="plot") +
  labs(subtitle=expression(bold("B")))


#Write to file - FIGURE 1
tiff(file=paste0("Figure_1-", Sys.Date(), ".tiff"), height=18, width=15, units="cm", res=300)
ggpubr::ggarrange(gg.seed, gg.accum, nrow=2, heights=c(1.5,1))
dev.off()


###################################
## FIGURE 2 - OTU CLASSIFICATION ##
###################################

## T-BAS VISUALISATION OF OTU ASSIGNMENT ##

#Make dataframe of OTU count per isolation method
method.count <- data.frame(otu=unique(df$OTU), direct=NA, culture=NA, clone=NA)

for (i in 1:length(unique(df$OTU))) {
  method.count$direct[i] <- sum(df$Direct.Sequencing.or.Culture[df$OTU == unique(df$OTU)[i]] == "Direct" & df$Clone[df$OTU == unique(df$OTU)[i]] != "Y")
  method.count$culture[i] <- sum(df$Direct.Sequencing.or.Culture[df$OTU == unique(df$OTU)[i]] == "Culture")
  method.count$clone[i] <- sum(df$Clone[df$OTU == unique(df$OTU)[i]] == "Y")
}

#Read in T-BAS tree
tbas <- read.tree("data/TBAS_asco_081220.tre")
#Make dataframe of isolation method
tbas.df <- data.frame(tip=tbas$tip.label, direct=NA, culture=NA, clone=NA)
#Add counts to T-BAS plot dataframe
tbas.df$direct <- method.count$direct[match(tbas.df$tip, method.count$otu)]
tbas.df$culture <- method.count$culture[match(tbas.df$tip, method.count$otu)]
tbas.df$clone <- method.count$clone[match(tbas.df$tip, method.count$otu)]
#Replace 0 with NA
tbas.df[tbas.df == 0] <- NA
#Read in node and colour data for T-BAS tree
cladelabels <- read.csv("data/TBAS_clades.csv")

#Plot tree
gg.tbas <- ggtree(tbas,
                  layout="circular",
                  cex=0.2,
                  branch.length="none")

for (i in which(!is.na(cladelabels$node))) {
  gg.tbas <- gg.tbas +
    geom_hilight(node=cladelabels[i,2], fill=cladelabels[i,3])
}

gg.tbas <- gg.tbas %<+% tbas.df +
  geom_tree(cex=0.2) +
  geom_tippoint(aes(x=x+1, size=culture, colour="Culturing"), na.rm=TRUE) +    
  geom_tippoint(aes(x=x+3, size=direct, colour="Direct extraction"), na.rm=TRUE) +
  geom_tippoint(aes(x=x+3, size=clone, colour="Cloning"), alpha=0.4, na.rm=TRUE) +
  scale_colour_manual(name="Sampling method", labels=c("Cloning", "Culturing","Direct sequencing"), values=c("black", "#E69F00", "#56B4E9")) +
  labs(size="OTU abundance") +
  guides(colour=guide_legend(override.aes=list(size=3))) +
  theme(legend.position="right",
        plot.margin=unit(c(-10,0,-10,0), "mm"),
        legend.box.margin=margin(c(0,0,0,0)),
        plot.title.position="plot")


## PIECHART VISUALISATION OF OTU CLASSIFICATION ##

#Make taxonomy dataframe
tmp <- data.frame(phylum=df$Phylum, class=df$Class, order=df$Order, genus=df$Genus, species=df$Species)
#Filter for Ascomycota
pie.df <- tmp[tmp$phylum == "Ascomycota",]
#Remove duplicate rows for genera
pie.df <- pie.df[!duplicated(pie.df$species),]

#Add column of total count for each genus
pie.df["num"] <- NA

for (i in 1:length(pie.df$species)) {
  pie.df$num[i] <- sum(df$Species == pie.df$species[i])
}

#Sort dataframe sequentially by class, order and genus
pie.df <- transform(pie.df, freq=ave(seq(nrow(pie.df)), class, FUN=length))
pie.df <- transform(pie.df, freq2=ave(seq(nrow(pie.df)), order, FUN=length))
pie.df <- transform(pie.df, freq3=ave(seq(nrow(pie.df)), genus, FUN=length))
pie.df <- pie.df[order(-pie.df["freq"], -pie.df["freq2"], -pie.df["freq3"], -pie.df["num"], pie.df["species"]), ]

#Add ymin and ymax columns for size of pie slice
ymin <- c(0,pie.df$num)

for (i in 2:length(ymin)) {
  ymin[i] <- ymin[i] + ymin[i-1]
}

pie.df["ymin"] <- ymin[1:length(pie.df$species)]
pie.df["ymax"] <- ymin[2:length(ymin)]

#Add column with midpoint position for genus label in pie slice
pie.df["lab.pos.sp"] <- (pie.df$ymax + pie.df$ymin) / 2


for (i in c("genus", "order", "class")) {
  
  #Make dataframe of genus label positions
  pos <- data.frame(taxon=unique(pie.df[,i]), lab.pos=NA)
  
  for (j in 1:length(unique(pie.df[,i]))){
    pos$lab.pos[pos$taxon == unique(pie.df[,i])[j]] <- pie.df$ymax[pie.df[,i] == unique(pie.df[,i])[j]][which.max(pie.df$ymax[pie.df[,i] == unique(pie.df[,i])[j]])]
  }
  
  ymin <- c(0, pos$lab.pos)
  pos["ymin"] <- ymin[1:length(pos$taxon)]
  pos["ymax"] <- ymin[2:length(ymin)]
  
  res <- vector(length=length(pos$lab.pos))
  
  for (k in 2:length(ymin)) {
    res[k-1] <- (ymin[k] + ymin[k-1]) / 2
  }
  
  pos$lab.pos <- res
  
  assign(paste0("pos.", i), pos)
  
}

#Add unclassified Ascomycota to the T-BAS clade colours
colours <- rbind(cladelabels[-2], data.frame(class="Ascomycota class.", colour="#818181", mid="#9b9b9b", low="#b5b5b5", lower="#b5b5b5", lowest="#c7c7c7"))
#Sort the T-BAS clade colours according to the count data
colours <- colours[match(unique(pie.df$class), colours$class),]

#Add colour columns to main dataframe
for (i in 1:length(unique(pie.df$class))) {
  
  #Class
  pie.df$colour.class[pie.df$class == unique(pie.df$class)[i]] <- as.character(colours$colour[i])
  
  #Order
  tmp <- subset(pie.df, class == unique(pie.df$class)[i])
  for (j in 1:length(unique(pie.df$order[pie.df$class == unique(pie.df$class)[i]]))) {
    pie.df$colour.order[pie.df$order == unique(tmp$order)[j]] <- colorRampPalette(c(as.character(colours$mid[i]),as.character(colours$low[i])))(length(unique(pie.df$order[pie.df$class == unique(pie.df$class)[i]])))[j]
  }
  tmp <- subset(pie.df, class == unique(pie.df$class)[i])
  for (j in 1:length(unique(pie.df$genus[pie.df$class == unique(pie.df$class)[i]]))) {
    pie.df$colour.genus[pie.df$genus == unique(tmp$genus)[j]] <- colorRampPalette(c(as.character(colours$low[i]),as.character(colours$lower[i])))(length(unique(pie.df$genus[pie.df$class == unique(pie.df$class)[i]])))[j]
  }
  
  #Species
  pie.df$colour.species[pie.df$class == unique(pie.df$class)[i]] <- colorRampPalette(c(as.character(colours$lower[i]),as.character(colours$lowest[i])))(length(unique(pie.df$species[pie.df$class == unique(pie.df$class)[i]])))
}

#Create a vector with names from all taxon levels and sort alphabetically
all.tax <- c(as.character(pie.df$class), as.character(pie.df$order), as.character(pie.df$genus), as.character(pie.df$species))
all.tax <- all.tax[order(all.tax)]
all.tax <- all.tax[!duplicated(all.tax)]

#Create a dataframe of colours for ggplot
pie.colours <- rbind(data.frame(name=pie.df$class[!duplicated(pie.df$class)],
                                colour=pie.df$colour.class[!duplicated(pie.df$class)]),
                     data.frame(name=pie.df$order[!duplicated(pie.df$order)],
                                colour=pie.df$colour.order[!duplicated(pie.df$order)]),
                     data.frame(name=pie.df$genus[!duplicated(pie.df$genus)],
                                colour=pie.df$colour.genus[!duplicated(pie.df$genus)]),
                     data.frame(name=pie.df$species[!duplicated(pie.df$species)],
                                colour=pie.df$colour.species[!duplicated(pie.df$species)]))
#Sorted alphabetically and by taxonomic level
pie.colours <- pie.colours[order(pie.colours$name),]
#Make into vector
pie.colours <- as.vector(pie.colours$colour[match(all.tax, pie.colours$name)])
#Edit Fusarium solani name
pie.df$species[grep("Fusarium solani", pie.df$species)] <- "Fusarium solani \u2020"

set.seed(1)

#Plot piechart
gg.pie <- ggplot(pie.df) +
  geom_rect(aes(fill=species, ymax=ymax, ymin=ymin, xmax=6, xmin=4.5), 
            colour="white",
            show.legend=FALSE) +
  geom_rect(data=pos.genus,
            aes(fill=taxon, ymax=ymax, ymin=ymin, xmax=4.5, xmin=3), 
            colour="white",
            show.legend=FALSE) +
  geom_rect(data=pos.order, 
            aes(fill=taxon, ymax=ymax, ymin=ymin, xmax=3, xmin=1.5), 
            colour="white",
            show.legend=FALSE) +
  geom_rect(data=pos.class,
            aes(fill=taxon, ymax=ymax, ymin=ymin, xmax=1.5, xmin=0), 
            colour="white",
            show.legend=FALSE) +
  annotate("text", x=0.75, y=max(pie.df$ymax)*0.75,
           label="Class", color="white", size=3, fontface="bold") +
  annotate("text", x=2.25, y=max(pie.df$ymax)*0.75,
           label="Order", color="white", size=3, fontface="bold") +
  annotate("text", x=3.75, y=max(pie.df$ymax)*0.75,
           label="Genus", color="white", size=3, fontface="bold") +
  annotate("text", x=5.25, y=max(pie.df$ymax)*0.75,
           label="Species", color="white", size=3, fontface="bold") +
  geom_label_repel(aes(x=6, y=lab.pos.sp, label=species, fill=species),
                   #direction="y",
                   fontface="italic", 
                   size=2, 
                   nudge_x=1,
                   show.legend=FALSE) +
  geom_label_repel(data=pos.order[-c(grep("ord", pos.order$taxon), 15),], 
                   aes(x=2.25, y=lab.pos, label=taxon, fill=taxon),
                   direction="y",
                   size=2,
                   nudge_x=0.5,
                   show.legend=FALSE) +
  geom_point(data=pos.class,
             aes(x=1, y=lab.pos, color=taxon),
             size=2,
             shape=15,
             alpha=0) +
  scale_fill_manual(values=pie.colours) +
  scale_color_manual(name="", labels=c("Unclassified Ascomycota",as.character(sort(pos.class$taxon)[-1])), values=as.vector(colours$colour[order(as.vector(colours$class))])) +
  xlim(c(0, 8)) +
  coord_polar(theta="y") +
  theme_void() +
  guides(color=guide_legend(override.aes=list(alpha=1, size=5))) +
  theme(legend.position="left",
        plot.margin=unit(c(-10,0,-10,0), "mm"),
        legend.box.margin=margin(c(0,-80,0,0)))


#Write to file - FIGURE 2
tiff(file=paste0("Figure_2-", Sys.Date(), ".tiff"), height=26, width=20, units="cm", res=300)
ggpubr::ggarrange(gg.tbas, gg.pie, heights=c(1,1.5), ncol=1)
dev.off()


################################################
## FIGURE 3 - COMPARISON OF DETECTION METHODS ##
################################################

#Test for difference in OTUs between methods
#Create OTU abundance matrix with accessions split into three methods
otu.count.methods <- matrix(nrow=length(unique(spec.df$Serial.No))*3,ncol=length(unique(df$OTU)))
colnames(otu.count.methods) <- unique(df$OTU)
rownames(otu.count.methods) <- c(paste0(unique(spec.df$Serial.No), ".culture"), 
                                 paste0(unique(spec.df$Serial.No), ".direct"), 
                                 paste0(unique(spec.df$Serial.No), ".clone"))

#For each row...
for (i in 1:length(rownames(otu.count.methods))) {
  
  #For each OTU...
  for (j in 1:length(colnames(otu.count.methods))) {
    
    #Grab the accession
    accession <- sub("\\..*", "", rownames(otu.count.methods)[i])
    
    #If the method is culturing, count the number of occurrences of the OTU and enter into matrix
    if (sub(".*\\.", "", rownames(otu.count.methods)[i]) == "culture") {
      otu.count.methods[i,j] <- length(df$OTU[grep(accession, df$Serial.No)[which(grepl(accession, df$Serial.No)) %in% which(df$Direct.Sequencing.or.Culture == "Culture")]][df$OTU[grep(accession, df$Serial.No)[which(grepl(accession, df$Serial.No)) %in% which(df$Direct.Sequencing.or.Culture == "Culture")]] == colnames(otu.count.methods)[j]])
    }
    #If the method is direct sequencing, count the number of occurrences of the OTU and enter into matrix
    if (sub(".*\\.", "", rownames(otu.count.methods)[i]) == "direct") {
      otu.count.methods[i,j] <- length(df$OTU[grep(accession, df$Serial.No)[which(grepl(accession, df$Serial.No)) %in% which(df$Direct.Sequencing.or.Culture == "Direct" & df$Clone != "Y")]][df$OTU[grep(accession, df$Serial.No)[which(grepl(accession, df$Serial.No)) %in% which(df$Direct.Sequencing.or.Culture == "Direct" & df$Clone != "Y")]] == colnames(otu.count.methods)[j]])
    }
    #If the method is cloning, count the number of occurrences of the OTU and enter into matrix
    if (sub(".*\\.", "", rownames(otu.count.methods)[i]) == "clone") {
      otu.count.methods[i,j] <- length(df$OTU[grep(accession, df$Serial.No)[which(grepl(accession, df$Serial.No)) %in% which(df$Clone == "Y")]][df$OTU[grep(accession, df$Serial.No)[which(grepl(accession, df$Serial.No)) %in% which(df$Clone == "Y")]] == colnames(otu.count.methods)[j]])
    }
  }
}

#Remove rows with no occurrences
otu.count.methods <- otu.count.methods[-which(rowSums(otu.count.methods) == 0),]

#betadisper by methods
anova(betadisper(vegdist(otu.count.methods, method="bray"), sub(".*\\.", "", rownames(otu.count.methods))))

#ANOSIM by methods
method.anosim <- anosim(otu.count.methods, sub(".*\\.", "", rownames(otu.count.methods)), distance="bray", permutations=999)

#Indices for OTUs found by different combinations of detection methods
#Direct
direct <- which(method.count$direct > 0 & method.count$culture == 0 & method.count$clone == 0)
#Culture
culture <- which(method.count$direct == 0 & method.count$culture > 0 & method.count$clone == 0)
#Clone
clone <- which(method.count$direct == 0 & method.count$culture == 0 & method.count$clone > 0)
#Direct&Culture
direct.culture <- which(method.count$direct > 0 & method.count$culture > 0 & method.count$clone == 0)
#Direct&Clone
direct.clone <- which(method.count$direct > 0 & method.count$culture == 0 & method.count$clone > 0)
#Culture&Clone
culture.clone <- which(method.count$direct == 0 & method.count$culture > 0 & method.count$clone > 0)
#All
direct.culture.clone <- which(method.count$direct > 0 & method.count$culture > 0 & method.count$clone > 0)

#Make matrix of number of OTUs for each detection method
method.otus <- c("Direct"=length(method.count$otu[direct]),
                 "Culture"=length(method.count$otu[culture]),
                 "Clone"=length(method.count$otu[clone]),
                 "Direct&Culture"=length(method.count$otu[direct.culture]),
                 "Direct&Clone"=length(method.count$otu[direct.clone]),
                 "Culture&Clone"=length(method.count$otu[culture.clone]),
                 "Direct&Culture&Clone"=length(method.count$otu[direct.culture.clone]))





#Create and plot euler diagram
set.seed(1)
method.euler <- euler(method.otus, shape="ellipse")
euler <- plot(method.euler,
              fills=list(fill=c("#56B4E9", "#E69F00", "black"), alpha=0.3),
              labels=FALSE,
              #quantities=list(col="white"),
              edges=FALSE)
#Convert to ggplot
gg.euler <- as.ggplot(euler)
#Create dataframe for OTU wordcloud labels
euler.df <- rbind(data.frame(otu=df$Species[match(method.count$otu[direct], df$OTU)],
                             count=rowSums(method.count[direct, 2:4]),
                             method="direct",
                             x=0.3, y=0.25),
                  data.frame(otu=df$Species[match(method.count$otu[culture], df$OTU)],
                             count=rowSums(method.count[culture, 2:4]),
                             method="culture",
                             x=0.7, y=0.25),
                  data.frame(otu=df$Species[match(method.count$otu[clone], df$OTU)],
                             count=rowSums(method.count[clone, 2:4]),
                             method="clone",
                             x=0.52, y=0.8),
                  data.frame(otu=df$Species[match(method.count$otu[direct.culture], df$OTU)],
                             count=rowSums(method.count[direct.culture, 2:4]),
                             method="direct.culture",
                             x=0.52, y=0.34),
                  data.frame(otu=df$Species[match(method.count$otu[direct.clone], df$OTU)],
                             count=rowSums(method.count[direct.clone, 2:4]),
                             method="direct.clone",
                             x=0.42, y=0.59),
                  data.frame(otu=df$Species[match(method.count$otu[culture.clone], df$OTU)],
                             count=rowSums(method.count[culture.clone, 2:4]),
                             method="culture.clone",
                             x=0.63, y=0.51),
                  data.frame(otu=df$Species[match(method.count$otu[direct.culture.clone], df$OTU)],
                             count=rowSums(method.count[direct.culture.clone, 2:4]),
                             method="direct.culture.clone",
                             x=0.5, y=0.49))

#Plot euler
gg.euler <- gg.euler +
  geom_text(x=0.1, y=0.5, size=4, fontface="bold", colour="#56B4E9", 
            label=paste0("Direct sequencing\n", length(method.count$otu[direct]))) +
  geom_text(x=0.85, y=0.5, colour="#E69F00", size=4, fontface="bold",
            label=paste0("Culture\n", length(method.count$otu[culture]))) +
  geom_text(x=0.3, y=0.9, colour="black", size=4, fontface="bold",
            label=paste0("Clone\n", length(method.count$otu[clone]))) +
  geom_text(x=0.38, y=0.5, colour="white",
            size=5, fontface="bold",
            label=length(method.count$otu[direct.clone])) +
  geom_text(x=0.515, y=0.285, colour="white",
            size=5, fontface="bold",
            label=length(method.count$otu[direct.culture])) +
  geom_text(x=0.6, y=0.58, colour="white",
            size=5, fontface="bold",
            label=length(method.count$otu[culture.clone])) +
  geom_text(x=0.5, y=0.59, colour="white",
            size=5, fontface="bold",
            label=length(method.count$otu[direct.culture.clone])) +
  geom_text(x=0.85, y=0.9, size=3, fontface="bold",
            label=paste0("ANOSIM")) +
  geom_text(x=0.85, y=0.85, size=3, 
            label=paste0("R=", round(method.anosim$statistic, 2), " p=", round(method.anosim$signif, 3))) +
  geom_text_wordcloud(data=euler.df,
                      eccentricity=1,
                      area_corr_power=1,
                      aes(label=otu, colour=method, size=count)) +
  scale_size_continuous(range=c(1, 3)) +
  scale_colour_manual(values=c("#000000", #clone
                               "#E69F00", #culture
                               hex(mixcolor(alpha=0.5, hex2RGB("#E69F00"), hex2RGB("#000000"))), #culture.clone
                               "#56B4E9", #direct
                               hex(mixcolor(alpha=0.5, hex2RGB("#56B4E9"), hex2RGB("#000000"))), #direct.clone
                               hex(mixcolor(alpha=0.5, hex2RGB("#56B4E9"), hex2RGB("#E69F00"))), #direct.culture
                               hex(mixcolor(alpha=0.5, (mixcolor(alpha=0.5, hex2RGB("#E69F00"), hex2RGB("#000000"))), hex2RGB("#56B4E9"))))) + #direct.culture.clone
  theme(legend.position="none")

#Write to file - FIGURE 3
tiff(file=paste0("Figure_3-", Sys.Date(), ".tiff"), height=12, width=18, units="cm", res=300)
gg.euler
dev.off()


################################################################
## FIGURE 4 - COMMUNITY COMPOSITION, DIVERSITY AND ABUNDANCE  ##
################################################################

## FIGURE 4A - NMDS ANALYSIS ##

#Subset OTU abundance matrix by Musa accessions for acuminata, balbisiana and itinerans
otu.count.nmds <- subset(otu.count, rownames(otu.count) %in% rownames(otu.count)[spec.df$MusaSpecies[match(rownames(otu.count), spec.df$Serial.No)] == musa[1] | spec.df$MusaSpecies[match(rownames(otu.count), spec.df$Serial.No)] == musa[2] | spec.df$MusaSpecies[match(rownames(otu.count), spec.df$Serial.No)] == musa[3]])
#Filter for most abundant OTUs
otu.count.nmds <- otu.count.nmds[,-which(colSums(otu.count.nmds) < 20)]
#Remove rows with no OTUs
otu.count.nmds <- otu.count.nmds[-which(rowSums(otu.count.nmds) == 0),]

#Function to perform NMDS for 1-10 dimensions to pick optimal number of dimensions
NMDS.scree <- function(x) {
  plot(rep(1, 10), replicate(10, metaMDS(x, autotransform=F, k=1)$stress), xlim=c(1, 10),ylim=c(0, 0.30), xlab="# of Dimensions", ylab="Stress")
  for (i in 1:10) {
    points(rep(i + 1,10),replicate(10, metaMDS(x, autotransform=F, k=i + 1)$stress))
  }
}

#Write to file scree and stress plot - SUPPLEMENTARY FIGURE 1
tiff(file=paste0("Supplementary_Figure_1-", Sys.Date(), ".tiff"), height=20, width=17, units="cm", res=300)
#Plot screeplot
par(mfrow=c(2,1), mar=c(4,4,2,2))
NMDS.scree(otu.count.nmds)

set.seed(2)
#NMDS with optimal dimensions from scree plot
NMDS <- metaMDS(otu.count.nmds, k=6, trymax=1000, trace=F)

#Plot stressplot
stressplot(NMDS)
dev.off()

#Make dataframes of NMDS site scores for ggplot (site=Musa accession)
metadata.nmds <- as.data.frame(scores(NMDS))
metadata.nmds$Serial.No <- rownames(metadata.nmds)
#Add Musa species
metadata.nmds$MusaSpecies <- spec.df$MusaSpecies[match(metadata.nmds$Serial.No, spec.df$Serial.No)]
#Add habitat
metadata.nmds$Habitat.summary <- metadata$Habitat.summary[match(metadata.nmds$Serial.No, metadata$Serial.No)]
#Add TZ
metadata.nmds$TZ <- metadata$TZ[match(metadata.nmds$Serial.No, metadata$Serial.No)]
#Add germination rate
metadata.nmds$Germination <- metadata$ER.shoots[match(metadata.nmds$Serial.No, metadata$Serial.No)] / metadata$ER.embryos[match(metadata.nmds$Serial.No, metadata$Serial.No)] * 100

#Make dataframe of NMDS species scores for ggplot (species=fungal OTUs)
otu.gg <- as.data.frame(scores(NMDS, "species"))
otu.gg$label <- df$Species[match(rownames(otu.gg), df$OTU)]
#Edit Fusarium solani name
otu.gg$label[grep("Fusarium solani", otu.gg$label)] <- "Fusarium solani \u2020"

#Fit TZ variable to NMDS
nmds.ordi.tz <- ordisurf(NMDS ~ metadata.nmds$TZ, plot=FALSE)
#Pull out coordinates for plotting
nmds.ordi.tz.gg <- expand.grid(x=nmds.ordi.tz$grid$x, y=nmds.ordi.tz$grid$y)
#Get z scores
nmds.ordi.tz.gg$z <- as.vector(nmds.ordi.tz$grid$z)
#Remova NAs
nmds.ordi.tz.gg <- data.frame(na.omit(nmds.ordi.tz.gg))

#Fit germination rate variable to NMDS
nmds.ordi.germ <- ordisurf(NMDS ~ metadata.nmds$Germination, plot=FALSE)
#Pull out coordinates for plotting
nmds.ordi.germ.gg <- expand.grid(x=nmds.ordi.germ$grid$x, y=nmds.ordi.germ$grid$y)
#Get z scores
nmds.ordi.germ.gg$z <- as.vector(nmds.ordi.germ$grid$z)
#Remova NAs
nmds.ordi.germ.gg <- data.frame(na.omit(nmds.ordi.germ.gg))


##PERMANOVA

#PERMANOVA on accessions with subset of common taxa used in NMDS

#Marginal PERMANOVA to test for unique variable impact
set.seed(1)
nmds.permanova.marginal <- as.data.frame(adonis2(otu.count.nmds ~ metadata.nmds$Habitat.summary + metadata.nmds$MusaSpecies, method="bray", perm=999, by="margin"))
#Main PERMANOVA
set.seed(1)
nmds.permanova <- as.data.frame(adonis(otu.count.nmds ~ metadata.nmds$Habitat.summary + metadata.nmds$MusaSpecies, method="bray", perm=999)$aov.tab)

#PERMANOVA on all accessions including rare taxa (excluding oil palm plantation and Kew accessions)

#Filter OTU counts for non-empty accessions
otu.count.rare <- otu.count[-which(rowSums(otu.count) == 0),]
#Exclude low-sampled habitat accessions
otu.count.rare <- subset(otu.count.rare, rownames(otu.count.rare) %in% metadata$Serial.No[-c(35, 36, 39, 45)])
#Filter metadata similarly
metadata.rare <- metadata[match(rownames(otu.count.rare), metadata$Serial.No),]
#Marginal PERMANOVA to test for unique variable impact
set.seed(1)
permanova.marginal <- as.data.frame(adonis2(otu.count.rare ~ metadata.rare$Habitat.summary + metadata.rare$MusaSpecies, method="bray", perm=999, by="margin"))
#Main PERMANOVA
set.seed(1)
permanova <- as.data.frame(adonis(otu.count.rare ~ metadata.rare$Habitat.summary + metadata.rare$MusaSpecies, method="bray", perm=999)$aov.tab)

#Combine PERMANOVA testing into dataframe
Table.1 <- merge(rbind(data.frame(variable=rownames(nmds.permanova.marginal[1:2,]),
                                  nmds.permanova.marginal[1:2, c(3,5)]),
                       data.frame(variable=rownames(permanova.marginal[1:2,]),
                                  permanova.marginal[1:2, c(3,5)])),
                 rbind(data.frame(variable=rownames(nmds.permanova[1:2,]),
                                  nmds.permanova[1:2, c(5,6)]),
                       data.frame(variable=rownames(permanova[1:2,]),
                                  permanova[1:2, c(5,6)])), by="variable", all=TRUE)


#PERMDISP

#For each of the data subsets (common and all)...
for (i in c("nmds", "rare")) {
  
  #For each variable...
  for (j in c("Habitat.summary", "MusaSpecies")) {
    
    #betadisper analysis
    disp <- betadisper(vegdist(get(paste0("otu.count.", i)), method="bray"), get(paste0("metadata.", i))[,j])
    #ANOVA for significance
    disp.anova <- anova(betadisper(vegdist(get(paste0("otu.count.", i)), method="bray"),
                                   get(paste0("metadata.", i))[,j]))
    #Print results based on p-value significance
    if (disp.anova$`Pr(>F)`[1] >= 0.05) {
      print(paste(i, j, ": even data dispersion, p=", signif(disp.anova$`Pr(>F)`[1], digits=3)))
    } else {
      print(paste(i, j, ": uneven data dispersion, p=", signif(disp.anova$`Pr(>F)`[1], digits=3)))
    }
    
    assign(paste0("disp.anova.", i, ".", j), disp.anova)
    
    #For the habitat variable...
    if (j == "Habitat.summary") {
      #Make a dataframe with the dispersion distances from centroids for each accession
      disp.df <- as.data.frame(disp$distances)
      colnames(disp.df) <- "distances"
      #Add accession
      disp.df$accession <- rownames(disp.df)
      #Add habitat
      disp.df$hab <- metadata$Habitat.summary[match(disp.df$accession, metadata$Serial.No)]
      #Add data subset
      disp.df$data <- i
      
      assign(paste0("disp.", i, ".", j), disp.df)
    }
  }
}

#Add PERMDISP results to table
Table.1$betadisper.p <- NA
Table.1$betadisper.p[1] <- disp.anova.nmds.Habitat.summary$`Pr(>F)`[1]
Table.1$betadisper.p[2] <- disp.anova.nmds.MusaSpecies$`Pr(>F)`[1]
Table.1$betadisper.p[3] <- disp.anova.rare.Habitat.summary$`Pr(>F)`[1]
Table.1$betadisper.p[4] <- disp.anova.rare.MusaSpecies$`Pr(>F)`[1]

#Write to file
write.csv(file="Table_1.csv", Table.1, row.names=FALSE)


## Supplementary Figure 2 ##
#Combine habitat dataframes
disp.df <- rbind(disp.nmds.Habitat.summary, disp.rare.Habitat.summary)
#Set order for plot
disp.df$hab <- factor(disp.df$hab, levels=c("Roadside", "Jungle edge", "Jungle buffer", "Ravines"))
#Make dataframe for plot labels
disp.lab <- rbind(cbind(as.data.frame(metadata.nmds %>% group_by(Habitat.summary) %>% dplyr::summarise(n=n())),
                        data="nmds"),
                  cbind(as.data.frame(metadata.rare %>% group_by(Habitat.summary) %>% dplyr::summarise(n=n())),
                        data="rare"))

#Plot boxplot
gg.disp <- ggplot(disp.df, aes(x=distances, y=hab)) +
  facet_wrap(. ~ data, labeller=labeller(data=c(`nmds`="Common OTUs", `rare`="All OTUs (including rare)"))) +
  geom_boxploth(aes(fill=hab)) +
  geom_text(data=disp.lab,
            aes(x=Inf, y=Habitat.summary, label=paste0("n=", n)),
            vjust=-0.5,
            hjust=-0.25,
            size=2.5,) +
  labs(x="Distance from centroid") +
  scale_fill_manual(values=c("#6FD4FF", "#F9FF00", "#79EC14", "#fb9a99"),
                    name="Habitat") +
  theme_minimal() +
  theme(axis.title.x=element_text(margin=unit(c(2, 0, 0, 0), "mm"), size=10),
        axis.title.y=element_blank(),
        axis.text=element_text(size=8),
        panel.spacing.x=unit(0.5, "in"),
        legend.position="top",
        legend.title=element_text(size=6),
        legend.text=element_text(size=5),
        plot.margin=unit(c(10, 35, 0, 10), "pt")) +
  coord_cartesian(clip="off")

#Write to file - SUPPLEMENTARY FIGURE 2
tiff(file=paste0("Supplementary_Figure_2-", Sys.Date(), ".tiff"), height=7, width=12, units="cm", res=300)
gg.disp
dev.off()


#Plot NMDS with TZ fit
gg.nmds1 <- ggplot() +
  geom_contour(data=nmds.ordi.tz.gg, 
               aes(x=x, y=y, z=z),
               colour="dimgrey",
               binwidth=5, 
               show.legend=FALSE,
               size=0.3) +
  geom_label_contour(data=nmds.ordi.tz.gg,
                     aes(x=x, y=y, z=z),
                     colour="dimgrey",
                     binwidth=5,
                     size=2,
                     label.size=NA) +
  stat_ellipse(data=metadata.nmds,
               aes(x=NMDS1, y=NMDS2, fill=Habitat.summary, colour=Habitat.summary),
               alpha=0.15,
               lty="dotted",
               geom="polygon",
               type='t',
               size=0.4,
               show.legend=FALSE) +
  geom_point(data=metadata.nmds,
             aes(x=NMDS1, y=NMDS2, colour=Habitat.summary, shape=MusaSpecies),
             size=2) +
  geom_text_repel(data=otu.gg,
                  aes(x=NMDS1, y=NMDS2, label=label),
                  fontface="italic",
                  size=2.5,
                  show.legend=FALSE) +
  geom_text(aes(x=-1.2, y=-3),
            label=paste("stress=",round(NMDS$stress, 2)),
            size=3) +
  geom_text(aes(x=-0.9, y=4.4),
            label="PERMANOVA (adonis)",
            fontface="bold",
            size=2.5) +
  geom_text(aes(x=-0.9, y=4.2),
            label=paste0("Habitat: R=", round(nmds.permanova$R2[1], 2),
                         " p=", round(nmds.permanova$`Pr(>F)`[1], 3)),
            size=2.5) +
  geom_text(aes(x=1, y=3),
            label="Seed viability",
            fontface="bold",
            size=2.5) +
  geom_text(aes(x=1.1, y=2.75),
            label=paste0("ordisurf: adj. R=", round(summary(nmds.ordi.tz)$r.sq, 2),
                         " p=", signif(summary(nmds.ordi.tz)$s.table[4], 3)),
            size=2) +
  scale_shape_manual(values=c(15:17), 
                     name=expression(paste(italic("Musa"), " species")),
                     labels=c(expression(italic("Musa acuminata")),
                              expression(italic("Musa balbisiana")),
                              expression(italic("Musa itinerans")))) +
  scale_color_manual(values=c("#79EC14","#F9FF00","#fb9a99","#6FD4FF"),
                     name="Habitat") +
  scale_fill_manual(values=c("#79EC14","#F9FF00","#fb9a99","#6FD4FF")) +
  guides(fill=FALSE,
         shape=guide_legend(order=1),
         colour=guide_legend(order=2)) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme_bw() +
  theme(legend.text.align=0,
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        aspect.ratio=1,
        panel.grid=element_blank(),
        legend.position="top",
        legend.box="vertical",
        plot.margin=margin(0, 5, 0, 5),
        legend.margin=margin(0, 0, 0, 0),
        legend.text=element_text(size=8),
        legend.title=element_text(size=9),
        plot.title.position="plot") +
  coord_fixed(ylim=c(-3.2, 3.2),
              clip="off")


#Plot NMDS with germination fit
gg.nmds2 <- ggplot() +
  geom_contour(data=nmds.ordi.germ.gg, 
               aes(x=x, y=y, z=z),
               colour="dimgrey",
               binwidth=5, 
               show.legend=FALSE,
               size=0.3) +
  geom_label_contour(data=nmds.ordi.germ.gg,
                     aes(x=x, y=y, z=z),
                     colour="dimgrey",
                     binwidth=5,
                     size=2,
                     label.size=NA) +
  stat_ellipse(data=metadata.nmds,
               aes(x=NMDS1, y=NMDS2, fill=Habitat.summary, colour=Habitat.summary),
               alpha=0.15,
               lty="dotted",
               geom="polygon",
               type='t',
               size=0.4,
               show.legend=FALSE) +
  geom_point(data=metadata.nmds,
             aes(x=NMDS1, y=NMDS2, colour=Habitat.summary, shape=MusaSpecies),
             size=2) +
  geom_text_repel(data=otu.gg,
                  aes(x=NMDS1, y=NMDS2, label=label),
                  fontface="italic",
                  size=2.5,
                  show.legend=FALSE) +
  geom_text(aes(x=-1.2, y=-3),
            label=paste("stress=",round(NMDS$stress, 2)),
            size=3) +
  geom_text(aes(x=1.2, y=3),
            label="Germination rate",
            fontface="bold",
            size=2.5) +
  geom_text(aes(x=1.2, y=2.75),
            label=paste0("ordisurf: adj. R=", round(summary(nmds.ordi.germ)$r.sq, 2),
                         " p=", round(summary(nmds.ordi.germ)$s.table[4], 3)),
            size=2) +
  scale_shape_manual(values=c(15:17), 
                     name=expression(paste(italic("Musa"), " species")),
                     labels=c(expression(italic("Musa acuminata")),
                              expression(italic("Musa balbisiana")),
                              expression(italic("Musa itinerans")))) +
  scale_color_manual(values=c("#79EC14","#F9FF00","#fb9a99","#6FD4FF"),
                     name="Habitat") +
  scale_fill_manual(values=c("#79EC14","#F9FF00","#fb9a99","#6FD4FF")) +
  guides(fill=FALSE,
         shape=guide_legend(order=1),
         colour=guide_legend(order=2)) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme_bw() +
  theme(legend.text.align=0,
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        aspect.ratio=1,
        panel.grid=element_blank(),
        legend.position="top",
        legend.box="vertical",
        plot.margin=margin(0, 5, 0, 5),
        legend.margin=margin(0, 0, 0, 0),
        legend.text=element_text(size=8),
        legend.title=element_text(size=9),
        plot.title.position="plot") +
  coord_fixed(clip="off")


## FIGURE 4B - pairwise PERMANOVA ##

#For each of the data subsets...
for (i in c("nmds", "rare")) {
  
  #Do pairwise PERMANOVA on habitat pairs
  set.seed(1)
  pw.permanova <- as.matrix(as.data.frame(pairwise.perm.manova(vegdist(get(paste0("otu.count.", i)), method="bray"),
                                                               get(paste0("metadata.", i))[,"Habitat.summary"],
                                                               p.method="BH",
                                                               nperm=999)[3]))
  #Make row and column labels consistent
  colnames(pw.permanova) <- sub("\\.", " ", sub("p.value.", "", colnames(pw.permanova)))
  #Convert to dataframe for plotting
  pw.permanova <- melt(pw.permanova, na.rm = TRUE)
  #Add data subset
  pw.permanova$data <- rep(i, length(pw.permanova$value))
  
  assign(paste0("pw.permanova.", i), pw.permanova)
}

#Combine the results dataframes
pw.permanova <- rbind(pw.permanova.nmds, pw.permanova.rare)

#Make custom annotations for coloured boxes behind habitat labels
#(https://stackoverflow.com/questions/45956883/set-text-background-to-ggplot-axis-text)

element_custom.y <- function(...) {
  structure(list(...), class = c("element_custom.y", "element_blank"))
}

element_grob.element_custom.y <- function(element, label, x, y, ...)  {
  tg <- textGrob(label, y=y, gp=gpar(col=element$colour, fontsize=5), rot=90)
  padding <- unit(0.3, "line")
  rg <- rectGrob(y=y,width=grobWidth(tg)+padding, height=unit(2,"line")+padding, 
                 gp=gpar(fill = element$fill, col=NA, alpha=1))
  gTree(children=gList(rg, tg), width=grobWidth(tg), cl="custom_axis.y")
}

widthDetails.custom_axis.y <- function(x) x$width + unit(2,"mm")

element_custom.x <- function(...) {
  structure(list(...), class = c("element_custom.x", "element_blank"))
}

element_grob.element_custom.x <- function(element, label, x, y, ...)  {
  tg <- textGrob(label, x=x, gp=gpar(col=element$colour, fontsize=5))
  padding <- unit(0.3, "line")
  rg <- rectGrob(x=x,height=grobHeight(tg)+padding, width=unit(2,"line")+padding, 
                 gp=gpar(fill = element$fill, col=NA, alpha=1))
  gTree(children=gList(rg, tg), width=grobWidth(tg), cl="custom_axis.x")
}

widthDetails.custom_axis.x <- function(x) x$width + unit(2,"mm")

#Plot grid
gg.pwperm <- ggplot(pw.permanova, aes(Var2, Var1, fill=value>0.05)) +
  facet_grid(. ~ data, labeller=labeller(data=c(`nmds`="Common OTUs", `rare`="All OTUs (including rare)"))) +
  geom_tile(color="grey", size=2, show.legend=FALSE) +
  geom_text(aes(label=value), size=3) +
  scale_fill_manual(values=c("white", "darkgrey")) +
  theme_minimal() + 
  theme(axis.text.x=element_custom.x(fill=c("#79EC14","#F9FF00", "#fb9a99")),
        axis.text.y=element_custom.y(fill=c("#F9FF00", "#fb9a99", "#6FD4FF")),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        strip.text=element_text(face="bold", size=9),
        panel.grid=element_blank(),
        plot.margin=unit(c(0, 10, 0, 10), "pt")) +
  coord_fixed()


## FIGURE 4C - DIVERSITY ##

#Subset matrix of OTU count, removing low sampled habitats
otu.count.hab.rare <- subset(otu.count, rownames(otu.count) %in% metadata$Serial.No[-c(35, 36, 39, 45)])

#Compute diversity indices
shannon.diversity <- diversity(otu.count.hab.rare, index="shannon")
simpson.diversity <- diversity(otu.count.hab.rare, index="simpson")

#Make dataframe for ggplot
div.hab.df <- data.frame(Serial.no=c(names(shannon.diversity), names(simpson.diversity)),
                         Diversity=c(shannon.diversity, simpson.diversity),
                         Index=rep(c("Shannon", "Simpson"),
                                   each=length(shannon.diversity)))
#Add habitat
div.hab.df$Habitat <- metadata$Habitat.summary[match(div.hab.df$Serial.no, metadata$Serial.No)]
#Set order in plot
div.hab.df$Habitat <- factor(div.hab.df$Habitat, levels=c("Roadside", "Jungle edge", "Jungle buffer", "Ravines"))

#For each diversity index...
for (i in c("Shannon", "Simpson")) {
  
  #TukeyHSD to test for significance between habitats
  div.hab.tukey <- TukeyHSD(aov(lm(Diversity ~ Habitat, div.hab.df[div.hab.df$Index == i,])))
  
  #Make dataframe for ggplot with tukey groups
  div.hab.labels <- data.frame(multcompLetters(div.hab.tukey[["Habitat"]][,4])["monospacedLetters"])
  
  #If only one statistical group found...
  if (length(div.hab.labels) == 0) {
    #Add tukey groups (not monospaced)
    div.hab.labels <- data.frame(multcompLetters(div.hab.tukey[["Habitat"]][,4])["Letters"])
    colnames(div.hab.labels)[1] <- "monospacedLetters"
  }
  
  #Move rownames into own column
  div.hab.labels <- cbind(Treatment=rownames(div.hab.labels), data.frame(div.hab.labels, row.names=NULL))
  #Add diversity index
  div.hab.labels$Index <- i
  
  #Add sample size and position for labelling
  div.hab.labels$n <- NA
  div.hab.labels$pos <- NA
  
  for (j in div.hab.labels$Treatment) {
    div.hab.labels[div.hab.labels$Treatment == j, "n"] <- length(div.hab.df$Habitat[div.hab.df$Habitat == j & div.hab.df$Index == i])
    div.hab.labels[div.hab.labels$Treatment == j, "pos"] <- quantile(div.hab.df$Diversity[div.hab.df$Habitat == j & div.hab.df$Index == i])[4]
  }
  
  div.hab.labels$n <- paste0("n=", div.hab.labels$n)
  
  assign(paste0("div.hab.tukey.", i), div.hab.tukey)
  assign(paste0("div.hab.labels.", i), div.hab.labels)
}

#Combine label dataframes
div.hab.labels <- rbind(div.hab.labels.Shannon, div.hab.labels.Simpson)
#Edit index names
div.hab.df$Index <- paste(div.hab.df$Index, "index")
div.hab.labels$Index <- paste(div.hab.labels$Index, "index")

#Plot boxplot
gg.div <- ggplot(div.hab.df, aes(x=Diversity, y=Habitat, fill=Habitat)) +
  geom_boxploth(show.legend=FALSE) +
  facet_wrap(. ~ Index, scales="free_x") +
  geom_text(data=div.hab.labels,
            aes(x=pos, y=Treatment, label=n),
            vjust=-0.5,
            hjust=-0.25,
            size=2,
            inherit.aes=FALSE) +
  geom_text(data=div.hab.labels,
            aes(x=Inf, y=Treatment, label=monospacedLetters),
            family="mono",
            hjust=-0.5,
            size=2.5,
            inherit.aes=FALSE) +
  scale_x_continuous(name="Diversity",
                     expand=expansion(mult=c(0, 0.1))) +
  scale_fill_manual(values=c("#6FD4FF", "#F9FF00", "#79EC14", "#fb9a99")) +
  theme_minimal() +
  theme(axis.title.x=element_text(margin=unit(c(2, 0, 0, 0), "mm"), size=8),
        axis.title.y=element_blank(),
        axis.text=element_text(size=6),
        plot.title.position="plot",
        panel.spacing.x=unit(0.5, "in"),
        plot.margin=unit(c(0, 40, 5, 10), "pt")) +
  coord_cartesian(clip="off")


## FIGURE 4D - ABUNDANCE ##

#Make dataframe of fungi per Musa individual
abund.hab.df <- data.frame(Serial.No=rownames(otu.count.hab.rare),
                           Habitat=metadata$Habitat.summary[match(rownames(otu.count.hab.rare), metadata$Serial.No)],
                           Count=rowSums(otu.count.hab.rare))
#Set order in plot
abund.hab.df$Habitat <- factor(abund.hab.df$Habitat, levels=c("Roadside", "Jungle edge", "Jungle buffer", "Ravines"))

#TukeyHSD to test for significance between habitats
hab.tukey <- TukeyHSD(aov(lm(Count ~ Habitat, data=abund.hab.df)))

#Make dataframe for ggplot with tukey groups
hab.labels <- data.frame(multcompLetters(hab.tukey[["Habitat"]][,4])["monospacedLetters"])
hab.labels <- cbind(Treatment=rownames(hab.labels), data.frame(hab.labels, row.names=NULL))

#Add sample size and position for labelling
hab.labels$n <- NA
hab.labels$pos <- NA

for (i in hab.labels$Treatment) {
  hab.labels$n[hab.labels$Treatment == i] <- length(abund.hab.df$Habitat[abund.hab.df$Habitat == i])
  hab.labels$pos[hab.labels$Treatment == i] <- quantile(abund.hab.df$Count[abund.hab.df$Habitat == i])[4]
}

hab.labels$n <- paste0("n=", hab.labels$n)

#Plot boxplot
gg.abund <- ggplot(abund.hab.df, aes(x=Count, y=Habitat, fill=Habitat)) +
  geom_boxploth(show.legend=FALSE) +
  geom_text(data=hab.labels,
            aes(x=pos, y=Treatment, label=n),
            vjust=-0.5,
            hjust=-0.25,
            size=2,
            inherit.aes=FALSE) +
  geom_text(data=hab.labels,
            aes(x=Inf, y=Treatment, label=monospacedLetters),
            family="mono",
            hjust=-0.5,
            size=2.5,
            inherit.aes=FALSE) +
  scale_x_continuous(name=expression(paste("Number of OTUs per ", italic("Musa")," accession")),
                     expand=c(0, 0)) +
  scale_fill_manual(values=c("#6FD4FF", "#F9FF00", "#79EC14", "#fb9a99")) +
  theme_minimal() +
  theme(axis.title.x=element_text(margin=unit(c(2, 0, 0, 0), "mm"), size=8),
        axis.title.y=element_blank(),
        axis.text=element_text(size=6),
        plot.title.position="plot",
        plot.margin=unit(c(5, 40, 0, 10), "pt")) +
  coord_cartesian(clip="off")


#Write to file - FIGURE 4
tiff(file=paste0("Figure_4-", Sys.Date(), ".tiff"), height=18, width=20, units="cm", res=300)
ggpubr::ggarrange(ggpubr::ggarrange(gg.nmds1, gg.nmds2, common.legend=TRUE),
                  ggpubr::ggarrange(gg.pwperm,
                                    egg::ggarrange(gg.div, gg.abund, ncol=1, labels=c("C", "D"),
                                                   label.args=list(gp=grid::gpar(font=2, cex=1.2))),
                                    ncol=2, labels=c("B", "")),
                  ncol=1, heights=c(2,1), labels=c("A", ""))
dev.off()


###################################
## FIGURE 5 - FUSARIUM PHYLOGENY ##
###################################

#Read in tree
phy <- read.tree("data/RAxML_bipartitions.Musa_fusarium_130T_3loci_introns")
#Root tree
outgroup <- c("Cylindrocarpon_cylindroides_22505", "Neonectria_coccinea_20485")
phy <- root(phy, outgroup, resolve.root=TRUE)
#Remove insignificant bootstrap values
BS <- as.numeric(phy$node.label)
BS[(BS<70)] <- ""
#Replace underscores in tip labels
phy$tip.label <- gsub("_", " ", phy$tip.label)
#Replace FIESC species names
names <- read.csv("data/tree_names.csv", header=FALSE)
phy$tip.label[match(names$V1, phy$tip.label)] <- names$V2 
#Create dataframe to signify OTUs from this study to colour
tips <- data.frame(species=phy$tip.label)
tips$colour <- ifelse(grepl("Otu", tips$species), "new", "")
#Create vector for tip colours
tip.colours <- c("black","dodgerblue3")

#Create dataframes for node collapsing, scaling and collapsing of species complexes for ggtree
clades <- data.frame(node=c(225,198,240,244,166,156,154,194,245,147,246,133,143), 
                     name=c("equiseti clade",
                            "Fusarium sambucinum SC",
                            "Fusarium tricinctum SC",
                            "Fusarium heterosporum SC",
                            "Fusarium nisikadoi SC",
                            "Fusarium oxysporum SC",
                            "Fusarium redolens SC",
                            "Fusarium concolor SC",
                            "Fusarium lateritium SC",
                            "Fusarium buharicum SC",
                            "Cyanonectria (=Fusarium buxicola SC)",
                            "Bisifusarium (=Fusarium dimerum SC)",
                            "Geejayessia (=Fusarium staphyleae SC)"), 
                     offset=c(0.029,0.15,0.05,0.038,0.035,0,0,0.15,0.11,0.17,0.085,0.4,0.15), 
                     scale=c(0.1,0.1,0.5,1,0.5,0.1,1,1,0.7,1,1,0.3,1))

highlight <- data.frame(node=c(170,210,250), 
                        name=c("Fusarium fujikuroi SC",
                               "Fusarium incarnatum-equiseti SC",
                               "Neocosmospora \n(=Fusarium solani SC)"), 
                        offset=c(0.7,0.8,1.52))

#Plot tree
gg.tree <- ggtree(phy, cex=0.1) +
  xlim(0, 3)

for (i in 1:length(clades$node)) {
  gg.tree <- scaleClade(gg.tree, node=clades$node[i], scale=clades$scale[i])
}

for (i in 1:length(highlight$node)) {
  gg.tree <- gg.tree +
    geom_cladelabel(node=highlight$node[i],
                    label=highlight$name[i],
                    offset=highlight$offset[i],
                    barsize=1,
                    fontsize=3,
                    fontface='bold.italic') +
    geom_hilight(node=highlight$node[i],
                 extend=highlight$offset[i],
                 alpha=0.3, fill="wheat2")
}

for (i in 1:length(clades$node)) {
  gg.tree <- collapse(gg.tree, node=clades$node[i], mode="max", colour="wheat2", fill="wheat2")
  gg.tree <- gg.tree +
    geom_cladelabel(node=clades$node[i],
                    label=clades$name[i],
                    barsize=0,
                    offset=clades$offset[i],
                    fontsize=3,
                    fontface='bold.italic')
}

gg.tree.1 <- gg.tree %<+% tips +
  geom_tree() + 
  geom_tiplab(aes(color=colour), size=3,  fontface='italic') +
  scale_colour_manual(values=tip.colours) +
  geom_nodelab(aes(x=branch), label=BS, size=2, vjust=-0.5) +
  geom_nodepoint(aes(subset=(node %in% c(146, 132))),
                 size=3, colour="darkred") +
  geom_nodelab(aes(subset=(node %in% c(146, 132))),
               label=c("O'Donnell et al. 2020", "Lombard et al. 2015"),
               hjust=1.1,
               vjust=-1.2,
               size=3,
               colour="darkred") +
  theme(legend.position="none")

#Write to file - FIGURE 5
tiff(file=paste0("Figure_5-", Sys.Date(), ".tiff"), height=30, width=20, units="cm", res=300)
plot(gg.tree.1)
dev.off()


#################################################################################################
## SUPPLEMENTARY Figure 2 - Abundance including oil palm plantation and botanic garden habitat ##
#################################################################################################

#Make dataframe of fungi per Musa individual
abund.hab.df.supp <- data.frame(Serial.No=rownames(otu.count),
                                Habitat=metadata$Habitat.summary[match(rownames(otu.count), metadata$Serial.No)],
                                Count=rowSums(otu.count))
#Set order in plot
abund.hab.df.supp$Habitat <- factor(abund.hab.df.supp$Habitat, levels=c("Oil palm plantation", "Botanical garden", "Roadside", "Jungle edge", "Jungle buffer", "Ravines"))


#TukeyHSD to test for significance between habitats
hab.tukey.supp <- TukeyHSD(aov(lm(Count ~ Habitat, data=abund.hab.df.supp)))

#Make dataframe for ggplot with tukey groups
hab.labels.supp <- data.frame(multcompLetters(hab.tukey.supp[["Habitat"]][,4])["monospacedLetters"])
hab.labels.supp <- cbind(Treatment=rownames(hab.labels.supp), data.frame(hab.labels.supp, row.names=NULL))

#Add sample size and position for labelling
hab.labels.supp$n <- NA
hab.labels.supp$pos <- NA

for (i in hab.labels.supp$Treatment) {
  hab.labels.supp$n[hab.labels.supp$Treatment == i] <- length(abund.hab.df.supp$Habitat[abund.hab.df.supp$Habitat == i])
  hab.labels.supp$pos[hab.labels.supp$Treatment == i] <- quantile(abund.hab.df.supp$Count[abund.hab.df.supp$Habitat == i])[4]
}

hab.labels.supp$n <- paste0("n=", hab.labels.supp$n)

#Plot boxplot
gg.abund.supp <- ggplot(abund.hab.df.supp, aes(x=Count, y=Habitat, fill=Habitat)) +
  geom_boxploth() +
  geom_text(data=hab.labels.supp,
            aes(x=pos, y=Treatment, label=n),
            vjust=-0.5,
            hjust=-0.25,
            size=2.5,
            inherit.aes=FALSE) +
  geom_text(data=hab.labels.supp,
            aes(x=Inf, y=Treatment, label=monospacedLetters),
            family="mono",
            hjust=-0.5,
            size=3,
            inherit.aes=FALSE) +
  scale_x_continuous(name=expression(paste("Number of OTUs per ", italic("Musa")," accession")),
                     expand=c(0, 0)) +
  scale_fill_manual(values=c("White", "grey", "#6FD4FF", "#F9FF00", "#79EC14", "#fb9a99")) +
  theme_minimal() +
  theme(axis.title.x=element_text(margin=unit(c(2, 0, 0, 0), "mm"), size=10),
        axis.title.y=element_blank(),
        axis.text=element_text(size=8),
        legend.position="top",
        legend.title=element_text(size=6),
        legend.text=element_text(size=5),
        plot.margin=unit(c(10, 35, 0, 10), "pt")) +
  coord_cartesian(clip = 'off')
plot(gg.abund.supp)

#Write to file - SUPPLEMENTARY FIGURE 3
tiff(file=paste0("Supplementary_Figure_3-", Sys.Date(), ".tiff"), height=7.5, width=10, units="cm", res=300)
gg.abund.supp
dev.off()


######################################################
## SUPPLEMENTARY TABLE 4 - OTUs FOUND OUTSIDE SEEDS ##
######################################################

#Add OTU to metadata dataframe
df.controls["OTU"] <- otus$OTU[match(df.controls$Serial.No, otus$Serial.No)]
#Add classification to metadata dataframe
df.controls[,c("Phylum","Class","Order","Genus","Species")] <- classification.df[,c("phy.","class.","ord.","gen.","sp.")][match(df.controls$OTU, classification.df$OTU),]
#Remove accessions which couldn't cluster into OTUs
df.controls <- subset(df.controls, OTU != "#N/A")
#Filter for imprints
df.controls <- subset(df.controls, df.controls$Direct.Sequencing.or.Culture ==  "Control")
#Remove unnecessary columns
df.controls <- df.controls[c(4, 9)]
#Add number of times found in imprints
df.controls$Count <- NA

for (i in 1:length(df.controls$OTU)) {
  df.controls$Count[i] <- length(grep(df.controls$OTU[i], df.controls$OTU))
}

#Remove duplicate rows
df.controls <- unique(df.controls)
#Write to file
write.csv(df.controls, "Supplementary_Table_4.csv", row.names=FALSE)


####################
## STATS FOR TEXT ##
####################

## Proportion of endophytes reported as other guilds in FUNguild ##

#Get function to download FUNGuild database from https://rdrr.io/github/vmikk/metagMisc/man/parse_funguild.html
parse_funguild <- function(url='http://www.stbates.org/funguild_db.php', tax_name=TRUE){
  
  # require(XML)
  # require(jsonlite)
  # require(RCurl)
  
  ## Parse data
  tmp <- XML::htmlParse(url)
  tmp <- XML::xpathSApply(doc=tmp, path="//body", fun=XML::xmlValue)
  
  ## Read url and convert to data.frame
  db <- jsonlite::fromJSON(txt=tmp)
  
  ## Remove IDs
  db$`_id` <- NULL
  
  if(tax_name == TRUE){
    
    ## Code legend
    ## Taxon Level: A numeral corresponding the correct taxonomic level for the taxon
    taxons <- c(
      "keyword",                                                       # 0
      "Phylum", "Subphylum", "Class", "Subclass", "Order", "Suborder", # 3:8
      "Family", "Subfamily", "Tribe", "Subtribe", "Genus",             # 9:13
      "Subgenus", "Section", "Subsection", "Series", "Subseries",      # 15:19
      "Species", "Subspecies", "Variety", "Subvariety", "Form",        # 20:24
      "Subform", "Form Species")
    
    ## Table with coding
    taxmatch <- data.frame(
      TaxID=c(0, 3:13, 15:26),
      Taxon=factor(taxons, levels=taxons))
    
    ## Match taxon codes
    db$taxonomicLevel <- taxmatch[match(x=db$taxonomicLevel, table=taxmatch$TaxID), "Taxon"]
  }
  
  # remove rows with missing data
  # which(
  # 	with(db, trophicMode == "NULL" & guild == "NULL" & growthForm == "NULL" & trait == "NULL" & notes == "NULL")
  # 	)
  
  ## Add database dump date as attributes to the result
  attr(db, "DownloadDate") <- date()
  
  return(db)
}

#Download FUNguild database
funguild <- parse_funguild(url="http://www.stbates.org/funguild_db.php", tax_name=TRUE)
#Filter for endophytes
funguild <- funguild[grep("Endophyt*", funguild$guild),]
#Filter for species level
funguild <- funguild[funguild$taxonomicLevel == "Species",]
#Make list of guilds for each entry
funguild.guilds <- str_split(as.list(funguild$guild), "-")
#For each guild...
for (i in 1:length(unique(unlist(funguild.guilds)))) {
  #Print number reported as both endophytes and another guild
  print(unique(unlist(funguild.guilds))[i])
  print(paste0(length(grep(unique(unlist(funguild.guilds))[i], funguild.guilds)), ", ", round(length(grep(unique(unlist(funguild.guilds))[i], funguild.guilds)) / length(funguild.guilds) * 100), "%"))
}


#Total number of seeds
print(sum(seed.df$y))
#Number of seeds containing endophytes
print(paste0(sum(seed.count[2:6]), ", ", round(sum(seed.count[2:6]) / sum(seed.df$y) * 100), "%"))

#Number of endophytes
paste0("Total: ", length(df$Serial.No),
       " (Culture: ", length(which(df$Direct.Sequencing.or.Culture == "Culture")),
       ", Direct: ", length(which(df$Direct.Sequencing.or.Culture == "Direct" & df$Clone == "")),
       ", Clones: ", length(which(df$Clone == "Y")), ")")
#Total number of OTUs
print(length(unique(df$OTU)))
#Number of singleton OTUs
singletons <- data.frame(otu=unique(df$OTU), count=NA)
for (i in 1:length(singletons$otu)) {
  singletons$count[i] <- length(grep(singletons$otu[i], df$OTU))
}
print(paste0(length(singletons$otu[singletons$count == 1]), ", ", round(length(singletons$otu[singletons$count == 1]) / length(singletons$otu) * 100), "%"))

#Number of Ascomycota, Basidio and unknown
for (i in c("Ascomycota", "Basidiomycota", "Fungi phy.")) {
  print(paste0(i, ": ", length(unique(df$OTU[df$Phylum == i])),
               " (", round(length(unique(df$OTU[df$Phylum == i])) / length(unique(df$OTU)) * 100), "%)"))
}


#Percentage of endophytes that were Lasiodiplodia, Aspergillus and Fusarium
for (i in c("Lasiodiplodia", "Fusarium", "Aspergillus")) {
  print(paste0(i, ", ", length(grep(i, df$Genus)), ", ", round(length(grep(i, df$Genus)) / length(df$Genus) * 100), "%"))
}
#Proportion of all endophytes belonging to those three genera
(length(grep("Lasiodiplodia", df$Genus)) + length(grep("Fusarium", df$Genus)) + length(grep("Aspergillus", df$Genus))) / length(df$Serial.No) * 100

#Percentage of significantly supported branches in tree
length(which(as.numeric(BS) >= 70)) / length(BS) * 100


########################
## GRAPHICAL ABSTRACT ##
########################

#Read in image of seeds and cultures
img <- readPNG("Seeds_cultures.png")

#Plot image as ggplot
gg.img <- ggplot() + 
  background_image(img) +
  coord_fixed(xlim=c(0, dim(img)[2]), 
              ylim=c(0, dim(img)[1]), 
              expand = FALSE) +
  ggtitle("Fungal endophytes from\nstored wild banana seeds") +
  theme_void() +
  theme(plot.title=element_text(size=15, face="bold", margin=margin(0, 0, 20, 0), hjust=0.5))

#Adapt NMDS plot
gg.nmds.abs <- ggplot() +
  geom_contour(data=nmds.ordi.germ.gg, 
               aes(x=x, y=y, z=z),
               colour="dimgrey",
               binwidth=5, 
               show.legend=FALSE,
               size=0.3) +
  geom_label_contour(data=nmds.ordi.germ.gg,
                     aes(x=x, y=y, z=z),
                     colour="dimgrey",
                     binwidth=5,
                     size=2,
                     label.size=NA) +
  stat_ellipse(data=metadata.nmds,
               aes(x=NMDS1, y=NMDS2, fill=Habitat.summary, colour=Habitat.summary),
               alpha=0.15,
               lty="dotted",
               geom="polygon",
               type='t',
               size=0.4,
               show.legend=FALSE) +
  geom_point(data=metadata.nmds,
             aes(x=NMDS1, y=NMDS2, colour=Habitat.summary, shape=MusaSpecies),
             size=2) +
  geom_text_repel(data=otu.gg,
                  aes(x=NMDS1, y=NMDS2, label=label),
                  fontface="italic",
                  size=3,
                  show.legend=FALSE) +
  geom_text(aes(x=-1.2, y=-3),
            label=paste("stress=",round(NMDS$stress, 2)),
            size=3) +
  geom_text(aes(x=1, y=3),
            label="Germination rate",
            fontface="bold",
            size=4) +
  geom_text(aes(x=1, y=2.75),
            label=paste0("ordisurf: adj. R=", round(summary(nmds.ordi.germ)$r.sq, 2),
                         " p=", round(summary(nmds.ordi.germ)$s.table[4], 3)),
            size=3) +
  geom_curve(aes(x=1, y=2.5, xend=0.3, yend=0.6),
             size=1,
             angle=90,
             curvature=0.1,
             arrow=arrow(length=unit(0.03, "npc"))) +
  scale_shape_manual(values=c(15:17)) +
  scale_color_manual(values=c("#79EC14","#F9FF00","#fb9a99","#6FD4FF"),
                     name="Habitat") +
  scale_fill_manual(values=c("#79EC14","#F9FF00","#fb9a99","#6FD4FF")) +
  guides(fill=FALSE,
         shape=FALSE,
         colour=guide_legend(override.aes=list(size=5))) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme_bw() +
  theme(legend.text.align=0,
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        aspect.ratio=1,
        panel.grid=element_blank(),
        legend.position="top",
        legend.box="vertical",
        plot.margin=margin(0, 5, 0, 5),
        legend.margin=margin(0, 0, 0, 0),
        legend.text=element_text(size=11),
        legend.title=element_text(size=14, face="bold"),
        plot.title.position="plot") +
  coord_fixed(clip="off")

#Adapt abundance plot
gg.abund.abs <- gg.abund +
  scale_x_continuous(name="Abundance",
                     expand=c(0, 0)) +
  theme(axis.title.x=element_text(margin=unit(c(2, 0, 0, 0), "mm"), size=15, face="bold"),
        axis.text=element_text(size=8),
        axis.text.y=element_blank(),
        plot.margin=unit(c(10, 40, 5, 10), "pt")) +
  coord_cartesian(clip="off")

#Adapt diversity plot
gg.div.abs <- gg.div +
  scale_x_continuous(name="Diversity",
                     expand=expansion(mult=c(0, 0.1))) +
  scale_fill_manual(values=c("#6FD4FF", "#F9FF00", "#79EC14", "#fb9a99")) +
  theme_minimal() +
  theme(axis.title.x=element_text(margin=unit(c(2, 0, 0, 0), "mm"), size=15, face="bold"),
        axis.title.y=element_blank(),
        axis.text=element_text(size=8),
        axis.text.y=element_blank(),
        panel.spacing.x=unit(0.5, "in"),
        plot.margin=unit(c(5, 40, 10, 10), "pt")) +
  coord_cartesian(clip="off")

#Write to file
png(file=paste0("Graphical-abstract-", Sys.Date(), ".png"), height=5, width=12.5, units="in", res=300)
ggarrange(gg.img,
          gg.nmds.abs,
          ggarrange(gg.div.abs, gg.abund.abs, ncol=1),
          ncol=3,
          widths=c(1, 1.5, 1))
dev.off()
