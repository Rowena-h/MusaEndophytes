#####################################################
#####################################################
#####################################################
######                                         ######
######    Fungal Endophytes in Musa CWR Seeds  ######
######                                         ######
#####################################################
#####################################################
#####################################################

setwd("~/Kew/Musa Endophytes/R analysis/")

library(stringr)
library(ggplot2)
library(scales)
library(vegan)
library(ape)
library(ggtree)
library(ggpubr)
library(ggrepel)
library(ggstance)
library(ggnewscale)
library(egg)
library(multcompView)
library(eulerr)
library(ggplotify)
library(colorspace)
library(ggwordcloud)
library(ggtree)


#############################################################
## FIGURE 1 - OTUs PER SEED AND SPECIES ACCUMULATION CURVE ##
#############################################################

## FIGURE 1A - UNIQUE OTUs PER SEED ##

#Read in all data
df.all <- read.csv("data/All_data.csv")
#Remove controls
df.all <- subset(df.all, df.all$Sterilised == "Y")
df.all <- subset(df.all, df.all$Whole.Seed.Half.Seed != "Seed Removed Control")
df.all <- df.all[-grep("CONTROL", df.all$Culture.Notes),]

#Read in ITS OTU mapping from USEARCH
otus <- read.table("data/030919_5.8SITS2LSU_trimmed_otus99map_relabel.txt", sep = '\t', col.names=c("Serial.No","OTU"))
#Read in OTU classification from T-BAS
tbas.class <- read.csv("data/TBAS_classification.csv")
#Read in metadata
metadata <- read.csv("data/metadata.csv")

#Make new OTU metadata dataframe
df <- df.all
#Add OTU to metadata dataframe
df["OTU"] <- otus$OTU[match(df$Serial.No, otus$Serial.No)]
#Add classification to metadata dataframe
df[,c("Phylum","Class","Order","Genus")] <- tbas.class[,c("Phylum","Class","Order","Genus")][match(df$OTU, tbas.class$OTU),]
#Remove accessions which couldn't cluster into OTUs
df <- subset(df, OTU != "#N/A")
#Remove empty rows
df <- df[!apply(is.na(df) | df == "", 1, all),]
#Remove duplicate clones
df <- subset(df, !duplicated(data.frame(Serial.No=sub("\\--.*", "", df$Serial.No), OTU=df$OTU)))
#Remove controls
df <- df[df$Sterilised == "Y",]

for (method in c("Culture", "Direct")) {
  
  #Extract vector of seeds
  seed <- substring(df$Serial.No[df$Direct.Sequencing.or.Culture == method],1,9)[!duplicated(substring(df$Serial.No[df$Direct.Sequencing.or.Culture == method],1,9))]
  
  seed.df <- data.frame(Serial.No=seed, count=NA)
  
  #Count number of fungi per seed (all Musa individuals)
  for (i in 1:length(seed)) {
    seed.df$count[i] <- length(grep(seed[i], df$Serial.No))
  }
  
  seed.count <- data.frame(x=sort(unique(seed.df$count)), y=NA)
  
  for (i in 1:length(seed.count$x)) {
    seed.count$y[i] <- length(seed.df$count[seed.df$count == seed.count$x[i]])
  }
  
  #Count number of fungi per seed (by Musa species)
  
  #Add Musa species column to seed dataframe
  seed.df["musa"] <- metadata$Species[match(gsub("-.*", "", seed.df$Serial.No), metadata$Serial.No)]
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
  df.tot[df.tot$method == "Culture" & df.tot$musa == i, "y"] <- length(unique(substring(df.all$Serial.No[grep(i, df.all$Species)][df.all$Direct.Sequencing.or.Culture[grep(i, df.all$Species)] == "Culture"], 1, 9)))
  df.tot[df.tot$method == "Direct" & df.tot$musa == i, "y"] <- length(unique(substring(df.all$Serial.No[grep(i, df.all$Species)][df.all$Direct.Sequencing.or.Culture[grep(i, df.all$Species)] == "Direct"], 1, 9)))
}

#Remove numbers of seeds with isolates from the total seed count
df.tot$y[df.tot$method == "Culture"] <- df.tot$y[df.tot$method == "Culture"] - colSums(seed.count.Culture[2:6])
df.tot$y[df.tot$method == "Direct"] <- df.tot$y[df.tot$method == "Direct"] - colSums(seed.count.Direct[2:6])

#Combine counts for all individuals and for separate Musa species into one dataframe for plotting
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
                     limits=c(0, 250),
                     breaks=seq(0, 250, 50),
                     expand=c(0, 5)) +
  theme_minimal() +
  theme(axis.title.x=element_text(margin=unit(c(3, 0, 0, 0), "mm")),
        axis.title.y=element_text(margin=unit(c(0, 3, 0, 0), "mm")),
        panel.spacing.y=unit(1, "lines"),
        strip.text.x=element_text(face="italic", size=7),
        strip.text.y=element_text(size=9),
        plot.title.position="plot") +
  labs(subtitle=expression(bold("A")))
plot(gg.seed)


## FIGURE 1B - SPECIES ACCUMULATION CURVE ##

#Make dataframe of species
spec.df <- metadata[c("Serial.No", "Species")]
#Group subspecies and varieties
spec.df$Species <- word(spec.df$Species, 1, 2)

#Create OTU abundance matrix
tree.count <- matrix(nrow=length(unique(spec.df$Serial.No)),ncol=length(unique(df$OTU)))
colnames(tree.count) <- unique(df$OTU)
rownames(tree.count) <- unique(spec.df$Serial.No)

for (i in 1:length(rownames(tree.count))) {
  for (j in 1:length(colnames(tree.count))) {
    tree.count[i,j] <- length(df$OTU[grep(rownames(tree.count)[i], df$Serial.No)][df$OTU[grep(rownames(tree.count)[i], df$Serial.No)] == colnames(tree.count)[j]])
  }
}

#Accumulation curves including singletons

#Create dataframe to collect rarefaction accumulation curve data for acuminata, balbisians and itinerans
spec.accum.df <- data.frame(individual=NA, richness=NA, error=NA, musa=sort(spec.df$Species[spec.df$Species == musa[1] | spec.df$Species == musa[2] | spec.df$Species == musa[3]]), rare="Including singletons")

#Rarefaction for acuminata, balbisians and itinerans
for (i in 1:length(musa[1:3])) {
  #Subset OTU abundance matrix for top three host species
  tmp.df <- subset(tree.count, rownames(tree.count) %in% rownames(tree.count)[spec.df$Species[match(rownames(tree.count), spec.df$Serial.No)] == musa[i]])
  #Species accumulation
  tmp.accum <- specaccum(tmp.df, method="rarefaction")
  #Add to results dataframe
  spec.accum.df[spec.accum.df$musa == musa[i],] <- data.frame(individual=tmp.accum$sites, richness=tmp.accum$richness, error=tmp.accum$sd, musa=musa[i], rare="Including singletons")
}

#Rarefaction for all Musa individuals
tmp.accum <- specaccum(tree.count, method="rarefaction")
tmp.accum.df <- data.frame(individual=tmp.accum$sites, richness=tmp.accum$richness, error=tmp.accum$sd, musa="All", rare="Including singletons")

#Combine species accumulation data including singletons
spec.accum.df <- rbind(spec.accum.df, tmp.accum.df)


#Accumulation curves excluding singletons

#Create dataframe to collect rarefaction accumulation curve data for acuminata, balbisians and itinerans
spec.accum.excl.df <- data.frame(individual=NA, richness=NA, error=NA, musa=sort(spec.df$Species[spec.df$Species == musa[1] | spec.df$Species == musa[2] | spec.df$Species == musa[3]]), rare="Excluding singletons")

#Remove singleton OTUs
tree.count.excl <- subset(tree.count, select=-which(colSums(tree.count) == 1))

#Rarefaction for acuminata, balbisians and itinerans
for (i in 1:length(musa[1:3])) {
  #Subset for top three host species
  tmp.df <- subset(tree.count.excl, rownames(tree.count.excl) %in% rownames(tree.count.excl)[spec.df$Species[match(rownames(tree.count.excl), spec.df$Serial.No)] == musa[i]])
  #Species accumulation
  tmp.accum <- specaccum(tmp.df, method="rarefaction")
  #Add to results dataframe
  spec.accum.excl.df[spec.accum.excl.df$musa == musa[i],] <- data.frame(individual=tmp.accum$sites, richness=tmp.accum$richness, error=tmp.accum$sd, musa=musa[i], rare="Excluding singletons")
}

#Rarefaction for all Musa individuals
tmp.accum <- specaccum(tree.count.excl, method="rarefaction")
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
  labs(x=expression(paste("Number of ",italic("Musa")," individuals sampled")), y="Number of unique OTUs") +
  theme_minimal() +
  theme(axis.title.x=element_text(margin=unit(c(3, 0, 0, 0), "mm")),
        axis.title.y=element_text(margin=unit(c(0, 3, 0, 0), "mm")),
        legend.text.align=0,
        plot.title.position="plot") +
  labs(subtitle=expression(bold("B")))
plot(gg.accum)

#Write to file
tiff(file=paste0("Figure_1-", Sys.Date(), ".tiff"), height=7, width=6, units="in", res=600)
ggpubr::ggarrange(gg.seed, gg.accum, nrow=2, heights=c(1.5,1))
dev.off()


###################################
## FIGURE 2 - OTU CLASSIFICATION ##
###################################

## T-BAS VISUALISATION OF OTU CLASSIFICATION ##

#Read in TBAS tree
tbas <- read.tree("data/TBAS_030919.tre")
#Make dataframe of isolation method
tbas.df <- data.frame(tip=tbas$tip.label, direct=NA, culture=NA, clone=NA)

#Make dataframe of OTU count per isolation method
method.count <- data.frame(otu=unique(df$OTU), direct=NA, culture=NA, clone=NA)

for (i in 1:length(unique(df$OTU))) {
  method.count$direct[i] <- sum(df$Direct.Sequencing.or.Culture[df$OTU == unique(df$OTU)[i]] == "Direct" & df$Clone[df$OTU == unique(df$OTU)[i]] != "Y")
  method.count$culture[i] <- sum(df$Direct.Sequencing.or.Culture[df$OTU == unique(df$OTU)[i]] == "Culture")
  method.count$clone[i] <- sum(df$Clone[df$OTU == unique(df$OTU)[i]] == "Y")
  if (sum(df$Sterilised[df$OTU == unique(df$OTU)[i]] == "N") > 0) {
    method.count$control[i] <- 1
  }
}

#Add counts to TBAS plot dataframe
tbas.df$direct <- method.count$direct[match(tbas.df$tip, method.count$otu)]
tbas.df$culture <- method.count$culture[match(tbas.df$tip, method.count$otu)]
tbas.df$clone <- method.count$clone[match(tbas.df$tip, method.count$otu)]
#Replace 0 with NA
tbas.df[tbas.df == 0] <- NA

cladelabels <- read.csv("data/TBAS_clades.csv")

#Plot tree
gg.tbas <- ggtree(tbas,
                  layout="circular",
                  cex=0.2,
                  branch.length="none")

for (i in 1:nrow(cladelabels)) {
  gg.tbas <- gg.tbas +
    geom_hilight(node=cladelabels[i,2], fill=cladelabels[i,3])
}

gg.tbas <- gg.tbas %<+% tbas.df +
  geom_tree(cex=0.2) +
  geom_tippoint(aes(x=x+1, size=culture, colour="Culturing"), na.rm=TRUE) +    
  geom_tippoint(aes(x=x+3, size=direct, colour="Direct extraction"), na.rm=TRUE) +
  geom_tippoint(aes(x=x+3, size=clone, colour="Cloning"), alpha=0.4, na.rm=TRUE) +
  scale_colour_manual(name="Isolation method", labels=c("Cloning", "Culturing","Direct extraction"), values=c("black", "#E69F00", "#56B4E9")) +
  labs(size="Number of isolates per OTU") +
  guides(colour=guide_legend(override.aes=list(size=3))) +
  theme(legend.position = "right",
        plot.margin=unit(c(-10,0,-10,0), "mm"),
        legend.box.margin=margin(c(0,0,0,0)),
        plot.title.position="plot")
plot(gg.tbas)


## PIECHART VISUALISATION OF OTU CLASSIFICATION ##

#Make taxonomy dataframe
pie.df <- data.frame(phylum=df$Phylum, class=df$Class, order=df$Order, genus=df$Genus)
#Remove duplicate rows for genera
pie.df <- pie.df[!duplicated(pie.df$genus),]

#Add column of total count for each genus
pie.df["num"] <- NA

for (i in 1:length(pie.df$genus)) {
  pie.df$num[i] <- sum(df$Genus == pie.df$genus[i])
}

#Sort dataframe sequentially by class, order and genus
pie.df <- transform(pie.df, freq=ave(seq(nrow(pie.df)), class, FUN=length))
pie.df <- transform(pie.df, freq2=ave(seq(nrow(pie.df)), order, FUN=length))
pie.df <- pie.df[order(-pie.df["freq"], -pie.df["freq2"], pie.df["order"], -pie.df["num"], pie.df["genus"]), ]

#Add ymin and ymax columns for size of pie slice
ymin <- c(0,pie.df$num)

for (i in 2:length(ymin)) {
  ymin[i] <- ymin[i] + ymin[i-1]
}

pie.df["ymin"] <- ymin[1:length(pie.df$genus)]
pie.df["ymax"] <- ymin[2:length(ymin)]

#Add column with midpoint position for genus label in pie slice
pie.df["lab.pos.gen"] <- (pie.df$ymax + pie.df$ymin) / 2

#Make dataframe of order label positions
pos.order <- data.frame(order=unique(pie.df$order), lab.pos.order=NA)

for (i in 1:length(unique(pie.df$order))){
  pos.order$lab.pos.order[pos.order$order == unique(pie.df$order)[i]] <- pie.df$ymax[pie.df$order == unique(pie.df$order)[i]][which.max(pie.df$ymax[pie.df$order == unique(pie.df$order)[i]])]
}

ymin <- c(0,pos.order$lab.pos.order)
pos.order["ymin"] <- ymin[1:length(pos.order$order)]
pos.order["ymax"] <- ymin[2:length(ymin)]

res <- vector(length=length(pos.order$lab.pos.order))

for (i in 2:length(ymin)) {
  res[i-1] <- (ymin[i] + ymin[i-1]) / 2
}

pos.order$lab.pos.order <- res

#Make dataframe of class label positions
pos.class <- data.frame(class=unique(pie.df$class), lab.pos.class=NA)

for (i in 1:length(unique(pie.df$class))){
  pos.class$lab.pos.class[pos.class$class == unique(pie.df$class)[i]] <- pie.df$ymax[pie.df$class == unique(pie.df$class)[i]][which.max(pie.df$ymax[pie.df$class == unique(pie.df$class)[i]])]
}

ymin <- c(0,pos.class$lab.pos.class)
pos.class["ymin"] <- ymin[1:length(pos.class$class)]
pos.class["ymax"] <- ymin[2:length(ymin)]

res <- vector(length=length(pos.class$lab.pos.class))

for (i in 2:length(ymin)) {
  res[i-1] <- (ymin[i] + ymin[i-1]) / 2
}

pos.class$lab.pos.class <- res

#Add unclassified Ascomycota to the T-BAS clade colours
colours <- rbind(cladelabels[-2], data.frame(class="Ascomycota class", colour="#818181", mid="#9b9b9b", low="#b5b5b5", lowest="#b5b5b5"))
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
  
  #Genus
  pie.df$colour.genus[pie.df$class == unique(pie.df$class)[i]] <- colorRampPalette(c(as.character(colours$low[i]),as.character(colours$lowest[i])))(length(unique(pie.df$genus[pie.df$class == unique(pie.df$class)[i]])))
}

#Create a vector with names from all taxon levels and sort alphabetically
all.tax <- c(as.character(pie.df$class), as.character(pie.df$order), as.character(pie.df$genus))
all.tax <- all.tax[order(all.tax)]
all.tax <- all.tax[!duplicated(all.tax)]

#Create a dataframe of colours for ggplot
pie.colours <- rbind(data.frame(name=pie.df$class[!duplicated(pie.df$class)], colour=pie.df$colour.class[!duplicated(pie.df$class)]), data.frame(name=pie.df$order[!duplicated(pie.df$order)], colour=pie.df$colour.order[!duplicated(pie.df$order)]), data.frame(name=pie.df$genus[!duplicated(pie.df$genus)], colour=pie.df$colour.genus[!duplicated(pie.df$genus)]))
#Sorted alphabetically and by taxonomic level
pie.colours <- pie.colours[order(pie.colours$name),]
#Make into vector
pie.colours <- as.vector(pie.colours$colour[match(all.tax, pie.colours$name)])

#Plot piechart
gg.pie <- ggplot(pie.df) +
  geom_rect(aes(fill=genus, ymax=ymax, ymin=ymin, xmax=6, xmin=4), 
            colour="white",
            show.legend = FALSE) +
  geom_rect(data=pos.order, 
            aes(fill=order, ymax=ymax, ymin=ymin, xmax=4, xmin=2), 
            colour="white",
            show.legend = FALSE) +
  geom_rect(data=pos.class,
            aes(fill=class, ymax=ymax, ymin=ymin, xmax=2, xmin=0), 
            colour="white",
            show.legend = FALSE) +
  annotate("text", x=1, y=439, vjust=-0.5, label="Class", color="white", size=4, fontface="bold", angle=45) +
  annotate("text", x=3, y=439, vjust=1, label="Order", color="white", size=4, fontface="bold", angle=45) +
  annotate("text", x=5, y=439, vjust=1.5, label="Genus", color="white", size=4, fontface="bold", angle=45) +
  geom_label_repel(aes(x=6, y=lab.pos.gen, label=genus, fill=genus),
                   direction="y",
                   fontface="italic", 
                   size=2.5, 
                   nudge_x=1,
                   show.legend = FALSE) +
  geom_label_repel(data=pos.order[1:9,], 
                   aes(x=3, y=lab.pos.order, label=order, fill=order),
                   direction="y",
                   size=2.5,
                   nudge_x=1,
                   show.legend = FALSE) +
  geom_point(data=pos.class,
             aes(x=1, y=lab.pos.class, color=class),
             size=2,
             shape=15,
             alpha=0) +
  scale_fill_manual(values=pie.colours) +
  scale_color_manual(name="", labels=c("Unclassified Ascomycota",as.character(sort(pos.class$class)[2:9])), values=as.vector(colours$colour[order(as.vector(colours$class))])) +
  xlim(c(0, 8)) +
  coord_polar(theta="y") +
  theme_void() +
  guides(color=guide_legend(override.aes=list(alpha=1, size=5))) +
  theme(legend.position="left",
        plot.margin=unit(c(-10,0,-10,0), "mm"),
        legend.box.margin=margin(c(-300,-80,0,0)))
plot(gg.pie)

#Write to file
tiff(file=paste0("Figure_2-", Sys.Date(), ".tiff"), height=11, width=9, units="in", res=600)
ggpubr::ggarrange(gg.tbas, gg.pie, heights=c(1,1.5), ncol=1)
dev.off()


################################################
## FIGURE 3 - COMPARISON OF DETECTION METHODS ##
################################################

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
euler.df <- rbind(data.frame(otu=df$Genus[match(method.count$otu[direct], df$OTU)],
                          count=rowSums(method.count[direct, 2:4]),
                          method="direct",
                          x=0.3, y=0.25),
               data.frame(otu=df$Genus[match(method.count$otu[culture], df$OTU)],
                          count=rowSums(method.count[culture, 2:4]),
                          method="culture",
                          x=0.7, y=0.25),
               data.frame(otu=df$Genus[match(method.count$otu[clone], df$OTU)],
                          count=rowSums(method.count[clone, 2:4]),
                          method="clone",
                          x=0.52, y=0.8),
               data.frame(otu=df$Genus[match(method.count$otu[direct.culture], df$OTU)],
                          count=rowSums(method.count[direct.culture, 2:4]),
                          method="direct.culture",
                          x=0.5, y=0.34),
               data.frame(otu=df$Genus[match(method.count$otu[direct.clone], df$OTU)],
                          count=rowSums(method.count[direct.clone, 2:4]),
                          method="direct.clone",
                          x=0.4, y=0.59),
               data.frame(otu=df$Genus[match(method.count$otu[culture.clone], df$OTU)],
                          count=rowSums(method.count[culture.clone, 2:4]),
                          method="culture.clone",
                          x=0.61, y=0.52),
               data.frame(otu=df$Genus[match(method.count$otu[direct.culture.clone], df$OTU)],
                          count=rowSums(method.count[direct.culture.clone, 2:4]),
                          method="direct.culture.clone",
                          x=0.5, y=0.49))

#Plot euler
gg.euler <- gg.euler +
  geom_text(x=0.1, y=0.5, size=5, fontface="bold", colour="#56B4E9", 
            label=paste0("Direct Extraction\n", length(method.count$otu[direct]))) +
  geom_text(x=0.85, y=0.5, colour="#E69F00", size=5, fontface="bold",
            label=paste0("Culture\n", length(method.count$otu[culture]))) +
  geom_text(x=0.3, y=0.9, colour="black", size=5, fontface="bold",
            label=paste0("Clone\n", length(method.count$otu[clone]))) +
  geom_text(x=0.39, y=0.49, colour="white",
            size=5, fontface="bold",
            label=length(method.count$otu[direct.clone])) +
  geom_text(x=0.515, y=0.285, colour="white",
            size=5, fontface="bold",
            label=length(method.count$otu[direct.culture])) +
  geom_text(x=0.58, y=0.58, colour="white",
            size=5, fontface="bold",
            label=length(method.count$otu[culture.clone])) +
  geom_text(x=0.5, y=0.59, colour="white",
            size=5, fontface="bold",
            label=length(method.count$otu[direct.culture.clone])) +
  geom_text_wordcloud(data=euler.df,
                      eccentricity=1,
                      area_corr_power=1,
                      aes(label=otu, colour=method, size=count)) +
  scale_size_continuous(range=c(1.5, 4)) +
  scale_colour_manual(values=c("#000000", #clone
                               "#E69F00", #culture
                               hex(mixcolor(alpha=0.5, hex2RGB("#E69F00"), hex2RGB("#000000"))), #culture.clone
                               "#56B4E9", #direct
                               hex(mixcolor(alpha=0.5, hex2RGB("#56B4E9"), hex2RGB("#000000"))), #direct.clone
                               hex(mixcolor(alpha=0.5, hex2RGB("#56B4E9"), hex2RGB("#E69F00"))), #direct.culture
                               hex(mixcolor(alpha=0.5, (mixcolor(alpha=0.5, hex2RGB("#E69F00"), hex2RGB("#000000"))), hex2RGB("#56B4E9"))))) + #direct.culture.clone
  theme(legend.position="none")
plot(gg.euler)

#Write to file
tiff(file=paste0("Figure_3-", Sys.Date(), ".tiff"), height=6, width=9, units="in", res=300)
gg.euler
dev.off()


#####################################################
## FIGURE 4 - COMMUNITY COMPOSITION AND ABUNDANCE  ##
#####################################################

## FIGURE 4A - NMDS ANALYSIS ##

#Subset OTU abundance matrix by Musa individuals for acuminata, balbisiana and itinerans
species.count.nmds <- subset(tree.count, rownames(tree.count) %in% rownames(tree.count)[spec.df$Species[match(rownames(tree.count), spec.df$Serial.No)] == musa[1] | spec.df$Species[match(rownames(tree.count), spec.df$Serial.No)] == musa[2] | spec.df$Species[match(rownames(tree.count), spec.df$Serial.No)] == musa[3]])
#Filter for most abundant OTUs
species.count.nmds <- species.count.nmds[,-which(colSums(species.count.nmds) < 20)]
#Remove rows with no OTUs
species.count.nmds <- species.count.nmds[-which(rowSums(species.count.nmds) == 0),]

#Perform NMDS for 1-10 dimensions to pick optimal number of dimensions
NMDS.scree <- function(x) {
  plot(rep(1, 10), replicate(10, metaMDS(x, autotransform = F, k = 1)$stress), xlim = c(1, 10),ylim = c(0, 0.30), xlab = "# of Dimensions", ylab = "Stress")
  for (i in 1:10) {
    points(rep(i + 1,10),replicate(10, metaMDS(x, autotransform = F, k = i + 1)$stress))
  }
}

## Supplementary Figure 1 ##
tiff(file=paste0("Supplementary_Figure_1-", Sys.Date(), ".tiff"), height=8, width=7, units="in", res=600)
#Plot screeplot
par(mfrow=c(2,1), mar=c(4,4,2,2))
NMDS.scree(species.count.nmds)

#NMDS with optimal dimensions
NMDS <- metaMDS(species.count.nmds, k=6, trymax=1000, trace=F)

#Plot stressplot
stressplot(NMDS)
dev.off()

#Make dataframes of NMDS site scores for ggplot (site=Musa individual)
nmds.gg <- as.data.frame(scores(NMDS))
nmds.gg$Serial.No <- rownames(nmds.gg)
#Add Musa species
nmds.gg$Species <- spec.df$Species[match(nmds.gg$Serial.No, spec.df$Serial.No)]
#Add habitat
nmds.gg$Habitat <- metadata$Habitat.summary[match(nmds.gg$Serial.No, metadata$Serial.No)]
#Make dataframe of NMDS species scores for ggplot (species=fungal OTUs)
otu.gg <- as.data.frame(scores(NMDS, "species"))
otu.gg$Genus <- df$Genus[match(rownames(otu.gg), df$OTU)]

#Make dataframe for fitting environmental variables
nmds.envdf <- data.frame(Serial.No=nmds.gg$Serial.No)
#Add seed viability (TZ)
nmds.envdf$TZ <- metadata$TZ[match(nmds.envdf$Serial.No, metadata$Serial.No)]
#Add germination rates
nmds.envdf$Germination <- metadata$ER.shoots[match(nmds.envdf$Serial.No, metadata$Serial.No)] / metadata$ER.embryos[match(nmds.envdf$Serial.No, metadata$Serial.No)]
#Add contamination rates
#nmds.envdf$Contamination <- metadata$ER.contamination[match(nmds.envdf$Serial.No, metadata$Serial.No)] / metadata$ER.seeds[match(nmds.envdf$Serial.No, metadata$Serial.No)]
#Add seed emptiness rates
#nmds.envdf$Emptiness <- metadata$ER.empty[match(nmds.envdf$Serial.No, metadata$Serial.No)] / metadata$ER.seeds[match(nmds.envdf$Serial.No, metadata$Serial.No)]
#Add seed emptiness rates
nmds.envdf$Altitude <- metadata$Altitude[match(nmds.envdf$Serial.No, metadata$Serial.No)]

set.seed(1)
#Fit environmental variables
NMDS.env <- envfit(NMDS, nmds.envdf[-1], type="t", permutations=1000, na.rm=TRUE)
#Make dataframe of TZ fit for ggplot
env.gg <- rbind(data.frame(NMDS.env$vectors$arrows*sqrt(NMDS.env$vectors$r), R=NMDS.env$vectors$r, p=NMDS.env$vectors$pvals), data.frame(NMDS.env$factors$centroids, R=NMDS.env$factors$r, p=NMDS.env$factors$pvals))
rownames(env.gg) <-sub("TZ", "Seed viability", rownames(env.gg))
#env.gg <- env.gg[-which(rownames(env.gg) == "Emptiness"),]

#ANOSIM to test for significant difference between habitats
nmds.anosim <- anosim(species.count.nmds, nmds.gg$Habitat, distance="bray", permutations=1000)
#ANOSIM including all (rare) OTUs
anosim(tree.count[-which(rowSums(tree.count) == 0),], metadata$Habitat.summary[match(rownames(tree.count[-which(rowSums(tree.count) == 0),]), metadata$Serial.No)], distance="bray", permutations=1000)
#ANOSIM to test for signficant difference between Musa species
anosim(species.count.nmds, nmds.gg$Species, distance="bray", permutations=1000)
anosim(tree.count[-which(rowSums(tree.count) == 0),], metadata$Species[match(rownames(tree.count[-which(rowSums(tree.count) == 0),]), metadata$Serial.No)], distance="bray", permutations=1000)

#Plot NMDS
gg.nmds.hab <- ggplot() +
  stat_ellipse(data=nmds.gg,
               aes(x=NMDS1, y=NMDS2, fill=Habitat),
               alpha=0.1,
               geom="polygon",
               type='t',
               size=1) +
  geom_point(data=nmds.gg,
             aes(x=NMDS1, y=NMDS2, colour=Habitat, shape=Species),
             size=3) +
  geom_text_repel(data=otu.gg,
                  aes(x=NMDS1, y=NMDS2, label=Genus, fontface="italic"),
                  size=4) +
  geom_segment(data=env.gg,
               aes(x=0, xend=NMDS1, y=0, yend=NMDS2, alpha=p <= 0.05),
               arrow=arrow(length=unit(0.2, "cm")),
               size=0.8,
               show.legend=FALSE) +
  geom_text(data=env.gg,
            aes(x=NMDS1*1.5, y=NMDS2, label=rownames(env.gg), fontface="bold", alpha=p <= 0.05),
            #direction="x",
            #nudge_x=,
            size=4,
            show.legend=FALSE) +
  geom_text(aes(x=-1.2, y=-3),
            label=paste("stress=",round(NMDS$stress, 2)),
            size=4) +
  geom_text(aes(x=1.5, y=3),
            label=expression(bold("ANOSIM")),
            size=4) +
  geom_text(aes(x=1.5, y=2.75),
            label=paste("R=", round(nmds.anosim$statistic, 1)),
            size=4) +
  geom_text(aes(x=1.5, y=2.5),
            label=paste("p=", round(nmds.anosim$signif, 3)),
            size=4) +
  scale_alpha_discrete(range=c(0.5, 1)) +
  scale_color_manual(values=c("#90E84B","#faff45","#fb9a99","#6FD4FF")) +
  scale_fill_manual(values=c("#90E84B","#faff45","#fb9a99","#6FD4FF")) +
  scale_shape_manual(values=c(15:17),
                     name=expression(paste(italic("Musa")," species")),
                     labels=c(expression(italic("Musa acuminata")), expression(italic("Musa balbisiana")), expression(italic("Musa itinerans")))) +
  guides(fill=FALSE) +
  theme_bw() +
  theme(legend.text.align=0,
        aspect.ratio=1,
        panel.grid=element_blank(),
        legend.position="top",
        legend.box="vertical",
        legend.margin=margin(0,0,0,0),
        plot.title.position="plot") +
  coord_fixed()
plot(gg.nmds.hab)


## FIGURE 4B - ABUNDANCE OF ENDOPHYTES IN HABITATS ##

#Make dataframe of fungi per Musa individual
hab.df <- data.frame(Serial.No=rownames(tree.count), Habitat=metadata$Habitat.summary[match(rownames(tree.count), metadata$Serial.No)], Count=rowSums(tree.count))

#TukeyHSD to test for significance between habitats
hab.tukey <- TukeyHSD(aov(lm(Count ~ Habitat, data=hab.df)))

#Make dataframe for ggplot with tukey groups
hab.labels <- data.frame(multcompLetters(hab.tukey[["Habitat"]][,4])["monospacedLetters"])
hab.labels <- cbind(Treatment=rownames(hab.labels), data.frame(hab.labels, row.names=NULL))

#Add sample size and position for labelling
hab.labels$n <- NA
hab.labels$pos <- NA

for (i in hab.labels$Treatment) {
  hab.labels$n[hab.labels$Treatment == i] <- length(hab.df$Habitat[hab.df$Habitat == i])
  hab.labels$pos[hab.labels$Treatment == i] <- quantile(hab.df$Count[hab.df$Habitat == i])[4]
}

hab.labels$n <- paste0("n=", hab.labels$n)

#Plot OTUs per individual bargraphs
gg.hab <- ggplot(hab.df, aes(x=Count, y=Habitat, fill=Habitat)) +
  geom_boxploth(show.legend=FALSE) +
  geom_text(data=hab.labels,
            aes(x=pos, y=Treatment, label=n),
            vjust=-0.5,
            hjust=-0.25,
            size=2.5,
            inherit.aes=FALSE) +
  geom_text(data=hab.labels,
            aes(x=Inf, y=Treatment, label=monospacedLetters),
            family="mono",
            hjust=-0.5,
            size=3,
            inherit.aes=FALSE) +
  scale_x_continuous(name=expression(paste("Number of OTUs per ", italic("Musa")," individual")),
                     expand=c(0, 0)) +
  scale_fill_manual(values=c("grey","#90E84B","#faff45","grey", "#fb9a99","#6FD4FF")) +
  theme_minimal() +
  theme(axis.title.x=element_text(margin=unit(c(2, 0, 0, 0), "mm"),
                                  size=10),
        axis.title.y=element_blank(),
        axis.text=element_text(size=8),
        plot.title.position="plot",
        plot.margin=unit(c(20,30,5.5,5.5), "pt")) +
  coord_cartesian(clip = 'off')
plot(gg.hab)

#Write to file
tiff(file=paste0("Figure_4-", Sys.Date(), ".tiff"), height=8, width=8, units="in", res=600)
egg::ggarrange(gg.nmds.hab, gg.hab, nrow=2, widths=c(4,1), labels=c("A", "B"), label.args=list(gp=grid::gpar(font=2, cex=1.2)))
dev.off()


###################################
## FIGURE 5 - FUSARIUM PHYLOGENY ##
###################################

#Read in tree
phy <- read.tree("data/RAxML_bipartitions.Musa_fusarium_138T_3loci")
#Root tree
outgroup <- c("Cylindrocarpon_cylindroides_22505", "Neonectria_coccinea_20485")
phy <- root(phy, outgroup, resolve.root=TRUE)
#Remove insignificant bootstrap values
BS <- as.numeric(phy$node.label)
BS[(BS<70)] <- ""
BS[c(1,2)] <- ""
#Replace underscores in tip labels
phy$tip.label <- gsub("_", " ", phy$tip.label)
#Replace FIESC species names
fiesc <- read.csv("data/FIESC_names.csv", header=FALSE)
phy$tip.label[match(fiesc$V1, phy$tip.label)] <- paste0("Fusarium ",fiesc$V2, " ", word(fiesc$V1,-1)) 
#Read in Fusarium EF OTU mapping
fusotus <- read.csv("data/Fus_ef_060120_otus99map2.txt", sep = '\t', col.names=c("Serial.No","OTU"))
#Summarise OTU counts
fusotus <- fusotus %>% group_by(OTU) %>% summarise(n=n())
#Add to tip labels
phy$tip.label[match(fusotus$OTU, phy$tip.label)] <- paste0(phy$tip.label[match(fusotus$OTU, phy$tip.label)], " n=", fusotus$n)
#Create dataframe to signify OTUs from this study to colour
tips <- data.frame(species=phy$tip.label)
tips$colour <- ifelse(grepl("Otu", tips$species), "new", "")
#Create vector for tip colours
tip.colours <- c("black","dodgerblue3")

#Create dataframes for node collapsing, scaling and collapsing of species complexes for ggtree
clades <- data.frame(node=c(229,244,256,208,207,165,191,162,159,257,259,260,141,148,274,272), 
                     name=c("equiseti clade","sambucinum","chlamydosporum","tricinctum","heterosporum","nisikadoi","oxysporum","redolens","concolor","lateritium","buharicum","buxicola","dimerum","albidum","staphyleae","decemcellulare"), 
                     offset=c(0.029,0.15,0.11,0.05,0.038,0.035,0,0,0.15,0.11,0.17,0.085,0.39,0.09,0.165,0.015), 
                     scale=c(0.1,0.1,1,0.5,1,0.5,0.1,1,1,0.7,1,1,0.3,1,1,1))

highlight <- data.frame(node=c(169,213,262), 
                        name=c("fujikuroi","incarnatum-equiseti","solani"), 
                        offset=c(0.55,0.6,0.65))

#Plot tree
gg <- ggtree(phy, cex=0.3) +
  xlim(0,2.5)

for (i in 1:length(clades$node)) {
  gg <- scaleClade(gg, node=clades$node[i], scale=clades$scale[i])
}

for (i in 1:length(clades$node)) {
  gg <- collapse(gg, node=clades$node[i], mode="max", colour="wheat2", fill="wheat2")
  gg <- gg +
    geom_cladelabel(node=clades$node[i],
                    label=clades$name[i],
                    barsize=0,
                    offset=clades$offset[i],
                    fontsize=7,
                    fontface='bold.italic')
}

for (i in 1:length(highlight$node)) {
  gg <- gg +
    geom_cladelabel(node=highlight$node[i],
                    label=highlight$name[i],
                    offset=highlight$offset[i],
                    barsize=1,
                    fontsize=7,
                    fontface='bold.italic') +
    geom_hilight(node=highlight$node[i],
                 extend=highlight$offset[i],
                 alpha=0.3, fill="wheat2")
}

gg1 <- gg %<+% tips +
  geom_tree() + 
  geom_tiplab(aes(color=colour), size=7,  fontface='italic') +
  scale_colour_manual(values=tip.colours) +
  geom_nodelab(aes(x=branch), label=BS, size=4, vjust=-0.5) +
  #geom_text2(aes(subset=!isTip, label=node), hjust=-.3) +
  theme(legend.position="none")
plot(gg1)

#Write to file
tiff(file=paste0("Figure_5-", Sys.Date(), ".tiff"), height=30, width=20, units="in", res=300)
plot(gg1)
dev.off()


######################################################
## SUPPLEMENTARY TABLE 3 - OTUs FOUND OUTSIDE SEEDS ##
######################################################

#Read in all data
df.cont <- read.csv("data/All_data.csv")
#Add OTU to metadata dataframe
df.cont["OTU"] <- otus$OTU[match(df.cont$Serial.No, otus$Serial.No)]
#Add classification to metadata dataframe
df.cont[,c("Phylum","Class","Order","Genus")] <- tbas.class[,c("Phylum","Class","Order","Genus")][match(df.cont$OTU, tbas.class$OTU),]
#Remove accessions which couldn't cluster into OTUs
df.cont <- subset(df.cont, OTU != "#N/A")
#Filter for imprints
df.cont <- subset(df.cont, df.cont$Sterilised != "Y")
#Remove unnecessary columns
df.cont <- df.cont[44:48]
#Add number of times found in imprints
df.cont$Count <- NA

for (i in 1:length(df.cont$OTU)) {
  df.cont$Count[i] <- length(grep(df.cont$OTU[i], df.cont$OTU))
}

#Remove duplicate rows
df.cont <- unique(df.cont)
write.csv(df.cont, "Supplementary_Table_3.csv", row.names=FALSE)


####################
## STATS FOR TEXT ##
####################

## Proportion of endophytes reported as other guilds in FUNguild ##

#Get function to download FUNGuild database from https://rdrr.io/github/vmikk/metagMisc/man/parse_funguild.html
parse_funguild <- function(url = 'http://www.stbates.org/funguild_db.php', tax_name = TRUE){
  
  # require(XML)
  # require(jsonlite)
  # require(RCurl)
  
  ## Parse data
  tmp <- XML::htmlParse(url)
  tmp <- XML::xpathSApply(doc = tmp, path = "//body", fun = XML::xmlValue)
  
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
      TaxID = c(0, 3:13, 15:26),
      Taxon = factor(taxons, levels = taxons))
    
    ## Match taxon codes
    db$taxonomicLevel <- taxmatch[match(x = db$taxonomicLevel, table = taxmatch$TaxID), "Taxon"]
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

#Total number of endophytes
print(length(df$Serial.No))
#Total number of OTUs
print(length(unique(df$OTU)))
#Number of singleton OTUs
singletons <- data.frame(otu=unique(df$OTU), count=NA)
for (i in 1:length(singletons$otu)) {
  singletons$count[i] <- length(grep(singletons$otu[i], df$OTU))
}
print(paste0(length(singletons$otu[singletons$count == 1]), ", ", round(length(singletons$otu[singletons$count == 1]) / length(singletons$otu) * 100), "%"))

#Percentage of endophytes that were Lasiodiplodia, Aspergillus and Fusarium
for (i in c("Lasiodiplodia", "Aspergillus", "Fusarium")) {
  print(paste0(i, ", ", length(grep(i, df$Genus)), ", ", round(length(grep(i, df$Genus)) / length(df$Genus) * 100), "%"))
}




genera <- data.frame(genus=unique(df$Genus), count=NA)

for (i in 1:length(genera$genus)) {
  genera$count[i] <- length(df$Genus[df$Genus == genera$genus[i]])
}

ggplot(genera, aes(label=genus, size=count)) +
  geom_text_wordcloud()


test <- spec.df
test$otus <- NA

for (i in 1:length(spec.df$Serial.No)) {
  test$otus[i] <- length(unique(df$OTU[grep(test$Serial.No[i], df$Serial.No)]))
}

test.labels <- test %>% group_by(Species) %>% summarise(value=max(otus), count=n())


ggplot(test, aes(x=Species, y=otus)) +
  geom_violin() +
  geom_boxplot(width=0.2) +
  geom_text(data=test.labels, nudge_y=1, aes(x=Species, y=value, label=paste0("n=", count)))