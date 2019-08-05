#################################################
##### SUPPLEMENTARY DATA 1 FROM MELLO ET AL. 2019
#################################################


##### Ecological Synthesis Lab (SintECO): https://marcomellolab.wordpress.com
##### Authors: Marco A. R. Mello, Gabriel M. Felix, Rafael B. P. Pinheiro, Renata L. Muylaert, Cullen Geiselman, Sharlene E. Santana, Marco Tschapka, Nastaran Lotfi, Francisco A. Rodrigues & Richard D. Stevens
##### E-mail: marmello@gmail.com 
##### How to cite: Mello et al. 2019. Assembly rules of a continent-wide multilayer network. Published as a preprint in bioRxiv in 2018, accepted by Nature Ecology and Evolution in 2019. Supplementary Data 1. 
##### Updated on August 5th, 2019 (English version).
##### Run in R version 3.6.0 (2019-04-26) -- "Planting of a Tree"

##### Disclaimer: You may use this script freely for non-comercial purposes at your own risk. We assume no responsibility or liability for the use of this software, convey no license or title under any patent, copyright, or mask work right to the product. We reserve the right to make changes in the software without notification. We also make no representation or warranty that such application will be suitable for the specified use without further testing or modification. If this script helps you produce any academic work (paper, book, chapter, dissertation etc.), please acknowledge the authors and cite the source.


##############################################
##### Set the working directory
##############################################


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


##############################################
##### Load the packages 
##############################################


library("igraph")
library("bipartite")
library("reshape2")
library("gdata")
library("Matrix")
library("paco")
library("ggplot2")
library("tcltk")
library("corrgram")
library("glm2")
library("lme4")
library("gridExtra")
library("dplyr")


##############################################
##### Create multilayer network for igraph 
##############################################


#Create nodes and edges
nodes <- read.csv("nodesc1.csv", header=T, as.is=T)
links <- read.csv("linksc1.csv", header=T, as.is=T)


#Inspect objects
class(nodes)
head(nodes)
class(links)
head(links)


#Identify and label dual links
multi_index=which(duplicated(links[,1:2], fromLast = T))
links[c(multi_index),][,4]="dual"
to_remove=which(duplicated(links[,1:2]))
links=links[-to_remove,]
links = arrange(links,type)
unique(links$type)


#Create multilayer network
multilayer <- graph_from_data_frame(d=links, vertices=nodes, directed=F)


#Add information on two-mode strucutre
V(multilayer)$type <- V(multilayer)$name %in% links[,1]
V(multilayer)$type = nodes[,2]


#Inspect multilayer network
class(multilayer)
multilayer

E(multilayer)$weight
E(multilayer)$type

V(multilayer)$name
V(multilayer)$taxon


#Import module membership from DIRT_LPA+
#You should run this analysis using the package biparite for R.
#Then you need to extract module membership and save it as 
#tab-delimited TXT file.
modules=read.table("partitions.txt", h=T)
head(modules)

length(modules$nodes)
length(V(multilayer)$name)

#######################################################################################
###### Draw the graphs
#######################################################################################


#Set the same layout for all graphs
l <- layout_nicely(multilayer)


#Set node shapes for all vertices
V(multilayer)$shape = V(multilayer)$taxon
V(multilayer)$shape = gsub("Bats","square",V(multilayer)$shape)
V(multilayer)$shape = gsub("Plants","circle",V(multilayer)$shape)


#######################################################################################
###### Draw the graphs
#######################################################################################


##### Draw the graph with node colors by taxon #####

V(multilayer)$color = V(multilayer)$taxon
V(multilayer)$color = gsub("Bats","#D3802B",V(multilayer)$color)
V(multilayer)$color = gsub("Plants","#2C8437",V(multilayer)$color)

E(multilayer)$color = E(multilayer)$type
E(multilayer)$color = gsub("dual","#A151A1",E(multilayer)$color)
E(multilayer)$color = gsub("frugivory","#B3DCF9",E(multilayer)$color)
E(multilayer)$color = gsub("nectarivory","#FFD26B",E(multilayer)$color)


E(multilayer)$arrow.mode = 0

E(multilayer)$width = E(multilayer)$type
E(multilayer)$width = gsub("dual", 6, E(multilayer)$width)
E(multilayer)$width = gsub("frugivory", 1, E(multilayer)$width)
E(multilayer)$width = gsub("nectarivory", 1, E(multilayer)$width)

png(filename= "multilayer.bip.png", res= 300, height= 3000, width= 4900)

par(mfrow=c(1,1),mar=c(1,1,3,17))
plot(multilayer, 
     vertex.color = V(multilayer)$color, 
     vertex.frame.color= V(multilayer)$color, 
     vertex.shape = V(multilayer)$shape, 
     vertex.size=3,
     vertex.label=V(multilayer)$name,
     vertex.label.color="white",
     vertex.label.cex=.1,
     edge.color = E(multilayer)$color, 
     edge.curved=0.3, 
     layout=l)

legend(x = 1.3, y = 0.9, title="Taxon",
       legend = c("Bats", "Plants"), pch = c(15,16),   
       text.col = "gray20", title.col = "black", 
       box.lwd = 0, cex = 2, col=c("#D3802B", "#2C8437"))
legend(x = 1.3,y = 0.1, title="Layers (interaction types)",
       legend = c("frugivory", "nectarivory", "dual"),
       fill = c("#B3DCF9", "#FFD26B", "#A151A1"), border = "white", 
       text.col = "gray20", title.col = "black", 
       box.lwd = 0, cex = 2, bty="n")
title(main = "Bat-plant multilayer network", cex.main=2)
par(mfrow=c(1,1))

dev.off()


#####  Draw the graph with node colors by module #####

#First, you need to analyze the modularity of the network using the package bipartite, and then extract infromation on module membership, and then saving the membership list as a tab-delimited TXT file, whcih should be imported before as the object "modules".
vertex_name_order = data.frame(nodes=V(multilayer)$name)
modules_order = merge(vertex_name_order, modules, by="nodes", sort = F)
#colrs = rainbow(12, alpha=1) #rainbow palette with 12 colors
colrs = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#B41572", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928") #colorblind-friendly palette with 12 colors.
V(multilayer)$color <- colrs[modules_order$modules]

E(multilayer)$color = E(multilayer)$type
E(multilayer)$color = gsub("dual","#A151A1",E(multilayer)$color)
E(multilayer)$color = gsub("frugivory","#B3DCF9", E(multilayer)$color)
E(multilayer)$color = gsub("nectarivory","#FFD26B", E(multilayer)$color)


E(multilayer)$arrow.mode = 0

E(multilayer)$width = E(multilayer)$type
E(multilayer)$width = gsub("dual", 6, E(multilayer)$width)
E(multilayer)$width = gsub("frugivory", 1, E(multilayer)$width)
E(multilayer)$width = gsub("nectarivory", 1, E(multilayer)$width)

png(filename= "multilayer.mod.png", res= 300, height= 3000, width= 4900)

par(mfrow=c(1,1),mar=c(1,1,3,17))
plot(multilayer, 
     vertex.color = V(multilayer)$color, 
     vertex.frame.color= V(multilayer)$color, 
     vertex.shape = V(multilayer)$shape, 
     vertex.size=3,
     vertex.label=V(multilayer)$name,
     vertex.label.color="white",
     vertex.label.cex=.2,
     edge.color = E(multilayer)$color, 
     edge.curved=0.3, 
     layout=l)

legend(x = 1.3,y = 0.9, title="Taxon",
       legend = c("Bats", "Plants"), pch = c(15,16),   
       text.col = "gray20", title.col = "black", 
       box.lwd = 0, cex = 2, col=c("grey", "grey"))
legend(x = 1.3,y = 0.1, title="Layers (interaction types)",
       legend = c("frugivory", "nectarivory", "dual"),
       fill = c("#B3DCF9", "#FFD26B", "#A151A1"), border = "white", 
       text.col = "gray20", title.col = "black", box.lwd = 0, cex = 2)
legend(x = 1.3,y = -0.6, title="Modules",
       legend = c("node colors = modules"),
       fill = c("grey"), border = "white", 
       text.col = "gray20", title.col = "black", box.lwd = 0, cex = 2)
title(main = "Bat-plant multilayer network", cex.main=2)
par(mfrow=c(1,1))

dev.off()
