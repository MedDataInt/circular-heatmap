###### Gene Expression Circular Heatmap in Wild Type and Knock-Out Mice 
# Load libraries 
library(tidyverse)
library(dendextend)
library(circlize)
library(ComplexHeatmap)

# load data 
gene_data <-read.csv("GeneList_Jie.csv", sep = ',' ,check.names = F, row.names = 1)
head(gene_data)

# select the data, convert rows to columns 
mat <- t(gene_data[,1:4])  
head(mat)

# normalize data 
mat <- scale(mat, center = TRUE, scale = TRUE) 
View(mat)

# calculate dendextend
dend <-as.dendrogram(hclust(dist(t(mat))))  

# Set the number of dend clusters 
dend <- dend %>% set("branches_k_color",k = 3)  

# plot dend
par(mar=c(7.5,3,1,0))   
plot(dend)

# order the gene to match with dend
mat2 <- mat[,order.dendrogram(dend)]  

# label and number row names and col names
label1 <- row.names(mat2)
label1

label2 <- colnames(mat2)
label2

nr <- nrow(mat2)
nr

nc <- ncol(mat2)
nc

# convert value to colors
col_fun <- colorRamp2(c(-1, 0, 1), c("skyblue", "white", "red"))   
col_mat <- col_fun(mat2)  

# view the color showing in each cells
col_mat[,1] 

# set a space for the circular heatmap, 
par(mar=c(0,0,0,0))

circos.par(canvas.xlim=c(-1.5,1.5),        
           canvas.ylim=c(-1.5,1.5),
           cell.padding=c(0,0,0,0),
           gap.degree=90) 

# initialize the figure  
factors = "a"
circos.initialize(factors,xlim=c(0,ncol(mat2))) 

# map heatmap data into color tiles 
circos.track(ylim=c(0,nr), bg.border = NA, track.height = 0.1*nr,
             panel.fun = function(x,y){
               for (i in 1:nr) {
                 circos.rect(xleft = 1:nc-1, 
                             ybottom = rep(nr-i, nc),
                             xright = 1:nc, 
                             ytop = rep(nr-i+1, nc),
                             border="white",
                             col=col_mat[i,])
                 
                 circos.text(x=nc,
                             y=4.4-i,
                             labels=label1[i],
                             facing = "downward", 
                             niceFacing = TRUE,
                             cex = 0.8,   
                             adj = c(-0.2, 0))
               }
             })
    

# add annotations 
for (i in 1:nc){
  circos.text(x = i-0.4,
              y = 5,   # circle for gene name , you can adjust to larger, if needed
              labels = label2[i],
              facing = "clockwise",
              niceFacing = TRUE,
              cex = 1,        
              adj = c(0,0)
              )
            }

# merge dendrogram cluster tree with circular map  
max_height<-max(attr(dend,"height"))
circos.track(ylim=c(0, max_height), 
             bg.border = NA, 
             track.height = 0.3,
             panel.fun=function(x,y){
               circos.dendrogram(dend = dend,
                                 max_height = max_height)
                                 }
            )

circos.clear() 

# add legend key 
lgd <- Legend(at = c(-2,-1,0,1,2), 
              col_fun = col_fun, # col_fun was set above 
              title_position = "topcenter", 
              title = "Z-score")
draw(lgd, x = unit(0.7,"npc"), y=unit(0.7,"npc"))
#---------------------------End---------------------------------#