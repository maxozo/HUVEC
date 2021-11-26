library(PCAtools)
library("factoextra")
library(ade4)
Experimental_grops_load = t(read.table('/Users/mo11/work/HUVEC/Data2/Data_Extracted_Peaks/Extracted_Peaks_3e.csv',, fill = TRUE,row.names = 1,header = TRUE,sep = ','))

res.pca <- dudi.pca(Experimental_grops_load,
                    scannf = FALSE,   # Hide scree plot
                    nf = 5            # Number of components kept in the results
)

fviz_eig(res.pca)
fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


means <- apply(Experimental_grops_load,2,mean)
sds <- apply(Experimental_grops_load,2,sd)
nor <- scale(Experimental_grops_load,center=means,scale=sds)

distance = dist(nor)
mydata.hclust = hclust(distance)
plot(mydata.hclust)



