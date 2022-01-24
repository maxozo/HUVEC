library(PCAtools)
library("factoextra")
library(ade4)
#/Users/mo11/work/HUVEC/Data3/min_max/Data_Extracted_Peaks/Metrics_Calculations.csv
Experimental_grops_load_peaks = t(read.table('/Users/mo11/work/HUVEC/Data3/FC_norm/Data_Extracted_Peaks/Extracted_Peaks.csv',, fill = TRUE,row.names = 1,header = TRUE,sep = ','))
Experimental_grops_load_features = read.table('/Users/mo11/work/HUVEC/Data3/FC_norm/Data_Extracted_Peaks/Metrics_Calculations.csv',, fill = TRUE,row.names = 1,header = TRUE,sep = ',')
Experimental_grops_load_features_normalised = read.table('/Users/mo11/work/HUVEC/Data3/FC_norm/Data_Extracted_Peaks/Norm_Metrics_Calculations.csv',, fill = TRUE,row.names = 1,header = TRUE,sep = ',')


Experimental_grops_load = Experimental_grops_load_features_normalised
rownames(Experimental_grops_load)=gsub('X', '',rownames(Experimental_grops_load))
Metadata = read.table('/Users/mo11/work/HUVEC/Data3/Thrombin_Data_metadata.tsv', fill = TRUE,row.names = 1,header = TRUE,sep = '\t')

rownames(Metadata)=gsub(' ', '.',rownames(Metadata))
t2 = transform(merge(Experimental_grops_load, Metadata, by=0,all.x = TRUE,), row.names=Row.names, Row.names=NULL)
d1 = dim(t2)[2]-dim(Metadata)[2]
Exp_data = t2[0:d1]

Exp_data[is.na(Exp_data)] <- 100000
res.pca <- dudi.pca(Exp_data,
                    scannf = FALSE,   # Hide scree plot
                    nf = 5            # Number of components kept in the results
)

#Colored based on he time
fviz_pca_ind(res.pca,axes = c(1, 2),col.ind = t2$Date_time, ,repel = TRUE) +
  scale_color_gradient2(low="white", mid="blue",
                        high="red", midpoint=mean(t2$Date_time))+ theme_minimal()

#Colored based on the experiment
fviz_pca_ind(res.pca,axes = c(1, 2),col.ind = "cos2", ,repel = TRUE, habillage=t2$Group)

#Colored based on the experiment type
fviz_pca_ind(res.pca,axes = c(1, 2),col.ind = "cos2", ,repel = TRUE, habillage=t2$exp_id)

#Colored based on the experiment type
fviz_pca_ind(res.pca,axes = c(1, 2),col.ind = "cos2", ,repel = TRUE, habillage=t2$label)

Experimental_grops_load[is.na(Exp_data)] <- 100000
res.pca <- prcomp(t(Experimental_grops_load),  scale = TRUE)

fviz_eig(res.pca)

means <- apply(Experimental_grops_load,2,mean)
sds <- apply(Experimental_grops_load,2,sd)
nor <- scale(Experimental_grops_load,center=means,scale=sds)
distance = dist(nor)
mydata.hclust = hclust(distance)
plot(mydata.hclust)



