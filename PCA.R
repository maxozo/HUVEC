library(PCAtools)
library("factoextra")
library(ade4)
#/Users/mo11/work/HUVEC/Data3/min_max/Data_Extracted_Peaks/Metrics_Calculations.csv

treatment='Thrombin'
measurement='Capacitance'
norm_method = 'FC_controls'
freq=4000

Experimental_grops_load_peaks_pre = read.table(paste('/Users/mo11/work/HUVEC/Data3/',norm_method,'/Data_Extracted_Peaks/',treatment,'_Extracted_Peaks_',measurement,'.csv',sep=''), fill = TRUE,row.names = 1,header = TRUE,sep = ',')
Experimental_grops_load_features_pre = read.table(paste('/Users/mo11/work/HUVEC/Data3/',norm_method,'/Data_Extracted_Peaks/',treatment,'_Metrics_Calculations_',measurement,'.csv',sep=''), fill = TRUE,row.names = 1,header = TRUE,sep = ',')

Experimental_grops_load_peaks = t(Experimental_grops_load_peaks_pre[which(Experimental_grops_load_peaks_pre['freq'] == freq),])
Experimental_grops_load_features = Experimental_grops_load_features_pre[which(Experimental_grops_load_features_pre['freq'] == freq),]
#Experimental_grops_load_features_normalised = read.table('/Users/mo11/work/HUVEC/Data3/FC_norm/Data_Extracted_Peaks/Norm_Metrics_Calculations.csv',, fill = TRUE,row.names = 1,header = TRUE,sep = ',')

#Experimental_grops_load_peaks = t(read.table('/Users/mo11/work/HUVEC/Data3/Thrombin_Data_metadata.tsv',, fill = TRUE,header = TRUE,sep = ','))
#ix <- which(c('X', 'freq') == rownames(Experimental_grops_load_peaks))
#Experimental_grops_load_peaks = Experimental_grops_load_peaks[-c(ix), ]   

Experimental_grops_load = Experimental_grops_load_features
rownames(Experimental_grops_load)=gsub('X', '',rownames(Experimental_grops_load))
Metadata = read.table(paste('/Users/mo11/work/HUVEC/Data3/',treatment,'_Data_metadata.tsv',sep=''), fill = TRUE,row.names = 1,header = TRUE,sep = '\t')

rownames(Metadata)=gsub(' ', '.',rownames(Metadata))
t2 = transform(merge(Experimental_grops_load, Metadata, by=0,all.x = TRUE,), row.names=Row.names, Row.names=NULL)
d1 = dim(t2)[2]-dim(Metadata)[2]
Exp_data = t2[0:d1]

Exp_data[is.na(Exp_data)] <- 100000
drops=c('flagged.as.outlier')
Exp_data2 = Exp_data[ , !(names(Exp_data) %in% drops)]
res.pca <- dudi.pca(Exp_data2,
                    scannf = FALSE,center = FALSE, scale = FALSE,  # Hide scree plot
                    nf = 5            # Number of components kept in the results
)
#Exp_data2[which(!is.finite(Exp_data2))] <- 0

res.pca <-princomp(Exp_data)
res.pca <-prcomp(Exp_data,scale=True)
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



