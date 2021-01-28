#load packages
library(Seurat) #needed to deal with seurat object


#parameter
celltype='Tcells28'
inputfile='Tcells28.data.R'
outputfile=paste0(celltype,'.csv')




#load inputfile
load(inputfile)
expr=GetAssayData(object=Tcells28, assay="RNA", slot="data")


#output as matrix
write.table(expr,outputfile,quote=F,sep=',')