#Illumina BeadArray data analysis
#select genes to make heatmap

#make a list of genes, whose whose LFC is greater than specified
#and adjusted expression p-value is less than specified
#by adjustment of p-values in eBayes-processed fit
#using p.adjust (although topTable does produce adjusted p-values, eBayes does NOT)
#IMPORTANT!!! Method acts on a VECTOR of p-values
#and its output DEPENDS on the WHOLE vector 

produce.heatmaps<-function(data,    #data.norm.filt
                           controls,
                           treatments,
                           fit2Eb,
                           contrastsMatrix,
                           LOGFC,
                           PVALUE,
                           forceLabelsON){
        
        contrastsMatrix.colnames<-colnames(contrastsMatrix)
        for (i in 1:length(contrastsMatrix.colnames)) {
                #adjust p-values from Bayes fit using FDR method
                adjPVAL<-p.adjust(fit2Eb$p.value[,contrastsMatrix.colnames[i]],method="BH")
                #select only features with p-values which are lower than preset 
                tempPVAL<-adjPVAL<PVALUE
                #log fold change values for a given coefficient
                logFoldChange<-fit2Eb$coefficients[,contrastsMatrix.colnames[i]]
                #select only features with abs LFC higher than preset
                tempLFC<-abs(logFoldChange)>LOGFC
                
                #make and populate logical vector that is TRUE for genes satisfying 
                #p-value cutoff and LFC cutoff
                #also store all adjusted p-values in one vector for all contrasts
                #also store all lfc values in one vector for all contrasts
                if (i==1){
                        selected<-tempPVAL & tempLFC               
                        selGeneAdjPval<-adjPVAL[tempPVAL]
                        selGeneLogFC<-logFoldChange[tempLFC]
                } else {
                        selected<-selected|(tempPVAL & tempLFC)
                        selGeneAdjPval<-c(selGeneAdjPval,adjPVAL[tempPVAL])
                        selGeneLogFC<-c(selGeneLogFC,logFoldChange[tempLFC])
                }
                
                #print(sum(selected))
        }
        
        phenData<-pData(data)
        samples.sel<-unique(c(controls,treatments))
        #NOTE:  "sectionNames" is automatically generated
        columns.sel<-as.character(phenData$sectionNames[which(phenData$Sample.ID %in% 
                                                                      colnames(designMatrix)[samples.sel])])
        data.sel<-data[selected,columns.sel]
        
        #generate a list of four entities: for no sorting, for sorting by pvalue, and by LOGFC
        dataList.for.heatmap<-prep.heatmap.sort.data(data.sel=data.sel,   #data.norm.filt.sel
                                                     selGeneAdjPval=selGeneAdjPval,
                                                     selGeneLogFC=selGeneLogFC,
                                                     forceLabelsON=forceLabelsON)
        #if class of dataList.for.heatmap is not matrix, abort
        if (class(dataList.for.heatmap)!="list") return()
        
        #########################PLOT HEATMAPS#########################
        #loop trough the dataLabel element of dataList.for.heatmap
        for (i in 1:(length(dataList.for.heatmap[[4]]))){
                #indicate how to sort the genes on heatmap
                SORTED<-ifelse(i==1,"NONE",ifelse(i==2, "negLogPVAL","LOGFC"))
                #indicate which row labels to use
                rowLabels<-ifelse(dataList.for.heatmap[[4]][[1]]==FALSE,
                                  FALSE,
                                  dataList.for.heatmap[[4]][[i]])
                
                #print(rowLabels)

                Rowv<-ifelse(i==1,TRUE,FALSE)
                dendrogram<-ifelse(i==1,"both","column")
                
                file.name=paste0("heatmap",
                                 "-sortedBY-",
                                 SORTED,
                                 ".pdf")
                plot.heatmap(file.name=file.name,
                             data=dataList.for.heatmap[[i]],
                             rowLabels=rowLabels,
                             Rowv=Rowv,
                             dendrogram=dendrogram)
        } 

              
                
}