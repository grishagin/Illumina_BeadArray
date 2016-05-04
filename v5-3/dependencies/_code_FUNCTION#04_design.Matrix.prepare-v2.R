#make design matrix for Illumina BeadArray

design.Matrix.prepare<-function(data){
        #making design matrix
        descrFac<-factor(pData(data)$Design) #define sample names
        #define sample names factor
        designMatrix<-model.matrix(~0+descrFac) #make the design matrix
        #and edit its column names
        #reorder unique Sample.ID's in the order of unique
        #descriptors of the experiment design
        #they basically prescribe the correct order in the design matrix
        colnames(designMatrix)<-
                unique(pData(data)$Sample.ID)[unique(descrFac)]
        return(designMatrix)
}