#make design matrix for Illumina BeadArray

design.Matrix.prepare<-function(data){
        #making design matrix
        descr<-pData(data)$Sample.ID #define sample names
        #define sample names factor
        #need to do that as simple levels(descr) sorts levels in alphabetical order
        descrFac<-factor(descr,levels=unique(descr)) 
        designMatrix<-model.matrix(~0+descrFac) #make the design matrix
        colnames(designMatrix)<-unique(descr) #and edit its column names
        return(designMatrix)
}