#Illumina BeadArray data
#fitting to a linear regression model and Bayes processing

fit.data<-function(eset,designMatrix,contMatrix){
        fit<-lmFit(exprs(eset),designMatrix)
        fit2<-contrasts.fit(fit,contMatrix)
        fit2Eb<-eBayes(fit2)
        fit2Eb$genes$ProbeID<-fData(eset)$ProbeID
        fit2Eb$genes$SYMBOL<-fData(eset)$SYMBOL
        fit2Eb$genes$GENENAME<-fData(eset)$GENENAME
        fit2Eb$genes$ENTREZID<-fData(eset)$ENTREZID
        
        return(fit2Eb)
}


