#analysis of Illumina BeadArray data
#read raw bead level data from *.idat files
#normalize, add annotations, filter out bad probes

read.norm.filt.data<-function(writeGEO=FALSE,softTemplate=NULL){
        # ALL *.idat files pertaining to experiment must be in one folder!
        dir<-choose.dir(default = getwd(),caption = "Choose directory with IDAT files.")
        files<-list.files(dir,full.names = TRUE,pattern=".idat")
        data<-readIdatFiles(idatFiles = files)
        
        # make description vectors and add them to phenoData
        phenDF<-read.table(file=choose.files(caption="Select a file with phenoData"),
                           sep="\t")
        phenoData(data)<-AnnotatedDataFrame(data=cbind(pData(data),phenDF))
        
        # normalize the data
        data.norm<-normaliseIllumina(data, method="quantile", transform= "log2")
        
        if (writeGEO) {
                #write the normalized and raw data as GEO submission files
                makeGEOSubmissionFiles(rawData = data, normData = data.norm,
                                       softTemplate=softTemplate)
        }
        
        # write feature names into a variable 
        fNamesILMN<-featureNames(data.norm)
        
        # add information from illuminaMousev2 to annotation data frame
        SYMBOL<-unlist(lookUp(fNamesILMN,
                              "illuminaMousev2","SYMBOLREANNOTATED")) 
        GENENAME<-unlist(lookUp(fNamesILMN,
                                "illuminaMousev2","GENENAME"))
        PROBEQUALITY<-unlist(lookUp(fNamesILMN,
                                    "illuminaMousev2","PROBEQUALITY"))
        ENTREZID<-unlist(lookUp(fNamesILMN,
                                "illuminaMousev2","ENTREZID"))
        annot<-as.data.frame(cbind(SYMBOL,GENENAME,ENTREZID,PROBEQUALITY))
        
        #add annotation table to the fData
        fData(data.norm)<-cbind(fData(data.norm),annot) 
        
        # filter out the bad probes
        qual<-fData(data.norm)$PROBEQUALITY
        #table(qual)
        removeProbes<-qual=="No match" | qual=="Bad"
        data.norm.filt<-data.norm[!removeProbes,]
        
        return(data.norm.filt)
}