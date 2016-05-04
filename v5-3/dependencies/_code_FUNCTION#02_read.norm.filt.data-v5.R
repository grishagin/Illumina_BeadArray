#analysis of Illumina BeadArray data
#read raw bead level data from *.idat files
#normalize, add annotations, filter out bad probes

read.norm.filt.data<-function(batchCorrect=FALSE,
                              writeGEO=FALSE){
        # ALL *.idat files pertaining to experiment must be in one folder!
        dir<-tk_choose.dir(default = getwd(),caption = "Choose directory with IDAT files.")
        files<-list.files(dir,full.names = TRUE,pattern=".idat")
        
        # make description vectors and add them to phenoData
        phenDF<-xlsx:::read.xlsx(file = tk_choose.files(caption="Select xlsx file with phenoData"),
                                 sheetIndex = 1,
                                 header = TRUE,
                                 stringsAsFactors=FALSE)
        #if need GEO, ask for a soft template
        if (writeGEO) {
                softTemplate<-
                        choose.files(caption = "Select SOFT file to annotate against (if available).")
                if (softTemplate==0) {softTemplate=NULL}
        }
        
        #load the idat files
        data<-readIdatFiles(idatFiles = files)
        
        #order phenDF such that it matches the existing pheno data
        row_order<-match(pData(data)$sectionNames,phenDF$Chip.ID)
        if(sum(is.na(row_order)>0)){
                print("phenoData does not seem to match chipIDs extracted from file names.")
                readline(prompt="Hit 'Enter' to abort.")
				stop()
        }
        phenDF<-phenDF[row_order,]
        
        #add custom pheno data to dataset
        phenoData(data)<-AnnotatedDataFrame(data=cbind(pData(data),phenDF))

        # normalize the data
        data.norm<-normaliseIllumina(data, method="quantile", transform= "log2")
        
        #perform batch correction
        if(batchCorrect){
                designM<-model.matrix(~0+Design, pData(data))
                exprs(data.norm)<-removeBatchEffect(exprs(data.norm),
                                               batch=pData(data.norm)$Batch,
                                               design=designM)
        }
        
        #if GEO files required, make them
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