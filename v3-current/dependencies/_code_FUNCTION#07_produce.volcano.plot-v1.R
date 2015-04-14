#produce volcano plots for a given Bayes fit

produce.volcano.plot<-function(fit2Eb,
                               PVALUE=1,
                               LOGFC=0,
                               markSHAPE=FALSE,
                               markCOLOR=FALSE){
        for (i in 1:ncol(fit2Eb)) { 
                #make a dataframe with log fold change values, and adjusted p-values
                tempLFC<-fit2Eb$coefficients[,i]
                tempAdjPval<-p.adjust(fit2Eb$p.value[,i],method="BH")
                tempLODS<-fit2Eb$lods[,i]
                tempNegLogAdjPval<- -log10(tempAdjPval)
                flag<-ifelse(tempAdjPval>PVALUE, 16, ifelse(abs(tempLFC)<LOGFC, 16,8))
                #define shapes for plot
                #16 - filled circle, 8 - star
                if(length(unique(flag))==1) {
                        shapes<-16 
                        clrs<-"gray24" 
                } else {
                        shapes<-c(8,16)
                        clrs<-c("red3","gray24")
                }
                
                #assemble dataframe for volcano plot
                tempDF<-as.data.frame(cbind(tempLFC,tempLODS,tempNegLogAdjPval,flag))
                
                for (j in 1:2){
                        subType<-ifelse(j==1,"LOdds.","negLogPval.")
                        y.data<-ifelse(j==1,"tempLODS","tempNegLogAdjPval")
                        file.name<-paste0("Volcano Plot.",subType,
                                          colnames(fit2Eb)[i],".pdf")
                        
                        make.volcano(file.name=file.name,
                                     dFrame=tempDF,
                                     x="tempLFC",
                                     y=y.data,
                                     shapes=shapes,
                                     clrs=clrs,
                                     markSHAPE=markSHAPE,
                                     markCOLOR=markCOLOR)
                }  
        }
}

