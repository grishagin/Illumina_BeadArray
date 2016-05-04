######################################
#2015-03-16
#fully tested
#WORKS
######################################
#######SET UP KEY PARAMETERS##########
        #analysis p-value and log fold change
PVALUEvector<-c(0.03)
LOGFCvector<-c(0.5)
        #numbers of columns in design matrix corresponding to controls
controls<-c(1) 
        #numbers of columns in design matrix corresponding to treatments
treatments<-list(c(2,3))
        #sort topTable by LFC or p-value?
        #specify "logFC" or "B" or "p"
sortBY<-"p" 
        #specify main folder for the code
main.code.folder<-"D:\\_\\Templates_and_Scripts\\R\\Illumina_BeadArray"
        #force labels on heatmaps with more than 25 genes?
forceLabelsON=TRUE
######################################

#source all functions for auxiliary files
source(paste0(main.code.folder,"\\_code_FUNCTION#0__source.all.functions-v1.R"))
source.all.functions()

#load all necessary packages
loadPackages()

#load and normalize the data form idat files
data.norm.filt<-read.norm.filt.data()

#produce QC on selected genes on filtered data
QC(data.norm.filt)

        #make design and contrasts matrices in automatic mode
#CONTROLS AND TREATMENTS ARE SET UP in THE BEGINNING of this code (above)!
desM.contM.newCont.list<-design.contrasts.Matrix.prepare(data=data.norm.filt,
                                                         controls=controls,
                                                         treatments=treatments)
designMatrix<-desM.contM.newCont.list[[1]]
contMatrix<-desM.contM.newCont.list[[2]]
newCont<-desM.contM.newCont.list[[3]]

#scan throu all p-values and log fold changes as specified above
#make directory for each combination to store analyses in
for (PVALUE in PVALUEvector){
        for(LOGFC in LOGFCvector){
                #create a working directory for a given p-value and lfc
                dirName<-paste0("p",PVALUE,"-lfc",LOGFC)
                dir.create(dirName, showWarnings = FALSE)
                setwd(dirName)
        
                #fitting: linear regression, then Bayes
                fit2Eb<-fit.data(eset=data.norm.filt,
                                  designMatrix=designMatrix,
                                  contMatrix=contMatrix)
                
                results<-as.data.frame(decideTests(fit2Eb,
                                                   adjust.method="BH",
                                                   p.value = PVALUE,
                                                   lfc =LOGFC))
                
                # make volcano plots - highlight top genes; provide list of gene symbols
                for (i in 1:ncol(fit2Eb)) { 
                        #make a dataframe with log fold change values, and adjusted p-values
                        tempLFC<-fit2Eb$coefficients[,i]
                        tempAdjPval<-p.adjust(fit2Eb$p.value[,i],method="BH")
                        tempLODS<-fit2Eb$lods[,i]
                        tempNegLogAdjPval<--log10(tempAdjPval)
                        tempDF<-as.data.frame(cbind(tempLFC,tempLODS,tempNegLogAdjPval))
                                #colFlag<-ifelse(tempAdjPval>PVALUE, 16, ifelse(abs(tempLFC)<LOGFC, 16,8))
                                #define shapes for plot 16 - filled circle
                                #if(length(unique(colFlag))<2) shapes<-16 else shapes<-c(16,8)
                        
                        for (j in 1:2){
                                subType<-ifelse(j==1,"LOdds.","negLogPval.")
                                y.data<-ifelse(j==1,"tempLODS","tempNegLogAdjPval")
                                file.name<-paste0("Volcano Plot.",subType,".",
                                                 colnames(fit2Eb)[i],".pdf")
                                
                                make.volcano(file.name=file.name,
                                             df=tempDF,
                                             x="tempLFC",
                                             y=y.data)
                        }  
                }
                
                #print venn diagrams for up-, down-, and both up- and down-regulated genes
                #into respective files
                vennD(controls=controls,
                      treatments=treatments,
                      designMatrix=designMatrix,
                      results=results)
                
                #make the lists of top genes
                top.tables(controls=controls,
                           fit2Eb=fit2Eb,
                           contMatrix=contMatrix,
                           newCont=newCont,
                           LOGFC=LOGFC,
                           PVALUE=PVALUE,
                           sortBY=sortBY)
                
                #make a list of genes 
                #whose adjusted expression p-value is less than specified
                #and lfc is greater than specified
                produce.heatmaps(data=data.norm.filt,
                                 controls=controls,
                                 treatments=treatments,
                                 fit2Eb=fit2Eb,
                                 contMatrix=contMatrix,
                                 newCont=newCont,
                                 LOGFC=LOGFC,
                                 PVALUE=PVALUE,
                                 forceLabelsON=forceLabelsON)
 
                #go out of current working directory
                setwd("../")
        }
}
        
