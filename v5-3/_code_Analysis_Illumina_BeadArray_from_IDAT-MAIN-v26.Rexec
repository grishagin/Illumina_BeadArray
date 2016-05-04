######################################
#2016-05-04

#specify main folder for the code
main.code.folder<-getwd()
######################################
if (!require(tcltk)){
  message("Package 'tcltk' is missing! Close the window to abort.")
  Sys.sleep(600)
}

#source all functions for auxiliary files
dir<-tk_choose.dir(default = main.code.folder,
                caption = "Choose directory with auxiliary R scripts.")
files<-list.files(dir,
                  full.names = TRUE,
                  pattern=".R")
sapply(files,FUN=function(a) source(a))

#load all necessary packages
loadPackages()

#get a dataframe with all parameters
params<-read.params()

#######SET UP KEY PARAMETERS################################################
PVALUEvector<-params$p.value[!(is.na(params$p.value))]
LOGFCvector<-params$logFC[!(is.na(params$logFC))]
#numbers of columns in design matrix corresponding to controls
controls<-strsplit(as.character(params$Controls[!(is.na(params$Controls))]),
                   split="[[:punct:]]")
controls<-lapply(controls,as.numeric)
#numbers of columns in design matrix corresponding to treatments
treatments<-strsplit(as.character(params$Treatments[!(is.na(params$Treatments))]),
                     split="[[:punct:]]")
treatments<-lapply(treatments,as.numeric)

for (index in 1:length(controls)){
    if (length(controls[[index]])!=length(treatments[[index]])){
        tkmessageBox(message = "Number of Controls is not equal to number of Treatments!
Click 'OK' to abort and check your Parameters file."
                     ,icon = "error"
                     ,type = "ok")
        stop()
    }
}

#sort topTable by LFC or p-value?
#specify "logFC" or "B" or "p"
sortBY<-params$Sort.heatmap.by[!(is.na(params$Sort.heatmap.by))]
#force labels on heatmaps with more than 25 genes?
forceLabelsON<-params$Force.labels.ON[!(is.na(params$Force.labels.ON))]
forceLabelsON<-ifelse(forceLabelsON=="Yes" || forceLabelsON=="yes",TRUE,FALSE)
#determine if QC is required: Yes or No
performQC<-params$perform.QC[!(is.na(params$perform.QC))]
#determine if GEO submission files are desired or not
writeGEO<-params$make.GEO.files[!(is.na(params$make.GEO.files))]
writeGEO<-ifelse(writeGEO=="Yes" || writeGEO=="yes",TRUE,FALSE)
#perform batchCorrect?
batchCorrect<-params$Batch.Correct[!(is.na(params$Batch.Correct))]
batchCorrect<-ifelse(batchCorrect=="Yes" || batchCorrect=="yes",TRUE,FALSE)
############################################################################

dir.create("analysis", showWarnings = FALSE)
setwd("analysis")

#load and normalize the data form idat files
#arguments are optional
#defaults pertain to MouseRef-8 v2
data.norm.filt<-read.norm.filt.data(batchCorrect=batchCorrect,
                                    writeGEO=writeGEO)

if(performQC=="Yes" || performQC=="yes"){
        #produce QC on selected genes on filtered data
        QC(eset=data.norm.filt)
}

#make design matrix
designMatrix<-design.Matrix.prepare(data.norm.filt)

#go through all controls (and corresponding treatments)
#for vector element of controls list there must be
#one vector in treatments list
for (controls.counter in 1:length(controls)){   

        #make contrasts matrix
        #again, CONTROLS AND TREATMENTS ARE SET UP in THE BEGINNING of this module (above)!
        contrastsMatrix<-contrasts.Matrix.prepare(designMatrix=designMatrix,
                                                  controls=controls[[controls.counter]],
                                                  treatments=treatments[[controls.counter]])
        
        #fitting: linear regression, then Bayes
        fit2Eb<-fit.data(eset=data.norm.filt,
                         designMatrix=designMatrix,
                         contrastsMatrix=contrastsMatrix)
        
        #create a directory for each control
        #name it as "Treatment1-control1, treatment2-control2", etc.
        dirName<-paste0(colnames(contrastsMatrix),collapse=", ")
        dir.create(dirName, showWarnings = FALSE)
        setwd(dirName)
        
        #make volcano plots with no markings for each control
        produce.volcano.plot(fit2Eb=fit2Eb,
                             markSHAPE=FALSE,
                             markCOLOR=FALSE)

        #scan throu all p-values and log fold changes as specified above
        #make directory for each combination to store analyses in
        for (PVALUE in PVALUEvector){
                for(LOGFC in LOGFCvector){
                        #create a working directory for a given p-value and lfc
                        dirName<-paste0("p",PVALUE,"-lfc",LOGFC)
                        dir.create(dirName, showWarnings = FALSE)
                        setwd(dirName)
                        
                        results<-as.data.frame(decideTests(fit2Eb,
                                                           adjust.method="BH",
                                                           p.value = PVALUE,
                                                           lfc =LOGFC))
                        
                        #make volcano plots
                        produce.volcano.plot(fit2Eb=fit2Eb,
                                             PVALUE=PVALUE,
                                             LOGFC=LOGFC,
                                             markSHAPE=TRUE,
                                             markCOLOR=TRUE)
                        
                        #print venn diagrams for up-, down-, and both up- and down-regulated genes
                        #into respective files
                        vennD(treatments=treatments[[controls.counter]],
                              designMatrix=designMatrix,
                              results=results)
                        
                        #make the lists of top genes
                        top.tables(fit2Eb=fit2Eb,
                                   contrastsMatrix=contrastsMatrix,
                                   LOGFC=LOGFC,
                                   PVALUE=PVALUE,
                                   sortBY=sortBY)
                        
                        #make a list of genes 
                        #whose adjusted expression p-value is less than specified
                        #and lfc is greater than specified
                        produce.heatmaps(data=data.norm.filt,
                                         controls=controls[[controls.counter]],
                                         treatments=treatments[[controls.counter]],
                                         fit2Eb=fit2Eb,
                                         contrastsMatrix=contrastsMatrix,
                                         LOGFC=LOGFC,
                                         PVALUE=PVALUE,
                                         forceLabelsON=forceLabelsON)
         
                        #go out of current working directory
                        setwd("../")
                }
        }
        #go out of current working directory
        setwd("../")
}               
tkmessageBox(message = "All done!", icon = "info", type = "ok")

