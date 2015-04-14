######################################
#2015-04-14
#fully tested
#WORKS
######################################
#######SET UP KEY PARAMETERS##########
        #analysis p-value and log fold change
PVALUEvector<-c(0.001,0.01,0.05,1)
LOGFCvector<-c(0.5,1,2)
        #numbers of columns in design matrix corresponding to controls
controls<-c(1,1) 
        #numbers of columns in design matrix corresponding to treatments
treatments<-list(c(2,3),4)
        #sort topTable by LFC or p-value?
        #specify "logFC" or "B" or "p"
sortBY<-"p" 
        #force labels on heatmaps with more than 25 genes?
forceLabelsON=TRUE

        #specify main folder for the code
main.code.folder<-getwd()
######################################

#source all functions for auxiliary files
dir<-choose.dir(default = main.code.folder,caption = "Choose directory with auxiliary R scripts.")
files<-list.files(dir,full.names = TRUE,pattern=".R")
sapply(files,FUN=function(a) source(a))

#load all necessary packages
loadPackages()

#load and normalize the data form idat files
#arguments are optional
#defaults pertain to MouseRef-8 v2
data.norm.filt<-read.norm.filt.data(writeGEO=FALSE,softTemplate=NULL)

#produce QC on selected genes on filtered data
QC(eset=data.norm.filt)

#make design matrix
designMatrix<-design.Matrix.prepare(data.norm.filt)

#go through all controls (and corresponding treatments)
#for each element of controls vector there must be
#one vector in treatments list
for (controls.counter in 1:length(controls)){   

        #make contrasts matrix
        #again, CONTROLS AND TREATMENTS ARE SET UP in THE BEGINNING of this module (above)!
        contrastsMatrix<-contrasts.Matrix.prepare(designMatrix=designMatrix,
                                                  control=controls[controls.counter],
                                                  treatments=treatments[[controls.counter]])

        #fitting: linear regression, then Bayes
        fit2Eb<-fit.data(eset=data.norm.filt,
                         designMatrix=designMatrix,
                         contrastsMatrix=contrastsMatrix)
        
        #create a directory for each control
        #name it as "Treatment1-control, treatment2-control", etc.
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
                        vennD(control=controls[controls.counter],
                              treatments=treatments[[controls.counter]],
                              designMatrix=designMatrix,
                              results=results)
                        
                        #make the lists of top genes
                        top.tables(control=controls[controls.counter],
                                   fit2Eb=fit2Eb,
                                   contrastsMatrix=contrastsMatrix,
                                   LOGFC=LOGFC,
                                   PVALUE=PVALUE,
                                   sortBY=sortBY)
                        
                        #make a list of genes 
                        #whose adjusted expression p-value is less than specified
                        #and lfc is greater than specified
                        produce.heatmaps(data=data.norm.filt,
                                         control=controls[controls.counter],
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
