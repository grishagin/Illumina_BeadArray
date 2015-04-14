#make design and contrasts matrices for Illumina BeadArray

#making contrasts matrix
design.contrasts.Matrix.prepare<-function(data,controls,treatments){
        #making design matrix
        descr<-pData(data)$Sample.ID #define sample names
        #define sample names factor
        #need to do that as simple levels(descr) sorts levels in alphabetical order
        descrFac<-factor(descr,levels=unique(descr)) 
        designMatrix<-model.matrix(~0+descrFac) #make the design matrix
        colnames(designMatrix)<-unique(descr) #and edit its column names
        
        #declare list of contrasts...
        cont<-list()
        #and their counter...
        cont.counter<-1 
        #and a string with arguments to make contrasts matrix later
        mArg<-"makeContrasts("
        #list of counters of contrasts for each given control
        newCont<-list()
        
        for (c in 1:length(controls)){
                
                for (t in 1:length(treatments[[c]])){
                        
                        cont[cont.counter]<-paste0(unique(descr)[treatments[[c]][t]],"-",
                                                   unique(descr)[controls[c]])
                        
                        if (cont.counter==1) {
                                mArg<-paste0(mArg,unlist(cont)[cont.counter])
                        } else {
                                mArg<-paste(mArg,unlist(cont)[cont.counter],sep=",")
                        }
                        
                        newCont[[c]]<-c(unlist(newCont[c]),cont.counter)
                        
                        cont.counter<-cont.counter+1
                }
        }
        mArg<-paste0(mArg,",levels=designMatrix)")
        contMatrix<-eval(parse(text=mArg))
        desM.contM.newCont.list<-list(designMatrix,contMatrix,newCont)
        return(desM.contM.newCont.list)
}