#make contrasts matrix for Illumina BeadArray

contrasts.Matrix.prepare<-function(designMatrix,
                                   controls,
                                   treatments){
        #design matrix columns/sample names
        sampleNames<-colnames(designMatrix)
        #declare list of contrasts...
        contrasts<-NULL
        #and their counter...
        contrasts.counter<-1 
        #and a string with arguments to make contrasts matrix later
        mArg<-"makeContrasts("
        
        for (t in 1:length(treatments)){
                #make expressions "Treatment-Control" for each treatment
                contrasts[t]<-paste0(sampleNames[treatments[t]],"-",
                                                                sampleNames[controls[t]])
        }
        #make expression to prepare contrasts matrix       
        mArg<-paste0(mArg,paste(contrasts,collapse=","),",levels=designMatrix)")
        #make contrasts matrix and return as function output
        contrastsMatrix<-eval(parse(text=mArg))
        return(contrastsMatrix)
}