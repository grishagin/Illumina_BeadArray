#Illumina BeadArray data analysis
#print venn diagrams for up-, down-, and both up- and down-regulated genes
#into respective files

vennD<-function(controls,
                treatments,
                designMatrix,
                results){
        for (c in 1:length(controls)){ 
                
                for (i in 1:3){
                        subType<-ifelse(i==1, "UP",ifelse(i==2, "DOWN","BOTH"))
                        
                        file.name=paste0("Venn Diagrams.-control#",c,"-",
                                         colnames(designMatrix)[controls[c]],
                                         "-",subType,".pdf")
                        pdf(file=file.name)
                        
                        vennDiagram.mod(results[,newCont[[c]]],
                                        include=subType,
                                        names=colnames(designMatrix)[treatments[[c]]],
                                        show.include=TRUE,
                                        cex=1.4,
                                        cexCounts=2)
                        dev.off() 
                        
                }
                
        }
}

