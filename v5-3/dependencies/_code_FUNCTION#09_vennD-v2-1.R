#Illumina BeadArray data analysis
#print venn diagrams for up-, down-, and both up- and down-regulated genes
#into respective files

vennD<-function(treatments,
                designMatrix,
                results){
        
        for (i in 1:3){
                subType<-ifelse(i==1, "UP",ifelse(i==2, "DOWN","BOTH"))
                
                file.name=paste0("Venn.Diagram-",
                                 subType,
                                 ".pdf")
                pdf(file=file.name)
                
                vennDiagram.mod(results,
                                include=subType,
                                names=colnames(designMatrix)[treatments],
                                show.include=TRUE,
                                cex=1.4,
                                cexCounts=2)
                dev.off() 
                
        }
}

