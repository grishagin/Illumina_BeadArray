#Illumina BeadArray data analysis
#plot heatmap

plot.heatmap<-function(file.name,
                       data,
                       rowLabels,
                       Rowv,
                       dendrogram){
        
        pdf(file=file.name)
        par(cex.main=1)                                  #set main title font size
        heatmap.2.mod(exprs(data),             #dataset - genes expressed differentially
                      Rowv=Rowv,                            #no row clustering
                      dendrogram=dendrogram,                     #dendrogram - only for columns
                      col=greenred(300),                      #color scheme - green-black-red
                      scale="row",                           #scale by row
                      #ColSideColors = sideColors,            #color clustering of the samples
                      main="Differentially Expressed Genes",   #title of the heatmap
                      margins = c(10,10),                    #margins around heatmap    
                      key=TRUE,                              #legend is present
                      key.title="none",                      #no title for legend
                      key.xlab="",                           #no label for legend axis
                      keysize=1.2,                           #size of legend strip
                      symkey=FALSE,                          #do not make the colors symmetrical around 0
                      density.info="none",                   #no density information
                      trace="none",                          #
                      cexCol=1.2,                            #column labels' size
                      #cexRow=,                               #column labels' size
                      labRow=rowLabels,                      #row labels
                      labCol=pData(data)$Sample.ID)          #column labels
        dev.off()
}
