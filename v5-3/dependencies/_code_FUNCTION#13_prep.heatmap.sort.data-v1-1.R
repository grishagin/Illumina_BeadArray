#Illumina BeadArray data analysis
#prepare gene lists for heatmaps sorted by PVAL and LFC
#and appropriate labels

prep.heatmap.sort.data<-function(data.sel,   #data.norm.filt.sel
                                 selGeneAdjPval,
                                 selGeneLogFC,
                                 forceLabelsON){
  
        #remove duplicate probe entries from the lists
        #of adjusted p-values and lfc's
        #by leaving only one with the smallest p-values and largest LFC's
        #which are amongst overall selected genes
        uniqueNamesPVAL<-unique(names(selGeneAdjPval))
        selGeneAdjPval<-sapply(uniqueNamesPVAL,
                               FUN=function(a) min(selGeneAdjPval[names(selGeneAdjPval)==a]))
        selGeneAdjPval<-selGeneAdjPval[rownames(data.sel)]
        
        uniqueNamesLogFC<-unique(names(selGeneLogFC))
        selGeneLogFC<-sapply(uniqueNamesLogFC,
                             FUN=function(a) max(selGeneLogFC[names(selGeneLogFC)==a]))
        selGeneLogFC<-selGeneLogFC[rownames(data.sel)]
        
        #if there's less than two genes - abort
        if (nrow(data.sel)<2) return("Error: cannot plot heatmap.")
        
        #ORDER selected genes by adjusted p-value
        data.for.heatmap.PVAL<-data.sel
        fData(data.for.heatmap.PVAL)<-fData(data.for.heatmap.PVAL)[order(selGeneAdjPval),]
        exprs(data.for.heatmap.PVAL)<-exprs(data.for.heatmap.PVAL)[order(selGeneAdjPval),]
        
        #ORDER selected genes by lfc
        data.for.heatmap.LFC<-data.sel
        fData(data.for.heatmap.LFC)<-fData(data.for.heatmap.LFC)[order(abs(selGeneLogFC),
                                                                       decreasing = TRUE),]
        exprs(data.for.heatmap.LFC)<-exprs(data.for.heatmap.LFC)[order(abs(selGeneLogFC),
                                                                       decreasing = TRUE),]
        
        #determine if gene labels are desired
        if ((length(fData(data.sel)$SYMBOL)>25) & (forceLabelsON==FALSE)) {
                geneLabels<-FALSE
        } else {
                geneLabels<-list()
                geneLabels[[1]]<-as.character(fData(data.sel)$SYMBOL)
                geneLabels[[2]]<-as.character(fData(data.for.heatmap.PVAL)$SYMBOL)
                geneLabels[[3]]<-as.character(fData(data.for.heatmap.LFC)$SYMBOL)
        }
        dataList.for.heatmap<-list(data.sel,
                                   data.for.heatmap.PVAL,
                                   data.for.heatmap.LFC,
                                   geneLabels)
        return(dataList.for.heatmap)
}
