#Illumina BeadArray data
#source all necessary auxiliaries

source.all.functions<-function(){
        source(paste0(main.code.folder,"\\_code_FUNCTION#1__loadPackages-v1.R"))
        
        source(paste0(main.code.folder,"\\_code_FUNCTION#2_read.norm.filt.data-v2.R"))
        source(paste0(main.code.folder,"\\_code_FUNCTION#3_QC.CDKN1A.PPIA.HPRT-v4.R"))
        
        source(paste0(main.code.folder,"\\_code_FUNCTION#4_designMatrix.contMatrix.prep-v1.R"))
        source(paste0(main.code.folder,"\\_code_FUNCTION#5_fit.data-v1.R"))
        
        source(paste0(main.code.folder,"\\_code_FUNCTION#6_volcano.plot-v1.R"))
        
        source(paste0(main.code.folder,"\\_code_FUNCTION#7_venn-v1.R"))
        source(paste0(main.code.folder,"\\_code_FUNCTION#8_vennDiagram.mod-v1.R"))
        
        source(paste0(main.code.folder,"\\_code_FUNCTION#9_top.tables-v1.R"))
        
        source(paste0(main.code.folder,"\\_code_FUNCTION#10_produce.heatmaps-v1.R")) 
        source(paste0(main.code.folder,"\\_code_FUNCTION#11_prep.heatmap.sort.data-v1.R"))
        source(paste0(main.code.folder,"\\_code_FUNCTION#12_plot.heatmap-v1-1.R"))
        source(paste0(main.code.folder,"\\_code_FUNCTION#13_heatmap.2.mod-v1-1.R")) 
}
