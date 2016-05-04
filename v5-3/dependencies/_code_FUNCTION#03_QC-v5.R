#Illumina BeadArray data analysis
#Looks for control (housekeeping, etc.) genes
#plots them and outputs into separate files (one file per gene)
#2015-04-14


QC<-function(eset){
        
        geneNames<-list("ppia","cdkn1a","hprt","hox")
        genes<-list()
        
        #create directory for quality control plots
        dir.create("./QC/", showWarnings = FALSE)
        
        #get unique chip ID's
        uniqueSId<-paste(unique(str_extract(colnames(eset),"[0-9]+")),collapse=".")
        
        #first, use plotMA to look for mean-average trends
        pdf(file = paste("./QC/quality_control.MAplots",uniqueSId,"pdf",sep="."),)
        for (i in 1:ncol(eset)) {
                limma::plotMA(exprs(eset),array=i)
        }
        dev.off()
        
        #go through the list of gene names and split every name by letter
        #then turn every letter x (or X) into [xX]; leave every number unaltered
        #for grep search; perform grep search, output into genes list
        for (i in 1:length(geneNames)){
                temp<-unlist(strsplit(geneNames[[i]],split=""))
                for(j in 1:length(temp)){
                        if (temp[j] %in% letters) {
                                num<-which(letters %in% temp[j])
                                temp[j]<-paste0("[",letters[num],LETTERS[num],"]")
                        } else if (temp[j] %in% LETTERS) {
                                num<-which(LETTERS %in% temp[j])
                                temp[j]<-paste0("[",letters[num],LETTERS[num],"]")
                        }

                }
                #make a regexp search query
                geneNames[[i]]<-paste0("^",paste0(temp,collapse=""))
                #find all unique gene symbols - they will make up new search query
                genes[[i]]<-unique(grep(geneNames[[i]],fData(eset)$SYMBOL,value=TRUE))                
        }
        #turn the list of genes into a vector
        genes<-unlist(genes)
        #search the eset for the genes again, this time using proper list of unique genes
        genes<-lapply(genes,FUN=function(a) grep(a,fData(eset)$SYMBOL))

        for (g in genes){
                
                #declare data to plot as data frame
                toPlot.data.N.F<-data.frame(stringsAsFactors = FALSE)
                
                for (i in 1:length(g)){
                        Value<-exprs(eset)[g[i],]
                        ProbeID<-as.character(fData(eset)[g[i],"IlluminaID"])
                        Gene<-as.character(fData(eset)$SYMBOL[g[i]])
                        SampleID<-names(Value)
                        
                        toPlot.data.N.F<-rbind(toPlot.data.N.F,
                                               data.frame(cbind(Value,ProbeID,Gene,SampleID),
                                                          stringsAsFactors = FALSE))
                }
                
                #compose a file name from chip IDs and gene name
                file.name<-paste("quality_control",uniqueSId,Gene,"jpg",sep=".")
                
                gg<-ggplot(data=toPlot.data.N.F,
                           mapping = aes(x=SampleID,
                                         y=as.numeric(Value),
                                         color=as.factor(ProbeID),
                                         fill=as.factor(ProbeID)))+
                        ggtitle(paste0("Normalized Expression of ",Gene,"\n"))+
                        geom_point(size=4,shape=21)+
                        geom_text(aes(label=round(as.numeric(Value),digits=1),hjust=0, vjust=1))+
                        scale_y_continuous("Log(Expression Value)\n")+
                        scale_x_discrete("Sample Name")+
                        facet_wrap(facets = ~ProbeID,nrow=1)+
                        theme_bw(base_size = 20)+
                        theme(axis.text.x = element_text(angle = 90, hjust = 1))+
                        guides(fill=FALSE, color=FALSE)
                
                        #+scale_shape_discrete(name="Probe ID")+
                        #scale_fill_discrete(name="Probe ID")+
                        #scale_color_discrete(name="Probe ID")
                
                ggsave(file=file.name,
                       path="./QC/",
                       plot=gg,
                       width=length(g)*12,
                       height=12,
                       dpi=300,
                       units = "in",
                       limitsize=FALSE)     
        } 
}