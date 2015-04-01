#Illumina BeadArray data analysis
#Looks for control (housekeeping, etc.) genes
#plots them and outputs into separate files (one file per gene)
#2015-03-18


QC<-function(data){
        
        geneNames<-list("ppia","cdkn1a","hprt")
        genes<-list()
        
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
                geneNames[[i]]<-paste0(temp,collapse="")
                genes[[i]]<-grep(geneNames[[i]],fData(data.norm.filt)$SYMBOL)
        }
        
        my.plots<-list()
        n=1
        
        for (g in genes){
                
                #declare data to plot as data frame
                toPlot.data.N.F<-data.frame(stringsAsFactors = FALSE)
                
                for (i in 1:length(g)){
                        Value<-exprs(data)[g[i],]
                        ProbeID<-as.character(fData(data)[g[i],"IlluminaID"])
                        Gene<-as.character(fData(data)$SYMBOL[g[i]])
                        SampleID<-names(Value)
                        
                        toPlot.data.N.F<-rbind(toPlot.data.N.F,
                                               data.frame(cbind(Value,ProbeID,Gene,SampleID),
                                                          stringsAsFactors = FALSE))
                }
                
                #get unique chip ID's
                uniqueSId<-paste(unique(str_extract(rownames(toPlot.data.N.F),"[0-9]+")),collapse=".")
                #compose a file name from chip IDs and gene name
                file.name<-paste("quality_control",uniqueSId,Gene,"jpg",sep=".")
                dir.create("./QC/", showWarnings = FALSE)
                
                
                gg<-ggplot(data=toPlot.data.N.F,
                           mapping = aes(x=SampleID,
                                         y=as.numeric(Value),
                                         color=as.factor(ProbeID),
                                         fill=as.factor(ProbeID)))+
                        ggtitle(paste0("Normalized Expression of ",Gene,"\n"))+
                        geom_point(size=5,shape=21)+
                        geom_text(aes(label=round(as.numeric(Value),digits=1),hjust=0, vjust=1))+
                        scale_y_continuous("Log(Expression Value)\n")+
                        scale_x_discrete("Sample Name")+
                        facet_wrap(facets = ~ProbeID)+
                        theme_bw(base_size = 20)+
                        theme(axis.text.x = element_text(angle = 90, hjust = 1))+
                        guides(fill=FALSE, color=FALSE)
                
                        #+scale_shape_discrete(name="Probe ID")+
                        #scale_fill_discrete(name="Probe ID")+
                        #scale_color_discrete(name="Probe ID")
                
                ggsave(file=file.name,
                       path="./QC/",
                       plot=gg,
                       width=length(g)*10,
                       height=10,
                       dpi=300,
                       units = "in")     
        } 
}