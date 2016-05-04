#function to make volcano plot for BeadArray

make.volcano<-function(file.name,
                       dFrame,
                       x,
                       y,
                       shapes=shapes,
                       clrs=clrs,
                       markSHAPE=FALSE,
                       markCOLOR=FALSE){
        
        if (y=="tempLODS") {
                y.label="Ln(Odds)\n"
        } else if (y=="tempNegLogAdjPval") {
                y.label = expression(atop(-log[10]*"(Adjusted P-value)"," "))
        } else {
                simpleError(message = "Wrong type!")
                stop()
        }
        
        
        if (markSHAPE & markCOLOR){
                g<-ggplot(data=dFrame, 
                          aes_string(x=x, y=y),environment=environment()) +
                        geom_point(size=3,
                                   aes(shape=factor(dFrame$flag),
                                       color=factor(dFrame$flag)))
        } else if (markSHAPE) {
                g<-ggplot(data=dFrame, 
                          aes_string(x=x, y=y)) +
                        geom_point(size=3,
                                   aes(shape=factor(dFrame$flag)))
        } else if(markCOLOR) {
                g<-ggplot(data=dFrame, 
                          aes_string(x=x, y=y),environment=environment()) +
                        geom_point(size=3,
                                   aes(color=factor(dFrame$flag)))
        } else {
                g<-ggplot(data=dFrame, 
                          aes_string(x=x, y=y),environment=environment()) +
                        geom_point(size=3)
        }
        
        g<-g+theme_bw(base_size = 20)+
                theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      panel.border=element_blank(),
                      axis.line = element_line(colour = "black"))+
                guides(fill=FALSE, color=FALSE,shape=FALSE)+
                xlab(label = expression(atop(" ",Log[2]*"(Fold Change)")))+
                ylab(label = y.label)+
                scale_shape_manual(values=shapes)+
                scale_color_manual(values=clrs)
        
        ggsave(filename = file.name,
               plot = g,
               width = 10,
               height = 10,
               units = "in",
               dpi=300)
}
