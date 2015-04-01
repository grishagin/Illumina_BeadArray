#function to make volcano plot for BeadArray

make.volcano<-function(file.name,df,x,y){
        
        if (y=="tempLODS") {
                y.label="Ln(Odds)\n"
        } else if (y=="tempNegLogAdjPval") {
                y.label = expression(atop(-log[10]*"(Adjusted P-value)"," "))
        } else {
                simpleError(message = "Wrong type!")
                stop()
        }
 
        g<-ggplot(data=df, 
                  aes_string(x=x, y=y)) +
                geom_point(size=3)+
                theme_bw(base_size = 20)+
                theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      panel.border=element_blank(),
                      axis.line = element_line(colour = "black"))+
                guides(fill=FALSE, color=FALSE,shape=FALSE)+
                xlab(label = expression(atop(" ",Log[2]*"(Fold Change)")))+
                ylab(label = y.label)+
                scale_shape_manual(values=shapes)
        
        ggsave(filename = file.name,
               plot = g,
               width = 10,
               height = 10,
               units = "in",
               dpi=300)
}
