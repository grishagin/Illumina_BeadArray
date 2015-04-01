#Illumina BeadArray data analysis
#Prepare Top Tables

top.tables<-function(controls,
                     fit2Eb,
                     contMatrix,
                     newCont,
                     LOGFC,
                     PVALUE,
                     sortBY){
        for (c in 1:length(controls)){ 
                #which columns from contrasts matrix to look at
                annotContMatrix<-colnames(contMatrix)[newCont[[c]]]
                for (i in 1:length(newCont[[c]])){
                        # store table in temp var
                        temp <- topTable(fit2Eb,
                                         coef=newCont[[c]][i],
                                         number=Inf,
                                         lfc=LOGFC,
                                         p.value=PVALUE,
                                         sort.by = sortBY)
                        # if there are no genes in top table,
                        # add one NA value to temp var to avoid error
                        if (length(temp)==0) temp[1,1]<-NA
                        print(length(temp$ID.SYMBOL))
                        
                        # write top table into a variable
                        assign(paste0("topExpr.",annotContMatrix[i]),
                               temp) 
                        
                        # write top table into a file 
                        write.xlsx(temp,
                                   paste0("topExpr.",annotContMatrix[i],".xlsx"))
                        
                        # append all of top table's first N genes into a file
                        if(i==1) {
                                allTopTables<-loadWorkbook(paste0("topExpr.ALLsample-",c,"-",
                                                                  colnames(designMatrix)[controls[c]],".xlsx"),
                                                           create=TRUE)
                                createSheet(allTopTables,name="allTopTables")
                        } 
                        writeWorksheet(allTopTables,paste0("topExpr.",annotContMatrix[i]),
                                       sheet="allTopTables",
                                       startRow=(8*(i-1)+1),
                                       header=FALSE)
                        writeWorksheet(allTopTables,temp[1:5,],
                                       sheet="allTopTables",
                                       startRow=(8*(i-1)+2))
                }     
                saveWorkbook(allTopTables)
        }
}