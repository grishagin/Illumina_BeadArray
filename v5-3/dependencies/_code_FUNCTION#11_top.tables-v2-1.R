#Illumina BeadArray data analysis
#Prepare Top Tables

top.tables<-function(fit2Eb,
                     contrastsMatrix,
                     LOGFC,
                     PVALUE,
                     sortBY){
        
        #which columns from contrasts matrix to look at
        contrastsMatrix.colnames<-colnames(contrastsMatrix)
        for (i in 1:length(contrastsMatrix.colnames)){
                # store table in temp var
                temp <- topTable(fit2Eb,
                                 coef=i,
                                 number=Inf,
                                 lfc=LOGFC,
                                 p.value=PVALUE,
                                 sort.by = sortBY)
                # if there are no genes in top table,
                # add one NA value to temp var to avoid error
                if (length(temp)==0) temp[1,1]<-NA
                #print(length(temp$ID.SYMBOL))
                
                # write top table into a variable
                assign(paste0("topExpr.",contrastsMatrix.colnames[i]),
                       temp) 
                
                # write top table into a file 
                write.xlsx(temp,
                           paste0("topExpr.",contrastsMatrix.colnames[i],".xlsx"))
                
                # append all of top table's first N genes into a file
                if(i==1) {
                        allTopTables<-loadWorkbook(paste0("topExpr_allContr_smallSample.xlsx"),
                                                   create=TRUE)
                        createSheet(allTopTables,name="allTopTables")
                } 
                writeWorksheet(allTopTables,paste0("topExpr.",contrastsMatrix.colnames[i]),
                               sheet="allTopTables",
                               startRow=(8*(i-1)+1),
                               header=FALSE)
                writeWorksheet(allTopTables,temp[1:5,],
                               sheet="allTopTables",
                               startRow=(8*(i-1)+2))
        }     
        saveWorkbook(allTopTables)
}
                