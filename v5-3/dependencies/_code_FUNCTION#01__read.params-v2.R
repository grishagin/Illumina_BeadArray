read.params<-function(){
        #read a file with parameters
        params<-xlsx:::read.xlsx(file = tk_choose.files(caption="Select a file with parameters."),
                                 sheetIndex = 1,
                                 header = TRUE,
                                 stringsAsFactors=FALSE)
        return(params)
}