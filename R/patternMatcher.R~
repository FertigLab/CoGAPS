#' PatternMatcher Shiny Ap
#'
#' @param PBySet list of matched set solutions for the Pmatrix from an NMF algorithm
#' @param out optional name for saving output
#' @param order optional vector indicating order of samples for plotting. Default is NULL.
#' @param sample.color optional vector of colors of same lenght as colnames. Default is NULL.
#'
#' @return either an index of selected sets' contributions or the editted \code{PBySet} object
#' @export
#'
#' @examples \dontrun{
#' patternMatcher(PBySet,out,order,sample.color)
#' }
#'
#'
patternMatcher<-function(PBySet=NULL,out=NULL,order=NULL, sample.color=NULL) {

runApp(list(
  ui = pageWithSidebar(
    # Application title
    headerPanel('NMF Pattern Matching'),
    # Side pannel with controls
    sidebarPanel(
      # to upload file
      fileInput('file1',
                'Choose .Rda File',
                accept=c('.RData','.Rda','R data object','.rda')
      ),
      #
      uiOutput("pickplot"),
      uiOutput("checkbs"),
      downloadButton('downloadData', 'Download'),
      actionButton("end", "Return")
    ),
    # Main panel with plots
    mainPanel(
      plotOutput('plot1')
    )
  ),

  server = function(input, output, session) {
    #load in the data
    df<-reactive({
      if(!is.null(PBySet)){
        df<-PBySet
        return(df)
      }
      inFile <- input$file1 # get the path to the input file on the server
      if (is.null(inFile)){return(NULL)}
      load(inFile$datapath) #load it
      df <- get(load(inFile$datapath))# get the name of the object that was loaded
      return(df)
    })

    # get data name
    datName<-reactive({
      if(!is.null(out)){
        datName<-paste(out,'.SelectedPatterns.Rda',sep="")
        return(datName)
      }
      inFile <- input$file1
      if (is.null(inFile) & is.null(out)){
        datName<-"SelectedPatterns.Rda"
        return(datName)
      }
      if (is.null(inFile)){return(NULL)}
      fn<-strsplit(inFile$name,"[.]")[[1]][1]
      datName<-paste(fn,'.SelectedPatterns.Rda',sep="")
      return(datName)
    })


    mdf=reactive({# use that to give options for subsetting, some formatting may need to be removed
      dfx=df()
      if (is.null(dfx)){return(NULL)}
      mdf=melt(dfx,stringsAsFactors=FALSE) # melt the elements of the list
      colnames(mdf)<-c("BySet","Samples","value","Patterns")
      mdf$BySet<-as.character(mdf$BySet) # change them to characters
      mdf$Samples<-as.character(mdf$Samples)
      mdf$value=as.numeric(mdf$value) #make sure value is numeric for plotting
      str(mdf)
      return(mdf)
    })


    output$pickplot <- renderUI({# menu to select which matrix to look at/edit
      if (is.null(df())){return(NULL)}
      mdf2=mdf()
      selectInput("whichplot", "Select the Pattern to Plot",choices=unique(mdf2$Patterns))
    })


    output$checkbs <- renderUI({# make the checkboxes for each row of each matrix
      if (is.null(df())){return(NULL)}
      mdf2=mdf()
      lapply(unique(mdf2$Patterns), function(i) {
        subss <- unique(mdf2$BySet[mdf2$Patterns==i]) # find the rows (after it has been melted)
        tmp=sprintf("input.whichplot ==  %g", i) # create the javascript code to make this a conditional panel
        conditionalPanel(
          condition = tmp,
          checkboxGroupInput(paste("subs",i,sep=""), i, choices=subss, selected=subss) # the actual checkboxes for each, subs1, subs2, subsn
        )
      })
    })


    output$plot1 <- renderPlot({#plot the data, subset to the desired columns
      # if there has not been an uploaded matrix yet, don't even try to make a plot
      if (is.null(df())){return(NULL)}
      if (is.null(input$whichplot)){return(NULL)}
      par(mar = c(5.1, 4.1, 0, 1))
      mdf2=mdf() # grab the melted data frame to use the ggplot2 plot
      x=input$whichplot # which matrix to show
      x=as.numeric(x)
      tmp=input[[paste("subs",x,sep="")]] # get the rows that have been selected
      mdf2x=mdf2[which(mdf2$BySet%in%tmp),]
      if (!is.null(order) & !is.null(sample.color)){
       ggplot(mdf2x, aes(x=Samples, y=value, col=BySet,group=BySet))+
        geom_line() + scale_x_discrete(limits=order) +
        theme(axis.text.x = element_text(angle=45,family="Helvetica-Narrow", hjust = 1,colour = sample.color))
      } else if(!is.null(sample.color) & is.null(order)) {
      ggplot(mdf2x, aes(x=Samples, y=value, col=BySet,group=BySet))+
          geom_line() +
          theme(axis.text.x = element_text(angle=45,family="Helvetica-Narrow", hjust = 1,colour = sample.color))
      } else if(!is.null(order) & is.null(sample.color) ) {
        ggplot(mdf2x, aes(x=Samples, y=value, col=BySet,group=BySet))+
          geom_line() + scale_x_discrete(limits=order) +
          theme(axis.text.x = element_text(angle=45,family="Helvetica-Narrow", hjust = 1))
      } else {
        ggplot(mdf2x, aes(x=Samples, y=value, col=BySet,group=BySet))+
          geom_line() +
          theme(axis.text.x = element_text(angle=45,family="Helvetica-Narrow", hjust = 1))
      }
      #pplot
      #browser()
    })

    # create and download the final result file
    output$downloadData <- downloadHandler(
      filename = datName(), # set the file name
      content = function(file) {
        PatternsSelect <- lapply(1:length(mdf()), function(i) {input[[paste("subs",i,sep="")]]})
        save(PatternsSelect, file=file) # generate the object to save
      }
    )
    #stop app and return to R
    observeEvent(input$end, {
      mdf2=mdf()
      PatternsSelect <- sapply(1:length(df()), function(i) {input[[paste("subs",i,sep="")]]})
      selectPBySet <- mdf2[which(mdf2$BySet%in%PatternsSelect),]
      stopApp(returnValue = selectPBySet)
    })



  }

)
)

}
