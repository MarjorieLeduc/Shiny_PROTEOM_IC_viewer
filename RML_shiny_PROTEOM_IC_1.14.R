library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(ggplot2)
library(plotly)
library(magrittr)
library(purrr)
library(RColorBrewer)
library(imputeLCMD)
library(ggalt)
library(ggrepel)
library(ComplexHeatmap)
library(svglite)
library(ragg)
library(stringr)

####UI----
options(shiny.maxRequestSize=100*1024^2)

ui <- dashboardPage(
  dashboardHeader(title = "Proteom'IC"),
  dashboardSidebar( width = 130,
                    sidebarMenu(id = "tabs",
                                menuItem("Import data", tabName = "ImportData"),
                                menuItem("Descriptives Plots", tabName = "DescPlot"),
                                menuItem("Statistical Plots", tabName = "StatPlot"),
                                menuItem("Proteins Plots",tabName="ProtPlot")
                    )
  ),
  dashboardBody(
    tabItems(
      ##### Read data ----
      tabItem(tabName = "ImportData",
              h1("Import Data"),
              fileInput("dataFile",label = NULL,
                        buttonLabel = "Browse...",
                        placeholder = "No file selected"),
              textOutput("TextNbProt"),
              br(),
              dataTableOutput("value")
              
      ),
      #####Descriptiveplot ----
      tabItem(tabName = "DescPlot",
              uiOutput("colors"),
              tabBox( width = 0,
                      ###### NbProt ----
                      tabPanel("Nbprot",
                               plotlyOutput('nbprot',width = "100%", height = "400"),
                               selectInput(
                                 inputId = "NbProtformat",
                                 label = "Format of downloaded plot",
                                 choices = c("svg", "png", "pdf"),
                                 selected = "svg"
                               ),
                               textInput("WidthNbProt", "Width of downloaded plot (cm)", "20"),
                               textInput("HeightNbProt", "Height of downloaded plot (cm)", "15"),
                               helpText("Plot is downloaded at 600dpi."),
                               downloadButton("downloadNbProtplot", "Download")
                               ),
                      ###### PCA ----
                      tabPanel("PCA",  
                               selectInput(
                                 inputId = "TypeDataPCA",
                                 label = "Type of data used to generate PCA",
                                 choices = c("Z-score(intensity)", "Log2(intensity)"),
                                 selected = "Z-score(intensity)"
                               ),
                               selectInput(
                                 inputId = "FilterType",
                                 label = "Type of filter used to generate PCA",
                                 choices = c("% of Valid value in at least one group", "% of Valid value in each group","% of Valid value in total"),
                                 selected = "% of Valid value in at least on group"
                               ),
                               textInput("PcFilter","% of valid values used to generate PCA","70"),
                               actionButton("DrawPCA","Draw PCA"),
                               helpText("WARRNING : PCA take a long time to be displayed."),
                               textOutput("NbProtPCA"),
                               plotOutput('PCA',width = "100%", height = "400"),
                               selectInput(
                                 inputId = "PCAformat",
                                 label = "Format of downloaded plot",
                                 choices = c("svg", "png", "pdf"),
                                 selected = "svg"
                               ),
                               textInput("WidthPCA", "Width of downloaded plot (cm)", "26"),
                               textInput("HeightPCA", "Height of downloaded plot (cm)", "17"),
                               helpText("Plot is downloaded at 600dpi."),
                               downloadButton("downloadPCAplot", "Download"),
                               
                               plotlyOutput('PCAprotInteractive',width = "100%", height = "400"),
                               selectizeInput(
                                 inputId = "GenePCA",
                                 label = "Select which genes to show",
                                 choices = NULL,
                                 multiple = TRUE
                               ),
                               textInput("ForcePCA", "Force to avoid overlap of gene names", "1"),
                               plotOutput('PCAprot',width = "100%", height = "400"),
                               selectInput(
                                 inputId = "PCAprotformat",
                                 label = "Format of downloaded plot",
                                 choices = c("svg", "png", "pdf"),
                                 selected = "svg"
                               ),
                               textInput("WidthPCAprot", "Width of downloaded plot (cm)", "26"),
                               textInput("HeightPCAprot", "Height of downloaded plot (cm)", "17"),
                               helpText("Plot is downloaded at 600dpi."),
                               downloadButton("downloadPCAprotplot", "Download")
                               ),
                      
                      ###### Pearson Correlations ----
                      tabPanel("Pearson correlations", 
                               checkboxInput(
                                 inputId = "PearsonCluster",
                                 label = "Add Sample clustering",
                                 value = FALSE
                               ),
                               plotOutput('Pearson',width = "100%", height = "400"),
                               selectInput(
                                   inputId = "Pearsonformat",
                                   label = "Format of downloaded plot",
                                   choices = c("svg", "png", "pdf"),
                                   selected = "svg"
                                 ),
                               textInput("WidthPearson", "Width of downloaded plot (cm)", "20"),
                               textInput("HeightPearson", "Height of downloaded plot (cm)", "15"),
                               helpText("Plot is downloaded at 600dpi."),
                               downloadButton("downloadPearsonplot", "Download"),
                               br(),
                               plotOutput('PearsonBox',width = "100%", height = "400"),
                               selectInput(
                                 inputId = "PearsonBoxformat",
                                 label = "Format of downloaded plot",
                                 choices = c("svg", "png", "pdf"),
                                 selected = "svg"
                               ),
                               textInput("WidthPearsonBox", "Width of downloaded plot (cm)", "20"),
                               textInput("HeightPearsonBox", "Height of downloaded plot (cm)", "15"),
                               helpText("Plot is downloaded at 600dpi."),
                               downloadButton("downloadPearsonBoxplot", "Download")
                               ),
                      
                      ###### Euclidean distances ----
                      tabPanel("Euclidean distances", 
                               checkboxInput(
                                 inputId = "EuclidCluster",
                                 label = "Add Sample clustering",
                                 value = FALSE
                               ),
                               plotOutput('Euclid',width = "100%", height = "400"),
                               selectInput(
                                 inputId = "Euclidformat",
                                 label = "Format of downloaded plot",
                                 choices = c("svg", "png", "pdf"),
                                 selected = "svg"
                               ),
                               textInput("WidthEuclid", "Width of downloaded plot (cm)", "20"),
                               textInput("HeightEuclid", "Height of downloaded plot (cm)", "15"),
                               helpText("Plot is downloaded at 600dpi."),
                               downloadButton("downloadEuclidplot", "Download"),
                               br(),
                               plotOutput('EuclidBox',width = "100%", height = "400"),
                               selectInput(
                                 inputId = "EuclidBoxformat",
                                 label = "Format of downloaded plot",
                                 choices = c("svg", "png", "pdf"),
                                 selected = "svg"
                               ),
                               textInput("WidthEuclidBox", "Width of downloaded plot (cm)", "20"),
                               textInput("HeightEuclidBox", "Height of downloaded plot (cm)", "15"),
                               helpText("Plot is downloaded at 600dpi."),
                               downloadButton("downloadEuclidBoxplot", "Download")
                               ),
                      ###### Heatmap Log2(intensity) ----
                      tabPanel("Heatmap Log2(intensity)", 
                               colourpicker::colourInput(
                                 inputId = "LowColor",
                                 label = "Choose the color of low expression proteins",
                                 value = "#00FF00"
                               ),
                               colourpicker::colourInput(
                                 inputId = "MiddleColor",
                                 label = "Choose the color of middle expression proteins",
                                 value = "#000000"
                               ),
                               colourpicker::colourInput(
                                 inputId = "HighColor",
                                 label = "Choose the color of high expression proteins",
                                 value = "#FF0000"
                               ),
                               checkboxInput(
                                 inputId = "SampleCluster",
                                 label = "Add Sample clustering",
                                 value = FALSE
                               ),
                               plotOutput('HeatmapLog2',width = "100%", height = "400"),
                               selectInput(
                                 inputId = "HeatmapLog2format",
                                 label = "Format of downloaded plot",
                                 choices = c("svg", "png", "pdf"),
                                 selected = "svg"
                               ),
                               textInput("WidthHeatmapLog2", "Width of downloaded plot (cm)", "15"),
                               textInput("HeightHeatmapLog2", "Height of downloaded plot (cm)", "15"),
                               helpText("Plot is downloaded at 600dpi."),
                               downloadButton("downloadHeatmapLog2plot", "Download")
                      ),
                      
                      ###### Scaterplot ----
                      tabPanel("Scater plot", 
                               selectInput(
                                 inputId = "Sample1",
                                 label = "Sample on X axis",
                                 choices = NULL,
                                 selected = NULL
                               ),
                               selectInput(
                                 inputId = "Sample2",
                                 label = "Sample on Y axis",
                                 choices = NULL,
                                 selected = NULL
                               ),
                               textOutput('CorScaterplot'),
                               actionButton("DrawScaterPlot","Draw Scaterplot"),
                               plotlyOutput('Scaterplotinteractive',width = "100%", height = "400"),
                               selectizeInput(
                                 inputId = "GenesScaterplot",
                                 label = "Select genes to show",
                                 choices = NULL,
                                 multiple = TRUE
                               ),
                               textInput("ForceScaterplot", "Force to avoid overlap of gene names", "1"),
                               plotOutput('Scaterplot',width = "100%", height = "400"),
                               selectInput(
                                 inputId = "ScaterPlotformat",
                                 label = "Format of downloaded plot",
                                 choices = c("svg", "png", "pdf"),
                                 selected = "svg"
                               ),
                               textInput("WidthScaterPlot", "Width of downloaded plot (cm)", "15"),
                               textInput("HeightScaterPlot", "Height of downloaded plot (cm)", "15"),
                               helpText("Plot is downloaded at 600dpi."),
                               downloadButton("downloadScaterPlot", "Download")
                      ),
                      ###### MeanScaterplot ----
                      tabPanel("Mean Scater plot", 
                               selectInput(
                                 inputId = "Cond1",
                                 label = "Condition on X axis",
                                 choices = NULL,
                                 selected = NULL
                               ),
                               selectInput(
                                 inputId = "Cond2",
                                 label = "Condition on y axis",
                                 choices = NULL,
                                 selected = NULL
                               ),
                               textOutput('CorMeanScaterplot'),
                               actionButton("DrawMeanScaterPlot","Draw Scaterplot"),
                               plotlyOutput('MeanScaterplotinteractive',width = "100%", height = "400"),
                               selectizeInput(
                                 inputId = "GenesMeanScaterplot",
                                 label = "Select genes to show",
                                 choices = NULL,
                                 multiple = TRUE
                               ),
                               textInput("ForceMeanScaterplot", "Force to avoid overlap of gene names", "1"),
                               plotOutput('MeanScaterplot',width = "100%", height = "400"),
                               selectInput(
                                 inputId = "MeanScaterPlotformat",
                                 label = "Format of downloaded plot",
                                 choices = c("svg", "png", "pdf"),
                                 selected = "svg"
                               ),
                               textInput("WidthMeanScaterPlot", "Width of downloaded plot (cm)", "15"),
                               textInput("HeightMeanScaterPlot", "Height of downloaded plot (cm)", "15"),
                               helpText("Plot is downloaded at 600dpi."),
                               downloadButton("downloadMeanScaterPlot", "Download")
                      )
              )
      ),
      # #####Statisticalplot ----
      tabItem(tabName = "StatPlot",
              selectInput(
                inputId = "ComparisonSelect",
                label = "Select the comparison",
                choices = NULL,
                selected = NULL
              ),
              selectInput(
                inputId = "ThreasholdType",
                label = "Select the type of threashold",
                choices = c("Pvalue","Qvalue"),
                selected = "Pvalue"
              ),
              textInput("SignifThreashold", "Significant threashold value", "0.05"),
              textInput("FCThreashold", "log2(Fold change) threashold value", "0"),
              textOutput("NbProtSignif"),
              tabBox( width = 0,

                      ###### VolcanoPlot ----
                      tabPanel("Volcanoplot",
                               fluidRow(
                                  box(title = "Colors",
                                    width = 4,
                                    colourpicker::colourInput(
                                       inputId = "UpColor",
                                       label = "Choose the color of up regulated proteins",
                                       value = "#FF0000"
                                     ),
                                     colourpicker::colourInput(
                                       inputId = "DownColor",
                                       label = "Choose the color of down regulated proteins",
                                       value = "limegreen"
                                     ),
                                     colourpicker::colourInput(
                                       inputId = "NSColor",
                                       label = "Choose the color of non significant proteins",
                                       value = "#C4C4C4"
                                     )    
                                 ),
                                 box(title="Line of threashold",
                                     width=4,
                                       checkboxInput(
                                         inputId = "ShowHline",
                                         label = "Add P/Qvalue threashhold line",
                                         value = FALSE
                                       ),
                                       checkboxInput(
                                       inputId = "ShowVline",
                                       label = "Add Log2(FC) threashhold line",
                                       value = FALSE
                                       ),
                                       selectInput(
                                         inputId = "LineType",
                                         label = "Select the type of the line",
                                         choices = c( "solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "1F", "F1", "4C88C488", "12345678"),
                                         selected = "solid"
                                       ),
                                       colourpicker::colourInput(
                                         inputId = "LineColor",
                                         label = "Choose the color of the line",
                                         value = "#000000"
                                       )
                                  )
                               ),
                               plotlyOutput('IntercativeVolcanoPlot',width = "100%", height = "400"),
                               fluidRow(
                                 box(title="Gene names",
                                        width=4,
                                   selectizeInput(
                                   inputId = "sel_gene_nm",
                                   label = "Select which gene names to show",
                                   choices = NULL,
                                   multiple = TRUE
                                 ),
                                 textInput("Force", "Force to avoid overlap of gene names", "1"),
                                 ),
                                 box(title="Gene points",
                                     width=4,
                                     selectizeInput(
                                       inputId = "GenePoints",
                                       label = "Select which gene points to show",
                                       choices = NULL,
                                       multiple = TRUE
                                     ),
                                     textInput("PointSize", "Points size", "1"),
                                     colourpicker::colourInput(
                                       inputId = "PointColor",
                                       label = "Points Color",
                                       value = "#000000"
                                     )
                                     
                                 )
                               ),
                               plotOutput('VolcanoPlot',width = "100%", height = "400"),
                               selectInput(
                                 inputId = "Volcanoformat",
                                 label = "Format of downloaded plot",
                                 choices = c("svg", "png", "pdf"),
                                 selected = "svg"
                               ),
                               textInput("WidthVolcano", "Width of downloaded plot (cm)", "15"),
                               textInput("HeightVolcano", "Height of downloaded plot (cm)", "15"),
                               helpText("Plot is downloaded at 600dpi."),
                               downloadButton("downloadVolcanoplot", "Download"),
                              
                              
                      ),
                      ###### Heatmap Significant----
                      tabPanel("Heatmap",
                               colourpicker::colourInput(
                                 inputId = "LowColorSignif",
                                 label = "Choose the color of low expression proteins",
                                 value = "#00FF00"
                               ),
                               colourpicker::colourInput(
                                 inputId = "MiddleColorSignif",
                                 label = "Choose the color of middle expression proteins",
                                 value = "#000000"
                               ),
                               colourpicker::colourInput(
                                 inputId = "HighColorSignif",
                                 label = "Choose the color of high expression proteins",
                                 value = "#FF0000"
                               ),
                               checkboxInput(
                                 inputId = "ShowGeneName",
                                 label = "Show gene names",
                                 value = FALSE
                               ),
                               helpText("WARRNING : You must clic on 'descriptive plot' before the heatmap to be display."),
                               plotOutput('HMSignif',width = "100%", height = "400"),
                               selectInput(
                                 inputId = "HMSignifformat",
                                 label = "Format of downloaded plot",
                                 choices = c("svg", "png", "pdf"),
                                 selected = "svg"
                               ),
                               textInput("WidthHMSignif", "Width of downloaded plot (cm)", "15"),
                               textInput("HeightHMSignif", "Height of downloaded plot (cm)", "15"),
                               helpText("Plot is downloaded at 600dpi."),
                               downloadButton("downloadHMSignif", "Download")
                      )
              )
      ),
      tabItem("ProtPlot",
              selectInput(
                inputId = "ProtPlot_GeneName",
                label = "Select the Gene to show",
                choices = NULL,
                selected = NULL
              ),
              tabBox( width = 0,
                      
                      ###### ProtPlot ----
                      tabPanel("Protein level",
                               textInput("ytitle", "y axis title", "Intensity"),
                               selectInput(
                                 inputId = "BarOrLine",
                                 label = "Select the type of plot display",
                                 choices = c("Bar plot", "Line plot"),
                                 selected = "Bar plot"
                               ),
                               plotlyOutput('ProtBarplot',width = "100%", height = "400"),
                               selectInput(
                                 inputId = "ProtBarPlotformat",
                                 label = "Format of downloaded plot",
                                 choices = c("svg", "png", "pdf"),
                                 selected = "svg"
                               ),
                               textInput("WidthProtBarPlot", "Width of downloaded plot (cm)", "10"),
                               textInput("HeightProtBarPlot", "Height of downloaded plot (cm)", "10"),
                               helpText("Plot is downloaded at 600dpi."),
                               downloadButton("downloadProtBarPlot", "Download")
                      ),
                      ###### PepPlot ----
                      tabPanel("Peptide level",
                               plotlyOutput('PepLinePlot',width = "100%", height = "400"),
                               selectInput(
                                 inputId = "PepLinePlotformat",
                                 label = "Format of downloaded plot",
                                 choices = c("svg", "png", "pdf"),
                                 selected = "svg"
                               ),
                               textInput("WidthPepLinePlot", "Width of downloaded plot (cm)", "30"),
                               textInput("HeightPepLinePlot", "Height of downloaded plot (cm)", "15"),
                               helpText("Plot is downloaded at 600dpi."),
                               downloadButton("downloadPepLinePlot", "Download")
                      )
              )
    )
  )
  )
)


####Server ----
server <- function(input, output, session) {

  output$value <- renderDataTable({
    out <- readRDS(input$dataFile$datapath)
    out[["edata"]]
  })

  edata<-eventReactive(input$dataFile,{
    out <- readRDS(input$dataFile$datapath)
    out[["edata"]]
  })
  
  
  output$TextNbProt<-renderText({
    paste(as.character(nrow(edata())),"proteins,",as.character(ncol(edata())),"samples.")
  })

  fdata<-eventReactive(input$dataFile,{
    out <- readRDS(input$dataFile$datapath)
    out<-out[["fdata"]]
    colnames(out)[which(colnames(out)=="Protein.Group"|colnames(out)=="Protein IDs")]<-"ProteinGroup"
    out$Gene_ProteinGroup<-paste(out$Genes,out$ProteinGroup)
    out
  })
  
  MyCond<-eventReactive(input$dataFile,{
    out <- readRDS(input$dataFile$datapath)
    out<-unique(out[["phenodata"]]$Condition)
    out
  })
  
  phenodata<-eventReactive(input$dataFile,{
    out <- readRDS(input$dataFile$datapath)
    out[["phenodata"]]
  })

  Pepdata<-eventReactive(input$dataFile,{
    out <- readRDS(input$dataFile$datapath)
    out<-out[["Pepdata"]]
    out[,3:ncol(out)]<-log2(out[,3:ncol(out)])
    out
  })

  observeEvent(input$dataFile$datapath, {
    out <- readRDS(input$dataFile$datapath)
    out<-out[["fdata"]]
    colnames(out)[which(colnames(out)=="Protein.Group"|colnames(out)=="Protein IDs")]<-"ProteinGroup"
    out$Gene_ProteinGroup<-paste(out$Genes,out$ProteinGroup)
    
    updateSelectizeInput(
      inputId = "sel_gene_nm",
      choices = out$Gene_ProteinGroup,
      server = TRUE,
      selected = NULL
    )
    updateSelectizeInput(
      inputId = "GenePoints",
      choices = out$Gene_ProteinGroup,
      server = TRUE,
      selected = NULL
    )
    updateSelectizeInput(
      inputId = "ProtPlot_GeneName",
      choices = out$Gene_ProteinGroup,
      server = TRUE,
      selected = NULL
    )
    updateSelectizeInput(
      inputId = "GenePCA",
      choices = out$Gene_ProteinGroup,
      server = TRUE,
      selected = NULL
    )
    updateSelectizeInput(
      inputId = "GenesScaterplot",
      choices = out$Gene_ProteinGroup,
      server = TRUE,
      selected = NULL
    )
    updateSelectizeInput(
      inputId = "GenesMeanScaterplot",
      choices = out$Gene_ProteinGroup,
      server = TRUE,
      selected = NULL
    )
  })
  #####Functions  ----  
  
  Exportggplot<-function(graph,basename,format,width,height){
    downloadHandler(
      filename = function() {
        paste0(basename,".", req(format))
      },
      content = function(file) {
        ggsave(file,
               device=req(format),
               plot=graph,
               dpi=600,
               width=as.numeric(width),
               height=as.numeric(height),
               units="cm")
      })
  }
  
  Exportheatmap<-function(graph,basename,format,width,height){
    downloadHandler(
    filename = function() {
      paste0(basename,".", req(format))
    },
    content = function(file) {
      if (format == "png") {
        ragg::agg_png(file,
                      res = 600,
                      width=as.numeric(width),
                      height=as.numeric(height), 
                      units = "cm"
        )
      } else if (format == "pdf") {
        pdf(file,
            width=as.numeric(width)*0.393701,
            height=as.numeric(height)*0.393701)
      } else if (format == "svg") {
        svglite::svglite(file,
                         width=as.numeric(width)/2,
                         height =as.numeric(height)/2)
      }
      ComplexHeatmap::draw(graph)
      dev.off()
    }
  )}
  
  Density<-function(x,y,points){
    #CalcDensityOnGrid
    MyData<-data.frame(x=x,y=y)
    MyData<-MyData[!is.na(MyData[,1])&!is.na(MyData[,2]),]
    
    xmin<-min(MyData[,1])
    xmax<-max(MyData[,1])
    ymin<-min(MyData[,2])
    ymax<-max(MyData[,2])
    
    #GetValuesOnGrid
    xvals<-MyData[,1]
    xStep<-(xmax - xmin) / points
    yvals<-MyData[,2]
    yStep<-(ymax - ymin) / points
    
    n = length(xvals)
    
    #CalcCovariance
    CalcCovariance<-function(data){
      n = nrow(data)
      p = ncol(data)
      
      means<-apply(data,2,mean)
      
      cov = matrix(0,nrow=p,ncol=p)
      for (i in 1:p){
        for (j in 1:i){
          cov[i, j] = sum((data[, i] - means[i]) * (data[, j] - means[j]))
          cov[i, j] = cov[i, j]/ n
          cov[j, i] = cov[i, j]
        }
      }
      return(cov)
    }
    
    cov<-CalcCovariance(MyData)#OK
    
    fact = n^(1 / 6)#OK
    hinv = fact/sqrt(abs(cov))
    hinv[c(1,2), 1] <-hinv[c(1,2),1]*xStep#OK
    hinv[c(1,2), 2] <-hinv[c(1,2),2]*yStep#OK
    
    dx = as.integer(1.0 / hinv[1, 1] * 5)#OK
    dy = as.integer(1.0 / hinv[2, 2] * 5)#OK
    
    values = matrix(0,ncol=points,nrow = points)
    
    xind = as.integer(floor((xvals - xmin) / xStep))#OK
    yind = as.integer(floor((yvals - ymin) / yStep))#OK
    
    for (i in 1:n){
      
      ii <-max(xind[i] - dx, 1): min(xind[i] + dx, points)
      jj <-max(yind[i] - dy, 1): min(yind[i] + dy, points)
      
      #MatrixTimesVector & StandardGaussian
      a1<-((hinv[1,1]*(ii - xind[i]))+(hinv[2,1]*(ii - xind[i])))^2
      a1<-matrix(rep(a1,length(jj)),nrow=length(ii))
      
      a2<-((hinv[1,2]*(jj - yind[i]))+(hinv[2,2]*(jj - yind[i])))^2
      a2<-matrix(rep(a2,each=length(ii)),ncol=length(jj))
      
      MySum<-a1+a2
      StandardGaussian<-exp(-0.5*MySum)/(2*pi)#OK
      
      values[ii, jj] <-values[ii, jj]+StandardGaussian
      
    }
    values<-as.data.frame(values)
    values[values == "NaN"]<-0 
    
    values<-values/max(values)#OK
    xmat <- seq(xmin, xmax, length.out = points)
    ymat <- seq(ymin, ymax, length.out = points)
    
    dvals = rep(0,length(x))
    for (i in 1:length(x)){
      xx = x[i]
      yy = y[i]
      if (!is.na(xx) && !is.na(yy)){
        xind = length(xmat[xmat<=xx])
        yind = length(ymat[ymat<=yy])
        dvals[i] = values[xind, yind]
      } else{
        dvals[i] = NA
      }
    }
    return(dvals)
  }
  
  DataScaterPlot<-function(x,y,edata,fdata){
    MyData<-data.frame(x=x,y=y)
    MyData$Density<-Density(x,y,300)
    MyData$edataRowNames<-row.names(edata)
    MyData<-merge(MyData,fdata,by.x="edataRowNames",by.y="RowNamesfdata",all.x=TRUE)
    return(MyData)
  }
  
  DrawScaterplot<-function(MyData,TitleX,TitleY){

    hmcol<- c(colorRampPalette(c("green","yellow"))(5),
              colorRampPalette(c("yellow","orange"))(5),
              colorRampPalette(c("orange","red"))(15),
              colorRampPalette(c("red","blue"))(25),
              colorRampPalette(c("blue","cyan"))(50))
    
    ggplot(data=MyData,aes_string(x=colnames(MyData)[2],y=colnames(MyData)[3],color="Density",label="Gene_ProteinGroup"))+
      geom_point(size=0.5)+
      scale_color_gradientn(colours = hmcol)+
      geom_abline(slope=1,intercept = 0,colour ="white",linewidth=1)+
      theme(panel.background = element_rect(fill = "black"),
            panel.grid.major=element_line(colour="grey40"),
            panel.grid.minor=element_line(colour="grey40"))+
      xlab(TitleX) +
      ylab(TitleY)
  }
  
  #####Descriptive plots  ----
    ######Choose Colors per condition----
  output$colors <- renderUI({
    req(MyCond())
    purrr::map2(
      MyCond(),
      colorRampPalette(brewer.pal(9, "Set1"))(if(length(MyCond())>9){length(MyCond())}else{9})[1:length(MyCond())],
      ~ colourpicker::colourInput(
        inputId = session$ns(.x),
        paste("Choose the color of : ", .x),
        value = .y
      )
    )
  })

  phenodata2<-eventReactive(
    magrittr::set_names(
    purrr::map_chr(
      MyCond(),
      ~ input[[.x]] %||% ""
    ),MyCond()),{
      color_by_level <- magrittr::set_names(
      purrr::map_chr(
        MyCond(),
        ~ input[[.x]] %||% ""
      ),MyCond())

    phenodata2<-phenodata()
    for(i in 1:length(MyCond())){
      phenodata2[phenodata2$Condition==MyCond()[i],"CondColor"]<- color_by_level[i]
    }
    phenodata2
  })
  
  ######NbProt----
  
  GraphNbProt<-reactive({
    req(edata(),phenodata(),MyCond())

    NbProt <- colSums(!is.na(edata()))
    
    x <- data.frame(Sample=colnames(edata()),
                    NbProt=NbProt,
                    Condition=phenodata2()$Condition,
                    CondColor=phenodata2()$CondColor)
    
    Conditions<-factor(x$Condition,levels=unique(x$Condition))
    MyColor<-unique(x$CondColor)
    MySample<-factor(x$Sample,levels=x$Sample)
    
      ggplot(x, aes(x=MySample, y=NbProt, color=Conditions,fill=Conditions)) +
      geom_bar(stat="identity")+
      ylim(0,NA)+
      scale_colour_manual(values = MyColor)+
      scale_fill_manual(values = MyColor)+
      theme_light()+
      theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.25),
            axis.title.x = element_blank())+
      ggtitle("Number of Protein per sample")

  }) 
  
  output$nbprot<-renderPlotly({ GraphNbProt() })

  output$downloadNbProtplot <-Exportggplot(graph=GraphNbProt(),
                                           basename="NbProt",
                                           format=input$NbProtformat,
                                           width=input$WidthNbProt,
                                           height=input$HeightNbProt)

  
  ######PCA----

  DataToShow <-eventReactive(input$DrawPCA,{
    req(edata(),phenodata())

    FilterGroup <- function(data, MinPc,type=c("1Group","EachGroup","Total")) {
      ValidProt<-rep(NA,nrow(data))
      
      
      if(type=="Total"){
        ValidProt<-100*rowSums(!is.na(edata()))/ncol(edata())
        ValidProt<-ValidProt>=MinPc
      }else{
        for(i in 1:nrow(data)){
          nbVV <- tapply(as.numeric(data[i,]),
                         INDEX=phenodata()[,"Condition"],
                         FUN=complete.cases)#show TRUE FALSE Valid Values for each sample per group
          nbVV <- lapply(nbVV,sum)#count number of Valid Values per group
          
          for(condition in levels(factor(phenodata()[,"Condition"]))) {
            nbVV[[condition]] <- nbVV[[condition]]/table(phenodata()[,"Condition"])[[condition]]*100#%Valid Values per condition.
          }
          
          if (type=="1Group"){
            ValidProt[i] <- any(nbVV >= MinPc)#return TRUE if there is at least one group with 70% of valid values.
          }else if(type=="EachGroup"){
            ValidProt[i] <- !any(nbVV < MinPc)
          }
        }
      }
      return(ValidProt)
    }
    
    if(input$FilterType=="% of Valid value in at least one group"){
      MyType<-"1Group"
    } else if (input$FilterType=="% of Valid value in each group"){
      MyType<-"EachGroup"
    } else if (input$FilterType=="% of Valid value in total"){
      MyType<-"Total"
    }
    
    edata()[FilterGroup(edata(), MinPc = as.numeric(input$PcFilter),type=MyType),]

  })
  
  output$NbProtPCA<-renderText({
    paste(as.character(nrow(DataToShow())),"/",as.character(nrow(edata())),"proteins used for PCA.")
  })
  
  my.prc <-eventReactive(input$DrawPCA,{
      req(edata(),phenodata(),DataToShow())

    DataToShow <- impute.MinProb(DataToShow(),q = 0.01,tune.sigma =1)
    
    if(input$TypeDataPCA=="Z-score(intensity)"){
      x<-colnames(DataToShow)
      DataToShow <-t(apply(DataToShow,1,scale))
      colnames(DataToShow)<-x  
      prcomp(t(DataToShow) , center=F, scale=F)
    }else if (input$TypeDataPCA=="Log2(intensity)"){
      prcomp(t(DataToShow) , center=T, scale=F)
    }


  })
  
  PCATitle <-eventReactive(input$DrawPCA,{
    req(edata(),phenodata(),DataToShow())
    paste("PCA", input$TypeDataPCA ,"of proteins with", input$PcFilter,input$FilterType)
  })
  
  GraphPCA<-reactive({
    req(edata(),phenodata(),MyCond())

    #PC1 versus PC2
      var_explained <- round(my.prc()$sdev^2/sum(my.prc()$sdev^2)*100,1)

      Conditions<-as.factor(phenodata2()$Condition)
      MyColors<-levels(Conditions)
      MyColors<-match(MyColors,phenodata2()$Condition)
      MyColors<-phenodata2()$CondColor[MyColors]
      
      x<-as.data.frame(my.prc()$x ) 
      x$Sample<-rownames(x)
      
        ggplot(x,aes(x=PC1,y=PC2,color=Conditions)) +
        geom_point(size=4) +
        scale_colour_manual(values = MyColors)+
        scale_fill_manual(values = MyColors)+
        geom_text_repel(data = x,
                        aes(label = Sample),
                        nudge_y = 2,size=5,
                        segment.size  = 0.5,
                        force=1)+
        geom_encircle(
          aes(group = Conditions,
              fill = Conditions),
          alpha = 0.1,
          size = 2,
          show.legend = FALSE,
          na.rm = TRUE)+
        theme_bw(base_size=17) + 
        labs(x=paste0("PC1: ",var_explained[1],"%"),
             y=paste0("PC2: ",var_explained[2],"%")) +
        theme(legend.position="right")+
        ggtitle(PCATitle())
  }) 
  
  output$PCA<-renderPlot({ GraphPCA() })

  output$downloadPCAplot <-Exportggplot(graph=GraphPCA(),
                                           basename="PCA",
                                           format=input$PCAformat,
                                           width=input$WidthPCA,
                                           height=input$HeightPCA)
  
  
  GraphPCAprot<-reactive({
    req(my.prc(), fdata())
    
    x<-as.data.frame(my.prc()$rotation ) 
    x$ProtNumber<-rownames(x)
    y<-fdata()[,c("Gene_ProteinGroup","Genes","RowNamesfdata")]
    x<-merge(x,y,by.x="ProtNumber",by.y="RowNamesfdata",all.x=TRUE, all.y=FALSE)

    ggplot(x,aes(x=PC1,y=PC2,label=Gene_ProteinGroup)) +
      geom_point(color="grey")+
      geom_text_repel(data = x[match(input$GenePCA,x$Gene_ProteinGroup),],
                      aes(label = Genes),
                      segment.size  = 0.5,
                      force=as.numeric(input$ForcePCA),
                      color = "black",
                      max.overlaps = Inf)+
      theme_bw()
  }) 
  
  output$PCAprotInteractive<-renderPlotly({ GraphPCAprot() })
  
  output$PCAprot<-renderPlot({ GraphPCAprot() })
  
  output$downloadPCAprotplot <-Exportggplot(graph=GraphPCAprot(),
                                        basename="PCAprot",
                                        format=input$PCAprotformat,
                                        width=input$WidthPCAprot,
                                        height=input$HeightPCAprot)
  
  ######Pearson correlations----
  ####### Heatmap ----
  GraphPearson<-reactive({
    req(edata(),phenodata(),MyCond())
    
    MyConditions<-data.frame(Conditions=phenodata2()$Condition,row.names=colnames(edata()))
    
    MycolorCond<-phenodata2()[,c(2,3)]
    MycolorCond<-distinct(MycolorCond)
    x<-MycolorCond[,2]
    names(x) <- MycolorCond[,1]
    MycolorCond<-list(x)
    names(MycolorCond)<-"Conditions"
    
    matDist <- cor(edata(),use = "pairwise.complete.obs")
    matDist[matDist == 1]<-NA 
    
    hmcol<- c(colorRampPalette(c("red","sienna1"))(250),
              colorRampPalette(c("sienna1","white"))(150),
              colorRampPalette(brewer.pal(9, 'GnBu'))(100))

    MyBreaks<-c(seq(from=0,to=0.4999,length.out=250),
                seq(from=0.5,to=0.7999,length.out=150),
                seq(from=0.8,to=1,length.out=100))

    pheatmap(matDist,  
                             color=hmcol,
                             breaks=MyBreaks,
                             border_color = NA,
                             cluster_rows=input$PearsonCluster,
                             cluster_cols = input$PearsonCluster,
                             legend_breaks=c(0, 0.5, 0.8, 0.9, 1),
                             annotation_col = MyConditions,
                             annotation_row = MyConditions,
                             annotation_colors=MycolorCond,
                             annotation_names_col = FALSE,
                             annotation_names_row =  FALSE,
                             show_rownames=F,
                             show_colnames=T,
                             main = "Pearson's correlation between log2(Intensity)of proteins",
                             na_col = "#000000",
                             fontsize = 13,
                             fontsize_col = 13,
                             angle_col = "90"
                          )
  }) 
  
  
  output$Pearson<-renderPlot({ 
    GraphPearson() 
    })
  
  output$downloadPearsonplot <-Exportheatmap(graph=GraphPearson(),
                                             basename="PearsonHeatmap",
                                             format=input$Pearsonformat,
                                             width=input$WidthPearson,
                                             height=input$HeightPearson)
  
  ####### Boxplot ----
  GraphPearsonBox<-reactive({
    req(edata(),phenodata(),MyCond())
    matDist <- cor(edata(),use = "pairwise.complete.obs")
    matDist[matDist == 1]<-NA
    matDist<-as.data.frame(matDist)
    matDist$Sample1<-row.names(matDist)
    matDist<-reshape(matDist,
                     idvar="Sample1",
                     varying = colnames(matDist)[1:(ncol(matDist)-1)],
                     times= colnames(matDist)[1:(ncol(matDist)-1)], 
                     direction = "long",
                     v.name="Pearson")
    matDist<-matDist[is.na(matDist$Pearson)==FALSE,]
    colnames(matDist)[2]<-"Sample2"
    matDist<-merge(x=matDist,
                   y=phenodata2(),
                   by.x="Sample1",
                   by.y="ShortName",
                   all.x=TRUE,
                   all.y=FALSE)
    colnames(matDist)[4:5]<-paste0(colnames(matDist)[4:5],"1")
    
    matDist<-merge(x=matDist,
                   y=phenodata2(),
                   by.x="Sample2",
                   by.y="ShortName",
                   all.x=TRUE,
                   all.y=FALSE)
    colnames(matDist)[6:7]<-paste0(colnames(matDist)[6:7],"2")
    
    for(i in 1:nrow(matDist)){
      matDist$Cond1_Cond2[i]<-paste(sort(c(matDist$Condition1[i],matDist$Condition2[i])),collapse = "/")
    }
    MyColor<-rep("lightgray",length(levels(factor(matDist$Cond1_Cond2))))
    
    MyNames<-boxplot(data=matDist,Pearson~Cond1_Cond2)$names
    MyNames<-str_split(MyNames,pattern="/")
    for(i in 1:length(MyNames)){
      if(MyNames[[i]][1]==MyNames[[i]][2]){
        #Find color condition
        x<-matDist[matDist$Condition1==MyNames[[i]][1],]$CondColor1
        MyColor[i]<-x
      }
    }
    
    matDist$Cond1_Cond2 <- as.factor(matDist$Cond1_Cond2)
    
    ggplot(matDist, aes(x=Cond1_Cond2, y=Pearson,fill=Cond1_Cond2)) + 
      geom_boxplot(outlier.shape = NA)+
      scale_fill_manual(values=MyColor)+
      geom_jitter(shape=16, position=position_jitter(0.2),size=0.5)+
      theme_light()+
      theme(axis.text.x = element_text(angle = 90,hjust=1),
            axis.title.x = element_blank())+
      theme(legend.position="none")+
      ggtitle("Boxplot of Pearson's correlation between samples") 
    
  }) 
  
  output$PearsonBox<-renderPlot({ GraphPearsonBox() })
  
  output$downloadPearsonBoxplot <-Exportggplot(graph=GraphPearsonBox(),
                                            basename="PearsonBox",
                                            format=input$PearsonBoxformat,
                                            width=input$WidthPearsonBox,
                                            height=input$HeightPearsonBox)
  
  
  ######Euclidean distances ----
  #######Heatmap ----
  GraphEuclid<-reactive({
    req(edata(),phenodata(),MyCond())

    MyConditions<-data.frame(Conditions=phenodata2()$Condition,row.names=colnames(edata()))
    
    MycolorCond<-phenodata2()[,c(2,3)]
    MycolorCond<-distinct(MycolorCond)
    x<-MycolorCond[,2]
    names(x) <- MycolorCond[,1]
    MycolorCond<-list(x)
    names(MycolorCond)<-"Conditions"
    
    matDist <- as.matrix(dist(t(edata())))
    matDist[matDist == 0]<-NA 
    
    hmcol<- rev(colorRampPalette(brewer.pal(9, 'GnBu'))(100))
    
    pheatmap(matDist,  
             color=hmcol,
             border_color = NA,
             cluster_rows=input$EuclidCluster,
             cluster_cols = input$EuclidCluster,
             annotation_col = MyConditions,
             annotation_row = MyConditions,
             annotation_colors=MycolorCond,
             annotation_names_col = FALSE,
             annotation_names_row =  FALSE,
             show_rownames=F,
             show_colnames=T,
             main = "Euclidean distance between log2(Intensity)of proteins",
             na_col = "#000000",
             fontsize = 13,
             fontsize_col = 13,
             angle_col = "90"
    )
  }) 
  
  
  output$Euclid<-renderPlot({ GraphEuclid() })
  
  output$downloadEuclidplot <-Exportheatmap(graph=GraphEuclid(),
                                             basename="EuclideanHeatmap",
                                             format=input$Euclidformat,
                                             width=input$WidthEuclid,
                                             height=input$HeightEuclid)
  
  #######BoxPlot ----
  GraphEuclidBox<-reactive({
    req(edata(),phenodata(),MyCond())
    matDist <- as.matrix(dist(t(edata())))
    matDist[matDist == 0]<-NA 
    matDist<-as.data.frame(matDist)
    matDist$Sample1<-row.names(matDist)
    matDist<-reshape(matDist,
                     idvar="Sample1",
                     varying = colnames(matDist)[1:(ncol(matDist)-1)],
                     times= colnames(matDist)[1:(ncol(matDist)-1)], 
                     direction = "long",
                     v.name="Euclidean.distance")
    matDist<-matDist[is.na(matDist$Euclidean.distance)==FALSE,]
    colnames(matDist)[2]<-"Sample2"
    matDist<-merge(x=matDist,
                   y=phenodata2(),
                   by.x="Sample1",
                   by.y="ShortName",
                   all.x=TRUE,
                   all.y=FALSE)
    colnames(matDist)[4:5]<-paste0(colnames(matDist)[4:5],"1")
    
    matDist<-merge(x=matDist,
                   y=phenodata2(),
                   by.x="Sample2",
                   by.y="ShortName",
                   all.x=TRUE,
                   all.y=FALSE)
    colnames(matDist)[6:7]<-paste0(colnames(matDist)[6:7],"2")
    
    for(i in 1:nrow(matDist)){
      matDist$Cond1_Cond2[i]<-paste(sort(c(matDist$Condition1[i],matDist$Condition2[i])),collapse = "/")
    }
    
    MyColor<-rep("lightgray",length(levels(factor(matDist$Cond1_Cond2))))
    
    MyNames<-boxplot(data=matDist,Euclidean.distance~Cond1_Cond2)$names
    MyNames<-str_split(MyNames,pattern="/")
    for(i in 1:length(MyNames)){
      if(MyNames[[i]][1]==MyNames[[i]][2]){
        #Find color of the condition
        MyColor[i]<-matDist[matDist$Condition1==MyNames[[i]][1],]$CondColor1
      }
    }
    
    matDist$Cond1_Cond2 <- as.factor(matDist$Cond1_Cond2)
    
    ggplot(matDist, aes(x=Cond1_Cond2, y=Euclidean.distance,fill=Cond1_Cond2)) + 
      geom_boxplot(outlier.shape = NA)+
      scale_fill_manual(values=MyColor)+
      geom_jitter(shape=16, position=position_jitter(0.2),size=0.5)+
      theme_light()+
      theme(axis.text.x = element_text(angle = 90,hjust=1),
            axis.title.x = element_blank())+
      theme(legend.position="none")+
      ggtitle("Boxplot of Euclidean distance between samples") 
    
    
  }) 
  
  output$EuclidBox<-renderPlot({ 
    GraphEuclidBox() 
  })
  
  output$downloadEuclidBoxplot <-Exportggplot(graph=GraphEuclidBox(),
                                               basename="EuclideanBox",
                                               format=input$EuclidBoxformat,
                                               width=input$WidthEuclidBox,
                                               height=input$HeightEuclidBox)
  
  ######Heatmap Log2 ----
  GraphLog2<-reactive({
    req(edata(),phenodata(),MyCond())
  
    # Select proteins with at least 1 valid value in total
    x <- edata()[rowSums(is.na(edata())) != ncol(edata()), ]
    
    # Mean intensity
    x$MeanLFQ<-apply(x,1,function(x){mean(x,na.rm = TRUE)})
    
    # Sort proteins by mean intensity
    x<-x[order(x$MeanLFQ,decreasing = TRUE),]
    
    # Delete mean columns
    x<-x %>% select(-'MeanLFQ')
    
    #select 2000 prot if nbprot>2000
    nbprot<-nrow(x)
    if (nbprot>2000){
      step<-ceiling(nbprot/2000)
      x<-x[seq(from=1, to=nbprot,by=step ),]
    }
    
    hmcol<- colorRampPalette(c(input$LowColor,input$MiddleColor,input$HighColor))(100)
    
    MyConditions<-data.frame(Conditions=phenodata2()$Condition,row.names=colnames(edata()))
    
    MycolorCond<-phenodata2()[,c(2,3)]
    MycolorCond<-distinct(MycolorCond)
    y<-MycolorCond[,2]
    names(y) <- MycolorCond[,1]
    MycolorCond<-list(y)
    names(MycolorCond)<-"Conditions"
    
    
    #heatmap with sample clustering
    pheatmap(x,
             cluster_rows=FALSE,
             cluster_cols = input$SampleCluster,
             clustering_distance_cols = "euclidean",
             clustering_method = "average",
             col=hmcol,
             annotation_col = MyConditions,
             annotation_colors = MycolorCond,
             na_col = "grey",
             show_rownames=F,
             show_colnames=T,
             annotation_names_col = FALSE,
             annotation_names_row =  FALSE,
             main = "Heatmap of Log2(Intensity)",
             fontsize = 8,
             fontsize_col = 8,
             border_color = NA,
             angle_col="90")

  }) 


  output$HeatmapLog2<-renderPlot({ GraphLog2() })

  output$downloadHeatmapLog2plot <-Exportheatmap(graph=GraphLog2(),
                                            basename="HeatmapLog2",
                                            format=input$HeatmapLog2format,
                                            width=input$WidthHeatmapLog2,
                                            height=input$HeightHeatmapLog2)
  
  
  ######Scaterplot ----
  observeEvent(input$dataFile, {
    out<-colnames(edata())
    out<-as.list(out)
    names(out)<-out
    
    updateSelectInput(
      inputId = "Sample1",
      choices = out,
      selected=out[[1]]
    )
    
    updateSelectInput(
      inputId = "Sample2",
      choices = out,
      selected=out[[2]]
    )
    
  })
  
  output$CorScaterplot<-renderText({
    x<-edata()[,which(colnames(edata())==input$Sample1)]
    y<-edata()[,which(colnames(edata())==input$Sample2)]
    paste("r =",round(cor(x,y,use="pairwise.complete.obs"),5))
  })
  
  
  DensityScater <-eventReactive(input$DrawScaterPlot,{
    req(edata(),fdata())
    
    x<-edata()[,which(colnames(edata())==input$Sample1)]
    y<-edata()[,which(colnames(edata())==input$Sample2)]
    
    DataScaterPlot(x,y,edata(),fdata())

  })
  
  DrawScater<-eventReactive(input$DrawScaterPlot,{
    req(DensityScater())
    DrawScaterplot(DensityScater(),TitleX=input$Sample1,TitleY=input$Sample2)
    
  })
  
  ScaterPlot<-reactive({
    req(DrawScater())
    
    g<-DrawScater()
    g+geom_text_repel(data = DensityScater()[match(input$GenesScaterplot,DensityScater()$Gene_ProteinGroup),],
                      aes(label = Genes),
                      segment.size  = 0.5,
                      force=as.numeric(input$ForceScaterplot),
                      color = "white",
                      max.overlaps = Inf)


  })
  
  output$Scaterplotinteractive<-renderPlotly({ ScaterPlot() })
  
  output$Scaterplot<-renderPlot({ ScaterPlot() })

  output$downloadScaterPlot <-Exportggplot(graph=ScaterPlot(),
                                              basename=paste0("ScaterPlot_",input$Sample1,"_",input$Sample2),
                                              format=input$ScaterPlotformat,
                                              width=input$WidthScaterPlot,
                                              height=input$HeightScaterPlot)
  
  ######MeanScaterplot ----
  observeEvent(input$dataFile, {
    out<-as.list(MyCond())
    names(out)<-out
    
    updateSelectInput(
      inputId = "Cond1",
      choices = out,
      selected=out[[1]]
    )
    
    updateSelectInput(
      inputId = "Cond2",
      choices = out,
      selected=out[[2]]
    )
    
  })
  
  output$CorMeanScaterplot<-renderText({
    x<-edata()[,which(phenodata2()$Condition == input$Cond1)]
    x<-apply(x,MARGIN=1,mean, na.rm = TRUE)
    y<-edata()[,which(phenodata2()$Condition == input$Cond2)]
    y<-apply(y,MARGIN=1,mean, na.rm = TRUE)
    paste("r =",round(cor(x,y,use="pairwise.complete.obs"),5))
  })

  
  DensityMeanScater  <-eventReactive(input$DrawMeanScaterPlot,{
    req(edata(),fdata())
    
    x<-edata()[,which(phenodata2()$Condition == input$Cond1)]
    x<-apply(x,MARGIN=1,mean, na.rm = TRUE)
    y<-edata()[,which(phenodata2()$Condition == input$Cond2)]
    y<-apply(y,MARGIN=1,mean, na.rm = TRUE)
    
    DataScaterPlot(x,y,edata(),fdata())
    
  })
  
  DrawMeanScater<-eventReactive(input$DrawMeanScaterPlot,{
    req(DensityMeanScater())
    DrawScaterplot(DensityMeanScater(),TitleX=input$Cond1,TitleY=input$Cond2)
    
  })

  MeanScaterPlot<-reactive({
    req(DrawMeanScater())
    
    g<-DrawMeanScater()
    g+geom_text_repel(data = DensityMeanScater()[match(input$GenesMeanScaterplot,DensityMeanScater()$Gene_ProteinGroup),],
                      aes(label = Genes),
                      segment.size  = 0.5,
                      force=as.numeric(input$ForceMeanScaterplot),
                      color = "white",
                      max.overlaps = Inf)
    
    
  })

  output$MeanScaterplotinteractive<-renderPlotly({ MeanScaterPlot() })
  
  output$MeanScaterplot<-renderPlot({ MeanScaterPlot() })

  output$downloadMeanScaterPlot <-Exportggplot(graph=MeanScaterPlot(),
                                           basename=paste0("MeanScaterPlot_",input$Cond1,"_",input$Cond2),
                                           format=input$MeanScaterPlotformat,
                                           width=input$WidthMeanScaterPlot,
                                           height=input$HeightMeanScaterPlot)
  
  #####Statistical plots  ----
  observeEvent(input$dataFile$datapath, {
    out <- readRDS(input$dataFile$datapath)
    out<-out[["fdata"]]
    outComp<-colnames(out)[str_detect(colnames(out),"Significant")]
    outComp<-gsub("Significant ","",outComp)
    out<-as.list(outComp)
    names(out)<-outComp

    updateSelectInput(
      inputId = "ComparisonSelect",
      choices = out
    )
  })
  
  output$NbProtSignif<-renderText({
    req(MyListeProt(),MyColNumber())
    
    Pval<-max(fdata()[which(MyListeProt()==TRUE),MyColNumber()])
    Qval<-max(fdata()[which(MyListeProt()==TRUE),(MyColNumber()+1)])
    paste(as.character(length(which(MyListeProt()==TRUE))),"/",as.character(length(MyListeProt())),"significant proteins.","\nPvalue=",Pval,",FDR=",Qval,".")
  })
  
  ######Volcanoplot ----

  GraphVolcano<-reactive({
    req(fdata(),MyColNumber())
    
    n<-if(input$ThreasholdType=="Pvalue"){0}else{1}
    
      x<-fdata()[is.na(fdata()[,MyColNumber()])==FALSE,]
      x<-x[,c(MyColNumber(),MyColNumber()+n,MyColNumber()+2,5,ncol(x))]
      x$Type<-input$NSColor
      x[x[,2]<as.numeric(input$SignifThreashold) & x[,3]>as.numeric(input$FCThreashold),"Type"]<-input$UpColor
      x[x[,2]<as.numeric(input$SignifThreashold) & x[,3]<(-as.numeric(input$FCThreashold)),"Type"]<-input$DownColor
      x[,2]<--log10(x[,2])
      colnames(x)<-gsub(pattern="\\s",replacement="_",x=colnames(x))
      colnames(x)<-gsub(pattern="\\/",replacement="_",x=colnames(x))
      
      g<-ggplot(x,aes_string(x=colnames(x)[3],y=colnames(x)[2],color="Type",label="Gene_ProteinGroup")) +
        geom_point(size=1) +
        scale_colour_manual(values = levels(as.factor(x$Type)))+
        geom_point(data=x[match(input$GenePoints,x$Gene_ProteinGroup),],
                   color=input$PointColor,
                   size=as.numeric(input$PointSize))+
        geom_text_repel(data = x[match(input$sel_gene_nm,x$Gene_ProteinGroup),],
                        aes_string(label = "Genes"),
                        segment.size  = 0.5,
                        force=as.numeric(input$Force),
                        color = "black",
                        max.overlaps = Inf)+
        theme_bw(base_size=10) + 
        labs(x=paste0("Log2(",input$ComparisonSelect,")"),
             y=paste0("-log10(",input$ThreasholdType,")")) +
        theme(legend.position="none")+
        ggtitle(paste("Volcano plot",input$ComparisonSelect))
      
      if(input$ShowHline==TRUE){
        g<-g+geom_hline(linetype=input$LineType,
                        color=input$LineColor,
                        yintercept=-log10(as.numeric(input$SignifThreashold))
                        )
      }
      
      if(input$ShowVline==TRUE){
        g<-g+geom_vline(linetype=input$LineType,
                        color=input$LineColor,
                        xintercept=as.numeric(input$FCThreashold)
                        )+
            geom_vline(linetype=input$LineType,
                        color=input$LineColor,
                        xintercept=-as.numeric(input$FCThreashold)
                        )
      }
    g
  })
  
  output$VolcanoPlot<-renderPlot({ GraphVolcano() })
  
  output$IntercativeVolcanoPlot<-renderPlotly({ GraphVolcano() })
  
  output$downloadVolcanoplot <-Exportggplot(graph=GraphVolcano(),
                                           basename=paste0("Volcanoplot_",input$ComparisonSelect),
                                           format=input$Volcanoformat,
                                           width=input$WidthVolcano,
                                           height=input$HeightVolcano)

  
  ######Heatmap significant ----
  MyColNumber <-reactive({
    req(fdata())
    MyColNumber<-which(str_detect(colnames(fdata()),paste("Pvalue",input$ComparisonSelect)))
    MyColNumber
  })
  
  MyListeProt <-reactive({
    req(fdata(),MyColNumber())
    n<-if(input$ThreasholdType=="Pvalue"){0}else{1}
    
    #filter significative proteins
    MyListeProt<-(fdata()[,(MyColNumber()+n)]<as.numeric(input$SignifThreashold)
                  & (fdata()[,(MyColNumber()+2)]<(-as.numeric(input$FCThreashold ))
                     |fdata()[,(MyColNumber()+2)]>as.numeric(input$FCThreashold)))
    MyListeProt
    })
  
  
  GraphHMSignif<-reactive({
    req(fdata(),phenodata2(),MyListeProt(),edata())

    #select columns of the comparison
    MyListX<-which(phenodata2()$Condition == str_split(input$ComparisonSelect,"/")[[1]][1])
    MyListY<-which(phenodata2()$Condition == str_split(input$ComparisonSelect,"/")[[1]][2])

    #filter significative proteins
    Genes<-fdata()[which(MyListeProt()==TRUE),c(5,2)]#"Genes","Protein.Group"
    Genes[is.na(Genes$Genes),"Genes"]<-Genes[is.na(Genes$Genes),2]#"Protein.Group"
    Genes<-Genes$Genes
    
    MyListeProt<-fdata()[which(MyListeProt()==TRUE),"RowNamesfdata"]
    MyListeProt<-data.frame(x1=MyListeProt,x2=MyListeProt)
    
    DataToShow<-cbind(edata(),RowNamesedata=row.names(edata()))
    DataToShow<-merge(DataToShow,MyListeProt,by.x="RowNamesedata",by.y="x1",all=FALSE)
    DataToShow<-DataToShow %>% select(-c('RowNamesedata','x2'))
    
    DataToShow<-DataToShow[,c(MyListX,MyListY)]
    
    if(nrow(MyListeProt)>1){  
      #Z-score
      MyZscore<-t(apply(DataToShow,1,scale))
      colnames(MyZscore)<-colnames(DataToShow)
      rownames(MyZscore)<-Genes
      
      hmcol<- colorRampPalette(c(input$LowColorSignif,input$MiddleColorSignif,input$HighColorSignif))(100)
      hmcol<-c(input$LowColorSignif,hmcol,input$HighColorSignif)
      
      MyBreaks<-seq(from=-2,to=2,length.out=100)
      MinBreaks<-if(min(MyZscore,na.rm=TRUE)<(-2)){min(MyZscore,na.rm=TRUE)}else{-2.0001}
      MaxBreaks<-if(max(MyZscore,na.rm=TRUE)>2){max(MyZscore,na.rm=TRUE)}else{2.0001}
      MyBreaks<-c(MinBreaks,MyBreaks,MaxBreaks)
      
      MyConditions<-data.frame(Conditions=phenodata2()$Condition[c(MyListX,MyListY)],row.names=colnames(MyZscore))
      
      MycolorCond<-phenodata2()[c(MyListX[1],MyListY[1]),c(2,3)]
      x<-MycolorCond[,2]
      names(x) <- MycolorCond[,1]
      MycolorCond<-list(x)
      names(MycolorCond)<-"Conditions"
      
      pheatmap(MyZscore,
               cluster_rows=T,
               cluster_cols = T,
               color=hmcol,
               breaks=MyBreaks,
               legend_breaks=c(if(MinBreaks==-2.0001){-2}else{c(round(MinBreaks,1),-2)},
                               if(MaxBreaks==2.0001){2}else{c(2,round(MaxBreaks,1))}),
               annotation_col = MyConditions,
               annotation_colors=MycolorCond,
               na_col = "grey",
               show_rownames=input$ShowGeneName,
               show_colnames=TRUE,
               annotation_names_col = FALSE,
               annotation_names_row =  FALSE,
               main = paste0("Heatmap of Zscore of differential proteins ",
                             input$ComparisonSelect,
                             " (",input$ThreasholdType,input$SignifThreashold,if(input$FCThreashold!=0){paste0(", FC>",2^as.numeric(input$FCThreashold))}else{""},
                             ")"
               ),
               border_color = NA,
               angle_col="90"
               )
    }
  })
  
  output$HMSignif<-renderPlot({ GraphHMSignif() })
  
  output$downloadHMSignif <-Exportheatmap(graph=GraphHMSignif(),
                                                 basename=paste0("HeatmapSignificant_",input$ComparisonSelect),
                                                 format=input$HMSignifformat,
                                                 width=input$WidthHMSignif,
                                                 height=input$HeightHMSignif)
  
  
  #####Proteins plots  ----
  ###### ProtPlot ----
  df<-eventReactive(input$dataFile,{
    Myedata<-2^edata()
    Myedata<-cbind(Myedata,RowNamesedata=row.names(edata()))
    Myfdata<-fdata()[,c("RowNamesfdata","Gene_ProteinGroup")]
    df<-merge(Myfdata,Myedata,by.x="RowNamesfdata",by.y="RowNamesedata")
    df<-df[,-1]
    df<-reshape(data=df,
                idvar="Gene_ProteinGroup",
                varying = colnames(df)[2:ncol(df)],
                times=colnames(df)[2:ncol(df)],
                v.names="Intensity",
                direction="long")
    
    #Ajoute conditions
    Myphenodata<-phenodata()[,1:2]
    df<-merge(df,Myphenodata,by.x="time",by.y="ShortName",all.x=TRUE)
    colnames(df)[1]<-"Sample"
    df<-df[!is.na(df$Intensity),]
    df
  })
  
  GraphProt<-reactive({
    req(df(),MyCond())
    i<-input$ProtPlot_GeneName
    x<-df()[df()$Gene_ProteinGroup==i,]
    
    for(j in MyCond()){
      if(length(which(x$Condition==j))==0){
        y<-data.frame(Sample=NA,
                      Gene_ProteinGroup=i,
                      Intensity=0,
                      Condition=j)
        x<-rbind(x,y)
      }
    }
    
    #Calcul SD et Mean
    x.summary <- x %>%
      group_by(Condition) %>%
      summarise(
        sd = sd(Intensity, na.rm = TRUE),
        Intensity = mean(Intensity)
      )
    x.summary$Sample<-NA
    
    x<-x[x$Intensity>0,]
    
    g<-ggplot(x, aes(x=Condition, y=Intensity, label=Sample)) +
      geom_errorbar( aes(ymin = Intensity-sd, ymax = Intensity+sd), 
                     data = x.summary, width = 0.3) +
      theme_light()+
      ylab(input$ytitle)+
      theme(axis.title.x=element_blank(),
            axis.text.x=element_text(angle=90, hjust=1),
            plot.title = element_text(size=8))+
      ggtitle(i)
    
    if(input$BarOrLine=="Bar plot"){
        g<-g+geom_col(data = x.summary, fill = "grey", color = "grey50") +
        geom_point(pch=16,size=4,alpha=0.4, color = "black")
      }else{
        g<-g+geom_line(data=x.summary,aes(x=Condition,y=Intensity,group=1))+
        geom_point(data=x.summary,pch=16,size=4)
      }
    g
  })
  
  output$ProtBarplot<-renderPlotly({ GraphProt() })

  output$downloadProtBarPlot <-Exportggplot(graph=GraphProt(),
                                            basename=input$ProtPlot_GeneName,
                                            format=input$ProtBarPlotformat,
                                            width=input$WidthProtBarPlot,
                                            height=input$HeightProtBarPlot)
  
  
  ###### PeptidePlot ----
  GraphPep<-reactive({
    req(Pepdata())
    i<-input$ProtPlot_GeneName
    x<-Pepdata()[Pepdata()$Gene_ProteinGroup==i,]
    
    x<-reshape(data=x,
               idvar=c("Gene_ProteinGroup","Peptide"),
               varying=colnames(x)[3:ncol(x)],
               times=colnames(x)[3:ncol(x)],
               v.names="Log2Intensity",
               direction="long")
    colnames(x)[3]<-"Sample"
    x$OxyM<-" Unmodified"
    x$OxyM[str_detect(string=x$Precursor.Id,pattern="(UniMod:35)")]<-"OxyM"
    
    #Pepide plot
    g<-ggplot(x,aes(x=Sample,y=Log2Intensity,group=Peptide,color=Peptide))+
      geom_line(aes(linetype=OxyM),linewidth=1.1)+
      geom_point(size=2)+
      theme_light()+
      theme(axis.title.x=element_blank(),
            axis.text.x=element_text(angle=90, hjust=1),
            plot.title = element_text(size=8),
            legend.position = "none")+
      ggtitle(i)
    
    #Add prot data
    Myfdata<-fdata()
    y<-Myfdata[Myfdata$Gene_ProteinGroup==i,"RowNamesfdata"]
    y<-edata()[row.names(edata())==y,]
    y<-as.data.frame(t(y))
    y$Sample<-row.names(y)
    colnames(y)[1]<-"Log2Intensity"
    y$Peptide<-"Protein"
    
    g<-g+geom_line(data=y,aes(x=Sample,y=Log2Intensity),color="black",linewidth=1.2,linetype="dotted")+
      geom_point(data=y,pch=1,size=3,color="black")
    g
    
  })
  
  output$PepLinePlot<-renderPlotly({ GraphPep() })

  output$downloadPepLinePlot <-Exportggplot(graph=GraphPep(),
                                            basename=paste0(input$ProtPlot_GeneName,"_peptides"),
                                            format=input$PepLinePlotformat,
                                            width=input$WidthPepLinePlot,
                                            height=input$HeightPepLinePlot)

}


####Run----
shinyApp(ui, server)