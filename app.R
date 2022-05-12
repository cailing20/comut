library(shiny)
library(data.table)
library(DT)
library(ggplot2)
library(plotly)
library(ggrepel)
library(RColorBrewer)
library(scales)
library(shinyjs)
library(shinythemes)
source('load_data.R')
# Define UI for application that draws a histogram
ui <- fluidPage(
    navbarPage(title = "Co-mutation check with AACR-GENIE v11.0 data",theme = shinytheme("paper"),
               tabPanel(title = "Review screening results",
                        radioButtons(inputId = 'screen.w.CNA',label = 'Count in copy number alterations',choices = c('no','yes'),inline = T),
                        fluidRow(
                            column(5,
                                   h5('Number of significant hits by gene pair'),
                                   DT::dataTableOutput('gp-dt'),downloadButton(outputId = 'gp-download',label = 'Download gene pair-specific results')
                                   ),
                            column(7,
                                   h5('Number of significant hits by cancer type'),
                                   DT::dataTableOutput('type-dt'),downloadButton(outputId = 'type-download',label = 'Download cancer type-specific results')
                            )
                        ),
                        fluidRow(
                            h5('All significant hits'),
                            DT::dataTableOutput('full-dt'),downloadButton(outputId = 'full-download',label = 'Download all significant results')
                        )
               ),
               
               tabPanel(title = "Adhoc analysis for a specific gene pair",
                        fluidRow(
                            column(2,selectInput(inputId = 'g1',label = "Select gene 1",choices = sort(g.ref$gene),selected = 'TP53')),
                            column(2,selectInput(inputId = 'g2',label = "Select gene 2",choices = sort(g.ref$gene),selected = 'RB1')),
                            column(3,radioButtons(inputId = 'CNA',label = "Count in copy number alterations",choices = c("no","yes"),inline = T)),
                            column(2,radioButtons(inputId = 'show',label = "Display",choices = c("table","scatter plot"),inline = T)),
                            column(1,actionButton(inputId = 'submit',label = 'Submit')),
                            column(2,downloadButton(outputId = 'result-download',label = 'Download result table'))
                        ),
                        fluidRow(
                            useShinyjs(),
                            DT::dataTableOutput('dt'),
                            plotlyOutput(outputId = 'scPlot',height = '900px')
                        )
               ))
)

# Define server logic required to draw a histogram
server <- function(input, output,session) {
    # screen results
    output$`gp-dt`<-DT::renderDataTable(DT::datatable(result.list[[ifelse(input$screen.w.CNA=='yes','w','wo')]][['comut']],options = list(pageLength=5,scrollX = TRUE),filter='top',rownames = F))
    output$`type-dt`<-DT::renderDataTable(DT::datatable(result.list[[ifelse(input$screen.w.CNA=='yes','w','wo')]][['type']],options = list(pageLength=5,scrollX = TRUE),filter='top',rownames = F))
    output$`full-dt`<-DT::renderDataTable(DT::datatable(result.list[[ifelse(input$screen.w.CNA=='yes','w','wo')]][['full']],options = list(pageLength=5),filter='top',rownames = F))
    
    output$`gp-download` <- downloadHandler(
        filename=function(){paste0('significant_counts_by_gene_pair_from_screening_',ifelse(input$screen.w.CNA=='yes','','not_'),'considering_copy_number_alterations.csv')},content = function(file){fwrite(result.list[[ifelse(input$screen.w.CNA=='yes','w','wo')]][['comut']],file)}
    )
    output$`type-download` <- downloadHandler(
        filename=function(){paste0('significant_counts_by_cancer_type_from_screening_',ifelse(input$screen.w.CNA=='yes','','not_'),'considering_copy_number_alterations.csv')},content = function(file){fwrite(result.list[[ifelse(input$screen.w.CNA=='yes','w','wo')]][['type']],file)}
    )
    output$`full-download` <- downloadHandler(
        filename=function(){paste0('significant_results_from_screening_',ifelse(input$screen.w.CNA=='yes','','not_'),'considering_copy_number_alterations.csv')},content = function(file){fwrite(result.list[[ifelse(input$screen.w.CNA=='yes','w','wo')]][['full']],file)}
    )
    # adhoc analyses
    output$dt<-NULL;output$scPlot<-NULL;shinyjs::hide('result-download')
    co.df<-reactive({
        test.comut(g1 = input$g1,g2=input$g2,include.CNA = (input$CNA=='yes'))
    })
    co.g<-reactive({
        plot.gene(g1 = input$g1,g2=input$g2,by_cancer.df = co.df())
    })
    observeEvent(list(input$g1,input$g2),{shinyjs::hide('result-download')})
    observeEvent(input$submit,{
        if(input$g1==input$g2) showNotification("Please choose two different genes.")else{
            output$dt<-DT::renderDataTable(DT::datatable({isolate(co.df())},options = list(pageLength=10),filter='top',rownames = F))
            output$scPlot<-renderPlotly(ggplotly(isolate(co.g())))
            if(input$show=='table'){
                shinyjs::show('dt')
                shinyjs::hide('scPlot')
                shinyjs::show('result-download')
            }else{
                shinyjs::hide('dt')
                shinyjs::show('scPlot')
                shinyjs::hide('result-download')
            }
        }
    })
    output$`result-download` <- downloadHandler(
        filename=function(){paste0(input$g1,'_',input$g2,'_comutation_result.csv')},content = function(file){fwrite(co.df(),file)}
    )
}

# Run the application 
shinyApp(ui = ui, server = server)
