library(shiny)
library(Rsamtools)
library(rlist)
library(tidyverse)
library(jsonlite)
library(docopulae)
library(shinythemes) 
library(htmlwidgets) 
library(shinyWidgets) 
library(shinydashboard)

#load("data/genes_table_npm1_nk.RData")
#load("data/genes_table_npm1_imp.RData")
source("functions.R")

runApp(list(
  ui = fluidPage(
    #headerPanel(""),
    fluidRow(style='padding:20px;',
             conditionalPanel("output.show_input",
             sidebarPanel(
               fileInput("file1", "Choose BAM File",
                         multiple = FALSE,
                         accept = c(".bam")),
               
               fileInput("file2", "Choose SAM File",
                         multiple = FALSE,
                         accept = c(".sam")),
               
               actionButton(inputId = "getFiles", label = "Get annotation")
             ))
             ,
             
    ),
    fluidRow(
      mainPanel(
        conditionalPanel("output.show_output",
          fluidRow(
            column(width = 6, offset = 0,style='padding:20px;',
                   div(
                     uiOutput("genes_list"),
                     uiOutput("region_list")
                   )    
            ),
            column(width = 6, offset = 0,style='padding:35px;',
                   div(
                     #uiOutput("count_table")
                     uiOutput("chrm_list"),
                     uiOutput("strand_list")
                   )
            )
          ),
          fluidRow(
            div(style='padding:20px;',
                verbatimTextOutput("placeholder"),
                actionButton("btn_back", "Clear")
            )
          )
        )
      )
    )
    
  ),
  
  server = function(input, output, session) {
    
    ShowInput <- reactiveVal(TRUE)
    output$show_input <- reactive({
      ShowInput() 
    })
    outputOptions(output, "show_input", suspendWhenHidden = FALSE)
    
    ShowOutput <- reactiveVal(FALSE)
    output$show_output <- reactive({
      ShowOutput() 
    })
    outputOptions(output, "show_output", suspendWhenHidden = FALSE)
    
    observeEvent(input$getFiles, {
      ShowInput(FALSE)
      ShowOutput(TRUE)
    })
    
    observeEvent(input$btn_back, {
      ShowInput(TRUE)
      ShowOutput(FALSE)
    })
    
    rv <- reactiveValues(data_table = NULL)
    
    observeEvent(input$getFiles, {
      
      req(input$file1)
      req(input$file2)
      
      tryCatch(
        {
          files_annot$bam <- scanBam(input$file1$datapath)
          
        },
        error = function(e) {
          # return a safeError if a parsing error occurs
          stop(safeError(e))
        }
      )
      
      tryCatch(
        {
          files_annot$sam <- read.csv(input$file2$datapath, comment.char = "@", sep = "\t", header = F)
          
        },
        error = function(e) {
          # return a safeError if a parsing error occurs
          stop(safeError(e))
        }
      )
      
      downloading_genes()
      
    })
    
    downloading_genes <- reactive({
      input$getFiles
      isolate({
        
        withProgress(message = 'Downloading annotation', value = 0, {
          genes_table <- mapping(files_annot$bam, files_annot$sam)
        })
        
        genes_table <- genes_table[which(genes_table$external_name != "NA"),]
        genes_table$strand <- factor(genes_table$strand, levels = c("-1", "1"), 
                                     labels = c("reverse", "forward"))
        genes_table <- genes_table[with(genes_table, order(strand, start_kmer)), ]
        
        rv$data_table <- genes_table
        
      })
    })
    
    files_annot <- reactiveValues(sam = NULL, bam = NULL, kmer = NULL, cond = NULL)
    
    output$genes_list <- renderUI({
        input_genes()
    })
    
    output$region_list <- renderUI({
        input_region()
    })
    
    output$chrm_list <- renderPrint({
      cat(paste0("Chrm: ",
                 unique(rv$data_table[rv$data_table$external_name %in% input$gene & rv$data_table$region %in% input$region, "chrm"])))
    })
    
    output$strand_list <- renderPrint({
      cat(paste0("Strand: ",
                 unique(rv$data_table[rv$data_table$external_name %in% input$gene & rv$data_table$region %in% input$region, "strand"])))
    })
    
    output$placeholder <- renderPrint({ 
        tryCatch(
          {
            printSeqs(rv$data_table[rv$data_table$external_name %in% input$gene & rv$data_table$region %in% input$region,])
          },
          error = function(e){
            cat(" ")
          }
        )
    })
    
    input_genes <- reactive({
      selectInput(inputId = "gene",
                  label = "Gene:",
                  choices = c( unique(as.character(rv$data_table$external_name))))
    })
    
    input_region <- reactive({
      selectInput(inputId = "region",
                  label = "Region:",
                  choices = c( unique(as.character(rv$data_table[rv$data_table$external_name %in% input$gene, "region"]))))
    })
    
    
    # output$count_table  <- renderPrint({
    #   files_annot$kmer <- read.csv("../../output.tsv", sep = "\t")
    #   files_annot$cond <- read.csv("../../samples_cond", sep = ",", header = F)
    #   
    #   #cat(nrow(files_annot$kmer))
    #   
    # })
    
  }
), launch.browser = TRUE)