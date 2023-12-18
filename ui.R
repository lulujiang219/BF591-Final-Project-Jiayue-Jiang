#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
options(shiny.maxRequestSize = 10 * 1024^200)
library(shiny)
library(readr)
library(ggplot2)
library(shinyWidgets)
library(DT)
library(colourpicker)



ui <- fluidPage(
  
  titlePanel('BF591 final project'),
  h4("Author: Jiayue(Lulu) Jiang"),
  h5("This is an R Shiny application to examine a variety of bioinformatics processes implemented in R."),
  
  
  tabsetPanel(
    
    tabPanel( #Sample Information Tab
      "Samples",
      h5("This section allows the user to load and examine a sample information matrix before conducting analysis of corresponding sample data."),
             
             sidebarLayout(
               sidebarPanel( #Inputs
                 fileInput(inputId = 'uploaded_file',label = 'Upload a file', accept=c('.csv','.tsv'))
               ),
              
               
               mainPanel( #Outputs
                 tabsetPanel( 
                   
                   # summary table
                   tabPanel("Summary Table", tableOutput("summarytable")),
                   
                   # a data table displaying the sample information, with sortable columns
                   tabPanel("Datatable", 
                            dataTableOutput("datatable")
                            ),
                   
                   # a density plot of continuous variables.
                   tabPanel("Density Plot",
                            selectInput(inputId = "x_col", 
                                        label = "Choose a numeric column for x axis!",
                                        choices = names(data)[sapply(data, is.numeric)]),
                            plotOutput("density_plot"),
                            conditionalPanel(
                              condition = "output.plot_loading === true",
                              tags$div(class = "loader", "Loading...")  # Customize as needed
                            )
                            )
               )
             )
    )),
    
    
    

    tabPanel( #Counts Matrix Exploration
      "Counts",
      h5("This section allows the user to choose different gene filtering thresholds and assess their effects using diagnostic plots of the input counts matrix."),
             sidebarLayout(
               
               sidebarPanel( #Inputs
          #Upload Normalized counts matrix, in CSV
                 fileInput(inputId = 'counts_file',
                           label = 'Upload normalized counts matrix', 
                           accept='.csv'),
          #Slider to include genes with at least X percentile of variance
                 sliderInput(inputId = "variance_percentile",
                 label = "Select variance percentile threshold",
                 min = 0,
                 max = 100,
                 value = 50),
          #Slider to include genes with at least X samples that are non-zero
                 sliderInput(inputId = "nonzero_samples",
                             label = "Select how many gene samples to include that are non-zero",
                             min = 0,
                             max = 100,
                             value = 20)
               ),
               
               
               mainPanel( #Outputs
                   tabsetPanel(
              #Tab with text or a table summarizing the effect of the filtering
                     tabPanel("Filtering Summary",verbatimTextOutput("filter_summary")),
              #Tab with diagnostic scatter plots, where genes passing filters are marked in a darker color, and genes filtered out are lighter
                     tabPanel("Diagnostic Scatterplot", 
                              plotOutput(outputId ="scatter_plot_median_vs_variance"),
                              plotOutput("scatter_plot_median_vs_zeros")),
              #Tab with a clustered heatmap of counts remaining after filtering
                     tabPanel("Clustered Heatmap", plotOutput(outputId = "heatmap", width = "100%", height = "1600px")),
              #Tab with a scatter plot of principal component analysis projections.
                     tabPanel("PCA Visualization",
                              selectInput("x_pca", "Select the first principal component (x-axis):", choices = ""), #will update in server.R
                              selectInput("y_pca", "Select the second principal component (y-axis):", choices = ""),
                              plotOutput("pca_scatter_plot")
                              )
                              
                   
               )))),
                 
              
    
    
    tabPanel( # Differential Expression Tab
      "Differential Expression Analysis",
      h5("This section allows the user to load and explore a differential expression dataset."),
             sidebarLayout(
               
               sidebarPanel( #Inputs
                 fileInput(inputId = "de_file", 
                           label = "Please upload the result of a differential expression analysis in CSV format.",
                           accept = ".csv")
               ),
               
               mainPanel(
                 tabsetPanel(
                   tabPanel("DE Result",
                            DTOutput("de_table")
                            ),
                   tabPanel("Volcano Plot",
                            sidebarPanel(
                              radioButtons(inputId="volc_plot_x",
                                           label = "Choose the x-axis column:",
                                           choices = c("baseMean","HD.mean","Control.mean","log2FoldChange",
                                                       "lfcSE","stat","pvalue","padj"),
                                           selected = "log2FoldChange"),
                              radioButtons(inputId = "volc_plot_y",
                                           label = "Choose the y-axis column:",
                                           choices = c("baseMean","HD.mean","Control.mean","log2FoldChange",
                                                       "lfcSE","stat","pvalue","padj"),
                                           selected = "padj"),
                              colourInput("volc_plot_base_col", "Base point color:", value = "#0563fa"),
                              colourInput("volc_plot_high_col", "Highlight point color:", value = "#f2aa6b"),
                              sliderInput("volc_plot_padj_slider", "Select the magnitude of the adjust p-value coloring threshold:",
                                          min = -35, max = 0, value = -15, step = 1)
                            ),
                            mainPanel(plotOutput("volcano_plot")))
                            
                 )
               )
               
             )),
    
    
    
    
    tabPanel("Individual Gene Expression Visualization",
             h5("This section allows counts from an arbitrary gene to be selected and visualized broken out by a desired sample information variable.

"),
             sidebarLayout(
               
               sidebarPanel(fileInput(inputId = "normalized_count_matrix", label = "Upload normalized counts matrix (CSV)", accept=".csv"),
                            fileInput(inputId = "sample_info_matrix", label = "Upload sample information matrix (CSV)", accept = ".csv"),
                            uiOutput(outputId = "choose_category"),
                            textInput(inputId ="choose_gene", label = "Type a gene name", placeholder = "Type a gene name"),
                            selectInput(inputId="choose_plot", label = "select plot type", choices = c("Bar Plot", "Boxplot", "Violin Plot", "Beeswarm Plot")),
                            tags$h4("Matching Genes:"),
                            verbatimTextOutput("gene_name_output"),
                            actionButton("go_button", "Generate Plot")
                            ),
               
               mainPanel(plotOutput(outputId="chosen_plot")
                         )
               
             )
    )
    
  )
)