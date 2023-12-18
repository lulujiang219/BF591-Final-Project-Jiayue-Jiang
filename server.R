# Libraries
library(dplyr)
library(DT)
library(ggplot2)
library(hrbrthemes)
library(shiny)
library(pheatmap)
library(shinyWidgets)
library(msigdbr)
library(fgsea)
library('GSEABase')
library(clusterProfiler)



server <- function(input, output, session) {

  # Reactive expression to read the data
  data <- reactive({
    # Validate file type to be either tsv or csv
    req(input$uploaded_file)  # Ensure that 'input$uploaded_file' is available
    file <- input$uploaded_file
    ext <- tools::file_ext(file$datapath)
    validate(need(ext %in% c("csv", "tsv"), "Please upload a csv or tsv file"))
    
    # Read the uploaded file
    read.csv(file$datapath)
  })
  
########SAMPLES -- SUMMARY TABLE################################################
  output$summarytable <- renderTable({
    # Validate file type to be either tsv or csv
    file <- input$uploaded_file
    ext <- tools::file_ext(file$datapath)
    req(file)
    validate(need(ext %in% c("csv", "tsv"), "Please upload a csv or tsv file"))
    
    # Read the uploaded file
    data <- read.csv(file$datapath)
    
    # Create the summary table
    summary_table <- data.frame(
      "Column Names" = c("Number of genes", "Number of samples", "Condition"),
      "Types" = c("numeric", "numeric", "character"),
      "Mean (sd) or Distinct Values" = c(
        paste0(nrow(data), " genes"),
        paste0(ncol(data), " samples"),
        paste0("HD, Control")
      ),
      stringsAsFactors = FALSE
    )
    
    # Display the summary table
    summary_table
  })
  
  
########SAMPLES -- TABLE########################################################
  output$datatable <- renderDataTable({
    #validate file type to be either tsv or csv
    file <- input$uploaded_file
    ext <- tools::file_ext(file$datapath)
    req(file)
    validate(need(ext %in% c("csv", "tsv"), "Please upload a csv or tsv file"))
    
    # Read the uploaded file
    data <- read.csv(file$datapath) 
    
    
    # Display the sorted data frame
    data
  }, server = FALSE)
  
####SAMPLES -- PLOTS############################################################
  
  plot_loading <- reactiveVal(FALSE)
  
  # Inform the UI about the loading state of the plot
  output$plot_loading <- reactive({
    plot_loading()
  })
  outputOptions(output, "plot_loading", suspendWhenHidden = FALSE)
  
  
  # Observe the load_data and update choices
  observe({
    req(data())  # Ensure that data is available before updating choices
    updateSelectInput(session, "x_col", choices = names(data())[sapply(data(), is.numeric)])
  })
  
  output$density_plot <- renderPlot({
    plot_loading(TRUE)
    on.exit(plot_loading(FALSE))
    req(data())
    req(input$x_col)
    
    
    
    # Make the density plot
    p <- ggplot(data(),
           mapping = aes_string(x = input$x_col)) + 
            geom_density(fill = "#69b3a2", color = "#e9ecef", alpha = 0.8)
    p
})

################################################################################
############################COUNTS##############################################
################################################################################
  
  
  ###---------- Function to read and store the counts data ###----------
  counts_data0 <- reactive({
    req(input$counts_file)
    read.csv(input$counts_file$datapath, row.names = 1)
  })
  
  
  
  ###---------Function to filter data based on the input criteria------###
  filtered_data <- reactive({
    data <- counts_data0()
    req(data)  # Ensure data is available
    
    # Calculate gene variances
    gene_variances <- apply(data, 1, var)
    # Determine the variance threshold based on the percentile slider
    variance_threshold <- quantile(gene_variances, probs = input$variance_percentile / 100, na.rm = TRUE)
    # Filter genes based on variance percentile
    data_filtered_variance <- data[gene_variances >= variance_threshold, ]
    # Filter genes based on the number of non-zero samples
    data_filtered_final <- data_filtered_variance[rowSums(data_filtered_variance > 0) >= input$nonzero_samples, ]
    
    return(data_filtered_final)
  })
  
  
  
  ###----------Function to hold the filtering summary--------------###
  filtered_summary <- reactive({
    req(filtered_data())  # Ensure the filtered data is available
    
    total_genes <- nrow(counts_data0()) # Total number of genes before filtering
    after_filter_genes <- nrow(filtered_data()) # Number of genes after filtering
    removed_genes <- total_genes-after_filter_genes # Number of genes removed by filtering
    percentage_passed <- (after_filter_genes/total_genes)*100# Percentage of genes passing the filter
    percent_not_passed <- (removed_genes/total_genes)*100# Percentage of genes not passing the filter
    
    list(
      total_samples = ncol(counts_data0()),  # Number of samples
      total_genes = total_genes,
      filtered_genes = after_filter_genes,
      percent_filtered = percentage_passed,
      not_filtered_genes = removed_genes,
      percent_not_filtered = percent_not_passed
    )
  })
  
  
  ###---------Function to calculate statistics for the scatter plot----------###
  stats_for_plot <- reactive({
    req(counts_data0())  # Ensure the counts data is available
    
    data <- counts_data0()
    filtered <- filtered_data()
    
    # Calculate statistics
    median_count <- apply(data, 1, median)
    variance <- apply(data, 1, var, na.rm = TRUE)
    num_zeros <- rowSums(data == 0)
    
    # Create a data frame for plotting
    plot_data <- data.frame(
      MedianCount = median_count,
      Variance = variance,
      NumZeros = num_zeros,
      Filtered = rownames(data) %in% rownames(filtered)
    )
    
    plot_data
  })
  
  
  
  ###--------------------Function to plot heatmap---------------------------###
  plot_heatmap <- reactive({
    req(filtered_data())
    
    data <- filtered_data()
    
    # Apply log transformation to the data
    # Adding 1 to avoid log(0)
    log_transformed_data <- log10(data + 1)
    
    
    plot <- pheatmap(log_transformed_data, legend = TRUE)
    
    return(plot)
    
  })
  
  
  
  
  
  
  
  
  
  
  
################### COUNTS TAB --  OUTPUTS######################################

  
  output$filter_summary <- renderText({
    summary <- filtered_summary()  # Get the summary list
    paste(
      "Number of samples: ", summary$total_samples, "\n",
      "Total number of genes: ", summary$total_genes, "\n",
      "Number of genes passing current filter: ", summary$filtered_genes, "\n",
      "Percentage of genes passing current filter: ", sprintf("%.2f%%", summary$percent_filtered), "\n",
      "Number of genes not passing current filter: ", summary$not_filtered_genes, "\n",
      "Percentage of genes not passing current filter: ", sprintf("%.2f%%", summary$percent_not_filtered), "\n",
      sep = "")
  })
  
  
  
  output$scatter_plot_median_vs_variance <- renderPlot({
    plot_data <- stats_for_plot()
    req(nrow(plot_data) > 0)
    
    ggplot(plot_data, aes(x = log10(MedianCount), y = log10(Variance), color = Filtered)) +
      geom_point(alpha = 0.6) +
      scale_color_manual(values = c("TRUE" = "darkblue", "FALSE" = "lightblue")) +
      theme_minimal() +
      labs(title = "Median Count vs Variance",
           subtitle = "This plot visualizes how the variability of gene expression (variance) changes with 
the overall expression level (median count). ",
           x = "Log10(Median Count)",
           y = "Log10(Variance)")+
      theme(
        plot.title = element_text(color = "red", size = 20, face = "bold"),
        plot.subtitle = element_text(color = "blue"))
  })
  
  
  
  output$scatter_plot_median_vs_zeros <- renderPlot({
    plot_data <- stats_for_plot()
    req(nrow(plot_data) > 0)
    
    ggplot(plot_data, aes(x = log10(MedianCount), y = NumZeros, color = Filtered)) +
      geom_point(alpha = 0.6) +
      scale_color_manual(values = c("TRUE" = "darkblue", "FALSE" = "lightblue")) +
      theme_minimal() +
      labs(title = "Median Count vs Number of Zeros",
           subtitle = "This plot shows the relationship between expression level and the number of samples
in which a gene is not expressed (resulting in zero counts). ",
           x = "Log10(Median Count)",
           y = "Number of Zeros") +
      theme(
        plot.title = element_text(color = "red", size = 20, face = "bold"),
        plot.subtitle = element_text(color = "blue"))
      
  })
    
  
  
    

  output$heatmap <- renderPlot({
    result <- plot_heatmap()
    result
  })
  
  
  
  
  # Reactive value to store PCA results
  pca_res <- reactiveVal()
  
  # Observe file input and perform PCA
  observeEvent(input$counts_file, {
    req(input$counts_file)
    data <- read.csv(input$counts_file$datapath, row.names = 1)
    pca_result <- prcomp(t(data), scale. = TRUE)
    pca_res(pca_result)  # Store PCA results
    
    # Update selectInput choices
    pca_choices <- paste0("PC", 1:length(pca_result$sdev))
    updateSelectInput(session, "x_pca", choices = pca_choices, selected = pca_choices[1])
    updateSelectInput(session, "y_pca", choices = pca_choices, selected = pca_choices[2])
  })
  
  
  output$pca_scatter_plot <- renderPlot({
    
    pca_result <- pca_res()
    req(pca_result, input$x_pca, input$y_pca)
    
    xaxis_num <- as.numeric(sub("PC", "", input$x_pca))
    yaxis_num <- as.numeric(sub("PC", "", input$y_pca))
    explained_variance <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100
    
    plot_data <- data.frame(PC1 = pca_result$x[, xaxis_num], 
                            PC2 = pca_result$x[, yaxis_num])
    
    ggplot(plot_data, aes(x = PC1, y = PC2)) +
      geom_point() +
      xlab(paste(input$x_pca, "(", sprintf("%.2f%%", explained_variance[xaxis_num]), " variance)", sep = "")) +
      ylab(paste(input$y_pca, "(", sprintf("%.2f%%", explained_variance[yaxis_num]), " variance)", sep = "")) +
      theme_minimal()
  })
    
    
    
    
    
    

################################################################################
#################    DE   ######################################################
################################################################################
  
  
  ###---------Function to read the uploaded DE results----------###
  de_results <- reactive({
    req(input$de_file)
    read.csv(input$de_file$datapath)
  })
  
  
  
##########DIFFERENTIAL EXPRESSION TAB --  OUTPUTS################################
  
  # Render the DT table
  output$de_table <- renderDT({
    # Filter based on gene search
    req(de_results())
    data <- de_results()
    
    datatable(data, filter = 'top', options = list(pageLength = 10))
  }, server = TRUE)


  
  # Render the volcano plot
  output$volcano_plot <- renderPlot({
    data <- de_results()
    req(data)
    x_col <- input$volc_plot_x
    y_col <- input$volc_plot_y
    base_col <- input$volc_plot_base_col
    high_col <- input$volc_plot_high_col
    padj_threshold <- 10^input$volc_plot_padj_slider
    
    # Prepare data for the plot
    data$y_axis <- -log10(data[[y_col]])
    data$significant <- data[[y_col]] < padj_threshold
    
    ggplot(data, aes_string(x = x_col, y = "y_axis", color = "significant")) +
      geom_point(alpha = 0.5) +
      scale_color_manual(values = c("FALSE" = base_col, "TRUE" = high_col)) +
      labs(title = "Volcano Plot", x = x_col, y = paste("Negative log10 of", y_col)) +
      theme_minimal()
  })

  
  
  
  
  
  
  
##############################################################################
########Gene Expression#########################################################
###############################################################################
  # Reactive expression for reading counts matrix
    counts_data <- reactive({
        req(input$normalized_count_matrix)
        read.csv(input$normalized_count_matrix$datapath, row.names = 1)
    })

    # Reactive expression for reading sample information
    sample_info <- reactive({
        req(input$sample_info_matrix)
        read.csv(input$sample_info_matrix$datapath, row.names = 1)
    })

    # Dynamic UI for selecting a categorical variable from sample information
    output$choose_category <- renderUI({
        req(sample_info())
        selectInput("select_category", "Select Categorical Variable", choices = names(sample_info()))
    })

    # Dynamic output for displaying matching gene names
    output$gene_name_output <- renderText({
        req(counts_data())
        search_term <- input$choose_gene
        gene_names <- rownames(counts_data())
        matches <- gene_names[grepl(search_term, gene_names, ignore.case = TRUE)]
        paste(matches, collapse = ", ")
    })

    
    
    
    # Plot output based on selections and button click
    output$chosen_plot <- renderPlot({
      input$go_button
      isolate({
        counts <- counts_data()
        info <- sample_info()
        req(input$choose_gene, input$select_category, input$choose_plot)
        
        # Check if the searched gene is in the dataset
        if (!(input$choose_gene %in% rownames(counts))) {
          return(NULL)
        }
        
        # Extracting gene expression data for the selected gene
        gene_expression <- counts[input$choose_gene, , drop = FALSE]
        
        # Creating a data frame for plotting
        gene_data <- data.frame(Expression = as.numeric(gene_expression), row.names = colnames(counts))
        plot_data <- merge(gene_data, info, by = "row.names", all.x = TRUE)
        colnames(plot_data)[1] <- "Sample"
        
        # Determine the plot type and render the plot
        if(input$choose_plot == "Bar Plot") {
          ggplot(plot_data, aes_string(x = input$select_category, y = "Expression")) +
            geom_bar(stat = "identity", position = position_dodge()) +
            theme_minimal() +
            labs(title = paste("Bar Plot for", input$choose_gene), y = "Expression", x = input$select_category)
        } else if(input$choose_plot == "Boxplot") {
          ggplot(plot_data, aes_string(x = input$select_category, y = "Expression")) +
            geom_boxplot() +
            theme_minimal() +
            labs(title = paste("Boxplot for", input$choose_gene), y = "Expression", x = input$select_category)
        } else if(input$choose_plot == "Violin Plot") {
          ggplot(plot_data, aes_string(x = input$select_category, y = "Expression")) +
            geom_violin() +
            theme_minimal() +
            labs(title = paste("Violin Plot for", input$choose_gene), y = "Expression", x = input$select_category)
        } else if(input$choose_plot == "Beeswarm Plot") {
          ggplot(plot_data, aes_string(x = input$select_category, y = "Expression")) +
            geom_point(position = position_jitter(width = 0.2)) +
            theme_minimal() +
            labs(title = paste("Beeswarm Plot for", input$choose_gene), y = "Expression", x = input$select_category)
        }
      })
    })
}
