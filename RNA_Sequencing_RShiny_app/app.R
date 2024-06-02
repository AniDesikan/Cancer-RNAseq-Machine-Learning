#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

######################
#       IMPORTS      #
######################

library(shiny)
library(bslib)
library(DT)
library(matrixStats)
library(tidyverse)
library(pheatmap)
library(colourpicker)
library(fgsea)
library(beeswarm)
library(DESeq2)
library(data.table)
library(GSEABase)
library(rlang)

# Change this number in order to change the max upload size for files
options(shiny.maxRequestSize = 200 * 1024^2)  # Set max upload size to 100 MB

###################################################
#                 FRONT END                       #
###################################################

# This section of the code contains all of the UI for the code
# Rshiny code is written in a similar format to HTML, but it's all stored in the UI variable.

####################
# Frontend Variables
####################
# Variables that are used in the dropdowns in the UI and need to be defined beforehand

deseq_choices <-
  c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")

####################
# UI
####################

ui <- fluidPage(
  # Title 
  titlePanel("Differential Expression Pipeline"),
  helpText("Insert metadata and raw counts data into the data preprocessing tab in order to perform DEseq2, or else input DEseq2 dataframe into DE tab for analysis."),
  # Beginning of all tabs
  tabsetPanel(
    tabPanel("Data Preprocessing",
    tabsetPanel(
    ####################################################
    # SAMPLES TAB: where you look at metadata
    ####################################################
    tabPanel("Samples",
      ####################################################
      # Samples Sidebar
      ####################################################
       sidebarPanel(
         helpText("Upload a CSV file containing sample information. The file should have a header row with column names."),
         fileInput("samples_file", "Choose a CSV file for sample information exploration:"),
         
         #selectInput("numeric_columns", "Select Numeric Column", choices =c(
         #   "age_of_death","AvgSpotLen","mrna.seq_reads","pmi","RIN", "age_of_onset", "cag", "duration", "h.v_cortical_score", "h.v_striatal_score", "vonsattel_grade")
         #), # TODO: Make these choices not hardcoded, have the same UI_output as below 

         helpText("Select columns to perform DEseq2 analysis on. Make sure that the first column is the id column which has the same sample names as the RNA sequencing data. The other columns will be used as the basis for the differential expression."),
         uiOutput("samples_columns"), # dropdown menu that shows all of the columns in the csv
         
         helpText("Select a numeric column to explore. This will be only be used for the histogram plot."),
         uiOutput("numeric_columns"),
         
         submitButton(
           text = "Submit",
           icon = icon("car-crash"),
           width = "100%"
         )
         
        ),
      ####################################################
      # Samples Main Panel
      ####################################################
       mainPanel(
         tabsetPanel(
           
           # Table tab has a Dataframe output of the samples table given
           tabPanel("Table", 
                    helpText("This is the full input table. Move to the next tab to see which columns are selected."),
                    DTOutput("samples_table"),
                    downloadButton("download_samples_csv", "Download CSV"),
                    downloadButton("download_samples_xlsx", "Download Excel")),
           
           # Select columns gives the table with the columns that you selected 
           tabPanel("Selected Columns", 
                    helpText("These are the columns that are selected, to have differential expression performed on them."),
                    DTOutput("filtered_samples_table"),
                    downloadButton("download_filtered_samples_csv", "Download CSV"),
                    downloadButton("download_filtered_samples_xlsx", "Download Excel")),
           
           # Summary tab has a summary of all the selected columns, their type (int, string, etc.) and most common entry
           tabPanel("Summary", 
                    helpText("This is a summary of the full metadata table, with the names of each column, what type they are, and what the most common entry in the column is."),
                    tableOutput("samples_summary"),
                    downloadButton("download_samples_summary_csv", "Download CSV"),
                    downloadButton("download_samples_summary_xlsx", "Download Excel")),
           
           # Plots numeric column vs selected column
           # TODO: Currently the numeric columns are hardcoded so you can't select anything else
           tabPanel("Plots", 
                    helpText("This is a histogram of the numeric column selected on the sidebar."),
                    plotOutput("selected_plot"),
                    downloadButton("download_plot_jpeg", "Download JPEG"),
                    downloadButton("download_plot_png", "Download PNG"),
                    downloadButton("download_plot_svg", "Download SVG"))
         )
       )
    ),
    ####################################################
    # End of Samples Tab
    ####################################################
    
    ####################################################
    # Start of Counts Tab: Tab for filtering Counts Matrix
    ####################################################
    tabPanel("Counts", 
       ####################################################
       # Counts Sidebar
       ####################################################
       sidebarPanel(
         
         # CSV Input
         fileInput("counts_file", "Choose a CSV file for counts matrix exploration:"),
         
         # Choose how much variance you want the row to have to be selected
         helpText("Choose how much variance you want each selected gene to have."),
         sliderInput("counts_variance", "Select the percentile of variance:", min = 0, max = 100, value = 0),
         
         # Choose how many non zero samples you want the row to have to be selected
         helpText("Choose how many non-zero samples you want each selected gene to have."),
         sliderInput("counts_non_zero", "Select the number of samples that are non-zero to be included:", min = 0, max = 100, value = 0),
         
         # Submit Button
         submitButton(
           text = "Submit",
           icon = icon("car-crash"),
           width = "100%"
         )
       ),
       ####################################################
       # Counts Main Panel
       ####################################################
       mainPanel(
         tabsetPanel(
           
           # Shows the results of filtering from the sidebar
           
           # TODO: Show DTOutput of the filtered table underneath the counts filtering table
           tabPanel("Filtering Effect", tableOutput("counts_filtering_table"),
                    helpText("Shows the results of filtering the counts table with the variables from the sidebar."),
                    downloadButton("download_counts_filtering_csv", "Download CSV"),
                    downloadButton("download_counts_filtering_xlsx", "Download Excel")),
           
           # Shows how many samples actually pass the filters and where they are in terms of variance and # of zeros
           
           tabPanel("Diagnostic Scatter Plots",
                    helpText("Shows how many samples actually pass the filters and where they are in terms of variance and # of zeros."),
                    plotOutput("plot_variance"),
                    downloadButton("download_plot_variance_jpeg", "Download JPEG"),
                    downloadButton("download_plot_variance_png", "Download PNG"),
                    downloadButton("download_plot_variance_svg", "Download SVG"),
                    plotOutput("plot_zeros"),
                    downloadButton("download_plot_zeros_jpeg", "Download JPEG"),
                    downloadButton("download_plot_zeros_png", "Download PNG"),
                    downloadButton("download_plot_zeros_svg", "Download SVG")),
           
           # Shows Samples vs Genes heatmap
           
           tabPanel("Clustered Heatmap", plotOutput("counts_heatmap"),
                    helpText("Shows Samples vs Genes heatmap"),
                    downloadButton("download_heatmap_jpeg", "Download JPEG"),
                    downloadButton("download_heatmap_png", "Download PNG"),
                    downloadButton("download_heatmap_svg", "Download SVG")),
           
           # Beeswarm of which components are most important from a PCA
           
           # TODO: Add PCA as well as beeswarm
           tabPanel("PCA", 
                    sliderInput("num_components", "Choose number of components for PCA",min = 0, max = 69, value = 2),
                    helpText("Beeswarm of which components are most important from a PCA"),
                    plotOutput("counts_pca"),
                    downloadButton("download_pca_jpeg", "Download JPEG"),
                    downloadButton("download_pca_png", "Download PNG"),
                    downloadButton("download_pca_svg", "Download SVG")
         )
       )
    )
    )
    )
    ),
    
    ####################################################
    # End of Counts Tab
    ####################################################
    
    ####################################################
    # Start of DE Tab
    ####################################################
    
    # This tab will either let you insert differential expression results from a csv, or take the sample and counts
    # matrices from before and do differential expression.
    
    # TODO: Change the tab structure so Samples and Counts are under "DE Preprocessing" tab, 
    # and DE tab becomes "DE Analysis".
    
    tabPanel("DE", 
       ####################################################
       # DE Analysis Sidebar
       ####################################################
       sidebarPanel(
         # Insert DE CSV file
         # TODO: Make it so that there isn't an error when this is empty
         
         
         # Radio buttons where the user can choose whether they want to use the new data or the old data
         # TODO: Automatically use old data unless they want to use new data
         # Possibly hide file input unless no is checked on this radio button?
         helpText("If you would like to insert your own DEseq2 data, please insert click the no button and insert the file below"),
         radioButtons(
           inputId = "deseq",
           label = "Run Deseq2 on the filtered counts data from previous data?",
           choices = c("yes", "no"),
           selected = "yes"
         ),
         
         fileInput("DE_file", "Choose a CSV file for differential expression:"),
         # The following are choices for the volcano plots
         # deseq_choices are at the top 
         # TODO: Make hovering the points give you the name of the gene
         # TODO: Make clicking the points take you to the uniprot page for the gene
         
         helpText("These are options to customize the volcano plot."),
         radioButtons(
           inputId = "x_axis",
           label = "Choose the column for the x-axis",
           choices = deseq_choices,
           selected = "log2FoldChange"
         ),
         radioButtons(
           inputId = "y_axis",
           label = "Choose the column for the y-axis",
           choices = deseq_choices,
           selected = "padj"
         ),
         colourInput(
           inputId = "base",
           label = "Base point color",
           value = "#22577A",
           closeOnClick = T
         ),
         colourInput(
           inputId = "highlight",
           label = "Highlight point color",
           value = "#FFCF56",
           closeOnClick = T
         ),
         sliderInput(
           "slider",
           "Select the magnitude of the p adjusted coloring:",
           min = -50,
           max = 0,
           value = -5,
           step = 1
         ),
         
         # Submit Button
         submitButton(
           text = "Plot",
           icon = icon("car-crash"),
           width = "100%"
         )
       ),
       
       ####################################################
       # DE Analysis Main Panel
       ####################################################
       mainPanel(
         tabsetPanel(
           # Show the DE Table
           tabPanel("Sortable Table", DTOutput("DE_table"),
                    helpText("Sortable table of differential expression results. This may take a minute to load."),
                    downloadButton("download_DE_table_csv", "Download CSV"),
                    downloadButton("download_DE_table_xlsx", "Download Excel")
                    ),
           
           # Show the volcano plot and the corresponding table``
           tabPanel("DE Analysis",
                    tabsetPanel(tabPanel("Volcano Plots", 
                                         plotOutput("volcano"),
                                         downloadButton("download_volcano_jpeg", "Download JPEG"),
                                         downloadButton("download_volcano_png", "Download PNG"),
                                         downloadButton("download_volcano_svg", "Download SVG")), 
                                tabPanel("Table", 
                                         DTOutput("table"),
                                         downloadButton("download_DE_analysis_table_csv", "Download CSV"),
                                         downloadButton("download_DE_analysis_table_xlsx", "Download Excel")
                                         )
                                )
                    )
         )
       )
    ),
    ####################################################
    # End of DE Tab
    ####################################################
    
    ####################################################
    # Start of EA Tab
    ####################################################
    tabPanel("Enrichment Analysis", 
       tabsetPanel(
         ####################################################
         # EA Table
         ####################################################
         # This tab creates a table of the enrichment analysis that you can sort through
         # You can also select which NES pathways to select.
         
         tabPanel("Sortable Data Table", 
            ####################################################
            # EA Table Sidebar
            ####################################################
            sidebarPanel(
              # Insert Differential Expression CSV
              # TODO: Have tab automatically detect previous differential expression tab
              helpText("Insert the differential expression dataset to be analyzed"),
              fileInput("EA_file", "Choose a CSV file for differential expression:"),
              
              helpText("Choose which species' gmt file to use"),
              radioButtons("species", "Select Species:", 
                           c("Human" = "human", "Mouse" = "mouse")),
              
              
              helpText("Slider to select max p adjusted for NES pathways"),
              sliderInput("EA_p_adj_table", "p adjusted to select pathways by:", min =-30, max = 0, value = 0),
              
              helpText("Select whether you want negative, positive or all pathways"),
              radioButtons("EA_NES_select", "Choose which NES pathways to select:",
                           c("positive", "negative", "all")),
              
            
              submitButton(
                text = "Update Table",
                icon = icon("car-crash"),
                width = "100%"
              )
            ),
            ####################################################
            # EA Table Main Panel
            ####################################################
            mainPanel(
              DTOutput("EA_Table"),
              downloadButton("download_EA_table_csv", "Download CSV"),
              downloadButton("download_EA_table_xlsx", "Download Excel")
            )
         ),
         
         ####################################################
         # EA Top Pathways
         ####################################################
         tabPanel("Top Pathways", 
            ####################################################
            # EA Top Pathways Sidebar
            ####################################################
            sidebarPanel(
              fileInput("EA_file", "Choose a CSV file for differential expression:"),
              helpText("Adjust which pathways are selected by changing the log 10 p adjusted value."),
              sliderInput("EA_p_adj", "p adjusted to select pathways by:", min = -30, max = 0, value = -10),
              submitButton(
                text = "Update Table",
                icon = icon("car-crash"),
                width = "100%"
              )
            ),
            ####################################################
            # EA Top Pathways Main Panel
            ####################################################
            mainPanel(
              plotOutput("EA_barplot"),
              downloadButton("download_EA_barplot_jpeg", "Download JPEG"),
              downloadButton("download_EA_barplot_png", "Download PNG"),
              downloadButton("download_EA_barplot_svg", "Download SVG")
            ),
                  
         ),
         ####################################################
         # EA Scatter
         ####################################################
         tabPanel("Scatter Plot of NES", 
            ####################################################
            # EA Scatter Sidebar
            ####################################################
            sidebarPanel(
              fileInput("EA_file", "Choose a CSV file for differential expression:"),
              sliderInput("EA_p_adj_plot", "p adjusted to select pathways by:", min = -30, max = 0, value = -10),
              submitButton(
                text = "Update Table",
                icon = icon("car-crash"),
                width = "100%"
              )
            ),
            ####################################################
            # EA Scatter Main Panel
            ####################################################
            mainPanel(
              plotOutput("EA_scatter_plot"),
              helpText("This is a scatter plot showing the normalized enrichment score of different pathways vs their adjusted p values."),
              downloadButton("download_EA_scatter_plot_jpeg", "Download JPEG"),
              downloadButton("download_EA_scatter_plot_png", "Download PNG"),
              downloadButton("download_EA_scatter_plot_svg", "Download SVG")
            )
         )
       )
    )
  )
)


###################################################
#                 BACK END                        #
###################################################

# All of the functions that are used in the frontend 
server <- function(input, output) {
  
  ####################################################
  # SAMPLES FUNCTIONS
  ####################################################
  
  # Function to read the file input into "samples_file" variable of input
  load_samples_data <- reactive({
    req(input$samples_file)
    
    # Try reading the file
    new_file <- tryCatch({
      fread(input$samples_file$datapath)
    }, error = function(e) {
      NULL
    })
    
    # Return the data as a tibble, or NULL if an error occurred
    if (is.null(new_file)) {
      return(NULL)
    } else {
      return(as.tibble(new_file))
    }
  })
  
  # Function to filter the data for the input sample columns
  column_filter <- reactive({
    selected_data <- load_samples_data() %>%
      dplyr::select(input$samples_columns)
    return(selected_data)
  })
  
  # Function to create the summary for the samples tab
  samples_summary <- function(dataf) {
    newdf <- tibble(
      Columns = colnames(dataf),
      type = sapply(dataf, typeof),
      most_common = sapply(dataf, function(col) {
        if (typeof(col) %in% c("character", "factor")) {
          names(table(col))[which.max(table(col))]
        } else {
          mean(na.omit(col))
        }
      })
    )
    return(newdf)
  }
  
  plot_numeric <- function() {
    data <- load_samples_data()
    vector <- input$numeric_columns
    ggplot(data, aes(x = !!sym(vector))) +
      geom_histogram(binwidth = 10) +  # You can adjust binwidth as needed
      labs(title = paste("Histogram of", vector), x = vector, y = "Frequency")
  }
  ####################################################
  # SAMPLES OUTPUTS
  ####################################################
  
  # Display the columns of the dataframe in the dropdown so the user can select them
  output$samples_columns <- renderUI({
    data <- load_samples_data()
    choices <- colnames(data)
    selectInput("samples_columns", "Select Columns", choices, multiple = TRUE)
  })
  
  # Display the numeric columns of the dataframe in the dropdown
  output$numeric_columns <- renderUI({
    data <- load_samples_data()
    numeric_cols <- colnames(data)[sapply(data, is.numeric)]
    selectInput("numeric_columns", "Select Numeric Columns", numeric_cols, multiple = TRUE)
  })
  
  # Render the summary table
  output$samples_summary <- renderTable({
    data <- load_samples_data()
    if (is.null(data)) {
      return(NULL)
    }
    samples_summary(data)
  })
  
  # Render the samples table
  output$samples_table <- renderDT({
    data <- load_samples_data()
    if (is.null(data)) {
      return(NULL)
    }
    datatable(data, options = list(
      pageLength = 10,       # Number of rows to display per page
      lengthMenu = c(10, 25, 50, 100), # Dropdown menu for number of rows to display
      ordering = TRUE
    ))
  })
  
  # Render the filtered samples table based on selected columns
  output$filtered_samples_table <- renderDT({
    data <- load_samples_data()
    if (is.null(data) || is.null(input$samples_columns)) {
      return(NULL)
    }
    selected_data <-dplyr::select(data, input$samples_columns)
    datatable(column_filter(), options = list(
      pageLength = 10,       # Number of rows to display per page
      lengthMenu = c(10, 25, 50, 100), # Dropdown menu for number of rows to display
      ordering = TRUE
    ))
  })
  
  # Render the plot for the selected numeric column
  output$selected_plot <- renderPlot({
    data <- load_samples_data()
    if (is.null(data) || is.null(input$numeric_columns)) {
      return(NULL)
    }
    vector <- input$numeric_columns
    hist(data[[vector]], main = paste("Histogram of", vector))
  })
  
  ####################################################
  # SAMPLES DOWNLOADS
  ####################################################
  
  output$download_samples_csv <- downloadHandler(
    filename = function() {
      paste("samples-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(load_samples_data(), file)
    }
  )
  
  output$download_samples_xlsx <- downloadHandler(
    filename = function() {
      paste("samples-", Sys.Date(), ".xlsx", sep = "")
    },
    content = function(file) {
      writexl::write_xlsx(load_samples_data(), file)
    }
  )
  
  output$download_filtered_samples_csv <- downloadHandler(
    filename = function() {
      paste("filtered-samples-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(column_filter(), file)
    }
  )
  
  output$download_filtered_samples_xlsx <- downloadHandler(
    filename = function() {
      paste("filtered-samples-", Sys.Date(), ".xlsx", sep = "")
    },
    content = function(file) {
      writexl::write_xlsx(column_filter(), file)
    }
  )
  
  output$download_samples_summary_csv <- downloadHandler(
    filename = function() {
      paste("samples-summary-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(samples_summary(load_samples_data()), file)
    }
  )
  
  output$download_samples_summary_xlsx <- downloadHandler(
    filename = function() {
      paste("samples-summary-", Sys.Date(), ".xlsx", sep = "")
    },
    content = function(file) {
      writexl::write_xlsx(samples_summary(load_samples_data()), file)
    }
  )
  
  # Download handlers for plots
  output$download_plot_jpeg <- downloadHandler(
    filename = function() {
      paste("plot-", Sys.Date(), ".jpeg", sep = "")
    },
    content = function(file) {
      
      ggsave(file, plot = plot_numeric(), device = "jpeg")
    }
  )
  
  output$download_plot_png <- downloadHandler(
    filename = function() {
      paste("plot-", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = plot_numeric(), device = "png")
    }
  )
  
  output$download_plot_svg <- downloadHandler(
    filename = function() {
      paste("plot-", Sys.Date(), ".svg", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = plot_numeric(), device = "svg")
    }
  )

  ####################################################
  # COUNTS FUNCTIONS
  ####################################################
  
  load_counts_data <- reactive({
    # Read the data using read.table
    new_file <- fread(input$counts_file$datapath)
    return(as.tibble(new_file))
  })
  
  counts_filtering <- function(variance, non_zero) {
    dataf <- load_counts_data()
    
    row_variances <- rowVars(as.matrix(dataf[,-1]), na.rm = TRUE)
    
    # Filter the data based on variance threshold
    filtered_counts <- dplyr::filter(dataf, row_variances >= variance)
    
    # Filter out genes with too few non-zero samples
    filtered_counts <- dplyr::filter(filtered_counts, rowSums(filtered_counts != 0) >= non_zero)
    
    return(filtered_counts)
  }
  
  # Make a table that summarizes the effect of filtering
  # Return a table with the number of samples, number of genes, number and % of genes passing filter and not passing filter
  counts_filtering_table <- function(dataf, variance, zero) {
    filtered_counts <- counts_filtering(variance, zero)
    
    # Create the summary table
    table <- tibble(
      `Number of samples` = ncol(dataf),
      `Total Number of Genes` = nrow(dataf),
      `Number of Genes passing the filter` = nrow(filtered_counts),
      `Percent of Genes passing the filter` = nrow(filtered_counts)/nrow(dataf) * 100,
      `Number of Genes not passing filter` = nrow(dataf) - nrow(filtered_counts),
      `Percent of Genes not passing the filter` =  ((nrow(dataf) - nrow(filtered_counts))/nrow(dataf))*100
    )
    
    # Return the summary table
    return(table)
  }
  
  #Make diagnostic scatter plots
  # Genes passing in a darker color, genes not passing in a lighter color
  # Do median count vs variance and median count vs number of zeros
  counts_scatter_plot <- function(variance, non_zero) {
    dataf <- load_counts_data()
    filtered_counts <- counts_filtering(variance, non_zero)
    df_var <- rowVars(as.matrix(dataf[, -1]), na.rm = TRUE)
    df_zeros <- rowSums(dataf == 0)
    
    diagnostic_data <- cbind(dataf, df_var, df_zeros)
    diagnostic_data$filtered <- ifelse(df_var >= variance & df_zeros>=non_zero, "Pass", "Fail")
    
    # Create scatter plot for median count vs variance
    plot_variance <- ggplot(diagnostic_data, aes(x = rank(rowMeans(dataf[, -1], na.rm = TRUE)), y = log2(df_var))) +
      geom_point(aes(color = filtered), alpha = 0.7) +
      labs(title = "Median Count vs Variance",
           x = "rank Median Count",
           y = "Log2 Variance") +
      theme_minimal()
    
    # Create scatter plot for median count vs number of zeros
    plot_zeros <- ggplot(diagnostic_data, aes(x = rank(rowMeans(dataf[, -1], na.rm = TRUE)), y = df_zeros)) +
      geom_point(aes(color = filtered), alpha = 0.7) +
      labs(title = "Median Count vs Number of Zeros",
           x = "rank Median Count",
           y = "Number of Zeros") +
      theme_minimal()
    
    # Return a list of plots
    return(list(plot_variance, plot_zeros))
  }
  
  #Make clustered Heatmap
  # enable log transforming counts
  # make sure there's a legend with color bar
  counts_heatmap <-
    function(variance, non_zero) {
      filtered_counts <- counts_filtering(variance, non_zero)
      filtered_counts <- filtered_counts[,-1]
      
      # print(dim(filtered_counts))
      # print(head(filtered_counts))
      log_transformed_counts <- log2(filtered_counts + 1)
      
      clustered_data <- as.matrix(log_transformed_counts)
      row_order <- order(rowMeans(clustered_data, na.rm = TRUE))
      col_order <- order(colMeans(clustered_data, na.rm = TRUE))
      
      # Create the clustered heatmap
      heatmap(
        clustered_data[row_order, col_order],
        Colv = FALSE,
        Rowv = FALSE,
        scale = "row",  # Scale rows (genes)
        col = colorRampPalette(c("white", "blue"))(100),
        margins = c(5, 10),
        xlab = "Samples",
        ylab = "Genes"
      )
      
      
      # Return the heatmap object
      return(heatmap)
    }
  
  counts_pca <-
    function(variance, non_zero, num_components) {
      filtered_counts <- counts_filtering(variance, non_zero)
      df <- as.data.frame(filtered_counts[, -1])  # Assuming the first column is sample/gene names
      
      # Transpose the data for PCA
      expr_mat <- prcomp(t(df), center = TRUE)
      num <- num_components
      # Create a data frame with PC1 and PC2 scores
      # pca_data <- data.frame(PC1 = expr_mat$x[, 1], PC2 = expr_mat$x[, 2])
      pca_data <- as.data.frame(expr_mat$x[, 1:num])
      
      # Create a PCA object without centering and scaling
      pca <- prcomp(
        pca_data,
        center = FALSE,
        scale = FALSE
      )
      
      # Create a summary table
      summary_table <- summary(pca)
      
      component_colors <- rainbow(num_components)
      # Create a scatter plot
      # final <- ggplot(pca_data, aes(x = PC1, y = PC2,)) +
      #   geom_point() +
      #   labs(
      #     title = "PCA Scatter Plot",
      #     x = paste("PC1 (", round(summary_table$sdev[1] / sum(summary_table$sdev) * 100, 2), "%)"),
      #     y = paste("PC2 (", round(summary_table$sdev[2] / sum(summary_table$sdev) * 100, 2), "%)")
      #   )+
      #   theme_minimal()
      
      # final <- ggplot(NULL, aes(x = 1:length(summary_table$sdev), y = summary_table$sdev^2 / sum(summary_table$sdev^2))) +
      #   geom_bar(stat = "identity", fill = "skyblue", width = 0.7) +
      #   labs(
      #     title = "Scree Plot",
      #     x = "Principal Components",
      #     y = "Proportion of Variance"
      #   ) +
      #   theme_minimal()
      
      final <- beeswarm(pca_data, pch = 16, col = component_colors, main = "Beeswarm Plot of Principal Components")
      
      legend("bottomright", legend = paste("PC", 1:num_components, ": ", round(summary(expr_mat)$sdev[1:num_components] / sum(summary(expr_mat)$sdev) * 100, 2), "%"), col = component_colors, pch = 16)
      
      return(final)
    }
  
  ####################################################
  # COUNTS OUTPUTS
  ####################################################
  
  output$counts <- renderDT(
    datatable(load_counts_data(), options = list(ordering = TRUE))
  )
  
  output$counts_filtering_table <- renderTable({
    req(input$counts_file)
    data <- load_counts_data()
    counts_filtering_table(data, input$counts_variance, input$counts_non_zero)
  })
  
  # Display the plots
  output$plot_variance <- renderPlot({
    req(input$counts_file)
    plots <- counts_scatter_plot(input$counts_variance, input$counts_non_zero)
    print(plots[[1]])
  })
  
  output$plot_zeros <- renderPlot({
    req(input$counts_file)
    plots <- counts_scatter_plot(input$counts_variance, input$counts_non_zero)
    print(plots[[2]])
  })
  
  output$counts_heatmap<- renderPlot({
    req(input$counts_file)
    counts_heatmap(input$counts_variance, input$counts_non_zero)
  })
  
  output$counts_pca<- renderPlot({
    req(input$counts_file)
    counts_pca(input$counts_variance, input$counts_non_zero, input$num_components)
  })
  
  ####################################################
  # COUNTS DOWNLOADS
  ####################################################

  output$download_counts_filtering_csv <- downloadHandler(
    filename = function() {
      paste("counts-filtered-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(counts_filtering(data, input$counts_variance, input$counts_non_zero), file)
    }
  )
  
  output$download_counts_filtering_xlsx <- downloadHandler(
    filename = function() {
      paste("counts-filtered-", Sys.Date(), ".xlsx", sep = "")
    },
    content = function(file) {
      writexl::write_xlsx(counts_filtering(data, input$counts_variance, input$counts_non_zero), file)
    }
  )

   output$download_plot_variance_jpeg <- downloadHandler(
     filename = function() {
       paste("plot-variance-", Sys.Date(), ".jpeg", sep = "")
     },
     content = function(file) {
       ggsave(file, plot = counts_scatter_plot(input$counts_variance, input$counts_non_zero)[[1]], device = "jpeg")
     }
   )
   
   output$download_plot_variance_png <- downloadHandler(
     filename = function() {
       paste("plot-variance-", Sys.Date(), ".png", sep = "")
     },
     content = function(file) {
       ggsave(file, plot = counts_scatter_plot(input$counts_variance, input$counts_non_zero)[[1]], device = "png")
     }
   )
   
   output$download_plot_variance_svg <- downloadHandler(
     filename = function() {
       paste("plot-variance-", Sys.Date(), ".svg", sep = "")
     },
     content = function(file) {
       ggsave(file, plot = counts_scatter_plot(input$counts_variance, input$counts_non_zero)[[1]], device = "svg")
     }
   )


  output$download_plot_zeros_jpeg <- downloadHandler(
    filename = function() {
      paste("plot-zeros-", Sys.Date(), ".jpeg", sep = "")
    },
    content = function(file) {
      
      ggsave(file, plot = plot_numeric(), device = "jpeg")
    }
  )
  
  output$download_plot_zeros_png <- downloadHandler(
    filename = function() {
      paste("plot-zeros-", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = plot_numeric(), device = "png")
    }
  )
  
  output$download_plot_zeros_svg <- downloadHandler(
    filename = function() {
      paste("plot-zeros-", Sys.Date(), ".svg", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = plot_numeric(), device = "svg")
    }
  )

  output$download_heatmap_jpeg <- downloadHandler(
    filename = function() {
      paste("heatmap-", Sys.Date(), ".jpeg", sep = "")
    },
    content = function(file) {
      
      ggsave(file, plot = plot_numeric(), device = "jpeg")
    }
  )
  
  output$download_heatmap_png <- downloadHandler(
    filename = function() {
      paste("heatmap-", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = plot_numeric(), device = "png")
    }
  )
  
  output$download_heatmap_svg <- downloadHandler(
    filename = function() {
      paste("heatmap-", Sys.Date(), ".svg", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = plot_numeric(), device = "svg")
    }
  )

  
  output$download_pca_jpeg <- downloadHandler(
    filename = function() {
      paste("pca-", Sys.Date(), ".jpeg", sep = "")
    },
    content = function(file) {
      
      ggsave(file, plot = plot_numeric(), device = "jpeg")
    }
  )
  
  output$download_pca_png <- downloadHandler(
    filename = function() {
      paste("pca-", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = plot_numeric(), device = "png")
    }
  )
  
  output$download_pca_svg <- downloadHandler(
    filename = function() {
      paste("pca-", Sys.Date(), ".svg", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = plot_numeric(), device = "svg")
    }
  )
  
  ####################################################
  # DE FUNCTIONS
  ####################################################
  
  load_DE_data <- reactive({
    # Read the data using read.table
    if(input$deseq == "yes"){
      new_file <- run_deseq(counts_filtering(input$counts_variance, input$counts_non_zero))
    }else{
      new_file <- fread(input$DE_file$datapath)
    }
    return(as.tibble(new_file))
  })
  
  
  run_deseq <- function(counts_data) {
    print("working")
    
    # turn the first column of the counts matrix into rownames
    counts_row_names <- counts_data[[1]]
    print(counts_data)
    counts_data <- as.data.frame(counts_data[,-1])
    print(counts_data)
    print(length(counts_row_names))
    print(nrow(counts_data))
    rownames(counts_data) <- counts_row_names
    
    # Make sure that the ids in the metadata match the ids in the counts data
    # Load and process the total_design
    total_design <- as.data.frame(column_filter())
    sample_ids <- colnames(counts_data) 
    print(sample_ids)
    filtered_design <- total_design[total_design[[1]] %in% sample_ids, ]
    
    # Make the ids of the metadata the rownames
    metadata_row_names <- filtered_design[, 1]
    metadata_col_names<- colnames(filtered_design)[-1]
    filtered_design <- as.data.frame(filtered_design[, -1])
    colnames(filtered_design) <- metadata_col_names
    rownames(filtered_design) <- metadata_row_names
    print(filtered_design)
    
    # Print the row names
    print(rownames(filtered_design))
    print(colnames(counts_data))
    
    # Check if counts_data has row names and column names
    if (is.null(rownames(counts_data))) {
      stop("countData must have row names (gene identifiers)")
    }
    
    if (is.null(colnames(counts_data))) {
      stop("countData must have column names (sample identifiers)")
    }
    
    # Ensure the sample names in counts_data match row names in total_design
    #if (!all(colnames(counts_data) %in% rownames(total_design))) {
    #  stop("Column names in counts matrix must match sample_id column in metadata!")
    #}
    
    design_formula <- as.formula(paste("~", metadata_col_names))
    print(design_formula)
    # Create DESeqDataSet
    dds <- DESeqDataSetFromMatrix(countData = counts_data,
                                  colData = filtered_design,  # Sample information
                                  design = design_formula)    # Design formula with all columns in total_design
    
    # Run DESeq2 analysis
    dds <- DESeq(dds)
    
    # Get results
    res <- results(dds)
    print(res)
    res$Gene <- counts_row_names
    res <- res[, c("Gene", setdiff(names(res), "Gene"))]
    
    return(res)
  }
  
  #Volcano Plot Function
  volcano_plot <-
    function(dataf, x_name, y_name, slider, color1, color2) {
      p <- ggplot(dataf, aes(x = !!sym(x_name),
                             y = -log10(!!sym(y_name)))) +
        geom_point(aes(color = !!sym(y_name) < 1 * 10 ^ (as.numeric(slider)))) +
        theme_bw() +
        scale_color_manual(values = c(color1, color2)) +
        theme(legend.position = "bottom") +
        labs(color = paste0(y_name, " < 1 Ã— 10^", slider))
      return(p)
    }
  
  draw_table <- function(dataf, slider) {
    if (!"padj" %in% colnames(dataf)) {
      stop("Column 'padj' not found in the data frame.")
    }
    new_table <- dplyr::filter(dataf, dataf$padj < 10^(slider))
    new_table$padj <- prettyNum(new_table$padj, digits = 5)
    return(new_table)
  }
  
  ####################################################
  # DE OUTPUTS
  ####################################################
  
  output$DE_table <- renderDT(
    datatable(load_DE_data(), options = list(ordering = TRUE))
  )
  
  output$volcano <- renderPlot({
    p <-volcano_plot(load_DE_data(),
                     input$x_axis,
                     input$y_axis,
                     input$slider,
                     input$base,
                     input$highlight)
    return(p)
  })
  
  output$table <- renderDT(
    datatable(draw_table(load_DE_data(),input$slider))
  )
  
  ####################################################
  # DE DOWNLOADS
  ####################################################
  
  
  output$download_DE_table_csv <- downloadHandler(
    filename = function() {
      paste("differential_expression-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(load_DE_data(), file)
    }
  )
  
  output$download_DE_table_xlsx <- downloadHandler(
    filename = function() {
      paste("differential-expression-", Sys.Date(), ".xlsx", sep = "")
    },
    content = function(file) {
      writexl::write_xlsx(load_DE_data(), file)
    }
  )


  
output$download_DE_analysis_csv <- downloadHandler(
  filename = function() {
    paste("DE-analysis-", Sys.Date(), ".csv", sep = "")
  },
  content = function(file) {
    write.csv(draw_table(load_DE_data(),input$slider), file)
  }
)

output$download_DE_analysis_xlsx <- downloadHandler(
  filename = function() {
    paste("DE-analysis-", Sys.Date(), ".xlsx", sep = "")
  },
  content = function(file) {
    writexl::write_xlsx(draw_table(load_DE_data(),input$slider), file)
  }
)            
                       
  ####################################################
  # EA FUNCTIONS
  ####################################################
  
  load_EA_data <- reactive({
    req(input$EA_file)
    new_file <- read.csv(input$EA_file$datapath)
    print(head(new_file))
    return(as_tibble(new_file))
  })
  
  make_ranked_log2fc <- function(dataf) {
    new_table <- dataf %>%
      mutate(rank = rank(log2FoldChange)) %>%
      arrange(rank) %>%
      dplyr::select(Gene, log2FoldChange) %>%
      filter(is.finite(log2FoldChange))
    new_vector <- pull(new_table, log2FoldChange, Gene)
    return(new_vector)
  }
  
  run_fgsea <- reactive({
    req(input$species)
    species <- input$species
    gmt_pathways <- if (species == "human") {
      gmtPathways("~/Documents/GitHub/591_files/R_Shiny_app/rnaseq_rewritten/gmt_files/h.all.v2023.2.Hs.symbols.gmt.txt")
    } else if (species == "mouse") {
      gmtPathways("~/Documents/GitHub/591_files/R_Shiny_app/rnaseq_rewritten/gmt_files/mh.all.v2023.2.Mm.symbols.gmt.txt")
    }
    
    rnk_list <- make_ranked_log2fc(load_EA_data())
    print(head(rnk_list))
    fgsea_results <- fgsea(pathways = gmt_pathways, stats = rnk_list, minSize = 15, maxSize = 500) %>%
      as_tibble()
    print(head(fgsea_results))
    return(fgsea_results)
  })
  
  fgsea_filter <- reactive ({
    ea_table <- run_fgsea
    if(select == "positive"){
      ea_table <- dplyr::filter(ea_table, NES > 0)
    } 
    if(select == "negative"){
      ea_table <- dplyr::filter(ea_table, NES < 0)
    } 
    ea_table <- dplyr::filter(ea_table, padj < 10^(p_adj))
    return(ea_table)
  })
  
  EA_barplot <- function(fgsea_results, slider){
    
    top_pos <- fgsea_results %>% filter((padj) < 10^slider) %>% pull(pathway)
    
    subset <- fgsea_results %>% 
      filter(pathway %in% c(top_pos)) %>%
      mutate(pathway = factor(pathway)) %>%
      mutate(plot_name = str_replace_all(pathway, '_', ' '))
    
    plot <- subset %>% 
      mutate(plot_name = forcats::fct_reorder(factor(plot_name), NES)) %>%
      ggplot() +
      geom_bar(aes(x=plot_name, y=NES, fill = NES > 0), stat='identity', show.legend = FALSE) +
      scale_fill_manual(values = c('TRUE' = 'red', 'FALSE' = 'blue')) + 
      theme_minimal(base_size = 8) +
      ggtitle('fgsea results for Hallmark MSigDB gene sets') +
      ylab('Normalized Enrichment Score (NES)') +
      xlab('') +
      scale_x_discrete(labels = function(x) str_wrap(x, width = 80)) +
      coord_flip()
    return(plot)
  }
  
  EA_Table <-
    function(p_adj, select) {
      ea_table <- run_fgsea()
      if(select == "positive"){
        ea_table <- dplyr::filter(ea_table, NES > 0)
      } 
      if(select == "negative"){
        ea_table <- dplyr::filter(ea_table, NES < 0)
      } 
      ea_table <- dplyr::filter(ea_table, padj < 10^(p_adj))
      return(ea_table)
    }
  
  scatter_plot_data <- reactive({
    ea_table <- run_fgsea()
    
    # Filter gene sets below p-value threshold
    ea_table <- dplyr::filter(ea_table, padj < 10^(-input$EA_p_adj_plot))
    
    return(ea_table)
  })
  
  
  ####################################################
  # EA OUTPUTS
  ####################################################
  
  output$EA_barplot  <- renderPlot(
    EA_barplot(run_fgsea(), input$EA_p_adj)
  )
  
  output$EA_Table  <- renderDT(
    datatable(EA_Table(input$EA_p_adj_table, input$EA_NES_select))
  )
  
  output$EA_scatter_plot <- renderPlot({
    scatter_data <- scatter_plot_data()
    
    # Create scatter plot using ggplot2
    ggplot(scatter_data, aes(x = NES, y = -log10(padj), color = factor(padj < 10^(input$EA_p_adj_plot)))) +
      geom_point(alpha = 0.7) +
      scale_color_manual(values = c("TRUE" = "grey", "FALSE" = "blue")) +
      labs(
        title = "Scatter Plot of NES vs. -log10 Adjusted p-value",
        x = "Normalized Enrichment Score (NES)",
        y = "-log10 Adjusted p-value"
      )
  })
  
  ####################################################
  # EA DOWNLOADS
  ####################################################

  output$download_EA_table_csv <- downloadHandler(
    filename = function() {
      paste("enrichment_analysis-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(EA_Table(input$EA_p_adj_table, input$EA_NES_select), file)
    }
  )
  
  output$download_EA_table_xlsx <- downloadHandler(
    filename = function() {
      paste("enrichment_analysis-", Sys.Date(), ".xlsx", sep = "")
    },
    content = function(file) {
      writexl::write_xlsx(EA_Table(input$EA_p_adj_table, input$EA_NES_select), file)
    }
  )
  
  output$download_EA_barplot_jpeg <- downloadHandler(
    filename = function() {
      paste("EA_barplot-", Sys.Date(), ".jpeg", sep = "")
    },
    content = function(file) {
      ggsave("EA_barplot.jpeg", plot = EA_barplot(run_fgsea(), input$EA_p_adj), device = "jpeg")
    }
  )
  
  output$download_EA_barplot_png <- downloadHandler(
    filename = function() {
      paste("EA_barplot-", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      ggsave("EA_barplot.png", plot = EA_barplot(run_fgsea(), input$EA_p_adj), device = "png")
    }
  )
  
  output$download_EA_barplot_svg <- downloadHandler(
    filename = function() {
      paste("EA_barplot-", Sys.Date(), ".svg", sep = "")
    },
    content = function(file) {
      ggsave("EA_barplot.svg", plot = EA_barplot(run_fgsea(), input$EA_p_adj), device = "svg")
    }
  )
  
  output$download_EA_scatter_plot_jpeg <- downloadHandler(
    filename = function() {
      paste("EA_scatter_plot-", Sys.Date(), ".jpeg", sep = "")
    },
    content = function(file) {
      ggsave("EA_scatter_plot.jpeg", plot = ggplot(scatter_plot_data(), aes(x = NES, y = -log10(padj), color = factor(padj < 10^(input$EA_p_adj_plot)))) + 
               geom_point(alpha = 0.7) +
               scale_color_manual(values = c("TRUE" = "grey", "FALSE" = "blue")) +
               labs(
                 title = "Scatter Plot of NES vs. -log10 Adjusted p-value",
                 x = "Normalized Enrichment Score (NES)",
                 y = "-log10 Adjusted p-value"
               ), device = "jpeg")
    }
  )
  
  output$download_EA_scatter_plot_png <- downloadHandler(
    filename = function() {
      paste("EA_scatter_plot-", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      ggsave("EA_scatter_plot.png", plot = ggplot(scatter_plot_data(), aes(x = NES, y = -log10(padj), color = factor(padj < 10^(input$EA_p_adj_plot)))) + 
               geom_point(alpha = 0.7) +
               scale_color_manual(values = c("TRUE" = "grey", "FALSE" = "blue")) +
               labs(
                 title = "Scatter Plot of NES vs. -log10 Adjusted p-value",
                 x = "Normalized Enrichment Score (NES)",
                 y = "-log10 Adjusted p-value"
               ), device = "png")
    }
  )
  
  output$download_EA_scatter_plot_svg <- downloadHandler(
    filename = function() {
      paste("EA_scatter_plot-", Sys.Date(), ".svg", sep = "")
    },
    content = function(file) {
      ggsave("EA_scatter_plot.svg", plot = ggplot(scatter_plot_data(), aes(x = NES, y = -log10(padj), color = factor(padj < 10^(input$EA_p_adj_plot)))) + 
               geom_point(alpha = 0.7) +
               scale_color_manual(values = c("TRUE" = "grey", "FALSE" = "blue")) +
               labs(
                 title = "Scatter Plot of NES vs. -log10 Adjusted p-value",
                 x = "Normalized Enrichment Score (NES)",
                 y = "-log10 Adjusted p-value"
               ), device = "svg")
    }
  )
}


###################################################
#                 RUN APP                         #
###################################################


# Run the application 
shinyApp(ui = ui, server = server)

