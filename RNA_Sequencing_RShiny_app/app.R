#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(DT)
library(matrixStats)
library(tidyverse)
library(pheatmap)
library(colourpicker)
library(fgsea)
library(beeswarm)
library(DESeq2)
library(data.table)
options(shiny.maxRequestSize = 200 * 1024^2)  # Set max upload size to 100 MB

deseq_choices <-
  c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")

ui <- fluidPage(
  titlePanel("Plotting of Normalized Counts Data"),
  
  tabsetPanel(
    tabPanel("Samples",
      sidebarPanel(
        fileInput("samples_file", "Choose a CSV file for sample information exploration:"),
        selectInput("numeric_columns", "Select Numeric Column", choices =c(
        "age_of_death","AvgSpotLen","mrna.seq_reads","pmi","RIN", "age_of_onset", "cag", "duration", "h.v_cortical_score", "h.v_striatal_score", "vonsattel_grade")
        ),
        uiOutput("samples_columns"),
        submitButton(
          text = "Submit",
          icon = icon("car-crash"),
          width = "100%"
        )
      )
      ,
      mainPanel(
        tabsetPanel(
          tabPanel("Table", 
                   DTOutput("samples_table")),
          tabPanel("Select Columns", 
                   DTOutput("filtered_samples_table")),
          tabPanel("Summary", 
                   tableOutput("samples_summary")),
          tabPanel("Plots", 
                   plotOutput("selected_plot"))
        )
      )
    ),
    
    tabPanel("Counts", 
      sidebarPanel(
        fileInput("counts_file", "Choose a CSV file for counts matrix exploration:"),
        sliderInput("counts_variance", "Select the percentile of variance:", min = 0, max = 100, value = 0),
        sliderInput("counts_non_zero", "Select the number of samples that are non-zero to be included:", min = 0, max = 100, value = 0),
        
        submitButton(
          text = "Submit",
          icon = icon("car-crash"),
          width = "100%"
        )
      ),
      mainPanel(
        tabsetPanel(
          tabPanel("Filtering Effect", tableOutput("counts_filtering_table")),
          tabPanel("Diagnostic Scatter Plots", plotOutput("plot_variance"), plotOutput("plot_zeros")),
          tabPanel("Clustered Heatmap", plotOutput("counts_heatmap")),
          tabPanel("PCA", sliderInput("num_components", "Choose number of components for PCA",min = 0, max = 69, value = 2 ),
                   plotOutput("counts_pca"))
        )
      )
    ),
    
    tabPanel("DE", 
      sidebarPanel(
        fileInput("DE_file", "Choose a CSV file for differential expression:"),
        radioButtons(
          inputId = "deseq",
          label = "Run Deseq2 on the filtered counts data from previous data?",
          choices = c("yes", "no"),
          selected = "no"
        ),
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
        submitButton(
          text = "Plot",
          icon = icon("car-crash"),
          width = "100%"
        )
      ),
      mainPanel(
        tabsetPanel(
          tabPanel("Sortable Table", DTOutput("DE_table")),
          tabPanel("DE Analysis",
           tabsetPanel(tabPanel("Volcano Plots", plotOutput("volcano")), tabPanel("Table", DTOutput("table"))))
          #   tabPanel("Volcano Plots", plotOutput("volcano")),
          #   tabPanel("Table", tableOutput("table"))
        )
      )
    ),
    
    tabPanel("Enrichment Analysis", 
      tabsetPanel(
        tabPanel("Sortable Data Table", 
                 sidebarPanel(
                   fileInput("EA_file", "Choose a CSV file for differential expression:"),
                   sliderInput("EA_p_adj_table", "p adjusted to select pathways by:", min =-30, max = 0, value = 0),
                   radioButtons("EA_NES_select", "Choose which NES pathways to select:",
                                c("positive", "negative", "all")),
                   submitButton(
                     text = "Update Table",
                     icon = icon("car-crash"),
                     width = "100%"
                   ),
                   p("Download the current filtered table:"),
                   downloadButton("EA_download", "Download Table")
                   
                 ),
                 mainPanel(
                   DTOutput("EA_Table")
                 )
        ),
        tabPanel("Top Pathways", 
          sidebarPanel(
            fileInput("EA_file", "Choose a CSV file for differential expression:"),
            sliderInput("EA_p_adj", "p adjusted to select pathways by:", min = -30, max = 0, value = -10),
            submitButton(
              text = "Update Table",
              icon = icon("car-crash"),
              width = "100%"
            )
          ),
          mainPanel(
            plotOutput("EA_barplot")
          ),
          
        ),
        
        tabPanel("Scatter Plot of NES", 
          sidebarPanel(
            fileInput("EA_file", "Choose a CSV file for differential expression:"),
            sliderInput("EA_p_adj_plot", "p adjusted to select pathways by:", min = -30, max = 0, value = -10),
            submitButton(
              text = "Update Table",
              icon = icon("car-crash"),
              width = "100%"
            )
          ),
          mainPanel(
            plotOutput("EA_scatter_plot")
          )
        )
      )
    )
  )
)






server <- function(input, output) {
  
  #HELPER FUNCTIONS:
  
  #Filter sample columns
  
  column_filter <- reactive({
    selected_data <- load_samples_data() %>%
      select(input$samples_columns)
    return(selected_data)
  })
  
  #deseq function to run deseq on normalized data if deseq data isn't available
  run_deseq <- function(counts_data) {
    total_design <- column_filter()
    filter <- total_design[,2]
    print(counts_data)
    # Assuming counts_data has columns with sample names and rows with gene counts
    # colnames(counts_data) <- make.names(colnames(counts_data))  # Ensure valid column names
    # Create DESeqDataSet
    dds <- DESeqDataSetFromMatrix(countData = counts_data,
                                  colData = colData,  # You might want to add sample information here
                                  design = ~1)    # Simple design, adjust as needed
    # Run DESeq2 analysis
    dds <- DESeq(dds)
    # Get results
    res <- results(dds)
    return(res)
  }
  
  # Filter counts
  counts_filtering <- function(variance, zero) {
    dataf <- load_counts_data()
    
    row_variances <- rowVars(as.matrix(dataf[,-1]), na.rm = TRUE)
    
    # Filter the data based on variance threshold
    filtered_counts <- dplyr::filter(dataf, row_variances >= variance)
    filtered_counts <- dplyr::filter(filtered_counts, rowSums(filtered_counts == 0) >= zero)
    return(filtered_counts)
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
  
  #Table drawing function
  draw_table <- function(dataf, slider) {
    if (!"padj" %in% colnames(dataf)) {
      stop("Column 'padj' not found in the data frame.")
    }
    new_table <- dplyr::filter(dataf, dataf$padj < 10^(slider))
    new_table$padj <- prettyNum(new_table$padj, digits = 5)
    return(new_table)
  }
  
  # Function to load data for each tab, will see if it's possible to make all of this one function later
  
  # load_samples_data <- reactive({
  #   new_file <- read.csv(input$samples_file$datapath)
  #   return(new_file)
  # })
  load_samples_data <- reactive({
    # con <- gzfile(input$samples_file$datapath, "rt")
    # Read the data using read.table
    # new_file <- read.table(con, header = TRUE, sep = "\t", fill = TRUE)
    new_file <- fread(input$samples_file$datapath)
    return(as.tibble(new_file))
  })
  load_counts_data <- reactive({
    # Read the data using read.table
    new_file <- fread(input$counts_file$datapath)
    return(as.tibble(new_file))
  })
  load_DE_data <- reactive({
    # Read the data using read.table
    
    if(input$deseq == "yes"){
      new_file <- run_deseq(counts_filtering(input$counts_variance, input$counts_non_zero))
    }else{
      new_file <- fread(input$DE_file$datapath)
    }
    
    return(as.tibble(new_file))
  })
  load_EA_data <- reactive({
    con <- gzfile(input$EA_file$datapath, "rt")
    # Read the data using read.table
    new_file <- read.table(con, header = TRUE, sep = "\t")
    print(new_file)
    return(as.tibble(new_file))
  })
  
  
  output$samples_columns <- renderUI({
    # Get the column names of the data frame
    choices <- colnames(load_samples_data())
    
    # Create a selectInput with dynamically set choices
    selectInput("samples_columns", "Select Columns", choices, multiple = TRUE)
  })
  
  #samples functions

  
  # Creates the summary for the samples tab
  samples_summary <-
    function(dataf) {
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
  
  # Load the samples summary into an output object for the UI
  output$samples_summary <- renderTable(
    samples_summary(load_samples_data())
  ) 
  
  # Create table for the samples tab

  # Load the samples summary into an output object for the UI
  output$samples_table <- renderDT(
    datatable(load_samples_data(), options = list(ordering = TRUE))
  )
  
  # In order to choose what's going to be in the violin plot, we're going to make another output which is a list of column names 
  # For the radio buttons to choose which things to plot
  output$filtered_samples_table <- renderDT({
    # Subset the data based on the selected columns
    # Render the table
    datatable(column_filter(), options = list(dom = 't',pageLength = nrow(column_filter())))
  })
  
  
  # Render the plot or display for the selected numeric column
  output$selected_plot <- renderPlot({
    # Example: Create a histogram for the selected numeric column
    data <- load_samples_data()
    vector <- input$numeric_columns
    # print(vector)
    hist(data[[vector]], main = paste("Histogram of", vector))
  })
  
  # Counts Functions
  
  output$counts <- renderDT(
    datatable(load_counts_data(), options = list(ordering = TRUE))
  )
    
  
  
  # Make a table that summarizes the effect of filtering
  # Return a table with the number of samples, number of genes, number and % of genes passing filter and not passing filter
  counts_filtering_table <- function(dataf, variance, zero) {
    # # Filter the data
    # row_variances <- rowVars(as.matrix(dataf[,-1]), na.rm = TRUE)
    # 
    # # Filter the data based on variance threshold
    # filtered_counts <- dplyr::filter(dataf, row_variances >= variance)
    # filtered_counts <- dplyr::filter(filtered_counts, rowSums(filtered_counts == 0) >= zero)
    
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
  
  output$counts_filtering_table <- renderTable({
    counts_filtering_table(load_counts_data(), input$counts_variance, input$counts_non_zero)
  })
  
  
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
  
  # Display the plots
  output$plot_variance <- renderPlot({
    plots <- counts_scatter_plot(input$counts_variance, input$counts_non_zero)
    print(plots[[1]])
  })
  
  output$plot_zeros <- renderPlot({
    plots <- counts_scatter_plot(input$counts_variance, input$counts_non_zero)
    print(plots[[2]])
  })
  
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
  
  output$counts_heatmap<- renderPlot(
    counts_heatmap(input$counts_variance, input$counts_non_zero)
  )
  
  # Principal Component Analysis
  # Allow the scatter plot to be plotted against whichever PCs they want
  # Implies we need a drop down box for how many PCs they want the PCA to be by
  
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
  
  output$counts_pca<- renderPlot(
    counts_pca(input$counts_variance, input$counts_non_zero, input$num_components)
  )
  
  # DE Analysis Functions
  
  # Make a sortable data table with DE results

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
  }, height = 700)
  
  

  
  output$table <- renderDT(
    datatable(draw_table(load_DE_data(),input$slider))
  )
  
  
  
  
  # FGSEA Functions
  
  
  # Get list of diff expressed genes
  
  
  make_ranked_log2fc <- function(dataf) {
    new_table <- mutate(dataf, rank = rank(log2FoldChange) )%>%
      dplyr::arrange(rank) %>%
      dplyr::select(symbol, log2FoldChange)%>%
      dplyr::filter(is.finite(log2FoldChange))
    new_vector <- pull(new_table, log2FoldChange, symbol)
    return(new_vector)
  }
  
  
  #Run FGSEA
  
  run_fgsea <- reactive ({
    rnk_list <- make_ranked_log2fc(load_EA_data())
    c2_pathways <- gmtPathways("~/591_Final_Project/data/h.all.v2023.2.Hs.symbols.gmt")
    
    fgsea_results <- fgsea(c2_pathways, 
                           rnk_list, 
                           minSize=15, 
                           maxSize=500) %>% 
      as_tibble()
    
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
  
  #Barplot of FGSEA NES for top pathways from slider
  
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
  
  output$EA_barplot  <- renderPlot(
    EA_barplot(run_fgsea(), input$EA_p_adj)
  )
    
  #Sortable data table of EA results
  
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
  
  output$EA_Table  <- renderDT(
    datatable(EA_Table(input$EA_p_adj_table, input$EA_NES_select))
  )
    
  #Scatter plot of NES on x-axis and -log10 adjusted p-value on y-axis, with gene sets below threshold in grey color
  
  scatter_plot_data <- reactive({
    ea_table <- run_fgsea()
    
    # Filter gene sets below p-value threshold
    ea_table <- dplyr::filter(ea_table, padj < 10^(-input$EA_p_adj_plot))
    
    return(ea_table)
  })
  
  # Render the scatter plot
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
}

# Run the application 
shinyApp(ui = ui, server = server)
