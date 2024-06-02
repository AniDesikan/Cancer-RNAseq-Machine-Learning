library('tidyverse')
library('RColorBrewer')

#' Read the expression data "csv" file as a dataframe, not tibble
#'
#' @param filename (str): the path of the file to read
#' @param delimiter (str): generalize the function so it can read in data with
#'   your choice of delimiter
#'
#' @return A dataframe containing the example intensity data with rows as probes
#'   and columns as samples
#' @export
#'
#' @examples
read_data <- function(intensity_data, delimiter) {
    return(read.delim(intensity_data, sep = delimiter))
}

#' Define a function to calculate the proportion of variance explained by each PC
#'
#' @param pca_results (obj): the results returned by `prcomp()`
#'
#' @return A vector containing the values of the variance explained by each PC
#' @export
#'
#' @examples
calculate_variance_explained <- function(pca_results) {
  # variance = sd squared
  #We need to take variance of each PC and divide it by the sum of variance of all PCs
  # use pca_results$dev, make a new tibble with column called variance
  prop_var <- tibble( variance = pca_results$sdev ** 2)
  # make a new variable, varsum, which is the values in the column added up
  varsum = sum(prop_var$variance)
  # mutate a new column propsum which is the variance divided by the varsum
  prop_var <- dplyr::mutate(prop_var, prop_variance = variance / varsum)
  # return the new propsum column
  return(prop_var$prop_variance)
}

#' Define a function that takes in the variance values and the PCA results to
#' make a tibble with PC names, variance explained by each PC, and the
#' cumulative sum of variance explained. These columns should be named 
#' "principal_components", "variance_explained", and "cumulative", respectively.
#' 
#'
#'
#' @param pca_ve (vector): the vector generated in the previous function with
#'   the variance explained values
#' @param pca_results (object): the results returned by `prcomp()`
#' @return A tibble that contains the names of the PCs, the individual variance
#'   explained, and the cumulative variance explained with names described above
#' @export
#' @examples 
make_variance_tibble <- function(pca_ve, pca_results) {
  x <- pca_results$x
  # Make a vector that has the cumulative variance
  Cum_vector <- vapply(seq_along(pca_ve), function(i) sum(pca_ve[1:i]), 1)
    #Make a tibble with three columns, (names of the PCs, variance of each pc, and the sum of variance(PC1 explains 35%, PC1 and PC2 together explain 42%, etc.))
  variance_tibble <- tibble(principal_components =  colnames(x), variance_explained = pca_ve, cumulative = Cum_vector) %>%
    dplyr::arrange(cumulative)
  return(variance_tibble)
}



#' Define a function to create a biplot of PC1 vs. PC2 labeled by
#' SixSubTypesClassification
#'
#' @param metadata (str): The path to the proj_metadata.csv file
#' @param pca_results (obj): The results returned by `prcomp()`
#'
#' @return A ggplot consisting of a scatter plot of PC1 vs PC2 labeled by
#'   SixSubTypesClassification found in the metadata
#' @export
#'
#' @examples
make_biplot <- function(metadata, pca_results) {
  x <- pca_results$x
  metadata_r <- read_csv(metadata)
  matched <- match(rownames(x), metadata_r$geo_accession)
  plot <- as_tibble(pca_results$x) %>%
    mutate(type= metadata_r$SixSubtypesClassification[matched]) %>%
    ggplot(aes(x=PC1,y=PC2, color = type)) +
    geom_point()
  return(plot)
}

#' Define a function to return a list of probeids filtered by signifiance
#'
#' @param diff_exp_tibble (tibble): A tibble containing the differential expression results
#' @param fdr_threshold (float): an appropriate FDR threshold, we will use a
#'   value of .01. This is the column "padj" in the tibble.
#'
#' @return A list with the names of the probeids passing the fdr_threshold
#' @export
#'
#' @examples
list_significant_probes <- function(diff_exp_tibble, fdr_threshold) {
  #This returns list of probe ids, which is the first column, filtered by if it meets the fdr threshold
  #So we can use the dplyr::filter to see if the padj is greater than fdr_threshold
  list <- dplyr::filter(diff_exp_tibble, padj < fdr_threshold)
    return(list$probeid)
}

#' Define a function that uses the list of significant probeids to return a
#' matrix with the intensity values for only those probeids.
#' @param intensity (dataframe): The dataframe of intensity data generated in
#'   part 1
#' @param sig_ids_list (list/vector): The list of differentially expressed
#'   probes generated in part 6
#'
#' @return A `matrix()` of the probe intensities for probes in the list of
#'   significant probes by FDR determined in the previous function.
#'
#' @export
#'
#' @examples
return_de_intensity <- function(intensity, sig_ids_list) {
  #First we find the indices of the matching probeids by using the match()
  matching <- match( sig_ids_list,rownames(intensity))
  intensity_matrix <- as.matrix(intensity[matching, ])
    return(intensity_matrix)
}

#' Define a function that takes the intensity values for significant probes and
#' creates a color-blind friendly heatmap
#'
#' @param de_intensity (matrix): The matrix of intensity values for significant
#'   differentially expressed probes returned in part 7
#' @param num_colors (int): The number of colors in a specificed RColorBrewer
#'   palette
#' @param palette (str): The name of the chosen RColorBrewer palette
#'
#' @return A heatmap displaying the intensity values for the differentially
#'   expressed probes
#' @export
#'
#' @examples
plot_heatmap <- function(de_intensity, num_colors, palette) {
  color <- brewer.pal(num_colors, palette)
  map <- heatmap(de_intensity, col = color)
    return(map)
}