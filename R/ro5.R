
#' Lipinskis Rule of 5 Filter
#'
#' @param compound_dataset The input compound dataset
#'
#' @return New filtered compound dataset based on Lipinskis Rule of 5.
#' @export
#'
#' @importFrom readr read_csv
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @examples
#' ro5filter(compound_dataset)
#' \dontrun{
#' ro5filter(compound_dataset)
#' }
ro5filter <- function(compound_dataset) {
  #read in csv file containing compound dataset
  mols <<- read_csv(file = compound_dataset)
  #get number of molecules in dataset
  num_mols <- nrow(mols)
  print(sprintf("%s, molecules in dataset.", num_mols))

  #filter dataset according to Lipinskis Rule of 5
  filtered_dataset <<- mols %>%
    filter(MW < 500, cLogP < 5, RotatableBonds < 10, HBA < 10, HBD <5) %>%
    select(MW, cLogP, RotatableBonds, HBA, HBD)

  #find number of compounds not meeeting Lipinskis requirements
  difference <- num_mols - nrow(filtered_dataset)
  print(sprintf("%s, molecules removed from dataset due to not meeting Lipinski's Rules.", difference))
  #statistical summary
  summary <- summary(filtered_dataset)
  print(summary)

}

#' Compound Dataset PCA
#'
#' @param filtered_dataset Filtered compound dataset based on Rule of 5.
#'
#' @return PCA plot of filtered compound dataset.
#' @export
#'
#' @importFrom stats prcomp
#' @import ggbiplot
#' @examples
#' pca(filtered_dataset)
#' \dontrun{
#' pca(filtered_dataset)
#' }
pca <- function(filtered_dataset) {
  mols.pca <<- prcomp(filtered_dataset, center = TRUE, scale. = TRUE)
  pca_summary <- summary(mols.pca)
  print(pca_summary)
  ggbiplot(mols.pca, ellipse = TRUE, labels = rownames(filtered_dataset))

}

#' Radar Chart Summarising Compound Dataset
#'
#' @param filtered_dataset
#'
#' @return A radar plot summarising Compound Dataset
#' @export
#'
#' @import fmsb
#'
#' @examples
#' radar_chart(filtered_dataset)
#' \dontrun{
#' radar_chart(filtered_dataset)
#' }
radar_chart <- function(filtered_dataset) {
  #Define the variable ranges: maximum and minimum
  max_min <<- data.frame(
    MW = c(500, 0), cLogP = c(5, 0),
    RotatableBonds = c(10, 0), HBA = c(10, 0), HBD = c(5, 0)
  )
  rownames(max_min) <- c("Max", "Min")

  #Bind the variable ranges to the data
  df_radar <<- rbind(max_min, filtered_dataset)

  radarchart(data = df_radar)

}
