
#' Lipinskis Rule of 5
#'
#' @param compound_dataset The input compound dataset
#'
#' @return New filtered compound dataset based on Lipinskis Rule of 5 plus PCA
#' and Radar Chart
#' @export
#'
#' @importFrom readr read_csv
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom stats prcomp
#' @import fmsb
#' @import ggbiplot
#' @examples
#'
lipro5 <- function(compound_dataset) {
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


  mols.pca <<- prcomp(filtered_dataset, center = TRUE, scale. = TRUE)
  pca_summary <- summary(mols.pca)
  print(pca_summary)
  ggbiplot(mols.pca, ellipse = TRUE, labels = rownames(filtered_dataset))

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
