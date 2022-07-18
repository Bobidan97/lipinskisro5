
#' Lipinskis Rule of 5 Function
#'
#' @param compound_dataset The input compound dataset
#'
#' @return New filtered compound dataset based on Rule of 5.
#' @export
#'
#' @examples
#' mol_table_construction(compound_dataset)
#' \dontrun{
#' mol_table_construction(compound_dataset)
#' }
mol_table_construction <- function(compound_dataset) {
  #read in csv file containing compound dataset
  mols <<- read_csv(file = compound_dataset)
  #get number of molecules in dataset
  num_mols <- nrow(mols)
  print(sprintf("%s, molecules in dataset.", num_mols))

  #filter dataset according to Lipinskis Rule of 5
  mols_filtered <<- mols %>%
    filter(MW < 500, cLogP < 5, TPSA < 140, RotatableBonds < 10, HBA < 10, HBD <5) %>%
    select(MW, cLogP, TPSA, RotatableBonds, HBA, HBD)

  #find number of compounds not meeeting Lipinskis requirements
  difference <- num_mols - nrow(mols_filtered)
  print(sprintf("%s, molecules removed from dataset due to not meeting Lipinski's Rules.", difference))
  #statistical summary
  summary <- summary(mols_filtered)
  print(summary)

}

#' Compound Dataset PCA
#'
#' @param mols_filtered Filtered compound dataset based on Rule of 5.
#'
#' @return PCA plot of filtered compound dataset.
#' @export
#'
#' @examples
#' pca(mols_filtered)
#' \dontrun{
#' pca(mols_filtered)
#' }
pca <- function(mols_filtered){
  mols.pca <<- prcomp(mols_filtered, center = TRUE, scale. = TRUE)
  summary(mols.pca)
  ggbiplot(mols.pca, ellipse = TRUE, labels = rownames(mols_filtered))

}

#' Radar Chart Summarising Compound Dataset
#'
#' @param mols_filtered
#'
#' @return A radar plot summarising Compound Dataset
#' @export
#'
#' @examples
#' radar_chart(mols_filtered)
#' \dontrun{
#' radar_chart(mols_filtered)
#' }
radar_chart <- function(mols_filtered) {
  #Define the variable ranges: maximum and minimum
  max_min <<- data.frame(
    MW = c(500, 0), cLogP = c(5, 0), TPSA = c(140, 0),
    RotatableBonds = c(10, 0), HBA = c(10, 0), HBD = c(5, 0)
  )
  rownames(max_min) <- c("Max", "Min")

  #Bind the variable ranges to the data
  df_radar <<- rbind(max_min, mols_filtered)

  radarchart(data = df_radar)

}
