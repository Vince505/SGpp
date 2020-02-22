// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org


#pragma once

#include <sgpp/datadriven/datamining/modules/dataSource/DataSource.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/visualization/VisualizerConfiguration.hpp>
#include <sgpp/datadriven/tools/CSVTools.hpp>
#include <string>

namespace sgpp {
namespace datadriven {

class Visualizer{
 public:
  Visualizer();
 /**
  * Virtual destructor
  */
  virtual ~Visualizer() = default;

  /**
   * Method to run the visualization process for a given batch and fold
   * @param model The model used to evaluate the visualization
   * @param dataSource The datasource from where the data points are obtained
   * @param epoch The current epoch running
   * @param fold The current fold being processed
   * @param batch The current batch being processed
   */
  virtual void runVisualization(ModelFittingBase &model, DataSource &dataSource, size_t epoch,
    size_t fold, size_t batch) = 0;

  /**
   * Method to run the visualization process when executing a post Process
   * @param model The model used to evaluate the visualization
   * @param dataSource The datasource from where the data points are obtained
   * @param fold The current fold being processed
   */
  virtual void runPostProcessingVisualization(ModelFittingBase &model, DataSource &dataSource,
    size_t fold = 0)= 0;

  /**
   * Runs the tsne algorithm to visualize high dimensional data in 2 dimensions
   * @param originalData Matrix with the original points in high dimensional space
   * @param compressedData Matrix which will contain the compressed points
   */
  void runTsne(DataMatrix &originalData, DataMatrix &compressedData);

  /**
   * Get the configuration of the visualizer object.
   * @return configuration of the visualizer object
   */
  const VisualizerConfiguration &getVisualizerConfiguration() const;

  /**
   * Method which starts the heatmap generation for Density Estimation Models
   * @param model The model used to evaluate the heatmap
   * @param currentDirectory The current directory to store the heatmap results
   * @param matrix The matrix containing the points to evaluate the heatmap
   * @param nDimensions Number of dimensions of the data to generate the Heatmap
 */
  void getHeatmap(ModelFittingBase &model, std::string currentDirectory, DataMatrix &matrix,
                  size_t &nDimensions);

  /**
   * Method which starts the linear cut generation for Density Estimation Models
   * @param model The model used to evaluate the linear cuts
   * @param currentDirectory The current directory to store the linear cuts results
   * @param matrix The matrix containing the points to evaluate the cuts
   * @param nDimensions Number of dimensions of the data to generate the linear cuts
   */
  void getLinearCuts(ModelFittingBase &model, std::string currentDirectory, DataMatrix &matrix,
                     size_t &nDimensions);

  void setResolution(size_t resolution);

 protected:
  /**
   * Method to create a folder based on the operating system of the user
   * @param folder_path absolute or relative path of the folder to be created
   */
  void createFolder(std::string folder_path);

  /**
   * Method which builds the matrices used to generate a heatmap of a given model
   * @param heatmapMatrix matrix to be initialized
   * @param nDimensions Number of dimensions of the data to generate the Heatmap
   */
  void initializeHeatmapMatrix(DataMatrix &heatmapMatrix, size_t &nDimensions);

  /**
   * Method which builds the matrices used to generate the linear Cuts of a given model
   * @param cutMatrix matrix to be initialized
   * @param nDimensions Number of dimensions of the data to generate the linear cuts
   */
  void initializeCutMatrix(DataMatrix &cutMatrix, size_t &nDimensions);

  /**
 * Method which generates the linear cuts graphs for models of 3 or more dimensions
 * @param model the model used to evaluate the linear cuts
 * @param currentDirectory The current directory to store the linear cuts results
 * @param matrix The matrix containing the points to evaluate the cuts
 */
  void getLinearCutsMore3D(ModelFittingBase &model, std::string currentDirectory,
                           DataMatrix &matrix);

  /**
   * Method which generates the linear cuts graphs for models of 1 dimension
   * @param model the model used to evaluate the linear cuts
   * @param currentDirectory The current directory to store the linear cuts results
   * @param matrix The matrix containing the points to evaluate the cuts
   */
  void getLinearCuts1D(ModelFittingBase &model, std::string currentDirectory, DataMatrix &matrix);

  /**
   * Method which generates the linear cuts graphs for models of 2 dimensions
   * @param model the model used to evaluate the linear cuts
   * @param currentDirectory The current directory to store the linear cuts results
   * @param matrix The matrix containing the points to evaluate the cuts
   */
  void getLinearCuts2D(ModelFittingBase &model, std::string currentDirectory, DataMatrix &matrix);

  /**
   * Method which generates the heatmap of models of 4 or more dimensions
   * @param model The model used to evaluate the heatmap
   * @param currentDirectory The current directory to store the heatmap results
   * @param matrix The matrix containing the points to evaluate the heatmap
   */
  void getHeatmapMore4D(ModelFittingBase &model, std::string currentDirectory, DataMatrix &matrix);

  /**
   * Method which generates the heatmap of models of 3 dimensions
   * @param model The model used to evaluate the heatmap
   * @param currentDirectory The current directory to store the heatmap results
   * @param matrix The matrix containing the points to evaluate the heatmap
   */
  void getHeatmap3D(ModelFittingBase &model, std::string currentDirectory, DataMatrix &matrix);

  /**
   * Method which generates the heatmap of models of 2 dimensions
   * @param model The model used to evaluate the heatmap
   * @param currentDirectory The current directory to store the heatmap results
   * @param matrix The matrix containing the points to evaluate the heatmap
   */
  void getHeatmap2D(ModelFittingBase &model, std::string currentDirectory, DataMatrix &matrix);

  /**
   * Method which shifts one position the columns of a matrix from left to right
   * in a circular fashion until the column given by the parameters maxColumns
   * @param matrix The matrix to be shifted
   * @param maxColumns The max number of columns used when shifting
   */
  void translateColumns(DataMatrix &matrix, size_t maxColumns);

  /**
   * Method which shifts the columns given by the vector indexes
   * of a matrix from left to right
   * in a circular fashion. If indexes are <1,3,6> Then column 1 will be shifted to
   * 3, 3 to 6 and 6 to 1.
   * @param matrix The matrix to be shifted
   * @param indexes Vector containing the columns to shift
   */
  void translateColumnsRight(DataMatrix &matrix, std::vector<size_t> indexes);

  /**
   * Method which shifts the columns given by the vector indexes
   * of a matrix from right to left
   * in a circular fashion. If indexes are <1,3,6> Then column 1 will be shifted to
   * 6, 6 to 3 and 1 to 6.
   * @param matrix The matrix to be shifted
   * @param indexes Vector containing the columns to shift
   */
  void translateColumnsLeft(DataMatrix &matrix, std::vector<size_t> indexes);

  /**
   * Method to update the columns indexes to be shifted when generating
   * the linear cuts
   * @param columnIndexes Vector to update
   * @param matrix The matrix whose columns are being shifted
   */
  void updateIndexesCuts(std::vector<size_t> &columnIndexes, DataMatrix &matrix);

  /**
   * Method to update the columns indexes to be shifted when generating
   * the heatmap
   * @param columnIndexes Vector to update
   * @param matrix The matrix whose columns are being shifted
   */
  void updateIndexesHeatmap(std::vector<size_t> &columnIndexes, DataMatrix &matrix);

  /**
   * Method to swap to columns of a matrix
   * @param matrix The matrix whose columns are to be swaped
   * @param col1 The column identifier to be swaped with the column identified by col2
   * @param col2 The colum identifier to be swaped with the column identified by col1
   */
  void swapColumns(DataMatrix &matrix, size_t col1, size_t col2);


  /**
   * Method to generate and store in json  format for the
   * plotly library the output of the data in a scatterplot.
   * @param matrix Matrix with the content to be stored
   * @param model Model used in the evaluation
   * @param currentDirectory The current directory to store the json file
   */
  virtual void storeScatterPlotJson(DataMatrix &matrix, ModelFittingBase &model,
                            std::string currentDirectory) = 0;

  /**
   * Method to generate and store in json  format for the
   * plotly library the output of the linear cuts for models of 2 or more dimensions
   * @param matrix Matrix with the content to be stored
   * @param indexes Vectors containing the dimensions used when generating these cuts
   * @param varDim THe dimension number varying and whose evaluation is shown in the model
   * @param filepath The current directory to store the json file
   */
  virtual void storeCutJson(DataMatrix &matrix,
                    std::vector<size_t> indexes, size_t &varDim, std::string filepath) = 0;

  /**
   * Method to generate and store in json  format for the
   * plotly library the output of the linear cuts for models of 1 dimension
   * @param matrix Matrix with the content to be stored
   * @param filepath The current directory to store the json file
   */
  virtual void storeCutJson(DataMatrix &matrix, std::string filepath) = 0;

  /**
   * Method to generate and store in json format for the
   * plotly library the output of the classification
   * heatmaps for models of 2 dimensions
   * @param matrix Matrix with the content to be stored
   * @param model The model used when evaluating the heatmaps
   * @param filepath The current directory to store the json file
   */
  virtual void storeHeatmapJson(DataMatrix &matrix, ModelFittingBase &model,
                                      std::string filepath)= 0;

  /**
   * Method to generate and store in json  format for the
   * plotly library the output of the classification heatmaps for models of 3 or more dimensions
   * @param matrix Matrix with the content to be stored
   * @param model The model used when evaluating the heatmaps
   * @param indexes Vectors containing the dimensions used when generating these heatmaps
   * @param varDim1 The first dimension number varying and whose evaluation
   * is shown in the model
   * @param varDim2 The second dimension number varying and whose evaluation
   * is shown in the model
   * @param filepath The current directory to store the json file
   */
  virtual void storeHeatmapJson(DataMatrix &matrix, ModelFittingBase &model,
                                      std::vector<size_t> indexes, size_t &varDim1,
                                      size_t &varDim2, std::string filepath) = 0;
  /**
   * Configuration object for the fitter.
   */
  VisualizerConfiguration config;

  /**
   * Matrix with reduced dimensional data
   */
  DataMatrix compressedData;
  /**
   * Resolution used in the graphs
   */
  size_t resolution;

  /**
   * Variable to store the current folder in which the visualizer is working on
   */
  std::string currentDirectory;
};

}  // namespace datadriven
}  // namespace sgpp
