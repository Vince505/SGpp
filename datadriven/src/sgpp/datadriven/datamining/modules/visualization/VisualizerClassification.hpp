// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/datamining/modules/dataSource/DataSource.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/visualization/VisualizerConfiguration.hpp>
#include <sgpp/datadriven/datamining/modules/visualization/VisualizerDensityEstimation.hpp>
#include <sgpp/datadriven/tools/CSVTools.hpp>
#include <vector>
#include <string>

namespace sgpp {
namespace datadriven {

class VisualizerClassification:public Visualizer {
 public:
  VisualizerClassification()= default;

  /**
   * Constructor given a configuration
   * @param config The VisualizerConfiguration object which contains
   * the configuration to run the visualization module
   */
  explicit VisualizerClassification(VisualizerConfiguration config);

  ~VisualizerClassification() override = default;

  /**
   * Method to run the visualization process for a given batch and fold
   * @param model The model used to evaluate the visualization
   * @param dataSource The datasource from where the data points are obtained
   * @param epoch The current epoch running
   * @param fold The current fold being processed
   * @param batch The current batch being processed
   */
  void runVisualization(ModelFittingBase &model, DataSource &dataSource, size_t epoch,
    size_t fold, size_t batch) override;

  /**
   * Method to run the visualization process when executing a post Process
   * @param model The model used to evaluate the visualization
   * @param dataSource The datasource from where the data points are obtained
   * @param fold The current fold being processed
   */
  void runPostProcessingVisualization(ModelFittingBase &model, DataSource &dataSource,
    size_t fold = 0) override;

 protected:
  /**
   * Method to generate and store in json  format for the
   * plotly library the output of the linear cuts for models of 2 or more dimensions
   * @param matrix Matrix with the content to be stored
   * @param indexes Vectors containing the dimensions used when generating these cuts
   * @param varDim THe dimension number varying and whose evaluation is shown in the model
   * @param filepath The current directory to store the json file
   */
  void storeCutJson(DataMatrix &matrix,
                    std::vector<size_t> indexes, size_t &varDim, std::string filepath) override;

  /**
   * Method to generate and store in json  format for the
   * plotly library the output of the linear cuts for models of 1 dimension
   * @param matrix Matrix with the content to be stored
   * @param filepath The current directory to store the json file
   */
  void storeCutJson(DataMatrix &matrix, std::string filepath) override;

  /**
   * Method to generate and store in json  format for the
   * plotly library the output of the data with its predicted labels in a scatterplot.
   * @param matrix Matrix with the content to be stored
   * @param model Model used in the evaluation
   * @param currentDirectory The current directory to store the json file
   */
  void storeScatterPlotJson(DataMatrix &matrix, ModelFittingBase &model,
                            std::string currentDirectory) override;

  /**
   * Method to generate and store in json format for the
   * plotly library the output of the classification
   * heatmaps for models of 2 dimensions
   * @param matrix Matrix with the content to be stored
   * @param model The model used when evaluating the heatmaps
   * @param filepath The current directory to store the json file
   */
  void storeHeatmapJson(DataMatrix &matrix, ModelFittingBase &model,
  std::string filepath) override;

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
  void storeHeatmapJson(DataMatrix &matrix, ModelFittingBase &model,
  std::vector<size_t> indexes, size_t &varDim1, size_t &varDim2, std::string filepath) override;


  /**
   * List of colors to add to the grids per class in high dimensional cases
   */
  std::vector<std::string> colors = {"red", "darkviolet", "orange", "palegreen",
                    "plum", "purple", "chocolate", "darkcyan", "gold", "tomato"};


  /**
   * Vector which contains the all of the class label values in the model
   */
  DataVector classes;

  VisualizerDensityEstimation* visualizerDensityEstimation;

  DataMatrix originalData;

};

}  // namespace datadriven
}  // namespace sgpp
