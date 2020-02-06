// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/datamining/modules/visualization/VisualizerConfiguration.hpp>
#include <sgpp/datadriven/datamining/modules/visualization/VisualizerClassification.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingClustering.hpp>
#include <sgpp/datadriven/datamining/tools/hierarchyTree/ClusterNode.hpp>
#include <sgpp/base/tools/json/JSON.hpp>

#include <vector>
#include <string>


namespace sgpp {
namespace datadriven {

class VisualizerClustering : public VisualizerClassification {
 public:
  VisualizerClustering() = default;
  /**
   * Constructor given a configuration
   * @param config The VisualizerConfiguration object which contains
   * the configuration to run the visualization module
   */
  explicit VisualizerClustering(VisualizerConfiguration config);

  ~VisualizerClustering() override = default;
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
   * Method to generate and store in json format for the
   * plotly library the output of the clustering in a scatterplot.
   * @param matrix Matrix with the content to be stored
   * @param model Model used in the evaluation
   * @param currentDirectory The current directory to store the json file
   */
  void storeScatterPlotJson(DataMatrix &matrix, ModelFittingBase &model,
                       std::string currentDirectory) override;
  VisualizerDensityEstimation* visualizerDensityEstimation;

 private:
  /**
   * Method to generate and store in json format for the plotly library the plot containing
   * the full graph with the graph and evaluated densities
   * @param matrix Data Points to be plotted
   * @param model THe clustering model trained in the pipeline
   * @param currentDirectory Directory to store the data
   */
  void getDensityGraphPlot(DataMatrix &matrix,
    ModelFittingClustering &model, std::string currentDirectory);

  void getHierarchyAnimation(DataMatrix &matrix,
                    ModelFittingClustering &model, std::string currentDirectory);

  /**
   * Method which separates a matrix of points into multiple matrices depending on their value
   * obtained by the model. This makes the generation of he json easier
   * @param points The matrix of points to separate
   * @param labels THe labels associated to the points
   * @param traces vector of matrices in which the points will be distributed
   */
  void separateClustersIntoTraces(sgpp::base::DataMatrix &points,
    sgpp::base::DataVector &labels, std::vector<sgpp::base::DataMatrix> &traces);

  /**
   * Stores the cluster labels of all of the dimensionally
   * compressed points per each label in a csv file
   * @param matrix Matrix containing the 2D embedding of the points
   * @param model Trained model
   * @param currentDirectory Directory to store the files
   */
  void storeHierarchyCsv(DataMatrix &matrix,
                         ModelFittingClustering &model, std::string currentDirectory);
};
}  // namespace datadriven
}  // namespace sgpp
