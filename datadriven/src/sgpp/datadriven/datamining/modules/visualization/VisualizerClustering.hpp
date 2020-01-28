// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/datamining/modules/visualization/VisualizerConfiguration.hpp>
#include <sgpp/datadriven/datamining/modules/visualization/VisualizerClassification.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingClustering.hpp>
#include <sgpp/datadriven/tools/hierarchyTree/ClusterNode.hpp>
#include <vector>
#include <string>

#include <sgpp/base/tools/json/JSON.hpp>
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
   */
  void runPostProcessingVisualization(ModelFittingBase &model, DataSource &dataSource) override;

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
  void getGraphPlot(DataMatrix &matrix,
    ModelFittingClustering &model, std::string currentDirectory);

  void getHierarchyAnimation(DataMatrix &matrix,
                    ModelFittingClustering &model, std::string currentDirectory);

  void separateClustersIntoTraces(sgpp::base::DataMatrix &points,
    sgpp::base::DataVector &labels, std::vector<sgpp::base::DataMatrix> &traces);
};
}  // namespace datadriven
}  // namespace sgpp
