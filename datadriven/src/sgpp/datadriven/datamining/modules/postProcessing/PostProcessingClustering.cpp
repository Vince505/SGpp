// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/postProcessing/PostProcessingClustering.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingClustering.hpp>
#include <sgpp/datadriven/datamining/tools/Graph.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>
#include <sgpp/datadriven/tools/CSVWriter.hpp>

#include<map>

using sgpp::datadriven::CSVWriter;
namespace sgpp {
namespace datadriven {

void PostProcessingClustering::postProcessing(DataSource &datasource,
  ModelFittingBase &model, Visualizer &visualizer, size_t fold) {
  std::cout << "Starting post processing" << std::endl;
  ModelFittingClustering *clusteringModel =
    dynamic_cast<ModelFittingClustering *>(&model);

  clusteringModel->updateVpTree(datasource.getAllSamples()->getData());
  clusteringModel->generateSimilarityGraph();

  clusteringModel->intializeHierarchyTree();
  // Obtaining the cluster hierarchy for this batch
  std::map<UndirectedGraph::vertex_descriptor, size_t> clusterMap;

  double minThreshold =
    clusteringModel->getFitterConfiguration().getClusteringConfig().minDensityThreshold;
  double maxThreshold =
    clusteringModel->getFitterConfiguration().getClusteringConfig().maxDensityThreshold;
  size_t steps = clusteringModel->getFitterConfiguration().getClusteringConfig().steps;

  if (minThreshold > maxThreshold) {
    maxThreshold = minThreshold;
  }
  double stepIncrease = (maxThreshold-minThreshold)/static_cast<double>(steps);

  for (double step = minThreshold;
       step <= maxThreshold; step+=stepIncrease) {
    std::cout << "============= Step " << step << " =============" << std::endl;
    clusterMap.clear();
    clusteringModel->copyPreviousGraphStep();
    clusteringModel->applyDensityThresholds(step);
    clusteringModel->detectComponentsAndLabel(clusterMap);
    clusteringModel->getHierarchy(clusterMap, step);

    if (stepIncrease == 0) {
      break;
    }
  }

  if (clusteringModel->getFitterConfiguration().
    getClusteringConfig().storeHierarchy) {
    // Storing the metadata of the hierarchical tree
    clusteringModel->storeHierarchyTree();

    // Storing the evaluated labels

    Dataset *allSamples = datasource.getAllSamples();
    DataVector &targets = allSamples->getTargets();

    clusteringModel->evaluate(allSamples->getData(), targets);

    CSVWriter::writeCSVToFile(clusteringModel->getFitterConfiguration().
      getClusteringConfig().outputDirectory+"/hierarchyLabels.csv", *allSamples);
  }

  visualizer.runPostProcessingVisualization(model, datasource, fold);
}

}  // namespace datadriven
}  // namespace sgpp
