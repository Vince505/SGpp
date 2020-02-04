// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/application_exception.hpp>

#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingClustering.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationCG.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationCombi.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationOnOff.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationOnOffParallel.hpp>

#include <map>
#include <iostream>
#include <ctime>
#include <vector>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::Grid;
using sgpp::base::application_exception;

namespace sgpp {
namespace datadriven {

ModelFittingClustering::ModelFittingClustering(
  const FitterConfigurationClustering &config)
  : ModelFittingBase(), refinementsPerformed{0} {
  this->verboseSolver = true;
  this->config = std::unique_ptr<FitterConfiguration>(
    std::make_unique<FitterConfigurationClustering>(config));

  this->densityEstimationModel = createNewDensityModel(
    dynamic_cast<sgpp::datadriven::FitterConfigurationDensityEstimation &>(*(this->config)));

  this->graph = nullptr;
  this->vpTree = nullptr;
  this->hierarchy = nullptr;

#ifdef USE_SCALAPACK
  auto& parallelConfig = this->config->getParallelConfig();
if (parallelConfig.scalapackEnabled_) {
  processGrid = std::make_shared<BlacsProcessGrid>(config.getParallelConfig().processRows_,
                                                   config.getParallelConfig().processCols_);
}
#endif
}

void ModelFittingClustering::fit(Dataset &newDataset) {
  reset();
  update(newDataset);
}

void ModelFittingClustering::update(Dataset &newDataset) {
  // Update the model
  dataset = &newDataset;
  generateDensityEstimationModel(newDataset);
}

double ModelFittingClustering::evaluate(const DataVector &sample) {
  std::cout << "Method not defined for Clustering Models" << std::endl;
  return 0.0;
}

void ModelFittingClustering::evaluate(DataMatrix &samples, DataVector &results) {
  if (hierarchy == nullptr) {
    densityEstimationModel->evaluate(samples, results);
  } else {
    DataVector evaluateTemp;
    hierarchy->evaluateClustering(evaluateTemp);
    results.resize(samples.getNrows());

    for (size_t index = 0; index < samples.getNrows(); index++) {
      DataVector currentRow(samples.getNcols());
      samples.getRow(index, currentRow);
      size_t vpIndex = vpTree->getIndexedKeyFromPoint(currentRow);
      if (vpIndex == getPoints().getNrows()) {
        std::cout << "Point "+ currentRow.toString() +" wasn't used for clustering"<< std::endl;
        results.set(index, -1.0);
      } else {
        results.set(index, evaluateTemp.get(vpIndex));
      }
    }
  }
}

bool ModelFittingClustering::refine() {
  densityEstimationModel->refine();
  return true;
}

void ModelFittingClustering::reset() {
  densityEstimationModel->reset();
  graph.reset();
  prunedGraphPreviousStep .reset();
  hierarchy.reset();
}

double ModelFittingClustering::computeResidual(DataMatrix &validationData) const  {
throw sgpp::base::not_implemented_exception(
  "ModelFittingClustering::computeResidual() is not implemented!");
}

void ModelFittingClustering::resetTraining() {
  densityEstimationModel->resetTraining();
}

void ModelFittingClustering::updateRegularization(double lambda) {
  densityEstimationModel->updateRegularization(lambda);
}

std::unique_ptr<ModelFittingDensityEstimation> ModelFittingClustering::createNewDensityModel(
    sgpp::datadriven::FitterConfigurationDensityEstimation& densityEstimationConfig) {
  if (densityEstimationConfig.getGridConfig().generalType_ ==
      base::GeneralGridType::ComponentGrid) {
    return std::make_unique<ModelFittingDensityEstimationCombi>(densityEstimationConfig);
  }
  switch (densityEstimationConfig.getDensityEstimationConfig().type_) {
    case DensityEstimationType::CG: {
      return std::make_unique<ModelFittingDensityEstimationCG>(densityEstimationConfig);
    }
    case DensityEstimationType::Decomposition: {
#ifdef USE_SCALAPACK
      if (densityEstimationConfig.getParallelConfig().scalapackEnabled_) {
    return std::make_unique<ModelFittingDensityEstimationOnOffParallel>(densityEstimationConfig,
                                                                        processGrid);
  }
#endif  // USE_SCALAPACK
      return std::make_unique<ModelFittingDensityEstimationOnOff>(densityEstimationConfig);
    }
  }
  throw application_exception("Unknown density estimation type");
}

void ModelFittingClustering::generateDensityEstimationModel(Dataset &dataset) {
  densityEstimationModel->update(dataset);
}

void ModelFittingClustering::generateSimilarityGraph() {
  graph = std::make_shared<Graph>(getPoints().getNrows());
  DataVector currrentRow(getPoints().getNcols());

  for (size_t index = 0; index < getPoints().getNrows(); index++) {
    getPoints().getRow(index, currrentRow);
    auto nearestNeighbors = vpTree->getNearestNeighbors(currrentRow,
        config->getClusteringConfig().noNearestNeighbors);
    graph->createEdges(index, nearestNeighbors);
  }
  std::cout << "Num of vertices: " << graph->getNumberVertices() << std::endl;
  std::cout << "Num of edges " << graph->getNumberEdges() << std::endl;
}

void ModelFittingClustering::updateVpTree(DataMatrix &newDataset) {
  clock_t start = std::clock();
  if (vpTree == nullptr) {
    this->vpTree = std::make_unique<VpTree>(newDataset);
    clock_t end = std::clock();
    std::cout << "VpTree updated in  " <<
              std::to_string(static_cast<double>(end - start) / CLOCKS_PER_SEC) << " seconds"
              << std::endl;
  }
}

void ModelFittingClustering::applyDensityThresholds(double densityThreshold) {
  DataMatrix points = getPoints();
  DataVector evaluation(points.getNrows());

  densityEstimationModel->evaluate(points, evaluation);

  double maxValue = evaluation.max();
  for (size_t index = 0; index < evaluation.size(); index++) {
    if (graph->containsVertex(index)) {
      if (evaluation.get(index) < densityThreshold * maxValue) {
        graph->removeVertex(index);
      }
    }
  }
  std::cout << "Remaining vertices: " <<
  boost::num_vertices(*(graph->getGraph())) <<std::endl;
}

void ModelFittingClustering::detectComponentsAndLabel(
  std::map<UndirectedGraph::vertex_descriptor, size_t> &clusterMap) {
  auto numberComponents = graph->getConnectedComponents(clusterMap);
  std::cout << "Number of found components: "<< numberComponents <<std::endl;
}

void ModelFittingClustering::getHierarchy(
  std::map<UndirectedGraph::vertex_descriptor, size_t> &clusterMap, double densityThreshold) {
  std::map<size_t, std::vector<size_t>> labelstoPointsMap;

  for (auto vertex : clusterMap) {
    labelstoPointsMap[vertex.second].push_back(graph->getIndex(vertex.first));
  }
  std::cout << "Building the hierarchy" << std::endl;
  std::vector<ClusterNode*> updatedClusters;
  std::vector<size_t> visitedLabels;

  for (auto cluster : labelstoPointsMap) {
    ClusterNode* newChild = new ClusterNode(static_cast<int>(cluster.first),
              cluster.second, densityThreshold);
    ClusterNode* parentCluster = hierarchy->getMostSpecificCluster(cluster.second.at(0));
    parentCluster->addChild(newChild);
    newChild->setLevel(parentCluster->getLevel()+1);
    if (std::find(visitedLabels.begin(), visitedLabels.end(),
      parentCluster->getClusterLabel()) == visitedLabels.end()) {
      updatedClusters.push_back(parentCluster);
      visitedLabels.push_back(parentCluster->getClusterLabel());
    }
  }
  for (auto parent : updatedClusters) {
    if (parent->getChildren().size() > 1) {
      if (parent->split(prunedGraphPreviousStep, densityThreshold)) {
        // Check if the node is the root. We cannot split the root
        if (parent->getParent() != nullptr) {
          // Assigning the level of the node to its children
          for (auto child : parent->getChildren()) {
            child->setLevel(parent->getLevel());
          }

          // Linking the children to the parent of the node
          parent->getParent()->addChildren(parent->getChildren());

          // Removing the node from the tree
          parent->getParent()->removeChild(parent);
        }
      }
    } else {
      parent->removeChildren();
    }
  }
  hierarchy->postProcessing();
}

void ModelFittingClustering::intializeHierarchyTree() {
  this->hierarchy = std::make_unique<HierarchyTree>(getPoints().getNrows());
}

void ModelFittingClustering::copyPreviousGraphStep() {
  prunedGraphPreviousStep = std::make_shared<Graph>(*graph);
}

std::unique_ptr<ModelFittingDensityEstimation>*
    ModelFittingClustering::getDensityEstimationModel() {
  return &densityEstimationModel;
}

DataMatrix& ModelFittingClustering::getPoints() const {
  return vpTree->getStoredItems();
}

std::unique_ptr<HierarchyTree>* ModelFittingClustering::getHierarchyTree() {
  return &hierarchy;
}

std::shared_ptr<Graph> ModelFittingClustering::getGraph() {
  return graph;
}
}  // namespace datadriven
}  // namespace sgpp
