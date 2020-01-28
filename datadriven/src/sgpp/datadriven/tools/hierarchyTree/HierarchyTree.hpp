// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include  <sgpp/datadriven/tools/hierarchyTree/ClusterNode.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

namespace sgpp {
namespace datadriven {


class HierarchyTree {

public:

  explicit HierarchyTree(size_t numberPoints);

  ~HierarchyTree() {
    delete root;
  }

  void printTree();

  ClusterNode* getRoot();

  ClusterNode* getMostSpecificCluster(size_t vertexIndex);

  ClusterNode* getMostSpecificClusterAtLevel(size_t vertexIndex, size_t level);

  size_t getNumberClusters();

  size_t getNumberLevels();

  void postProcessing();

  void evaluateClustering(sgpp::base::DataVector &results);

  void evaluateClusteringAtLevel(sgpp::base::DataVector &results, size_t level);
private:

  ClusterNode* root;
  ClusterNode* findNode(ClusterNode* node, size_t vertexIndex);

  ClusterNode* findNodeAtLevel(ClusterNode* node, size_t vertexIndex, size_t level);

  void countCluster(ClusterNode* node, size_t &numberClusters);

  void countLevel(ClusterNode* node, size_t &maxLevel);

  void printNode(ClusterNode* node);
  void setPostProcessingLabel(ClusterNode* node, int &clusterLabel);

};

}  // namespace datadriven
}  // namespace sgpp
