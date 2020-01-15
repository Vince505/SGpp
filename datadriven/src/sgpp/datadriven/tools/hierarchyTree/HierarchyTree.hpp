// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include  <sgpp/datadriven/tools/hierarchyTree/ClusterNode.hpp>

namespace sgpp {
namespace datadriven {


class HierarchyTree {

public:

  explicit HierarchyTree(std::vector<size_t> vertexIndexes);

  ~HierarchyTree() {
    delete root;
  }

  void printTree();

  ClusterNode* getRoot();

  ClusterNode* getMostSpecificCluster(size_t vertexIndex);

  void postProcessing();
private:

  ClusterNode* root;
  ClusterNode* findNode(ClusterNode* node, size_t vertexIndex);
  void printNode(ClusterNode* node);
  void setPostProcessingLabel(ClusterNode* node, int &clusterLabel);

};

}  // namespace datadriven
}  // namespace sgpp
