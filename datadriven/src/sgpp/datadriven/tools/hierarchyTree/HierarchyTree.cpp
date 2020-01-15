// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org


#include  <sgpp/datadriven/tools/hierarchyTree/HierarchyTree.hpp>

#include <iostream>

namespace sgpp {
namespace datadriven {

HierarchyTree::HierarchyTree(std::vector<size_t> vertexIndexes) {
  this->root = new ClusterNode(-1, vertexIndexes, 0.0);
}

void HierarchyTree::printTree() {
  auto cluster = root->getVertexIndexes();
  std::cout << "=======================HIERACHY TREE=========================" << std::endl;
  std::cout << "Cluster Label: " << root->getClusterLabel() << std::endl;
  std::cout << "Indexes: ";
  for (auto vertexIndex: cluster) {
    std::cout << vertexIndex << ", ";
  }
  std::cout << std::endl;

  printNode(root);
}
void HierarchyTree::printNode(ClusterNode* node) {
  for(auto child: node->getChildren()) {
    auto cluster = child->getVertexIndexes();
    std::cout << "================================================"<< std::endl;
    std::cout << "Parent: " << child->getParent()->getClusterLabel();
    std::cout << std::endl;
    std::cout << "Cluster Label: " << child->getClusterLabel() << std::endl;
    std::cout << "Indexes: ";
    for (auto vertexIndex: cluster) {
      std::cout << vertexIndex << ", ";
    }
    std::cout << std::endl;
  }

  for(auto child: node->getChildren()) {
    printNode(child);
  }
}

ClusterNode* HierarchyTree::getRoot() {
  return this->root;
}

ClusterNode* HierarchyTree::getMostSpecificCluster(size_t vertexIndex) {
  return findNode(root,  vertexIndex);
}

ClusterNode* HierarchyTree::findNode(ClusterNode* node, size_t vertexIndex) {
  for (auto index: node->getVertexIndexes()) {
    if (vertexIndex == index) {
      if (node->getChildren().size() > 0) {
        for (auto children: node->getChildren()) {
          auto foundNode = findNode(children, vertexIndex);
          if (foundNode != nullptr) {
            return foundNode;
          }
        }
        return node;
      } else {
        return node;
      }
    }
  }
  return nullptr;
}

void HierarchyTree::postProcessing() {
  int labelCount = -1;
  root->setClusterLabel(labelCount++);
  setPostProcessingLabel(root, labelCount);
}

void HierarchyTree::setPostProcessingLabel(ClusterNode* node, int &label) {
  for (auto child: node->getChildren()) {
    child->setClusterLabel(label++);;
  }

  for (auto child: node->getChildren()) {
    setPostProcessingLabel(child, label);
  }
}
}  // namespace datadriven
}  // namespace sgpp