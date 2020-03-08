// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org


#include  <sgpp/datadriven/datamining/tools/hierarchyTree/HierarchyTree.hpp>

#include <iostream>
#include <vector>
#include <string>

using json::JSON;

using sgpp::base::DataVector;
namespace sgpp {
namespace datadriven {

HierarchyTree::HierarchyTree(size_t numberPoints) {
  std::vector<size_t> allVertexIndexes;
  for (size_t index = 0; index < numberPoints; index++) {
    allVertexIndexes.push_back(index);
  }
  this->root = new ClusterNode(-1, allVertexIndexes, 0.0);
}

void HierarchyTree::printTree() {
  auto cluster = root->getVertexIndexes();
  std::cout << "=======================HIERACHY TREE=========================" << std::endl;
  std::cout << "Cluster Label: " << root->getClusterLabel() << std::endl;
  std::cout << "Level: " << root->getLevel() << std::endl;
  std::cout << "Indexes: ";
  for (auto vertexIndex : cluster) {
    std::cout << vertexIndex << ", ";
  }
  std::cout << std::endl;

  printNode(root);
}
void HierarchyTree::printNode(ClusterNode* node) {
  for (auto child : node->getChildren()) {
    auto cluster = child->getVertexIndexes();
    std::cout << "================================================"<< std::endl;
    std::cout << "Parent: " << child->getParent()->getClusterLabel();
    std::cout << std::endl;
    std::cout << "Cluster Label: " << child->getClusterLabel() << std::endl;
    std::cout << "Level: " << child->getLevel() << std::endl;
    std::cout << "Indexes: ";
    for (auto vertexIndex : cluster) {
      std::cout << vertexIndex << ", ";
    }
    std::cout << std::endl;
  }

  for (auto child : node->getChildren()) {
    printNode(child);
  }
}

ClusterNode* HierarchyTree::getRoot() {
  return this->root;
}

size_t HierarchyTree::getMostSpecificLevel(size_t vertexIndex) {
  return findNode(root,  vertexIndex)->getLevel();
}
ClusterNode* HierarchyTree::getMostSpecificCluster(size_t vertexIndex) {
  return findNode(root,  vertexIndex);
}

ClusterNode* HierarchyTree::findNode(ClusterNode* node, size_t vertexIndex) {
  for (auto index : node->getVertexIndexes()) {
    if (vertexIndex == index) {
      if (node->getChildren().size() > 0) {
        for (auto children : node->getChildren()) {
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

ClusterNode* HierarchyTree::getMostSpecificClusterAtLevel(size_t vertexIndex, size_t level) {
  return findNodeAtLevel(root,  vertexIndex, level);
}

ClusterNode* HierarchyTree::findNodeAtLevel(ClusterNode* node, size_t vertexIndex, size_t level) {
  if (node->getLevel() > level) {
    return node->getParent();
  }
  for (auto index : node->getVertexIndexes()) {
    if (vertexIndex == index) {
      if (node->getChildren().size() > 0) {
        for (auto children : node->getChildren()) {
          auto foundNode = findNodeAtLevel(children, vertexIndex, level);
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

size_t HierarchyTree::getNumberClusters() {
  size_t count = 0;
  countCluster(root, count);

  return count;
}

size_t HierarchyTree::getNumberLevels() {
  size_t maxLevel = 0;
  countLevel(root, maxLevel);

  return maxLevel;
}
void HierarchyTree::countCluster(ClusterNode* node, size_t &numberClusters) {
  numberClusters += node->getChildren().size();

  for (auto child : node->getChildren()) {
    countCluster(child, numberClusters);
  }
}

void HierarchyTree::countLevel(ClusterNode* node, size_t &maxLevel) {
  if (node->getLevel() > maxLevel) {
    maxLevel = node->getLevel();
  }

  for (auto child : node->getChildren()) {
    countLevel(child, maxLevel);
  }
}

void HierarchyTree::postProcessing() {
  int labelCount = -1;
  root->setClusterLabel(labelCount++);
  setPostProcessingLabel(root, labelCount);
}

void HierarchyTree::setPostProcessingLabel(ClusterNode* node, int &label) {
  for (auto child : node->getChildren()) {
    child->setClusterLabel(label++);;
  }

  for (auto child : node->getChildren()) {
    setPostProcessingLabel(child, label);
  }
}

void HierarchyTree::evaluateClustering(DataVector &results) {
  results.resize(root->getVertexIndexes().size());
  for (auto vertex : root->getVertexIndexes()) {
    results.set(vertex, getMostSpecificCluster(vertex)->getClusterLabel());
  }
}

void HierarchyTree::evaluateClusteringAtLevel(DataVector &results, size_t level) {
  results.resize(root->getVertexIndexes().size());
  for (auto vertex : root->getVertexIndexes()) {
    results.set(vertex, getMostSpecificClusterAtLevel(vertex, level)->getClusterLabel());
  }
}

void HierarchyTree::storeHierarchy(std::string outputDirectory) {
  JSON output;
  size_t n_levels = getNumberLevels();
  for (size_t level = 1; level <= n_levels; level++) {
    output.addDictAttr("Level "+std::to_string(level));
  }
  storeHierarchy(root, output);

  std::cout << "Stroring hierarchy info at " << outputDirectory+"/hierarchyInfo.json" << std::endl;
  output.serialize(outputDirectory+"/hierarchyInfo.json");
}

void HierarchyTree::storeHierarchy(ClusterNode* node, json::JSON &output) {
  if (node->getClusterLabel() != -1) {  // We skip the root or the noise cluster
    output["Level "+std::to_string(node->getLevel())]
    .addDictAttr("Cluster " + std::to_string(node->getClusterLabel()));

    output["Level "+std::to_string(node->getLevel())]
    ["Cluster " + std::to_string(node->getClusterLabel())]
      .addIDAttr("Density Threshold", node->getDensityThreshold());

    output["Level "+std::to_string(node->getLevel())]
    ["Cluster " + std::to_string(node->getClusterLabel())]
      .addIDAttr("No. points", node->getVertexIndexes().size());

    output["Level "+std::to_string(node->getLevel())]
    ["Cluster " + std::to_string(node->getClusterLabel())].addListAttr("Children");

    if (node->getChildren().size() > 0) {
      for (auto child : node->getChildren()) {
        output["Level "+std::to_string(node->getLevel())]["Cluster " +
               std::to_string(node->getClusterLabel())]["Children"]
          .addIdValue("\"Cluster "+std::to_string(child->getClusterLabel())+"\"");
      }
    }
  }
  for (auto child : node->getChildren()) {
    storeHierarchy(child, output);
  }
}
}  // namespace datadriven
}  // namespace sgpp
