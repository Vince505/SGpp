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
  /**
   * Construnctor
   * @param numberPoints Number of points used for the clustering
   */
  explicit HierarchyTree(size_t numberPoints);

  /**
   * Destructor
   */
  ~HierarchyTree() {
    delete root;
  }

  /**
   * Prints the tree on the console for debugiing purposes
   */
  void printTree();

  /**
   * Gets the root of the tree
   * @return root of the tree
   */
  ClusterNode* getRoot();

  /**
   * Gets the deepest node (max level in the tree) in which a point is found
   * @param vertexIndex The index of the point
   * @return The deepest node in which a point is found
   */
  ClusterNode* getMostSpecificCluster(size_t vertexIndex);

  /**
   * Gets the deepest node (max level in the tree) in which a point is found without exceeding
   * a certain level.
   * @param vertexIndex The index of the point
   * @param level Maximum level to look into
   * @return THe depest node at the given level in which the point is found
   */
  ClusterNode* getMostSpecificClusterAtLevel(size_t vertexIndex, size_t level);

  /**
   * Gets the number of clusters in the tree
   * @return number of clusters
   */
  size_t getNumberClusters();

  /**
   * Gets the number of levels in the tree
   * @return total number of levels in the tree
   */
  size_t getNumberLevels();

  /**
   * Renames the labels of the nodes so that they follow a sequential order
   */
  void postProcessing();

  /**
   * Gets the labels at the maximum possible level of all the points
   * @param results list to store the labels of the points at the maximum possible level
   */
  void evaluateClustering(sgpp::base::DataVector &results);

  /**
   * Gets the labels of the points up to a certain level
   * @param results list to store the labels of the point up to a certain level
   * @param level Level limiting th evaluation
   */
  void evaluateClusteringAtLevel(sgpp::base::DataVector &results, size_t level);

  /**
   * Gets the most specific level in which a point is found in the tree
   * @param vertexIndex Index of the point
   * @return THe most specific level in which the point is found in the tree
   */
  size_t getMostSpecificLevel(size_t vertexIndex);
private:
  /**
   * Root of the tree
   */
  ClusterNode* root;

  /**
   * Does a recursive search to find the deepest node containing a point
   * @param node Node to process
   * @param vertexIndex Index of the point
   * @return Reference to the deepest node containing the point
   */
  ClusterNode* findNode(ClusterNode* node, size_t vertexIndex);

  /**
   * Does a recursive search to find the deepest node containing a point until a certain level
   * @param node Node to process
   * @param vertexIndex Index of the point
   * @param level Maximum level to look into
   * @return Reference to the deepest node containing the point at the given level
   */
  ClusterNode* findNodeAtLevel(ClusterNode* node, size_t vertexIndex, size_t level);

  /**
   * Recursive method to count the number of clusters in the tree
   * @param node Node to process
   * @param numberClusters Variable to keept track of the count of clusters
   */
  void countCluster(ClusterNode* node, size_t &numberClusters);

  /**
   * Recursive method to count the number of levels in the tree
   * @param node current node being processed
   * @param maxLevel Max level found so far.
   */
  void countLevel(ClusterNode* node, size_t &maxLevel);

  /**
   * Recursive method to print a node
   * @param node Node to process
   */
  void printNode(ClusterNode* node);

  /**
   * Recursive method to set the label in a node during the postprocessing of the tree
   * @param node Node to process
   * @param clusterLabel Label to assign
   */
  void setPostProcessingLabel(ClusterNode* node, int &clusterLabel);
};

}  // namespace datadriven
}  // namespace sgpp
