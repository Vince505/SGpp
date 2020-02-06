// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include  <sgpp/datadriven/datamining/tools/Graph.hpp>
#include <vector>

namespace sgpp {
namespace datadriven {
/**
 * Class defining the nodes of the hierarchy tree. Each node represents a cluster in the given level
 */
class ClusterNode {
 public:
  /**
   * Default constructor
   */
  ClusterNode() = default;
  /**
   * Constructor
   * @param clusterLabel The label given to the node
   * @param vertexIndexes List of indexes indicating which points belong to this cluster
   * @param density The density threshold in which the cluster appeared
   */
  ClusterNode(int clusterLabel, std::vector<size_t> vertexIndexes, double density);

  /**
   * Default destructor
   */
  ~ClusterNode() = default;

  /**
   * Removes a child from the node
   * @param child THe reference pointing to the child to remove
   */
  void removeChild(ClusterNode* child);

  /**
   * Removes all of the children of this node
   */
  void removeChildren();

  /**
   * Adds a new child to this node
   * @param child Refrence to the child to be added
   */
  void addChild(ClusterNode* child);
  /**
   * Adds the a list of children to this node
   * @param children Children to be added as this new node's children
   */
  void addChildren(std::vector<ClusterNode*> children);

  /**
   *
   * @param graph
   * @param densityThreshold
   * @return
   */
  bool split(std::shared_ptr<Graph> graph, double densityThreshold);

  /**
   *
   * @param child
   * @param graph
   * @param densityThreshold
   * @return
   */
  bool splitChild(ClusterNode* child, std::shared_ptr<Graph>  graph, double densityThreshold);

  /**
   * Gets the list of indexes pointing to the poinst contained in this node
   * @return List of indexes
   */
  std::vector<size_t> getVertexIndexes();

  /**
   * Gets all of the children of this node
   * @return List containing the references to all of the children of thsi node
   */
  std::vector<ClusterNode*> getChildren();

  /**
   * Gets the parent of this node
   * @return Reference pointning to the parent of this node
   */
  ClusterNode* getParent();

  /**
   * Gets the cluster label of this node
   * @return The cluster label of this node
   */
  int getClusterLabel();

  /**
   * Replaces the cluster label of this node with a new one
   * @param clusterLabel The new label of this node
   */
  void setClusterLabel(int clusterLabel);

  /**
   * Sets the reference to the parent node
   * @param parent New reference to the parent node
   */
  void setParent(ClusterNode* parent);

  /**
   * Gets the density threshold in which this cluster was created
   * @return  density threshold at which this node was created
   */

  double getDensityThreshold();

  /**
   * Gets the level of the hierarchy tree in which this node was placed
   * @return Level of the hierarchy tree in which this node was placed
   */
  size_t getLevel();

  /**
   * Sets the level of this node
   * @param level The new level of this node
   */
  void setLevel(size_t level);

 private:
  /**
   * Number indicating the label of the cluster stored in this node. Note: -1 is always used for
   * noise
   */
  int clusterLabel;
  /**
   * Number indicating at which density threshold this node was created
   */
  double density;
  /**
   * Number indicating the level of the tree in which this node was placed
   */
  size_t level;
  /**
   * List containing the indexes of the points that belong to this cluster.
   * These indexes are the ones
   * given by the VpTree
   */
  std::vector<size_t> vertexIndexes;
  /**
   * List containing the references to the children of this node
   */
  std::vector<ClusterNode*> children;
  /**
   * Rerence to the parent of this node.
   */
  ClusterNode* parent;
};
}  // namespace datadriven
}  // namespace sgpp
