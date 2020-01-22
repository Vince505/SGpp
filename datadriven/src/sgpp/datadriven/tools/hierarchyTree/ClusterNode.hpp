// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include  <sgpp/datadriven/tools/Graph.hpp>
#include <vector>

namespace sgpp {
namespace datadriven {

class ClusterNode {

public:

  ClusterNode() = default;

  ClusterNode(int clusterLabel, std::vector<size_t> vertexIndexes, double density);

  ~ClusterNode() = default;
  void removeChild(ClusterNode* child);

  void removeChildren();

  void addChildren(std::vector<ClusterNode*> children);

  bool split(std::shared_ptr<Graph> graph, double densityThreshold);

  bool splitChild(ClusterNode* child, std::shared_ptr<Graph>  graph, double densityThreshold);

  std::vector<size_t> getVertexIndexes();

  std::vector<ClusterNode*> getChildren();

  ClusterNode* getParent();

  int getClusterLabel();

  void setClusterLabel(int clusterLabel);

  void setParent(ClusterNode* parent);

  void addChild(ClusterNode* child);

  double getDensityThreshold();

  size_t getLevel();

  void setLevel(size_t level);

private:
  int clusterLabel;
  double density;
  size_t level;
  std::vector<size_t> vertexIndexes;
  std::vector<ClusterNode*> children;
  ClusterNode* parent;
};

}  // namespace datadriven
}  // namespace sgpp