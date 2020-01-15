// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include  <sgpp/datadriven/tools/hierarchyTree/ClusterNode.hpp>

#include <algorithm>
#include <iostream>

namespace sgpp {
namespace datadriven {

ClusterNode::ClusterNode(int clusterLabel, std::vector<size_t> vertexIndexes, double density) {
  this->clusterLabel = clusterLabel;
  this->vertexIndexes = vertexIndexes;
  this->parent = nullptr;
  this->density = density;
}

void ClusterNode::removeChild(ClusterNode* child) {
  children.erase(std::remove(children.begin(), children.end(), child),children.end());
  delete child;
}

void ClusterNode::removeChildren() {

  for (auto child: children) {
    delete child;
  }
  children.clear();
}

ClusterNode* ClusterNode::getParent() {
  return parent;
}

void ClusterNode::addChildren(std::vector<ClusterNode*> children) {
  for (auto child: children) {
    this->addChild(child);
  }
}

bool ClusterNode::splitChild(ClusterNode* child, std::shared_ptr<Graph> graph,
  double densityThreshold) {

  std::cout << "Processing Node: " << this->getClusterLabel() << std::endl;
  std::cout << "Processing Child: " << child->getClusterLabel() << std::endl;

  std::vector<size_t> parentVertexIndexes = this->getVertexIndexes();
  std::vector<size_t> childVertexIndexes = child->getVertexIndexes();

  std::cout <<"Parent size: "<< parentVertexIndexes.size()<<std::endl;
  std::cout <<"Child size: "<< childVertexIndexes.size()<<std::endl;
  size_t maxConnectionsParent = (parentVertexIndexes.size()*(parentVertexIndexes.size()-1)/2);
  size_t maxConnectionsChildParent =
    childVertexIndexes.size() *(parentVertexIndexes.size() - childVertexIndexes.size());

  size_t connectionsParent = 0;
  size_t connectionsChildParent = 0;
  std::vector<size_t> visitedVertex;

  for(auto vertex: parentVertexIndexes) {
    if (graph->containsVertex(vertex)) {
      auto adjacentVertices = graph->getAdjacentVertices(vertex);
      for (auto adjacentVertex : adjacentVertices) {
        // Check that the adjacent vertex is in the parent node and checking that we have not
        // procesed it in a previous iteration
        if (std::find(visitedVertex.begin(), visitedVertex.end(), adjacentVertex)
            == visitedVertex.end()) {
          if (std::find(parentVertexIndexes.begin(), parentVertexIndexes.end(), adjacentVertex)
              != parentVertexIndexes.end()) {
            connectionsParent++;
          }
        }

      }
      visitedVertex.push_back(vertex);
    }
  }

  visitedVertex.clear();
  for(auto vertex: childVertexIndexes) {
    if (graph->containsVertex(vertex)) {
      auto adjacentVertices = graph->getAdjacentVertices(vertex);
      for (auto adjacentVertex : adjacentVertices) {
        if (std::find(visitedVertex.begin(), visitedVertex.end(), adjacentVertex)
            == visitedVertex.end()) {
          // Check that the adjacent vertex is in the parent node and checking that we have not
          // procesed it in a previous iteration
          if (std::find(parentVertexIndexes.begin(), parentVertexIndexes.end(), adjacentVertex)
              != parentVertexIndexes.end()) {
            connectionsChildParent++;
          }
        }
      }
      visitedVertex.push_back(vertex);
    }
  }
  std::cout << "Connections in parent " << connectionsParent << std::endl;
  std::cout <<" Max connections parent " << maxConnectionsParent << std::endl;
  std::cout << "Connections child parent " << connectionsChildParent << std::endl;
  std::cout <<" Max connections child parent " << maxConnectionsChildParent << std::endl;

  double connectivityParent = connectionsParent/ static_cast<double>(maxConnectionsParent);
  double connectivityChildParent = connectionsChildParent/static_cast<double>(maxConnectionsChildParent);

  std::cout << "Connectivty Parent" << connectivityParent << std::endl;
  std::cout << "Connectiviy Child Parent" << connectivityChildParent << std::endl;

  double compare = connectivityChildParent/connectivityParent;

  std::cout << "Final value " << compare << std::endl;
  return compare <= densityThreshold;
}

bool ClusterNode::split(std::shared_ptr<Graph>  graph,  double densityThreshold) {
  for(auto child:children) {
    if (splitChild(child, graph, densityThreshold)) {
      return true;
    }
  }
  return false;
}

std::vector<size_t> ClusterNode::getVertexIndexes() {
  return this->vertexIndexes;
}

void ClusterNode::setParent(ClusterNode* parent) {
  this->parent = parent;
}

void ClusterNode::addChild(ClusterNode* child) {
  this->children.push_back(child);

  child->setParent(this);
}

int ClusterNode::getClusterLabel() {
  return this->clusterLabel;
}

void ClusterNode::setClusterLabel(int clusterLabel) {
  this->clusterLabel = clusterLabel;
}

std::vector<ClusterNode*> ClusterNode::getChildren() {
  return this->children;
}

double ClusterNode::getDensityThreshold() {
  return this->density;
}

}  // namespace datadriven
}  // namespace sgpp