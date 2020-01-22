// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/tools/Graph.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/graph_traits.hpp>

#include <iostream>
#include <ctime>
#include <queue>
#include <map>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;

namespace sgpp {
namespace datadriven {

Graph::Graph(size_t vertices) {
  this->graph = new UndirectedGraph(vertices);

  boost::graph_traits<UndirectedGraph>::vertex_iterator vi, vend;
  size_t cnt = 0;
  for (boost::tie(vi, vend) = boost::vertices(*graph); vi != vend; ++vi) {
    indexToPointer[cnt] = *vi;
    pointerToIndex[*vi] = cnt;
    cnt++;
  }
  maxIndex = cnt;
}

void Graph::addVertex() {
  auto vertex = boost::add_vertex(*graph);
  indexToPointer[maxIndex] = vertex;
  pointerToIndex[vertex] = maxIndex;
  maxIndex++;
}

void Graph::removeVertex(size_t vertex) {
    boost::clear_vertex(indexToPointer[vertex], *graph);
    boost::remove_vertex(indexToPointer[vertex], *graph);
    pointerToIndex.erase(indexToPointer[vertex]);
    indexToPointer.erase(vertex);
    deletedVertices.push_back(vertex);
}

size_t Graph::getConnectedComponents(
    std::map<UndirectedGraph::vertex_descriptor, size_t> &componentMap) {

  fillIndexMap();
  auto numberComponents = boost::connected_components(*graph,
      boost::make_assoc_property_map(componentMap));

  return numberComponents;
}

void Graph::createEdges(size_t vertex, std::priority_queue<VpHeapItem> nearestNeighbors) {
  while (!nearestNeighbors.empty()) {
    addEdge(vertex, nearestNeighbors.top().index);
    nearestNeighbors.pop();
  }
}

void Graph::addEdge(size_t vertex1, size_t vertex2) {
  boost::add_edge(indexToPointer[vertex1], indexToPointer[vertex2], *graph);
}

void Graph::deleteEdge(size_t vertex1, size_t vertex2) {
  boost::remove_edge(indexToPointer[vertex1], indexToPointer[vertex2], *graph);
}

void Graph::fillIndexMap() {
  /*
  * This index has to be created in order for the connected components to work
  * while using lists as a vertex conatiner
  */
  boost::property_map<UndirectedGraph, boost::vertex_index_t>::type
    index = get(boost::vertex_index, *graph);

  boost::graph_traits<UndirectedGraph>::vertex_iterator vi, vend;
  size_t cnt = 0;
  for (boost::tie(vi, vend) = boost::vertices(*graph); vi != vend; ++vi) {
    boost::put(index, *vi, cnt++);
  }
}

std::vector<size_t> Graph::getAdjacentVertices(size_t vertex) {
  fillIndexMap();
  std::vector<size_t> indexes;
  boost::graph_traits<UndirectedGraph>::vertex_iterator begin, end;
  // Getting the Iterator
  auto neighborsIterator = boost::adjacent_vertices(indexToPointer[vertex], *graph);
  for (auto vd : boost::make_iterator_range(neighborsIterator)) {
    // Translating to the indexes
    if(pointerToIndex.find(vd) != pointerToIndex.end()) {
      indexes.push_back(pointerToIndex[vd]);
    }
  }

  return indexes;
}

size_t Graph::getIndex(UndirectedGraph::vertex_descriptor vertexDescriptor) {
  return pointerToIndex[vertexDescriptor];
}

UndirectedGraph::vertex_descriptor Graph::getVertexDescriptor(size_t vertex) {
  return indexToPointer[vertex];
}

bool Graph::containsVertex(size_t vertex) {
  if (indexToPointer.find(vertex) != indexToPointer.end()) {
    return true;
  } else {
    return false;
  }
}

UndirectedGraph* Graph::getGraph() {
  return graph;
}
}  // namespace datadriven
}  // namespace sgpp