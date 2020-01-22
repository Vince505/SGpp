// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <sgpp/datadriven/tools/vpTree/VpTree.hpp>

#include <map>
#include <queue>
#include <iostream>

namespace sgpp {
namespace datadriven {
/**
 * Definition of the type UndirecetdGraph
 */
typedef boost::adjacency_list<boost::setS, boost::listS, boost::undirectedS,
boost::property<boost::vertex_index_t, size_t>, boost::no_property> UndirectedGraph;

/**
* @brief Class that encapsulates all methods and properties of an unidrected Graph
*/
class Graph {
 public:
 /**
  * Constructs with no edges and a certan number of vertices
  * @param vertices Number of vertices contained in the graph
  */
  explicit Graph(size_t vertices);

  /**
   * Copy constructor
   * @param rhs copy object
   */
  Graph(const Graph &rhs) {
    this->graph = new UndirectedGraph(*(rhs.graph));
    this->deletedVertices = rhs.deletedVertices;
    this->maxIndex = rhs.maxIndex;
    boost::graph_traits<UndirectedGraph>::vertex_iterator vi, vend;
    boost::tie(vi, vend) = boost::vertices(*(this->graph));

    for (size_t cnt = 0; cnt <= this->maxIndex; cnt++) {
      if (std::find(this->deletedVertices.begin(), this->deletedVertices.end(), cnt)
          == this->deletedVertices.end()) {
        this->indexToPointer[cnt] = *vi;
        this->pointerToIndex[*vi] = cnt;
        ++vi;
      } else {
        cnt++;
      }
    }
  }

  /**
   * Destructor
   */
  ~Graph() {
    delete this->graph;
  };

  /**
   * Method to add an additional vertex to the graph
   */
  void addVertex();

  /**
   * Removes a given vertex
   * @param vertex The index used to identify the vertex to remove
   */
  void removeVertex(size_t vertex);

  /**
   * Creates the edges of a vertex given its nearest neighbors
   * @param vertex The index used to identify the vertex to remove
   * @param nearestNeighbors priority queue obtaiend from a VP Tree with the indexes of the vertex
   * which are the nearest neighbors
   */
  void createEdges(size_t vertex, std::priority_queue<VpHeapItem> nearestNeighbors);

  /**
   * Adds an edge between to vertices
   * @param vertex1 The index used to identify the source vertex
   * @param vertex2 The index used to identify the sink vertex
   */
  void addEdge(size_t vertex1, size_t vertex2);

  /**
   * Deletes and edge between to vertices
   * @param vertex1 The index used to identify the source vertex
   * @param vertex2 The index used to identify the sink vertex
   */
  void deleteEdge(size_t vertex1, size_t vertex2);

  /**
   * Gets the underlying boost graph structure
   * @return Pointer to the boost graph data structure
   */
  UndirectedGraph* getGraph();

  /**
   * Obtains the number of connected components and stores the labels in a given map
   * @param componentMap Map which contains the mapping from vertices to assigned labels
   * @return  Number of connected components detected
   */
  size_t getConnectedComponents(std::map<UndirectedGraph::vertex_descriptor, size_t> &componentMap);

  /**
   * Obtains the indexes of the vertices connected to the given vertex
   * @param vertex The index used to identify the vertex
   * @return A list of the indexes of the vertices connceted to the given vertex
   */
  std::vector<size_t> getAdjacentVertices(size_t vertex);

  size_t getIndex (UndirectedGraph::vertex_descriptor vertexDescriptor);

  UndirectedGraph::vertex_descriptor getVertexDescriptor(size_t vertex);

  bool containsVertex(size_t vertex);

 private:

    void fillIndexMap();
    /**
     * Boost graph data structure
     */
    UndirectedGraph* graph;

    std::map<UndirectedGraph::vertex_descriptor, size_t> pointerToIndex;

    std::map<size_t, UndirectedGraph::vertex_descriptor> indexToPointer;

    std::vector<size_t> deletedVertices;

    size_t maxIndex;

};
}  // namespace datadriven
}  // namspace sgpp
