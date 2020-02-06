// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/tools/Graph.hpp>
#include <boost/test/unit_test.hpp>

#include <vector>
#include <map>
#include <algorithm>

using sgpp::datadriven::Graph;
using sgpp::datadriven::UndirectedGraph;

BOOST_AUTO_TEST_SUITE(test_graph)

Graph testGraph(11);

BOOST_AUTO_TEST_CASE(graphCreation) {
  // Test all vertices are in the graph
  BOOST_CHECK_EQUAL(testGraph.getNumberVertices(), 11);

  // Tests that edges already created are not duplicated when giving the vertices in inverse order
  testGraph.addEdge(0, 1);
  testGraph.addEdge(0, 2);
  testGraph.addEdge(0, 3);
  testGraph.addEdge(1, 0);
  testGraph.addEdge(1, 2);
  testGraph.addEdge(2, 0);
  testGraph.addEdge(2, 1);
  testGraph.addEdge(2, 4);
  testGraph.addEdge(2, 5);
  testGraph.addEdge(3, 0);
  testGraph.addEdge(3, 4);
  testGraph.addEdge(4, 2);
  testGraph.addEdge(4, 3);
  testGraph.addEdge(5, 2);
  testGraph.addEdge(5, 6);
  testGraph.addEdge(6, 5);
  testGraph.addEdge(6, 7);
  testGraph.addEdge(6, 10);
  testGraph.addEdge(7, 6);
  testGraph.addEdge(7, 8);
  testGraph.addEdge(8, 9);
  testGraph.addEdge(8, 10);
  testGraph.addEdge(9, 8);
  testGraph.addEdge(10, 6);
  testGraph.addEdge(10, 8);

  BOOST_CHECK_EQUAL(testGraph.getNumberEdges(), 13);
}

BOOST_AUTO_TEST_CASE(adjacentVertices) {
  // Tests that the adjacent vertices function works corrrectly
  std::vector<size_t> realAdjacentVertices = {0, 1, 4, 5};

  auto predictedAdjacentVertices = testGraph.getAdjacentVertices(2);
  std::sort(predictedAdjacentVertices.begin(), predictedAdjacentVertices.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(realAdjacentVertices.begin(), realAdjacentVertices.end(),
  predictedAdjacentVertices.begin(), predictedAdjacentVertices.end());
}

BOOST_AUTO_TEST_CASE(connectedComponents) {
  // Tests that the number of connected components is correct
  std::map<UndirectedGraph::vertex_descriptor, size_t> clusterMap;
  BOOST_CHECK_EQUAL(testGraph.getConnectedComponents(clusterMap), 1);
}

BOOST_AUTO_TEST_CASE(deletion) {
  testGraph.removeVertex(5);

  std::vector<size_t> realAdjacentVertices = {0, 1, 4};
  std::map<UndirectedGraph::vertex_descriptor, size_t> clusterMap;

  // After deletion vertex must not be there anymore
  BOOST_CHECK_EQUAL(testGraph.containsVertex(5), false);

  auto predictedAdjacentVertices = testGraph.getAdjacentVertices(2);
  std::sort(predictedAdjacentVertices.begin(), predictedAdjacentVertices.end());
  // 5 must not appear anymore as the adjacent vertex of 2
  BOOST_CHECK_EQUAL_COLLECTIONS(realAdjacentVertices.begin(), realAdjacentVertices.end(),
  predictedAdjacentVertices.begin(), predictedAdjacentVertices.end());

  // Connected componentes increase since the vertex 5 was previously linking to subgraphs
  BOOST_CHECK_EQUAL(testGraph.getConnectedComponents(clusterMap), 2);

  // Number of vertices and edges must be updated accordingly
  BOOST_CHECK_EQUAL(testGraph.getNumberVertices(), 10);
  BOOST_CHECK_EQUAL(testGraph.getNumberEdges(), 11);
}

BOOST_AUTO_TEST_CASE(copy) {
  Graph copyGraph = testGraph;

  // Both number of vertices and edges must be the same
  BOOST_CHECK_EQUAL(testGraph.getNumberEdges(), 11);
  BOOST_CHECK_EQUAL(copyGraph.getNumberEdges(), 11);

  // Both adjacent vertices of the same vertex must be the same
  std::vector<size_t> realAdjacentVertices = {6, 8};
  auto predictedAdjacentVertices = testGraph.getAdjacentVertices(7);
  auto copyAdjacentVertices = copyGraph.getAdjacentVertices(7);

  std::sort(predictedAdjacentVertices.begin(), predictedAdjacentVertices.end());
  std::sort(copyAdjacentVertices.begin(), copyAdjacentVertices.end());

  BOOST_CHECK_EQUAL_COLLECTIONS(realAdjacentVertices.begin(), realAdjacentVertices.end(),
  predictedAdjacentVertices.begin(), predictedAdjacentVertices.end());

  BOOST_CHECK_EQUAL_COLLECTIONS(realAdjacentVertices.begin(),realAdjacentVertices.end(),
  copyAdjacentVertices.begin(), copyAdjacentVertices.end());

  // Check if the vertex we removed in the previous test is removed for both graphs
  BOOST_CHECK_EQUAL(testGraph.containsVertex(5), false);
  BOOST_CHECK_EQUAL(copyGraph.containsVertex(5), false);

  // Number of components must be the same
  std::map<UndirectedGraph::vertex_descriptor, size_t> clusterMap;
  BOOST_CHECK_EQUAL(testGraph.getConnectedComponents(clusterMap), 2);
  BOOST_CHECK_EQUAL(copyGraph.getConnectedComponents(clusterMap), 2);

  // Test. When removing a vertex only one of the graphs is affected.
  copyGraph.removeVertex(8);

  // The original graph should still contain the vertex 8, the copy must not.
  BOOST_CHECK_EQUAL(testGraph.containsVertex(8), true);
  BOOST_CHECK_EQUAL(copyGraph.containsVertex(8), false);

  // Copy must have less vertices and edges
  BOOST_CHECK_EQUAL(testGraph.getNumberVertices(), 10);
  BOOST_CHECK_EQUAL(copyGraph.getNumberVertices(), 9);

  BOOST_CHECK_EQUAL(testGraph.getNumberEdges(), 11);
  BOOST_CHECK_EQUAL(copyGraph.getNumberEdges(), 8);

  // Vertex 8 does not appear anymore as an adjacent vertex in the copy graph.
  std::vector<size_t> copyRealAdjacentVertices = {6};
  copyAdjacentVertices = copyGraph.getAdjacentVertices(7);
  BOOST_CHECK_EQUAL_COLLECTIONS(realAdjacentVertices.begin(), realAdjacentVertices.end(),
  predictedAdjacentVertices.begin(), predictedAdjacentVertices.end());

  BOOST_CHECK_EQUAL_COLLECTIONS(copyRealAdjacentVertices.begin(), copyRealAdjacentVertices.end(),
  copyAdjacentVertices.begin(), copyAdjacentVertices.end());


  // Copy graph must have more connected components
  BOOST_CHECK_EQUAL(testGraph.getConnectedComponents(clusterMap), 2);
  BOOST_CHECK_EQUAL(copyGraph.getConnectedComponents(clusterMap), 3);
}


BOOST_AUTO_TEST_SUITE_END()
