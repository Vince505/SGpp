// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/tools/vpTree/VpNode.hpp>
#include <sgpp/datadriven/tools/vpTree/VpHeapItem.hpp>

#include <vector>
#include<utility>
#include <queue>
#include <iostream>
namespace sgpp {
namespace datadriven {
/**
 * Class for the Vantage Point Tree.
 * Based on the code by Steve Hanov's tutorial
 * at http://stevehanov.ca/blog/index.php?id=130
 */
class VpTree {
 public:
  /**
   * Constructor
   * @param matrix
   * Matrix of points to use
   */
  explicit VpTree(sgpp::base::DataMatrix matrix);

  // Destructor
  ~VpTree() {
    delete root;
  }

  /**
   * Gets the nearest neighbors of a given points
   * @param target The point whose nearest neighbors are to be searched
   * @param noNearestNeighbors The number of nearest neighbors to deliver
   * @return Priority queue with the index and coordinates of the neearest neighbors
   */
  std::priority_queue<VpHeapItem>  getNearestNeighbors(sgpp::base::DataVector &target,
      size_t noNearestNeighbors);

  /**
   * Updates the vpTree with a new set of points
   * @param matrix Matrix with the new points
   */
  void update(sgpp::base::DataMatrix &matrix);

  /**
   * Gets all of the points inside the vptTree
   * @return Matrix with the coordinates of all the points inside the vptTree
   */
  sgpp::base::DataMatrix &getStoredItems();

  /**
  * Defines the euclidean distance metric between two points
  */
  static double euclideanDistance(sgpp::base::DataVector point1, sgpp::base::DataVector point2);

  /**
   * Obtains the index given by the tree to a point
   * @param point Point whose index is being looked for
   * @return Index of the point
   */
  size_t getIndexedKeyFromPoint(sgpp::base::DataVector &point);

  /**
   * Prints the vpTRee in a preorder fashion
   */
  void printPreorder();

 private:
  /**
   * Node that defines the root of the tree
   */
  VpNode* root;
  /**
   * Value that keeps track the distance of the latest found nearest neighbor
   */
  double tau;

  /**
   * Matrix to keep track of all stored points when building the tree
   */
  sgpp::base::DataMatrix storedItems;

  /**
   * Method that builds a VpTree recursevily
   * @param startIndex Index of the matrix of the points to store
   * @param endIndex End index of the matrix of the poinst to store
   * @return The node corresponding to the roort of the tree
   */
  VpNode* buildRecursively(size_t startIndex, size_t endIndex);

  /**
   * Inserts new node in the tree
   * @param indexNewPoint THe index of the new point to be added in the node
   */
  void insertNewNode(size_t indexNewPoint);

  /**
   * Realizes a recursive search to find the node to append a new one
   * @param node Node to be processed
   * @param newPoint The coordinates of the point to be inserted in the tree
   * @return The reference to the node where the new node will be appended
   */
  VpNode* findInsertionNode(VpNode* node, sgpp::base::DataVector &newPoint);

  /**
   * Does a recursive search to find the index of a given point
   * @param node Node to be processed
   * @param point Point whose index is being looked for
   * @return The index of the given point
   */
  size_t findIndex(VpNode* node, sgpp::base::DataVector &point);

  /**
   * Recursive search to find the nearest neighbors of a point
   * @param node Node being processed
   * @param target Point whose nearest neighbors are being searched
   * @param noNearestNeighbors Number of nearest neighbors to deliver
   * @param heap The heap contaning the info of the currently found nearest neighbors
   */
  void searchRecursively(VpNode* &node, sgpp::base::DataVector &target,
      size_t noNearestNeighbors, std::priority_queue<VpHeapItem> &heap);

  /**
   * Swaps to rows in the storage matrix
   * @param index1 Index of the point to be swapped for the one given by index2
   * @param index2 Index of the point to be swapped for the one given by index1
   */
  void swap(size_t index1, size_t index2);

  /**
   * Calculates the distances of the points in the matrix to the vantage point whose indexes lie
   * between the startIndex and the endIndex
   * @param startIndex Index indicating the first row in the matrix to be processed
   * @param endIndex Indicating the last row in the matrix to be processed
   * @return A vetor containing pairs mapping the index of the point processed and its distance to
   * the vantage point
   */
  std::vector<std::pair <size_t, double>>  getDistances(size_t startIndex, size_t endIndex);

  /**
   * Sorts the rows of the storage matrix, whose indexes lie in the range given by the index
   * based on their distance to the vantage point
   * @param index1 indicating the first row in the matrix to be processed
   * @param index2 Indicating the last row in the matrix to be processed
   */
  void sortByDistances(size_t index1, size_t index2);

  /**
   * Prints recursively the information of the tree
   * @param node Node being processed
   */
  void printPreorder(VpNode* node);
};
}  // namespace datadriven
}  // namespace sgpp
