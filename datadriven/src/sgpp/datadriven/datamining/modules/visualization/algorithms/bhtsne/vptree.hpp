/*
 *
 * Copyright (c) 2014, Laurens van der Maaten (Delft University of Technology)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. All advertising materials mentioning features or use of this software
 *    must display the following acknowledgement:
 *    This product includes software developed by the Delft University of Technology.
 * 4. Neither the name of the Delft University of Technology nor the names of 
 *    its contributors may be used to endorse or promote products derived from 
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY LAURENS VAN DER MAATEN ''AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO 
 * EVENT SHALL LAURENS VAN DER MAATEN BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR 
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING 
 * IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY 
 * OF SUCH DAMAGE.
 *
 */


/* This code was adopted with minor modifications from Steve Hanov's great tutorial
 * at http://stevehanov.ca/blog/index.php?id=130 */

/**
 * Code originally taken from https://lvdmaaten.github.io/tsne/
 * It has been modified in order to be adapted to the SG++ datamining
 * pipeline structure and has been parallelized
 */

#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include <queue>
#include <limits>
#include <cmath>
#include <cfloat>


namespace sgpp {
namespace datadriven {

/**
 * Reference class to the define DataPoint structures for the Vp-Tree for the t-SNE algorithm
 */
class DataPoint {
  /**
   * Index idetifiying the point
   */
  size_t _ind;

 public:
  double* _x;
  size_t _D;
  /**
   * Constructor
   */
  DataPoint() {
    _D = 1;
    _ind = -1;
    _x = NULL;
  }
  /**
   * Constructor
   * @param D Dimensionality of the data
   * @param ind Index idetifiying the point
   * @param x Array containing all of the coordinates of a given point
   */
  DataPoint(size_t D, int ind, double* x) {
    _D = D;
    _ind = ind;
    _x = new double[_D];
    for (size_t d = 0; d < _D; d++) {
      _x[d] = x[d];
    }
  }
  /**
   * Move copy operator
   * @param other copy
   */
  DataPoint(const DataPoint& other) {  // this makes a deep copy -- should not free anything
    if (this != &other) {
      _D = other.dimensionality();
      _ind = other.index();
      _x = new double[_D];
      for (size_t d = 0; d < _D; d++) {
        _x[d] = other.x(d);
      }
    }
  }
  /**
   * Destructor
   */
  ~DataPoint() {
    if (_x != NULL) {
      delete[ ]_x;
    }
  }
  /**
   * Move operator
   * @param other old object
   * @return
   */
  DataPoint& operator= (const DataPoint& other) {         // asignment should free old object
    if (this != &other) {
      if (_x != NULL) {
        free(_x);
      }
      _D = other.dimensionality();
      _ind = other.index();
      _x = new double[_D];
      for (size_t d = 0; d < _D; d++) {
        _x[d] = other.x(d);
      }
    }
    return *this;
  }
  /**
   * Gets index of the point
   * @return Index of a point
   */
  size_t index() const {
    return _ind;
  }
  /**
   * Gets dimensionality of the point
   * @return Dimensionality of the point
   */
  size_t dimensionality() const {
    return _D;
  }
  /**
   * Returns the dth coordinate of a point
   * @param d The position of the coordinate
   * @return The dth coordinate of a point
   */
  double x(size_t d) const {
    return _x[d];
  }

  /**
   * Calculates the euclidiean distance between 2 datapoints
   * @param t1 Point 1
   * @param t2 Point 2
   * @return Euclidian distance of 2 Datapoints
   */
  static double euclidean_distance(const DataPoint &t1, const DataPoint &t2) {
    double dd = .0;
    double* x1 = t1._x;
    double* x2 = t2._x;
    double diff;
    for (size_t d = 0; d < t1._D; d++) {
        diff = (x1[d] - x2[d]);
        dd += diff * diff;
    }
    return sqrt(dd);
  }
};


template<typename T, double (*distance)( const T&, const T& )>
class VpTree {
 public:
  // Default constructor
  VpTree() : _root(0) {}

  // Destructor
  ~VpTree() {
      delete _root;
  }

  /**
   * Creates a vpTree based on an item lists
   * @param items Structure containing the elments used to build the vpTree
   */
  void create(const std::vector<T>& items) {
    delete _root;
    _items = items;
    _root = buildFromPoints(0, items.size());
  }

  /**
   * Search the nearest neighbors of a given targe
   * @param target Target whose nearest neighbors are to be found
   * @param k The number of nearest neighbors to seek
   * @param results Vector containing the nearest neighbors
   * @param distances Vector containing the distances to the nearest neighbors
   */
  void search(const T& target, size_t k, std::vector<T>* results,
    std::vector<double>* distances) {
    // Use a priority queue to store intermediate results on
    std::priority_queue<HeapItem> heap;

    // Variable that tracks the distance to the farthest point in our results
    _tau = DBL_MAX;

    // Perform the search
    search(_root, target, k, heap);

    // Gather final results
    results->clear(); distances->clear();
    while (!heap.empty()) {
        results->push_back(_items[heap.top().index]);
        distances->push_back(heap.top().dist);
        heap.pop();
    }

    // Results are in reverse order
    std::reverse(results->begin(), results->end());
    std::reverse(distances->begin(), distances->end());
  }

 private:
  /**
   * Vector containing the items of the VpTree
   */
  std::vector<T> _items;

  /**
   * Variable to store the distance of the currently found farthest nearest neighbor
   */
  double _tau;

  // Single node of a VP tree (has a point and radius;
  // left children are closer to point than the radius)
  struct Node {
    size_t index;              // index of point in node
    double threshold;       // radius(?)
    Node* left;             // points closer by than threshold
    Node* right;            // points farther away than threshold

    /**
     * Default constructor
     */
    Node() :
    index(0), threshold(0.), left(0), right(0) {}

    /**
     * Destructor
     */
    ~Node() {
      delete left;
      delete right;
    }
  }
  /**
   * Root of the tree
   */
  * _root;


  // An item on the intermediate result queue
  struct HeapItem {
    HeapItem(size_t index, double dist) :
    index(index), dist(dist) {}
    size_t index;
    double dist;
    bool operator<(const HeapItem& o) const {
        return dist < o.dist;
    }
  };

  // Distance comparator for use in std::nth_element
  struct DistanceComparator {
    const T& item;
    explicit DistanceComparator(const T& item) : item(item) {}
    bool operator()(const T& a, const T& b) {
        return distance(item, a) < distance(item, b);
    }
  };

  /**
   *
   * @param lower Start index of the vector of the elements to store
   * @param upper End index of the matrix of the elements to store
   * @return Reference to the created node
   */
  Node* buildFromPoints(size_t lower, size_t upper) {
    if (upper == lower) {     // indicates that we're done here!
      return NULL;
    }

    // Lower index is center of current node
    Node* node = new Node();
    node->index = lower;

    if (upper - lower > 1) {      // if we did not arrive at leaf yet
      // Choose an arbitrary point and move it to the start
      int i = static_cast<int> ((static_cast<double>(rand()) /
        static_cast<double>(RAND_MAX * (upper - lower - 1)) +
        static_cast<double>(lower)));
      std::swap(_items[lower], _items[i]);

      // Partition around the median distance
      size_t median = (upper + lower) / 2;
      std::nth_element(_items.begin() + lower + 1,
                       _items.begin() + median,
                       _items.begin() + upper,
                       DistanceComparator(_items[lower]));

      // Threshold of the new node will be the distance to the median
      node->threshold = distance(_items[lower], _items[median]);

      // Recursively build tree
      node->index = lower;
      node->left = buildFromPoints(lower + 1, median);
      node->right = buildFromPoints(median, upper);
    }

    // Return result
    return node;
  }

  /**
   * Searchs recursively for the k neearest neighbors
   * @param node Node being processed
   * @param target Target whose nearest neighbors are to be found
   * @param k The number of nearest neighbors to seek
   * @param heap Priority queue which keeps track of the currently found nearest neighbors
   */
  void search(Node* node, const T& target, size_t k, std::priority_queue<HeapItem>& heap) {
    if (node == NULL) {
      return;     // indicates that we're done here
    }
    // Compute distance between target and current node
    double dist = distance(_items[node->index], target);

    // If current node within radius tau
    if (dist < _tau) {
      if (heap.size() == k) {
        heap.pop();  // remove furthest node from result list (if we already have k results)
      }
      heap.push(HeapItem(node->index, dist));  // add current node to result list
      if (heap.size() == k) {
        _tau = heap.top().dist;  // update value of tau (farthest point in result list)
      }
    }

    // Return if we arrived at a leaf
    if (node->left == NULL && node->right == NULL) {
      return;
    }

    // If the target lies within the radius of ball
    if (dist < node->threshold) {
      if (dist - _tau <= node->threshold) {
        // if there can still be neighbors inside the ball, recursively search left child first
        search(node->left, target, k, heap);
      }

      if (dist + _tau >= node->threshold) {
        // if there can still be neighbors outside the ball, recursively search right child
        search(node->right, target, k, heap);
      }
    // If the target lies outsize the radius of the ball
    } else {
      if (dist + _tau >= node->threshold) {
        // if there can still be neighbors outside the ball, recursively search right child first
        search(node->right, target, k, heap);
      }

      if (dist - _tau <= node->threshold) {
        // if there can still be neighbors inside the ball, recursively search left child
        search(node->left, target, k, heap);
      }
    }
}
};
}  // namespace datadriven
}  // namespace sgpp
