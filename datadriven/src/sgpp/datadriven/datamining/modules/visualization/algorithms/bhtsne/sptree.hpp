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

/**
 * Code originally taken from https://lvdmaaten.github.io/tsne/
 * It has been modified in order to be adapted to the SG++ datamining
 * pipeline structure and has been parallelized
 */

#include <iostream>

namespace sgpp {
namespace datadriven {
/**
 * Reference class to represent a Cell in the quadtree
 */
class Cell {
  /**
   * Dimension of the cell
   */
  size_t dimension;
  /**
   * Stores the coordinates where the lowest right corner of the cell is found
   */
  double* corner;
  /**
   * Stores the width of each cell per dimension
   */
  double* width;

 public:
  /**
   * Constructor
   * @param inp_dimension Dimension of the cell
   */
  explicit Cell(size_t inp_dimension);
  /**
   * Constructor
   * @param inp_dimension Dimension of the cell
   * @param inp_corner Coordinates where the lowest right corner of the cell is found
   * @param inp_width Array containing the width per dimension
   */
  Cell(size_t inp_dimension, double* inp_corner, double* inp_width);
  /**
   * Destructor
   */
  ~Cell();

  /**
   * Returns the dth component of the corner
   * @param d Component of the corner to return
   * @return Dth component of the corner
   */
  double getCorner(size_t d);
  /**
   * Returns the dth component of the width
   * @param d Component of the width to return
   * @return Dth component of the corner
   */
  double getWidth(size_t d);
  /**
   * Sets the dth component of the corner
   * @param d dth component to set
   * @param val Value to set
   */
  void setCorner(size_t d, double val);
  /**
   * Sets the dth component of the width
   * @param d dth component to set
   * @param val Value to set
   */
  void setWidth(size_t d, double val);
  /**
   * Asks if a cell contains a given point
   * @param point Array containing the coordinates of the point to search
   * @return True if the cell contains the point, False otherwise
   */
  bool containsPoint(double point[]);
};

/**
 * Reference class to represent a quad or octree to approximate the t-SNE gradient
 */
class SPTree {
  // Fixed constants
  static const size_t QT_NODE_CAPACITY = 1;

  // A buffer we use when doing force computations
  double* buff;

  // Pointer to the parent
  SPTree* parent;
  // Dimension of the three
  size_t dimension;
  // Flag to determine if a node is a leaf
  bool is_leaf;
  // Size of the whole tree
  size_t size;
  // Size of the tree starting this node
  size_t cum_size;

  // Axis-aligned bounding box stored as a center with half-dimensions to
  // represent the boundaries of this quad tree
  Cell* boundary;

  // Indices in this space-partitioning tree node, corresponding center-of-mass,
  // and list of all children
  // Data inside this node
  double* data;
  // Centroid of this node.
  double* center_of_mass;
  size_t index[QT_NODE_CAPACITY];

  // Array of Children
  SPTree** children;
  // Number of children
  size_t no_children;

 public:
  /**
   * Default Constructor. Builds the a tree node and a tree
   * @param D Dimensionality of the data
   * @param inp_data Array containing all coordinates of all data points
   * @param N Number of data points
   */
  SPTree(size_t D, double* inp_data, size_t N);
  /**
   * Creates a tree node with a particular size.
   * @param D Dimensionality of the data
   * @param inp_data Array containing all coordinates of all data points
   * @param inp_corner Coordinates where the lowest right corner of the cell is found
   * @param inp_width Array containing the width per dimension of the cell
   */
  SPTree(size_t D, double* inp_data, double* inp_corner, double* inp_width);
  /**
   * Creates a tree node with a particular size.
   * @param D Dimensionality of the data
   * @param inp_data Array containing all coordinates of all data points
   * @param N Number of data points
   * @param inp_corner Coordinates where the lowest right corner of the cell is found
   * @param inp_width Array containing the width per dimension of the cell
   */
  SPTree(size_t D, double* inp_data, size_t N, double* inp_corner,
    double* inp_width);
  /**
   * Creates a tree node with a particular size and with a particular parent.
   * @param inp_parent Reference to the parent of this node
   * @param D Dimensionality of the data
   * @param inp_data Array containing all coordinates of all data points
   * @param N Number of data points
   * @param inp_corner Coordinates where the lowest right corner of the cell is found
   * @param inp_width Array containing the width per dimension of the cell
   */
  SPTree(SPTree* inp_parent, size_t D, double* inp_data, size_t N,
    double* inp_corner, double* inp_width);
  /**
   * Creates a tree node with a particular size.
   * @param inp_parent Reference to the parent of this node
   * @param D Dimensionality of the data
   * @param inp_data Array containing all coordinates of all data points
   * @param inp_corner Coordinates where the lowest right corner of the cell is found
   * @param inp_width Array containing the width per dimension of the cell
   */
  SPTree(SPTree* inp_parent, size_t D, double* inp_data,
    double* inp_corner, double* inp_width);
  /**
   * Default destructor
   */
  ~SPTree();
  /**
   * Sets the data of this node
   * @param inp_data Array containing all coordinates of the data points
   */
  void setData(double* inp_data);
  /**
   * Returns reference to the parent of this node
   * @return Reference to the parent of this node
   */
  SPTree* getParent();
  void construct(Cell boundary);
  /**
   * Inserts a new point in the tree
   * @param new_index Index referencing the point
   * @return True if point was inserted succesfully, false if otherwise
   */
  bool insert(size_t new_index);
  /**
   * Subdivides a tree when no more points fit in a node.
   */
  void subdivide();
  /**
   * Check if the tree starting from this node is correct
   * @return True if tree is correct, false otherwise
   */
  bool isCorrect();
  /**
   * Rebuilds the tree
   */
  void rebuildTree();
  /**
   * Gets all of the indices of the inserted points
   * @param indices Array to store the indices list
   */
  void getAllIndices(size_t* indices);
  /**
   * Gets the maximum depth of the tree
   * @return Maximum depth of he tree
   */
  size_t getDepth();
  /**
   * Computes the negative forces part of the gradient for t-SNE
   * @param point_index Index of the point we are processing
   * @param theta Theta parameter for gradient approximation
   * @param neg_f Array where the positive forces part of the gradient will be stored
   * @param sum_Q Stores the cumulative sum of all negative forces
   */
  void computeNonEdgeForces(size_t point_index, double theta, double neg_f[], double* sum_Q);
  /**
   * Computes the positive forces part of the gradient for t-SNE
   * @param row_P Array with pointers to the array inp_col_P indicating,
   * where the probability values for a given row lie in the inp_val_P
   * @param col_P  Array with pointers to the array val_P
   * @param val_P Array with the input probability values
   * @param N Total number of points
   * @param pos_f Array where the positive forces part of the gradient will be stored
   */
  void computeEdgeForces(size_t* row_P, size_t* col_P,
    double* val_P, size_t N, double* pos_f);

 private:
  /**
   * Initialize a Node and builds the tree
   * @param inp_parent Reference to the parent of this node
   * @param D Dimensionality of the data
   * @param inp_data Array containing all coordinates of all data points
   * @param inp_corner Coordinates where the lowest right corner of the cell is found
   * @param inp_width Array containing the width per dimension of the cell
   */
  void init(SPTree* inp_parent, size_t D, double* inp_data, double*
    inp_corner, double* inp_width);
  /**
   * Fills the tree
   * @param N Total number of points
   */
  void fill(size_t N);
  /**
   * Build a list of all indexes in the SPTree
   * @param indices Array to store the indices
   * @param loc The start location in the indices array for this node
   * @return The end location in the indicies array for this node
   */
  size_t getAllIndices(size_t* indices, size_t loc);
};

}  // namespace datadriven
}  // namespace sgpp
