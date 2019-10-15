// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/combigrid/LevelIndexTypes.hpp>
#include <sgpp/combigrid/basis/HeterogeneousBasis.hpp>

namespace sgpp {
namespace combigrid {

/**
 * Full grid essentially represented by its level and a HeterogeneousBasis.
 */
class FullGrid {
 public:
  /**
   * Default constructor, corresponds to the zero-dimensional case.
   */
  FullGrid() : level(), hasBoundary_(true), basis() {
  }

  /**
   * Constructor.
   *
   * @param level         level of the full grid
   * @param basis         type of basis functions for evaluating on the full grid
   * @param hasBoundary   whether the full grid has points on the boundary
   */
  FullGrid(const LevelVector& level, const HeterogeneousBasis& basis, bool hasBoundary = true) :
      level(level), hasBoundary_(hasBoundary), basis(basis) {
  }

  /**
   * @param other   other FullGrid instance
   * @return whether both instances are equal
   */
  bool operator==(const FullGrid& other) const {
    return (level == other.level) && (hasBoundary_ == other.hasBoundary_) && (basis == other.basis);
  }

  /**
   * @param other   other FullGrid instance
   * @return whether both instances are equal
   */
  bool operator!=(const FullGrid& other) const {
    return !(*this == other);
  }

  /**
   * @return level of the full grid
   */
  const LevelVector& getLevel() const {
    return level;
  }

  /**
   * @return level of the full grid
   */
  LevelVector& getLevel() {
    return level;
  }

  /**
   * @param d   dimension
   * @return level of the \f$d\f$-th dimension
   */
  size_t getLevel(size_t d) const {
    return level[d];
  }

  /**
   * @return level of the full grid
   */
  void setLevel(const LevelVector& level) {
    this->level = level;
  }

  /**
   * Minimum 1D index of grid points.
   *
   * @param d   dimension
   * @return minimum index in the \f$d\f$-th dimension
   */
  index_t getMinIndex(size_t d) const {
    return (hasBoundary_ ? 0 : 1);
  }

  /**
   * Minimum index of grid points.
   *
   * @param[out] index  minimum index as multi-index
   */
  void getMinIndex(IndexVector& index) const {
    index.resize(level.size());

    for (size_t d = 0; d < level.size(); d++) {
      index[d] = getMinIndex(d);
    }
  }

  /**
   * Maximum 1D index of grid points.
   *
   * @param d   dimension
   * @return maximum index in the \f$d\f$-th dimension
   */
  index_t getMaxIndex(size_t d) const {
    return (static_cast<index_t>(1) << level[d]) - getMinIndex(d);
  }

  /**
   * Maximum index of grid points.
   *
   * @param[out] index  maximum index as multi-index
   */
  void getMaxIndex(IndexVector& index) const {
    index.resize(level.size());

    for (size_t d = 0; d < level.size(); d++) {
      index[d] = getMaxIndex(d);
    }
  }

  /**
   * Number of index vectors (grid points) in 1D.
   *
   * @param d   dimension
   * @return number of index vectors in the \f$d\f$-th dimension
   */
  index_t getNumberOfIndexVectors(size_t d) const {
    return getMaxIndex(d) - getMinIndex(d) + 1;
  }

  /**
   * Number of index vectors (grid points) in 1D for all dimensions.
   *
   * @param[out] number   vector, \f$d\f$-th entry is number of index vectors in the
   *                      \f$d\f$-th dimension
   */
  void getNumberOfIndexVectors(IndexVector& number) const {
    number.resize(level.size());

    for (size_t d = 0; d < level.size(); d++) {
      number[d] = getNumberOfIndexVectors(d);
    }
  }

  /**
   * Total number of index vectors (grid points).
   *
   * @return number of index vectors
   */
  index_t getNumberOfIndexVectors() const {
    index_t result = 1;

    for (size_t d = 0; d < level.size(); d++) {
      result *= getNumberOfIndexVectors(d);
    }

    return result;
  }

  /**
   * @return dimensionality
   */
  size_t getDimension() const {
    return level.size();
  }

  /**
   * @return whether the full grid has points on the boundary
   */
  bool hasBoundary() const {
    return hasBoundary_;
  }

  /**
   * @param hasBoundary   whether the full grid has points on the boundary
   */
  void setHasBoundary(bool hasBoundary) {
    this->hasBoundary_ = hasBoundary;
  }

  /**
   * @return type of basis functions for evaluating on the full grid
   */
  const HeterogeneousBasis& getBasis() const {
    return basis;
  }

  /**
   * @param basis   type of basis functions for evaluating on the full grid
   */
  void setBasis(const HeterogeneousBasis& basis) {
    this->basis = basis;
  }

 protected:
  /// level of the full grid
  LevelVector level;
  /// whether the full grid has points on the boundary
  bool hasBoundary_;
  /// type of basis functions for evaluating on the full grid
  HeterogeneousBasis basis;
};

}  // namespace combigrid
}  // namespace sgpp
