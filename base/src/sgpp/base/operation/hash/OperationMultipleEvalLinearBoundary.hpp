// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONMULTIPLEEVALLINEARBOUNDARY_HPP
#define OPERATIONMULTIPLEEVALLINEARBOUNDARY_HPP

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

/**
 * This class implements OperationMultipleEval for a grids with linear basis ansatzfunctions
 * with boundaries
 *
 */
class OperationMultipleEvalLinearBoundary : public OperationMultipleEval {
 public:
  /**
   * Constructor
   *
   * @param grid grid
   * @param dataset the dataset the should be evaluated
   */
  OperationMultipleEvalLinearBoundary(Grid& grid, DataMatrix& dataset)
      : OperationMultipleEval(grid, dataset), storage(grid.getStorage()) {}

  /**
   * Destructor
   */
  ~OperationMultipleEvalLinearBoundary() override {}

  void mult(DataVector& alpha, DataVector& result) override;
  void multTranspose(DataVector& source, DataVector& result) override;

  double getDuration() override;

 protected:
  /// Pointer to GridStorage object
  GridStorage& storage;
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATIONMULTIPLEEVALLINEARBOUNDARY_HPP */
