// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearClenshawCurtisBasis.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

class OperationMultipleEvalLinearClenshawCurtisNaive : public OperationMultipleEval {
 public:
  OperationMultipleEvalLinearClenshawCurtisNaive(Grid& grid, DataMatrix& dataset)
      : OperationMultipleEval(grid, dataset), storage(grid.getStorage()) {}

  ~OperationMultipleEvalLinearClenshawCurtisNaive() override {}

  void mult(DataVector& alpha, DataVector& result) override;
  void multTranspose(DataVector& source, DataVector& result) override;

  double getDuration() override;

 protected:
  /// storage of the sparse grid
  GridStorage& storage;
  /// 1D B-spline basis
  SLinearClenshawCurtisBase base;
  /// untransformed evaluation point (temporary vector)
  DataMatrix pointsInUnitCube;
};

}  // namespace base
}  // namespace sgpp
