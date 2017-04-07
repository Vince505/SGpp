/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DBMatOfflineIChol.cpp
 *
 *  Created on: Feb 27, 2017
 *      Author: Michael Lettrich
 */

#include <sgpp/datadriven/algorithm/DBMatOfflineIChol.hpp>

#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/datadriven/algorithm/IChol.hpp>
#include <sgpp/datadriven/algorithm/SparseDataMatrix.hpp>

namespace sgpp {
namespace datadriven {

using sgpp::base::algorithm_exception;

DBMatOfflineIChol::DBMatOfflineIChol(const DBMatDensityConfiguration& oc) : DBMatOfflineChol(oc) {}

void DBMatOfflineIChol::decomposeMatrix() {
  if (isConstructed) {
    if (isDecomposed) {
      return;
    } else {
      SparseDataMatrix sparseLHS(lhsMatrix);
      DataVector norms{sparseLHS.getNrows()};
      IChol::normToUnitDiagonal(sparseLHS, norms);
      IChol::decompose(sparseLHS, 3);
      IChol::reaplyDiagonal(sparseLHS, norms);
      SparseDataMatrix::toDataMatrix(sparseLHS, lhsMatrix);
    }

    isDecomposed = true;

  } else {
    throw algorithm_exception("Matrix has to be constructed before it can be decomposed");
  }
}

} /* namespace datadriven */
} /* namespace sgpp */
