// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org


#include <sgpp/datadriven/datamining/modules/scoring/FowlkesMallows.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

#include <map>
#include <cmath>
#include <iostream>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
namespace sgpp {
namespace datadriven {

Metric *FowlkesMallows::clone() const { return new FowlkesMallows(*this); }

double FowlkesMallows::measurePostProcessing(const DataVector &predictedValues,
  const DataVector &trueValues, const ModelFittingBase &model, Dataset &testDataset) const {
  DataMatrix countMatrix(0,0,0.0);

  size_t nTrueValues = 0;
  size_t nPredictedValues = 0;

  std::map<int, size_t> trueLabelsMap;
  std::map<int, size_t> predLabelsMap;


  for (size_t index = 0; index < predictedValues.size() ; index++) {
    auto trueValue = trueValues.get(index);
    auto predictedValue = predictedValues.get(index);

    if (trueLabelsMap.find(trueValue) == trueLabelsMap.end()) {
      trueLabelsMap[trueValue] = nTrueValues++;
      countMatrix.resizeRows(nTrueValues);
    }

    if (predLabelsMap.find(predictedValue) == predLabelsMap.end()) {
      predLabelsMap[predictedValue] = nPredictedValues++;
      countMatrix.resizeRowsCols(nTrueValues, nPredictedValues);
    }

    countMatrix.set(trueLabelsMap[trueValue], predLabelsMap[predictedValue],
      countMatrix.get(trueLabelsMap[trueValue], predLabelsMap[predictedValue])+1);
  }

  DataMatrix squareMatrix (countMatrix);

  squareMatrix.sqr();

  double tk = squareMatrix.sum() - trueValues.size();

  DataVector rowVector(countMatrix.getNcols());

  double pk = 0;
  for (size_t row = 0 ; row < countMatrix.getNrows() ; row++) {
    countMatrix.getRow(row, rowVector);
    pk += pow(rowVector.sum(), 2);
  }

  pk = pk - trueValues.size();

  countMatrix.transpose();
  DataVector colVector(countMatrix.getNcols());
  double qk = 0;

  for (size_t col = 0 ; col < countMatrix.getNrows(); col++) {
    countMatrix.getRow(col, colVector);
    qk += pow(colVector.sum(), 2);
  }

  qk = qk - trueValues.size();

  double score = tk/sqrt(pk*qk);

  return score;
}

} //  namespace datadriven
} //  namespace sgpp