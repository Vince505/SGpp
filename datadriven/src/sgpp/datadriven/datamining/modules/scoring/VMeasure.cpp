// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org


#include <sgpp/datadriven/datamining/modules/scoring/VMeasure.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

#include <map>
#include <cmath>
#include <iostream>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
namespace sgpp {
namespace datadriven {

Metric *VMeasure::clone() const { return new VMeasure(*this); }

double VMeasure::measurePostProcessing(const DataVector &predictedValues,
  const DataVector &trueValues, const ModelFittingBase &model, Dataset &testDataset) const {
  // Obtaining the count of points belonging to a pair of classes

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
  // CALCULATING HOMOGENEITY
  //  Entropy of classes

  double classConditionalEntropy = 0.0;
  double classEntropy = 0.0;

  DataMatrix transposeCountMatrix(countMatrix);

  transposeCountMatrix.transpose();
  DataVector predictedVector (transposeCountMatrix.getNcols());

  DataVector trueVector(countMatrix.getNcols());

  double n = trueValues.size();

  for (auto &trueLabel: trueLabelsMap) {

    // Calculating class conditional entropy
    for (auto &predLabel: predLabelsMap) {
      double nck = countMatrix.get(trueLabel.second, predLabel.second);
      transposeCountMatrix.getRow(predLabel.second, predictedVector);
      double nk = predictedVector.sum();

      if (nck != 0) {
        classConditionalEntropy = classConditionalEntropy +
                                  (nck / n) * log(nck / nk);
      }
    }

    countMatrix.getRow(trueLabel.second, trueVector);

    double nc = trueVector.sum();

    // Calculating class entropy
    classEntropy = classEntropy + (nc/n) * log(nc/n);
  }

  classConditionalEntropy = -1*classConditionalEntropy;

  classEntropy = -1*classEntropy;

  double h;
  if (classEntropy != 0) {
    h = 1 - (classConditionalEntropy / classEntropy);
  } else {
    h = 1;
  }


  // CALCULATING COMPLETENESS
  //  Entropy of clusters

  double clusterConditionalEntropy = 0.0;
  double clusterEntropy = 0.0;


  for (auto &predLabel: predLabelsMap) {
    // Calculating cluster conditional entropy
      for (auto &trueLabel: trueLabelsMap) {
      double nck = countMatrix.get(trueLabel.second, predLabel.second);
      countMatrix.getRow(trueLabel.second, trueVector);

      double nc = trueVector.sum();
      if (nck != 0) {
        clusterConditionalEntropy = clusterConditionalEntropy +
                                    (nck / n) * log(nck / nc);
      }
    }

    transposeCountMatrix.getRow(predLabel.second, predictedVector);
    double nk = predictedVector.sum();

    // Calculating cluster entropy
    clusterEntropy = clusterEntropy + (nk/n) * log2(nk/n);
  }

  clusterConditionalEntropy = -1*clusterConditionalEntropy;

  clusterEntropy = -1*clusterEntropy;

  // Calculating Completeness
  double c;
  if (clusterEntropy != 0) {
    c = 1 - (clusterConditionalEntropy / clusterEntropy);
  } else {
    c = 1;
  }

  // Calculating VMeasure
  double vMeasure = 2 *((h*c)/(h+c));

  std::cout << "Homogeneity: " << h << std::endl;
  std::cout << "Completeness: " << c << std::endl;
  return vMeasure;
}

} //  namespace datadriven
} //  namespace sgpp