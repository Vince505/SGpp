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

double VMeasure::measurePostProcessing(ModelFittingBase &model, DataSource &datasource) const {
  // Obtaining the count of points belonging to a pair of classes

  DataMatrix countMatrix(0, 0, 0.0);

  size_t nTrueValues = 0;
  size_t nPredictedValues = 0;

  std::map<int, size_t> trueLabelsMap;
  std::map<int, size_t> predLabelsMap;

  Dataset* testDataset = datasource.getAllSamples();
  DataVector trueValues = testDataset->getTargets();

  DataVector predictedValues(testDataset->getNumberInstances());
  model.evaluate(testDataset->getData(), predictedValues);

  for (size_t index = 0; index < predictedValues.size() ; index++) {
    auto trueValue = static_cast<int>(trueValues.get(index));
    auto predictedValue = static_cast<int>(predictedValues.get(index));

    if (trueLabelsMap.find(trueValue) == trueLabelsMap.end()) {
      trueLabelsMap[trueValue] = nTrueValues++;
    }

    if (predLabelsMap.find(predictedValue) == predLabelsMap.end()) {
      predLabelsMap[predictedValue] = nPredictedValues++;
    }
  }
  countMatrix.resizeRowsCols(nTrueValues, nPredictedValues);

  for (size_t index = 0; index < predictedValues.size() ; index++) {
    auto trueValue = static_cast<int>(trueValues.get(index));
    auto predictedValue = static_cast<int>(predictedValues.get(index));
    auto previousValue = countMatrix.get(trueLabelsMap[trueValue], predLabelsMap[predictedValue]);
    countMatrix.set(trueLabelsMap[trueValue], predLabelsMap[predictedValue], (previousValue+1));
  }
  // CALCULATING HOMOGENEITY
  //  Entropy of classes

  double classConditionalEntropy = 0.0;
  double classEntropy = 0.0;

  DataVector predictedVector(countMatrix.getNrows());

  DataVector trueVector(countMatrix.getNcols());

  double n = static_cast<double>(trueValues.size());

  for (auto &trueLabel : trueLabelsMap) {
    // Calculating class conditional entropy
    for (auto &predLabel : predLabelsMap) {
      double nck = countMatrix.get(trueLabel.second, predLabel.second);
      countMatrix.getColumn(predLabel.second, predictedVector);
      double nk = predictedVector.sum();

      if (nck != 0) {
        classConditionalEntropy = classConditionalEntropy +
                                  (nck / n) * log(nck / nk);
      }
    }

    countMatrix.getRow(trueLabel.second, trueVector);

    double nc = trueVector.sum();

    // Calculating class entropy
    if (nc > 0) {
      classEntropy = classEntropy + (nc / n) * log(nc / n);
    }
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


  for (auto &predLabel : predLabelsMap) {
    // Calculating cluster conditional entropy
      for (auto &trueLabel : trueLabelsMap) {
      double nck = countMatrix.get(trueLabel.second, predLabel.second);
      countMatrix.getRow(trueLabel.second, trueVector);

      double nc = trueVector.sum();
      if (nck != 0) {
        clusterConditionalEntropy = clusterConditionalEntropy +
                                    (nck / n) * log(nck / nc);
      }
    }

    countMatrix.getColumn(predLabel.second, predictedVector);
    double nk = predictedVector.sum();

    // Calculating cluster entropy
    if (nk > 0) {
      clusterEntropy = clusterEntropy + (nk / n) * log(nk / n);
    }
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

  std::cout << "V Measure Score is " << vMeasure << std::endl;
  return vMeasure;
}
}  // namespace datadriven
}  // namespace sgpp
