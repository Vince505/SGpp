// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org


#include <sgpp/datadriven/datamining/modules/scoring/ClusteringMetric.hpp>

namespace sgpp {
namespace datadriven {

double ClusteringMetric::measure(const DataVector &predictedValues, const DataVector &trueValues,
                             const ModelFittingBase &model, Dataset &testDataset) const {
  return mseMetric.measure(predictedValues, trueValues, model, testDataset);
}

double ClusteringMetric::measureLowerIsBetter(const DataVector &predictedValues,
  const DataVector &trueValues, const ModelFittingBase &model, Dataset &testDataset) const {
  return measure(predictedValues, trueValues, model, testDataset);
}

}  // namespace datadriven
}  // namespace sgpp

