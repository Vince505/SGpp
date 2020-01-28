// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/datamining/modules/scoring/ClusteringMetric.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

namespace sgpp {
namespace datadriven {
/**
 * Metric that delivers David-Bouldin Index
 * that determine the quality of a clustering if the data doesn't contain labels
 * to compare against
 */
class DavidBouldin: public ClusteringMetric {
public:
  Metric *clone() const override;
  /**
   * Quantify the David-Bouldin score
   * of a clustering after postprocessing
   *
   * @param predictedValues ignored
   * @param trueValues ignored
   * @param model reference to the model
   * @param testDataset ignored
   * @return the David-Bouldin score of a clustering
   */
  double measurePostProcessing(const DataVector &predictedValues, const DataVector &trueValues,
                               const ModelFittingBase &model, Dataset &testDataset) const override;


};

} //  namespace datadriven
} //  namespace sgpp


