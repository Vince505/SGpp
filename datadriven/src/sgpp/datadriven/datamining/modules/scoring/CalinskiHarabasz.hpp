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
 * Metric that delivers the Calinski-Harabasz score
 * to determine the quality of a clustering if the data doesn't contain labels
 * to compare against
 */
class CalinskiHarabasz: public ClusteringMetric {
 public:
  Metric *clone() const override;

  /**
   * Quantify the Calinski-Harabasz score
   * of a clustering after postprocessing
   * @params model The fitted model
   * @params datasource The source pointing the data
   * @return the Calinski-Harabasz score of a clustering
   */
  double measurePostProcessing(ModelFittingBase &model, DataSource &datasource) const override;
};

}  // namespace datadriven
}  // namespace sgpp


