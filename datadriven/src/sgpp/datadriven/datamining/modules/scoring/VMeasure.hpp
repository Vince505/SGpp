// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/datamining/modules/scoring/ClusteringMetric.hpp>

namespace sgpp {
namespace datadriven {
/**
 * Metric that delivers a V measure score metric that determine the quality of a clustering if
 * the data contains labels
 * to compare against
 */
class VMeasure: public ClusteringMetric {
 public:
  Metric *clone() const override;
  /**
   * Quantify the V Measure score
   * of a clustering
   * @param model The fitted model
   * @param datasource The source pointing the data
   * @return he v Measure score of a clustering
   */
  double measurePostProcessing(ModelFittingBase &model, DataSource &datasource) const override;
};
}  //  namespace datadriven
}  //  namespace sgpp
