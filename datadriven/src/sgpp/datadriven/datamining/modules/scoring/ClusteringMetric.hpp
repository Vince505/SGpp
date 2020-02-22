// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/datamining/modules/scoring/Metric.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/MSE.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

namespace sgpp {
namespace datadriven {
class ClusteringMetric: public Metric {
 public:
  /**
   * Quantify the score for the currently trained density estimation model using the mean error
   * squared metric
   *
   * @param predictedValues values calculated by the model for testing data
   * @param trueValues actual values as taken from the dataset.
   * @param model reference to the model
   * @param testDataset dataset with test data
   * @return mean squared error (MSE) - strictly positive such that smaller values are better.
   */
  double measure(const DataVector &predictedValues, const DataVector &trueValues,
                 const ModelFittingBase &model, Dataset &testDataset) const override;

  /**
   * Quantify the score for the currently trained density estimation model using the mean error
   * squared metric
   *
   * @param predictedValues values calculated by the model for testing data
   * @param trueValues actual values as taken from the dataset.
   * @param model reference to the model
   * @param testDataset dataset with test data
   * @return mean squared error (MSE) - strictly positive such that smaller values are better.
   */
  double measureLowerIsBetter(const DataVector &predictedValues, const DataVector &trueValues,
                              const ModelFittingBase &model, Dataset &testDataset) const override;


  /**
   * Quantify the quality of a cluster using a specific metric
   * @param model The fitted model
   * @param datasource The source pointing the data
   * @return Cluster quality score Value
   */
  virtual double measurePostProcessing(ModelFittingBase &model, DataSource &datasource) const = 0;

 private:
  MSE mseMetric;
};

}  // namespace datadriven
}  // namespace sgpp

