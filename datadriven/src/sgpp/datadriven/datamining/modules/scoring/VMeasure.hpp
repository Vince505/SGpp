// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/datamining/modules/scoring/Metric.hpp>

namespace sgpp {
namespace datadriven {
/**
 * Metric that delivers a V measure score metric that determine the quality of a clustering if
 * the data contains labels
 * to compare against
 */
class VMeasure: public Metric {
public:
  Metric *clone() const override;

  /**
    * Quantify the V Measure score
    * of a clustering
    *
    * @param predictedValues clustering given  by the model for testing data
    * @param trueValues trueLabels of the clustering
    * @param model reference to the model
    * @param testDataset dataset with test data
    * @return the v Measure score of a clustering
    */
  double measure(const DataVector &predictedValues, const DataVector &trueValues,
                 const ModelFittingBase &model, Dataset &testDataset) const override;


  /**
    * Quantify the V Measure score
      * of a clustering
    *
    * @param predictedValues clustering given  by the model for testing data
    * @param trueValues trueLabels of the clustering
    * @param model reference to the model
    * @param testDataset dataset with test data
    * @return the v Measure score of a clustering
    */
  double measureLowerIsBetter(const DataVector &predictedValues, const DataVector &trueValues,
                              const ModelFittingBase &model, Dataset &testDataset) const override;

};

} //  namespace datadriven
} //  namespace sgpp
