// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/visualization/Visualizer.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSource.hpp>

namespace sgpp {
namespace datadriven {
/**
 * Class used to run post processing procedures after the training within the Pipeline
 */
class PostProcessingBase {
 public:
  /**
   * Default constructor
   */
  PostProcessingBase() = default;

  /**
   * Default destructor
   */
  ~PostProcessingBase() = default;

  /**
   * Method which executes post processing processes of the fitter and the visualizer
   * @param datasource Datasource pointing to the data
   * @param model Fitter containing the model
   * @param visualizer Visualizer to generate the viualization part
   * @param fold Fold being processed in case of cross validation
   */
  virtual void postProcessing(DataSource &datasource,
    ModelFittingBase &model, Visualizer &visualizer, size_t fold = 0) = 0;
};
}  // namespace datadriven
}  // namespace sgpp
