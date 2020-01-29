// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/postProcessing/PostProcessingDensityEstimation.hpp>

#include <iostream>


namespace sgpp {
namespace datadriven {

void PostProcessingDensityEstimation::postProcessing(DataSource &datasource,
  ModelFittingBase &model, Visualizer &visualizer, size_t fold) {
  visualizer.runPostProcessingVisualization(model, datasource, fold);
}
}  // namespace datadriven
}  // namespace sgpp
