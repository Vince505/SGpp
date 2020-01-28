// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/postProcessing/PostProcessingClassification.hpp>

#include <iostream>


namespace sgpp {
namespace datadriven {

void PostProcessingClassification::postProcessing(DataSource &datasource,
  ModelFittingBase &model, Visualizer &visualizer) {
  std::cout << "For Classification Models no post processing is required" << std::endl;
}

} // namespace datadriven
} // namespace sgpp