// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

# pragma once

#include <sgpp/datadriven/datamining/modules/postProcessing/PostProcessingBase.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/visualization/Visualizer.hpp>

namespace sgpp {
namespace datadriven {
class PostProcessingDensityEstimation : public PostProcessingBase {
public:
  void postProcessing(DataSource &datasource,
    ModelFittingBase &model, Visualizer &visualizer) override;
};

} // namespace datadriven
} // namespace sgpp

