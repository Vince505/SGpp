// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/visualization/VisualizerDummy.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>


using sgpp::datadriven::ModelFittingBase;

namespace sgpp {
namespace datadriven {

void VisualizerDummy::runVisualization(ModelFittingBase &model, DataSource &dataSource, size_t epoch,
  size_t fold, size_t batch) {
  return;
}

void VisualizerDummy::runPostProcessingVisualization(ModelFittingBase &model,
  DataSource &dataSource) {
  return;
}

void VisualizerDummy::storeScatterPlotJson(DataMatrix &matrix, ModelFittingBase &model,
                          std::string currentDirectory) {
  return;
}

void VisualizerDummy::storeHeatmapJson(DataMatrix &matrix, ModelFittingBase &model,
  std::vector<size_t> indexes, size_t &varDim1, size_t &varDim2, std::string filepath) {
  return;
}

void VisualizerDummy::storeHeatmapJson(DataMatrix &matrix,
  ModelFittingBase &model, std::string filepath) {
  return;
}

void VisualizerDummy::storeCutJson(DataMatrix &matrix, std::vector<size_t> indexes,
                                            size_t &varDim, std::string filepath) {
  return;
}

void VisualizerDummy::storeCutJson(DataMatrix &matrix, std::string filepath) {
  return;
}

}  // namespace datadriven
}  // namespace sgpp
