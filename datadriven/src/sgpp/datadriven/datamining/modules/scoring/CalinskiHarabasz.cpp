// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/scoring/CalinskiHarabasz.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingClustering.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

#include <iostream>
#include <map>
#include <vector>
#include <cmath>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;

namespace sgpp {
namespace datadriven {

Metric *CalinskiHarabasz::clone() const { return new CalinskiHarabasz(*this);}

double CalinskiHarabasz::measurePostProcessing(ModelFittingBase &model, DataSource &datasource)
const {
  auto clusteringModel = dynamic_cast<ModelFittingClustering*>(&model);

  DataMatrix samples = clusteringModel->getPoints();

  auto hierarchy = clusteringModel->getHierarchyTree();

  DataVector allLabels;

  (*hierarchy)->evaluateClustering(allLabels);

  // Stores the labels
  std::vector<int> clusterLabels;
  // stores the number of points per cluster labels
  std::map<int, size_t> pointsPerLabel;
  // Stores the cluster centroids of a given label
  std::map<int, DataVector> clusterCentroids;

  DataVector globalCentroid(samples.getNcols());
  // Getting the number of labels and the centroids
  for (size_t index = 0; index < allLabels.size() ; index++) {
    auto value = static_cast<int>(allLabels.get(index));
    if (std::find(clusterLabels.begin(), clusterLabels.end(), value) == clusterLabels.end()) {
        clusterLabels.push_back(value);
        DataVector centroid(samples.getNcols(), 0);
        clusterCentroids[value] = centroid;
      }
      DataVector row(samples.getNcols());
      samples.getRow(index, row);
      clusterCentroids[value].add(row);
      pointsPerLabel[value]++;
      globalCentroid.add(row);
  }
  if (clusterLabels.size() == 1) {
    std::cout << "To use the Calinski-Harabasz score, "
                 "at least 2 clusters are required" << std::endl;
    return 0.0;
  }
  // Getting the global centroid
  globalCentroid.mult(1/ static_cast<double>(allLabels.size()));

  // Between group dispersion
  double betweenClusterDispersion = 0;
  for (auto &label : clusterLabels) {
    // Getting the centroid of the cluster
    clusterCentroids[label].mult(1/ static_cast<double>(pointsPerLabel[label]));
    DataVector temp(clusterCentroids[label]);
    temp.sub(globalCentroid);

    betweenClusterDispersion+= (static_cast<double>(pointsPerLabel[label])*pow(temp.l2Norm(), 2));
  }

  // Within group dispersion
  double withinClusterDispersion = 0;
  for (size_t index = 0; index < allLabels.size() ; index++) {
    auto value = static_cast<int>(allLabels.get(index));
    DataVector row(samples.getNcols());
    samples.getRow(index, row);
    row.sub(clusterCentroids[value]);
    withinClusterDispersion += pow(row.l2Norm(), 2);
  }
  // Calculating final score

  double score = (static_cast<double>(betweenClusterDispersion)/
    static_cast<double>(withinClusterDispersion))*
    (static_cast<double>(samples.getNrows()-clusterLabels.size())/
    static_cast<double>(clusterLabels.size()-1));

  std::cout << "Calinksi Harabasz Score is " << score << std::endl;
  return score;
}

}  // namespace datadriven
}  // namespace sgpp
