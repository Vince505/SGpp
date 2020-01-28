// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/scoring/DavidBouldin.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingClustering.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

#include <iostream>
#include <map>
#include <cmath>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;

namespace sgpp {
namespace datadriven {

Metric *DavidBouldin::clone() const { return new DavidBouldin(*this);}

double DavidBouldin::measurePostProcessing(const DataVector &predictedValues,
  const DataVector &trueValues,
  const ModelFittingBase &model, Dataset &testDataset) const {
  auto clusteringModel = dynamic_cast<const ModelFittingClustering*>(&model);

  DataMatrix samples = clusteringModel->getPoints();

  std::vector<size_t> clusterLabels;
  std::map<size_t, size_t> pointsPerLabel;
  std::map<size_t, DataVector> clusterCentroids;
  std::map<size_t, double> averageDistances;

  size_t numberNoise = 0;

  // Getting the number of labels and the centroids
  for (size_t index = 0; index < predictedValues.size() ; index++) {
    auto value = predictedValues.get(index);
    if (value != -1) {  // Skipping all noisy data
      if (std::find(clusterLabels.begin(), clusterLabels.end(), value) == clusterLabels.end()) {
        clusterLabels.push_back(value);
        DataVector centroid(samples.getNcols(), 0);
        clusterCentroids[value] = centroid;
      }
      DataVector row(samples.getNcols());
      samples.getRow(index, row);
      clusterCentroids[value].add(row);
      pointsPerLabel[value]++;
    } else {
      numberNoise++;
    }
  }
  for (auto &label: clusterLabels) {
    clusterCentroids[label].mult(1/ static_cast<double>(pointsPerLabel[label]));
  }

  // Calculating the average distances to the cluster centroids
  for (size_t index = 0; index < predictedValues.size() ; index++) {
    auto value = predictedValues.get(index);
    if (value != -1) {  // Skipping all noisy data
      DataVector row(samples.getNcols());
      samples.getRow(index, row);

      row.sub(clusterCentroids[value]);

      averageDistances[value] = averageDistances[value] + row.l2Norm();
    }
  }
  for (auto &label: clusterLabels) {
    averageDistances[label] = averageDistances[label]/static_cast<double>(pointsPerLabel[label]);
  }

  double davidBouldinIndex = 0;
  for (auto labeli : clusterLabels) {
    double maxIndex = 0;
    DataVector tempCentroid(clusterCentroids[labeli]);

    for (auto labelj : clusterLabels) {
      if (labeli == labelj) {
        continue;
      } else {
        // to calculate distance between centroids
        tempCentroid.sub(clusterCentroids[labelj]);
        double index =
          (averageDistances[labeli] + averageDistances[labelj])/tempCentroid.l2Norm();

        if (index > maxIndex) {
          maxIndex = index;
        }
      }
    }
    davidBouldinIndex+=maxIndex;
  }
  davidBouldinIndex = davidBouldinIndex/clusterLabels.size()+numberNoise/4.0;

  return davidBouldinIndex;
}


}//  namespace datadriven
} //  namespace sgpp
