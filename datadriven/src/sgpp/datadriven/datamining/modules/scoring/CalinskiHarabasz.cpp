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
#include <cmath>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;

namespace sgpp {
namespace datadriven {

Metric *CalinskiHarabasz::clone() const { return new CalinskiHarabasz(*this);}

double CalinskiHarabasz::measure(const DataVector &predictedValues, const DataVector &trueValues,
                                 const ModelFittingBase &model, Dataset &testDataset) const {
    auto clusteringModel = dynamic_cast<const ModelFittingClustering*>(&model);
    DataMatrix samples = clusteringModel->getPoints();

    DataVector values = clusteringModel->getLabels();

    if (predictedValues.size() != 0) {
      if (samples.getNrows() % testDataset.getData().getNrows() != 0) {
        for (size_t index = 0; index < predictedValues.size(); index++) {
          DataVector tempRow(testDataset.getData().getNcols());
          testDataset.getData().getRow(index, tempRow);
          samples.appendRow(tempRow);
          values.append(predictedValues.get(index));
        }
      }
    }


    std::map<size_t, size_t> clusterLabels;
    std::map<size_t, size_t> pointsPerLabel;
    std::map<size_t, DataVector> clusterCentroids;

    size_t N = 0;
    size_t vectorPosition = 0;
    DataVector globalCentroid(samples.getNcols());
    // Getting the number of labels and the centroids
    for (size_t index = 0; index < values.size() ; index++) {
      auto value = values.get(index);
        if (clusterLabels.find(value) == clusterLabels.end()) {
          clusterLabels[value] = vectorPosition++;
          DataVector centroid(samples.getNcols(), 0);
          clusterCentroids[value] = centroid;
        }

        DataVector row(samples.getNcols());
        samples.getRow(index, row);
        clusterCentroids[value].add(row);
        pointsPerLabel[value]++;
        globalCentroid.add(row);
        N++;
    }

    if (clusterLabels.size() == 1) {
      std::cout << "To use the Calinski-Harabasz score, "
                   "at least 2 clusters are required" << std::endl;
      return 0.0;
    }
    // Getting the global centroid
    globalCentroid.mult(1/ static_cast<double>(values.size()));

    DataVector globalSumSquareVector(clusterLabels.size());
    // Between group dispersion
    DataVector betweenClusterDispersionVector(clusterLabels.size());
    for (auto &label: clusterLabels) {
      clusterCentroids[label.first].mult(1/ static_cast<double>(pointsPerLabel[label.first]));

      DataVector temp(clusterCentroids[label.first]);
      temp.sub(globalCentroid);

      betweenClusterDispersionVector.set(label.second, pointsPerLabel[label.first]*pow(temp.l2Norm(), 2));
    }

   double betweenClusterDispersion = betweenClusterDispersionVector.sum();

   // Within group dispersion
   DataVector withinClusterDispersionVector(clusterLabels.size());

   for (size_t index = 0; index < values.size() ; index++) {
     auto value = values.get(index);
       DataVector row(samples.getNcols());
       samples.getRow(index, row);
       row.sub(clusterCentroids[value]);
       withinClusterDispersionVector.set(clusterLabels[value],
         withinClusterDispersionVector.get(clusterLabels[value]) + pow(row.l2Norm(), 2));
   }

   double withinClusterDispersion = withinClusterDispersionVector.sum();

   std::cout << "Number noisy points " << values.size() - N<<std::endl;
   double score = (betweenClusterDispersion/withinClusterDispersion)*
     ((N-clusterLabels.size())/(clusterLabels.size()-1));

   return score;
}

double CalinskiHarabasz::measureLowerIsBetter(const DataVector &predictedValues,
  const DataVector &trueValues,
  const ModelFittingBase &model, Dataset &testDataset) const {
  return measure(predictedValues, trueValues, model, testDataset);
}
}//  namespace datadriven
} //  namespace sgpp