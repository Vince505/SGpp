// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/ScorerMetricTypeParser.hpp>

#include <algorithm>
#include <map>
#include <string>

namespace sgpp {
namespace datadriven {

using sgpp::base::data_exception;

const ScorerMetricTypeParser::MetricTypeMap_t ScorerMetricTypeParser::metricTypeMap = []() {
  return MetricTypeMap_t{std::make_pair(ScorerMetricType::mse, "MSE"),
                         std::make_pair(ScorerMetricType::nll, "NLL"),
                         std::make_pair(ScorerMetricType::accuracy, "Accuracy"),
                         std::make_pair(ScorerMetricType::fowlkes_mallows, "Fowlkes-Mallows"),
                         std::make_pair(ScorerMetricType::v_measure, "V-Measure"),
                         std::make_pair(ScorerMetricType::david_bouldin, "David-Bouldin"),
                         std::make_pair(ScorerMetricType::calinski_harabasz, "Calinski-Harabasz")};
}();

const ScorerMetricTypeParser::RegularizationTypeMap_t
    ScorerMetricTypeParser::regularizationTypeMap = []() {
      return RegularizationTypeMap_t{
          std::make_pair(RegularizationMetricType::mse, "MSE"),
          std::make_pair(RegularizationMetricType::nll, "NLL"),
          std::make_pair(RegularizationMetricType::accuracy, "Accuracy"),
          std::make_pair(RegularizationMetricType::residual, "Residual")};
    }();

const std::string& ScorerMetricTypeParser::toString(ScorerMetricType type) {
  return metricTypeMap.at(type);
}
ScorerMetricType ScorerMetricTypeParser::parse(const std::string& input) {
  auto inputLower = input;
  std::transform(inputLower.begin(), inputLower.end(), inputLower.begin(), ::tolower);

  if (inputLower == "mse") {
    return ScorerMetricType::mse;
  } else if (inputLower == "nll") {
    return ScorerMetricType::nll;
  } else if (inputLower == "accuracy") {
    return ScorerMetricType::accuracy;
  } else if (inputLower == "fowlkes-mallows") {
    return ScorerMetricType::fowlkes_mallows;
  } else if (inputLower == "v-measure") {
    return ScorerMetricType::v_measure;
  } else if (inputLower == "calinski-harabasz") {
    return ScorerMetricType::calinski_harabasz;
  } else if (inputLower == "david-bouldin") {
    return ScorerMetricType::david_bouldin;
  } else {
    const auto errorMsg = "Failed to convert string \"" + input + "\" to any known ScorerMetric";
    throw data_exception(errorMsg.c_str());
  }
}

RegularizationMetricType ScorerMetricTypeParser::parseRegularizationMetric(
    const std::string& input) {
  auto inputLower = input;
  std::transform(inputLower.begin(), inputLower.end(), inputLower.begin(), ::tolower);

  if (inputLower == "mse") {
    return RegularizationMetricType::mse;
  } else if (inputLower == "nll") {
    return RegularizationMetricType::nll;
  } else if (inputLower == "accuracy") {
    return RegularizationMetricType::accuracy;
  } else if (inputLower == "residual") {
    return RegularizationMetricType::residual;
  } else {
    const auto errorMsg =
        "Failed to convert string \"" + input + "\" to any known RegularizationMetricType";
    throw data_exception(errorMsg.c_str());
  }
}

const std::string& ScorerMetricTypeParser::regularizationMetricToString(
    RegularizationMetricType type) {
  return regularizationTypeMap.at(type);
}

} /* namespace datadriven */
} /* namespace sgpp */
