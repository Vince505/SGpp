// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include<cstddef>

namespace sgpp {
namespace datadriven {

struct VisualizationParameters {
  /**
   * The perplexity to use in case tsne is the selected algorithm
   */
  double perplexity = 30;

  /**
   * The theta parameter to use in case tsne is the selected algorithm
   */
  double theta = 0.5;

  /*
   * The random seed to initialize the selected algorithm
   */
  std::size_t seed = 100;

  /*
   * The maximum number of iteration to run on the gradient descent of a selected
   * algorithm
   */
  std::size_t maxNumberIterations = 1000;
};
}  // namespace datadriven
}  // namespace sgpp
