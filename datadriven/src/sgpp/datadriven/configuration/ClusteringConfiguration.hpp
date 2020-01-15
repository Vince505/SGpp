// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

/**
 * Struct that stores all the configuration information
 * for an offline object for density based clustering
 */
namespace sgpp {
namespace datadriven {

struct ClusteringConfiguration {
    // nearest neighbors parameter to build the similiarity graph
    size_t noNearestNeighbors = 5;

    // density threshold minimum prune the graph
    double minDensityThreshold = 0.1;

    // density threshold minimum prune the graph
    double maxDensityThreshold = 0.9;

    // Number of steps to run the hierarchy
    size_t steps = 0;
};
}  // namespace datadriven
}  // namespace sgpp