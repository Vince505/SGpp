// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <boost/test/test_tools.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <sgpp/datadriven/algorithm/DBMatObjectStore.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineOrthoAdapt.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDEFactory.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDEOrthoAdapt.hpp>
#include <sgpp/datadriven/algorithm/DBMatPermutationFactory.hpp>
#include <sgpp/datadriven/algorithm/GridFactory.hpp>

#include <vector>

BOOST_AUTO_TEST_SUITE(PermutationBlowUpTest)

BOOST_AUTO_TEST_CASE(ComponentGridOrthoTest) {
  sgpp::base::GeneralGridConfiguration baseGridConfig;
  baseGridConfig.generalType_ = sgpp::base::GeneralGridType::ComponentGrid;
  baseGridConfig.type_ = sgpp::base::GridType::Linear;
  baseGridConfig.levelVector_ = std::vector<size_t>{3, 2, 2};
  baseGridConfig.dim_ = 3;

  sgpp::base::GeneralGridConfiguration desiredGridConfig;
  desiredGridConfig.generalType_ = sgpp::base::GeneralGridType::ComponentGrid;
  desiredGridConfig.type_ = sgpp::base::GridType::Linear;
  desiredGridConfig.levelVector_ = std::vector<size_t>{2, 3, 2, 1, 1};
  desiredGridConfig.dim_ = 5;

  sgpp::base::AdaptivityConfiguration adaptivityConfig;
  adaptivityConfig.numRefinements_ = 0;

  sgpp::datadriven::RegularizationConfiguration regConfig;
  regConfig.lambda_ = 0;

  sgpp::datadriven::DensityEstimationConfiguration densConfig;
  densConfig.type_ = sgpp::datadriven::DensityEstimationType::Decomposition;
  densConfig.decomposition_ = sgpp::datadriven::MatrixDecompositionType::OrthoAdapt;

  // build grids
  std::cout << "Build Grid" << std::endl;
  sgpp::datadriven::GridFactory gridFactory;
  std::vector<std::vector<size_t>> interactions = std::vector<std::vector<size_t>>();
  std::unique_ptr<sgpp::base::Grid> desiredGrid{
      gridFactory.createGrid(desiredGridConfig, interactions)};

  // build offline objects
  std::cout << "Instanciate Offline" << std::endl;

  std::unique_ptr<sgpp::datadriven::DBMatOfflineOrthoAdapt> desiredOff{
      dynamic_cast<sgpp::datadriven::DBMatOfflineOrthoAdapt *>(
          sgpp::datadriven::DBMatOfflineFactory::buildOfflineObject(
              desiredGridConfig, adaptivityConfig, regConfig, densConfig))};

  // object store
  std::shared_ptr<sgpp::datadriven::DBMatObjectStore> store =
      std::make_shared<sgpp::datadriven::DBMatObjectStore>();

  // build and decompose
  std::cout << "Build and decompose" << std::endl;
  // sgpp::datadriven::DBMatBaseObjectStore store2(adaptivityConfig, regConfig, densConfig);
  // Add base to store
  sgpp::datadriven::GeometryConfiguration gc;
  gc.stencilType = sgpp::datadriven::StencilType::None;
  sgpp::datadriven::DBMatPermutationFactory factory(store);
  // put base object into store
  factory.getPermutedObject(baseGridConfig, gc, adaptivityConfig, regConfig, densConfig);
  // build and decompose desired offline object
  desiredOff->buildMatrix(desiredGrid.get(), regConfig);
  desiredOff->decomposeMatrix(regConfig, densConfig);
  // get permuted base object from factory
  std::unique_ptr<sgpp::datadriven::DBMatOfflineOrthoAdapt> permOff(
      (sgpp::datadriven::DBMatOfflineOrthoAdapt *)factory.getPermutedObject(
          desiredGridConfig, gc, adaptivityConfig, regConfig, densConfig));

  // build online objects
  std::unique_ptr<sgpp::datadriven::DBMatOnlineDEOrthoAdapt> online1{
      dynamic_cast<sgpp::datadriven::DBMatOnlineDEOrthoAdapt *>(
          sgpp::datadriven::DBMatOnlineDEFactory::buildDBMatOnlineDE(
              *permOff, *desiredGrid, regConfig.lambda_, 0,
              sgpp::datadriven::MatrixDecompositionType::OrthoAdapt))};
  std::unique_ptr<sgpp::datadriven::DBMatOnlineDEOrthoAdapt> online2{
      dynamic_cast<sgpp::datadriven::DBMatOnlineDEOrthoAdapt *>(
          sgpp::datadriven::DBMatOnlineDEFactory::buildDBMatOnlineDE(
              *desiredOff, *desiredGrid, regConfig.lambda_, 0,
              sgpp::datadriven::MatrixDecompositionType::OrthoAdapt))};

  // Generate sample dataset
  std::cout << "Generate samples" << std::endl;
  sgpp::base::DataMatrix samples(100, desiredGridConfig.dim_);
  for (size_t i = 0; i < 100; i++) {
    sgpp::base::DataVector vec(desiredGridConfig.dim_);
    for (size_t j = 0; j < desiredGridConfig.dim_; j++) {
      vec.at(j) = (double)std::rand() / RAND_MAX;
    }
    samples.setRow(i, vec);
  }

  // alphas
  sgpp::base::DataVector alpha1(permOff->getGridSize());
  sgpp::base::DataVector alpha2(desiredOff->getGridSize());

  // compute alphas
  online1->computeDensityFunction(alpha1, samples, *desiredGrid, densConfig, false);
  online2->computeDensityFunction(alpha2, samples, *desiredGrid, densConfig, false);

  for (size_t i = 0; i < alpha1.getSize(); i++) {
    BOOST_CHECK(std::abs(alpha1[i] - alpha2[i]) / std::abs(alpha2[i]) < 0.001);
  }
}

BOOST_AUTO_TEST_CASE(FullCombiSchemeOrthoTest) {}

BOOST_AUTO_TEST_SUITE_END()
