// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/base/grid/generation/functors/RefinementFunctor.hpp>

#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimation.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingClassification.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationClustering.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationDensityEstimation.hpp>
#include <sgpp/datadriven/datamining/tools/Graph.hpp>
#include <sgpp/datadriven/datamining/tools/hierarchyTree/HierarchyTree.hpp>
#include <sgpp/datadriven/datamining/tools/vpTree/VpTree.hpp>

#include <map>
#include <vector>

namespace sgpp {
namespace datadriven {

/**
 * Fitter object that encapsulates density based classification using instances of
 * ModelFittingDensityEstimation for each class.
 */
class ModelFittingClustering : public ModelFittingBase {
 public:
  /**
   * Constructor
   *
   * @param config configuration object that specifies grid, refinement, and regularization
   */
  explicit ModelFittingClustering(const FitterConfigurationClustering& config);

  ~ModelFittingClustering() = default;

  /**
   * Runs the series of stpes to obtain the clustering model of the dataset based on the
   * algorithm using a density estimation model.
   * Requires only data samples and no targets (since those are irrelevant for the clustering)
   * @param newDataset the training dataset that is used to fit the model.
   */
  void fit(Dataset& newDataset) override;

  /**
   * Updates the internal the density estimation model with new samples
   * @param newDataset Training datasetused to update the model
   */
  void update(Dataset& newDataset) override;


  double evaluate(const DataVector& sample) override;

  /**
   * Evaluates a series of samples and assign them a cluster label if the point
   * was used in the clustering. Otherwise it will be labeled as noise
   * @param samples Samples to evaluate
   * @param results Labels assigned to the samples
   */
  void evaluate(DataMatrix& samples, DataVector& results) override;

  /**
   * Performs a refinement given the new grid size and the points to coarsened
   * @return if the grid was refined (true)
   */
  bool refine() override;

  /**
   * Clears the model
   */
  void reset() override;

  /**
   * Resets any trained representations of the model, but does not reset the entire state.
   */
  void resetTraining() override;

  /**
  * Should compute some kind of Residual to evaluate the fit of the model.
  *
  * In the case of density estimation, this is
  * || R * alpha_lambda - b_val ||_2
  *
  * This is useful for unsupervised learning models, where normal evaluation cannot be used as
  * there are no targets.
  *
  * @param validationData Matrix for validation data
  *
  * @returns the residual score
  */
  double computeResidual(DataMatrix &validationData) const override;

  /**
   * Updates the regularization parameter lambda of the underlying model.
   *
   * @param lambda the new lambda parameter
   */
  void updateRegularization(double lambda) override;

  /**
   * Returns the pointer to the pointer pointing to the density estimation model
   * @return Pointer of the unique pointer of the density estimation model
   */
  std::unique_ptr<ModelFittingDensityEstimation>* getDensityEstimationModel();

  /**
   * Returns the pointer to the graph
   * @return Pointer to the graph
   */
  std::shared_ptr<Graph> getGraph();

  /**
   * Returns the pointer to the pointer pointing to the hierarchy tree data structure
   * @return Pointer of the unique pointer to the hierarchy tree
   */
  std::unique_ptr<HierarchyTree>* getHierarchyTree();

  /**
   * Initializes the instances of the hierarchy tree
   */
  void intializeHierarchyTree();

  /**
   * Saves a copy of the current graph in the variable prunedGraphPreviousStep
   */
  void copyPreviousGraphStep();

  /**
   * Method which updates the vpTree with new samples .
   * @param newDataset New dataset to be added
   */
  void updateVpTree(DataMatrix &newDataset);

  /**
   * Method that generates the similarity graph of a given dataset based on the nearest
   * neighbors hyperparameter given in the configuration
   *
   **/
  void generateSimilarityGraph();

  /**
   * Method that deletes from the graph all of the points that do not have the minimum density value
   * @params densityThreshold minimum density threshold used for deletion
   */
  void applyDensityThresholds(double densityThreshold);

  /**
   * Method that detects all of the disconected components from the graph and assigns a label
   * to each of its corresponding points
   * @params clusterMap A map mapping the descriptor of a vertex to a certain component
   */
  void detectComponentsAndLabel(std::map<UndirectedGraph::vertex_descriptor, size_t> &clusterMap);

  /**
   * Method which  updates the hierachy tree after obtaining the connected components
   * @param clusterMap A map mapping the descriptor of a vertex to a certain component
   * @param densityThreshold  densityThreshold minimum density threshold used for splitting a child
   */
  void getHierarchy(std::map<UndirectedGraph::vertex_descriptor, size_t> &clusterMap,
    double densityThreshold);

  /**
   * Stores the info of the hierarchy
   * @param outputDirectory Directory to store the file with the hierarchy info
   */
  void storeHierarchy(std::string outputDirectory);

  /**
   * Method to obtain the points used in the clustering
   * @return Matrix with all of the points used in the clustering
   */
  DataMatrix& getPoints() const;

 protected:
  /**
   * Count the amount of refinement operations performed on the current dataset.
   */
  size_t refinementsPerformed;

 private:
  /**
  * Density Estimation Model
  */
  std::unique_ptr<ModelFittingDensityEstimation> densityEstimationModel;

  /**
   * VpTree to do nearest neighbors queries
   */
  std::unique_ptr<VpTree> vpTree;

  /**
   * Nearest Neighbors' Graph
   */
  std::shared_ptr<Graph> graph;

  /**
   * Graph that keeps a copy of the nearest neighbors graph during the hierarchy generation
   */
  std::shared_ptr<Graph> prunedGraphPreviousStep;

  /**
   * Hierarchy Tree containing the hierarchy of clusters
   */
  std::unique_ptr<HierarchyTree> hierarchy;

  /**
   * Method that generates the density estimation model of the unlabeled dataset
   * @params dataset the training dataset that is used to fit the model.
   */
  void generateDensityEstimationModel(Dataset &dataset);

  /**
   * Creates a density estimation model that fits the model settings.
   * @param densityEstimationConfig configuration for the density estimation
   * @return a new density estimation model
   */
  std::unique_ptr<ModelFittingDensityEstimation> createNewDensityModel(
      sgpp::datadriven::FitterConfigurationDensityEstimation& densityEstimationConfig);
};
}  // namespace datadriven
}  // namespace sgpp
