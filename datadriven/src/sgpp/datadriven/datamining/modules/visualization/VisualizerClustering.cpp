// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/visualization/VisualizerClustering.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingClustering.hpp>

#include <sgpp/base/tools/json/JSON.hpp>
#include <sgpp/base/tools/json/ListNode.hpp>
#include <sgpp/base/tools/json/DictNode.hpp>

#include <omp.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <queue>

using json::JSON;
namespace sgpp {
namespace datadriven {
  VisualizerClustering::VisualizerClustering(VisualizerConfiguration config) {
    this->config = config;
    this->visualizerDensityEstimation = new VisualizerDensityEstimation(config);

  }

  void VisualizerClustering::runVisualization(ModelFittingBase &model, DataSource &dataSource,
      size_t epoch, size_t fold, size_t batch) {

    if (batch % config.getGeneralConfig().numBatches != 0 ||
        !config.getGeneralConfig().execute) {
      return;
    }

    size_t nDimensions = model.getDataset()->getDimension();
    if (epoch == 0 && fold == 0 && batch == 0) {
      //originalData = dataSource.getAllSamples()->getData();
      resolution = static_cast<size_t>(pow(2,
                                           model.getFitterConfiguration().getGridConfig().level_+2));
      visualizerDensityEstimation->setResolution(resolution);
    }

    createFolder(config.getGeneralConfig().
      targetDirectory);

    // Creating the output directory
    if (config.getGeneralConfig().crossValidation) {
      currentDirectory = config.getGeneralConfig().
        targetDirectory+"/Epoch_" + std::to_string(epoch);
      createFolder(currentDirectory);
      currentDirectory = config.getGeneralConfig().
        targetDirectory+"/Epoch_" + std::to_string(epoch)+"/Fold_" + std::to_string(fold);
      createFolder(currentDirectory);
      currentDirectory = config.getGeneralConfig().
        targetDirectory+"/Epoch" + std::to_string(epoch)+
                         "/Fold_" + std::to_string(fold) + "/Batch_" + std::to_string(batch);
      createFolder(currentDirectory);

    } else {
      currentDirectory = config.getGeneralConfig().
        targetDirectory+"/Epoch_" + std::to_string(epoch);
      createFolder(currentDirectory);
      currentDirectory = config.getGeneralConfig().
        targetDirectory+"/Epoch_" + std::to_string(epoch)+"/Batch_" + std::to_string(batch);
      createFolder(currentDirectory);
    }

    std::cout << "Creating output directory " << config.getGeneralConfig().targetDirectory
              << std::endl;

    createFolder(currentDirectory+"/Clustering");

    omp_set_num_threads(static_cast<int> (config.getVisualizationParameters().numberCores));

    ModelFittingClustering *clusteringModel =
      dynamic_cast<ModelFittingClustering *>(&model);

    tsneCompressedData = clusteringModel->getPoints();
    #pragma omp parallel sections
    {
      #pragma omp section
      {
        storeTsneJson(tsneCompressedData,
          model, currentDirectory+"/Clustering");
      }
      #pragma omp section
      {
        if (std::find(config.getGeneralConfig().algorithm.begin(),
                      config.getGeneralConfig().algorithm.end(), "linearcuts")
            != config.getGeneralConfig().algorithm.end()) {
          createFolder(currentDirectory+"/DensityEstimation");
          DataMatrix cutMatrix;
          visualizerDensityEstimation->getLinearCuts
            (**(clusteringModel->getDensityEstimationModel()),
              currentDirectory+"/DensityEstimation", cutMatrix, nDimensions);
        }
      }
      #pragma omp section
      {
        if (std::find(config.getGeneralConfig().algorithm.begin(),
                      config.getGeneralConfig().algorithm.end(), "heatmaps")
            != config.getGeneralConfig().algorithm.end()) {
          createFolder(currentDirectory+"/DensityEstimation");
          DataMatrix heatMapMatrix;
          visualizerDensityEstimation->getHeatmap
            (**(clusteringModel->getDensityEstimationModel()),
              currentDirectory+"/DensityEstimation", heatMapMatrix, nDimensions);
        }
      }
    }
      /* To be added after tsne is functioning again
       * #pragma omp section
      {
        if (config.getGeneralConfig().algorithm == "tsne") {
          if (fold == 0 && batch == 0) {
            runTsne(model);
          }
          if (originalData.getNcols() >= 1) {
            DataVector evaluation(originalData.getNrows());
            model.evaluate(originalData, evaluation);
            tsneCompressedData.setColumn(tsneCompressedData.getNcols()-1, evaluation);
            if (config.getGeneralConfig().targetFileType == VisualizationFileType::CSV) {
              CSVTools::writeMatrixToCSVFile(currentDirectory +
                "/tsneCompression", tsneCompressedData);
            } else if (config.getGeneralConfig().targetFileType == VisualizationFileType::json) {
              if (config.getVisualizationParameters().targetDimension != 2) {
                std::cout << "A json output is only available for compressions in 2 dimensions"
                "Storing the CSV instead" << std::endl;
                CSVTools::writeMatrixToCSVFile(currentDirectory +
                  "/tsneCompression", tsneCompressedData);
              }
                storeTsneJson(tsneCompressedData, model, currentDirectory);
            }
          }
        }
      }
    } */
  }

  void VisualizerClustering::storeTsneJson(DataMatrix &matrix, ModelFittingBase &model,
      std::string currentDirectory) {
    JSON jsonOutput;

    ModelFittingClustering *clusteringModel =
        dynamic_cast<ModelFittingClustering *>(&model);

    jsonOutput.addListAttr("data");

    // Trace for the data points
    jsonOutput["data"].addDictValue();
    jsonOutput["data"][0].addIDAttr("type", "\"scatter\"");
    jsonOutput["data"][0].addIDAttr("mode", "\"markers\"");

    DataVector xCol(matrix.getNrows());

    matrix.getColumn(0, xCol);

    jsonOutput["data"][0].addIDAttr("x", xCol.toString());

    DataVector yCol(matrix.getNrows());

    matrix.getColumn(1, yCol);
    jsonOutput["data"][0].addIDAttr("y", yCol.toString());

    /*// Trace for the edges
    jsonOutput["data"].addDictValue();
    jsonOutput["data"][1].addIDAttr("type", "\"scatter\"");
    jsonOutput["data"][1].addIDAttr("mode", "\"lines\"");
    jsonOutput["data"][1].addDictAttr("line");
    jsonOutput["data"][1]["line"].addIDAttr("width", 0.5);


    DataVector source(matrix.getNcols());
    DataVector sink(matrix.getNcols());

    jsonOutput["data"][1].addListAttr("x");
    jsonOutput["data"][1].addListAttr("y");*/

    // Layout of the plot
    jsonOutput.addDictAttr("layout");

    jsonOutput["layout"].addDictAttr("title");

    jsonOutput["layout"]["title"].addIDAttr("text", "\"Hierarchichal Clustering\"");
    jsonOutput["layout"]["title"].addIDAttr("x", 0.5);

    jsonOutput["layout"].addListAttr("updatemenus");
    jsonOutput["layout"]["updatemenus"].addDictValue();
    jsonOutput["layout"]["updatemenus"][0].addIDAttr("type","\"buttons\"");

    jsonOutput["layout"]["updatemenus"][0].addListAttr("buttons");
    jsonOutput["layout"]["updatemenus"][0]["buttons"].addDictValue();
    jsonOutput["layout"]["updatemenus"][0]["buttons"][0].addIDAttr("label","\"Play and try\"");
    jsonOutput["layout"]["updatemenus"][0]["buttons"][0].addIDAttr("method","\"animate\"");

    jsonOutput["layout"]["updatemenus"][0].addIDAttr("showactive", true);

    // FRAMES THIS IS FOR THE ANIMATION
    jsonOutput.addListAttr("frames");
    processHierarchyTree(matrix, jsonOutput, *clusteringModel);
    jsonOutput.serialize(currentDirectory + "/Hierarchy.json");
  }

  void VisualizerClustering::processHierarchyTree(DataMatrix &matrix,
    JSON &jsonOutput, ModelFittingClustering &model) {
    auto tree = model.getHierarchyTree();

    size_t frame = 0;

    double stepIncrease =
      ( model.getFitterConfiguration().getClusteringConfig().maxDensityThreshold -
        model.getFitterConfiguration().getClusteringConfig().minDensityThreshold)
      / model.getFitterConfiguration().getClusteringConfig().steps;

    for (double step = model.getFitterConfiguration().getClusteringConfig().minDensityThreshold;
         step <= model.getFitterConfiguration().getClusteringConfig().maxDensityThreshold;
         step+=stepIncrease) {
      std::cout << "Pocessing step "<<step<<std::endl;
      DataMatrix framePoints(0, matrix.getNcols() + 1);
      jsonOutput["frames"].addDictValue();
      jsonOutput["frames"][frame].addListAttr("data");
      jsonOutput["frames"][frame]["data"].addDictValue();
      jsonOutput["frames"][frame]["data"][0].addIDAttr("mode", "\"markers\"");

      processHierarchyNode(matrix, framePoints, (*(tree))->getRoot(), step);

      DataVector xCol(framePoints.getNrows());
      DataVector yCol(framePoints.getNrows());
      DataVector zCol(framePoints.getNrows());

      framePoints.getColumn(0, xCol);
      framePoints.getColumn(1, yCol);
      framePoints.getColumn(2, zCol);

      jsonOutput["frames"][frame]["data"][0].addIDAttr("x", xCol.toString());
      jsonOutput["frames"][frame]["data"][0].addIDAttr("y", yCol.toString());

      jsonOutput["frames"][frame]["data"][0].addDictAttr("marker");
      jsonOutput["frames"][frame]["data"][0]["marker"].addIDAttr("color", zCol.toString());

      jsonOutput["frames"][frame]["data"][0]["marker"].addIDAttr("colorscale", "\"Viridis\"");


      frame++;
    }

  }

  void VisualizerClustering::processHierarchyNode(DataMatrix &fullMatrix,
    DataMatrix &framePoints, ClusterNode* node, double step) {
    if (node->getDensityThreshold() >= step) {
      for (auto index: node->getVertexIndexes()) {
        DataVector point(fullMatrix.getNcols());
        fullMatrix.getRow(index, point);
        if (node->getDensityThreshold() == step) {
          point.append(node->getClusterLabel());
        } else {
          point.append(node->getParent()->getClusterLabel());
        }
        framePoints.appendRow(point);
      }
    } else {
      for (auto child: node->getChildren()) {
        processHierarchyNode(fullMatrix, framePoints, child, step);
      }
    }
  }

}  // namespace datadriven
}  // namespace sgpp
