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

        getGraphPlot(tsneCompressedData, *clusteringModel, currentDirectory+"/Clustering");


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

    auto tree = clusteringModel->getHierarchyTree();

    size_t numberClusters = (*tree)->getNumberClusters();

    size_t maxLevel = (*tree)->getNumberLevels();

    // Add empty data traces that can be identified by name in the frames. Otherwise animation
    // won't work. Each cluster, including the noise, is considered as a data trace.
    jsonOutput.addListAttr("data");
    for (size_t cluster = 0; cluster <= numberClusters; cluster++) {

    }
    // Layout of the plot
    jsonOutput.addDictAttr("layout");

    jsonOutput["layout"].addDictAttr("title");

    jsonOutput["layout"]["title"].addIDAttr("text", "\"Hierarchichal Clustering\"");
    jsonOutput["layout"]["title"].addIDAttr("x", 0.5);

    // ANIMATION SLIDERS
    jsonOutput["layout"].addListAttr("sliders");
    jsonOutput["layout"]["sliders"].addDictValue();
    jsonOutput["layout"]["sliders"][0].addIDAttr("visible", true);
    jsonOutput["layout"]["sliders"][0].addIDAttr("active", static_cast<size_t>(0));
    jsonOutput["layout"]["sliders"][0].addIDAttr("x", 0.0);
    jsonOutput["layout"]["sliders"][0].addIDAttr("y", static_cast<size_t>(0));
    jsonOutput["layout"]["sliders"][0].addIDAttr("xanchor", "\"left\"");
    jsonOutput["layout"]["sliders"][0].addIDAttr("yanchor", "\"top\"");

    jsonOutput["layout"]["sliders"][0].addDictAttr("pad");
    jsonOutput["layout"]["sliders"][0]["pad"].addIDAttr("t", static_cast<size_t>(50));
    jsonOutput["layout"]["sliders"][0]["pad"].addIDAttr("r", static_cast<size_t>(0));
    jsonOutput["layout"]["sliders"][0]["pad"].addIDAttr("b", static_cast<size_t>(10));
    jsonOutput["layout"]["sliders"][0]["pad"].addIDAttr("l", static_cast<size_t>(0));

    jsonOutput["layout"]["sliders"][0].addListAttr("steps");

    // FRAMES THIS IS FOR THE ANIMATION
    jsonOutput.addListAttr("frames");

    std::vector<DataMatrix> traces(numberClusters+1);

    // Each frame is a level
    size_t traceIndex = 0;
    for (size_t frame = 0; frame <= maxLevel; frame++) {
      jsonOutput["layout"]["sliders"][0]["steps"].addDictValue();
      jsonOutput["layout"]["sliders"][0]["steps"][frame].addIDAttr("visible", true);
      jsonOutput["layout"]["sliders"][0]["steps"][frame].addIDAttr("method", "\"restyle\"");
      jsonOutput["layout"]["sliders"][0]["steps"][frame].addIDAttr("label",
                                                                   "\"Level: "+std::to_string(frame)+"\"");

      jsonOutput["layout"]["sliders"][0]["steps"][frame].addListAttr("args");
      jsonOutput["layout"]["sliders"][0]["steps"][frame]["args"].addIdValue(
        "\"visible\"");

      jsonOutput["layout"]["sliders"][0]["steps"][frame]["args"].addListValue();

      // This is to indicate plotlyy in an array which traces are to be activated when
      // we move the slider to a certain value, all of the previously processed traces of
      // previous are set to false
      for (size_t index = 0; index <traceIndex; index++) {
        jsonOutput["layout"]["sliders"][0]["steps"][frame]["args"][1].addIdValue(false);
      }

      DataVector zCol(matrix.getNcols());

      (*tree)->evaluateClusteringAtLevel(zCol, frame);
      separateClustersIntoTraces(matrix, zCol, traces);

      // Per frame Each cluster is processed as a separate component in a trace
      for (size_t cluster = 0; cluster <= numberClusters; cluster++) {

        jsonOutput["data"].addDictValue();
        jsonOutput["data"][traceIndex].addIDAttr("mode", "\"markers\"");

        // We only show at the beginning the first level of the tree
        if (frame == 0) {
          jsonOutput["data"][traceIndex].addIDAttr("visible", true);
        } else {
          jsonOutput["data"][traceIndex].addIDAttr("visible", false);
        }

        if(cluster == 0) {
          if (frame == 0) {
            jsonOutput["data"][traceIndex].addIDAttr("name", "\"Unclustered\"");
          } else {
            jsonOutput["data"][traceIndex].addIDAttr("name", "\"Noise\"");
          }

        } else {
          jsonOutput["data"][traceIndex].addIDAttr("name", "\"Cluster: " +
                                                           std::to_string(cluster - 1) +
                                                           "\"");
        }
        jsonOutput["data"][traceIndex].addIDAttr("showlegend", true);
        DataVector xCol(traces[cluster].getNrows());
        traces[cluster].getColumn(0, xCol);
        jsonOutput["data"][traceIndex].addIDAttr("x", xCol.toString());

        DataVector yCol(traces[cluster].getNrows());
        traces[cluster].getColumn(1, yCol);
        jsonOutput["data"][traceIndex].addIDAttr("y", yCol.toString());

        jsonOutput["layout"]["sliders"][0]["steps"][frame]["args"][1].addIdValue(true);
        traceIndex++;
      }


      for (size_t index = traceIndex; index <(maxLevel+1)*(numberClusters+1); index++) {
        jsonOutput["layout"]["sliders"][0]["steps"][frame]["args"][1].addIdValue(false);
      }

    }

    jsonOutput.serialize(currentDirectory + "/Hierarchy.json");
    std::cout << "Writing file " << currentDirectory + "/Hierarchy.json" << std::endl;
  }

  void VisualizerClustering::getGraphPlot(DataMatrix &matrix, ModelFittingClustering &model,
    std::string currentDirectory) {
    JSON jsonOutput;

    jsonOutput.addListAttr("data");
    // Layout of the plot
    jsonOutput.addDictAttr("layout");

    jsonOutput["layout"].addDictAttr("title");

    jsonOutput["layout"]["title"].addIDAttr("text", "\"Graph and Clustering\"");
    jsonOutput["layout"]["title"].addIDAttr("x", 0.5);

    // Trace for the edges of the graph
    jsonOutput["data"].addDictValue();
    jsonOutput["data"][0].addIDAttr("type", "\"scatter\"");
    jsonOutput["data"][0].addIDAttr("mode", "\"lines\"");
    jsonOutput["data"][0].addIDAttr("name", "\"graph\"");
    jsonOutput["data"][0].addDictAttr("line");
    jsonOutput["data"][0]["line"].addIDAttr("width", 0.5);

    DataVector source(matrix.getNcols());
    DataVector sink(matrix.getNcols());

    jsonOutput["data"][0].addListAttr("x");
    jsonOutput["data"][0].addListAttr("y");

    auto graph = model.getGraph();

    for (size_t index = 0; index < matrix.getNrows(); index++) {

      matrix.getRow(index, source);
      if (graph->containsVertex(index)) {
        auto neighbors = graph->getAdjacentVertices(index);
        for(size_t neighbor: neighbors) {
          jsonOutput["data"][0]["x"].addIdValue(source.get(0));
          jsonOutput["data"][0]["y"].addIdValue(source.get(1));
          matrix.getRow(neighbor, sink);
          jsonOutput["data"][0]["x"].addIdValue(sink.get(0));
          jsonOutput["data"][0]["y"].addIdValue(sink.get(1));
          jsonOutput["data"][0]["x"].addIdValue("\"None\"\n");
          jsonOutput["data"][0]["y"].addIdValue("\"None\"\n");
        }

     }
   }
    // Adding the data points

    auto tree = model.getHierarchyTree();

    size_t numberClusters = (*tree)->getNumberClusters();

    std::vector<DataMatrix> traces(numberClusters+1);

    DataVector zCol(matrix.getNcols());
    model.evaluate(matrix, zCol);

    separateClustersIntoTraces(matrix, zCol, traces);

    for(size_t clusterNumber = 0 ; clusterNumber <= numberClusters; clusterNumber++) {
      jsonOutput["data"].addDictValue();
      jsonOutput["data"][clusterNumber + 1].addIDAttr("type", "\"scatter\"");
      jsonOutput["data"][clusterNumber + 1].addIDAttr("mode", "\"markers\"");

      if (clusterNumber == 0) {
        jsonOutput["data"][clusterNumber + 1].addIDAttr("name", "\"Noise\"");
      } else {
        jsonOutput["data"][clusterNumber + 1].addIDAttr("name", "\"Cluster: " +
                                                                std::to_string(clusterNumber - 1) +
                                                                "\"");
      }

      DataVector xCol(traces[clusterNumber].getNrows());
      traces[clusterNumber].getColumn(0, xCol);
      jsonOutput["data"][clusterNumber + 1].addIDAttr("x", xCol.toString());

      DataVector yCol(traces[clusterNumber].getNrows());
      traces[clusterNumber].getColumn(1, yCol);
      jsonOutput["data"][clusterNumber + 1].addIDAttr("y", yCol.toString());
    }
    jsonOutput.serialize(currentDirectory + "/Graph.json");

    std::cout << "Writing file " << currentDirectory + "/Graph.json" << std::endl;
  }


void VisualizerClustering::separateClustersIntoTraces(DataMatrix &points,
  DataVector &labels, std::vector<DataMatrix> &traces) {

    for(size_t cluster = 0 ; cluster < traces.size(); cluster++) {
      traces[cluster].resize(0, points.getNcols());
    }

    for(size_t index = 0 ; index<labels.size(); index++) {
      DataVector traceRow(points.getNcols());
      points.getRow(index, traceRow);
      traces[static_cast<int>(labels.get(index))+1].appendRow(traceRow);
    }

  }
}  // namespace datadriven
}  // namespace sgpp
