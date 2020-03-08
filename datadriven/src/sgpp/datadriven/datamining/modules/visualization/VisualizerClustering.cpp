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
      resolution = static_cast<size_t>(pow(2,
        model.getFitterConfiguration().getGridConfig().level_+2));
      visualizerDensityEstimation->setResolution(resolution);
    }

    createFolder(config.getGeneralConfig().
      targetDirectory);

    // Creating the output directory
    if (config.getGeneralConfig().crossValidation) {
      currentDirectory = config.getGeneralConfig().
        targetDirectory+"/Fold_" + std::to_string(fold);
      createFolder(currentDirectory);
      currentDirectory = config.getGeneralConfig().
        targetDirectory+"/Fold_" + std::to_string(fold)+"/Epoch_" + std::to_string(epoch);
      createFolder(currentDirectory);
      currentDirectory = config.getGeneralConfig().
        targetDirectory+"/Fold_" + std::to_string(fold)+"/Epoch_" + std::to_string(epoch)+
                         + "/Batch_" + std::to_string(batch);
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


    omp_set_num_threads(static_cast<int> (config.getGeneralConfig().numberCores));

    ModelFittingClustering *clusteringModel =
      dynamic_cast<ModelFittingClustering *>(&model);

    #pragma omp parallel sections
    {
      #pragma omp section
      {
        if (std::find(config.getGeneralConfig().plots.begin(),
                      config.getGeneralConfig().plots.end(), "linearcuts")
            != config.getGeneralConfig().plots.end()) {
          createFolder(currentDirectory+"/DensityEstimation");
          DataMatrix cutMatrix;
          visualizerDensityEstimation->getLinearCuts
            (**(clusteringModel->getDensityEstimationModel()),
              currentDirectory+"/DensityEstimation", cutMatrix, nDimensions);
        }
      }
      #pragma omp section
      {
        if (std::find(config.getGeneralConfig().plots.begin(),
                      config.getGeneralConfig().plots.end(), "heatmaps")
            != config.getGeneralConfig().plots.end()) {
          createFolder(currentDirectory+"/DensityEstimation");
          DataMatrix heatMapMatrix;
          visualizerDensityEstimation->getHeatmap
            (**(clusteringModel->getDensityEstimationModel()),
              currentDirectory+"/DensityEstimation", heatMapMatrix, nDimensions);
        }
      }
    }
  }

  void VisualizerClustering::runPostProcessingVisualization(ModelFittingBase &model,
        DataSource &dataSource, size_t fold) {
    if (std::find(config.getGeneralConfig().plots.begin(),
                  config.getGeneralConfig().plots.end(), "scatterplots")
        != config.getGeneralConfig().plots.end() && config.getGeneralConfig().execute ) {
      auto currentDirectory = config.getGeneralConfig().targetDirectory;

      if (config.getGeneralConfig().crossValidation) {
        currentDirectory = currentDirectory
                           +"/Fold_" + std::to_string(fold)+"/Clustering";
      } else {
        currentDirectory = currentDirectory+"/Clustering";
      }

      createFolder(currentDirectory);

      ModelFittingClustering *clusteringModel =
        dynamic_cast<ModelFittingClustering *>(&model);

      DataMatrix originalData = clusteringModel->getPoints();
      // TSNE just needs to be run once, in the remaining folds we just append the
      // evaluation of the model
      if (config.getGeneralConfig().algorithm == "tsne" && fold == 0) {
        runTsne(originalData, compressedData);
      }
      clusteringModel->generateSimilarityGraph();
      auto densityEstimationModel = clusteringModel->getDensityEstimationModel();
      if (config.getGeneralConfig().targetFileType == VisualizationFileType::json) {
        DataMatrix toJsonMatrix(compressedData);
        DataVector densityValues(toJsonMatrix.getNrows());
        (*densityEstimationModel)->evaluate(originalData, densityValues);
        toJsonMatrix.appendCol(densityValues);

        storeScatterPlotJson(toJsonMatrix,
                      model, currentDirectory);
      } else {
        storeHierarchyCsv(compressedData, *clusteringModel, currentDirectory);

        DataMatrix toCsvMatrix(compressedData);
        DataVector densityValues(toCsvMatrix.getNrows());
        (*densityEstimationModel)->evaluate(originalData, densityValues);
        toCsvMatrix.appendCol(densityValues);
        CSVTools::writeMatrixToCSVFile(currentDirectory +
                                       "/densityValues", toCsvMatrix);
      }
    }
  }

  void VisualizerClustering::storeScatterPlotJson(DataMatrix &matrix, ModelFittingBase &model,
                                                 std::string currentDirectory) {
    ModelFittingClustering *clusteringModel =
      dynamic_cast<ModelFittingClustering *>(&model);
    #pragma omp parallel sections
    {
      #pragma omp section
      {
        getHierarchyAnimation(matrix,
                              *clusteringModel, currentDirectory);
      }
      #pragma omp section
      {
        getDensityGraphPlot(matrix, *clusteringModel, currentDirectory);
      }
    }
  }

  void VisualizerClustering::getHierarchyAnimation(DataMatrix &matrix,
    ModelFittingClustering &model, std::string currentDirectory) {
    JSON jsonOutput;

    auto tree = model.getHierarchyTree();

    size_t numberClusters = (*tree)->getNumberClusters();

    size_t maxLevel = (*tree)->getNumberLevels();

    // Add empty data traces that can be identified by name in the frames. Otherwise animation
    // won't work. Each cluster, including the noise, is considered as a data trace.
    jsonOutput.addListAttr("data");

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

      separateClustersIntoTraces(matrix, model, zCol, traces, frame);

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

        if (cluster == 0) {
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

    // To keep track where the graph traces start
    size_t graphTraceIndex = traceIndex;
    // Adding the graphs per level
    for (size_t frame = 0; frame <= maxLevel; frame++) {
      for (size_t index = 0; index < frame; index++) {
        jsonOutput["layout"]["sliders"][0]["steps"][frame]["args"][1].addIdValue(false);
      }
      jsonOutput["data"].addDictValue();
      jsonOutput["data"][traceIndex].addIDAttr("mode", "\"lines\"");
      jsonOutput["data"][traceIndex].addIDAttr("name", "\"Graph\"");
      jsonOutput["data"][traceIndex].addDictAttr("line");
      jsonOutput["data"][traceIndex]["line"].addIDAttr("width", 1.0);
      jsonOutput["data"][traceIndex]["line"].addIDAttr("color", "\"navy\"");
      jsonOutput["data"][traceIndex].addListAttr("x");
      jsonOutput["data"][traceIndex].addListAttr("y");
      if (frame == 0) {
        jsonOutput["data"][traceIndex].addIDAttr("visible", true);
      } else {
        jsonOutput["data"][traceIndex].addIDAttr("visible", false);
      }
      jsonOutput["layout"]["sliders"][0]["steps"][frame]["args"][1].addIdValue(true);
      traceIndex++;

      for (size_t index = frame+1; index <=maxLevel; index++) {
        jsonOutput["layout"]["sliders"][0]["steps"][frame]["args"][1].addIdValue(false);
      }
    }

    auto graph = model.getGraph();
    DataVector source(matrix.getNcols());
    DataVector sink(matrix.getNcols());
    std::vector<size_t> visitedVertex;
    for (size_t index = 0; index < matrix.getNrows(); index++) {
      visitedVertex.push_back(index);
      size_t pointLevel =  (*tree)->getMostSpecificLevel(index);
      matrix.getRow(index, source);
      auto neighbors = graph->getAdjacentVertices(index);
      for (size_t neighbor : neighbors) {
        size_t neighborLevel = (*tree)->getMostSpecificLevel(neighbor);
        for (size_t level = 0; level <= pointLevel; level++) {
          if (neighborLevel >= pointLevel) {
            jsonOutput["data"][graphTraceIndex + level]["x"].addIdValue(source.get(0));
            jsonOutput["data"][graphTraceIndex + level]["y"].addIdValue(source.get(1));
            matrix.getRow(neighbor, sink);
            jsonOutput["data"][graphTraceIndex + level]["x"].addIdValue(sink.get(0));
            jsonOutput["data"][graphTraceIndex + level]["y"].addIdValue(sink.get(1));
            jsonOutput["data"][graphTraceIndex + level]["x"].addIdValue("\"None\"\n");
            jsonOutput["data"][graphTraceIndex + level]["y"].addIdValue("\"None\"\n");
          }
        }
      }
    }

    std::cout << "Writing file " << currentDirectory + "/Hierarchy.json" << std::endl;
    jsonOutput.serialize(currentDirectory + "/Hierarchy.json");
  }

  void VisualizerClustering::getDensityGraphPlot(DataMatrix &matrix, ModelFittingClustering &model,
    std::string currentDirectory) {
    JSON jsonOutput;

    jsonOutput.addListAttr("data");
    // Layout of the plot
    jsonOutput.addDictAttr("layout");

    jsonOutput["layout"].addDictAttr("title");

    jsonOutput["layout"]["title"].addIDAttr("text", "\"Graph and Densities\"");
    jsonOutput["layout"]["title"].addIDAttr("x", 0.5);

    // Trace for the edges of the graph
    jsonOutput["data"].addDictValue();
    jsonOutput["data"][0].addIDAttr("type", "\"scatter\"");
    jsonOutput["data"][0].addIDAttr("mode", "\"lines\"");
    jsonOutput["data"][0].addIDAttr("name", "\"Graph\"");
    jsonOutput["data"][0].addDictAttr("line");
    jsonOutput["data"][0]["line"].addIDAttr("width", 0.5);

    DataVector source(matrix.getNcols());
    DataVector sink(matrix.getNcols());

    jsonOutput["data"][0].addListAttr("x");
    jsonOutput["data"][0].addListAttr("y");

    auto graph = model.getGraph();
    std::vector<size_t> visitedVertex;
    for (size_t index = 0; index < matrix.getNrows(); index++) {
      visitedVertex.push_back(index);
      matrix.getRow(index, source);
      if (graph->containsVertex(index)) {
        auto neighbors = graph->getAdjacentVertices(index);
        for (size_t neighbor : neighbors) {
          // We skip points which were already linked previously to avoid edges repetitions
          if (std::find(visitedVertex.begin(), visitedVertex.end(), neighbor)
              != visitedVertex.end()) {
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
    }

    // Adding the data points¿
    jsonOutput["data"].addDictValue();
    jsonOutput["data"][1].addIDAttr("type", "\"scatter\"");
    jsonOutput["data"][1].addIDAttr("mode", "\"markers\"");
    jsonOutput["data"][1].addIDAttr("showlegend", false);
    DataVector xCol(matrix.getNrows());
    matrix.getColumn(0, xCol);
    jsonOutput["data"].addDictValue();
    jsonOutput["data"][1].addIDAttr("x", xCol.toString());

    DataVector yCol(matrix.getNrows());
    matrix.getColumn(1, yCol);
    jsonOutput["data"][1].addIDAttr("y", yCol.toString());

    DataVector zCol(matrix.getNrows());

    matrix.getColumn(2, zCol);

    jsonOutput["data"][1].addDictAttr("marker");
    jsonOutput["data"][1].addIDAttr("hovertext", zCol.toString());
    jsonOutput["data"][1].addIDAttr("hoverinfo", "\"x+y+text\"");
    jsonOutput["data"][1]["marker"].addIDAttr("color", zCol.toString());
    jsonOutput["data"][1]["marker"].addIDAttr("colorscale", "\"Viridis\"");
    jsonOutput["data"][1]["marker"].addIDAttr("opacity", 0.8);
    jsonOutput["data"][1]["marker"].addIDAttr("showscale", true);
    jsonOutput["data"][1]["marker"].addDictAttr("colorbar");
    jsonOutput["data"][1]["marker"]["colorbar"].addDictAttr("title");
    jsonOutput["data"][1]["marker"]["colorbar"]["title"].addIDAttr("text", "\"Density value\"");

    // Places graph legend on the left side of the plot
    jsonOutput["layout"].addDictAttr("legend");
    jsonOutput["layout"]["legend"].addIDAttr("x", -0.1);
    jsonOutput["layout"]["legend"].addIDAttr("y", 1.0);

    std::cout << "Writing file " << currentDirectory + "/Graph.json" << std::endl;
    jsonOutput.serialize(currentDirectory + "/Graph.json");
  }


void VisualizerClustering::separateClustersIntoTraces(DataMatrix &points,
  ModelFittingClustering &model, DataVector &labels,
  std::vector<DataMatrix> &traces, size_t level) {

    auto tree = model.getHierarchyTree();
    for (size_t cluster = 0 ; cluster < traces.size(); cluster++) {
      traces[cluster].resize(0, points.getNcols());
    }

    for (size_t index = 0 ; index < labels.size(); index++) {
      if (level >= 2) {
        size_t pointLevel = (*tree)->getMostSpecificLevel(index);
        if (pointLevel <= level - 2) {
          continue;
        }
      }
      DataVector traceRow(points.getNcols());
      points.getRow(index, traceRow);
      traces[static_cast<int>(labels.get(index))+1].appendRow(traceRow);
    }
  }

void VisualizerClustering::storeHierarchyCsv(DataMatrix &matrix,
  ModelFittingClustering &model, std::string currentDirectory) {
    auto tree = model.getHierarchyTree();
    size_t maxLevel = (*tree)->getNumberLevels();
    DataMatrix toCsvMatrix(matrix);
    DataVector labels(toCsvMatrix.getNrows(), 0.0);
    toCsvMatrix.appendCol(labels);
    for (size_t level = 0; level <= maxLevel; level++) {
      (*tree)->evaluateClusteringAtLevel(labels, level);
      toCsvMatrix.setColumn(2, labels);
      CSVTools::writeMatrixToCSVFile(currentDirectory +
                                     "/clustering_level_"+std::to_string(level), toCsvMatrix);
    }
  }
}  // namespace datadriven
}  // namespace sgpp
