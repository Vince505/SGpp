// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/visualization/VisualizerDensityEstimation.hpp>
#include <sgpp/base/tools/json/ListNode.hpp>
#include <sgpp/base/tools/json/DictNode.hpp>
#include <sgpp/base/tools/json/JSON.hpp>
#include <omp.h>
#include <iostream>
#include <memory>
#include <vector>
#include <algorithm>
#include <string>

namespace sgpp {
namespace datadriven {
VisualizerDensityEstimation::VisualizerDensityEstimation(VisualizerConfiguration config) {
  this->config = config;
}

void VisualizerDensityEstimation::runVisualization(ModelFittingBase &model, DataSource &dataSource,
  size_t fold, size_t batch) {
  if (batch % config.getGeneralConfig().numBatches != 0 ||
    !config.getGeneralConfig().execute) {
    return;
  }
  std::cout << "Creating output directory " << config.getGeneralConfig().targetDirectory
     << std::endl;
  createFolder(config.getGeneralConfig().
      targetDirectory);
  // Creating the output directory
  if (config.getGeneralConfig().crossValidation) {
    currentDirectory = config.getGeneralConfig().
    targetDirectory+"/Fold_" + std::to_string(fold);
    createFolder(currentDirectory);
    currentDirectory = config.getGeneralConfig().
        targetDirectory+"/Fold_" + std::to_string(fold) + "/Batch_" + std::to_string(batch);
    createFolder(currentDirectory);

  } else {
    currentDirectory = config.getGeneralConfig().
    targetDirectory+"/Batch_" + std::to_string(batch);
    createFolder(currentDirectory);
  }

  // If it's the first time executing obtaining the data from the datasource and
  // assign the resolution
  if (fold == 0 && batch == 0) {
    originalData = dataSource.getAllSamples()->getData();
    resolution = static_cast<size_t>(pow(2,
      model.getFitterConfiguration().getGridConfig().level_+2));
  }

  size_t nDimensions = model.getDataset()->getDimension();
  omp_set_num_threads(static_cast<int> (config.getVisualizationParameters().numberCores));

  #pragma omp parallel sections
  {
    #pragma omp section
    {
      if (std::find(config.getGeneralConfig().algorithm.begin(),
          config.getGeneralConfig().algorithm.end(), "linearcuts") !=
          config.getGeneralConfig().algorithm.end()) {
        DataMatrix cutMatrix;
        getLinearCuts(model, currentDirectory, cutMatrix, nDimensions);
      }
    }

    #pragma omp section
    {
      if (std::find(config.getGeneralConfig().algorithm.begin(),
                   config.getGeneralConfig().algorithm.end(), "heatmaps") !=
                   config.getGeneralConfig().algorithm.end()) {
        DataMatrix heatMapMatrix;
        getHeatmap(model, currentDirectory, heatMapMatrix, nDimensions);
      }
    }

    /*#pragma omp section
    {
      if (config.getGeneralConfig().algorithm == "tsne") {
        // Just run tsne if it's the first time the visualization module is executed
        if (fold == 0 && batch == 0) {
          runTsne(model);
        }
        // Evaluate the model on the original data and append the vector to the
        // compressed data.
        if (originalData.getNcols() > 1) {
          DataVector evaluation(originalData.getNrows());
          model.evaluate(originalData, evaluation);
          tsneCompressedData.setColumn(tsneCompressedData.getNcols()-1, evaluation);
          if ( config.getGeneralConfig().targetFileType == VisualizationFileType::CSV ) {
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
    }*/

    #pragma omp section
    {
      // Just store the grid if the output is CSV.
      if (config.getGeneralConfig().targetFileType == VisualizationFileType::CSV) {
        storeGrid(model, currentDirectory);
      }
    }
  }
}

void VisualizerDensityEstimation::storeGrid(ModelFittingBase &model,
  std::string currentDirectory) {
  ModelFittingBaseSingleGrid* gridModel = dynamic_cast<ModelFittingBaseSingleGrid*>(&model);

  auto grid = gridModel->getGrid().clone();

  DataMatrix gridMatrix;

  grid->getStorage().getCoordinateArrays(gridMatrix);

  CSVTools::writeMatrixToCSVFile(currentDirectory + "/grid", gridMatrix);
}

/*void VisualizerDensityEstimation::runTsne(ModelFittingBase &model) {
    if ( originalData.getNcols() == 1 ) {
      std::cout << "The tsne algorithm can only be applied if "
      "the dimension is greater than 1" << std::endl;
      return;
    }

    size_t N = originalData.getNrows();
    size_t D = originalData.getNcols();

    std::cout << "Rows: "<< std::to_string(N) << std::endl;
    std::cout << "Columns: " << std::to_string(D) << std::endl;

    // For the weird case in which no data can be found TSNE won't run
    if (N == 0 || D == 0) {
      std::cout << "No data found. TSNE wo't run" << std::endl;
      return;
    }
    std::unique_ptr<double[]> input (new double[N*D]);

    std::copy(originalData.data(), originalData.data()+N*D,
      input.get());

    std::unique_ptr<double[]> output(new double[N* config.getVisualizationParameters().
                                 targetDimension]);

    if ( D > config.getVisualizationParameters().targetDimension ) {
      std::cout << "Compressing with tsne to " <<
      std::to_string(config.getVisualizationParameters().targetDimension)
      << " dimensions" << std::endl;

      TSNE tsne;
      tsne.run(input, N, D , output, config.getVisualizationParameters().targetDimension,
      config.getVisualizationParameters().perplexity, config.getVisualizationParameters().theta,
      config.getVisualizationParameters().seed, false,
      config.getVisualizationParameters().maxNumberIterations);
      D = config.getVisualizationParameters().targetDimension;
    } else {
     std::copy(output.get(), output.get()+N*D,
           input.get());
    }
    DataVector evaluation(originalData.getNrows());
    model.evaluate(originalData, evaluation);
    tsneCompressedData = DataMatrix(output.get(), N, D);
    tsneCompressedData.appendCol(evaluation);

}*/




void VisualizerDensityEstimation::storeTsneJson(DataMatrix &matrix, ModelFittingBase &model,
  std::string currentDirectory) {
  json::JSON jsonOutput;

  jsonOutput.addListAttr("data");

  jsonOutput["data"].addDictValue();
  jsonOutput["data"][0].addIDAttr("type", "\"scatter\"");
  jsonOutput["data"][0].addIDAttr("mode", "\"markers\"");

  DataVector xCol(matrix.getNrows());

  matrix.getColumn(0, xCol);

  jsonOutput["data"][0].addIDAttr("x", xCol.toString());

  DataVector yCol(matrix.getNrows());

  matrix.getColumn(1, yCol);
  jsonOutput["data"][0].addIDAttr("y", yCol.toString());

  jsonOutput["data"][0].addDictAttr("marker");

  DataVector zCol(matrix.getNrows());

  matrix.getColumn(2, zCol);

  jsonOutput["data"][0]["marker"].addIDAttr("color", zCol.toString());

  jsonOutput["data"][0]["marker"].addIDAttr("colorscale", "\"Viridis\"");

  jsonOutput["data"][0]["marker"].addIDAttr("opacity", 0.8);

  jsonOutput["data"][0]["marker"].addIDAttr("showscale", true);

  jsonOutput["data"][0]["marker"].addDictAttr("colorbar");

  jsonOutput["data"][0]["marker"]["colorbar"].addDictAttr("title");
  jsonOutput["data"][0]["marker"]["colorbar"]["title"].addIDAttr("text", "\"Density value\"");

  jsonOutput.addDictAttr("layout");

  jsonOutput["layout"].addDictAttr("title");

  jsonOutput["layout"]["title"].addIDAttr("text", "\"TSNE Compression\"");
  jsonOutput["layout"]["title"].addIDAttr("x", 0.5);

  jsonOutput.serialize(currentDirectory + "/tsneCompression.json");
}

void VisualizerDensityEstimation::storeCutJson(DataMatrix &matrix, std::vector<size_t> indexes,
size_t &varDim, std::string filepath) {
  indexes.erase(std::find(indexes.begin(), indexes.end(), varDim));


  json::JSON jsonOutput;
  jsonOutput.addListAttr("data");

  jsonOutput.addDictAttr("layout");

  jsonOutput["layout"].addDictAttr("title");

  jsonOutput["layout"].addListAttr("annotations");

  jsonOutput["layout"]["title"].addIDAttr("text", "\"Density Estimation: Linear Cuts "
                                                  "Variable dimension " +
  std::to_string(varDim + 1) + "\"");
  jsonOutput["layout"]["title"].addIDAttr("x", 0.5);

  size_t totalGraphs = ((matrix.getNcols() == 3)?5:25);


  size_t rowsPerGraph = matrix.getNrows()/totalGraphs;

  for (size_t graphNumber=0; graphNumber < totalGraphs; graphNumber++) {
    DataMatrix temp(matrix.getNrows(), matrix.getNcols());

    temp.copyFrom(matrix);

    size_t beginIndex = graphNumber*rowsPerGraph+1;

    temp.resizeToSubMatrix(beginIndex, 1, beginIndex + rowsPerGraph-1, matrix.getNcols());

    // Adding data trace
    jsonOutput["data"].addDictValue();

    jsonOutput["data"][graphNumber].addIDAttr("type", "\"scatter\"");
    jsonOutput["data"][graphNumber].addIDAttr("mode", "\"lines\"");

    std::string xAxis("x" + ((graphNumber == 0)?"":std::to_string(graphNumber+1)));
    jsonOutput["data"][graphNumber].addIDAttr("xaxis",
      "\"x" + ((graphNumber == 0)?"":std::to_string(graphNumber+1))+"\"");

    std::string yAxis("y"+((graphNumber == 0)?"":std::to_string(graphNumber+1)));
    jsonOutput["data"][graphNumber].addIDAttr("yaxis",
      "\"y" + ((graphNumber == 0)?"":std::to_string(graphNumber + 1)) + "\"");

    DataVector xCol(temp.getNrows());

    temp.getColumn(varDim, xCol);

    jsonOutput["data"][graphNumber].addIDAttr("x", xCol.toString());

    DataVector yCol(temp.getNrows());

    temp.getColumn(temp.getNcols()-1, yCol);
    jsonOutput["data"][graphNumber].addIDAttr("y", yCol.toString());
    jsonOutput["data"][graphNumber].addIDAttr("showlegend", false);
    jsonOutput["data"][graphNumber].addIDAttr("hoverinfo", "\"x+y\"");

    // Layout part of the graph
    std::string xAxisName("xaxis" + ((graphNumber == 0)?"":std::to_string(graphNumber+1)));
    std::string yAxisName("yaxis" + ((graphNumber == 0)?"":std::to_string(graphNumber+1)));


    jsonOutput["layout"].addDictAttr(xAxisName);
    jsonOutput["layout"][xAxisName].addIDAttr("anchor", "\"" + yAxis + "\"");
    jsonOutput["layout"][xAxisName].addIDAttr("type", "\"linear\"");
    jsonOutput["layout"][xAxisName].addListAttr("domain");
    jsonOutput["layout"][xAxisName]["domain"].addIdValue(0.2*
      static_cast<double>(graphNumber%5));
    jsonOutput["layout"][xAxisName]["domain"].addIdValue(0.2*
      (static_cast<double>(graphNumber%5)+1)-0.05);

    jsonOutput["layout"].addDictAttr(yAxisName);
    jsonOutput["layout"][yAxisName].addIDAttr("anchor", "\""+xAxis+"\"");
    jsonOutput["layout"][yAxisName].addIDAttr("type", "\"linear\"");

    if (matrix.getNcols() > 3) {
     jsonOutput["layout"][yAxisName].addListAttr("domain");
     jsonOutput["layout"][yAxisName]["domain"].addIdValue(0.2*
       static_cast<double>(graphNumber/5));
     jsonOutput["layout"][yAxisName]["domain"].addIdValue(0.2*
       (static_cast<double>(graphNumber/5)+1)-0.1);
    }

    // Adding titles to subplots
    DataVector firstRow(temp.getNcols());
    temp.getRow(0, firstRow);

    // Adding titles to subplots
    std::string dim1Text(std::to_string(indexes.at(0) + 1));

    std::string dim1ValueText(std::to_string(firstRow.get(indexes.at(0))));
    dim1ValueText.erase(dim1ValueText.find_last_not_of('0') + 2, std::string::npos);
    dim1Text = "\"Dim " + dim1Text + "=" + dim1ValueText;

    std::string dim2Text = "";
    std::string dim2ValueText = "\"";

    if (matrix.getNcols() > 3) {
      dim2Text = std::to_string(indexes.at(1) + 1);
      dim2ValueText = std::to_string(firstRow.get(indexes.at(1)));
      dim2ValueText.erase(dim2ValueText.find_last_not_of('0') + 2, std::string::npos);
      dim2Text = ", Dim " + dim2Text + "= "+dim2ValueText + "\"";
    } else {
      dim2Text = dim2Text + dim2ValueText;
    }

    std::string subplot_title(dim1Text + dim2Text);

    jsonOutput["layout"]["annotations"].addDictValue();
    jsonOutput["layout"]["annotations"][graphNumber].addIDAttr("x",
    (std::stod(jsonOutput["layout"][xAxisName]["domain"][0].get()) +
    std::stod(jsonOutput["layout"][xAxisName]["domain"][1].get()))/2);

    if (matrix.getNcols() > 3) {
      jsonOutput["layout"]["annotations"][graphNumber].addIDAttr("y",
      (0.9-std::stod(jsonOutput["layout"][yAxisName]["domain"][0].get())));
    } else {
      jsonOutput["layout"]["annotations"][graphNumber].addIDAttr("y", 1.0);
    }

    jsonOutput["layout"]["annotations"][graphNumber].addIDAttr("showarrow", false);
    jsonOutput["layout"]["annotations"][graphNumber].addIDAttr("xanchor", "\"center\"");
    jsonOutput["layout"]["annotations"][graphNumber].addIDAttr("yanchor", "\"bottom\"");
    jsonOutput["layout"]["annotations"][graphNumber].addIDAttr("xref", "\"paper\"");
    jsonOutput["layout"]["annotations"][graphNumber].addIDAttr("yref", "\"paper\"");
    jsonOutput["layout"]["annotations"][graphNumber].addIDAttr("text", subplot_title);
  }
  std::cout << "Writing file " << filepath+".json" << std::endl;
  jsonOutput.serialize(filepath + ".json");
}

void VisualizerDensityEstimation::storeCutJson(DataMatrix &matrix, std::string filepath) {
  json::JSON jsonOutput;
  jsonOutput.addListAttr("data");

  jsonOutput.addDictAttr("layout");

  jsonOutput["layout"].addDictAttr("title");

  jsonOutput["layout"]["title"].addIDAttr("text", "\"Density Estimation: Fitted Model \"");
  jsonOutput["layout"]["title"].addIDAttr("x", 0.5);

  // Adding data trace
  jsonOutput["data"].addDictValue();

  jsonOutput["data"][0].addIDAttr("type", "\"scatter\"");
  jsonOutput["data"][0].addIDAttr("mode", "\"lines\"");

  DataVector xCol(matrix.getNrows());

  matrix.getColumn(0, xCol);

  jsonOutput["data"][0].addIDAttr("x", xCol.toString());

  DataVector yCol(matrix.getNrows());

  matrix.getColumn(1, yCol);
  jsonOutput["data"][0].addIDAttr("y", yCol.toString());
  jsonOutput["data"][0].addIDAttr("showlegend", false);
  jsonOutput["data"][0].addIDAttr("hoverinfo", "\"x+y\"");

  std::cout << "Writing file " << filepath+".json" << std::endl;
  jsonOutput.serialize(filepath + ".json");
}

void VisualizerDensityEstimation::storeHeatmapJson(DataMatrix &matrix, ModelFittingBase &model,
std::vector<size_t> indexes, size_t &varDim1, size_t &varDim2, std::string filepath) {
  bool showLegendGroup = true;

  indexes.erase(std::find(indexes.begin(), indexes.end(), varDim1));

  indexes.erase(std::find(indexes.begin(), indexes.end(), varDim2));

  ModelFittingBaseSingleGrid* gridModel = dynamic_cast<ModelFittingBaseSingleGrid*>(&model);

  auto grid = gridModel->getGrid().clone();

  DataMatrix gridMatrix;

  grid->getStorage().getCoordinateArrays(gridMatrix);

  double maxValue = matrix.max(matrix.getNcols()-1);
  double minValue = matrix.min(matrix.getNcols()-1);

  json::JSON jsonOutput;
  jsonOutput.addListAttr("data");
  jsonOutput.addDictAttr("layout");
  jsonOutput["layout"].addDictAttr("title");

  if (gridMatrix.getNcols() >= 4) {
    jsonOutput["layout"].addIDAttr("height", "1500");
  }
  jsonOutput["layout"].addListAttr("annotations");
  jsonOutput["layout"]["title"].addIDAttr("text", "\"Density Estimation: Heatmaps "
                                                  "Variable dimensions " +
  std::to_string(varDim1 + 1) + " and " + std::to_string(varDim2 + 1) + "\"");
  jsonOutput["layout"]["title"].addIDAttr("x", 0.5);

  unsigned int graphIndex = 0;

  size_t totalGraphs = ((matrix.getNcols() == 4)?5:25);

  size_t rowsPerGraph = matrix.getNrows()/totalGraphs;

  for (size_t graphNumber = 0; graphNumber < totalGraphs; graphNumber++) {
    DataMatrix temp(matrix.getNrows(), matrix.getNcols());

    temp.copyFrom(matrix);

    size_t beginIndex = graphNumber*rowsPerGraph+1;

    temp.resizeToSubMatrix(beginIndex, 1, beginIndex+rowsPerGraph-1, matrix.getNcols());

    // Adding data trace
    jsonOutput["data"].addDictValue();

    jsonOutput["data"][graphIndex].addIDAttr("type", "\"contour\"");

    std::string xAxis("x" + ((graphNumber == 0)?"":std::to_string(graphNumber+1)));
    jsonOutput["data"][graphIndex].addIDAttr("xaxis",
    "\"x"+((graphNumber == 0)?"":std::to_string(graphNumber+1))+"\"");

    std::string yAxis("y" + ((graphNumber == 0)?"":std::to_string(graphNumber+1)));
    jsonOutput["data"][graphIndex].addIDAttr("yaxis",
    "\"y"+((graphNumber == 0)?"":std::to_string(graphNumber+1))+"\"");

    DataVector xCol(temp.getNrows());

    temp.getColumn(varDim1, xCol);

    jsonOutput["data"][graphIndex].addIDAttr("x", xCol.toString());

    DataVector yCol(temp.getNrows());

    temp.getColumn(varDim2, yCol);

    jsonOutput["data"][graphIndex].addIDAttr("y", yCol.toString());

    DataVector zCol(temp.getNrows());

    temp.getColumn(temp.getNcols()-1, zCol);

    jsonOutput["data"][graphIndex].addIDAttr("z", zCol.toString());
    jsonOutput["data"][graphIndex].addIDAttr("showlegend", false);
    jsonOutput["data"][graphIndex].addIDAttr("hoverinfo", "\"x+y+z\"");

    jsonOutput["data"][graphIndex].addIDAttr("zmin", minValue);
    jsonOutput["data"][graphIndex].addIDAttr("zmax", maxValue);
    jsonOutput["data"][graphIndex].addIDAttr("colorscale", "\"Viridis\"");

    if (graphIndex != 0) {
      jsonOutput["data"][graphIndex].addIDAttr("showscale", false);
    }

    graphIndex++;

    // Adding the grid
    DataVector firstRow(temp.getNcols());
    temp.getRow(0, firstRow);
    DataMatrix tempGrid(0, gridMatrix.getNcols());

    for (size_t index = 0; index < gridMatrix.getNrows(); index++) {
      DataVector row(gridMatrix.getNcols());
      gridMatrix.getRow(index, row);
      if (gridMatrix.getNcols() >= 4) {
        if ((row.get(indexes.at(0)) == firstRow.get(indexes.at(0)))&
          (row.get(indexes.at(1)) == firstRow.get(indexes.at(1)))) {
          tempGrid.appendRow(row);
        }
      } else {
        if (row.get(indexes.at(0)) == firstRow.get(indexes.at(0))) {
          tempGrid.appendRow(row);
        }
      }
    }

    jsonOutput["data"].addDictValue();
    jsonOutput["data"][graphIndex].addIDAttr("type", "\"scatter\"");
    jsonOutput["data"][graphIndex].addIDAttr("mode", "\"markers\"");
    jsonOutput["data"][graphIndex].addIDAttr("xaxis",
      "\"x"+((graphNumber == 0)?"":std::to_string(graphNumber+1))+"\"");
    jsonOutput["data"][graphIndex].addIDAttr("yaxis",
      "\"y"+((graphNumber == 0)?"":std::to_string(graphNumber+1))+"\"");
    jsonOutput["data"][graphIndex].addDictAttr("marker");
    jsonOutput["data"][graphIndex]["marker"].addIDAttr("color", "\"red\"");

    DataVector xColGrid(tempGrid.getNrows());

    tempGrid.getColumn(varDim1, xColGrid);

    jsonOutput["data"][graphIndex].addIDAttr("x", xColGrid.toString());

    DataVector yColGrid(tempGrid.getNrows());

    tempGrid.getColumn(varDim2, yColGrid);

    jsonOutput["data"][graphIndex].addIDAttr("y", yColGrid.toString());
    jsonOutput["data"][graphIndex].addIDAttr("legendgroup", static_cast<size_t>(1));
    jsonOutput["data"][graphIndex].addIDAttr("hoverinfo", "\"x+y\"");

    if (showLegendGroup && tempGrid.getNrows() > 0) {
           jsonOutput["data"][graphIndex].addIDAttr("showlegend", true);
           jsonOutput["data"][graphIndex].addIDAttr("name", "\"Grid\"");
           showLegendGroup = false;
    } else {
      jsonOutput["data"][graphIndex].addIDAttr("showlegend", false);
    }

    graphIndex++;

    // Layout part of the graph
    std::string xAxisName("xaxis"+((graphNumber == 0)?"":std::to_string(graphNumber+1)));
    std::string yAxisName("yaxis"+((graphNumber == 0)?"":std::to_string(graphNumber+1)));
    jsonOutput["layout"].addDictAttr(xAxisName);
    jsonOutput["layout"][xAxisName].addIDAttr("anchor", "\"" + yAxis + "\"");
    jsonOutput["layout"][xAxisName].addIDAttr("type", "\"linear\"");
    jsonOutput["layout"][xAxisName].addListAttr("domain");
    jsonOutput["layout"][xAxisName]["domain"].addIdValue(0.2*
      static_cast<double>(graphNumber%5));
    jsonOutput["layout"][xAxisName]["domain"].addIdValue(0.2*
      (static_cast<double>(graphNumber%5)+1)-0.05);

    jsonOutput["layout"].addDictAttr(yAxisName);
    jsonOutput["layout"][yAxisName].addIDAttr("anchor", "\""+xAxis+"\"");
    jsonOutput["layout"][yAxisName].addIDAttr("type", "\"linear\"");

    if (matrix.getNcols() > 4) {
      jsonOutput["layout"][yAxisName].addListAttr("domain");
      jsonOutput["layout"][yAxisName]["domain"].addIdValue(0.2*
        static_cast<double>(graphNumber/5));
      jsonOutput["layout"][yAxisName]["domain"].addIdValue(0.2*
        (static_cast<double>(graphNumber/5)+1)-0.1);
    }
    // Adding titles to subplots
    std::string dim1Text(std::to_string(indexes.at(0)+1));

    std::string dim1ValueText(std::to_string(firstRow.get(indexes.at(0))));
    dim1ValueText.erase(dim1ValueText.find_last_not_of('0') + 2, std::string::npos);

    dim1Text = "\"Dim " +
    dim1Text + "=" + dim1ValueText;

    std::string dim2Text = "";
    std::string dim2ValueText = "\"";

    if (gridMatrix.getNcols() >= 4) {
      dim2Text = std::to_string(indexes.at(1) + 1);
      dim2ValueText = std::to_string(firstRow.get(indexes.at(1)));
      dim2ValueText.erase(dim2ValueText.find_last_not_of('0') + 2, std::string::npos);
      dim2Text = ", Dim " + dim2Text+"= " + dim2ValueText + "\"";
    } else {
      dim2Text = dim2Text + dim2ValueText;
    }

    std::string subplot_title(dim1Text + dim2Text);

    jsonOutput["layout"]["annotations"].addDictValue();
    jsonOutput["layout"]["annotations"][graphNumber].addIDAttr("x",
    (std::stod(jsonOutput["layout"][xAxisName]["domain"][0].get()) +
    std::stod(jsonOutput["layout"][xAxisName]["domain"][1].get()))/2);

    if (gridMatrix.getNcols() >= 4) {
      jsonOutput["layout"]["annotations"][graphNumber].addIDAttr("y",
      (0.9-std::stod(jsonOutput["layout"][yAxisName]["domain"][0].get())));
    } else {
      jsonOutput["layout"]["annotations"][graphNumber].addIDAttr("y", 1.0);
    }

    jsonOutput["layout"]["annotations"][graphNumber].addIDAttr("showarrow", false);
    jsonOutput["layout"]["annotations"][graphNumber].addIDAttr("xanchor", "\"center\"");
    jsonOutput["layout"]["annotations"][graphNumber].addIDAttr("yanchor", "\"bottom\"");
    jsonOutput["layout"]["annotations"][graphNumber].addIDAttr("xref", "\"paper\"");
    jsonOutput["layout"]["annotations"][graphNumber].addIDAttr("yref", "\"paper\"");
    jsonOutput["layout"]["annotations"][graphNumber].addIDAttr("text", subplot_title);
  }

  jsonOutput["data"][0].addDictAttr("colorbar");
  jsonOutput["data"][0]["colorbar"].addIDAttr("title", "\"Density Value\"");

  jsonOutput["layout"].addDictAttr("legend");
  jsonOutput["layout"]["legend"].addIDAttr("x", -0.1);
  jsonOutput["layout"]["legend"].addIDAttr("y", 1.0);

  std::cout << "Writing file " << filepath+".json" << std::endl;
  jsonOutput.serialize(filepath+".json");
}

void VisualizerDensityEstimation::storeHeatmapJson(DataMatrix &matrix, ModelFittingBase &model,
std::string filepath) {
  ModelFittingBaseSingleGrid* gridModel = dynamic_cast<ModelFittingBaseSingleGrid*>(&model);

  auto grid = gridModel->getGrid().clone();

  DataMatrix gridMatrix;

  grid->getStorage().getCoordinateArrays(gridMatrix);

  double maxValue = matrix.max();
  double minValue = matrix.min();

  json::JSON jsonOutput;
  jsonOutput.addListAttr("data");

  jsonOutput.addDictAttr("layout");

  jsonOutput["layout"].addDictAttr("title");

  jsonOutput["layout"].addListAttr("annotations");

  jsonOutput["layout"]["title"].addIDAttr("text", "\"Density Estimation: 2D Fitted Model\"");
  jsonOutput["layout"]["title"].addIDAttr("x", 0.5);

  // Adding data trace
  jsonOutput["data"].addDictValue();

  jsonOutput["data"][0].addIDAttr("type", "\"contour\"");

  DataVector xCol(matrix.getNrows());

  matrix.getColumn(0, xCol);

  jsonOutput["data"][0].addIDAttr("x", xCol.toString());

  DataVector yCol(matrix.getNrows());

  matrix.getColumn(1, yCol);

  jsonOutput["data"][0].addIDAttr("y", yCol.toString());

  DataVector zCol(matrix.getNrows());

  matrix.getColumn(2, zCol);

  jsonOutput["data"][0].addIDAttr("z", zCol.toString());
  jsonOutput["data"][0].addIDAttr("showlegend", false);
  jsonOutput["data"][0].addIDAttr("hoverinfo", "\"x+y+z\"");
  jsonOutput["data"][0].addIDAttr("zmin", minValue);
  jsonOutput["data"][0].addIDAttr("zmax", maxValue);
  jsonOutput["data"][0].addIDAttr("colorscale", "\"Viridis\"");
  jsonOutput["data"][0].addDictAttr("colorbar");
  jsonOutput["data"][0]["colorbar"].addIDAttr("title", "\"Density Value\"");

  // Adding the grid
  jsonOutput["data"].addDictValue();

  jsonOutput["data"][1].addIDAttr("type", "\"scatter\"");

  jsonOutput["data"][1].addIDAttr("mode", "\"markers\"");

  jsonOutput["data"][1].addDictAttr("marker");

  jsonOutput["data"][1]["marker"].addIDAttr("color", "\"red\"");

  DataVector xColGrid(gridMatrix.getNrows());

  gridMatrix.getColumn(0, xColGrid);

  jsonOutput["data"][1].addIDAttr("x", xColGrid.toString());

  DataVector yColGrid(gridMatrix.getNrows());

  gridMatrix.getColumn(1, yColGrid);

  jsonOutput["data"][1].addIDAttr("y", yColGrid.toString());
  jsonOutput["data"][1].addIDAttr("showlegend", true);
  jsonOutput["data"][1].addIDAttr("name", "\"Grid\"");
  jsonOutput["data"][1].addIDAttr("hoverinfo", "\"x+y\"");

  jsonOutput["layout"].addDictAttr("legend");
  jsonOutput["layout"]["legend"].addIDAttr("x", -0.15);
  jsonOutput["layout"]["legend"].addIDAttr("y", 1.0);


  std::cout << "Writing file " << filepath + ".json" << std::endl;
  jsonOutput.serialize(filepath + ".json");
}

}  // namespace datadriven
}  // namespace sgpp
