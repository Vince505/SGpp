// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org


#include <sgpp/datadriven/datamining/modules/visualization/Visualizer.hpp>
#include <sgpp/datadriven/datamining/modules/visualization/algorithms/bhtsne/tsne.hpp>
#include <string>
#include<iostream>
#ifdef _WIN32
#include <direct.h>
#elif defined __linux__
#include <sys/stat.h>
#include <sys/types.h>
#endif

using sgpp::datadriven::TSNE;
namespace sgpp {
namespace datadriven {

Visualizer::Visualizer() {
}

const VisualizerConfiguration &Visualizer::getVisualizerConfiguration() const {
  return config;
}

void Visualizer::createFolder(std::string folder_path) {
  #ifdef _WIN32
  _mkdir(folder_path.c_str());
  #elif __linux__
  mkdir(folder_path.c_str(), S_IRWXU | S_IRWXU | S_IROTH);
  #endif
}
void Visualizer::runTsne(DataMatrix &originalData,
  DataMatrix &compressedData) {
    if (originalData.getNcols() <= 2) {
      std::cout << "The tsne algorithm can only be applied if "
      "the dimension is greater than 2" << std::endl;
      compressedData = originalData;
      return;
    }

    size_t N = originalData.getNrows();
    size_t D = originalData.getNcols();

    std::unique_ptr<double[]> input (new double[N*D]);

    std::copy(originalData.data(), originalData.data()+N*D,
      input.get());

    std::unique_ptr<double[]> output(new double[N*2]);

    std::cout << "Compressing with tsne to 2"<<" dimensions" << std::endl;

    TSNE tsne;
    tsne.run(input, N, D , output, 2,
    config.getVisualizationParameters().perplexity, config.getVisualizationParameters().theta,
    config.getVisualizationParameters().seed, false,
    config.getVisualizationParameters().maxNumberIterations);

    compressedData = DataMatrix(output.get(), N, 2);
}

void Visualizer::initializeHeatmapMatrix(DataMatrix &heatmapMatrix, size_t &nDimensions) {
  double step = 1.0/static_cast<double>(resolution);
  heatmapMatrix.resize(0, nDimensions);
  if (nDimensions >=3) {
    if (nDimensions >= 4) {
      for (double dim1 = 0; dim1 <= 1; dim1+=0.25) {
        for (double dim2 = 0; dim2 <= 1; dim2+=0.25) {
          for (double dim3 = 0; dim3 <= 1; dim3+= step) {
            for (double dim4 = 0; dim4 <= 1; dim4+= step) {
              DataVector row(heatmapMatrix.getNcols(), 0.5);
              row.set(0, dim3);
              row.set(1, dim4);
              row.set(2, dim1);
              row.set(3, dim2);
              heatmapMatrix.appendRow(row);
            }
          }
        }
      }
    } else {
      for (double dim1 = 0; dim1 <= 1; dim1+=0.25) {
        for (double dim2 = 0; dim2 <= 1; dim2 += step) {
          for (double dim3 = 0; dim3 <= 1; dim3 += step) {
            DataVector row(3, 0.5);
            row.set(0, dim3);
            row.set(1, dim2);
            row.set(2, dim1);
            heatmapMatrix.appendRow(row);
          }
        }
      }
    }
  } else {
    for (double dim1 = 0; dim1 <= 1; dim1 += step) {
      for (double dim2 = 0; dim2 <= 1; dim2 += step) {
        DataVector row(2, 0.5);
        row.set(0, dim1);
        row.set(1, dim2);
        heatmapMatrix.appendRow(row);
      }
    }
  }
}

void Visualizer::initializeCutMatrix(DataMatrix &cutMatrix, size_t &nDimensions) {
  double step = 1.0 / static_cast<double>(resolution);
  cutMatrix.resize(0, nDimensions);
  if (nDimensions >= 3) {
    for (double dim1 = 0; dim1 <= 1; dim1 += 0.25) {
      for (double dim2 = 0; dim2 <= 1; dim2 += 0.25) {
        for (double dim3 = 0; dim3 <= 1; dim3 += step) {
          DataVector row(cutMatrix.getNcols(), 0.5);
          row.set(2, dim2);
          row.set(0, dim3);
          row.set(1, dim1);
          cutMatrix.appendRow(row);
        }
      }
    }
  } else if (nDimensions == 2) {
    for (double dim1 = 0; dim1 <= 1; dim1 += 0.25) {
      for (double dim2 = 0; dim2 <= 1; dim2 += step) {
        DataVector row(2, 0.5);
        row.set(0, dim2);
        row.set(1, dim1);
        cutMatrix.appendRow(row);
      }
    }
  } else {
    for (double dim1 = 0; dim1 <= 1; dim1 += step) {
      DataVector row(1);
      row.set(0, dim1);
      cutMatrix.appendRow(row);
    }
  }
}

void Visualizer::getLinearCuts(ModelFittingBase &model,
  std::string currentDirectory, DataMatrix &cutMatrix, size_t &nDimensions) {
  initializeCutMatrix(cutMatrix, nDimensions);
  // Here it decide which method to execute based on the dimensionality of the
  // model
  if ( nDimensions >= 2 ) {
    createFolder(currentDirectory+"/LinearCuts");
    if (nDimensions >=3) {
      getLinearCutsMore3D(model, currentDirectory, cutMatrix);
    } else {
      getLinearCuts2D(model, currentDirectory, cutMatrix);
    }
  } else {
    getLinearCuts1D(model, currentDirectory, cutMatrix);
  }
}

void Visualizer::getHeatmap(ModelFittingBase &model,
  std::string currentDirectory, DataMatrix &heatMapMatrix, size_t &nDimensions) {
  std::cout << "Generating the heatmaps" << std::endl;
  // Here it decide which method to execute based on the dimensionality of the
  // model

  if ( nDimensions == 1 ) {
    std::cout << "Heatmap generation is not available for models of 1 dimension" <<std::endl;
    return;
  } else {
    initializeHeatmapMatrix(heatMapMatrix, nDimensions);
  }
  if ( nDimensions >=3 ) {
    createFolder(currentDirectory+"/Heatmaps");
    if ( nDimensions >= 4 ) {
      getHeatmapMore4D(model, currentDirectory, heatMapMatrix);
    } else if ( nDimensions == 3 ) {
      getHeatmap3D(model, currentDirectory, heatMapMatrix);
    }
  } else {
    getHeatmap2D(model, currentDirectory, heatMapMatrix);
  }
}

void Visualizer::translateColumns(DataMatrix &matrix, size_t maxColumns) {
  DataMatrix temp(matrix);

  DataVector column(matrix.getNrows());

  for (size_t dimension=0; dimension < maxColumns-1; dimension++) {
    matrix.getColumn(dimension, column);

    temp.setColumn(dimension+1, column);
  }
  matrix.getColumn(maxColumns-1, column);
  temp.setColumn(0, column);
  matrix.copyFrom(temp);
}

void Visualizer::translateColumnsRight(DataMatrix &matrix,
                                                        std::vector <size_t> indexes) {
  DataMatrix temp(matrix);

  DataVector column(matrix.getNrows());

  for (size_t dimension = 0; dimension < indexes.size()-1; dimension++) {
    matrix.getColumn(indexes.at(dimension), column);
    temp.setColumn(indexes.at(dimension+1), column);
  }

  matrix.getColumn(indexes.back(), column);
  temp.setColumn(indexes.front(), column);
  matrix.copyFrom(temp);
}

void Visualizer::translateColumnsLeft(DataMatrix &matrix,
                                                       std::vector <size_t> indexes) {
  DataMatrix temp(matrix);

  DataVector column(matrix.getNrows());

  for (size_t dimension = indexes.size()-1; dimension > 0; dimension--) {
    DataVector column(matrix.getNrows());
    matrix.getColumn(indexes.at(dimension), column);
    temp.setColumn(indexes.at(dimension-1), column);
  }

  matrix.getColumn(indexes.front(), column);
  temp.setColumn(indexes.back(), column);
  matrix.copyFrom(temp);
}

void Visualizer::updateIndexesCuts(std::vector <size_t> &indexes,
                                                    DataMatrix &matrix) {
  if (indexes.at(2) < matrix.getNcols()-1) {
    indexes.at(2)++;
    swapColumns(matrix, indexes.at(2)-1, indexes.at(2));
  } else {
    if (indexes.at(1) < matrix.getNcols()-2) {
      indexes.at(1)++;
      indexes.at(2) = indexes.at(1)+1;

      swapColumns(matrix, indexes.at(2)-1, indexes.at(2));
      swapColumns(matrix, indexes.at(1)-1, indexes.at(1));

    } else {
      indexes.at(0)++;
      indexes.at(1) = indexes.at(0)+1;
      indexes.at(2) = indexes.at(1)+1;

      swapColumns(matrix, matrix.getNcols()-1, indexes.at(2));
      swapColumns(matrix, matrix.getNcols()-2, indexes.at(1));
      swapColumns(matrix, indexes.at(0)-1, indexes.at(0));
    }
  }
}

void Visualizer::updateIndexesHeatmap(std::vector <size_t> &indexes,
                                                       DataMatrix &matrix) {
  if (indexes.at(3) < matrix.getNcols()-1) {
    indexes.at(3)++;
    swapColumns(matrix, indexes.at(3)-1, indexes.at(3));
  } else {
    if (indexes.at(2) < matrix.getNcols()-2) {
      indexes.at(2)++;
      indexes.at(3) = indexes.at(2)+1;

      swapColumns(matrix, indexes.at(2)-1, indexes.at(2));
      swapColumns(matrix, matrix.getNcols()-1, indexes.at(3));
    } else {
      if (indexes.at(1) < matrix.getNcols()-3) {
        indexes.at(1)++;
        indexes.at(2) = indexes.at(1)+1;
        indexes.at(3) = indexes.at(2)+1;

        swapColumns(matrix, matrix.getNcols()-1, indexes.at(3));
        swapColumns(matrix, matrix.getNcols()-2, indexes.at(2));
        swapColumns(matrix, indexes.at(1)-1, indexes.at(1));
      } else {
        indexes.at(0)++;
        indexes.at(1) = indexes.at(0)+1;
        indexes.at(2) = indexes.at(1)+1;
        indexes.at(2) = indexes.at(1)+1;
        indexes.at(3) = indexes.at(2)+1;

        swapColumns(matrix, matrix.getNcols()-1, indexes.at(3));
        swapColumns(matrix, matrix.getNcols()-2, indexes.at(2));
        swapColumns(matrix, matrix.getNcols()-3, indexes.at(1));
        swapColumns(matrix, indexes.at(0)-1, indexes.at(0));
      }
    }
  }
}

void Visualizer::swapColumns(DataMatrix &matrix, size_t col1, size_t col2) {
  DataVector temp1(matrix.getNrows());

  DataVector temp2(matrix.getNrows());

  matrix.getColumn(col1, temp1);
  matrix.getColumn(col2, temp2);
  matrix.setColumn(col2, temp1);
  matrix.setColumn(col1, temp2);
}

// The algorithm used for the cuts is this:
// 1° Evaluate the initial matrix
// 2° Shift one position to the right based on the indexes to be used
// 3° Store current evaluation
// 4° Update the indexes
// 5° Repeat until the last column index exceeds the number of dimensions
void Visualizer::getLinearCutsMore3D(
  ModelFittingBase &model, std::string currentDirectory, DataMatrix &cutMatrix) {
  std::string outputDir(currentDirectory + "/LinearCuts/");

  std::vector <size_t> variableColumnIndexes = {0, 1, 2};

  while (variableColumnIndexes.at(2) < cutMatrix.getNcols()) {
    std::string subfolder(outputDir+"dimensions_"+
                          std::to_string(variableColumnIndexes.at(0)+1)+"_"+
                          std::to_string(variableColumnIndexes.at(1)+1)+"_"+
                          std::to_string(variableColumnIndexes.at(2)+1));

    createFolder(subfolder);

    for (size_t combination = 0; combination < 3; combination++) {
      DataMatrix cutResults(cutMatrix);
      DataVector evaluation(cutMatrix.getNrows());

      model.evaluate(cutMatrix, evaluation);

      cutResults.appendCol(evaluation);

      translateColumnsRight(cutMatrix, variableColumnIndexes);
      if (config.getGeneralConfig().targetFileType == VisualizationFileType::CSV) {
        CSVTools::writeMatrixToCSVFile(subfolder + "/Cut_var_dimension_" +
                                       std::to_string(variableColumnIndexes.at(combination)+1),
                                       cutResults);
      } else if (config.getGeneralConfig().targetFileType == VisualizationFileType::json) {
        storeCutJson(cutResults, variableColumnIndexes, variableColumnIndexes.at(combination),
                     subfolder + "/Cut_var_dimension_"+
                     std::to_string(variableColumnIndexes.at(combination)+1));
      }
    }
    updateIndexesCuts(variableColumnIndexes, cutMatrix);
  }
}

// The algorithm used for the cuts is this:
// 1° Evaluate the initial matrix
// 2° Shift one position to the right
// 3° Store current evaluation
// 4° Repeat until we have made D shifts
void Visualizer::getLinearCuts2D(
  ModelFittingBase &model, std::string currentDirectory, DataMatrix &cutMatrix) {
  std::string outputDir(currentDirectory + "/LinearCuts/");

  std::vector <size_t> variableColumnIndexes = {0, 1, 2};

  for (size_t combination = 0; combination < 2; combination++) {
    DataMatrix cutResults(cutMatrix);
    DataVector evaluation(cutMatrix.getNrows());

    model.evaluate(cutMatrix, evaluation);

    cutResults.appendCol(evaluation);

    translateColumns(cutMatrix, cutMatrix.getNcols());

    if (config.getGeneralConfig().targetFileType == VisualizationFileType::CSV) {
      CSVTools::writeMatrixToCSVFile(outputDir + "Cut_dimensions_1_2_variable_dimension" +
                                     std::to_string(combination+1),
                                     cutResults);
    } else if (config.getGeneralConfig().targetFileType == VisualizationFileType::json) {
      storeCutJson(cutResults, variableColumnIndexes, variableColumnIndexes.at(combination),
                   outputDir + "Cut_dimensions_1_2_variable_dimension" +
                   std::to_string(combination+1));
    }
  }
}

void Visualizer::getLinearCuts1D(ModelFittingBase &model,
                                                  std::string currentDirectory, DataMatrix &cutMatrix) {
  std::string outputDir(currentDirectory+"/");

  DataMatrix cutResults(cutMatrix);
  DataVector evaluation(cutMatrix.getNrows());

  model.evaluate(cutMatrix, evaluation);

  cutResults.appendCol(evaluation);
  if (config.getGeneralConfig().targetFileType == VisualizationFileType::CSV) {
    CSVTools::writeMatrixToCSVFile(outputDir+"1D_DensityEstimationLinearCut", cutResults);
  } else if (config.getGeneralConfig().targetFileType == VisualizationFileType::json) {
    storeCutJson(cutResults, outputDir+"1D_DensityEstimationLinearCut");
  }
}

// The algorithm used for the heatmaps is this:
// 1° Evaluate the initial matrix
// 2° Store current evaluation
// 3° Shift to the right based on the column indexes
// 4° Repeat steps 1 to 3 two more times
// 5° Repeat the steps 1 to 4 but shifting to the left instead
// 6° Shift to the right if its the first time you reach this step
// 7° Shift to the left
// 8° Update the indexes
// 9° Repeat until the last column index exceeds the number of dimensions
void Visualizer::getHeatmapMore4D(
  ModelFittingBase &model, std::string currentDirectory, DataMatrix &heatMapMatrix) {
  std::string outputDir(currentDirectory+"/Heatmaps/");

  std::vector <size_t> variableColumnIndexes = {0, 1, 2, 3};

  std::vector<size_t> workingIndexes = {0, 0, 0};

  while (variableColumnIndexes.at(3) < heatMapMatrix.getNcols()) {
    std::string subfolder(outputDir+"dimensions_"+
                          std::to_string(variableColumnIndexes.at(0)+1)+"_"+
                          std::to_string(variableColumnIndexes.at(1)+1)+"_"+
                          std::to_string(variableColumnIndexes.at(2)+1)+"_"+
                          std::to_string(variableColumnIndexes.at(3)+1));
    createFolder(subfolder);
    std::copy(variableColumnIndexes.begin()+1, variableColumnIndexes.end(),
              workingIndexes.begin());

    for (size_t iteration = 0; iteration < 2; iteration++) {
      for (size_t combination = 0; combination < 3; combination++) {
        DataMatrix heatMapResults(heatMapMatrix);
        DataVector evaluation(heatMapMatrix.getNrows());
        model.evaluate(heatMapMatrix, evaluation);
        heatMapResults.appendCol(evaluation);

        if (iteration == 0) {
          if (config.getGeneralConfig().targetFileType == VisualizationFileType::CSV) {
            CSVTools::writeMatrixToCSVFile(subfolder+"/Heatmap_var_dimensions_"
                                           +std::to_string(variableColumnIndexes.at(0)+1)+"_"+
                                           std::to_string(variableColumnIndexes.at(combination+1)+1), heatMapResults);
          } else if (config.getGeneralConfig().targetFileType == VisualizationFileType::json) {
            storeHeatmapJson(heatMapResults, model,
                             variableColumnIndexes, variableColumnIndexes.at(0),
                             variableColumnIndexes.at(combination+1),
                             subfolder+"/Heatmap_var_dimensions_"
                             +std::to_string(variableColumnIndexes.at(0)+1)+"_"+
                             std::to_string(variableColumnIndexes.at(combination+1)+1));
          }
          translateColumnsRight(heatMapMatrix, workingIndexes);
        } else {
          if (config.getGeneralConfig().targetFileType == VisualizationFileType::CSV) {
            CSVTools::writeMatrixToCSVFile(subfolder+"/Heatmap_var_dimensions_"
                                           +std::to_string(variableColumnIndexes.at((combination < 2)?1:2)+1)+"_"+
                                           std::to_string(variableColumnIndexes.at((combination < 1)?2:3)+1), heatMapResults);
          } else if (config.getGeneralConfig().targetFileType == VisualizationFileType::json) {
            storeHeatmapJson(heatMapResults, model,
                             variableColumnIndexes, variableColumnIndexes.at((combination < 2)?1:2),
                             variableColumnIndexes.at((combination < 1)?2:3),
                             subfolder+"/Heatmap_var_dimensions_"
                             +std::to_string(variableColumnIndexes.at((combination < 2)?1:2)+1)+"_"+
                             std::to_string(variableColumnIndexes.at((combination < 1)?2:3)+1));
          }
          translateColumnsLeft(heatMapMatrix, workingIndexes);
        }
      }
      if (iteration == 0) {
        translateColumnsRight(heatMapMatrix, variableColumnIndexes);
      }
    }
    translateColumnsLeft(heatMapMatrix, variableColumnIndexes);

    updateIndexesHeatmap(variableColumnIndexes, heatMapMatrix);
  }
}

// The algorithm used for the heatmaps is this:
// 1° Evaluate the initial matrix
// 2° Store current evaluation
// 3° Shift to the right
// 4° Repeat steps 1 to 3 two more times
void Visualizer::getHeatmap3D(ModelFittingBase &model,
  std::string currentDirectory, DataMatrix &heatMapMatrix) {
  std::string outputDir(currentDirectory+"/Heatmaps/");

  // Dummy to reutilize the storejson method
  std::vector <size_t> variableColumnIndexes = {0, 1, 2};

  std::cout << "Resolution " << std::to_string(resolution) << std::endl;

  for (size_t combination = 0; combination < 3; combination++) {
    DataMatrix heatMapResults(heatMapMatrix);
    DataVector evaluation(heatMapMatrix.getNrows());

    model.evaluate(heatMapMatrix, evaluation);

    heatMapResults.appendCol(evaluation);

    translateColumns(heatMapMatrix, heatMapMatrix.getNcols());
    if (config.getGeneralConfig().targetFileType == VisualizationFileType::CSV) {
      CSVTools::writeMatrixToCSVFile(outputDir+"Heatmap_var_dimensions_"
       +std::to_string(combination+1)
       +"_"+((combination < 2)?std::to_string(combination+2):std::to_string(1)), heatMapResults);

    } else if (config.getGeneralConfig().targetFileType == VisualizationFileType::json) {
      storeHeatmapJson(heatMapResults,
                       model,
                       variableColumnIndexes,
                       variableColumnIndexes.at(combination),
                       variableColumnIndexes.at((combination < 2)?combination+1:0),
                       outputDir+"Heatmap_var_dimensions_"+std::to_string(combination+1)+"_"+
                       ((combination < 2)?std::to_string(combination+2):std::to_string(1)));
    }
  }
}

void Visualizer::getHeatmap2D(
  ModelFittingBase &model, std::string currentDirectory, DataMatrix &heatMapMatrix) {
  std::string outputDir(currentDirectory+"/");

  DataMatrix heatMapResults(heatMapMatrix);
  DataVector evaluation(heatMapMatrix.getNrows());

  model.evaluate(heatMapMatrix, evaluation);

  heatMapResults.appendCol(evaluation);

  if (config.getGeneralConfig().targetFileType == VisualizationFileType::CSV) {
    CSVTools::writeMatrixToCSVFile(outputDir+"2D_Heatmap", heatMapResults);
  } else if (config.getGeneralConfig().targetFileType == VisualizationFileType::json) {
    storeHeatmapJson(heatMapResults,
                     model, outputDir+"2D_Heatmap");
  }
}

void Visualizer::setResolution(size_t resolution) {
  this->resolution = resolution;
}
}  // namespace datadriven
}  // namespace sgpp
