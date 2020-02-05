/*
 *
 * Copyright (c) 2014, Laurens van der Maaten (Delft University of Technology)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. All advertising materials mentioning features or use of this software
 *    must display the following acknowledgement:
 *    This product includes software developed by the Delft University of Technology.
 * 4. Neither the name of the Delft University of Technology nor the names of
 *    its contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY LAURENS VAN DER MAATEN ''AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
 * EVENT SHALL LAURENS VAN DER MAATEN BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
 * IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
 * OF SUCH DAMAGE.
 *
 */

/**
 * Code originally taken from https://lvdmaaten.github.io/tsne/
 * It has been modified in order to be adapted to the SG++ datamining
 * pipeline structure and has been parallelized
 */

#include <memory>

namespace sgpp {
namespace datadriven {
/**
 * Returns the mathematical sign of a given number
 * @param x Number to obtain the sign from
 * @return -1 if x is negative, 1 if positive and 0 if 0
 */
static inline double sign(double x) { return (x == .0 ? .0 : (x < .0 ? -1.0 : 1.0)); }
class TSNE {
 public:
  /**
   * Default constructor
   */
  TSNE();
  /**
   * Runs the Barnes-Huts version of the t-SNE algorithm
   * @param X Pointer to a double array with the input data
   * @param N Number of rows of the input data
   * @param D Dimension of the input data
   * @param Y Pointer to a double array where the output will be stored
   * @param no_dims Dimension to which the data will be compressed
   * @param perplexity Perplexity parameter
   * @param theta Theta parameter for gradient approximation
   * @param rand_seed Seed for random intitialization
   * @param skip_random_init bool to determine if we skip the random initizialization
   * @param max_iter Maximum number of iterations
   * @param mom_switch_iter Iteration in which to change the momentum value of the gradient descent
   */
  void run(std::unique_ptr<double[]> &X, size_t N, size_t D,
    std::unique_ptr<double[]> &Y, size_t no_dims, double perplexity,
    double theta, size_t rand_seed,
    bool skip_random_init, size_t max_iter, size_t mom_switch_iter =250);

 private:
  /**
   * Functions that computes the approximate gradient of the Kullback-leiber divergence
   * @param inp_row_P Array with pointers to the array inp_col_P indicating,
   * where the probability values for a given row lie in the inp_val_P
   * @param inp_col_P Array with pointers to the array inp_val_P
   * @param inp_val_P Array with the input probability values
   * @param Y Pointer to a double array where the result of the previous interation was stored
   * @param N Number of total input points
   * @param D Dimension of the input points
   * @param dC Array storing the approximate value of the gradient
   * @param theta Theta parameter for gradient approximation
   */
  void computeGradient(std::unique_ptr<size_t[]>  &inp_row_P,
    std::unique_ptr<size_t[]>  &inp_col_P,
    std::unique_ptr<double[]>  &inp_val_P, std::unique_ptr<double[]>  &Y, size_t N, size_t D,
    std::unique_ptr<double[]>  &dC, double theta);

  /**
   * Evaluates the current value of the Kullback-leiber divergence
   * @param row_P Array with pointers to the array inp_col_P indicating,
   * where the probability values for a given row lie in the inp_val_P
   * @param col_P Array with pointers to the array val_P
   * @param val_P Array with the input probability values
   * @param Y Pointer to a double array where the result of the previous interation was stored
   * @param N Number of total input points
   * @param D Dimension of the input points
   * @param theta Theta parameter for gradient approximation
   * @return Current value of the Kullback-leiber divergence. This the function we are trying
   * to optimize
   */
  double evaluateError(std::unique_ptr<size_t[]>  &row_P,
    std::unique_ptr<size_t[]>  &col_P,
    std::unique_ptr<double[]>  &val_P,
    std::unique_ptr<double[]>  &Y, size_t N, size_t D, double theta);

  /**
   * Translate all of the data points to a zero mean distribution
   * @param X Double array containing all of the coordinates of the input points
   * @param N Total number of input points
   * @param D Dimensionality of the input
   */
  void zeroMean(double* X, size_t N, size_t D);

  /**
   * Method that calculates the probability that a point in a high dimensional space will choose
   * another one as its neighbor if we assume a distribution is a gaussian centered in the former.
   * The variance of the gaussian is given by the perplexity value and the probabilities are only
   * calculated for the K nearest neighbors of any given point. Results are stored in a sparse-matrix
   * representation since most of the probabilities will be 0.
   * @param X Pointer to a double array with the input data
   * @param N Number of total input points
   * @param D Dimension of the input points
   * @param row_P Array with pointers to the array inp_col_P indicating,
   * where the probability values for a given row lie in the inp_val_P
   * @param col_P Array with pointers to the array val_P
   * @param val_P  Array where the probability values will be stored
   * @param perplexity Perplexity parameter which determines the variance of our
   * gaussian distribution
   * @param K The number of nearest neighbors to consider in the calculations
   */
  void computeGaussianPerplexity(std::unique_ptr<double[]>  &X,
    size_t N, size_t D, std::unique_ptr<size_t[]> &row_P,
    std::unique_ptr<size_t[]> &col_P, std::unique_ptr<double[]>  &val_P,
    double perplexity, size_t K);

  /**
   * Symmetrizes the sparse matrix with the probabilities
   * @param row_P Array with pointers to the array col_P indicating,
   * which columns of the latter have non zero values
   * @param col_P Array with pointers to the array val_P
   * @param val_P Array where the probability values will be stored
   * @param N Number of total input points
   */
  void symmetrizeMatrix(std::unique_ptr<size_t[]> &row_P,
    std::unique_ptr<size_t[]> &col_P,
    std::unique_ptr<double[]> &val_P,
    size_t N);

  /**
   * Deliver a random numbers
   * @return Random number
   */
  double randn();
};

}  // namespace datadriven
}  // namespace sgpp
