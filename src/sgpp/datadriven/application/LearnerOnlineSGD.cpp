#include <limits>

#include "LearnerOnlineSGD.hpp"

#include "parallel/datadriven/algorithm/DMSystemMatrixVectorizedIdentity.hpp"
#include "parallel/tools/TypesParallel.hpp"
#include "parallel/operation/ParallelOpFactory.hpp"
#include "base/grid/generation/hashmap/HashRefinementInconsistent.hpp"
#include "parallel/datadriven/basis/common/X86SimdKernelBase.hpp"

namespace sg
{

namespace datadriven
{

const sg::parallel::VectorizationType LearnerOnlineSGD::vecType_ = sg::parallel::VectorizationType::X86SIMD;


LearnerOnlineSGD::LearnerOnlineSGD(
    sg::datadriven::LearnerRegularizationType& regularization,
    const bool isRegression, const bool isVerbose) :
        Learner(regularization, isRegression, isVerbose),
        mainTrainDataset(NULL), mainClasses(NULL),
        testTrainDataset(NULL), testClasses(NULL),
        SGDCurrentIndex(0),
        numMainData(0), numMainDim(0),
        mainError(NULL),
        minibatchTrainDataset(NULL),
        minibatchClasses(NULL),
        minibatchError(NULL),
        errorPerBasisFunction(NULL),
        alphaAvg(NULL),
        currentGamma(0),
        currentRunIterations(0)

{}

void LearnerOnlineSGD::train(sg::base::DataMatrix& mainTrainDataset_,
                             sg::base::DataVector& mainClasses_,
                             sg::base::DataMatrix& testTrainDataset_,
                             sg::base::DataVector& testClasses_,
                             sg::base::RegularGridConfiguration& gridConfig,
                             sg::datadriven::LearnerOnlineSGDConfiguration& config_
                            )
{
    /*
     * Initialization
     */

    using namespace sg::base;

    if (alpha_ != NULL)
        delete alpha_;

    if (grid_ != NULL)
        delete grid_;

    if (mainError != NULL)
        delete mainError;

    if (minibatchTrainDataset != NULL)
        delete minibatchTrainDataset;

    if (minibatchClasses != NULL)
        delete minibatchClasses;

    if (minibatchError != NULL)
        delete minibatchError;

    if (errorPerBasisFunction != NULL)
        delete errorPerBasisFunction;

    if (alphaAvg != NULL)
        delete alphaAvg;

    if (isTrained_ == true)
        isTrained_ = false;

    InitializeGrid(gridConfig);
    if (grid_ == NULL)
    {
        return;
    }

    /*
     * Members
     */

    size_t chunkSize = sg::parallel::X86SimdKernelBase::getChunkDataPoints();
    size_t newSize = (mainTrainDataset_.getNrows() / chunkSize) * chunkSize;
    mainTrainDataset_.resize(newSize);
    mainClasses_.resize(newSize);

    mainTrainDataset = &mainTrainDataset_;
    DataMatrix mainTrainDatasetT(*mainTrainDataset);
    mainTrainDatasetT.transpose();

    mainClasses = &mainClasses_;

    testTrainDataset = &testTrainDataset_;
    DataMatrix testDatasetT(*testTrainDataset);
    testDatasetT.transpose();

    testClasses = &testClasses_;

    config = config_;

    numMainData = mainTrainDataset->getNrows();
    numMainDim = mainTrainDataset->getNcols();

    mainError = new DataVector(numMainData);
    mainError->setAll(0.0);

    minibatchTrainDataset = new sg::base::DataMatrix(0, numMainDim);
    minibatchTrainDataset->addSize(config.minibatchSize);
    minibatchTrainDataset->setAll(0.0);

    minibatchClasses = new DataVector(config.minibatchSize);
    minibatchClasses->setAll(0.0);

    minibatchError = new DataVector(config.minibatchSize);
    minibatchError->setAll(0.0);

    errorPerBasisFunction = new DataVector(grid_->getSize());
    errorPerBasisFunction->setAll(0.0);

    alphaAvg = new DataVector(grid_->getSize());
    alphaAvg->setAll(0.0);

    currentGamma = config.gamma;

    /*
     * File descriptors
     */

    std::fstream ferr0, ferr1, ferr2, fgrid, fcoor;
    std::string dir = config.experimentDir;
    ferr0.open((dir + std::string("/ferr0")).c_str(),
               std::ios::out | std::ios::trunc);
    ferr1.open((dir + std::string("/ferr1")).c_str(),
               std::ios::out | std::ios::trunc);
    ferr2.open((dir + std::string("/ferr2")).c_str(),
               std::ios::out | std::ios::trunc);
    fgrid.open((dir + std::string("/fgrid")).c_str(),
               std::ios::out | std::ios::trunc);
    fcoor.open((dir + std::string("/fcoor")).c_str(),
               std::ios::out | std::ios::trunc);

    /*
     * SGD Order
     */

    for (size_t i = 0; i < numMainData; i++)
    {
        SGDIndexOrder.push_back(i);
    }

    std::random_shuffle(SGDIndexOrder.begin(), SGDIndexOrder.end());
    SGDCurrentIndex = 0;

    /*
     * Refinement Functor
     */

    RefinementFunctor *functor = NULL;

    if(config.refinementType == "SURPLUS")
    {
        functor = new SurplusRefinementFunctor(alpha_,
                                               config.refinementNumPoints, 0.0);
    }

    else if(config.refinementType == "MSE")
    {
        /*functor = new SurplusRefinementFunctor(alpha_,
        		        RefineConfig.refinementNumPoints, 0.0);*/
        /*functor = new SurplusRefinementFunctor(alpha_,
                RefineConfig.refinementNumPoints, 0.0);*/
        functor = new SurplusRefinementFunctor(errorPerBasisFunction,
                                               config.refinementNumPoints, 0.0);
    }

    else if(config.refinementType == "WEIGHTED_ERROR_MINIBATCH")
    {
        functor = new WeightedErrorRefinementFunctor(alpha_, grid_,
                  config.refinementNumPoints,
                  -std::numeric_limits<double>::infinity());
        WeightedErrorRefinementFunctor* wfunctor =
            (WeightedErrorRefinementFunctor*) functor;
        wfunctor->setTrainDataset(minibatchTrainDataset);
        wfunctor->setClasses(minibatchClasses);
        wfunctor->setErrors(minibatchError);
    }

    else if(config.refinementType == "WEIGHTED_ERROR_ALL")
    {
        // FIXME: this case is not accounted for
        /*functor = new WeightedErrorRefinementFunctor(alpha_, grid_,
         RefineConfig.refinementNumPoints, 0.0);
         WeightedErrorRefinementFunctor* wfunctor =
         (WeightedErrorRefinementFunctor*) functor;

         wfunctor->setTrainDataset(mainTrainDataset);
         wfunctor->setClasses(mainClasses);*/
    }

    else if(config.refinementType == "PERSISTENT_ERROR")
    {
        functor = new PersistentErrorRefinementFunctor(alphaAvg, grid_,
                  config.refinementNumPoints,
                  -std::numeric_limits<double>::infinity());
        PersistentErrorRefinementFunctor* pfunctor =
            (PersistentErrorRefinementFunctor*) functor;
        pfunctor->setTrainDataset(minibatchTrainDataset);
        pfunctor->setClasses(minibatchClasses);
        pfunctor->setErrors(minibatchError);
    }

    else if(config.refinementType == "CLASSIFICATION")
    {

        functor = new ClassificationRefinementFunctor(alpha_, grid_,
                  config.refinementNumPoints, 0.0);
        ClassificationRefinementFunctor* cfunctor =
            (ClassificationRefinementFunctor*) functor;
        cfunctor->setTrainDataset(mainTrainDataset);
        cfunctor->setClasses(mainClasses);
    }

    else if(config.refinementType == "PREDICTIVE_REFINEMENT_DIMENSION")
    {

        functor = new PredictiveRefinementDimensionIndicator(grid_,
                  minibatchTrainDataset, minibatchError, config.refinementNumPoints);
    }

    if (functor == NULL)
    {
        throw base::application_exception("Invalid refinement type");
    }

    /*
     * Hash Refinement
     */

    if(config.hashRefinementType == "ONLINE_PREDICTIVE_REFINEMENT_DIMENSION")
    {
        if (!(config.refinementType == "PREDICTIVE_REFINEMENT_DIMENSION"))
        {
            throw base::application_exception("Online predictive refinement decorator supports only the corresponding ONLINE_PREDICTIVE_DIMENSION indicator");
        }
        HashRefinementInconsistent* hashRef = new HashRefinementInconsistent();
        online_refinement = new sg::base::OnlinePredictiveRefinementDimension(hashRef);
        online_refinement->setTrainDataset(minibatchTrainDataset);
        online_refinement->setErrors(minibatchError);
    }

    if(config.hashRefinementType == "ONLINE_PREDICTIVE_REFINEMENT_DIMENSION_OLD")
    {
        if (!(config.refinementType == "PREDICTIVE_REFINEMENT_DIMENSION"))
        {
            throw base::application_exception("Online predictive refinement decorator supports only the corresponding ONLINE_PREDICTIVE_DIMENSION indicator");
        }
        HashRefinementInconsistent* hashRef = new HashRefinementInconsistent();
        online_refinement_old = new sg::base::OnlinePredictiveRefinementDimensionOld(hashRef);
    }

    if (config.hashRefinementType == "HASH_REFINEMENT")
    {
        hash_refinement = new HashRefinement();
    }

    if (hash_refinement == NULL && online_refinement == NULL && online_refinement_old == NULL)
    {
        throw base::application_exception("Invalid hash refinement type");
    }

    /*
     * Perform runs
     */

    for (int countRun = 0; countRun < config.numRuns; countRun++)
    {

        currentRunIterations = 0;

        std::cout << "Run: " << countRun + 1 << std::endl;

        while (currentRunIterations < numMainData)
        {
            /*
             * Perform SGD
             */

            if(config.refinementCondition == "FIXED_NUMBER")
            {
                for (size_t i = 0;
                        i < config.numIterations;
                        i++)
                {
                    currentGamma = config.gamma*pow(1+config.gamma*config.lambda*(double)currentRunIterations, -2.0/3);
                    performSGDStep();
                    currentRunIterations++;
                }

            }

            if(config.refinementCondition == "SMOOTHED_ERROR_DECLINE")
            {

                double oldErrorSum = 0;
                double oldErrorLast = 0;

                double currentMinibatchError = 0;
                double ratio = 1;

                do
                {
                    // Run SGD
                    performSGDStep();
                    currentRunIterations++;

                    // Calculate smoothed error
                    currentMinibatchError = getError(minibatchTrainDataset, minibatchClasses, config.errorType, NULL, false);
                    if (smoothedErrorDeclineBuffer.size() >= config.smoothedErrorDeclineBufferSize)
                    {

                        // Calculate average of old minibatch errors
                        for (
                            std::list<double>::iterator it = smoothedErrorDeclineBuffer.begin();
                            it != smoothedErrorDeclineBuffer.end();
                            ++it
                        )
                        {
                            oldErrorSum += *it;
                        }

                        oldErrorLast = smoothedErrorDeclineBuffer.back();

                        // Update errorOnMinibatch
                        smoothedErrorDeclineBuffer.pop_back();
                        smoothedErrorDeclineBuffer.push_front(currentMinibatchError);

                        // Update ratio
                        ratio = (oldErrorLast - currentMinibatchError) / oldErrorSum;
                        ratio *= 1.0 / (double) config.smoothedErrorDeclineBufferSize;

                    }
                    else
                    {
                        smoothedErrorDeclineBuffer.push_front(currentMinibatchError);
                    }
                }
                while (ratio > config.smoothedErrorDecline);
            }

            /*
             * Error vectors
             */

            if (config.refinementType == "PERSISTENT_ERROR" ||
                    config.refinementType == "WEIGHTED_ERROR_MINIBATCH" || config.errorType == "ACCURACY"
                    || config.hashRefinementType == "ONLINE_PREDICTIVE_REFINEMENT_DIMENSION" ||
                    config.hashRefinementType == "ONLINE_PREDICTIVE_REFINEMENT_DIMENSION_OLD"
               )
            {
                double accuracy = getError(minibatchTrainDataset,
                                           minibatchClasses,
                                           config.errorType,
                                           minibatchError, false);
                if (accuracy == 1.0)
                {
                    continue;
                }
            }

            // Update main error
            getError(&mainTrainDatasetT, mainClasses, config.errorType, mainError, true);

            if (config.refinementType == "MSE")
            {
            	/*
                //getError(mainTrainDataset, mainClasses, config.errorType, mainError, false);
                getError(&mainTrainDatasetT, mainClasses, config.errorType, mainError, true);
                mainError->sqr();
                */
                OperationMultipleEval* eval = sg::op_factory::createOperationMultipleEval(*grid_, mainTrainDataset);
                eval->mult(*mainError, *errorPerBasisFunction);
                //sg::parallel::OperationMultipleEvalVectorized* eval = sg::op_factory::createOperationMultipleEvalVectorized(*grid_,
                //        vecType_, &mainTrainDatasetT);
                //eval->multTransposeVectorized(*mainError, *errorPerBasisFunction);
                delete eval;

            }

            /*
             * Output
             */

            size_t totalIterations = currentRunIterations + countRun * numMainData;

            //double err0 = getError(minibatchTrainDataset, minibatchClasses, config.errorType, NULL, false);
            double err1 = getError(&mainTrainDatasetT, mainClasses, config.errorType, NULL, true);
            double err2 = getError(&testDatasetT, testClasses, config.errorType, NULL, true);
            // double err0 = getError(minibatchTrainDataset, minibatchClasses, config.errorType, NULL, false);
            // double err1 = getError(mainTrainDataset, mainClasses, config.errorType, NULL, false);
            // double err2 = getError(testTrainDataset, testClasses, config.errorType, NULL, false);

            /*
            ferr0 << totalIterations << "," << err0 << std::endl;
            ferr1 << totalIterations << "," << err1 << std::endl;
            ferr2 << totalIterations << "," << err2 << std::endl;
            */

            size_t numGridPoints = grid_->getStorage()->size();
            //ferr0 << numGridPoints << "," << err0 << std::endl;
            ferr1 << numGridPoints << "," << err1 << std::endl;
            ferr2 << numGridPoints << "," << err2 << std::endl;

            std::string grid_str;
            grid_->serialize(grid_str);
            fgrid << grid_str << std::endl;

            double percent = (double) totalIterations / ((int) numMainData * config.numRuns);
            percent *= 100;
            if (percent > 100)
            {
                percent = 100;
            }

            if(config.errorType == "MSE")
            {
                printf("MSE: %2.10f (at %2.2f%%)\n", err1, percent);
                fflush(stdout);
            }
            if(config.errorType == "ACCURACY")
            {
                printf("Accuracy: %2.2f%% (at %2.2f%%)\n", err1 * 100, percent);
                fflush(stdout);
            }


            /*
             * Refinement
             */

            if(config.hashRefinementType == "ONLINE_PREDICTIVE_REFINEMENT_DIMENSION")
            {
                online_refinement->free_refine(grid_->getStorage(), dynamic_cast<PredictiveRefinementDimensionIndicator*>(functor));
            }

            if(config.hashRefinementType == "ONLINE_PREDICTIVE_REFINEMENT_DIMENSION_OLD")
            {
                online_refinement_old->free_refine(grid_->getStorage(), dynamic_cast<PredictiveRefinementDimensionIndicator*>(functor));
            }

            if (config.hashRefinementType == "HASH_REFINEMENT")
            {
                hash_refinement->free_refine(grid_->getStorage(), functor);
            }

            alpha_->resizeZero(grid_->getSize());
            alphaAvg->resizeZero(grid_->getSize());
        }
    }

    /*
     * Perform CG
     */

    /*
    std::cout << "Error before CG (ACCURACY): " <<
    getError(mainTrainDataset, mainClasses, "ACCURACY", NULL, false) << std::endl;
    */
    std::cout << "Error before CG (ACCURACY): " <<
    getError(&mainTrainDatasetT, mainClasses, "ACCURACY", NULL, true) << std::endl;

    /*
    std::cout << "Error before CG (MSE): " <<
    getError(mainTrainDataset, mainClasses, "MSE", NULL, false) << std::endl;
     */
    std::cout << "Error before CG (MSE): " <<
    getError(&mainTrainDatasetT, mainClasses, "MSE", NULL, true) << std::endl;

    sg::solver::ConjugateGradients *cg = new sg::solver::ConjugateGradients(
                                             config.CG_max, config.CG_eps);

    //    sg::base::OperationMatrix *C_ = sg::op_factory::createOperationIdentity(
    //                                        *this->grid_);
    //    sg::datadriven::DMSystemMatrix matrix(*grid_, *mainTrainDataset, *C_,
    //                                          config.lambda);

    sg::parallel::DMSystemMatrixVectorizedIdentity matrix(*grid_, *mainTrainDataset,
            config.lambda, vecType_);

    sg::base::DataVector b(alpha_->getSize());
    matrix.generateb(*mainClasses, b);

    cg->solve(matrix, *alpha_, b, true, false);
    *alphaAvg = *alpha_;

    /*
    std::cout << "Error after CG (ACCURACY): " <<
    getError(mainTrainDataset, mainClasses, "ACCURACY", NULL, false) << std::endl;
    std::cout << "Error after CG (MSE): " <<
    getError(mainTrainDataset, mainClasses, "MSE", NULL, false) << std::endl;
    */
    std::cout << "Error after CG (ACCURACY): " <<
    getError(&mainTrainDatasetT, mainClasses, "ACCURACY", NULL, true) << std::endl;
    std::cout << "Error after CG (MSE): " <<
    getError(&mainTrainDatasetT, mainClasses, "MSE", NULL, true) << std::endl;

    std::cout << "Error on test (ACCURACY): " <<
    getError(&testDatasetT, mainClasses, "ACCURACY", NULL, true) << std::endl;
    std::cout << "Error on test (MSE): " <<
    getError(&testDatasetT, mainClasses, "MSE", NULL, true) << std::endl;

    /*
     * Clean up
     */

    isTrained_ = true;

    ferr0.close();
    ferr1.close();
    ferr2.close();
    fgrid.close();
    fcoor.close();

    //    delete C_;
    delete cg;
    delete functor;
}

double LearnerOnlineSGD::getError(sg::base::DataMatrix* trainDataset,
                                  sg::base::DataVector* classes, std::string errorType, sg::base::DataVector* error, bool useEvalVectorized)
{
    using namespace sg::base;

    size_t numData;

    if( useEvalVectorized ) {
    	numData = trainDataset->getNcols();
    } else {
    	numData = trainDataset->getNrows();
    }


    bool cleanup = false;
    if( error == NULL )
    {
        error = new DataVector(numData);
        cleanup = true;
    }

    DataVector result(numData);

    if( useEvalVectorized )
    {
        sg::parallel::OperationMultipleEvalVectorized* eval = sg::op_factory::createOperationMultipleEvalVectorized(*grid_,
                vecType_, trainDataset);
        eval->multVectorized(*alphaAvg, result);

        delete eval;
    }
    else
    {
        OperationMultipleEval* eval = sg::op_factory::createOperationMultipleEval(*grid_, trainDataset);
        eval->mult(*alphaAvg, result);

        delete eval;
    }

    double res = -1.0;

    if(errorType == "MSE")
    {
        for (unsigned int i = 0; i < numData; i++)
        {
            error->set(i, classes->get(i) - result.get(i));
        }

        // Error
        double sum = 0;
        for (unsigned int i = 0; i < numData; i++)
        {
            sum += error->get(i) * error->get(i);
        }

        res = (sum / (double) numData);

    }
    if(errorType == "ACCURACY")
    {
        unsigned int correct = 0;
        for (unsigned int i = 0; i < numData; i++)
        {
            correct += (result.get(i) < 0) == (classes->get(i) < 0) ? 1 : 0;
            error->set(i, classes->get(i) - result.get(i));

        }
        res = static_cast<double>(correct) / static_cast<double>(numData);

    }

    if (cleanup)
    {
        delete error;
    }

    return res;
}

void LearnerOnlineSGD::performSGDStep()
{
    using namespace sg::base;

    size_t numCoeff = grid_->getStorage()->size();

    // Get x and y pair
    DataVector x(numMainDim);
    mainTrainDataset->getRow(SGDIndexOrder[SGDCurrentIndex], x);
    double y = mainClasses->get(SGDIndexOrder[SGDCurrentIndex]);

    // Store in minibatch
    pushMinibatch(x, y);


    // Update SGDCurrentIndex
    if (SGDCurrentIndex == SGDIndexOrder.size() - 1)
    {
        std::random_shuffle(SGDIndexOrder.begin(), SGDIndexOrder.end());
        SGDCurrentIndex = 0;
    }
    else
    {
        SGDCurrentIndex++;
    }

    // Calculate delta^n according to [Maier BA, 5.10]:
    // b_k^T * alpha^n - y_k
    double tmp1 = grid_->eval(*alpha_, x) - y;

    // delta^n = 2 * gamma * (b_k * tmp1 + lambda * a^n)
    DataVector delta(*alpha_);
    DataVector unit_alpha(numCoeff);
    unit_alpha.setAll(0.0);

    DataVector singleAlpha(1);
    singleAlpha[0] = 1.0;

    DataMatrix dm(x.getPointer(), 1, x.getSize());
    OperationMultipleEval *multEval = sg::op_factory::createOperationMultipleEval(*grid_, &dm);
    multEval->multTranspose(singleAlpha, delta);
    delete multEval;
    alpha_->mult(1-currentGamma * config.lambda);
    alpha_->axpy(-currentGamma * tmp1, delta);

    //double mu = 0.1; // boring exponential smoothing

    // L. Bottou exciting smoothing
    size_t t1 = (currentRunIterations > numMainDim+1) ? currentRunIterations - numMainDim : 1;
    size_t t2 = (currentRunIterations > numMainData+1) ? currentRunIterations - numMainData: 1;
    double mu = (t1>t2) ? static_cast<double>(t1) : static_cast<double>(t2);
    mu = 1.0/mu;

    alphaAvg->mult(1-mu);
    alphaAvg->axpy(mu, *alpha_);

    /*for (size_t i = 0; i < numCoeff; i++) {
    	unit_alpha[i] = 1;
    	delta[i] = grid_->eval(unit_alpha, x) * tmp1;
    	delta[i] += lambda * (*alpha_)[i];
    	delta[i] *= 2 * gamma;
    	unit_alpha[i] = 0;
}

    // update alpha
    // a^{n+1} = a^n - delta^n
    for (size_t i = 0; i < numCoeff; i++) {
    	(*alpha_)[i] = (*alpha_)[i] - delta[i];
}*/
}

void LearnerOnlineSGD::pushMinibatch(sg::base::DataVector& x, double y)
{
    static size_t next_idx = 0;
    if (minibatchTrainDataset->getUnused() > 0)
    {
        minibatchTrainDataset->appendRow(x);
        (*minibatchClasses)[next_idx] = y;
    }
    else
    {
        minibatchTrainDataset->setRow(next_idx, x);
        (*minibatchClasses)[next_idx] = y;
    }
    next_idx = (next_idx + 1) % config.minibatchSize;
}


LearnerOnlineSGD::~LearnerOnlineSGD()
{
    if (mainError != NULL)
        delete mainError;

    if (minibatchTrainDataset != NULL)
        delete minibatchTrainDataset;

    if (minibatchClasses != NULL)
        delete minibatchClasses;

    if (minibatchError != NULL)
        delete minibatchError;

    if (errorPerBasisFunction != NULL)
        delete errorPerBasisFunction;

    if (alphaAvg != NULL)
        delete alphaAvg;
}

}

}
