/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de), Joerg Blank (blankj@in.tum.de), Alexander Heinecke (alexander.heinecke@mytum.de)

%module(directors="1") finance
%feature("docstring");

%include "../../../base/src/sgpp/globaldef.hpp"

%include "stl.i"
%include "std_vector.i"
%include "std_pair.i"
%include "std_complex.i"
%include "std_map.i"

%include "cpointer.i" 
%include "typemaps.i"

%include "exception.i"

%{
#define SWIG_FILE_WITH_INIT
%}

%exception {
  try {
    $action
  } catch (const std::exception& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  }
}

%include "carrays.i"
%array_class(unsigned int, unsignedIntArray);
%array_class(bool,BoolArray);
%array_class(int, IntArray);


namespace std {
	%template(IntVector) vector<int>;
	%template(IntIntVector) vector< vector<int> >; 
	%template(BoolVector) vector<bool>;
	%template(DoubleVector) vector<double>;
	%template(IndexValPair) pair<size_t, double>;
        %template(IndexValVector) vector<pair<size_t, double> >;
        // For OnlinePredictiveRefinementDimension
        %template(refinement_key) std::pair<size_t, unsigned int>;
        %template(refinement_map) std::map<std::pair<size_t, unsigned int>, double>;

}

// This should include all necessary header files
%{
#include "src/sgpp_finance.hpp"
using namespace SGPP;
%}

// The Good, i.e. without any modifications
%include "finance/src/sgpp/finance/algorithm/BlackScholesParabolicPDESolverSystem.hpp"
%include "finance/src/sgpp/finance/algorithm/BlackScholesParabolicPDESolverSystemEuroAmer.hpp"
%include "finance/src/sgpp/finance/algorithm/BlackScholesParabolicPDESolverSystemEuroAmerParallelOMP.hpp"

%include "finance/src/sgpp/finance/application/BlackScholesSolver.hpp"
%include "finance/src/sgpp/finance/application/BlackScholesSolverWithStretching.hpp"

%include "finance/src/sgpp/finance/tools/VariableDiscountFactor.hpp"

%include "finance/src/sgpp/finance/operation/FinanceOpFactory.hpp"

%apply std::string *INPUT { std::string& istr };

%apply unsigned int *OUTPUT { unsigned int& l, unsigned int& i };

%apply std::vector<std::pair<size_t, double> > *OUTPUT { std::vector<std::pair<size_t, double> >& result };
%apply std::vector<double> *INPUT { std::vector<double>& point }; 