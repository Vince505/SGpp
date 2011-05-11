/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Chao qi(qic@in.tum.de)

#include "basis/linear/noboundary/operation/pde/financeHW1D/OperationLDLinear.hpp"

#include "basis/linear/noboundary/algorithm_sweep/XPhiPhiDownBBLinear.hpp"
#include "basis/linear/noboundary/algorithm_sweep/XPhiPhiUpBBLinear.hpp"

#include "algorithm/common/sweep.hpp"
using namespace sg::pde;
using namespace sg::base;

namespace sg
{
namespace finance
{

OperationLDLinear::OperationLDLinear(GridStorage* storage) : StdUpDown(storage)
{
}

OperationLDLinear::~OperationLDLinear()
{
}

void OperationLDLinear::up(DataVector& alpha, DataVector& result, size_t dim)
{
	// X * phi * phi
	XPhiPhiUpBBLinear func(this->storage);
	sweep<XPhiPhiUpBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationLDLinear::down(DataVector& alpha, DataVector& result, size_t dim)
{
	// X * phi * phi
	XPhiPhiDownBBLinear func(this->storage);
	sweep<XPhiPhiDownBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

}
}