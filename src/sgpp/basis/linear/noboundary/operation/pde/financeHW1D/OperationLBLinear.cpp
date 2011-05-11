/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Chao qi(qic@in.tum.de)

#include "basis/linear/noboundary/operation/pde/financeHW1D/OperationLBLinear.hpp"

#include "basis/linear/noboundary/algorithm_sweep/DPhiPhiDownBBLinear.hpp"
#include "basis/linear/noboundary/algorithm_sweep/DPhiPhiUpBBLinear.hpp"

#include "algorithm/common/sweep.hpp"
using namespace sg::pde;
using namespace sg::base;

namespace sg
{
namespace finance
{

OperationLBLinear::OperationLBLinear(GridStorage* storage) : StdUpDown(storage)
{
}

OperationLBLinear::~OperationLBLinear()
{
}

void OperationLBLinear::up(DataVector& alpha, DataVector& result, size_t dim)
{
	// Dphi * phi
	DPhiPhiUpBBLinear func(this->storage);
	sweep<DPhiPhiUpBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationLBLinear::down(DataVector& alpha, DataVector& result, size_t dim)
{
	// Dphi * phi
	DPhiPhiDownBBLinear func(this->storage);
	sweep<DPhiPhiDownBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

}
}