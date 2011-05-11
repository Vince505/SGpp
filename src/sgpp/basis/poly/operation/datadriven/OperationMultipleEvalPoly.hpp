/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONMULTIPLEEVALPOLY_HPP
#define OPERATIONMULTIPLEEVALPOLY_HPP

#include "operation/datadriven/OperationMultipleEval.hpp"
#include "grid/GridStorage.hpp"

#include "sgpp.hpp"

namespace sg
{
namespace base
{

/**
 * This class implements OperationMultipleEval for a grids with poly basis ansatzfunctions
 */
class OperationMultipleEvalPoly : public OperationMultipleEval
{
public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's GridStorage object
	 * @param degree the polynom's max. degree
	 */
	OperationMultipleEvalPoly(GridStorage* storage, size_t degree, DataMatrix* dataset) : OperationMultipleEval(dataset), base(degree) {
		this->storage = storage;
	}

	/**
	 * Destructor
	 */
	virtual ~OperationMultipleEvalPoly() {}

	virtual void mult(DataVector& alpha, DataVector& result);
	virtual void multTranspose(DataVector& source, DataVector& result);

protected:
	/// Pointer to GridStorage object
	GridStorage* storage;
	/// Poly Basis object
	SPolyBase base;
};

}
}

#endif /* OPERATIONMULTIPLEEVALPOLY_HPP */