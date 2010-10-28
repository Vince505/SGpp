/*
 * TruncatedTrapezoidGrid.cpp
 *
 *  Created on: Aug 4, 2010
 *      Author: Aliz Nagy
 */

/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "grid/Grid.hpp"
#include "grid/type/TruncatedTrapezoidGrid.hpp"

#include "grid/generation/TruncatedTrapezoidGridGenerator.hpp"
#include "basis/linear/boundary/operation/common/OperationEvalLinearBoundary.hpp"
#include "basis/linear/boundary/operation/common/OperationHierarchisationLinearBoundary.hpp"
#include "sgpp.hpp"

#include <iostream>

namespace sg
{

TruncatedTrapezoidGrid::TruncatedTrapezoidGrid(std::istream& istr) : Grid(istr)
{

}

TruncatedTrapezoidGrid::TruncatedTrapezoidGrid(size_t dim)
{
	this->storage = new GridStorage(dim);
}

TruncatedTrapezoidGrid::TruncatedTrapezoidGrid(BoundingBox& BB)
{
	this->storage = new GridStorage(BB);
}

TruncatedTrapezoidGrid::~TruncatedTrapezoidGrid()
{
}

const char* TruncatedTrapezoidGrid::getType()
{
	return "TruncatedTrapezoid";
}
Grid* TruncatedTrapezoidGrid::unserialize(std::istream& istr)
{
	return new TruncatedTrapezoidGrid(istr);
}
/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator* TruncatedTrapezoidGrid::createGridGenerator()
{
	return new TruncatedTrapezoidGridGenerator(this->storage);
}
OperationHierarchisation* TruncatedTrapezoidGrid::createOperationHierarchisation()
{
	return new OperationHierarchisationLinearBoundary(this->storage);
}
OperationEval* TruncatedTrapezoidGrid::createOperationEval()
{
	return new OperationEvalLinearBoundary(this->storage);
}

OperationConvert* TruncatedTrapezoidGrid::createOperationConvert()
{
	throw factory_exception("Unsupported operation");
}

}

