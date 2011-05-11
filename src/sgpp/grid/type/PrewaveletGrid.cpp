/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2008 Jörg Blank (blankj@in.tum.de)                          */
/* Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
/*                                                                           */
/* sgpp is free software; you can redistribute it and/or modify              */
/* it under the terms of the GNU General Public License as published by      */
/* the Free Software Foundation; either version 3 of the License, or         */
/* (at your option) any later version.                                       */
/*                                                                           */
/* sgpp is distributed in the hope that it will be useful,                   */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU General Public License         */
/* along with sgpp; if not, write to the Free Software                       */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/

#include "grid/Grid.hpp"
#include "grid/type/PrewaveletGrid.hpp"

#include "grid/generation/PrewaveletGridGenerator.hpp"

// Include all operations on the prewavelet grid
//#include "basis/prewavelet/operation/datadriven/OperationMultipleEvalPrewavelet.hpp"
//#include "basis/prewavelet/operation/datadriven/OperationTestPrewavelet.hpp"
//#include "basis/prewavelet/operation/datadriven/OperationLaplacePrewavelet.hpp"
//#include "basis/prewavelet/operation/common/OperationEvalPrewavelet.hpp"
//#include "basis/prewavelet/operation/common/OperationHierarchisationPrewavelet.hpp"
//#include "basis/prewavelet/operation/common/OperationConvertPrewavelet.hpp"

#include "exception/factory_exception.hpp"

#include "sgpp.hpp"



namespace sg
{
namespace base
{

PrewaveletGrid::PrewaveletGrid(std::istream& istr) : Grid(istr)
{
}

PrewaveletGrid::PrewaveletGrid(size_t dim)
{
	this->storage = new GridStorage(dim);
	this->shadowStorage = new GridStorage(dim);
}

PrewaveletGrid::~PrewaveletGrid()
{
}

const char* PrewaveletGrid::getType()
{
	return "prewavelet";
}

Grid* PrewaveletGrid::unserialize(std::istream& istr)
{
	return new PrewaveletGrid(istr);
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator* PrewaveletGrid::createGridGenerator()
{
	return new PrewaveletGridGenerator(this->storage, this->shadowStorage);
}

/*OperationMultipleEval* PrewaveletGrid::createOperationMultipleEval(DataMatrix* dataset)
{
	return new OperationMultipleEvalPrewavelet(this->storage, dataset);
}*/

//OperationMultipleEvalVectorized* PrewaveletGrid::createOperationMultipleEvalVectorized(const std::string& VecType, DataMatrix* dataset)
//{
//	throw factory_exception("Unsupported operation");
//}
//
//OperationMultipleEvalVectorizedSP* PrewaveletGrid::createOperationMultipleEvalVectorizedSP(const std::string& VecType, DataMatrixSP* dataset)
//{
//	throw factory_exception("Unsupported operation");
//}

/*OperationMatrix* PrewaveletGrid::createOperationLaplace()
{
	return new OperationLaplacePrewavelet(this->storage, this->shadowStorage);
}*/

/*OperationEval* PrewaveletGrid::createOperationEval()
{
	return new OperationEvalPrewavelet(this->storage);
}*/

//OperationHierarchisation* PrewaveletGrid::createOperationHierarchisation()
//{
//	return new OperationHierarchisationPrewavelet(this->storage, this->shadowStorage);
//}
/*
OperationMatrix* PrewaveletGrid::createOperationLTwoDotProduct()
{
	throw factory_exception("Unsupported operation");
}*/

/*OperationMatrix* PrewaveletGrid::createOperationLE()
{
	throw factory_exception("Unsupported operation");
}

OperationMatrix* PrewaveletGrid::createOperationLB()
{
	throw factory_exception("Unsupported operation");
}

OperationMatrix* PrewaveletGrid::createOperationLF()
{
	throw factory_exception("Unsupported operation");
}

OperationMatrix* PrewaveletGrid::createOperationLD()
{
	throw factory_exception("Unsupported operation");
}*/

//OperationTest* PrewaveletGrid::createOperationTest()
//{
//	return new OperationTestPrewavelet(this->storage);
//}

// finance operations
/////////////////////
//OperationMatrix* PrewaveletGrid::createOperationDelta(DataVector& coef)
//{
//	throw factory_exception("Unsupported operation");
//}

/*OperationMatrix* PrewaveletGrid::createOperationGamma(DataMatrix& coef)
{
	throw factory_exception("Unsupported operation");
}*/

//OperationMatrix* PrewaveletGrid::createOperationDeltaLog(DataVector& coef)
//{
//	throw factory_exception("Unsupported operation");
//}

/*OperationMatrix* PrewaveletGrid::createOperationGammaLog(DataMatrix& coef)
{
	throw factory_exception("Unsupported operation");
}*/

//OperationConvert* PrewaveletGrid::createOperationConvert()
//{
//	return new OperationConvertPrewavelet(this->storage,this->shadowStorage);
//}

GridStorage* PrewaveletGrid::getShadowStorage()
{
	return this->shadowStorage;
}

}
}