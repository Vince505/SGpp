/*
 * TruncatedTrapezoidGridGenerator.cpp
 *
 *  Created on: Aug 4, 2010
 *      Author: Aliz Nagy
 */
/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/


#include "base/grid/generation/TruncatedTrapezoidGridGenerator.hpp"
#include "base/grid/GridStorage.hpp"

#include "base/grid/generation/hashmap/HashGenerator.hpp"

namespace sg
{
namespace base
{

TruncatedTrapezoidGridGenerator::TruncatedTrapezoidGridGenerator(GridStorage* storage) : storage(storage)
{
}

TruncatedTrapezoidGridGenerator::~TruncatedTrapezoidGridGenerator()
{
}

void TruncatedTrapezoidGridGenerator::regular(size_t level)
{
	this->truncated( level,1);
}

void TruncatedTrapezoidGridGenerator::truncated(size_t level,size_t l_user)
{
	HashGenerator gen;
	gen.truncated(this->storage, level,l_user);
}

//void BoundaryGridGenerator::refine(RefinementFunctor* func)
//{
//	HashRefinementBoundaries refine;
//	refine.free_refine(this->storage, func);
//}
//
//size_t BoundaryGridGenerator::getNumberOfRefinablePoints()
//{
//	HashRefinementBoundaries refine;
//	return refine.getNumberOfRefinablePoints(this->storage);
//}
//
//void BoundaryGridGenerator::coarsen(CoarseningFunctor* func, DataVector* alpha)
//{
//	HashCoarsening coarsen;
//	coarsen.free_coarsen(this->storage, func, alpha);
//}
//
//size_t BoundaryGridGenerator::getNumberOfRemoveablePoints()
//{
//	HashCoarsening coarsen;
//	return coarsen.getNumberOfRemovablePoints(this->storage);
//}
//
//void BoundaryGridGenerator::refineMaxLevel(RefinementFunctor* func, unsigned int maxLevel)
//{
//	HashRefinementBoundariesMaxLevel refine;
//	refine.refineToMaxLevel(this->storage, func, maxLevel);
//}
//
//size_t BoundaryGridGenerator::getNumberOfRefinablePointsToMaxLevel(unsigned int maxLevel)
//{
//	HashRefinementBoundariesMaxLevel refine;
//	return refine.getNumberOfRefinablePointsToMaxLevel(this->storage, maxLevel);
//}

}
}

