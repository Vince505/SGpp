/* ****************************************************************************
 * Copyright (C) 2009 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 **************************************************************************** */
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de), Dirk Pflueger (pflueged@in.tum.de)

#ifndef HASHREFINEMENT_HPP
#define HASHREFINEMENT_HPP

#include "base/grid/GridStorage.hpp"
#include "base/grid/generation/functors/RefinementFunctor.hpp"
#include "base/grid/generation/hashmap/AbstractRefinement.hpp"
#include <sstream>
#include <map>

namespace sg {
namespace base {

//@TODO (lettrich, high) documentation
struct compareLevels
{
	bool operator()(AbstractRefinement::index_type gridPointA, AbstractRefinement::index_type gridPointB)
	{
		size_t dim = gridPointA.dim();
		//std::cout<< "comparing " << gridPointA.toString() << " & " << gridPointB.toString();
		//iterate over all dimensions
		for(size_t i = 0; i< dim; ++i)
		{
			//compare level in each dimension
			//std::cout << gridPointA.getLevel(i) << "<" << gridPointB.getLevel(i) << "? , ";
			if (gridPointA.getLevel(i) < gridPointB.getLevel(i)){
				//std::cout << " - SMALLER\n";
				return true;
			}else if (gridPointA.getLevel(i) > gridPointB.getLevel(i)) {
				//std::cout << " - LARGER\n";
				return false;
			}
		}
		//std::cout << " - LARGER\n";
		return false;
	}
};

class errorContainer
{

public:

	errorContainer(RefinementFunctor::value_type error, size_t contributionCounter)
{
		this->error = error;
		this->contributionCounter = contributionCounter;
}

	errorContainer(RefinementFunctor::value_type error)
	{
		this->error = error;
		this->contributionCounter = 1;
	}

	errorContainer()
	{
		error = 0;
		contributionCounter = 0;
	}

	RefinementFunctor::value_type getError()
	{
		return error;
	}

	size_t getContributionCounter()
	{
		return contributionCounter;
	}

	void operator= (RefinementFunctor::value_type newError)
	{
		error = newError;
		contributionCounter = 1;
	}

	errorContainer operator+ (RefinementFunctor::value_type newError)
	{
		errorContainer result = *this;
		result.error += newError;
		result.contributionCounter++;

		return result;
	}

	errorContainer operator+ (const errorContainer& container)
	{
		errorContainer result = *this;
		result.error += container.error;
		result.contributionCounter += container.contributionCounter;

		return result;
	}

	errorContainer operator+= (RefinementFunctor::value_type newError)
   			{
		this->error += newError;
		this->contributionCounter ++;
		return *this;
   			}

	errorContainer& operator+= (const errorContainer& container)
   			{
		this->error += container.error;
		this->contributionCounter += container.contributionCounter;
		return *this;
   			}


	errorContainer operator- (RefinementFunctor::value_type newError)
	{
		errorContainer result = *this;
		result.error -= newError;
		result.contributionCounter++;

		return result;
	}

	errorContainer operator- (const errorContainer& container)
	{
		errorContainer result = *this;
		result.error -= container.error;
		result.contributionCounter += container.contributionCounter;

		return result;
	}

	bool operator== (RefinementFunctor::value_type otherValue)
				{
		return error == otherValue ? true : false ;
				}

	bool operator== (errorContainer& container)
				{
		return error == container.error ? true : false ;
				}

	bool operator> (RefinementFunctor::value_type otherValue)
	{
		return error > otherValue ? true : false ;
	}

	bool operator> (errorContainer& container)
	{
		return error > container.error ? true : false ;
	}

	bool operator< (RefinementFunctor::value_type otherValue)
	{
		return error < otherValue ? true : false ;
	}

	bool operator< (errorContainer& container)
	{
		return error < container.error ? true : false ;
	}

	RefinementFunctor::value_type getContribPerPoint()
	{
		return error/static_cast<double>(contributionCounter);
	}

	std::string toString()
	{
		std::ostringstream tmp;
		tmp << "error: " << error << ", iterations: " << contributionCounter << ",contrib: " << getContribPerPoint() <<  "\n";
		return  tmp.str();
	}

private:

	RefinementFunctor::value_type error;
	size_t contributionCounter;


};




/**
 * Free refinement class for sparse grids
 */
class HashRefinement: public AbstractRefinement {

public:

	typedef std::map<index_type,errorContainer,compareLevels> SubspaceError;

	/**
	 * Refines a grid according to a RefinementFunctor provided.
	 * Refines up to RefinementFunctor::getRefinementsNum() grid points if
	 * possible, and if their refinement value is larger than RefinementFunctor::start()
	 * and their absolute value is larger or equal than RefinementFunctor::getRefinementThreshold()
	 *
	 * @param storage hashmap that stores the grid points
	 * @param functor a RefinementFunctor specifying the refinement criteria
	 */
	void free_refine(GridStorage* storage, RefinementFunctor* functor);

	/**
	 * Refines a grid by adding additional Subspaces according to a RefinementFunctor provided.
	 * Refines up to RefinementFunctor::getRefinementsNum() grid points if
	 * possible, and if their refinement value is larger than RefinementFunctor::start()
	 * and their absolute value is larger or equal than RefinementFunctor::getRefinementThreshold()
	 *
	 * @param storage hashmap that stores the grid points
	 * @param functor a RefinementFunctor specifying the refinement criteria
	 */
	void freeRefineSubspace(GridStorage* storage, RefinementFunctor* functor);

	/**
	 * Computes and returns the number of grid points, which can be refined.
	 * This is the number of grid points that have at least one child missing.
	 *
	 * @param storage hashmap that stores the grid points
	 * @return The number of grid points that can be refined
	 */
	size_t getNumberOfRefinablePoints(GridStorage* storage);

	/**
	 * Refine one grid point along a single direction
	 * @param storage hashmap that stores the grid points
	 * @param index point to refine
	 * @param d direction
	 */
	void refineGridpoint1D(GridStorage* storage, index_type& index, size_t d);

	void refineGridpoint1D(GridStorage* storage, HashGridIndex< unsigned int, unsigned int >* index, size_t d) {
		refineGridpoint1D(storage, *index, d);
	}

protected:
	/**
	 * This method refines a grid point by generating the children in every dimension
	 * of the grid and all their missing ancestors by calling create_gridpoint().
	 *
	 * @param storage hashmap that stores the gridpoints
	 * @param refine_index The index in the hashmap of the point that should be refined
	 */
	void refineGridpoint(GridStorage* storage, size_t refine_index);

	/**
	 * This method creates a new point on the grid. It checks if some parents or
	 * children are needed in other dimensions.
	 *
	 * @param storage hashmap that stores the gridpoints
	 * @param index The point that should be inserted
	 */
	void createGridpoint(GridStorage* storage, index_type& index);

	/**
	 * This method creates a new subspace in the grid. It checks if some parents or
	 * children are needed in other dimensions.
	 *
	 * @param storage hashmap that stores the gridpoints
	 * @param index - the index containing the level vector for the new subspace.
	 * the index vector of the object has to be set to 1 for all dimensions!.
	 */
	void createSubspace(GridStorage* storage, index_type& index);

	/**
	 * Examines the grid points and stores the indices those that can be refined
	 * and have maximal indicator values.
	 *
	 * @param storage hashmap that stores the grid points
	 * @param functor a RefinementFunctor specifying the refinement criteria
	 * @param refinements_num number of points to refine
	 * @param max_indices the array where the point indices should be stored
	 * @param max_values the array where the corresponding indicator values
	 * should be stored
	 */
	virtual void collectRefinablePoints(GridStorage* storage,
			RefinementFunctor* functor, size_t refinements_num, size_t* max_indices,
			RefinementFunctor::value_type* max_values);

	/**
	 * Examines the grid points, finds which ones are refinable and adds
	 * the error indicators of all points which belong to the same subspace.
	 * It returns an unsroted array with the refinements_num subspaces with the highest error indicator. (not sorted)
	 *
	 * @param storage hashmap that stores the grid points
	 * @param functor a RefinementFunctor specifying the refinement criteria
	 * @param refinements_num number of points to refine
	 * @param maxErrorSubspaces the array of hash grid indices, containing the level vector of the subspace
	 * @param maxErrorValues the array with the indicator values corresponding to the level of the subspace.
	 */
	virtual void collectRefinableSubspaces(GridStorage* storage,
			RefinementFunctor* functor,
			size_t refinements_num,
			AbstractRefinement::index_type* maxSubspaces,
			RefinementFunctor::value_type* maxErrors);

	/**
	 * Refines the collection of points.
	 *
	 * @param storage hashmap that stores the grid points
	 * @param functor a RefinementFunctor specifying the refinement criteria
	 * @param refinements_num number of points to refine
	 * @param max_indices the array with the indices of points that should be refined
	 * @param max_values the array with the corresponding indicator values
	 */
	virtual void refineGridpointsCollection(GridStorage* storage,
			RefinementFunctor* functor, size_t refinements_num, size_t* max_indices,
			RefinementFunctor::value_type* max_values);
	/**
	 * Creates all subspaces that are passed in the maxErrorSubspaces array according to
	 * the specifications in the refinement functor.
	 *
	 * @param storage hashmap that stores the grid points
	 * @param functor a RefinementFunctor specifying the refinement criteria
	 * @param refinements_num number of points to refine
	 * @param maxErrorSubspaces the array with the indices containing the level vectors
	 *  on which the subspaces should be created
	 * @param maxErrorValues the array with the corresponding indicator values
	 */
	virtual void refineSubspaceCollection(GridStorage* storage,
			RefinementFunctor* functor,
			size_t refinements_num,
			index_type* maxErrorSubspaces,
			RefinementFunctor::value_type* maxErrorValues);
private:

	/**
	 * recursive function to create all points on a subspace.
	 *
	 * @param storage
	 * @param index
	 * @param dim
	 *
	 */
	void createSubspaceHelper(GridStorage* storage, index_type& index, size_t dim);

	/**
	 * Sets all the elements of the index vector of a grid point to 1.
	 * @param gridPoint pointer to the grid point with the index array to be changed
	 */
	void resetIndexVector(index_type* gridPoint);

	void testAdmissibility(SubspaceError* subspaceError);
};
}
}

#endif /* HASHREFINEMENT_HPP */
