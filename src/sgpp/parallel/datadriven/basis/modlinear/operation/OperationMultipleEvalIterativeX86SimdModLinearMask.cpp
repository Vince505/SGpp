/******************************************************************************
* Copyright (C) 2013 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "parallel/datadriven/basis/modlinear/operation/OperationMultipleEvalIterativeX86SimdModLinearMask.hpp"
#include "parallel/datadriven/basis/modlinear/operation/impl/X86SimdModLinearMask.hpp"
#include "parallel/tools/PartitioningTool.hpp"

namespace sg
{
namespace parallel
{

OperationMultipleEvalIterativeX86SimdModLinearMask::OperationMultipleEvalIterativeX86SimdModLinearMask(
		sg::base::GridStorage* storage, sg::base::DataMatrix* dataset) :
	sg::parallel::OperationMultipleEvalVectorized(storage, dataset)
{
	this->level_ = new sg::base::DataMatrix(storage_->size(), storage_->dim());
	this->index_ = new sg::base::DataMatrix(storage_->size(), storage_->dim());
	this->mask_ = new sg::base::DataMatrix(storage_->size(), storage_->dim());
	this->offset_ = new sg::base::DataMatrix(storage_->size(), storage_->dim());
	
	storage_->getLevelIndexMaskArraysForModEval(*(this->level_), *(this->index_), *(this->mask_), *(this->offset_));
}

void OperationMultipleEvalIterativeX86SimdModLinearMask::rebuildLevelAndIndex()
{
	delete this->level_;
	delete this->index_;
	delete this->mask_;
	delete this->offset_;

	this->level_ = new sg::base::DataMatrix(storage_->size(), storage_->dim());
	this->index_ = new sg::base::DataMatrix(storage_->size(), storage_->dim());
	this->mask_ = new sg::base::DataMatrix(storage_->size(), storage_->dim());
	this->offset_ = new sg::base::DataMatrix(storage_->size(), storage_->dim());

	storage_->getLevelIndexMaskArraysForModEval(*(this->level_), *(this->index_), *(this->mask_), *(this->offset_));
}

double OperationMultipleEvalIterativeX86SimdModLinearMask::multTransposeVectorized(sg::base::DataVector& source, sg::base::DataVector& result)
{
	myTimer_->start();
	result.setAll(0.0);

	#pragma omp parallel
	{
		size_t start;
		size_t end;
		PartitioningTool::getOpenMPPartitionSegment(m_gridFrom, m_gridTo, &start, &end, 1);

		X86SimdModLinearMask::multTranspose(
			level_, index_, mask_, offset_, dataset_, source, result, start, end, 0, this->dataset_->getNcols());
	}

	return myTimer_->stop();
}

double OperationMultipleEvalIterativeX86SimdModLinearMask::multVectorized(sg::base::DataVector& alpha, sg::base::DataVector& result)
{
	myTimer_->start();
	result.setAll(0.0);

	#pragma omp parallel
	{
		size_t start;
		size_t end;
		PartitioningTool::getOpenMPPartitionSegment(m_datasetFrom, m_datasetTo, &start, &end, X86SimdModLinearMask::getChunkDataPoints());

		X86SimdModLinearMask::mult(
			level_, index_, mask_, offset_, dataset_, alpha, result, 0, alpha.getSize(), start, end);
	}

	return myTimer_->stop();
}

}
}
