// ==========================================================================
//                               Vaquita
// ==========================================================================
// Copyright (c) 2016, Jongkyu Kim, MPI-MolGen/FU-Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Jongkyu Kim or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL JONGKYU KIM OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Jongkyu Kim <j.kim@fu-berlin.de>
// ==========================================================================
#ifndef APP_INTERVALINDEX_H_
#define APP_INTERVALINDEX_H_

#include <utility>
#include <limits>
#include <vector>
#include <map>
#include <set>
#include <stdint.h>

	/* DEBUG
	#include "bin-interval-index.hpp"
	bool testfunc(int* a, int* b)
	{
	    std::cerr << "hi\n";
	    return true;
	}

    int a = 5;
    int b = 6;
    int c = 7;
    int d = 3;

    BinIntervalIndex<int*> test(100, testfunc);
    test.add(10000,20000, &a);
    test.add(19000,22000, &b);
    test.add(21000,25000, &c);

    std::set<int*> testset;
    test.find(testset, 10000, 25000, &d);
    std::cerr << "1 start\n";
    for(auto it = testset.begin(); it != testset.end(); ++it)
        std::cerr << *(*it) << "\n";

    test.remove(20000, 22000, &b);
    test.find(testset, 20000, 22000, &d);
    std::cerr << "2 start\n";
    for(auto it = testset.begin(); it != testset.end(); ++it)
        std::cerr << *(*it) << "\n";

    test.find(testset, 22001, 22000, &d);
    std::cerr << "3 start\n";
    for(auto it = testset.begin(); it != testset.end(); ++it)
        std::cerr << *(*it) << "\n";

    test.find(testset, 22001, 26000, &d);
    std::cerr << "4 start\n";
    for(auto it = testset.begin(); it != testset.end(); ++it)
        std::cerr << *(*it) << "\n";

    std::cerr << "hi\n";
    return 0;
    // DEBUG */

template <typename TRawData>
class IntervalIndex
{
    typedef uint32_t                           TBin;
    typedef uint32_t                           TPosition;
	typedef std::pair<TPosition, TPosition>    TInterval;
	typedef std::pair<TInterval, TRawData>     TRawDataWithInterval;
	typedef std::vector<TRawDataWithInterval>  TDataVector;
	typedef std::set<TRawData>                 TRawDataSet;

	private :
		static TBin const INVALID_POS = std::numeric_limits<TBin>::max();
		TPosition binSize;
        TPosition adjacency;
		std::map<TBin, TDataVector*> index;

		inline TBin getBin(TPosition p) { return (p / this->binSize); }
		bool (*isMatchedFunc)(TRawData, TRawData);

	public :
		IntervalIndex(TPosition binSize, TPosition adjacency, bool (*isMatchedFunc)(TRawData, TRawData));
		~IntervalIndex();

		bool add(const TPosition, const TPosition, const TRawData);
		bool remove(TPosition, TPosition, TRawData);
		bool find(TRawDataSet&, const TPosition, const TPosition, const TRawData);
		
        void getPositionWithAdj(TPosition &left, TPosition &right);
		static inline bool isValidInterval(TPosition, TPosition);
		static inline bool isOverlap(TPosition, TPosition, TPosition, TPosition);
        static inline bool isAdjacent(TPosition, TPosition, TPosition);
        inline bool isAdjacent(TPosition, TPosition, TPosition, TPosition);
        inline bool isAdjacent(TPosition, TPosition);
};

template <typename TRawData>
bool IntervalIndex<TRawData>::isOverlap(TPosition beginPos1, TPosition endPos1, TPosition beginPos2, TPosition endPos2)
{
    return (beginPos1 <= endPos2 && beginPos2 <= endPos1);
}

template <typename TRawData>
bool IntervalIndex<TRawData>::isAdjacent(TPosition beginPos1, TPosition endPos1, TPosition beginPos2, TPosition endPos2)
{
    return ( (endPos1 <= beginPos2 and ((beginPos2 - endPos1) <= this->adjacency)) or \
             (endPos2 <= beginPos1 and ((beginPos1 - endPos2) <= this->adjacency)) );
}
template <typename TRawData>
bool IntervalIndex<TRawData>::isAdjacent(TPosition p1, TPosition p2, TPosition tol)
{
    return (((p2 >= p1) && (p2 - p1) <= tol) || ((p1 > p2) && (p1 - p2) <= tol));
}

template <typename TRawData>
bool IntervalIndex<TRawData>::isAdjacent(TPosition p1, TPosition p2)
{
    return (((p2 >= p1) && (p2 - p1) <= this->adjacency) || ((p1 > p2) && (p1 - p2) <= this->adjacency));
}

template <typename TRawData>
inline bool IntervalIndex<TRawData>::isValidInterval(TPosition begin, TPosition end)
{
	return (begin != IntervalIndex::INVALID_POS && \
		end != IntervalIndex::INVALID_POS && \
		(begin <= end));
}

template <typename TRawData>
bool IntervalIndex<TRawData>::add(const TPosition begin, const TPosition end, const TRawData data)
{
	if (isValidInterval(begin, end) == false)
		return false;

	TBin minBin, maxBin;
	minBin = this->getBin(begin);
	maxBin = this->getBin(end);

    for (TBin bin=minBin; bin <= maxBin; ++bin)
    {
        // new bin
        if (this->index.find(bin) == this->index.end())
            this->index.insert( std::make_pair(bin, new TDataVector) );

        // add index        
        this->index[bin]->push_back( std::make_pair(std::make_pair(begin, end), data) );
    }

    return true;
}

template <typename TRawData>
bool IntervalIndex<TRawData>::remove(const TPosition begin, const TPosition end, const TRawData data)
{
	if (isValidInterval(begin, end) == false)
		return false;

	TBin minBin, maxBin;
	minBin = this->getBin(begin);
	maxBin = this->getBin(end);

	for (TBin bin = minBin; bin <= maxBin; ++bin)
    {
        if (this->index.find(bin) == this->index.end())
        	continue;

        auto it = this->index[bin]->begin();
        while (it != this->index[bin]->end())
        {
            if ( it->second == data )
                it = this->index[bin]->erase(it);
            else
                ++it;
        }
    }

    return true;
}

template <typename TRawData>
void IntervalIndex<TRawData>::getPositionWithAdj(TPosition &left, TPosition &right)
{
    if (left < this->adjacency)
        left = 0;
    else
        left -= this->adjacency;

    if ( right > IntervalIndex::INVALID_POS - this->adjacency)
        right = IntervalIndex::INVALID_POS;
    else
        right += this->adjacency;
}

template <typename TRawData>
bool IntervalIndex<TRawData>::find(TRawDataSet& outVector, const TPosition begin, const TPosition end, const TRawData data)
{
	if (isValidInterval(begin, end) == false)
		return false;
  
    TPosition beginWithAdj = begin;
    TPosition endWithAdj = end;
    this->getPositionWithAdj(beginWithAdj, endWithAdj);

	TBin minBin, maxBin;
	minBin = this->getBin(beginWithAdj);
	maxBin = this->getBin(endWithAdj);

	for (TBin bin = minBin; bin <= maxBin; ++bin)
    {
        if (this->index.find(bin) == this->index.end())
        	continue;

		auto it = this->index[bin]->begin();
      	while (it != this->index[bin]->end())
       	{
       		if ( isOverlap(beginWithAdj, endWithAdj, it->first.first, it->first.second) && (this->isMatchedFunc)(data, it->second) == true)
            {
                outVector.insert(it->second);
            }
            /*
            else if ( isAdjacent(begin, end, it->first.first, it->first.second) && (this->isMatchedFunc)(data, it->second) == true)
                outVector.insert(it->second);
            */
            ++it;
        }
    }

    return true;
}

template <typename TRawData>
IntervalIndex<TRawData>::IntervalIndex(TPosition binSize, TPosition adj, bool (*isMatchedFunc)(TRawData, TRawData))
{
	this->binSize = binSize;
	this->adjacency = adj;
    this->isMatchedFunc = isMatchedFunc;
}

template <typename TRawData>
IntervalIndex<TRawData>::~IntervalIndex()
{
	for(auto it = index.begin(); it != index.end(); ++it)
		delete it->second;
}
#endif // APP_INTERVALINDEX_H_