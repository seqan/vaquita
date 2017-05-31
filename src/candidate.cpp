// ==========================================================================
//                               Vaquita
// ==========================================================================
// Copyright (c) 2017, Jongkyu Kim, MPI-MolGen/FU-Berlin
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
#include <algorithm>
#include <iterator>
#include "misc.hpp"
#include "candidate.hpp"

void BreakpointCandidate::parseCIGAR(AlignmentInfo &alnInfo, TReadID &readID, BamAlignmentRecord &record, \
                                     bool checkIndel, bool checkClip, TPosition minSVSize, TPosition minClippedSeqSize)
{
    // reference position 0-based, [begin, end)
    alnInfo.refSegment.templateID = record.rID;
    alnInfo.refSegment.beginPos = record.beginPos;
    alnInfo.refSegment.isReverse = hasFlagRC(record);
    alnInfo.querySegment.templateID = record._qID; // TODO : SeqAn seems not working properly here. (qID)
    alnInfo.querySegment.isReverse = alnInfo.refSegment.isReverse;

    // can't be decided  yet
    TPosition refBeginPos = record.beginPos; 
    TPosition refEndPos = record.beginPos; 
    TPosition queryBeginPos = BreakpointEvidence::INVALID_POS;
    TPosition queryEndPos = BreakpointEvidence::INVALID_POS; 
    TPosition querySize = 0;

    // parse CIGAR string
    BreakpointEvidence be;
    for (unsigned i=0; i != length(record.cigar); ++i)
    { 
        switch (record.cigar[i].operation)
        {
            case 'M' : // alignment match
            case '=' : // sequence match
            case 'X' : // sequence mismatch
                refEndPos += record.cigar[i].count;
                if(queryBeginPos == BreakpointEvidence::INVALID_POS) // start position in query
                {
                    queryBeginPos = 0;
                    queryEndPos = record.cigar[i].count;
                }
                else
                    queryEndPos += record.cigar[i].count;
                querySize += record.cigar[i].count;
                break;
            case 'S' : // soft clip 
                if ( checkClip == true )
                {
                    if (record.cigar[i].count >= minClippedSeqSize)
                    {
                        be.leftSegment.isReverse = hasFlagRC(record);
                        be.rightSegment.isReverse = be.leftSegment.isReverse;
                        be.suppRead = 1;

                        if (i == 0) // left side clip
                        {
                            be.leftSegment.templateID = BreakpointEvidence::NOVEL_TEMPLATE;
                            be.leftSegment.beginPos = BreakpointEvidence::INVALID_POS;
                            be.leftSegment.endPos = BreakpointEvidence::INVALID_POS;
                            //be.leftSegment.templateID = record.rID;
                            //be.leftSegment.beginPos = refBeginPos;
                            //be.leftSegment.endPos = be.rightSegment.beginPos + 1;
                            be.rightSegment.templateID = record.rID;
                            be.rightSegment.beginPos = refBeginPos;
                            be.rightSegment.endPos = be.rightSegment.beginPos + 1;
                            be.rightSegment.isReverse = alnInfo.refSegment.isReverse;
                            be.sequence = prefix(record.seq, record.cigar[i].count);
                            be.orientation = BreakpointEvidence::ORIENTATION::CLIPPED;
                            //be.orientation = BreakpointEvidence::ORIENTATION::NOT_DECIDED;
                            alnInfo.clippedList.push_back(be);
                        }
                        
                        if (i == (length(record.cigar) - 1)) // right side clip
                        {
                            be.leftSegment.templateID = record.rID;
                            be.leftSegment.endPos = refEndPos;
                            be.leftSegment.beginPos = be.leftSegment.endPos - 1;    
                            be.leftSegment.isReverse = alnInfo.refSegment.isReverse;    
                            //be.rightSegment.templateID = record.rID;
                            //be.rightSegment.endPos = refEndPos;
                            //be.rightSegment.beginPos = be.leftSegment.endPos - 1;                   
                            be.rightSegment.templateID = BreakpointEvidence::NOVEL_TEMPLATE;
                            be.rightSegment.beginPos = BreakpointEvidence::INVALID_POS;
                            be.rightSegment.endPos = BreakpointEvidence::INVALID_POS;                               
                            be.sequence = suffix(record.seq, length(record.seq) - record.cigar[i].count);
                            be.orientation = BreakpointEvidence::ORIENTATION::CLIPPED;
                            //be.orientation = BreakpointEvidence::ORIENTATION::NOT_DECIDED;
                            alnInfo.clippedList.push_back(be);
                        }
                    }
                }
            case 'H' : // hard clip
                if (queryBeginPos == BreakpointEvidence::INVALID_POS) // start position in query
                {
                    queryBeginPos = record.cigar[i].count;
                    queryEndPos = queryBeginPos;
                }
                querySize += record.cigar[i].count;
                break;
            case 'P' : // padding (like insertion)
            case 'I' : // TODO : insertion
                queryEndPos += record.cigar[i].count;
                if (queryBeginPos == BreakpointEvidence::INVALID_POS) // start position in query
                    queryBeginPos = 0;
                querySize += record.cigar[i].count;
                break;
            case 'D' : // (reference) deletion
            case 'N' : // (reference) skip
                if (queryBeginPos == BreakpointEvidence::INVALID_POS) // start position in query
                {
                    queryBeginPos = record.cigar[i].count;
                    queryEndPos = queryBeginPos;
                }
            
                // find indels equal or larger than minSVSize (consider deletions only)
                if ( checkIndel == true && record.cigar[i].count >= minSVSize )
                {
                    TPosition delRefBeginPos = refEndPos;
                    refEndPos += record.cigar[i].count;

                    be.leftSegment.templateID = record.rID;
                    be.leftSegment.beginPos = delRefBeginPos;
                    be.leftSegment.endPos = be.leftSegment.beginPos + 1;
                    be.leftSegment.isReverse = hasFlagRC(record);

                    be.rightSegment.templateID = record.rID;
                    be.rightSegment.beginPos = refEndPos;
                    be.rightSegment.endPos = be.rightSegment.beginPos + 1;
                    be.rightSegment.isReverse = be.leftSegment.isReverse;

                    be.orientation = BreakpointEvidence::ORIENTATION::PROPERLY_ORIENTED_LARGE;
                    be.suppRead = 1;
                    alnInfo.indelList.push_back(be);
                }
                refEndPos += record.cigar[i].count;
        }
    } // end parse CIGAR
    
    // reference position 0-based, [begin, end)
    alnInfo.refSegment.endPos = refEndPos;
    alnInfo.querySize = querySize;
    alnInfo.querySegment.beginPos = queryBeginPos;
    alnInfo.querySegment.endPos = queryEndPos;
}

TBreakpointSet::iterator BreakpointCandidate::removeBreakpoint(Breakpoint* bp)
{
    return this->removeBreakpoint(this->breakpoints.find(bp));
}

TBreakpointSet::iterator BreakpointCandidate::removeBreakpoint(TBreakpointSet::iterator bpIt)
{
    // get indices
    TBreakpointIntervalIndex *leftIndex = getIndex(BreakpointEvidence::LEFT, (*bpIt)->leftTemplateID);
    TBreakpointIntervalIndex *rightIndex = getIndex(BreakpointEvidence::RIGHT, (*bpIt)->rightTemplateID);

    // remove indices
    leftIndex->remove((*bpIt)->minLeftPos, (*bpIt)->maxLeftPos, *bpIt);
    rightIndex->remove((*bpIt)->minRightPos, (*bpIt)->maxRightPos, *bpIt);
    
    // destroy breakpoint instance
    delete *bpIt;

    // delete in the set & return iterator
    return this->breakpoints.erase(bpIt);
}

void BreakpointCandidate::copyBreakpoint(Breakpoint& dest, Breakpoint& src)
{
    dest.leftTemplateID = src.leftTemplateID;
    dest.rightTemplateID = src.rightTemplateID;
    dest.leftReverseFlag = src.leftReverseFlag;
    dest.rightReverseFlag = src.rightReverseFlag;
    dest.orientation = src.orientation;

    dest.minLeftPos = src.minLeftPos;
    dest.maxLeftPos = src.maxLeftPos;
    dest.minRightPos = src.minRightPos;
    dest.maxRightPos = src.maxRightPos;

    dest.needLeftIndexUpdate = src.needLeftIndexUpdate;
    dest.needRightIndexUpdate = src.needRightIndexUpdate;
    dest.bFoundExactPosition = src.bFoundExactPosition;
    dest.clippedConsensusSequenceSize = src.clippedConsensusSequenceSize;

    dest.leftPos.insert(dest.leftPos.end(), src.leftPos.begin(), src.leftPos.end());
    dest.rightPos.insert(dest.rightPos.end(), src.rightPos.begin(), src.rightPos.end());
    dest.suppReads = src.suppReads;
    dest.clippedSequences.insert(dest.clippedSequences.end(), src.clippedSequences.begin(), src.clippedSequences.end());
}

void BreakpointCandidate::moveBreakpoint(Breakpoint& dest, Breakpoint& src)
{
    dest.leftTemplateID = src.leftTemplateID;
    dest.rightTemplateID = src.rightTemplateID;
    dest.leftReverseFlag = src.leftReverseFlag;
    dest.rightReverseFlag = src.rightReverseFlag;
    dest.orientation = src.orientation;

    dest.minLeftPos = src.minLeftPos;
    dest.maxLeftPos = src.maxLeftPos;
    dest.minRightPos = src.minRightPos;
    dest.maxRightPos = src.maxRightPos;

    dest.needLeftIndexUpdate = src.needLeftIndexUpdate;
    dest.needRightIndexUpdate = src.needRightIndexUpdate;
    dest.bFoundExactPosition = src.bFoundExactPosition;
    dest.clippedConsensusSequenceSize = src.clippedConsensusSequenceSize;

    std::move(src.leftPos.begin(), src.leftPos.end(), std::back_inserter(dest.leftPos));
    std::move(src.rightPos.begin(), src.rightPos.end(), std::back_inserter(dest.rightPos));
    dest.suppReads = src.suppReads;
    std::move(src.clippedSequences.begin(), src.clippedSequences.end(), std::back_inserter(dest.clippedSequences));
}

Breakpoint* BreakpointCandidate::copyAndUpdateBreakpoint(Breakpoint* bp, bool& isNew)
{
    Breakpoint* newBp = new Breakpoint;
    this->copyBreakpoint(*newBp, *bp);

    return this->updateBreakpoint(newBp, isNew);
}

Breakpoint* BreakpointCandidate::moveAndUpdateBreakpoint(Breakpoint* bp, bool& isNew)
{
    Breakpoint* newBp = new Breakpoint;
    this->moveBreakpoint(*newBp, *bp);

    return this->updateBreakpoint(newBp, isNew);
}


void BreakpointCandidate::findBreakpoint(TBreakpointSet& leftSideFound, TBreakpointSet& rightSideFound, TBreakpointSet& bothSideFound, Breakpoint* bp)
{
    // find breakpoints
    if (bp->leftTemplateID == BreakpointEvidence::NOVEL_TEMPLATE)
    {
        TBreakpointIntervalIndex *index = getIndex(BreakpointEvidence::RIGHT, bp->rightTemplateID);
        index->find(rightSideFound, bp->minRightPos, bp->maxRightPos, bp);
        bothSideFound = rightSideFound;
    }
    else if (bp->rightTemplateID == BreakpointEvidence::NOVEL_TEMPLATE)
    {
        // use left-matched set only
        TBreakpointIntervalIndex *index = getIndex(BreakpointEvidence::LEFT, bp->leftTemplateID);
        index->find(leftSideFound, bp->minLeftPos, bp->maxLeftPos, bp);
        bothSideFound = leftSideFound;
    }
    else
    {
        TBreakpointIntervalIndex *leftIndex = getIndex(BreakpointEvidence::LEFT, bp->leftTemplateID);
        TBreakpointIntervalIndex *rightIndex = getIndex(BreakpointEvidence::RIGHT, bp->rightTemplateID);
        leftIndex->find(leftSideFound, bp->minLeftPos, bp->maxLeftPos, bp);
        rightIndex->find(rightSideFound, bp->minRightPos, bp->maxRightPos, bp);

        // both side found
        for (auto it = leftSideFound.begin(); it != leftSideFound.end(); ++it)
            if (rightSideFound.find(*it) != rightSideFound.end())
                bothSideFound.insert(*it);
    }
}


Breakpoint* BreakpointCandidate::updateBreakpoint(Breakpoint* bp, bool& isNew)
{
    // find matched breakpoints
    TBreakpointSet leftSideFound, rightSideFound, bothSideFound;
    this->findBreakpoint(leftSideFound, rightSideFound, bothSideFound, bp);

    // merge found matches
    if (bothSideFound.size() > 0)
    {
        // Select one as the reprentative breakpoint(RB)
        auto itBothSide = bothSideFound.begin();
        Breakpoint* destBp = *(itBothSide++);

        // merge new bp to B
        mergeBreakpoint(destBp, bp);
        delete bp;

        for (; itBothSide != bothSideFound.end(); ++itBothSide) // merge others to RB
            mergeBreakpoint(destBp, *itBothSide);

        // update index
        updateBreakpointIndex(destBp);
        isNew = false;
        return destBp;
    } 
    else // or add (new breakpoint)
    {
        this->addNewBreakpoint(bp);
        isNew = true;
        return bp;
    }    
}

Breakpoint* BreakpointCandidate::updateBreakpoint(BreakpointEvidence& be, bool& isNew)
{
    // a breakpoint candidate. [begin, end]
    Breakpoint* bp = new Breakpoint;

    // left side
    if ( (be.leftSegment.endPos - be.leftSegment.beginPos) > 1)
    {          
        // segment contains multiple candidate positions
        bp->minLeftPos = be.leftSegment.beginPos;
        bp->maxLeftPos = be.leftSegment.endPos - 1; 
        bp->leftPos.push_back(bp->minLeftPos);
        bp->leftPos.push_back(bp->maxLeftPos);
    }
    else
    { 
        bp->minLeftPos = be.leftSegment.beginPos;
        bp->maxLeftPos = be.leftSegment.beginPos;
        bp->leftPos.push_back(bp->minLeftPos);

    }
    bp->leftTemplateID = be.leftSegment.templateID; 
    bp->leftReverseFlag = be.leftSegment.isReverse;
    bp->needLeftIndexUpdate = false;

    // right side
    if ( (be.rightSegment.endPos - be.rightSegment.beginPos) > 1)
    {   
        // segment contains multiple candidate positions
        bp->minRightPos = be.rightSegment.beginPos;
        bp->maxRightPos = be.rightSegment.endPos - 1; 
        bp->rightPos.push_back(bp->minRightPos);
        bp->rightPos.push_back(bp->maxRightPos);
    }
    else
    { 
        bp->minRightPos = be.rightSegment.beginPos;
        bp->maxRightPos = be.rightSegment.beginPos;
        bp->rightPos.push_back(bp->minRightPos);
    }
    bp->rightTemplateID = be.rightSegment.templateID;
    bp->rightReverseFlag = be.rightSegment.isReverse;
    bp->needRightIndexUpdate = false;

    // orientation
    bp->orientation = be.orientation;

    // supporting reads
    bp->suppReads = 1;

    // clipped sequences
    if(bp->leftTemplateID == BreakpointEvidence::NOVEL_TEMPLATE)
        bp->clippedSequences.push_back( std::move(std::make_pair(be.rightSegment.beginPos, be.sequence)) );
    if(bp->rightTemplateID == BreakpointEvidence::NOVEL_TEMPLATE)
        bp->clippedSequences.push_back( std::move(std::make_pair(be.leftSegment.endPos-1, be.sequence)) );

    return this->updateBreakpoint(bp, isNew);
}

void BreakpointCandidate::clearPosInfo(Breakpoint *bp)
{
    bp->leftPos.clear();
    bp->rightPos.clear();
    bp->minLeftPos = BreakpointEvidence::INVALID_POS;
    bp->maxLeftPos = 0;
    bp->minRightPos = BreakpointEvidence::INVALID_POS;
    bp->maxRightPos = 0;
}

void BreakpointCandidate::updateMinMaxPos(Breakpoint* bp)
{
    updateLeftMinMaxPos(bp);
    updateRightMinMaxPos(bp);
}

void BreakpointCandidate::updateLeftMinMaxPos(Breakpoint* bp)
{
    sort(bp->leftPos.begin(), bp->leftPos.end());
    GET_MIN_MAX_ELEMENT(bp->minLeftPos, bp->maxLeftPos, bp->leftPos);
}

void BreakpointCandidate::updateRightMinMaxPos(Breakpoint* bp)
{
    sort(bp->rightPos.begin(), bp->rightPos.end());
    GET_MIN_MAX_ELEMENT(bp->minRightPos, bp->maxRightPos, bp->rightPos);
}

void BreakpointCandidate::updateBreakpointIndex(Breakpoint* bp)
{
    if (bp->needLeftIndexUpdate == true)
    {
        TBreakpointIntervalIndex *index = getIndex(BreakpointEvidence::LEFT, bp->leftTemplateID);

        // remove
        index->remove(bp->minLeftPos, bp->maxLeftPos, bp);

        // update & add
        this->updateLeftMinMaxPos(bp);
        bp->needLeftIndexUpdate = false;

        index->add(bp->minLeftPos, bp->maxLeftPos, bp);
        
    }

    if (bp->needRightIndexUpdate == true)
    {
        TBreakpointIntervalIndex *index = getIndex(BreakpointEvidence::RIGHT, bp->rightTemplateID);
        index->remove(bp->minRightPos, bp->maxRightPos, bp);

        this->updateRightMinMaxPos(bp);
        bp->needRightIndexUpdate = false;

        index->add(bp->minRightPos, bp->maxRightPos, bp);
    }
}

void BreakpointCandidate::addNewBreakpoint(Breakpoint* bp)
{
    // register
    this->breakpoints.insert(bp);
 
    // update index
    TBreakpointIntervalIndex *leftIndex = this->getIndex(BreakpointEvidence::LEFT, bp->leftTemplateID);
    TBreakpointIntervalIndex *rightIndex = this->getIndex(BreakpointEvidence::RIGHT, bp->rightTemplateID);
    this->updateMinMaxPos(bp);
    leftIndex->add(bp->minLeftPos, bp->maxLeftPos, bp);
    rightIndex->add(bp->minRightPos, bp->maxRightPos, bp);

    // lazy updating
    bp->needLeftIndexUpdate = false;
    bp->needRightIndexUpdate = false;
}

TBreakpointSet::iterator BreakpointCandidate::mergeBreakpoint(Breakpoint* destBp, Breakpoint* srcBp)
{
    // destBp update
    std::move(srcBp->leftPos.begin(), srcBp->leftPos.end(), std::back_inserter(destBp->leftPos));
    if ( destBp->needLeftIndexUpdate == false && (destBp->minLeftPos > srcBp->minLeftPos || destBp->maxLeftPos < srcBp->maxLeftPos) )
        destBp->needLeftIndexUpdate = true;

    std::move(srcBp->rightPos.begin(), srcBp->rightPos.end(), std::back_inserter(destBp->rightPos));
    if ( destBp->needRightIndexUpdate == false && (destBp->minRightPos > srcBp->minRightPos || destBp->maxRightPos < srcBp->maxRightPos) )
        destBp->needRightIndexUpdate = true;

    // supporting reads
    destBp->suppReads += srcBp->suppReads;

    // clipped sequences
    if (this->op->isUsingAssembly())
    {
        if (destBp->clippedConsensusSequenceSize < srcBp->clippedConsensusSequenceSize)
            destBp->clippedConsensusSequenceSize = srcBp->clippedConsensusSequenceSize;
        std::move(srcBp->clippedSequences.begin(), srcBp->clippedSequences.end(), std::back_inserter(destBp->clippedSequences));
    }
    else
    {
        if (destBp->clippedConsensusSequenceSize < srcBp->clippedConsensusSequenceSize)
        {
            destBp->clippedConsensusSequenceSize = srcBp->clippedConsensusSequenceSize;
            destBp->clippedSequences.clear();
            std::move(srcBp->clippedSequences.begin(), srcBp->clippedSequences.end(), std::back_inserter(destBp->clippedSequences));
        }
    }
    destBp->bFoundExactPosition = (destBp->bFoundExactPosition || srcBp->bFoundExactPosition);
    
    // do additional jobs
    doAdditionalJobAfterMerge(destBp, srcBp);

    // srcBp removal
    auto it =  this->breakpoints.find(srcBp);
    if( it != this->breakpoints.end())
        return this->removeBreakpoint(it);
    else
        return this->breakpoints.end();
}

TBreakpointIntervalIndex* BreakpointCandidate::getIndex(BreakpointEvidence::SIDE side, TTemplateID id)
{
    if (side == BreakpointEvidence::LEFT)
    {
        if (leftIndexMap.find(id) == leftIndexMap.end())
        {
            leftIndexMap.insert(std::make_pair(id, \
                                new IntervalIndex<Breakpoint*> (GENOMIC_BIN_SIZE, op->getAdjTol(), this->isMatchedBreakpoint)));
        }
        
        return leftIndexMap[id];
    }
    else if (side == BreakpointEvidence::RIGHT)
    {
        if (rightIndexMap.find(id) == rightIndexMap.end())
        {
            rightIndexMap.insert(std::make_pair(id, \
                                new IntervalIndex<Breakpoint*> (GENOMIC_BIN_SIZE, op->getAdjTol(), this->isMatchedBreakpoint)));
        }

        return rightIndexMap[id];
    }
    else
    {
        printTimeMessage("Requesting invalid side/id: " + std::to_string((int)side) + "/" + std::to_string((int)id));
        exit(1);
    }
}

bool BreakpointCandidate::addIndex(TBreakpointIntervalIndex* index, TPosition minPos, TPosition maxPos, Breakpoint* bp)
{
    return index->add(minPos, maxPos, bp);
}

// 0-based, [begin, end)
bool BreakpointCandidate::isOverlap(TPosition beginPos1, TPosition endPos1, TPosition beginPos2, TPosition endPos2)
{
    return TBreakpointIntervalIndex::isOverlap(beginPos1, endPos1, beginPos2, endPos2);

}

bool BreakpointCandidate::isOverlap(SequenceSegment& s1, SequenceSegment& s2)
{
    return TBreakpointIntervalIndex::isOverlap(s1.beginPos, s1.endPos, s2.beginPos, s2.endPos);
}

// 0-based, [begin, end)
bool BreakpointCandidate::isAdjacent(TPosition beginPos1, TPosition endPos1, TPosition beginPos2, TPosition endPos2, TPosition tol)
{
    if (endPos1 <= beginPos2)
        return TBreakpointIntervalIndex::isAdjacent(beginPos2, endPos1, tol);
    else if (endPos2 <= beginPos1)
        return TBreakpointIntervalIndex::isAdjacent(beginPos1, endPos2, tol);
    return false;
}

bool BreakpointCandidate::isAdjacent(TPosition pos1, TPosition pos2, TPosition tol)
{
    return TBreakpointIntervalIndex::isAdjacent(pos1, pos2, tol);
}

bool BreakpointCandidate::isAdjacent(SequenceSegment& s1, SequenceSegment& s2)
{    
    return BreakpointCandidate::isAdjacent(s1.beginPos, s1.endPos, s2.beginPos, s2.endPos, this->getPositionalAdj()); 
}

void BreakpointCandidate::setPositionWithAdj(TPosition &left, TPosition &right, TPosition adj)
{
    if ( left < adj)
        left = 0;
    else
        left -= adj;

    if ( right > BreakpointEvidence::INVALID_POS - adj)
        right = BreakpointEvidence::INVALID_POS;
    else
        right += adj;
}

void BreakpointCandidate::setPositionWithAdj(TPosition &left, TPosition &right)
{
    setPositionWithAdj(left, right, this->getPositionalAdj());
    //setPositionWithAdj(left, right, this->getOptionManager()->getAdjTol());
}

void BreakpointCandidate::printBreakpoint(Breakpoint& bp)
{
    std::cerr << "<LEFT>" << "\n";
    std::cerr << "ID, LFLAG, NEED_UPDATE:" << "\t" << bp.leftTemplateID << "," << bp.leftReverseFlag << "," << bp.needLeftIndexUpdate << "\n\tPOS: ";
    for(unsigned int i=0; i < bp.leftPos.size(); ++i)
        std::cerr << bp.leftPos[i] << ",";
    std::cerr << "\nMIN, MAX:\t" << bp.minLeftPos << "," << bp.maxLeftPos << "\n";
    std::cerr << "\n";

    std::cerr << "<RIGHT>" << "\n";
    std::cerr << "ID, RFLAG, NEED_UPDATE:" << "\t" << bp.rightTemplateID << "," << bp.rightReverseFlag << "," << bp.needRightIndexUpdate << "\n\tPOS: ";
    for(unsigned int i=0; i < bp.rightPos.size(); ++i)
        std::cerr << bp.rightPos[i] << ",";
    std::cerr << "\nMIN, MAX:\t" << bp.minRightPos << "," << bp.maxRightPos << "\n";
    std::cerr << "\n";

    std::cerr << "<ORIENTATION>" << "\n\t";
    std::cerr << bp.orientation << "\n\n";

    std::cerr << "<READS>" << "\n\t";
    std::cerr << bp.suppReads << "\n\n";
    /*
    std::set<TReadID>::iterator it;
    for(it = bp.suppReads.begin(); it != bp.suppReads.end(); ++it)
        std::cerr << *it << ",";
    */
    std::cerr << "\n\n";
}

void BreakpointCandidate::printBreakpoint(Breakpoint *bp)
{
    printBreakpoint(*bp);
}

bool BreakpointCandidate::compareByQueryPos(AlignmentInfo& a1, AlignmentInfo& a2)
{
    return a1.querySegment.beginPos < a2.querySegment.beginPos;
}

bool BreakpointCandidate::compareByChrmAndPos(Breakpoint& bp1, Breakpoint& bp2)
{
    if (bp1.leftTemplateID != bp2.leftTemplateID)
        return bp1.leftTemplateID < bp2.leftTemplateID;

    if (bp1.minLeftPos != bp2.minLeftPos)
        return bp1.minLeftPos < bp2.minLeftPos;

    return bp1.maxRightPos < bp2.maxRightPos;
}

void BreakpointCandidate::setInsertionInfo(double median, double stdev, double min, double max)
{
    this->insertMedian = median;
    this->insertDev = stdev;
    this->minAbInsSize = min;
    this->maxAbInsSize = max;
}

BreakpointCandidate::~BreakpointCandidate()
{
    // release memory   
    for(auto it = this->breakpoints.begin(); it != this->breakpoints.end(); ++it)
        delete *it;
    breakpoints.clear();

    for(auto it = this->leftIndexMap.begin(); it != this->leftIndexMap.end(); ++it)
        delete it->second;
    leftIndexMap.clear();

    for(auto it = this->rightIndexMap.begin(); it != this->rightIndexMap.end(); ++it)
        delete it->second;
    rightIndexMap.clear();
} 

BreakpointCandidate::BreakpointCandidate(OptionManager* op)
{
    setOptionManager(op);
}

bool BreakpointCandidate::isMatchedBreakpoint(Breakpoint* a, Breakpoint* b)
{
    if (a->orientation == BreakpointEvidence::ORIENTATION::CLIPPED)
        return true;
    else
        return (a->orientation == b->orientation);
}
