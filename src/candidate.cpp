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
                            //be.orientation = BreakpointEvidence::ORIENTATION::LEFT_CLIP;
                            be.orientation = BreakpointEvidence::ORIENTATION::NOT_DECIDED;
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
                            //be.orientation = BreakpointEvidence::ORIENTATION::RIGHT_CLIP;
                            be.orientation = BreakpointEvidence::ORIENTATION::NOT_DECIDED;
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
    // delete indicies
    TBreakpointBinIndex *leftIndex, *rightIndex;
    this->getIndex(&leftIndex, &rightIndex, *bpIt);
    this->removeIndex(leftIndex, (*bpIt)->minLeftPos, (*bpIt)->maxLeftPos, *bpIt);
    this->removeIndex(rightIndex, (*bpIt)->minRightPos, (*bpIt)->maxRightPos, *bpIt);
    
    delete *bpIt;
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

    return this->updateBreakpoint(newBp, true, isNew);
}

Breakpoint* BreakpointCandidate::moveAndUpdateBreakpoint(Breakpoint* bp, bool& isNew)
{
    Breakpoint* newBp = new Breakpoint;
    this->moveBreakpoint(*newBp, *bp);

    return this->updateBreakpoint(newBp, true, isNew);
}


void BreakpointCandidate::findMatchedBreakpoint(TBreakpointSet& leftMatched, TBreakpointSet& rightMatched, Breakpoint* bp, bool checkOrientation)
{    
    TBreakpointBinIndex *leftIndex, *rightIndex;
    TPosition minPos, maxPos;

    // get index & search interval tress
    this->getIndex(&leftIndex, &rightIndex, bp);

    // same strand
    bool isSameStrand = (bp->leftReverseFlag == bp->rightReverseFlag);

    if (bp->leftTemplateID != BreakpointEvidence::NOVEL_TEMPLATE)
    {
        minPos = bp->minLeftPos;
        maxPos = bp->maxLeftPos;
        this->setPositionWithAdj(minPos, maxPos);
        this->findOverlap(leftMatched, leftIndex, minPos, maxPos, bp->orientation, checkOrientation);
    }

    if (bp->rightTemplateID != BreakpointEvidence::NOVEL_TEMPLATE)
    {
        minPos = bp->minRightPos;
        maxPos = bp->maxRightPos;
        this->setPositionWithAdj(minPos, maxPos);
        this->findOverlap(rightMatched, rightIndex, minPos, maxPos, bp->orientation, checkOrientation);
    }
}


void BreakpointCandidate::findMatchedBreakpoint(TBreakpointSet& leftMatched, TBreakpointSet& rightMatched, TBreakpointSet& bpEquivalent, Breakpoint* bp, bool checkBothSide)
{
    if (bp->orientation == BreakpointEvidence::ORIENTATION::NOT_DECIDED)
    {
        // do not check the orientation (last boolean)
        this->findMatchedBreakpoint(leftMatched, rightMatched, bp, false);
    }
    else
    {
        // check the coordinates and orientation (last boolean)
        this->findMatchedBreakpoint(leftMatched, rightMatched, bp, true);
    }

    // find equivalent breakpoin
    if (bp->leftTemplateID == BreakpointEvidence::NOVEL_TEMPLATE)
    {
        // use right-matched set only
        bpEquivalent = rightMatched;
    }
    else if (bp->rightTemplateID == BreakpointEvidence::NOVEL_TEMPLATE)
    {
        // use left-matched set only
        bpEquivalent = leftMatched;
    }
    else
    {
        // check both sides
        /*
        TBreakpointSet::iterator itLeft, itRight;
        for (itLeft = leftMatched.begin(); itLeft != leftMatched.end(); ++itLeft)
            for (itRight = rightMatched.begin(); itRight != rightMatched.end(); ++itRight)
                if ( *itLeft == *itRight)
                    bpEquivalent.insert(*itLeft);
        */
        /*
        for (auto it = leftMatched.begin(); it != leftMatched.end(); ++it)
            bpEquivalent.insert(*it);
        for (auto it = rightMatched.begin(); it != rightMatched.end(); ++it)
            bpEquivalent.insert(*it);
        */

        if (checkBothSide == true)
        {
            for (auto it = leftMatched.begin(); it != leftMatched.end(); ++it)
                if (rightMatched.find(*it) != rightMatched.end())
                    bpEquivalent.insert(*it);
        }
        else 
        {
            for (auto it = leftMatched.begin(); it != leftMatched.end(); ++it)
                bpEquivalent.insert(*it);
            for (auto it = rightMatched.begin(); it != rightMatched.end(); ++it)
                bpEquivalent.insert(*it);
        }
    }
}

Breakpoint* BreakpointCandidate::updateBreakpoint(Breakpoint* bp, bool checkBothSide, bool& isNew)
{
    // find matched breakpoints
    TBreakpointSet leftMatched, rightMatched, bpEquivalent;
    this->findMatchedBreakpoint(leftMatched, rightMatched, bpEquivalent, bp, checkBothSide);

    // merge found matches
    if (bpEquivalent.size() > 0)
    {
        // Select one as the reprentative breakpoint(RB)
        auto itEq = bpEquivalent.begin();
        Breakpoint* destBp = *(itEq++);

        // merge new bp to B
        mergeBreakpoint(destBp, bp);
        delete bp;

        for (; itEq != bpEquivalent.end(); ++itEq) // merge others to RB
            mergeBreakpoint(destBp, *itEq);

        // update index
        updateBreakpointIndex(destBp);
        isNew = false;
        return destBp;
    } 
    else // or add (new breakpoint)
    {
        this->addBreakpoint(bp);
        isNew = true;
        return bp;
    }    
}

Breakpoint* BreakpointCandidate::updateBreakpoint(BreakpointEvidence& be, bool checkBothSide, bool& isNew)
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

    return this->updateBreakpoint(bp, checkBothSide, isNew);
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
    TBreakpointBinIndex *leftIndex, *rightIndex;
    this->getIndex(&leftIndex, &rightIndex, bp);

    TPosition minPos, maxPos;
    if (bp->needLeftIndexUpdate == true)
    {
        // remove
        minPos = bp->minLeftPos;
        maxPos = bp->maxLeftPos;
        this->removeIndex(leftIndex, minPos, maxPos, bp);

        // update & add
        this->updateLeftMinMaxPos(bp);
        this->addIndex(leftIndex, bp->minLeftPos, bp->maxLeftPos, bp);
        bp->needLeftIndexUpdate = false;
    }

    if (bp->needRightIndexUpdate == true)
    {
        minPos = bp->minRightPos;
        maxPos = bp->maxRightPos;
        this->removeIndex(rightIndex, minPos, maxPos, bp);

        this->updateRightMinMaxPos(bp);
        this->addIndex(rightIndex, bp->minRightPos, bp->maxRightPos, bp);
        bp->needRightIndexUpdate = false;
    }
}

void BreakpointCandidate::addBreakpoint(Breakpoint* bp)
{
    // register
    this->breakpoints.insert(bp);
 
    // update index
    TBreakpointBinIndex *leftIndex, *rightIndex;
    this->getIndex(&leftIndex, &rightIndex, bp);

    this->updateLeftMinMaxPos(bp);
    this->addIndex(leftIndex, bp->minLeftPos, bp->maxLeftPos, bp);

    this->updateRightMinMaxPos(bp);
    this->addIndex(rightIndex, bp->minRightPos, bp->maxRightPos, bp);

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

void BreakpointCandidate::getIndex(TBreakpointBinIndex** leftIndex, TBreakpointBinIndex** rightIndex, Breakpoint* bp)
{
    *leftIndex = &(this->leftBpIndex[bp->leftTemplateID]);
    *rightIndex = &(this->rightBpIndex[bp->rightTemplateID]);
}

void BreakpointCandidate::addIndex(TBreakpointBinIndex* index, TPosition minPos, TPosition maxPos, Breakpoint* bp)
{
    if( minPos == BreakpointEvidence::INVALID_POS && maxPos == BreakpointEvidence::INVALID_POS)
        return;

    if (minPos == BreakpointEvidence::INVALID_POS)
        minPos = maxPos;
    else if (maxPos == BreakpointEvidence::INVALID_POS)
        maxPos = minPos;

    // select a range of positional bin
    int32_t minPosBin, maxPosBin;
    minPosBin = GENOMIC_BIN(minPos, BreakpointCandidate::GENOMIC_BIN_SIZE);
    maxPosBin = GENOMIC_BIN(maxPos, BreakpointCandidate::GENOMIC_BIN_SIZE);

    for (int32_t bin = minPosBin; bin <= maxPosBin; ++bin)
    {
        // make a new bin
        if (index->find(bin) == index->end())
            index->insert( std::make_pair(bin, new TPosBreakpointVector) );

        // add index        
        (*index)[bin]->push_back( std::make_pair(std::make_pair(minPos, maxPos), bp) );
    }
}

void BreakpointCandidate::removeIndex(TBreakpointBinIndex* index, TPosition minPos, TPosition maxPos, Breakpoint* bp)
{
    if( minPos == BreakpointEvidence::INVALID_POS && maxPos == BreakpointEvidence::INVALID_POS)
        return;

    if (minPos == BreakpointEvidence::INVALID_POS)
        minPos = maxPos;
    else if (maxPos == BreakpointEvidence::INVALID_POS)
        maxPos = minPos;

    int32_t minPosBin, maxPosBin;
    minPosBin = GENOMIC_BIN(minPos, BreakpointCandidate::GENOMIC_BIN_SIZE);
    maxPosBin = GENOMIC_BIN(maxPos, BreakpointCandidate::GENOMIC_BIN_SIZE);
 
    for (int32_t bin = minPosBin; bin <= maxPosBin; ++bin)
    {
        if (index->find(bin) != index->end())
        {
            TPosBreakpointVector::iterator it = (*index)[bin]->begin();
            while( it != (*index)[bin]->end() )
            {
                if ( it->second == bp )
                    it = (*index)[bin]->erase(it);
                else
                    ++it;
            }
        }
    }
}

void BreakpointCandidate::findOverlap(TBreakpointSet& foundSet, TBreakpointBinIndex* index, TPosition minPos, TPosition maxPos, BreakpointEvidence::ORIENTATION orientation, bool checkOrientation)
{
    if( minPos == BreakpointEvidence::INVALID_POS && maxPos == BreakpointEvidence::INVALID_POS)
        return;

    if (minPos == BreakpointEvidence::INVALID_POS)
        minPos = maxPos;
    else if (maxPos == BreakpointEvidence::INVALID_POS)
        maxPos = minPos;

    int32_t minPosBin, maxPosBin, adj;
    minPosBin = GENOMIC_BIN(minPos, BreakpointCandidate::GENOMIC_BIN_SIZE) - 0;
    maxPosBin = GENOMIC_BIN(maxPos, BreakpointCandidate::GENOMIC_BIN_SIZE) + 0;
    for (int32_t bin = minPosBin; bin <= maxPosBin; ++bin)
    {
        if (index->find(bin) != index->end())
        {
            TPosBreakpointVector::iterator it = (*index)[bin]->begin();
            while( it != (*index)[bin]->end() )
            {
                // it->first.first : begin pos
                // it->first.second : end pos
                // it->second : breakpoint pointer
                // 0-based, [begin, end)
                
                
                /*
                if ( isOverlap(it->first.first, it->first.second, minPos, maxPos) || \
                     (IS_ADJACENT(it->first.first, minPos, adj) || \
                      IS_ADJACENT(it->first.second, maxPos, adj)))
                */
                if (IS_OVERLAP(it->first.first, it->first.second, minPos, maxPos))
                {
                    Breakpoint* bp = it->second;
                    if (checkOrientation == true)
                    {
                        if ((bp->orientation) == orientation)
                            foundSet.insert(bp);
                    }
                    else
                        foundSet.insert(bp);
                }
                ++it;
            }
        }
    }
}

// 0-based, [begin, end)
bool BreakpointCandidate::isOverlap(TPosition beginPos1, TPosition endPos1, TPosition beginPos2, TPosition endPos2)
{
    return ((beginPos1 <= beginPos2 and beginPos2 < endPos1) or \
           (beginPos2 <= beginPos1 and beginPos1 < endPos2) or \
           (beginPos1 <= beginPos2 and endPos2 < endPos1) or \
           (beginPos2 <= beginPos1 and endPos1 < endPos2));
}

bool BreakpointCandidate::isOverlap(SequenceSegment& s1, SequenceSegment& s2)
{
    /*
    return (s1.isReverse == s2.isReverse) and (s1.templateID == s2.templateID) and \
           BreakpointCandidate::isOverlap(s1.beginPos, s1.endPos, s2.beginPos, s2.endPos);
    */
    return BreakpointCandidate::isOverlap(s1.beginPos, s1.endPos, s2.beginPos, s2.endPos);
}

// 0-based, [begin, end)
bool BreakpointCandidate::isAdjacent(TPosition beginPos1, TPosition endPos1, TPosition beginPos2, TPosition endPos2, TPosition tol)
{
    return ( (endPos1 <= beginPos2 and ((beginPos2 - endPos1) <= tol)) or \
             (endPos2 <= beginPos1 and ((beginPos1 - endPos2) <= tol)) );
}

bool BreakpointCandidate::isAdjacent(TPosition pos1, TPosition pos2, TPosition tol)
{
    return ( abs(pos2-pos1) <= tol );
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
    // TODO : release memory!
    //(breakpoints)
    std::set<Breakpoint*>::iterator it;
    for(it = this->breakpoints.begin(); it != this->breakpoints.end(); ++it)
        delete *it;
} 

BreakpointCandidate::BreakpointCandidate()
{
    
}

bool BreakpointCandidate::isMatchedBreakpoint(Breakpoint*, Breakpoint*)
{
    return true;
}