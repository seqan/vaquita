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
#include "pairedend.hpp"

void PairedEndRead::getCandidateRegion(TPosition &min, TPosition &max, BreakpointEvidence::SIDE toward)
{
    TPosition temp;

    if (toward == BreakpointEvidence::SIDE::BOTH)
    {
        BreakpointCandidate::setPositionWithAdj( min, max, this->getOptionManager()->getPairedEndSearchSize() );        
    }
    else if (toward == BreakpointEvidence::SIDE::RIGHT)
    {
        temp = min;
        BreakpointCandidate::setPositionWithAdj( min, temp, this->getOptionManager()->getMinSVSize() );
        temp = max;
        BreakpointCandidate::setPositionWithAdj( temp, max, this->getOptionManager()->getPairedEndSearchSize() );
    }
    else
    {
        temp = min;
        BreakpointCandidate::setPositionWithAdj( min, temp, this->getOptionManager()->getPairedEndSearchSize() );        
        temp = max;
        BreakpointCandidate::setPositionWithAdj( temp, max, this->getOptionManager()->getMinSVSize() );
    }
}

bool PairedEndRead::analyze(void)
{
    // change overlapped regions to breakpoint candidates
    //TBreakpointSet* breakpoints = this->tempBreakpoints.getCandidateSet(); 
    TBreakpointSet* breakpoints = this->getCandidateSet(); 
    auto it= breakpoints->begin();
    while ( it != breakpoints->end() )
    {
        Breakpoint* bp = *it;
        TPosition leftPos, rightPos, leftMin, leftMax, rightMin, rightMax;

        // get candidate region
        if (bp->orientation == BreakpointEvidence::ORIENTATION::PROPERLY_ORIENTED_LARGE || \
            bp->orientation == BreakpointEvidence::ORIENTATION::PROPERLY_ORIENTED_SMALL)
        {
            leftPos = bp->maxLeftPos;
            rightPos = bp->minRightPos;
            leftMin = leftPos;
            leftMax = leftPos;
            rightMin = rightPos;
            rightMax = rightPos;
            this->getCandidateRegion(leftMin, leftMax, BreakpointEvidence::SIDE::RIGHT);
            this->getCandidateRegion(rightMin, rightMax, BreakpointEvidence::SIDE::LEFT);
        }
        else if (bp->orientation == BreakpointEvidence::ORIENTATION::INVERTED)
        {
            if (bp->leftReverseFlag) // <-- <--
            {
                leftPos = bp->minLeftPos;
                rightPos = bp->minRightPos;
                leftMin = leftPos;
                leftMax = leftPos;
                rightMin = rightPos;
                rightMax = rightPos;
                this->getCandidateRegion(leftMin, leftMax, BreakpointEvidence::SIDE::BOTH);
                this->getCandidateRegion(rightMin, rightMax, BreakpointEvidence::SIDE::BOTH);
            }
            else // --> -->
            {
                leftPos = bp->maxLeftPos;
                rightPos = bp->maxRightPos;
                leftMin = leftPos;
                leftMax = leftPos;
                rightMin = rightPos;
                rightMax = rightPos;
                this->getCandidateRegion(leftMin, leftMax, BreakpointEvidence::SIDE::BOTH);
                this->getCandidateRegion(rightMin, rightMax, BreakpointEvidence::SIDE::BOTH);
           }
        }
        else if (bp->orientation == BreakpointEvidence::ORIENTATION::SWAPPED)
        {
            leftPos = bp->minLeftPos;
            rightPos = bp->maxRightPos;
            leftMin = leftPos;
            leftMax = leftPos;
            rightMin = rightPos;
            rightMax = rightPos;
            this->getCandidateRegion(leftMin, leftMax, BreakpointEvidence::SIDE::LEFT);
            this->getCandidateRegion(rightMin, rightMax, BreakpointEvidence::SIDE::RIGHT);
        }
        else // BreakpointEvidence::ORIENTATION::NOT_DECIDED
        {
            leftPos = bp->maxLeftPos;
            rightPos = bp->minRightPos;
            leftMin = leftPos;
            leftMax = leftPos;
            rightMin = rightPos;
            rightMax = rightPos;
            this->getCandidateRegion(leftMin, leftMax, BreakpointEvidence::SIDE::RIGHT);
            this->getCandidateRegion(rightMin, rightMax, BreakpointEvidence::SIDE::LEFT);
        }

        // they are in the same strand (used for NOT_DECIDED cases)
        if (bp->leftReverseFlag != bp->rightReverseFlag)
            bp->rightReverseFlag = bp->leftReverseFlag;   
        else
            bp->rightReverseFlag = (bp->leftReverseFlag == true) ? false : true;

        // avoid overlap
        if (leftMax > rightPos)
            leftMax = rightPos - 1;

        if (rightMin < leftPos)
            rightMin = leftPos + 1;

        // update region
        bp->leftPos.clear();
        bp->leftPos.push_back(leftMin);
        bp->leftPos.push_back(leftMax);
        bp->rightPos.clear();
        bp->rightPos.push_back(rightMin);
        bp->rightPos.push_back(rightMax);
        bp->needLeftIndexUpdate = true;
        bp->needRightIndexUpdate = true;
        this->updateBreakpointIndex(bp);

        // update position wihtout updating index
        bp->leftPos.clear();
        bp->rightPos.clear();
        bp->leftPos.push_back(leftPos);
        bp->rightPos.push_back(rightPos);

        // init.
        bp->bFoundExactPosition = false;
        this->setBreakpointUsed(bp, false);

        // to the next
        ++it;
    }

    return true;
}

bool PairedEndRead::isNew(TReadName &readName)
{
    return (this->pairInfo.find(readName) == this->pairInfo.end());
}

void PairedEndRead::parseReadRecord(TReadName &readName, seqan::BamAlignmentRecord &record)
{
    // can't specify the positions (should be filtered in record reading step)
    if (record.rID == seqan::BamAlignmentRecord::INVALID_REFID || record.rNextId == seqan::BamAlignmentRecord::INVALID_REFID)
        return;

    bool isNew = (this->pairInfo.find(readName) == this->pairInfo.end());
    BreakpointEvidence& be = this->pairInfo[readName]; // [start, end)
    unsigned readSize = getAlignmentLengthInRef(record);

    if (isNew)
    {
        be.leftSegment.templateID = BreakpointEvidence::NOVEL_TEMPLATE;
        be.rightSegment.templateID = BreakpointEvidence::NOVEL_TEMPLATE;
    }

    if (record.tLen >= 0)
    {
        be.leftSegment.beginPos = record.beginPos;
        be.leftSegment.endPos = record.beginPos + getAlignmentLengthInRef(record);
        be.leftSegment.templateID = record.rID;
        be.leftSegment.isReverse = hasFlagRC(record);
    }
    else
    {
        be.rightSegment.beginPos = record.beginPos;
        be.rightSegment.endPos = record.beginPos + getAlignmentLengthInRef(record);
        be.rightSegment.templateID = record.rID;
        be.rightSegment.isReverse = hasFlagRC(record);
    }

    // pair-information is completed
    if ( be.leftSegment.templateID != BreakpointEvidence::NOVEL_TEMPLATE && \
         be.rightSegment.templateID != BreakpointEvidence::NOVEL_TEMPLATE )
    {
        be.suppRead = this->getNextReadID(); // this must be unique
        if (be.leftSegment.templateID == be.rightSegment.templateID)
        {
            if (be.leftSegment.isReverse != be.rightSegment.isReverse)
            {
                if (be.leftSegment.isReverse)
                    be.orientation = BreakpointEvidence::ORIENTATION::SWAPPED; // <-- -->
                else
                {
                    // --> <--
                    if (abs(record.tLen) > this->getMaxAbInsSize())
                        be.orientation = BreakpointEvidence::ORIENTATION::PROPERLY_ORIENTED_LARGE;
                    else
                        be.orientation = BreakpointEvidence::ORIENTATION::PROPERLY_ORIENTED_SMALL;
                }
            }
            else
            {
                be.orientation = BreakpointEvidence::ORIENTATION::INVERTED; // --> --> or <-- <--
            }
        }
        else
        {
            be.orientation = BreakpointEvidence::ORIENTATION::NOT_DECIDED;
        }

        Breakpoint* bp = this->updateBreakpoint(be, isNew);
        this->pairInfo.erase(readName); // erase
    }
}