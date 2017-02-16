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
#include <fstream>
#include <math.h>  
#include <algorithm>
#include <iterator>
#include <map>
#include <seqan/find.h>
#include "misc.hpp"
#include "breakpoint.hpp"

bool BreakpointManager::merge(void)
{
    bool result = false;

    RUN(result,"From split-read evidences.", mergeSplitRead());    
    if (this->optionManager->doPairedEndAnalysis())
        RUN(result,"From read-pair evidences.", mergePairedEnd());

    if (this->optionManager->doClippedReadAnalysis())
        RUN(result,"From soft-clipped evidences.", mergeClippedRead());

    // 2) add imprecise breakpoints
    if (this->optionManager->doPairedEndAnalysis())
        RUN(result,"Add imprecise breakpoints to the merged set.", addImpreciseBreakpoints());

    // 3) get final positions
    RUN(result,"Get representative coordinatese.", findFinalBreakpoints());

    // 4) add read-depth based information (needs final breakpoints)
    if (this->optionManager->doReadDepthAnalysis())
        RUN(result,"Calcuate read-depth information", calculateReadDepth());
}

bool BreakpointManager::mergeSplitRead(void)
{
    // split-read
    TBreakpointSet* breakpoints = this->splitReadBreakpoints.getCandidateSet();
    auto itBreakpoint = breakpoints->begin();;
    while ( itBreakpoint != breakpoints->end() )
    {
        bool isNew;
        Breakpoint* currentBp = this->mergedBreakpoints.moveAndUpdateBreakpoint(*itBreakpoint, isNew);
        ReadSupportInfo* readInfo = this->mergedBreakpoints.getReadSupport(currentBp);

        readInfo->splitReadSupport = 0;
        readInfo->pairedEndSupport = 0;
        readInfo->clippedReadSupport = 0;
        readInfo->splitReadSupport = currentBp->suppReads;

        itBreakpoint = this->splitReadBreakpoints.removeBreakpoint(*itBreakpoint);
    }

    return true;
}

bool BreakpointManager::mergePairedEnd(void)
{    
    // paired-end
    TBreakpointSet* breakpoints = this->pairedEndBreakpoints.getCandidateSet();
    auto itBreakpoint = breakpoints->begin();;
 
    while ( itBreakpoint != breakpoints->end() )
    {
        Breakpoint* currentBp = *itBreakpoint;

        // find matches
        TBreakpointSet leftMatched, rightMatched, bpEquivalent;
        this->mergedBreakpoints.findMatchedBreakpoint(leftMatched, rightMatched, currentBp);
        for (auto itLeft = leftMatched.begin(); itLeft != leftMatched.end(); ++itLeft)
            for (auto itRight = rightMatched.begin(); itRight != rightMatched.end(); ++itRight)
                if ( *itLeft == *itRight)
                    bpEquivalent.insert(*itLeft);

        // found matches
        if (bpEquivalent.size() > 0)
        {
            // representative bp
            auto itEq = bpEquivalent.begin();
            Breakpoint* destBp = *(itEq++);
            ReadSupportInfo *readInfo = this->mergedBreakpoints.getReadSupport(destBp);

            // merge curret paired-end bp
            this->mergedBreakpoints.mergeBreakpoint(destBp, currentBp);
            readInfo->pairedEndSupport += currentBp->suppReads;

            // merge others that are linked by this paired-end bp
            for (; itEq != bpEquivalent.end(); ++itEq)
                this->mergedBreakpoints.mergeBreakpoint(destBp, *itEq);
            /*
            for (auto itEq=bpEquivalent.begin(); itEq != bpEquivalent.end(); ++itEq)
            {
                Breakpoint* destBp = *itEq;
                ReadSupportInfo* readInfo = this->mergedBreakpoints.getReadSupport(destBp);

                destBp->suppReads += currentBp->suppReads;
                readInfo->pairedEndSupport += currentBp->suppReads;
            }
            */
            this->pairedEndBreakpoints.setBreakpointUsed(currentBp, true);
            itBreakpoint = this->pairedEndBreakpoints.removeBreakpoint(currentBp);
        }
        else // This breakpoint is supported by paired-end reads only.
            ++itBreakpoint;
    }

    return true;
}

bool BreakpointManager::mergeClippedRead(void)
{
    // 1) clipped-read
    TBreakpointSet leftMatched, rightMatched;
    TBreakpointSet* breakpoints = this->clippedBreakpoints.getCandidateSet();
    auto itBreakpoint= breakpoints->begin();
    while ( itBreakpoint != breakpoints->end() )
    {
        Breakpoint* currentBp = *itBreakpoint;
        uint32_t clippedReadsCount = currentBp->suppReads;
        bool isLeftClip = (currentBp->leftTemplateID == BreakpointEvidence::NOVEL_TEMPLATE);

        leftMatched.clear();
        rightMatched.clear();

        // find matches from merged & paired-end set
        // last 'false' : do not check the orientations
        this->mergedBreakpoints.findMatchedBreakpoint(leftMatched, rightMatched, currentBp, false);        
        this->pairedEndBreakpoints.findMatchedBreakpoint(leftMatched, rightMatched, currentBp, false);

        CharString clipSeq = "";
        CharString revCompClipSeq = "";
        int32_t bestMatchScore;
        TFoundPosition foundPositions;
        foundPositions.clear();

        // search for candidate regions from paired-end information
        if ( (isLeftClip == true) && (rightMatched.size() > 0) || (isLeftClip == false) && (leftMatched.size() > 0))
        {               
            this->clippedBreakpoints.getConsensusSequence(clipSeq, currentBp);
            bestMatchScore = -(this->optionManager->getClippedSeqErrorRate() * length(clipSeq));

            if (isLeftClip == true)
            {
                for (auto itBp = rightMatched.begin(); itBp != rightMatched.end(); ++itBp )
                {
                    if ((*itBp)->orientation == BreakpointEvidence::ORIENTATION::INVERSED)
                    {
                        if (revCompClipSeq == "")
                        {
                            revCompClipSeq = clipSeq;
                            reverseComplement(revCompClipSeq);
                        }
                        this->clippedBreakpoints.searchPairRegion(foundPositions, *itBp, bestMatchScore, revCompClipSeq, ClippedRead::SIDE::LEFT, false, true);
                    }
                    else
                        this->clippedBreakpoints.searchPairRegion(foundPositions, *itBp, bestMatchScore, clipSeq, ClippedRead::SIDE::LEFT, false, false);
                }
            }
            else
            {
                for (auto itBp = leftMatched.begin(); itBp != leftMatched.end(); ++itBp )
                {
                    if ((*itBp)->orientation == BreakpointEvidence::ORIENTATION::INVERSED)
                    {
                        if (revCompClipSeq == "")
                        {
                            revCompClipSeq = clipSeq;
                            reverseComplement(revCompClipSeq);
                        }
                        this->clippedBreakpoints.searchPairRegion(foundPositions, *itBp, bestMatchScore, revCompClipSeq, ClippedRead::SIDE::RIGHT, false, true);
                    }
                    else
                        this->clippedBreakpoints.searchPairRegion(foundPositions, *itBp, bestMatchScore, clipSeq, ClippedRead::SIDE::RIGHT, false, false);
                }
            }
        }
 
        // search for twilight zone with clipped sequence
        if (foundPositions.size() == 0)
        {
            if (clipSeq == "")
            {
                this->clippedBreakpoints.getConsensusSequence(clipSeq, currentBp);
                bestMatchScore = -(this->optionManager->getClippedSeqErrorRate() * length(clipSeq));
            }

            if (isLeftClip == true)
                this->clippedBreakpoints.searchTwilightZone(foundPositions, currentBp, bestMatchScore, clipSeq, ClippedRead::SIDE::LEFT, false);
            else
                this->clippedBreakpoints.searchTwilightZone(foundPositions, currentBp, bestMatchScore, clipSeq, ClippedRead::SIDE::RIGHT, false);
        }

        // second try with reverse-complement
        if (foundPositions.size() == 0)
        {
            if (revCompClipSeq == "")
            {
                revCompClipSeq = clipSeq;
                reverseComplement(revCompClipSeq);
            }

            if (isLeftClip == true)
                this->clippedBreakpoints.searchTwilightZone(foundPositions, currentBp, bestMatchScore, revCompClipSeq, ClippedRead::SIDE::LEFT, true);
            else
                this->clippedBreakpoints.searchTwilightZone(foundPositions, currentBp, bestMatchScore, revCompClipSeq, ClippedRead::SIDE::RIGHT, true);
        }

        // found unique position
        if (foundPositions.size() > 0)
            addNewPositionsByClippedSequence(foundPositions, currentBp, isLeftClip);

        itBreakpoint = this->clippedBreakpoints.removeBreakpoint(currentBp);
    }
    return true;
}

void BreakpointManager::addNewPositionsByClippedSequence(TFoundPosition& foundPositions, Breakpoint* clippedBp, bool isLeftClip)
{
    // add new positions if any
    for (auto itFound = foundPositions.begin(); itFound != foundPositions.end(); ++itFound)
    {
        Breakpoint* matchedBp = itFound->matchedBp;
        bool isReverseComplemented = itFound->isReverseComplemented;

        // copy
        Breakpoint* newBp = new Breakpoint;
        BreakpointCandidate::copyBreakpoint(*newBp, *clippedBp);
        //BreakpointCandidate::moveBreakpoint(*newBp, *clippedBp);

        // get new exact position
        int32_t newPos = 0, svSize = 0;
        if (isLeftClip)
        {
            // check orientation
            if (isReverseComplemented)
            {
                newBp->orientation = BreakpointEvidence::ORIENTATION::INVERSED;
                newBp->leftReverseFlag = !newBp->rightReverseFlag;
                newPos = itFound->sequenceSegment.beginPos;
            }
            else
            {
                if (itFound->sequenceSegment.endPos < newBp->minRightPos)
                {
                    newBp->orientation = BreakpointEvidence::ORIENTATION::PROPERLY_ORIENTED;
                    newPos = itFound->sequenceSegment.endPos;
                }
                else
                {
                    newBp->orientation = BreakpointEvidence::ORIENTATION::SWAPPED;
                    newPos = itFound->sequenceSegment.endPos - 1;
                }
            }            

            // fill information
            if (matchedBp->leftTemplateID == BreakpointEvidence::NOVEL_TEMPLATE)
                newBp->leftTemplateID = newBp->rightTemplateID; // by twilight zone
            else
                newBp->leftTemplateID = matchedBp->leftTemplateID; // by paired-end region

            newBp->leftPos.clear();
            newBp->leftPos.push_back(newPos);            
            BreakpointCandidate::updateLeftMinMaxPos(newBp);

            // update depth
            if (this->optionManager->doReadDepthAnalysis())
            {
                this->readDepthBreakpoints.addUniformDepth(newBp->leftTemplateID, \
                                                           newPos, \
                                                           newBp->clippedConsensusSequenceSize,
                                                           newBp->suppReads);
            }
        }
        else
        {
            if (isReverseComplemented)
            {
                newBp->orientation = BreakpointEvidence::ORIENTATION::INVERSED;
                newBp->rightReverseFlag = !newBp->leftReverseFlag;
                newPos = itFound->sequenceSegment.endPos - 1;
            }
            else
            {
                if (newBp->maxLeftPos < itFound->sequenceSegment.beginPos)
                {
                    newBp->orientation = BreakpointEvidence::ORIENTATION::PROPERLY_ORIENTED;
                    newPos = itFound->sequenceSegment.beginPos - 1;
                }
                else
                {
                    newBp->orientation = BreakpointEvidence::ORIENTATION::SWAPPED;
                    newPos = itFound->sequenceSegment.beginPos;
                }
            }

            if (matchedBp->rightTemplateID == BreakpointEvidence::NOVEL_TEMPLATE)
                newBp->rightTemplateID = newBp->leftTemplateID;
            else
                newBp->rightTemplateID = matchedBp->rightTemplateID;


            newBp->rightTemplateID = newBp->leftTemplateID;
            newBp->rightPos.clear();
            newBp->rightPos.push_back(newPos);
            BreakpointCandidate::updateRightMinMaxPos(newBp);

            // update depth 
            if (this->optionManager->doReadDepthAnalysis())
            {
                this->readDepthBreakpoints.addUniformDepth(newBp->leftTemplateID, \
                                                           newPos, \
                                                           newBp->clippedConsensusSequenceSize,
                                                           newBp->suppReads);
            }
        }

        // add to the merged set
        bool isNewBp = false;
        newBp->needLeftIndexUpdate = true;
        newBp->needRightIndexUpdate = true;
        newBp = this->mergedBreakpoints.updateBreakpoint(newBp, isNewBp);

        // update supporting read information
        ReadSupportInfo* newReadInfo = this->mergedBreakpoints.getReadSupport(newBp);
        newReadInfo->clippedReadSupport += clippedBp->suppReads;
        if (isNewBp == true)
        {
            // supported by paired-end only
            if (matchedBp->leftTemplateID != BreakpointEvidence::NOVEL_TEMPLATE && \
                matchedBp->rightTemplateID != BreakpointEvidence::NOVEL_TEMPLATE)
            {
                this->pairedEndBreakpoints.setBreakpointUsed(matchedBp, true);
                newBp->suppReads += matchedBp->suppReads;
                newReadInfo->pairedEndSupport += matchedBp->suppReads;
            }
        }
   }
}

bool BreakpointManager::calculateReadDepth()
{
    TBreakpointSet* candidateSet = this->mergedBreakpoints.getCandidateSet();
    double avgReadDepth = 0;
    for(auto it = candidateSet->begin(); it != candidateSet->end(); ++it)
    {
        Breakpoint* bp = *it;
        ReadSupportInfo* info = this->mergedBreakpoints.getReadSupport(bp);
        FinalBreakpointInfo* finalBreakpoint = &this->finalBreakpoints[bp];

        int32_t breakpointSize = abs(finalBreakpoint->rightPosition - finalBreakpoint->leftPosition) + 1;
        int32_t windowSize = this->optionManager->getReadDepthWindowSize();

        // restrict the search region's size
        if (breakpointSize > windowSize)
            breakpointSize = windowSize;

        // select min or max value accoding to the orientation
        this->readDepthBreakpoints.getReadDepthDiffScore(info->leftReadDepth, \
                                                      info->rightReadDepth, \
                                                      info->leftReadDepthDiffScore, \
                                                      info->rightReadDepthDiffScore, \
                                                      finalBreakpoint->leftTemplateID, \
                                                      finalBreakpoint->leftPosition, \
                                                      finalBreakpoint->rightTemplateID, \
                                                      finalBreakpoint->rightPosition, \
                                                      windowSize, \
                                                      breakpointSize);

        if (info->leftReadDepthDiffScore > info->rightReadDepthDiffScore)
            info->leftReadDepthSelected = true;
        else
            info->leftReadDepthSelected = false;

        info->avgReadDepth = (info->leftReadDepth + info->rightReadDepth) / 2;
        avgReadDepth += info->avgReadDepth;
    }
    avgReadDepth = avgReadDepth / candidateSet->size();
    printTimeMessage("Average read-depth of the sample :" + std::to_string(avgReadDepth));
    this->optionManager->setAverageReadDepth(avgReadDepth);
}

bool BreakpointManager::find(void)
{
    bool result;

    // split-read analysis
    RUN(result,"Split-read analysis", this->splitReadBreakpoints.analyze());
    if (result == false) 
        return false;

    // paired-end analysis
    if (this->optionManager->doPairedEndAnalysis())
    {
        RUN(result,"Paired-end analysis", this->pairedEndBreakpoints.analyze());
        if (result == false) 
            return false;
    }
    
    // clipped-read analysis
    if (this->optionManager->doClippedReadAnalysis())
    {
        RUN(result,"Clipped-read analysis", this->clippedBreakpoints.analyze());
        if (result == false) 
            return false;
    }

    return true;
}

bool BreakpointManager::findFinalBreakpoints(void)
{
    TBreakpointSet* candidateSet = this->mergedBreakpoints.getCandidateSet();
    for(auto it = candidateSet->begin(); it != candidateSet->end(); ++it)
    {
        Breakpoint* bp = *it;
        ReadSupportInfo* info = this->mergedBreakpoints.getReadSupport(bp);

        FinalBreakpointInfo finalBreakpoint;
        finalBreakpoint.leftTemplateID = bp->leftTemplateID;
        finalBreakpoint.rightTemplateID = bp->rightTemplateID;

        // calculate positions
        uint32_t splitReadCnt = (info->splitReadSupport + info->clippedReadSupport);
        if (splitReadCnt > 0)
        {
            sort(bp->leftPos.begin(), bp->leftPos.end());
            sort(bp->rightPos.begin(), bp->rightPos.end());
            finalBreakpoint.leftPosition = MID_ELEMENT(bp->leftPos);
            finalBreakpoint.rightPosition = MID_ELEMENT(bp->rightPos);
        }
        else // imprecise cases (only supported by read-pairs)
        {
            finalBreakpoint.leftPosition = bp->minLeftPos;
            finalBreakpoint.rightPosition = bp->maxRightPos;
        }
        finalBreakpoints[bp] = finalBreakpoint;
    }

    return true;
}

bool BreakpointManager::addImpreciseBreakpoints(void)
{   
    // imprecise SVs
    TBreakpointSet* candidateSet = this->pairedEndBreakpoints.getCandidateSet();
    TPosition minSVSize = this->optionManager->getMinSVSize();

    for(auto it = candidateSet->begin(); it != candidateSet->end(); ++it)
    {
        if( this->pairedEndBreakpoints.isBreakpointUsed(*it) == false)
        {
            Breakpoint* newBp = new Breakpoint;
            TPosition leftPos, rightPos, tempPos;
            this->pairedEndBreakpoints.copyBreakpoint(*newBp, *(*it));

            sort(newBp->leftPos.begin(), newBp->leftPos.end());
            sort(newBp->rightPos.begin(), newBp->rightPos.end());
            leftPos = MID_ELEMENT(newBp->leftPos);
            rightPos = MID_ELEMENT(newBp->rightPos);

            /*
            if (newBp->orientation == BreakpointEvidence::ORIENTATION::PROPERLY_ORIENTED)
            {
                // shrink to inside
                leftPos = newBp->minLeftPos;
                rightPos = newBp->maxRightPos;

                // expand to outside : rightPos(-), leftPos(+)
                BreakpointCandidate::setPositionWithAdj(rightPos, leftPos, minSVSize);
            }
            else if (newBp->orientation == BreakpointEvidence::ORIENTATION::SWAPPED)
            {
                // shrink to inside
                leftPos = newBp->maxLeftPos;
                rightPos = newBp->minRightPos;

                // expand to outside : leftPos(-), rightPos(+)
                BreakpointCandidate::setPositionWithAdj(leftPos, rightPos, minSVSize);
            }
            else if (newBp->orientation == BreakpointEvidence::ORIENTATION::INVERSED)
            {
                if (newBp->leftReverseFlag) // <-- <--
                {
                    // shrink to right side & expand to left side
                    leftPos = newBp->maxLeftPos;
                    tempPos = leftPos;
                    BreakpointCandidate::setPositionWithAdj(leftPos, tempPos, minSVSize);

                    rightPos = newBp->maxRightPos;
                    tempPos = rightPos;
                    BreakpointCandidate::setPositionWithAdj(rightPos, tempPos, minSVSize);
                }
                else // -> ->
                {
                    // shrink to left side & expand to right side
                    leftPos = newBp->minLeftPos;
                    tempPos = leftPos;
                    BreakpointCandidate::setPositionWithAdj(tempPos, leftPos, minSVSize);

                    rightPos = newBp->minRightPos;
                    tempPos = rightPos;
                    BreakpointCandidate::setPositionWithAdj(tempPos, rightPos, minSVSize);
                }
            }
            else
                continue;
            */

            // left position
            newBp->leftPos.clear();
            newBp->leftPos.push_back(leftPos);
            newBp->minLeftPos = leftPos;
            newBp->maxLeftPos = leftPos;
            newBp->needLeftIndexUpdate = true;

            // right position
            newBp->rightPos.clear();
            newBp->rightPos.push_back(rightPos);
            newBp->minRightPos = rightPos;
            newBp->maxRightPos = rightPos;
            newBp->needRightIndexUpdate = true;

            // add
            bool isNew = false;
            newBp = this->mergedBreakpoints.updateBreakpoint(newBp, isNew);

            // update read-info
            ReadSupportInfo* newReadInfo;
            newReadInfo = this->mergedBreakpoints.getReadSupport(newBp);
            newReadInfo->pairedEndSupport += newBp->suppReads;
        }        
    }

    return true;
}

bool BreakpointManager::applyFilter(void)
{
    return this->filterByVote();
}

bool BreakpointManager::filterByVote(void)
{
    TBreakpointSet* candidateSet = this->mergedBreakpoints.getCandidateSet();
    auto itBreakpoint = candidateSet->begin(); 

    int minVote = this->optionManager->getMinVote();
    double lowerVoteBound = -1;
    double upperVoteBound = MaxValue<double>::VALUE;

    if (this->optionManager->doReadDepthAnalysis())
    {
        lowerVoteBound = this->optionManager->getAverageReadDepth() * (1-this->optionManager->getVoteBound());
        upperVoteBound = this->optionManager->getAverageReadDepth() * (1+this->optionManager->getVoteBound());
    }
    printTimeMessage("Upper & lower bounds for categorization: " + std::to_string(upperVoteBound) + "," + std::to_string(lowerVoteBound));

    while (itBreakpoint != candidateSet->end())
    {
        Breakpoint* bp = *itBreakpoint;
        ReadSupportInfo* info = this->mergedBreakpoints.getReadSupport(bp);
        FinalBreakpointInfo* finalBreakpoint = &this->finalBreakpoints[bp];
        unsigned voteCount = 0;

        // thresholds
        double seTH = 0.0;
        double peTH = 0.0;
        double reTHHigh = 0.0;
        double reTHMid = 0.0;
        double reTHLow = 0.0;
        if (this->optionManager->useGlobalThreshold() == true)
        {
            seTH = this->optionManager->getMinSplitReadSupport();
            peTH = this->optionManager->getMinPairSupport();
        }
        else
        {
            seTH = this->optionManager->getMinSplitReadSupport();
            peTH = this->optionManager->getMinPairSupport();
        }

        // limit the maximum
        if (info->splitReadSupport > ReadSupportInfo::MAX_VAL)
            info->splitReadSupport = ReadSupportInfo::MAX_VAL;
        if (info->clippedReadSupport > ReadSupportInfo::MAX_VAL)
            info->clippedReadSupport = ReadSupportInfo::MAX_VAL;
        if (info->pairedEndSupport > ReadSupportInfo::MAX_VAL)
            info->pairedEndSupport = ReadSupportInfo::MAX_VAL;

        // vote 1 : split-read + clipped-read
        if (this->optionManager->doClippedReadAnalysis()) // if clipped-read considered
        {
            if ((info->clippedReadSupport+info->splitReadSupport) >= seTH)
                voteCount += 1;
        }
        else
        {
            if (info->splitReadSupport >= seTH)
                voteCount += 1;
        }

        // vote 2 : pairend-read
        if (this->optionManager->doPairedEndAnalysis() && info->pairedEndSupport >= peTH)
            voteCount += 1;

        // vote 3 : read-depth
        if (this->optionManager->doReadDepthAnalysis())
        {   
            if (bp->orientation == BreakpointEvidence::PROPERLY_ORIENTED)
            {
                // potental deletions
                if (info->leftReadDepthDiffScore >= (info->leftReadDepth * this->optionManager->getDDSHigh()) || \
                    info->rightReadDepthDiffScore >= (info->rightReadDepth * this->optionManager->getDDSHigh()) )
                {
                    voteCount += 1;
                }
                else
                {
                    // potential duplications, translocations
                    if (info->leftReadDepthDiffScore >= (info->leftReadDepth * this->optionManager->getDDSMid()) || \
                        info->rightReadDepthDiffScore >= (info->rightReadDepth * this->optionManager->getDDSMid()) )
                    {
                        // this will not be reported as a deletion
                        info->pseudoDeletion = true;
                        voteCount += 1;
                    }
                }
            }
            else if (bp->orientation == BreakpointEvidence::SWAPPED)
            {
                if (info->leftReadDepthDiffScore >= (info->leftReadDepth * this->optionManager->getDDSMid()) || \
                    info->rightReadDepthDiffScore >= (info->rightReadDepth * this->optionManager->getDDSMid()) )
                {
                    voteCount += 1;
                }
            }
            else if (bp->orientation == BreakpointEvidence::INVERSED)
            {
                if (info->leftReadDepthDiffScore >= (info->leftReadDepth * this->optionManager->getDDSLow()) || \
                    info->rightReadDepthDiffScore >= (info->rightReadDepth * this->optionManager->getDDSLow()) )
                {
                    voteCount += 1;
                }
            }
        }
        
        // decision
        bool filterOut = false;
        if (info->avgReadDepth <= lowerVoteBound && voteCount < (minVote-1))
            filterOut = true;
        if (info->avgReadDepth >= upperVoteBound && voteCount < (minVote+1))
            filterOut = true;
        if (info->avgReadDepth > lowerVoteBound && info->avgReadDepth < upperVoteBound && voteCount < minVote)
            filterOut = true;

        if (filterOut == true)
        {
            this->finalBreakpoints.erase(bp);
            itBreakpoint = this->mergedBreakpoints.removeBreakpoint(bp);
        }
        else
        {
            finalBreakpoint->score = info->splitReadSupport + info->clippedReadSupport + info->pairedEndSupport;
            ++itBreakpoint;
        }
    }
    //std::cerr << "after filter : " << candidateSet->size() << "\n";
    return true;
}

bool BreakpointManager::writeBreakpoint()
{
    MergedCandidate* mergedBreakpoint = this->getMergedBreakpoint();
    TBreakpointSet* candidateSet = mergedBreakpoint->getCandidateSet();
    auto itBreakpoint = candidateSet->begin();
    int32_t nID = 1;

    std::ofstream outfile;
    outfile.open("breakpoints.tsv");

    outfile << "leftChr\tleftPos\trightChr\trightPos\t";
    outfile << "orientation\tdepth\tSR\tPE\tCR\tRD\t";
    outfile << "leftStart\tleftEnd\trightStart\trightEnd\n";

    while (itBreakpoint != candidateSet->end())
    {
        Breakpoint* bp = *itBreakpoint;
        ReadSupportInfo* info =  mergedBreakpoint->getReadSupport(bp);
        FinalBreakpointInfo* finalBreakpoint = this->getFinalBreakpointInfo(bp);
        
        outfile << finalBreakpoint->leftTemplateID << "\t";
        outfile << finalBreakpoint->leftPosition << "\t";
        outfile << finalBreakpoint->rightTemplateID << "\t";
        outfile << finalBreakpoint->rightPosition << "\t";

        //PROPERLY_ORIENTED, SWAPPED, INVERSED, NOT_DECIDED
        outfile << bp->orientation << "\t";
        outfile << std::to_string(info->avgReadDepth) << "\t";
        outfile << std::to_string(info->splitReadSupport)<< "\t";
        outfile << std::to_string(info->pairedEndSupport) << "\t";
        outfile << std::to_string(info->clippedReadSupport) << "\t";
        outfile << std::to_string(info->leftReadDepthDiffScore) << "\t";
        outfile << std::to_string(info->rightReadDepthDiffScore) << "\t";

        outfile << bp->minLeftPos << "\t";
        outfile << bp->maxLeftPos << "\t";
        outfile << bp->minRightPos << "\t";
        outfile << bp->maxRightPos << "\n";

        ++itBreakpoint;
    }
    outfile.close();

    return true;
}

void BreakpointManager::init(AlignmentManager& aln)
{ 
    this->optionManager = aln.getOptionManager();
    this->alignmentManager = &aln;

    this->splitReadBreakpoints.setOptionManager(this->optionManager);            
    this->pairedEndBreakpoints.setOptionManager(this->optionManager);   
    this->clippedBreakpoints.setOptionManager(this->optionManager);
    this->readDepthBreakpoints.setOptionManager(this->optionManager);
    this->mergedBreakpoints.setOptionManager(this->optionManager);

    // register objects to alignment manager
    this->alignmentManager->setBreakpointCandidate(&splitReadBreakpoints, &pairedEndBreakpoints, &clippedBreakpoints, &readDepthBreakpoints);
}