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
    printTimeMessage("Total breakpoints: " + std::to_string(this->mergedBreakpoints.getBreakpointCount()));
    if (this->optionManager->doPairedEndAnalysis())
        RUN(result,"From read-pair evidences.", mergePairedEnd());
    printTimeMessage("Total breakpoints: " + std::to_string(this->mergedBreakpoints.getBreakpointCount()));

    if (this->optionManager->doClippedReadAnalysis())
        RUN(result,"From soft-clipped evidences.", mergeClippedRead());
    printTimeMessage("Total breakpoints: " + std::to_string(this->mergedBreakpoints.getBreakpointCount()));

    // 2) add imprecise breakpoints
    if (this->optionManager->doPairedEndAnalysis())
        RUN(result,"Add imprecise breakpoints to the merged set.", addImpreciseBreakpoints());
    printTimeMessage("Total breakpoints: " + std::to_string(this->mergedBreakpoints.getBreakpointCount()));

    // 3) get final positions
    RUN(result,"Get final breakpoints.", findFinalBreakpoints());
    printTimeMessage("Total breakpoints: " + std::to_string(this->mergedBreakpoints.getBreakpointCount()));

    // 4) add read-depth based information (needs final breakpoints)
    if (this->optionManager->doReadDepthAnalysis())
        RUN(result,"Calculate read-depth information", calculateReadDepth());
    
    // 5) sequence featre
    if (this->optionManager->doClippedReadAnalysis())
        RUN(result,"Get sequence information.", getSequenceFeature());

    printTimeMessage("Total breakpoints: " + std::to_string(this->mergedBreakpoints.getBreakpointCount()));
}

bool BreakpointManager::mergeSplitRead(void)
{
    // split-read
    TBreakpointSet* breakpoints = this->splitReadBreakpoints.getCandidateSet();
    auto itBreakpoint = breakpoints->begin();;
    while ( itBreakpoint != breakpoints->end() )
    {
        bool isNew;
        Breakpoint* mergedBp = this->mergedBreakpoints.moveAndUpdateBreakpoint(*itBreakpoint, isNew);
        ReadSupportInfo* mergedReadInfo = this->mergedBreakpoints.getReadSupport(mergedBp);

        mergedReadInfo->splitReadSupport = mergedBp->suppReads;
        mergedReadInfo->pairedEndSupport = 0;
        mergedReadInfo->clippedReadSupport = 0;

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
   
        TBreakpointSet leftMatched, rightMatched, bpEquivalent;
        this->mergedBreakpoints.findMatchedBreakpoint(leftMatched, rightMatched, bpEquivalent, currentBp, true);

        // found matches
        if (bpEquivalent.size() > 0)
        {
            // representative bp
            auto itEq = bpEquivalent.begin();
            Breakpoint* destBp = *(itEq++);
            ReadSupportInfo* mergedReadInfo = this->mergedBreakpoints.getReadSupport(destBp);

            // merge curret paired-end bp without positional info.
            BreakpointCandidate::clearPosInfo(currentBp);
            this->mergedBreakpoints.mergeBreakpoint(destBp, currentBp);
            mergedReadInfo->pairedEndSupport += currentBp->suppReads;

            // merge others that are linked by this paired-end bp
            for (; itEq != bpEquivalent.end(); ++itEq)
                this->mergedBreakpoints.mergeBreakpoint(destBp, *itEq);

            // update index
            this->mergedBreakpoints.updateBreakpointIndex(destBp);

            // marking this BP 
            this->pairedEndBreakpoints.setBreakpointUsed(currentBp, true);
            itBreakpoint = this->pairedEndBreakpoints.removeBreakpoint(currentBp);
        }
        else // This BP will get a second chance in "addImpreciseBreakpoints"
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

        // find matches from merged & paired-end set
        CharString clipSeq = "";
        CharString revCompClipSeq = "";
        int32_t bestMatchScore;

        TFoundPosition foundPositions;
        foundPositions.clear();

        leftMatched.clear();
        rightMatched.clear();
        this->mergedBreakpoints.findMatchedBreakpoint(leftMatched, rightMatched, currentBp, false);
        this->pairedEndBreakpoints.findMatchedBreakpoint(leftMatched, rightMatched, currentBp, false);

        // 1. check pre-defined regions
        if ( rightMatched.size() > 0 || leftMatched.size() > 0)
        {              
            this->clippedBreakpoints.getConsensusSequence(clipSeq, currentBp);
            bestMatchScore = -(this->optionManager->getClippedSeqErrorRate() * length(clipSeq));

            if (isLeftClip == true)
            {
                // right match
                for (auto itBp = rightMatched.begin(); itBp != rightMatched.end(); ++itBp )
                {
                    if ((*itBp)->orientation == BreakpointEvidence::ORIENTATION::PROPERLY_ORIENTED_LARGE)
                    {
                        this->clippedBreakpoints.searchPairRegion( foundPositions, \
                                                                   *itBp, \
                                                                   bestMatchScore, \
                                                                   clipSeq, \
                                                                   ClippedRead::SIDE::LEFT, \
                                                                   false, \
                                                                   false, \
                                                                   BreakpointEvidence::ORIENTATION::PROPERLY_ORIENTED_LARGE);
                    }
                    else if ((*itBp)->orientation == BreakpointEvidence::ORIENTATION::INVERTED)
                    {
                        if (revCompClipSeq == "")
                        {
                            revCompClipSeq = clipSeq;
                            reverseComplement(revCompClipSeq);
                        }
                        this->clippedBreakpoints.searchPairRegion( foundPositions, \
                                                                   *itBp, \
                                                                   bestMatchScore, \
                                                                   revCompClipSeq, \
                                                                   ClippedRead::SIDE::LEFT, \
                                                                   false, \
                                                                   true, \
                                                                   BreakpointEvidence::ORIENTATION::INVERTED);
                    }
                }
            }
            else
            {
                // left match
                for (auto itBp = leftMatched.begin(); itBp != leftMatched.end(); ++itBp )
                {
                    if ((*itBp)->orientation == BreakpointEvidence::ORIENTATION::PROPERLY_ORIENTED_LARGE)
                    {
                        this->clippedBreakpoints.searchPairRegion( foundPositions, \
                                                                   *itBp, \
                                                                   bestMatchScore, \
                                                                   clipSeq, \
                                                                   ClippedRead::SIDE::RIGHT, \
                                                                   false, \
                                                                   false, \
                                                                   BreakpointEvidence::ORIENTATION::PROPERLY_ORIENTED_LARGE);
                    }
                    else if ((*itBp)->orientation == BreakpointEvidence::ORIENTATION::INVERTED)
                    {
                        if (revCompClipSeq == "")
                        {
                            revCompClipSeq = clipSeq;
                            reverseComplement(revCompClipSeq);
                        }
                        this->clippedBreakpoints.searchPairRegion( foundPositions, \
                                                                   *itBp, \
                                                                   bestMatchScore, \
                                                                   revCompClipSeq, \
                                                                   ClippedRead::SIDE::RIGHT, \
                                                                   false, \
                                                                   true, \
                                                                   BreakpointEvidence::ORIENTATION::INVERTED);
                    }
                }        
            }
        }
        // 2. search for twilight zone  (potential deletion)
        if (foundPositions.size() == 0)
        {
            if (clipSeq == "")
            {
                this->clippedBreakpoints.getConsensusSequence(clipSeq, currentBp);
                bestMatchScore = -(this->optionManager->getClippedSeqErrorRate() * length(clipSeq));
            }

            if (isLeftClip == true)
                this->clippedBreakpoints.searchTwilightZone(foundPositions, currentBp, bestMatchScore, clipSeq, ClippedRead::SIDE::LEFT, false, BreakpointEvidence::ORIENTATION::PROPERLY_ORIENTED_LARGE);
            else
                this->clippedBreakpoints.searchTwilightZone(foundPositions, currentBp, bestMatchScore, clipSeq, ClippedRead::SIDE::RIGHT, false, BreakpointEvidence::ORIENTATION::PROPERLY_ORIENTED_LARGE);
        }

        // 3. swap and search
        if (foundPositions.size() == 0)
        {
            // swap left and right
            std::swap(currentBp->leftTemplateID, currentBp->rightTemplateID);
            std::swap(currentBp->minLeftPos, currentBp->minRightPos);
            std::swap(currentBp->maxLeftPos, currentBp->maxRightPos);
            std::swap(currentBp->leftPos, currentBp->rightPos);

            leftMatched.clear();
            rightMatched.clear();
            this->mergedBreakpoints.findMatchedBreakpoint(leftMatched, rightMatched, currentBp, false);        
            this->pairedEndBreakpoints.findMatchedBreakpoint(leftMatched, rightMatched, currentBp, false);
            
            // left clip & left match (swap)
            bool foundInvBySwap = false;
            if (isLeftClip == true && leftMatched.size() > 0)
            {
                for (auto itBp = leftMatched.begin(); itBp != leftMatched.end(); ++itBp )
                {
                    if ((*itBp)->orientation == BreakpointEvidence::ORIENTATION::SWAPPED)
                    {
                        this->clippedBreakpoints.searchPairRegion( foundPositions, \
                                                                   *itBp, \
                                                                   bestMatchScore, \
                                                                   clipSeq, \
                                                                   ClippedRead::SIDE::RIGHT, \
                                                                   false, \
                                                                   false, \
                                                                   BreakpointEvidence::ORIENTATION::SWAPPED);
                    }
                    else if ((*itBp)->orientation == BreakpointEvidence::ORIENTATION::INVERTED)
                    {
                        if (revCompClipSeq == "")
                        {
                            revCompClipSeq = clipSeq;
                            reverseComplement(revCompClipSeq);
                        }
                        foundInvBySwap |= this->clippedBreakpoints.searchPairRegion( foundPositions, \
                                                                   *itBp, \
                                                                   bestMatchScore, \
                                                                   revCompClipSeq, \
                                                                   ClippedRead::SIDE::RIGHT, \
                                                                   false, \
                                                                   true, \
                                                                   BreakpointEvidence::ORIENTATION::INVERTED);
                    }
                }
            }
            else if (isLeftClip == false && rightMatched.size() > 0)
            {
                for (auto itBp = rightMatched.begin(); itBp != rightMatched.end(); ++itBp )
                {
                    if ((*itBp)->orientation == BreakpointEvidence::ORIENTATION::SWAPPED)
                    {
                        this->clippedBreakpoints.searchPairRegion( foundPositions, \
                                                                   *itBp, \
                                                                   bestMatchScore, \
                                                                   clipSeq, \
                                                                   ClippedRead::SIDE::LEFT, \
                                                                   false, \
                                                                   false, \
                                                                   BreakpointEvidence::ORIENTATION::SWAPPED);
                    }
                    else if ((*itBp)->orientation == BreakpointEvidence::ORIENTATION::INVERTED)
                    {
                        if (revCompClipSeq == "")
                        {
                            revCompClipSeq = clipSeq;
                            reverseComplement(revCompClipSeq);
                        }
                        foundInvBySwap |= this->clippedBreakpoints.searchPairRegion( foundPositions, \
                                                                   *itBp, \
                                                                   bestMatchScore, \
                                                                   revCompClipSeq, \
                                                                   ClippedRead::SIDE::LEFT, \
                                                                   false, \
                                                                   true, \
                                                                   BreakpointEvidence::ORIENTATION::INVERTED);
                    }
                }
            }

            // found inversion by swapping
            if (foundInvBySwap)
                isLeftClip = !isLeftClip;
        }

        // found positions
        if (foundPositions.size() > 0)
            addNewPositionsByClippedSequence(foundPositions, currentBp, isLeftClip);

        itBreakpoint = this->clippedBreakpoints.removeBreakpoint(currentBp);
    }
    return true;
}

void BreakpointManager::addNewPositionsByClippedSequence(TFoundPosition& foundPositions, Breakpoint* clippedBp, bool isLeftClip)
{
    // add new positions
    for (auto itFound = foundPositions.begin(); itFound != foundPositions.end(); ++itFound)
    {
        Breakpoint* matchedBp = itFound->matchedBp;

        // copy
        Breakpoint* newBp = new Breakpoint;
        BreakpointCandidate::copyBreakpoint(*newBp, *clippedBp);
        //BreakpointCandidate::moveBreakpoint(*newBp, *clippedBp);
        newBp->orientation = itFound->orientation;
       
        TPosition newPos = BreakpointEvidence::INVALID_POS, svSize = 0;
        if (itFound->orientation == BreakpointEvidence::ORIENTATION::SWAPPED)
        {
            if (isLeftClip)
            {
                newPos = itFound->sequenceSegment.endPos - 1;

                // fill information
                if (matchedBp->rightTemplateID == BreakpointEvidence::NOVEL_TEMPLATE)
                    newBp->rightTemplateID = newBp->rightTemplateID; // by twilight zone
                else
                    newBp->rightTemplateID = matchedBp->rightTemplateID; // by paired-end region

                // update positions
                newBp->rightPos.clear();
                newBp->rightPos.push_back(newPos);            
                BreakpointCandidate::updateRightMinMaxPos(newBp);

                // update depth
                if (this->optionManager->doReadDepthAnalysis())
                {
                    this->readDepthBreakpoints.addUniformDepth(newBp->rightTemplateID, \
                                                               newPos, \
                                                               newBp->clippedConsensusSequenceSize,
                                                               newBp->suppReads);
                }
            }
            else
            {
                newPos = itFound->sequenceSegment.beginPos;

                // fill information
                if (matchedBp->leftTemplateID == BreakpointEvidence::NOVEL_TEMPLATE)
                    newBp->leftTemplateID = newBp->rightTemplateID; // by twilight zone
                else
                    newBp->leftTemplateID = matchedBp->leftTemplateID; // by paired-end region

                // update positions
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

        }
        else
        {
            if (isLeftClip)
            {
                if (newBp->orientation == BreakpointEvidence::ORIENTATION::INVERTED)
                {
                    newPos = itFound->sequenceSegment.beginPos;
                    newBp->leftReverseFlag = !newBp->rightReverseFlag;

                }
                else if (newBp->orientation == BreakpointEvidence::ORIENTATION::PROPERLY_ORIENTED_LARGE)
                {
                    newPos = itFound->sequenceSegment.endPos;
                }

                if (newPos != BreakpointEvidence::INVALID_POS)
                {
                    // fill information
                    if (matchedBp->leftTemplateID == BreakpointEvidence::NOVEL_TEMPLATE)
                        newBp->leftTemplateID = newBp->rightTemplateID; // by twilight zone
                    else
                        newBp->leftTemplateID = matchedBp->leftTemplateID; // by paired-end region

                    // update positions
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
            }
            else
            {
                if (newBp->orientation== BreakpointEvidence::ORIENTATION::INVERTED)
                {
                    // original strand => outside inversion
                    newPos = itFound->sequenceSegment.endPos - 1;
                    newBp->rightReverseFlag = !newBp->leftReverseFlag;
                }
                else if (newBp->orientation== BreakpointEvidence::ORIENTATION::PROPERLY_ORIENTED_LARGE)
                {
                    newPos = itFound->sequenceSegment.beginPos - 1;
                }

                if (newPos != BreakpointEvidence::INVALID_POS)
                {
                    // fill information
                    if (matchedBp->rightTemplateID == BreakpointEvidence::NOVEL_TEMPLATE)
                        newBp->rightTemplateID = newBp->leftTemplateID; // by twilight zone
                    else
                        newBp->rightTemplateID = matchedBp->rightTemplateID; // by paired-end region

                    // update positions
                    newBp->rightPos.clear();
                    newBp->rightPos.push_back(newPos);            
                    BreakpointCandidate::updateRightMinMaxPos(newBp);
                
                    // update depth
                    if (this->optionManager->doReadDepthAnalysis())
                    {
                        this->readDepthBreakpoints.addUniformDepth(newBp->rightTemplateID, \
                                                                   newPos, \
                                                                   newBp->clippedConsensusSequenceSize,
                                                                   newBp->suppReads);
                    }
                }
            }

        }

        // update
        if (newPos != BreakpointEvidence::INVALID_POS)
        {
            // add to the merged set
            bool isNewBp = false;
            newBp->needLeftIndexUpdate = true;
            newBp->needRightIndexUpdate = true;
            newBp = this->mergedBreakpoints.updateBreakpoint(newBp, true, isNewBp);

            // update supporting read information
            ReadSupportInfo* newReadInfo = this->mergedBreakpoints.getReadSupport(newBp);
            newReadInfo->clippedReadSupport += clippedBp->suppReads;

            // new bp : based on read-pairs or twilight zone (one-side has to be the novel_template)
            if (isNewBp == true)
            {
                // based on read-pairs
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
}

bool BreakpointManager::calculateReadDepth()
{
    TBreakpointSet* candidateSet = this->mergedBreakpoints.getCandidateSet();
    std::map<TTemplateID, unsigned> svCntByTemplate;
    for(auto it = candidateSet->begin(); it != candidateSet->end(); ++it)
    {
        Breakpoint* bp = *it;
        ReadSupportInfo* info = this->mergedBreakpoints.getReadSupport(bp);
        FinalBreakpointInfo* finalBreakpoint = &this->finalBreakpoints[bp];

        int32_t breakpointSize = abs(finalBreakpoint->rightPosition - finalBreakpoint->leftPosition) + 1;
        int32_t windowSize = this->optionManager->getReadDepthWindowSize();

        // restrict the search region's size
        //if (breakpointSize > windowSize)
        //    breakpointSize = windowSize;
        if (windowSize > breakpointSize)
            windowSize = breakpointSize;

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
                                                      windowSize);

        if (info->leftReadDepthDiffScore > info->rightReadDepthDiffScore)
            info->leftReadDepthSelected = true;
        else
            info->leftReadDepthSelected = false;
        info->avgReadDepth = (info->leftReadDepth + info->rightReadDepth) / 2.0;

        // count
        if (svCntByTemplate.find(bp->leftTemplateID) == svCntByTemplate.end())
          svCntByTemplate[bp->leftTemplateID] = 0;
        ++svCntByTemplate[bp->leftTemplateID];
    }
    this->readDepthBreakpoints.calculateReadDepthStat(svCntByTemplate, candidateSet->size());
}

void BreakpointManager::getNTCount(CharString& sequence, unsigned& a, unsigned& t, unsigned& g, unsigned& c)
{
    for(unsigned i=0; i < length(sequence); ++i)
    {
        switch (sequence[i])
        {
            case 'A':
            case 'a':
                ++a;
                break;
            case 'T':
            case 't':
                ++t;
                break;
            case 'G':
            case 'g':
                ++g;
                break;
            case 'C':
            case 'c':
                ++c;
                break;
        }
    }
}


bool BreakpointManager::getSequenceFeature(void)
{
    //double GCContent = 0.0;
    //double sequenceComplexity = 0.0;
    unsigned kmerSize = 5;
    unsigned windowSize = this->optionManager->getReadDepthWindowSize();
    windowSize = 1000;
    TBreakpointSet* candidateSet = this->mergedBreakpoints.getCandidateSet();
    for(auto it = candidateSet->begin(); it != candidateSet->end(); ++it)
    {
        FinalBreakpointInfo& finalBreakpoint = this->finalBreakpoints[*it];
        CharString ref;
        unsigned gc = 0, a = 0, t = 0, g = 0, c = 0;
        
        TPosition start = 0;
        TPosition end = BreakpointEvidence::INVALID_POS;
        TPosition refSize = 0;

        // left
        end = finalBreakpoint.leftPosition;
        start = std::max(end-windowSize, (unsigned) 0);
        this->clippedBreakpoints.getReferenceSequence(ref, finalBreakpoint.leftTemplateID, start, end);
        refSize += length(ref);
        this->getNTCount(ref, a, t, g, c);

        // right
        start = finalBreakpoint.rightPosition;
        end = std::min(start + windowSize, this->readDepthBreakpoints.getRefSize(finalBreakpoint.rightTemplateID));
        this->clippedBreakpoints.getReferenceSequence(ref, finalBreakpoint.rightTemplateID, start, end);
        refSize += length(ref);
        this->getNTCount(ref, a, t, g, c);

        finalBreakpoint.gcContent = (double) (g+c) / refSize;
        double entropy = 0;
        entropy  = -(((double)a / refSize) * log2((double)a / refSize));
        entropy += -(((double)t / refSize) * log2((double)t / refSize));
        entropy += -(((double)g / refSize) * log2((double)g / refSize));
        entropy += -(((double)c / refSize) * log2((double)c / refSize));
        finalBreakpoint.sequenceComplexity = entropy;
    }

    return true;
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
        finalBreakpoint.isLeftReverse = bp->leftReverseFlag;
        finalBreakpoint.rightTemplateID = bp->rightTemplateID;
        finalBreakpoint.isRightReverse = bp->rightReverseFlag;

        // calculate positions
        uint32_t splitReadCnt = (info->splitReadSupport + info->clippedReadSupport);
        if (splitReadCnt > 0)
        {
            sort(bp->leftPos.begin(), bp->leftPos.end());
            sort(bp->rightPos.begin(), bp->rightPos.end());
            finalBreakpoint.leftPosition = MID_ELEMENT(bp->leftPos);
            finalBreakpoint.rightPosition = MID_ELEMENT(bp->rightPos);
            finalBreakpoint.imprecise = false;
        }
        else // imprecise cases (only supported by read-pairs)
        {
            finalBreakpoint.leftPosition = bp->minLeftPos;
            finalBreakpoint.rightPosition = bp->maxRightPos;
            finalBreakpoint.imprecise = true;
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
    //double cutoff = this->readDepthBreakpoints.getIQR() * this->optionManager->getDepthOutlier();
    //printTimeMessage("Depth outlier: " + std::to_string(cutoff));

    for(auto itBp = candidateSet->begin(); itBp != candidateSet->end(); ++itBp)
    {
        if( this->pairedEndBreakpoints.isBreakpointUsed(*itBp) == false)
        {
            ReadSupportInfo* info = this->mergedBreakpoints.getReadSupport(*itBp);
            double normalizedScore = info->pairedEndSupport;

            //normalizedScore += BreakpointCandidate::PREVENT_DIV_BY_ZERO();
            //normalizedScore /= (info->avgReadDepth + BreakpointCandidate::PREVENT_DIV_BY_ZERO());
            //normalizedScore *= this->readDepthBreakpoints.getMedianDepth();
            //if (normalizedScore >= (double) cutoff)

            //if ( info->avgReadDepth < cutoff)
            {
                //TBreakpointSet leftMatched, rightMatched, bpEquivalent;
                //this->mergedBreakpoints.findMatchedBreakpoint(leftMatched, rightMatched, bpEquivalent, *itBp, true);
                //if (leftMatched.size() == 0 && rightMatched.size() == 0)

                Breakpoint* newBp = new Breakpoint;
                TPosition leftPos, rightPos, tempPos;
                this->pairedEndBreakpoints.copyBreakpoint(*newBp, *(*itBp));

                sort(newBp->leftPos.begin(), newBp->leftPos.end());
                sort(newBp->rightPos.begin(), newBp->rightPos.end());
                leftPos = MID_ELEMENT(newBp->leftPos);
                rightPos = MID_ELEMENT(newBp->rightPos);

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
                newBp = this->mergedBreakpoints.updateBreakpoint(newBp, true, isNew);

                // update read-info
                ReadSupportInfo* newReadInfo;
                newReadInfo = this->mergedBreakpoints.getReadSupport(newBp);
                newReadInfo->pairedEndSupport += newBp->suppReads;
            }
            /*
            else
            {
                for (auto it = bpEquivalent.begin(); it != bpEquivalent.end(); ++it)
                {
                    ReadSupportInfo* newReadInfo = this->mergedBreakpoints.getReadSupport(*it);
                    newReadInfo->pairedEndSupport += (*itBp)->suppReads;
                }
            }
            */
        }        
    }

    return true;
}

bool BreakpointManager::applyFilter(void)
{
  //if (this->optionManager->getUseRankAggregation())
  //  this->applyNormalization();
  this->filterByEvidenceSumAndVote();
}

bool BreakpointManager::priByEvidenceSum(void)
{
    TBreakpointSet* candidateSet = this->mergedBreakpoints.getCandidateSet();
    auto itBreakpoint = candidateSet->begin();
    int cutoff = this->optionManager->getCutoff();

    while (itBreakpoint != candidateSet->end())
    {
        if (this->finalBreakpoints[*itBreakpoint].filtered == false)
        {
            FinalBreakpointInfo* finalBreakpointInfo = &this->finalBreakpoints[*itBreakpoint];
            ReadSupportInfo* info = this->mergedBreakpoints.getReadSupport(*itBreakpoint);
            finalBreakpointInfo->score =  (info->splitReadSupport + info->pairedEndSupport + info->clippedReadSupport);
        }
        ++itBreakpoint;
    }

    return true;
}

bool BreakpointManager::applyNormalization(void)
{
    TBreakpointSet* candidateSet = this->mergedBreakpoints.getCandidateSet();
    auto itBreakpoint = candidateSet->begin();
    int cutoff = this->optionManager->getCutoff();

    unsigned filteredBreakpointCount = 0;
    while (itBreakpoint != candidateSet->end())
    {
        FinalBreakpointInfo* finalBreakpointInfo = &this->finalBreakpoints[*itBreakpoint];
        ReadSupportInfo* info = this->mergedBreakpoints.getReadSupport(*itBreakpoint);

        //se = (se + cc) / (rd + cc) * 30.0
        info->splitReadSupport += BreakpointCandidate::PREVENT_DIV_BY_ZERO();
        info->pairedEndSupport += BreakpointCandidate::PREVENT_DIV_BY_ZERO();
        info->clippedReadSupport += BreakpointCandidate::PREVENT_DIV_BY_ZERO();
        info->leftReadDepthDiffScore += BreakpointCandidate::PREVENT_DIV_BY_ZERO();
        info->rightReadDepthDiffScore += BreakpointCandidate::PREVENT_DIV_BY_ZERO();

        info->splitReadSupport /= (info->avgReadDepth + BreakpointCandidate::PREVENT_DIV_BY_ZERO());
        info->pairedEndSupport /= (info->avgReadDepth + BreakpointCandidate::PREVENT_DIV_BY_ZERO());
        info->clippedReadSupport /= (info->avgReadDepth + BreakpointCandidate::PREVENT_DIV_BY_ZERO());
        info->leftReadDepthDiffScore /= (info->avgReadDepth + BreakpointCandidate::PREVENT_DIV_BY_ZERO());
        info->rightReadDepthDiffScore /= (info->avgReadDepth + BreakpointCandidate::PREVENT_DIV_BY_ZERO());

        info->splitReadSupport *= this->readDepthBreakpoints.getMedianDepth();
        info->pairedEndSupport *= this->readDepthBreakpoints.getMedianDepth();
        info->clippedReadSupport *= this->readDepthBreakpoints.getMedianDepth();
        info->leftReadDepthDiffScore *= this->readDepthBreakpoints.getMedianDepth();
        info->rightReadDepthDiffScore *= this->readDepthBreakpoints.getMedianDepth();

        ++itBreakpoint;
    }

    return true;
}

bool BreakpointManager::filterByEvidenceSumAndVote(void)
{
    double seTH = this->optionManager->getMinSplitReadSupport();
    double peTH = this->optionManager->getMinPairSupport();
    double reTH = this->readDepthBreakpoints.getReTH();
    unsigned cutoff = this->optionManager->getCutoff();
    unsigned minVote = this->optionManager->getMinVote();
    if (minVote < 0)
    {
        minVote = 1;
        if (this->optionManager->doPairedEndAnalysis())
            ++minVote;
        if (this->optionManager->doReadDepthAnalysis())
            ++minVote;
        this->optionManager->setMinVote(minVote);
    }

    TBreakpointSet* candidateSet = this->mergedBreakpoints.getCandidateSet();
    auto itBreakpoint = candidateSet->begin();
    unsigned filteredBreakpointCount = 0;
    while (itBreakpoint != candidateSet->end())
    {
        FinalBreakpointInfo* finalBreakpointInfo = &this->finalBreakpoints[*itBreakpoint];
        ReadSupportInfo* info = this->mergedBreakpoints.getReadSupport(*itBreakpoint);

        // Evidence sum
        double score = (info->splitReadSupport + info->pairedEndSupport + info->clippedReadSupport);
        finalBreakpointInfo->score = score;

        // Voting
        unsigned vote = 0;
        if ((info->splitReadSupport + info->clippedReadSupport) >= seTH)
            vote += 1;
        if (this->optionManager->doPairedEndAnalysis() && info->pairedEndSupport >= peTH)
            vote += 1;
        if (this->optionManager->doReadDepthAnalysis())
        {
            double re = 0.0;
            if (info->leftReadDepthDiffScore > info->rightReadDepthDiffScore)
                re = info->leftReadDepthDiffScore;
            else
                re = info->rightReadDepthDiffScore;
            if (re >= reTH)
                vote += 1;
        }
        finalBreakpointInfo->vote = vote;

        //std::cerr << score << "\t" << cutoff << "\t" << minVote << "\t" << vote << "\n";
        if (score >= (double)cutoff || vote >= minVote)
        {
          finalBreakpointInfo->filtered = false;
        }
        else
        {
          ++filteredBreakpointCount;
          finalBreakpointInfo->filtered = true;
        }
        ++itBreakpoint;
    }
    this->mergedBreakpoints.setFilteredBreakpointCount(filteredBreakpointCount);

    return true;
}

bool BreakpointManager::rescueByCombinedEvidence(void)
{
    double seTH = this->optionManager->getMinSplitReadSupport();
    double peTH = this->optionManager->getMinPairSupport();
    double reTH = this->readDepthBreakpoints.getReTH();
   
    // Min. vote for rescue
    int minVote = this->optionManager->getMinVote();
    if (minVote < 0)
    {
        minVote = 1;
        if (this->optionManager->doPairedEndAnalysis())
            ++minVote;
        if (this->optionManager->doReadDepthAnalysis())
            ++minVote;
        this->optionManager->setMinVote(minVote);
    }

    // for all breakpoints
    unsigned filteredBreakpointCount = this->getMergedBreakpoint()->getFilteredBreakpointCount();
    TBreakpointSet* candidateSet = this->mergedBreakpoints.getCandidateSet();
    auto itBreakpoint = candidateSet->begin();   
    while (itBreakpoint != candidateSet->end())
    {
        FinalBreakpointInfo* finalBreakpointInfo = &this->finalBreakpoints[*itBreakpoint];
        ReadSupportInfo* info = this->mergedBreakpoints.getReadSupport(*itBreakpoint);

        // vote
        unsigned vote = 0;
        if ((info->splitReadSupport + info->clippedReadSupport) >= seTH)
            vote += 1;
        if (this->optionManager->doPairedEndAnalysis() && info->pairedEndSupport >= peTH)
            vote += 1;
        if (this->optionManager->doReadDepthAnalysis())
        {
            double re = 0.0;
            if (info->leftReadDepthDiffScore > info->rightReadDepthDiffScore)
                re = info->leftReadDepthDiffScore;
            else
                re = info->rightReadDepthDiffScore;
            if (re >= reTH)
                vote += 1;
        }
        finalBreakpointInfo->vote = vote;

        // rescue
        if (vote >= minVote && finalBreakpointInfo->filtered == true)
        {
            --filteredBreakpointCount;
            finalBreakpointInfo->filtered = false;
        }
        ++itBreakpoint;
    }
    this->mergedBreakpoints.setFilteredBreakpointCount(filteredBreakpointCount);
}

void BreakpointManager::writeBreakpoint()
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