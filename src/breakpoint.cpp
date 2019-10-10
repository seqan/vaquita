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
    printTimeMessage("Total breakpoints: " + std::to_string(this->getMergedBreakpoint()->getBreakpointCount()));
    if (this->optionManager->doPairedEndAnalysis())
    {
        RUN(result,"From read-pair evidences.", mergePairedEnd());
        printTimeMessage("Total breakpoints: " + std::to_string(this->getMergedBreakpoint()->getBreakpointCount()));
    }

    if (this->optionManager->doClippedReadAnalysis())
    {
        RUN(result,"From soft-clipped evidences.", mergeClippedRead());
        printTimeMessage("Total breakpoints: " + std::to_string(this->getMergedBreakpoint()->getBreakpointCount()));
    }

    // 2) add imprecise breakpoints
    if (this->optionManager->doPairedEndAnalysis())
    {
        RUN(result,"Add imprecise breakpoints to the merged set.", addImpreciseBreakpoints());
        printTimeMessage("Total breakpoints: " + std::to_string(this->getMergedBreakpoint()->getBreakpointCount()));
    }

    // 3) get final positions
    RUN(result,"Get final breakpoints.", findFinalBreakpoints());
    printTimeMessage("Total breakpoints: " + std::to_string(this->getMergedBreakpoint()->getBreakpointCount()));

    // 4) add read-depth based information (needs final breakpoints)
    if (this->optionManager->doReadDepthAnalysis())
        RUN(result,"Calculate read-depth information", calculateReadDepth());

    // 5) sequence featre
    if (this->optionManager->doClippedReadAnalysis())
        RUN(result,"Get sequence information.", getSequenceFeature());

    printTimeMessage("Total breakpoints: " + std::to_string(this->getMergedBreakpoint()->getBreakpointCount()));

    return true;
}

bool BreakpointManager::mergeSplitRead(void)
{
    // split-read
    TBreakpointSet* breakpoints = this->getSplitRead()->getCandidateSet();
    auto itBreakpoint = breakpoints->begin();;
    while ( itBreakpoint != breakpoints->end() )
    {
        bool isNew;
        Breakpoint* mergedBp = this->getMergedBreakpoint()->moveAndUpdateBreakpoint(*itBreakpoint, isNew);
        ReadSupportInfo* mergedReadInfo = this->getMergedBreakpoint()->initReadSupport(mergedBp);
        mergedReadInfo->splitReadSupport = mergedBp->suppReads;

        itBreakpoint = this->getSplitRead()->removeBreakpoint(*itBreakpoint);
    }

    return true;
}

bool BreakpointManager::mergePairedEnd(void)
{
    // paired-end
    TBreakpointSet* breakpoints = this->getPairedEndRead()->getCandidateSet();
    auto itBreakpoint = breakpoints->begin();;
    while ( itBreakpoint != breakpoints->end() )
    {
        Breakpoint* currentBp = *itBreakpoint;

        TBreakpointSet leftMatched, rightMatched, bpEquivalent;
        this->getMergedBreakpoint()->findBreakpoint(leftMatched, rightMatched, bpEquivalent, currentBp);

        // found matches
        if (bpEquivalent.size() > 0)
        {
            // representative bp
            auto itEq = bpEquivalent.begin();
            Breakpoint* destBp = *(itEq++);
            ReadSupportInfo* mergedReadInfo = this->getMergedBreakpoint()->getReadSupport(destBp);

            // merge curret paired-end bp without positional info.
            destBp->suppReads += currentBp->suppReads;
            mergedReadInfo->pairedEndSupport += currentBp->suppReads;

            // merge others that are linked by this paired-end bp
            for (; itEq != bpEquivalent.end(); ++itEq)
                this->getMergedBreakpoint()->mergeBreakpoint(destBp, *itEq);

            // update index
            this->getMergedBreakpoint()->updateBreakpointIndex(destBp);

            // marking this BP
            itBreakpoint = this->getPairedEndRead()->removeBreakpoint(currentBp);
        }
        else // This BP will get a second chance in "addImpreciseBreakpoints"
            ++itBreakpoint;
    }

    return true;
}

bool BreakpointManager::mergeClippedRead(void)
{
    // 1) clipped-read
    TBreakpointSet leftMatched, rightMatched, bothMatched;
    TBreakpointSet* breakpoints = this->getClippedRead()->getCandidateSet();
    auto itBreakpoint= breakpoints->begin();
    while ( itBreakpoint != breakpoints->end() )
    {
        Breakpoint* currentBp = *itBreakpoint;
        uint32_t clippedReadsCount = currentBp->suppReads;
        bool isLeftClip = (currentBp->leftTemplateID == BreakpointEvidence::NOVEL_TEMPLATE);

        // find matches from merged & paired-end set
        seqan::CharString clipSeq = "";
        seqan::CharString revCompClipSeq = "";
        int32_t bestMatchScore;

        TFoundPosition foundPositions;
        foundPositions.clear();

        leftMatched.clear();
        rightMatched.clear();
        bothMatched.clear();
        this->getMergedBreakpoint()->findBreakpoint(leftMatched, rightMatched, bothMatched, currentBp);
        this->getPairedEndRead()->findBreakpoint(leftMatched, rightMatched, bothMatched, currentBp);

        // 1. check pre-defined regions
        if ( leftMatched.size() > 0 || rightMatched.size() > 0)
        {
            this->getClippedRead()->getConsensusSequence(clipSeq, currentBp);
            bestMatchScore = -(this->optionManager->getClippedSeqErrorRate() * length(clipSeq));

            if (isLeftClip == true)
            {
                // right match
                for (auto itBp = rightMatched.begin(); itBp != rightMatched.end(); ++itBp )
                {
                    if ((*itBp)->orientation == BreakpointEvidence::ORIENTATION::PROPERLY_ORIENTED_LARGE)
                    {
                        this->getClippedRead()->searchPairRegion( foundPositions, \
                                                                   *itBp, \
                                                                   bestMatchScore, \
                                                                   clipSeq, \
                                                                   BreakpointEvidence::SIDE::LEFT, \
                                                                   false, \
                                                                   false, \
                                                                   BreakpointEvidence::ORIENTATION::PROPERLY_ORIENTED_LARGE);
                    }
                    else if ((*itBp)->orientation == BreakpointEvidence::ORIENTATION::INVERTED)
                    {
                        if (revCompClipSeq == "")
                        {
                            revCompClipSeq = clipSeq;
                            seqan::reverseComplement(revCompClipSeq);
                        }
                        this->getClippedRead()->searchPairRegion( foundPositions, \
                                                                   *itBp, \
                                                                   bestMatchScore, \
                                                                   revCompClipSeq, \
                                                                   BreakpointEvidence::SIDE::LEFT, \
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
                        this->getClippedRead()->searchPairRegion( foundPositions, \
                                                                   *itBp, \
                                                                   bestMatchScore, \
                                                                   clipSeq, \
                                                                   BreakpointEvidence::SIDE::RIGHT, \
                                                                   false, \
                                                                   false, \
                                                                   BreakpointEvidence::ORIENTATION::PROPERLY_ORIENTED_LARGE);
                    }
                    else if ((*itBp)->orientation == BreakpointEvidence::ORIENTATION::INVERTED)
                    {
                        if (revCompClipSeq == "")
                        {
                            revCompClipSeq = clipSeq;
                            seqan::reverseComplement(revCompClipSeq);
                        }
                        this->getClippedRead()->searchPairRegion( foundPositions, \
                                                                   *itBp, \
                                                                   bestMatchScore, \
                                                                   revCompClipSeq, \
                                                                   BreakpointEvidence::SIDE::RIGHT, \
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
                this->getClippedRead()->getConsensusSequence(clipSeq, currentBp);
                bestMatchScore = -(this->optionManager->getClippedSeqErrorRate() * length(clipSeq));
            }

            if (isLeftClip == true)
            {
                this->getClippedRead()->searchTwilightZone(foundPositions, \
                                                           currentBp, \
                                                           bestMatchScore, \
                                                           clipSeq, \
                                                           BreakpointEvidence::SIDE::LEFT, \
                                                           false,
                                                           BreakpointEvidence::ORIENTATION::PROPERLY_ORIENTED_LARGE);
            }
            else
            {
                this->getClippedRead()->searchTwilightZone(foundPositions, \
                                                           currentBp, \
                                                           bestMatchScore, \
                                                           clipSeq, \
                                                           BreakpointEvidence::SIDE::RIGHT, \
                                                           false, \
                                                           BreakpointEvidence::ORIENTATION::PROPERLY_ORIENTED_LARGE);
            }
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
            bothMatched.clear();
            this->getMergedBreakpoint()->findBreakpoint(leftMatched, rightMatched, bothMatched, currentBp);
            this->getPairedEndRead()->findBreakpoint(leftMatched, rightMatched, bothMatched, currentBp);

            // left clip & left match (swap)
            bool foundInvBySwap = false;
            if (isLeftClip == true && leftMatched.size() > 0)
            {
                for (auto itBp = leftMatched.begin(); itBp != leftMatched.end(); ++itBp )
                {
                    if ((*itBp)->orientation == BreakpointEvidence::ORIENTATION::SWAPPED)
                    {
                        this->getClippedRead()->searchPairRegion( foundPositions, \
                                                                   *itBp, \
                                                                   bestMatchScore, \
                                                                   clipSeq, \
                                                                   BreakpointEvidence::SIDE::RIGHT, \
                                                                   false, \
                                                                   false, \
                                                                   BreakpointEvidence::ORIENTATION::SWAPPED);
                    }
                    else if ((*itBp)->orientation == BreakpointEvidence::ORIENTATION::INVERTED)
                    {
                        if (revCompClipSeq == "")
                        {
                            revCompClipSeq = clipSeq;
                            seqan::reverseComplement(revCompClipSeq);
                        }
                        foundInvBySwap |= this->getClippedRead()->searchPairRegion( foundPositions, \
                                                                   *itBp, \
                                                                   bestMatchScore, \
                                                                   revCompClipSeq, \
                                                                   BreakpointEvidence::SIDE::RIGHT, \
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
                        this->getClippedRead()->searchPairRegion( foundPositions, \
                                                                   *itBp, \
                                                                   bestMatchScore, \
                                                                   clipSeq, \
                                                                   BreakpointEvidence::SIDE::LEFT, \
                                                                   false, \
                                                                   false, \
                                                                   BreakpointEvidence::ORIENTATION::SWAPPED);
                    }
                    else if ((*itBp)->orientation == BreakpointEvidence::ORIENTATION::INVERTED)
                    {
                        if (revCompClipSeq == "")
                        {
                            revCompClipSeq = clipSeq;
                            seqan::reverseComplement(revCompClipSeq);
                        }
                        foundInvBySwap |= this->getClippedRead()->searchPairRegion( foundPositions, \
                                                                   *itBp, \
                                                                   bestMatchScore, \
                                                                   revCompClipSeq, \
                                                                   BreakpointEvidence::SIDE::LEFT, \
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

        itBreakpoint = this->getClippedRead()->removeBreakpoint(currentBp);
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
                if (itFound->matchedRightTemplateID == BreakpointEvidence::NOVEL_TEMPLATE)
                    newBp->rightTemplateID = newBp->rightTemplateID; // by twilight zone
                else
                    newBp->rightTemplateID = itFound->matchedRightTemplateID; // by paired-end region

                // update positions
                newBp->rightPos.clear();
                newBp->rightPos.push_back(newPos);
                BreakpointCandidate::updateRightMinMaxPos(newBp);

                // update depth
                if (this->optionManager->doReadDepthAnalysis())
                {
                    this->getReadDepth()->addUniformDepth(newBp->rightTemplateID, \
                                                               newPos, \
                                                               newBp->clippedConsensusSequenceSize,
                                                               newBp->suppReads);
                }
            }
            else
            {
                newPos = itFound->sequenceSegment.beginPos;

                // fill information
                if (itFound->matchedLeftTemplateID == BreakpointEvidence::NOVEL_TEMPLATE)
                    newBp->leftTemplateID = newBp->rightTemplateID; // by twilight zone
                else
                    newBp->leftTemplateID = itFound->matchedLeftTemplateID; // by paired-end region

                // update positions
                newBp->leftPos.clear();
                newBp->leftPos.push_back(newPos);
                BreakpointCandidate::updateLeftMinMaxPos(newBp);

                // update depth
                if (this->optionManager->doReadDepthAnalysis())
                {
                    this->getReadDepth()->addUniformDepth(newBp->leftTemplateID, \
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
                    if (itFound->matchedLeftTemplateID == BreakpointEvidence::NOVEL_TEMPLATE)
                        newBp->leftTemplateID = newBp->rightTemplateID; // by twilight zone
                    else
                        newBp->leftTemplateID = itFound->matchedLeftTemplateID; // by paired-end region

                    // update positions
                    newBp->leftPos.clear();
                    newBp->leftPos.push_back(newPos);
                    BreakpointCandidate::updateLeftMinMaxPos(newBp);

                    // update depth
                    if (this->optionManager->doReadDepthAnalysis())
                    {
                        this->getReadDepth()->addUniformDepth(newBp->leftTemplateID, \
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
                    if (itFound->matchedRightTemplateID == BreakpointEvidence::NOVEL_TEMPLATE)
                        newBp->rightTemplateID = newBp->leftTemplateID; // by twilight zone
                    else
                        newBp->rightTemplateID = itFound->matchedRightTemplateID; // by paired-end region

                    // update positions
                    newBp->rightPos.clear();
                    newBp->rightPos.push_back(newPos);
                    BreakpointCandidate::updateRightMinMaxPos(newBp);

                    // update depth
                    if (this->optionManager->doReadDepthAnalysis())
                    {
                        this->getReadDepth()->addUniformDepth(newBp->rightTemplateID, \
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

            // update supporting read information
            ReadSupportInfo* newReadInfo = this->getMergedBreakpoint()->initReadSupport(newBp);
            newReadInfo->clippedReadSupport += clippedBp->suppReads;

            // add to the merged set
            bool isNewBp = false;
            newBp->needLeftIndexUpdate = true;
            newBp->needRightIndexUpdate = true;
            newBp = this->getMergedBreakpoint()->updateBreakpoint(newBp, isNewBp);

            // new bp : based on read-pairs or twilight zone
            if (isNewBp == true)
            {
                // based on read-pairs
                // - twilight zone: one-side has to be the novel_template
                if (itFound->matchedLeftTemplateID != BreakpointEvidence::NOVEL_TEMPLATE && \
                    itFound->matchedRightTemplateID != BreakpointEvidence::NOVEL_TEMPLATE)
                {
                    this->getPairedEndRead()->setBreakpointUsed(matchedBp, true);
                    newBp->suppReads += matchedBp->suppReads;
                    newReadInfo->pairedEndSupport += matchedBp->suppReads;
                }
            }
        }
   }
}

bool BreakpointManager::calculateReadDepth()
{
    TBreakpointSet* candidateSet = this->getMergedBreakpoint()->getCandidateSet();
    std::map<TTemplateID, unsigned> svCntByTemplate;
    for(auto it = candidateSet->begin(); it != candidateSet->end(); ++it)
    {
        Breakpoint* bp = *it;
        ReadSupportInfo* info = this->getMergedBreakpoint()->getReadSupport(bp);
        FinalBreakpointInfo* finalBreakpoint = &this->finalBreakpoints[bp];

        uint32_t windowSize = this->getOptionManager()->getReadDepthWindowSize();
        uint32_t breakpointSize = 1;

        // get breakpoint size if two positions are on a same template
        if (finalBreakpoint->leftTemplateID == finalBreakpoint->rightTemplateID)
        {
          if (finalBreakpoint->rightPosition > finalBreakpoint->leftPosition)
            breakpointSize += (finalBreakpoint->rightPosition - finalBreakpoint->leftPosition);
          else
            breakpointSize += (finalBreakpoint->leftPosition - finalBreakpoint->rightPosition);
        }

        // restrict the search region's size
        if (windowSize > breakpointSize && breakpointSize > 0)
            windowSize = breakpointSize;

        this->getReadDepth()->getReadDepthDiffScore(info->leftReadDepth, \
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
    this->getReadDepth()->calculateReadDepthStat(svCntByTemplate, candidateSet->size());

    return true;
}

void BreakpointManager::getNTCount(seqan::CharString& sequence, unsigned& a, unsigned& t, unsigned& g, unsigned& c)
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
    unsigned windowSize = this->getOptionManager()->getReadDepthWindowSize();
    windowSize = 1000;
    TBreakpointSet* candidateSet = this->getMergedBreakpoint()->getCandidateSet();
    for(auto it = candidateSet->begin(); it != candidateSet->end(); ++it)
    {
        FinalBreakpointInfo& finalBreakpoint = this->finalBreakpoints[*it];
        seqan::CharString ref;
        unsigned gc = 0, a = 0, t = 0, g = 0, c = 0;

        TPosition start = 0;
        TPosition end = BreakpointEvidence::INVALID_POS;
        TPosition refSize = 0;

        // left
        end = finalBreakpoint.leftPosition;
        start = std::max(end-windowSize, (unsigned) 0);
        this->getClippedRead()->getReferenceSequence(ref, finalBreakpoint.leftTemplateID, start, end);
        refSize += length(ref);
        this->getNTCount(ref, a, t, g, c);

        // right
        start = finalBreakpoint.rightPosition;
        end = std::min(start + windowSize, this->getReadDepth()->getRefSize(finalBreakpoint.rightTemplateID));
        this->getClippedRead()->getReferenceSequence(ref, finalBreakpoint.rightTemplateID, start, end);
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
    RUN(result,"Split-read analysis", this->getSplitRead()->analyze());
    if (result == false)
        return false;

    // paired-end analysis
    if (this->getOptionManager()->doPairedEndAnalysis())
    {
        RUN(result,"Paired-end analysis", this->getPairedEndRead()->analyze());
        if (result == false)
            return false;
    }

    // clipped-read analysis
    if (this->getOptionManager()->doClippedReadAnalysis())
    {
        RUN(result,"Clipped-read analysis", this->getMergedBreakpoint()->analyze());
        if (result == false)
            return false;
    }

    return true;
}

bool BreakpointManager::findFinalBreakpoints(void)
{
    TBreakpointSet* candidateSet = this->getMergedBreakpoint()->getCandidateSet();
    for(auto it = candidateSet->begin(); it != candidateSet->end(); ++it)
    {
        Breakpoint* bp = *it;
        ReadSupportInfo* info = this->getMergedBreakpoint()->getReadSupport(bp);

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
    TBreakpointSet* candidateSet = this->getPairedEndRead()->getCandidateSet();
    TPosition minSVSize = this->getOptionManager()->getMinSVSize();

    for(auto itBp = candidateSet->begin(); itBp != candidateSet->end(); ++itBp)
    {
        if( this->getPairedEndRead()->isBreakpointUsed(*itBp) == false)
        {
            Breakpoint* newBp = new Breakpoint;
            TPosition leftPos, rightPos, tempPos;
            this->getPairedEndRead()->copyBreakpoint(*newBp, *(*itBp));

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
            newBp = this->getMergedBreakpoint()->updateBreakpoint(newBp, isNew);

            // update read-info
            ReadSupportInfo* newReadInfo;
            newReadInfo = this->getMergedBreakpoint()->initReadSupport(newBp);
            newReadInfo->pairedEndSupport += newBp->suppReads;
        }
    }

    return true;
}

bool BreakpointManager::applyFilter(void)
{
  //if (this->optionManager->getUseRankAggregation())
  //  this->applyNormalization();
  this->filterByEvidenceSumAndVote();
  return true;
}

bool BreakpointManager::priByEvidenceSum(void)
{
    TBreakpointSet* candidateSet = this->getMergedBreakpoint()->getCandidateSet();
    auto itBreakpoint = candidateSet->begin();
    int cutoff = this->optionManager->getCutoff();

    while (itBreakpoint != candidateSet->end())
    {
        if (this->finalBreakpoints[*itBreakpoint].filtered == false)
        {
            FinalBreakpointInfo* finalBreakpointInfo = &this->finalBreakpoints[*itBreakpoint];
            ReadSupportInfo* info = this->getMergedBreakpoint()->getReadSupport(*itBreakpoint);
            finalBreakpointInfo->score =  (info->splitReadSupport + info->pairedEndSupport + info->clippedReadSupport);
        }
        ++itBreakpoint;
    }

    return true;
}


bool BreakpointManager::filterByEvidenceSumAndVote(void)
{
    double seTH = this->getOptionManager()->getMinSplitReadSupport();
    double peTH = this->getOptionManager()->getMinPairSupport();
    double reTH = this->getReadDepth()->getReTH();
    int32_t cutoff = this->getOptionManager()->getCutoff();
    int32_t minVote = this->getOptionManager()->getMinVote();
    if (minVote < 0)
    {
        minVote = 1;
        if (this->getOptionManager()->doPairedEndAnalysis())
            ++minVote;
        if (this->getOptionManager()->doReadDepthAnalysis())
            ++minVote;
        this->getOptionManager()->setMinVote(minVote);
    }

    if (minVote == 1)
    {
        printTimeMessage("Voting based resecuing is disabled automatically.");
        minVote = seqan::MaxValue<int32_t>::VALUE;
    }

    TBreakpointSet* candidateSet = this->getMergedBreakpoint()->getCandidateSet();
    auto itBreakpoint = candidateSet->begin();
    unsigned filteredBreakpointCount = 0;
    while (itBreakpoint != candidateSet->end())
    {
        FinalBreakpointInfo* finalBreakpointInfo = &this->finalBreakpoints[*itBreakpoint];
        ReadSupportInfo* info = this->getMergedBreakpoint()->getReadSupport(*itBreakpoint);

        // Evidence sum
        double score = (info->splitReadSupport + info->pairedEndSupport + info->clippedReadSupport);
        finalBreakpointInfo->score = score;
        finalBreakpointInfo->filtered = true;

        // Voting
        int32_t vote = 0;
        if ((info->splitReadSupport + info->clippedReadSupport) >= seTH)
            vote += 1;
        if (this->getOptionManager()->doPairedEndAnalysis() && info->pairedEndSupport >= peTH)
            vote += 1;
        if (this->getOptionManager()->doReadDepthAnalysis())
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

        // by evidence sum
        if (score >= (double)cutoff)
          finalBreakpointInfo->filtered = false;

        // by voting
        if ( (*itBreakpoint)->orientation != BreakpointEvidence::INVERTED )
        {
          if (vote >= minVote)
            finalBreakpointInfo->filtered = false;
        }
        else if( (*itBreakpoint)->orientation == BreakpointEvidence::INVERTED)
        {
          if (this->getOptionManager()->getUseREforBalancedSV() && vote >= minVote)
            finalBreakpointInfo->filtered = false;
        }

        if (finalBreakpointInfo->filtered == true)
          ++filteredBreakpointCount;

        ++itBreakpoint;
    }
    this->getMergedBreakpoint()->setFilteredBreakpointCount(filteredBreakpointCount);

    return true;
}

bool BreakpointManager::rescueByCombinedEvidence(void)
{
    double seTH = this->getOptionManager()->getMinSplitReadSupport();
    double peTH = this->getOptionManager()->getMinPairSupport();
    double reTH = this->getReadDepth()->getReTH();

    // Min. vote for rescue
    int minVote = this->getOptionManager()->getMinVote();
    if (minVote < 0)
    {
        minVote = 1;
        if (this->getOptionManager()->doPairedEndAnalysis())
            ++minVote;
        if (this->getOptionManager()->doReadDepthAnalysis())
            ++minVote;
        this->getOptionManager()->setMinVote(minVote);
    }

    if (minVote == 1)
    {
        printTimeMessage("Voting based resecuing is disabled automatically.");
        return true;
    }

    // for all breakpoints
    unsigned filteredBreakpointCount = this->getMergedBreakpoint()->getFilteredBreakpointCount();
    TBreakpointSet* candidateSet = this->getMergedBreakpoint()->getCandidateSet();
    auto itBreakpoint = candidateSet->begin();
    while (itBreakpoint != candidateSet->end())
    {
        FinalBreakpointInfo* finalBreakpointInfo = &this->finalBreakpoints[*itBreakpoint];
        ReadSupportInfo* info = this->getMergedBreakpoint()->getReadSupport(*itBreakpoint);

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
    this->mergedBreakpoints->setFilteredBreakpointCount(filteredBreakpointCount);

    return true;
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

    this->splitReadBreakpoints = new SplitRead(this->optionManager);
    this->pairedEndBreakpoints = new PairedEndRead(this->optionManager);
    this->clippedBreakpoints = new ClippedRead(this->optionManager);
    this->readDepthBreakpoints = new ReadDepth(this->optionManager);
    this->mergedBreakpoints = new MergedCandidate(this->optionManager);

    // register objects to alignment manager
    this->alignmentManager->setBreakpointCandidate(splitReadBreakpoints, \
                                                   pairedEndBreakpoints, \
                                                   clippedBreakpoints, \
                                                   readDepthBreakpoints);
}

BreakpointManager::~BreakpointManager()
{
    delete this->splitReadBreakpoints;
    delete this->pairedEndBreakpoints;
    delete this->clippedBreakpoints;
    delete this->readDepthBreakpoints;
    delete this->mergedBreakpoints;
}
