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
#include <seqan/align.h>
#include <seqan/index.h>
#include <seqan/find.h>
#include <seqan/store.h>
#include <seqan/consensus.h>
#include <seqan/align_split.h>
#include "misc.hpp"
#include "clippedread.hpp"

void ClippedRead::prepAfterHeaderParsing(seqan::BamHeader& header, seqan::BamFileIn& fileIn, bool isLongRead)
{
    for (unsigned i=0; i < length(header); ++i)
    {
        if (header[i].type == seqan::BamHeaderRecordType::BAM_HEADER_REFERENCE)
        {
            seqan::CharString templateName, templateLengthStr;
            seqan::getTagValue(templateName, "SN", header[i]);
            seqan::getTagValue(templateLengthStr, "LN", header[i]);
            seqan::toCString(templateLengthStr);

            TTemplateID rID = BreakpointEvidence::NOVEL_TEMPLATE;
            TPosition templateLength = BreakpointEvidence::INVALID_POS;
            getIdByName(rID, contigNamesCache(seqan::context(fileIn)), templateName);
            templateLength = std::stol( CharStringToStdString(templateLengthStr) );

            this->templateSize[rID] = templateLength;
        }
    }

    if ( this->getOptionManager()->doClippedReadAnalysis() == true )
    {
        // load fasta file
        seqan::CharString pathToFaFile = this->getOptionManager()->getReferenceGenome();
        seqan::CharString pathToFaiFile = pathToFaFile;
        seqan::append(pathToFaiFile, ".fai");
        if (!open(this->faiIndex, seqan::toCString(pathToFaFile)))
        {
            printTimeMessage("ClippedRead : Can't load '" + CharStringToStdString(pathToFaiFile) + "', try to build it.");
            if (!build(this->faiIndex, seqan::toCString(pathToFaFile)))
            {
                printMessage("[ERROR] ClippedRead : Can't build FASTA index with '" + CharStringToStdString(pathToFaFile) + \
                             "'. Please check 'referenceGenome' option with '--help'.");
                exit(1);
            }
            else
            {
                if (!save(this->faiIndex, seqan::toCString(pathToFaiFile)))
                {
                    printMessage("ClippedRead : Can't save FASTA index to '" + CharStringToStdString(pathToFaiFile) + "'.");
                    exit(1);
                }

                printTimeMessage("ClippedRead : FASTA index was built.");
            }
        }
    }
}

void ClippedRead::setSearchRegionByOrientation(const BreakpointEvidence::ORIENTATION orientation, const BreakpointEvidence::SIDE side, Breakpoint &bp, TTemplateID &rID, TPosition &begin, TPosition &end)
{
    if (orientation == BreakpointEvidence::ORIENTATION::PROPERLY_ORIENTED_LARGE)
    {
        if (side == BreakpointEvidence::SIDE::RIGHT)
        {
            begin = bp.minRightPos;
            end = bp.maxRightPos;
            rID = bp.rightTemplateID;
        }
        else
        {
            begin = bp.minLeftPos;
            end = bp.maxLeftPos;
            rID = bp.leftTemplateID;
        }
    }
    else if (orientation == BreakpointEvidence::ORIENTATION::INVERTED)
    {
        if (side == BreakpointEvidence::SIDE::RIGHT)
        {
            begin = bp.minRightPos;
            end = bp.maxRightPos;
            rID = bp.rightTemplateID;
        }
        else
        {
            begin = bp.minLeftPos;
            end = bp.maxLeftPos;
            rID = bp.leftTemplateID;
        }
    }
    else if (orientation == BreakpointEvidence::ORIENTATION::SWAPPED)
    {
        if (side == BreakpointEvidence::SIDE::RIGHT)
        {
            begin = bp.minRightPos;
            end = bp.maxRightPos;
            rID = bp.rightTemplateID;
        }
        else
        {
            begin = bp.minLeftPos;
            end = bp.maxLeftPos;
            rID = bp.leftTemplateID;
        }
    }

    // extend search region
    int32_t PESearchSize = this->getOptionManager()->getPairedEndSearchSize();
    int32_t searchSize = (end - begin + 1);
    if ( searchSize < PESearchSize)
        BreakpointCandidate::setPositionWithAdj( begin, end, (PESearchSize - searchSize) / 2.0);
}

bool ClippedRead::searchPairRegion (TFoundPosition& foundPositions, \
                                    Breakpoint* bp, \
                                    int32_t &bestScore, \
                                    seqan::CharString& query, \
                                    BreakpointEvidence::SIDE searchSide, \
                                    bool useLocalAlignment, \
                                    bool useReverseComplement, \
                                    BreakpointEvidence::ORIENTATION orientation)
{
    // get search region
    TTemplateID refID;
    TPosition refBeginPosInGenome, refEndPosInGenome;
    setSearchRegionByOrientation(orientation, searchSide, *bp, refID, refBeginPosInGenome, refEndPosInGenome);
    if(refBeginPosInGenome > refEndPosInGenome)
        return false;

    // get reference sequence
    seqan::CharString ref;
    this->getReferenceSequence(ref, refID, refBeginPosInGenome, refEndPosInGenome);
    if (ref == "")
        return false;

    // search
    TFoundPosition newFoundPositions;
    int32_t newBestScore = bestScore;

    bool foundBestMatch = false;
    if (useLocalAlignment == true)
        foundBestMatch = this->alignByLocal(newFoundPositions, ref, query, newBestScore);
    else
        foundBestMatch = this->alignByMyersBitVector(newFoundPositions, ref, query, newBestScore);

    if (foundBestMatch == true)
    {
        // found best match
        if (newBestScore >= bestScore)
        {
            // new best score
            if (newBestScore > bestScore)
            {
                bestScore = newBestScore;
                foundPositions.clear();
            }

            // add to the list
            for (auto itNewFound = newFoundPositions.begin(); itNewFound != newFoundPositions.end(); ++itNewFound)
            {
                itNewFound->sequenceSegment.templateID = refID;
                itNewFound->sequenceSegment.beginPos += refBeginPosInGenome;
                itNewFound->sequenceSegment.endPos += refBeginPosInGenome;
                itNewFound->isReverseComplemented = useReverseComplement;
                itNewFound->orientation = orientation;
                itNewFound->matchedBp = bp;
                itNewFound->matchedRightTemplateID = bp->rightTemplateID;
                itNewFound->matchedLeftTemplateID = bp->leftTemplateID;

                // check duplicates
                bool duplicated = false;
                for (auto itFound = foundPositions.begin(); itFound != foundPositions.end(); ++itFound)
                {
                    if ( (itNewFound->sequenceSegment.beginPos == itFound->sequenceSegment.beginPos) && \
                         (itNewFound->sequenceSegment.endPos == itFound->sequenceSegment.endPos))
                    {
                        duplicated = true;
                        break;
                    }
                }
                if (duplicated == false)
                    foundPositions.push_back(*itNewFound);
            }
        }
        return true;
    }
    return false;
}

bool ClippedRead::searchTwilightZone (TFoundPosition& foundPositions, \
                                      Breakpoint* bp, \
                                      int32_t& bestScore, \
                                      seqan::CharString& query, \
                                      BreakpointEvidence::SIDE searchSide, \
                                      bool useReverseComplement, \
                                      BreakpointEvidence::ORIENTATION orientation)
{
    // get search region
    int32_t minSVSize = this->getOptionManager()->getMinSVSize();
    int32_t searchSize = this->getMaxAbInsSize();
    TTemplateID refID;
    TPosition refBeginPosInGenome, refEndPosInGenome;
    if (searchSide == BreakpointEvidence::SIDE::RIGHT)
    {
        // TODO : duplication? (search for both-side)
        refBeginPosInGenome = bp->maxLeftPos + minSVSize;
        refEndPosInGenome = refBeginPosInGenome + searchSize;
        refID = bp->leftTemplateID;
    }
    else
    {
        refEndPosInGenome = bp->minRightPos - minSVSize;
        refBeginPosInGenome = refEndPosInGenome - searchSize;
        refID = bp->rightTemplateID;
    }
    if(refBeginPosInGenome > refEndPosInGenome)
        return false;

    // get reference sequence
    seqan::CharString ref;
    this->getReferenceSequence(ref, refID, refBeginPosInGenome, refEndPosInGenome);
    if (ref == "")
        return false;

    // search
    TFoundPosition newFoundPositions;
    int32_t newBestScore = bestScore;
    if (this->alignByMyersBitVector(newFoundPositions, ref, query, newBestScore) == true)
    {
        // found best match
        if (newBestScore >= bestScore)
        {
            // new best score
            if (newBestScore > bestScore)
            {
                bestScore = newBestScore;
                foundPositions.clear();
            }

            // add to the list
            for (auto itNewFound = newFoundPositions.begin(); itNewFound != newFoundPositions.end(); ++itNewFound)
            {
                itNewFound->sequenceSegment.templateID = refID;
                itNewFound->sequenceSegment.beginPos += refBeginPosInGenome;
                itNewFound->sequenceSegment.endPos += refBeginPosInGenome;
                itNewFound->isReverseComplemented = useReverseComplement;
                itNewFound->orientation = orientation;
                itNewFound->matchedBp = bp;
                itNewFound->matchedRightTemplateID = bp->rightTemplateID;
                itNewFound->matchedLeftTemplateID = bp->leftTemplateID;

                // check duplicates
                bool duplicated = false;
                for (auto itFound = foundPositions.begin(); itFound != foundPositions.end(); ++itFound)
                {
                    if ( (itNewFound->sequenceSegment.beginPos == itFound->sequenceSegment.beginPos) && \
                         (itNewFound->sequenceSegment.endPos == itFound->sequenceSegment.endPos))
                    {
                        duplicated = true;
                        break;
                    }
                }
                if (duplicated == false)
                    foundPositions.push_back(*itNewFound);
            }
        }
        return true;
    }
    return false;
}

bool ClippedRead::onlineSearchBySegment(TFoundPosition& foundPositions, seqan::CharString &ref, seqan::CharString &query, int32_t &bestScore)
{
    // search
    int32_t currentScore = 0;

    int32_t k = 20;
    int32_t editDistanceInKmer = this->getOptionManager()->getClippedSeqErrorRate() * k;
    int32_t editDistanceInTotal = this->getOptionManager()->getClippedSeqErrorRate() * seqan::length(query);
    std::map<seqan::CharString, std::vector<TPosition> > kmers;
    for(unsigned i=0; i < seqan::length(query) - k + 1; ++i)
        kmers[seqan::infix(query, i, i+k)].push_back(i);

    // search minimal match with Myer's bitvector
    seqan::Finder<seqan::CharString> finder(ref);
    seqan::Pattern<seqan::CharString, seqan::Myers<> > patternQuery(query);
    for (auto itKmer = kmers.begin(); itKmer != kmers.end(); ++itKmer)
    {
        seqan::CharString kmer = itKmer->first;
        seqan::Pattern<seqan::CharString, seqan::Myers<> > patternKmer(kmer);

        goBegin(finder);
        seqan::clear(finder);
        while (find(finder, patternKmer, -editDistanceInKmer))
        {
            int32_t kmerScore = seqan::getScore(patternKmer);
            int32_t editDistanceRemains = editDistanceInTotal + kmerScore;
            while (seqan::findBegin(finder, patternKmer, kmerScore))
            {
                // extend search region
                int32_t subjectBeginPos = seqan::beginPosition(finder);
                for (auto itPos = itKmer->second.begin(); itPos != itKmer->second.end(); ++itPos)
                {
                    int32_t queryBeginPos = *itPos;
                    int32_t queryEndPos = queryBeginPos + k;
                    int32_t subjectEndPos = subjectBeginPos + k;

                    // search region
                    int32_t subjectSearchBeginPos = subjectBeginPos - (queryBeginPos + editDistanceRemains);
                    int32_t subjectSearchEndPos = subjectEndPos + ((seqan::length(query) - queryEndPos) + editDistanceRemains);
                    if (subjectSearchBeginPos < 0)
                        subjectSearchBeginPos = 0;
                    if (subjectSearchEndPos > length(ref))
                        subjectSearchEndPos = length(ref);

                    // semi-global alignment (allow begin/end gaps in query)
                    seqan::Align<seqan::CharString> align;
                    resize(seqan::rows(align), 2);
                    assignSource(row(align, 0), seqan::infix(ref, subjectSearchBeginPos, subjectSearchEndPos));
                    assignSource(row(align, 1), query);
                    currentScore = globalAlignment(align, seqan::Score<int, seqan::Simple>(0, -1, -1), seqan::AlignConfig<false, true, true, false>());

                    // got new best alignment
                    if (currentScore >= bestScore)
                    {
                        if (currentScore > bestScore)
                        {
                            foundPositions.clear();
                            bestScore = currentScore;
                        }

                        ClippedSequenceSegment s;
                        s.sequenceSegment.beginPos = clippedBeginPosition(row(align, 0)) + subjectSearchBeginPos;
                        s.sequenceSegment.endPos = clippedEndPosition(row(align, 0)) - 1 + subjectSearchBeginPos;

                        bool duplicated = false;
                        for (auto itFound = foundPositions.begin(); itFound != foundPositions.end(); ++itFound)
                        {
                            if ( (itFound->sequenceSegment.beginPos == s.sequenceSegment.beginPos) && \
                                 (itFound->sequenceSegment.endPos == s.sequenceSegment.endPos))
                            {
                                duplicated = true;
                                break;
                            }
                        }

                        if (duplicated == false)
                           foundPositions.push_back(s);
                    }
                }

            }
        }
    }

    return false;
}


bool ClippedRead::alignByLocal(TFoundPosition& foundPositions, seqan::CharString& ref, seqan::CharString& query, int32_t& bestScore)
{
    seqan::Align<seqan::CharString> align;
    resize(seqan::rows(align), 2);
    assignSource(row(align, 0), ref);
    assignSource(row(align, 1), query);

    int32_t currentScore = localAlignment(align, seqan::Score<int>(1, -1, -1, -2));
    //int32_t currentScore = globalAlignment(align, seqan::Score<int, seqan::Simple>(3, -3, -2, -2), seqan::AlignConfig<true, true, true, true>());
    if (currentScore >= bestScore)
    {
        if (currentScore > bestScore)
        {
            foundPositions.clear();
            bestScore = currentScore;
        }

        ClippedSequenceSegment s;
        s.sequenceSegment.beginPos = clippedBeginPosition(row(align, 0));
        s.sequenceSegment.endPos = clippedEndPosition(row(align, 0));
        s.querySegment.beginPos = seqan::beginPosition(row(align, 1));
        s.querySegment.endPos = seqan::endPosition(row(align, 1));

        /*
        std::cerr << align << "\n";
        std::cerr << s.sequenceSegment.beginPos << "\n";
        std::cerr << s.sequenceSegment.endPos << "\n";
        std::cerr << s.querySegment.beginPos << "\n";
        std::cerr << s.querySegment.endPos << "\n";
        std::cerr << seqan::beginPosition(row(align, 1)) << "\n";
        std::cerr << seqan::endPosition(row(align, 1)) << "\n";
        std::cerr << query << "\n";
        */

        foundPositions.push_back(s);
    }

    return (foundPositions.size() > 0);
}


bool ClippedRead::alignByMyersBitVector(TFoundPosition& foundPositions, seqan::CharString& ref, seqan::CharString& query, int32_t& bestScore)
{
    int32_t currentScore = 0;
    int32_t currentAlignSize = 0;
    int32_t bestAlignSize = 0;

    seqan::Finder<seqan::CharString> finder(ref);
    seqan::Pattern<seqan::CharString, seqan::Myers<> > pattern(query);
    while (find(finder, pattern, bestScore))
    {
        while (seqan::findBegin(finder, pattern, seqan::getScore(pattern)))
        {
            currentScore = seqan::getScore(pattern);
            currentAlignSize = seqan::endPosition(finder) - seqan::beginPosition(finder);

            if (currentScore >= bestScore && currentAlignSize >= bestAlignSize)
            {
                if (currentScore > bestScore || currentAlignSize > bestAlignSize)
                {
                    foundPositions.clear();
                    bestScore = currentScore;
                    bestAlignSize = currentAlignSize;
                }

                ClippedSequenceSegment s;
                s.sequenceSegment.beginPos = seqan::beginPosition(finder);
                s.sequenceSegment.endPos = seqan::endPosition(finder);
                foundPositions.push_back(s);
            }
        }
    }

    return (foundPositions.size() > 0);
}

void ClippedRead::getReferenceSequence(seqan::CharString& seq, TTemplateID templateID, TPosition start, TPosition end)
{
    seqan::clear(seq);
    if (templateID == BreakpointEvidence::NOVEL_TEMPLATE || end >= templateSize[templateID])
        return;

    try {
        seqan::readRegion(seq, this->faiIndex, templateID, start, end);
    }
    catch (...)
    {
        ;
    }
}

void ClippedRead::getReferenceSequence(seqan::CharString& seq, seqan::CharString chr, TPosition start, TPosition end)
{
    unsigned templateID;
    if (!getIdByName(templateID, this->faiIndex, chr))
    {
        seqan::clear(seq);
        return;
    }
    else
        getReferenceSequence(seq, templateID, start, end);
}

void ClippedRead::parseReadRecord(TReadName &readName, seqan::BamAlignmentRecord &record)
{
    if (record.mapQ < this->getOptionManager()->getMinMapQual())
        return;

    // parse CIGAR
    AlignmentInfo alnInfo;
    TReadID readID = this->getNextReadID();
    parseCIGAR(alnInfo, readID, record, false, true, this->getOptionManager()->getMinSVSize(), this->getOptionManager()->getMinClippedSeqSize());

    // for all clipped segments
    bool isNew;
    auto itClip = alnInfo.clippedList.begin();
    while (itClip != alnInfo.clippedList.end())
    {
        BreakpointEvidence& be = *(itClip++);
        be.orientation = BreakpointEvidence::ORIENTATION::CLIPPED;

        Breakpoint* bp = this->updateBreakpoint(be, isNew);
        if (isNew == true)
        {
            bp->bFoundExactPosition = false;
            bp->clippedConsensusSequenceSize = 0;
        }
    }
}

void ClippedRead::getConsensusSequence(seqan::CharString& seq, Breakpoint* bp)
{
    seqan::clear(seq);
    if (bp->clippedSequences.size() == 1)
    {
        seq = (bp->clippedSequences.begin())->second;
    }
    else if (bp->clippedSequences.size() == 2)
    {
        if (length((bp->clippedSequences.begin())->second) > \
            length((bp->clippedSequences.begin()+1)->second))
            seq = (bp->clippedSequences.begin())->second;
        else
            seq = (bp->clippedSequences.begin()+1)->second;
    }
    else
    {
        unsigned maxSeqIdx=0;
        unsigned maxSeqSize=0;

        // do assembly
        if (this->getOptionManager()->isUsingAssembly())
        {
            seqan::FragmentStore<> store;
            for (auto it = bp->clippedSequences.begin(); it != bp->clippedSequences.end(); ++it)
                appendRead(store, (*it).second);

            seqan::ConsensusAlignmentOptions options;
            options.useContigID = false;
            options.useGlobalAlignment = true;
            seqan::consensusAlignment(store, options);

            for( unsigned i=0; i < length(store.contigStore); ++i)
            {
                if (length(store.contigStore[i].seq) > maxSeqSize)
                {
                    maxSeqSize = length(store.contigStore[i].seq);
                    maxSeqIdx = i;
                }
            }
            seq = store.contigStore[maxSeqIdx].seq;
        }
        else // select the longest sequence
        {
            for(int i=0; i < bp->clippedSequences.size(); ++i)
            {
                if (length(bp->clippedSequences[i].second) > maxSeqSize)
                {
                    maxSeqSize = length(bp->clippedSequences[i].second);
                    maxSeqIdx = i;
                }
            }
            seq = bp->clippedSequences[maxSeqIdx].second;
        }
    }
    bp->clippedConsensusSequenceSize = length(seq);
}
