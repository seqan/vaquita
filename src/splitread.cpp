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
#include "misc.hpp"
#include "splitread.hpp"
#include <seqan/bam_io.h>

bool SplitRead::updateBreakpointByIndels(std::vector<BreakpointEvidence>& beList)
{
    bool isNew;
    for (unsigned i=0; i < beList.size(); ++i) // for each 
    {
        BreakpointEvidence& be = beList[i];
        Breakpoint* bp = this->updateBreakpoint(be, isNew);

        // every breakpoints based on split-read have exact positions.
        bp->bFoundExactPosition = true; 
    }
}

bool SplitRead::analyzeRead(TReadName& readName)
{
    auto itRead = this->recordByRead.find(readName);
    if (itRead == this->recordByRead.end())
        return false;

    if (itRead->second.size() > 2) // not confident (too-many segments)
        return false;
    
    // get readID
    TReadID readID = this->getNextReadID(); // this must be unique
    
    // extract information form local alignments
    uint32_t minQuality = 999999;
    std::vector<AlignmentInfo> localAligns;
    TPosition querySize = 0;
    for(auto it = itRead->second.begin(); it != itRead->second.end(); ++it) // for each segment
    {
        AlignmentInfo alnInfo;
        parseCIGAR(alnInfo, readID, *it, true, false, this->getOptionManager()->getMinSVSize(), 0);
        localAligns.push_back(alnInfo);

        // min. quality segments
        if (it->mapQ < minQuality)
            minQuality = it->mapQ;

        // size of the read
        if (alnInfo.querySize > querySize)
            querySize = alnInfo.querySize;
    }

    // quailiy check
    if (minQuality < this->getOptionManager()->getMinMapQual())
        return false;

    ///////////////////////////////////////////////////////////////////////////////////
    // TODO: These addtional quailiy controls must be shown to users (provide options)
    unsigned maxOverlapSize = 20;
    ///////////////////////////////////////////////////////////////////////////////////

    // sort according to the start position in the query
    sort(localAligns.begin(), localAligns.end(), BreakpointCandidate::compareByQueryPos);
  
    // check indels in the the first segment
    this->updateBreakpointByIndels( FIRST_ELEMENT(localAligns).indelList );

    // for each alignment
    for(unsigned int i=1; i < localAligns.size(); ++i)
    {
        // check indels
        this->updateBreakpointByIndels(localAligns[i].indelList);

        // check adjacency and overlaps of two consecutive segments
        SequenceSegment& rSeg1 = localAligns[i-1].refSegment;
        SequenceSegment& qSeg1 = localAligns[i-1].querySegment;
        SequenceSegment& rSeg2 = localAligns[i].refSegment;
        SequenceSegment& qSeg2 = localAligns[i].querySegment;

        bool skipThis = true; // this segment is not useuful
        unsigned overlapSize = maxOverlapSize + 1; // will be filtered by default. (must be rescued)

        // same strand
        if (qSeg1.isReverse == qSeg2.isReverse) 
        {
            if (this->isAdjacent(qSeg1, qSeg2)) // adjacent
                skipThis = false; 
            else if (this->isOverlap(qSeg1, qSeg2)) // overlap
                overlapSize = abs(qSeg2.beginPos-qSeg1.endPos);
        }
        else // different strand
        {
            // change coordinates according to the one strand
            SequenceSegment _qSeg2;
            _qSeg2.beginPos = querySize - qSeg2.endPos;
            _qSeg2.endPos = querySize - qSeg2.beginPos;

            if (this->isAdjacent(qSeg1, _qSeg2)) // adjacent
                skipThis = false; 
            else if (this->isOverlap(qSeg1, _qSeg2)) // overalp
            {
                if (qSeg1.beginPos <= _qSeg2.beginPos)
                    overlapSize = abs(_qSeg2.beginPos - qSeg1.endPos);
                else
                    overlapSize = abs(qSeg1.beginPos - _qSeg2.endPos);
            }
        }

        // rescued by overlap
        if ( (skipThis == true) && overlapSize <= maxOverlapSize )
            skipThis = false;

        // find breakpoint
        if (skipThis == false)
        {   
            SequenceSegment* leftRefSeg; 
            SequenceSegment* rightRefSeg;

            BreakpointEvidence be;
            be.suppRead = 1;

            // decided by positions in reference
            if (rSeg1.templateID == rSeg2.templateID)
            {
                if (rSeg1.isReverse == rSeg2.isReverse) // same strand
                {
                    if (rSeg1.beginPos < rSeg2.beginPos)
                    {
                        leftRefSeg = &rSeg1;
                        rightRefSeg = &rSeg2;
                        be.orientation = BreakpointEvidence::ORIENTATION::PROPERLY_ORIENTED_LARGE; // - ->
                        be.leftSegment.beginPos = leftRefSeg->endPos;
                        be.rightSegment.beginPos = rightRefSeg->beginPos - 1;
                    }
                    else
                    {
                        leftRefSeg = &rSeg2;
                        rightRefSeg = &rSeg1;
                        be.orientation = BreakpointEvidence::ORIENTATION::SWAPPED; // -> -
                        be.leftSegment.beginPos = leftRefSeg->beginPos;   
                        be.rightSegment.beginPos = rightRefSeg->endPos - 1;
                    }
                }
                else // different strand
                {
                    be.orientation = BreakpointEvidence::ORIENTATION::INVERTED;

                    bool orderKept, firstSegIsReverse;
                    if (rSeg1.beginPos <= rSeg2.beginPos)
                    {
                        leftRefSeg = &rSeg1;
                        rightRefSeg = &rSeg2;
                        orderKept = true;
                    }
                    else
                    {
                        leftRefSeg = &rSeg2;
                        rightRefSeg = &rSeg1;
                        orderKept = false;
                    }
                    firstSegIsReverse = leftRefSeg->isReverse;

                    if (orderKept)
                    {
                        // <-| -|
                        if (firstSegIsReverse)
                        {
                            be.leftSegment.beginPos = leftRefSeg->endPos;
                            be.rightSegment.beginPos = rightRefSeg->endPos - 1;
                        }
                        else // |-> |-
                        {
                            be.leftSegment.beginPos = leftRefSeg->beginPos;
                            be.rightSegment.beginPos = rightRefSeg->beginPos - 1;
                        }
                    }
                    else
                    {
                        // -| <-|
                        if (firstSegIsReverse)
                        {
                            be.leftSegment.beginPos = leftRefSeg->endPos;
                            be.rightSegment.beginPos = rightRefSeg->endPos - 1;
                        }
                        else // |- |->
                        {
                            be.leftSegment.beginPos = leftRefSeg->beginPos;
                            be.rightSegment.beginPos = rightRefSeg->beginPos - 1;
                        }   
                    }
                }

                be.leftSegment.templateID = leftRefSeg->templateID;
                be.leftSegment.isReverse = leftRefSeg->isReverse;
                be.leftSegment.endPos = be.leftSegment.beginPos + 1;
                be.rightSegment.templateID = rightRefSeg->templateID;
                be.rightSegment.isReverse = rightRefSeg->isReverse;
                be.rightSegment.endPos = be.rightSegment.beginPos + 1;
            }
            else // decided by positions in query
            {
                SequenceSegment* leftQuerySeg; 
                SequenceSegment* rightQuerySeg;

                leftQuerySeg = &qSeg1;
                leftRefSeg = &rSeg1;
                rightQuerySeg = &qSeg2;
                rightRefSeg = &rSeg2;

                be.orientation = BreakpointEvidence::ORIENTATION::NOT_DECIDED;

                be.leftSegment.templateID = leftRefSeg->templateID;
                be.leftSegment.isReverse = leftRefSeg->isReverse;
                be.leftSegment.beginPos = leftRefSeg->endPos;
                be.leftSegment.endPos = be.leftSegment.beginPos + 1;

                be.rightSegment.templateID = rightRefSeg->templateID;
                be.rightSegment.isReverse = rightRefSeg->isReverse;               
                be.rightSegment.beginPos = rightRefSeg->beginPos - 1;
                be.rightSegment.endPos = be.rightSegment.beginPos + 1;
            }

            bool isNew;
            Breakpoint *newBp = this->updateBreakpoint(be, isNew);
            newBp->bFoundExactPosition = true;                
        }
    }

    return true;
}

bool SplitRead::analyzeCurrentChromosome()
{
    auto itRead = this->recordByRead.begin();
    while (itRead != this->recordByRead.end()) // for each reads
    {
        // skip interchoromosomal reads
        if (this->multChrmRead[itRead->first] == true)
        {
            ++itRead;
            continue;
        }
        else
        {
            if (itRead->second.size() > 0)
                analyzeRead(const_cast<TReadName&>(itRead->first));

            // release memory
            this->multChrmRead.erase( this->multChrmRead.find(const_cast<TReadName&>(itRead->first)) );
            itRead = this->recordByRead.erase(itRead);
        }
    }

    return true;
}

bool SplitRead::analyze()
{
    // read all remaining records
    auto itRead = this->recordByRead.begin();
    while (itRead != this->recordByRead.end()) // for each reads
    {
        if (itRead->second.size() > 0)
            analyzeRead(const_cast<TReadName&>(itRead->first));
        itRead = this->recordByRead.erase(itRead);
    }
    this->multChrmRead.clear();

    return true;
}

void SplitRead::checkReadRecord(TReadName &readName, BamAlignmentRecord &record)
{
    TReadName nextPairName = readName;

    if ( hasFlagFirst(record) )
    {
        readName += "/1";
        nextPairName += "/2";
    }
    else // ( hasFlagLast(record) )
    {
        readName += "/2";        
        nextPairName += "/1";
    }    

    // add read
    if (this->recordByRead.find(readName) == this->recordByRead.end())
    {
        this->recordByRead[readName] = TBamRecords();
        this->multChrmRead[readName] = false;
    }
}

void SplitRead::parseReadRecord(TReadName &readName, BamAlignmentRecord &record)
{
    this->checkReadRecord(readName, record);

    // seq/qual/tag informations are removed after this.
    record.seq = "";
    record.qual = "";
    record.tags = "";
    record._buffer = "";
    resize(record.seq,0);
    resize(record.qual,0);
    resize(record.tags,0);
    resize(record._buffer,0);

    // first record
    if (lastChrmId == BreakpointEvidence::NOVEL_TEMPLATE)
        lastChrmId = record.rID;

    // inter-chromosomal read
    if ( record.rID != record.rNextId )
        this->multChrmRead[readName] = true;

    // chromosomal change (records are sorted by position)
    if (lastChrmId != record.rID)
    {
        bool result;
        RUN(result,"Analyze chromosome " + std::string(toCString(contigNames(context(*this->fileIn))[lastChrmId])), this->analyzeCurrentChromosome());
        lastChrmId = record.rID;
    }

    // add read
    this->recordByRead[readName].push_back(std::move(record));
}

void SplitRead::prepAfterHeaderParsing(BamHeader& header, BamFileIn& fileIn)
{
    this->fileIn = &fileIn;
}