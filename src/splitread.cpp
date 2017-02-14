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
#include "misc.hpp"
#include "splitread.hpp"
#include <seqan/bam_io.h>

bool SplitRead::analyzeRead(TReadName& readName)
{
    bool isNew;
    auto itRead = this->recordByRead.find(readName);
    if (itRead == this->recordByRead.end())
        return false;

    // not confident
    if (itRead->second.size() > 2)
        return false;

    // extract information
    TReadID readID = this->getAndUpdateCurrentReadID(); // this must be unique
    std::vector<AlignmentInfo> localAligns;
    TPosition querySize = 0;
    uint32_t minQuality = 999999;
    for(auto it = itRead->second.begin(); it != itRead->second.end(); ++it) // for each segment
    {
        AlignmentInfo alnInfo;
        parseCIGAR(alnInfo, readID, *it, true, false, this->getOptionManager()->getMinSVSize(), 0);
        localAligns.push_back(alnInfo);
        if (alnInfo.querySize > querySize)
            querySize = alnInfo.querySize;
        if (it->mapQ < minQuality)
            minQuality = it->mapQ;
    }

    if (minQuality < this->getOptionManager()->getMinMapQual())
        return false;

    // sort according to the start position in the query
    sort(localAligns.begin(), localAligns.end(), BreakpointCandidate::compareByQueryPos);

    // check indels for the first segment
    this->updateBreakpoint( FIRST_ELEMENT(localAligns).indelList, true );

    // for each alignment
    //unsigned maxOverlapSize = this->getOptionManager()->getAdjTol();
    unsigned maxOverlapSize = 20;
    //unsigned maxOverlapSize = 9999;
    for(unsigned int i=1; i < localAligns.size(); ++i)
    {
        // check indels
        this->updateBreakpoint(localAligns[i].indelList, true);

        SequenceSegment& rSeg1 = localAligns[i-1].refSegment;
        SequenceSegment& qSeg1 = localAligns[i-1].querySegment;
        SequenceSegment& rSeg2 = localAligns[i].refSegment;
        SequenceSegment& qSeg2 = localAligns[i].querySegment;

        bool skipThis = true;
        unsigned overlapSize = maxOverlapSize + 1;
        if (qSeg1.isReverse == qSeg2.isReverse) // same strand
        {
            if (this->isAdjacent(qSeg1, qSeg2))
                skipThis = false; // adjacent
            else if (this->isOverlap(qSeg1, qSeg2))
                overlapSize = abs(qSeg2.beginPos-qSeg1.endPos);
        }
        else
        {
            // change coordinates according to the one strand
            SequenceSegment _qSeg2;
            _qSeg2.beginPos = querySize - qSeg2.endPos;
            _qSeg2.endPos = querySize - qSeg2.beginPos;

            if (this->isAdjacent(qSeg1, _qSeg2))
                skipThis = false; // adjacent
            else if (this->isOverlap(qSeg1, _qSeg2))
            {
                if (qSeg1.beginPos <= _qSeg2.beginPos)
                    overlapSize = abs(_qSeg2.beginPos - qSeg1.endPos);
                else
                    overlapSize = abs(qSeg1.beginPos - _qSeg2.endPos);
            }
        }
        if ( (skipThis == true) && overlapSize <= maxOverlapSize )
            skipThis = false;

        if (skipThis == false)
        {   
            SequenceSegment* leftRefSeg; 
            SequenceSegment* rightRefSeg;

            BreakpointEvidence be;
            be.suppRead = 1;

            if (rSeg1.isReverse == rSeg2.isReverse) // same strand
            {
                if (rSeg1.beginPos < rSeg2.beginPos)
                {
                    leftRefSeg = &rSeg1;
                    rightRefSeg = &rSeg2;
                    be.orientation = BreakpointEvidence::ORIENTATION::PROPERLY_ORIENTED; // - ->
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
            else
            {
                be.orientation = BreakpointEvidence::ORIENTATION::INVERSED;

                if (rSeg1.beginPos < rSeg2.beginPos)
                {
                    leftRefSeg = &rSeg1;
                    rightRefSeg = &rSeg2;

                    if (leftRefSeg->isReverse == false) // - |    <-|
                    {
                        be.leftSegment.beginPos = leftRefSeg->endPos;
                        be.rightSegment.beginPos = rightRefSeg->endPos - 1;
                    }
                    else // |-    |->
                    {
                        be.leftSegment.beginPos = leftRefSeg->beginPos;
                        be.rightSegment.beginPos = rightRefSeg->beginPos - 1;
                    }
                }
                else
                {
                    leftRefSeg = &rSeg2;
                    rightRefSeg = &rSeg1;

                    if (leftRefSeg->isReverse == false) // |->    |-
                    {
                        be.leftSegment.beginPos = leftRefSeg->beginPos;
                        be.rightSegment.beginPos = rightRefSeg->beginPos - 1;
                    }
                    else // <-|   -|
                    {
                        be.leftSegment.beginPos = leftRefSeg->endPos;
                        be.rightSegment.beginPos = rightRefSeg->endPos - 1;
                    }
                }
            }

            be.leftSegment.templateID = leftRefSeg->templateID;
            be.leftSegment.isReverse = leftRefSeg->isReverse;
            be.leftSegment.endPos = be.leftSegment.beginPos + 1;
            be.rightSegment.templateID = rightRefSeg->templateID;
            be.rightSegment.isReverse = rightRefSeg->isReverse;
            be.rightSegment.endPos = be.rightSegment.beginPos + 1;

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
        analyzeRead(const_cast<TReadName&>(itRead->first));

        // release memory
        this->multChrmRead.erase( this->multChrmRead.find(const_cast<TReadName&>(itRead->first)) );
        itRead = this->recordByRead.erase(itRead);
    }

    return true;
}

void SplitRead::parseReadRecord(TReadName &readName, BamAlignmentRecord &record)
{
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
    if (this->recordByRead.find(readName) == this->recordByRead.end())
    {
        this->recordByRead[readName] = TBamRecords();
        this->multChrmRead[readName] = false;
    }
    this->recordByRead[readName].push_back(std::move(record));
}

void SplitRead::prepAfterHeaderParsing(BamHeader& header, BamFileIn& fileIn)
{
    this->fileIn = &fileIn;
}