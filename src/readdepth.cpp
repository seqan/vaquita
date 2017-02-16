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
#include "readdepth.hpp"

void ReadDepth::parseReadRecord(CharString &id, BamAlignmentRecord &record)
{  
    std::transform(readDepthInfo[record.rID].begin() + record.beginPos, \
                   readDepthInfo[record.rID].begin() + record.beginPos + getAlignmentLengthInRef(record),
                   readDepthInfo[record.rID].begin() + record.beginPos,
                   std::bind2nd(std::plus<TPosition>(), 1));
}


void ReadDepth::prepAfterHeaderParsing(BamHeader& header, BamFileIn& fileIn)
{
    for (unsigned i=0; i < length(header); ++i)
    {
        if (header[i].type == BamHeaderRecordType::BAM_HEADER_REFERENCE)
        {
            CharString templateName, templateLengthStr;
            getTagValue(templateName, "SN", header[i]);
            getTagValue(templateLengthStr, "LN", header[i]);
            toCString(templateLengthStr);
               
            TTemplateID rID = BreakpointEvidence::NOVEL_TEMPLATE;
            TPosition templateLength = BreakpointEvidence::INVALID_POS;
            getIdByName(rID, contigNamesCache(context(fileIn)), templateName);
            templateLength = std::stol( std::string(toCString(templateLengthStr)) );

            readDepthInfo[rID].resize(templateLength);
            std::fill(readDepthInfo[rID].begin(), readDepthInfo[rID].end(), 0);
        }
    }
}

void ReadDepth::addUniformDepth(TTemplateID id, TPosition beginPos, TPosition size, unsigned depth)
{  
    int endPos = beginPos + size;
    if ( beginPos > readDepthInfo[id].size() )
        beginPos = readDepthInfo[id].size() - 1;
    if ( endPos > readDepthInfo[id].size() )
        endPos = readDepthInfo[id].size() - 1;
   
    std::transform(readDepthInfo[id].begin() + beginPos, \
                   readDepthInfo[id].begin() + endPos,
                   readDepthInfo[id].begin() + beginPos,
                   std::bind2nd(std::plus<TPosition>(), depth));

}

double ReadDepth::printDepth(TTemplateID templateID, TPosition position, TPosition windowSize)
{
    for (TPosition pos = position; pos < position+windowSize; ++pos)
    {
        std::cerr << pos << "\t";
        std::cerr <<  unsigned(readDepthInfo[templateID][pos]) << "\n";
    }
}

void ReadDepth::getReadDepthDiffScore(double& leftOuterDepthAvg, double& rightOuterDepthAvg, double& leftScore, double& rightScore, TTemplateID t1, TPosition p1, TTemplateID t2, TPosition p2, TPosition windowSize, TPosition breakpointSize)
{
    double innerDepthAvg, outerDepthAvg, leftDepthAvg, rightDepthAvg;

    // left
    getAvgReadDepth(leftDepthAvg, rightDepthAvg, t1, p1, windowSize, breakpointSize, BreakpointCandidate::SIDE::LEFT);

    leftOuterDepthAvg = leftDepthAvg;
    outerDepthAvg = leftDepthAvg;
    innerDepthAvg = rightDepthAvg;
    if (leftDepthAvg > rightDepthAvg)
        leftScore = getKLScore(leftDepthAvg, rightDepthAvg);
    else
        leftScore = getKLScore(rightDepthAvg, leftDepthAvg);

    // right
    getAvgReadDepth(leftDepthAvg, rightDepthAvg, t2, p2, windowSize, breakpointSize, BreakpointCandidate::SIDE::RIGHT);
    rightOuterDepthAvg = rightDepthAvg;
    innerDepthAvg += leftDepthAvg;
    outerDepthAvg += rightDepthAvg;
    if (leftDepthAvg > rightDepthAvg)
        rightScore = getKLScore(leftDepthAvg, rightDepthAvg); // out -> in (del)
    else 
        rightScore = getKLScore(rightDepthAvg, leftDepthAvg); // in -> out (dup)

    // avg
    innerDepthAvg = innerDepthAvg / 2.0;
    outerDepthAvg = outerDepthAvg / 2.0;
}

void ReadDepth::getAvgReadDepth(double& leftAvg, double& rightAvg, TTemplateID templateID, TPosition position, TPosition windowSize, TPosition breakpointSize, BreakpointCandidate::SIDE side)
{
    if (templateID == BreakpointEvidence::NOVEL_TEMPLATE || readDepthInfo[templateID].size() == 0)
    {
       leftAvg = 0.0;
       rightAvg = 0.0;
       return;
    }

    if (position > readDepthInfo[templateID].size() - 1)
        position = readDepthInfo[templateID].size() - 1;

    TPosition begin, end;
    if (side == BreakpointCandidate::SIDE::LEFT)    
    {
        if (position > windowSize)
            begin = position - windowSize;
        else
            begin = 0;

        if (position < readDepthInfo[templateID].size() - breakpointSize - 1)
            end = position + breakpointSize + 1;
        else
            end = readDepthInfo[templateID].size() - 1;
    }
    else
    {
        if (position > breakpointSize)
            begin = position - breakpointSize;
        else
            begin = 0;

        if (position < readDepthInfo[templateID].size() - windowSize - 1)
            end = position + windowSize + 1;
        else
            end = readDepthInfo[templateID].size() - 1;
    }

    // sum
    leftAvg = std::accumulate(readDepthInfo[templateID].begin() + begin, readDepthInfo[templateID].begin() + position, 0.0);
    rightAvg = std::accumulate(readDepthInfo[templateID].begin() + position, readDepthInfo[templateID].begin() + end, 0.0);

    // avg
    leftAvg = leftAvg / ((position - begin) + PREVENT_DIV_BY_ZERO());
    rightAvg = rightAvg / ((end - position) + PREVENT_DIV_BY_ZERO());
}

double ReadDepth::getKLScore(double p0, double p1)
{
    // KL from p0 to p1
    if ( (KLTable.find(p0) == KLTable.end()) || (KLTable[p0].find(p1) == KLTable[p0].end()) )
    {
        double _p0 = p0 + PREVENT_DIV_BY_ZERO();
        double _p1 = p1 + PREVENT_DIV_BY_ZERO();
        KLTable[p0][p1] = _p0 - _p1 + (_p1 * log(_p1/_p0));
    }

    return KLTable[p0][p1];
}
