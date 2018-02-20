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
#include "readdepth.hpp"
#include "misc.hpp"
#include <cmath>

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

void ReadDepth::setRandomSeed(int seed)
{
    srand(seed);
}

void ReadDepth::getRandomPos(TTemplateID& t, TPosition &p)
{
    // select chromosome
    t = rand() % readDepthInfo.size();
    p = rand() % readDepthInfo[t].size();
}

bool ReadDepth::getRandomPosByTemplate(TTemplateID t, TPosition &p)
{
    if (t < readDepthInfo.size())
    {
        p = rand() % readDepthInfo[t].size();
        return true;
    }
    else
    {
        p = BreakpointEvidence::INVALID_POS;
        return false;
    }
}

void ReadDepth::calculateReadDepthStat(std::map<TTemplateID, unsigned>& svCntByTemplate, unsigned totalCnt)
{
    // get RE distribution
    std::vector<double> backgroundRE;
    std::vector<double> backgroundDepth;
    this->setRandomSeed(0);

    for (auto itSVCntByTemp = svCntByTemplate.begin(); itSVCntByTemp != svCntByTemplate.end(); ++itSVCntByTemp)
    {
        TTemplateID t = itSVCntByTemp->first;

        if (this->readDepthInfo.find(t) != this->readDepthInfo.end())
        {
            unsigned sampleCnt = ceil((double)this->getOptionManager()->getSamplingNum() * ((double) itSVCntByTemp->second / (double)totalCnt));
            for(unsigned i=0; i < sampleCnt; ++i)
            {
                TPosition p;
                double kl = 0.0, depth = 0.0;

                // get kl and depth at random position
                while (depth <= 0.0)
                {
                    if (this->getRandomPosByTemplate(t, p))
                        this->getReadDepthDiffScore(kl, depth, t, p, this->getOptionManager()->getReadDepthWindowSize());
                    else
                        break;
                }
                backgroundRE.push_back(kl);
                backgroundDepth.push_back(depth);
            }
        }
    }
    std::sort(backgroundRE.begin(), backgroundRE.end(), [](auto &left, auto &right) {return left < right;} );
    std::sort(backgroundDepth.begin(), backgroundDepth.end(), [](auto &left, auto &right) {return left < right;} );

    // Median depth
    double medianDepth = MID_ELEMENT(backgroundDepth);
    double q1 = backgroundDepth[backgroundDepth.size()*0.25];
    double q3 = backgroundDepth[backgroundDepth.size()*0.75];
    this->depthTH = q3 + (double)(q3-q1) * (double)this->getOptionManager()->getDepthOutlier();

    this->setMedianDepth( medianDepth );
    printTimeMessage("Depth median: " + std::to_string(medianDepth));
    printTimeMessage("Depth Q1: " + std::to_string(q1));
    printTimeMessage("Depth Q3: " + std::to_string(q3));
    printTimeMessage("Depth IQR: " + std::to_string(q3-q1));
    printTimeMessage("Depth threshold: " + std::to_string(this->depthTH));
    /*
    for (auto it = backgroundDepth.begin(); it != backgroundDepth.end(); ++it)
    {
        *it -= medianDepth;
        if (*it < 0.0) *it *= -1.0;
    }
    std::sort(backgroundDepth.begin(), backgroundDepth.end(), [](auto &left, auto &right) {return left < right;} );
    double madDepth = MID_ELEMENT(backgroundDepth);   
    //printTimeMessage("Depth median absolute deviation: " + std::to_string(madDepth));
    */
    
    // Read-depth evidence
    double medianRE = 0.0, madRE = 0.0;
    medianRE = MID_ELEMENT(backgroundRE);
    q1 = backgroundRE[backgroundRE.size()*0.25];
    q3 = backgroundRE[backgroundRE.size()*0.75];
    this->reTH = q3 + (double)(q3-q1) * (double)this->getOptionManager()->getReOutlierCutoff();
    printTimeMessage("RE median: " + std::to_string(medianRE));
    printTimeMessage("RE Q1: " + std::to_string(q1));
    printTimeMessage("RE Q3: " + std::to_string(q3));
    printTimeMessage("RE IQR: " + std::to_string(q3-q1));
    printTimeMessage("RE threshold : " + std::to_string(this->reTH));
    /*
    for (auto it = backgroundRE.begin(); it != backgroundRE.end(); ++it)
    {
        *it -= medianRE;
        if (*it < 0.0) *it *= -1.0;
    }
    std::sort(backgroundRE.begin(), backgroundRE.end(), [](auto &left, auto &right) {return left < right;} );
    madRE = MID_ELEMENT(backgroundRE);
    //this->reTH = medianRE + madRE * this->getOptionManager()->getReOutlierCutoff();
    //printTimeMessage("RE median absoulte deviation: " + std::to_string(madRE));
    */


    // KL figure
    /*
    std::cerr << "KL figure start\n";
    for (auto it = backgroundRE.begin(); it != backgroundRE.end(); ++it )
    {
        std::cout << "RANDOM\t" << *it << std::endl;
    }
    std::ifstream file("tp.txt");
    std::string sv_type;
    int chrm, start, end;
    while(file >> sv_type >> chrm >> start >> end)
    {
        double l_kl, r_kl, kl, l_depth, r_depth, depth;
        int32_t breakpointSize = abs(end - start) + 1;
        int32_t windowSize = this->getOptionManager()->getReadDepthWindowSize();
        if (windowSize > breakpointSize)
            windowSize = breakpointSize;

        this->getReadDepthDiffScore(l_depth, r_depth, l_kl, r_kl, chrm, start, chrm, end, windowSize, windowSize);
        kl = std::max(l_kl, r_kl);
        depth = (l_depth + r_depth) / 2.0;
        std::cout << sv_type << "\t" << start << "\t" << kl << "\t" << depth << "\n";
    }
    std::cerr << "KL figure end\n";
    exit(1);
    */
}

void ReadDepth::addUniformDepth(TTemplateID id, TPosition beginPos, TPosition size, unsigned depth)
{  
    int endPos = beginPos + size;
    if ( beginPos >= readDepthInfo[id].size() )
        beginPos = readDepthInfo[id].size() - 1;
    if ( endPos >= readDepthInfo[id].size() )
        endPos = readDepthInfo[id].size() - 1;
   
    std::transform(readDepthInfo[id].begin() + beginPos, \
                   readDepthInfo[id].begin() + endPos,
                   readDepthInfo[id].begin() + beginPos,
                   std::bind2nd(std::plus<TPosition>(), depth));

    this->baseCount += (size * depth);
}

void ReadDepth::printDepth(TTemplateID templateID, TPosition position, TPosition windowSize)
{
    for (TPosition pos = position; pos < position+windowSize; ++pos)
    {
        std::cerr << pos << "\t";
        std::cerr <<  unsigned(readDepthInfo[templateID][pos]) << "\n";
    }
}

void ReadDepth::getReadDepthDiffScore(double& score, double& depth, TTemplateID t, TPosition p, TPosition windowSize)
{
    double leftDepthAvg, rightDepthAvg;

    getAvgReadDepth(leftDepthAvg, rightDepthAvg, t, p, windowSize, windowSize, BreakpointEvidence::SIDE::LEFT);
    if (leftDepthAvg > rightDepthAvg)
        score = getKLScore(leftDepthAvg, rightDepthAvg); // / leftDepthAvg;
    else
        score = getKLScore(rightDepthAvg, leftDepthAvg); // / rightDepthAvg;
    depth = (leftDepthAvg + rightDepthAvg) / 2.0; // used to estimate the sample depth
}

void ReadDepth::getReadDepthDiffScore(double& leftOuterDepthAvg, double& rightOuterDepthAvg, double& leftScore, double& rightScore, TTemplateID t1, TPosition p1, TTemplateID t2, TPosition p2, TPosition windowSize, TPosition breakpointSize)
{
    double outerDepthAvg, leftDepthAvg, rightDepthAvg;

    // left
    getAvgReadDepth(leftDepthAvg, rightDepthAvg, t1, p1, windowSize, breakpointSize, BreakpointEvidence::SIDE::LEFT);
    leftOuterDepthAvg = leftDepthAvg;
    if (leftDepthAvg > rightDepthAvg)
        leftScore = getKLScore(leftDepthAvg, rightDepthAvg); // / leftDepthAvg;
    else
        leftScore = getKLScore(rightDepthAvg, leftDepthAvg); // / rightDepthAvg;

    // right
    getAvgReadDepth(leftDepthAvg, rightDepthAvg, t2, p2, windowSize, breakpointSize, BreakpointEvidence::SIDE::RIGHT);
    rightOuterDepthAvg = rightDepthAvg;
    if (leftDepthAvg > rightDepthAvg)
        rightScore = getKLScore(leftDepthAvg, rightDepthAvg); // / leftDepthAvg;
    else
        rightScore = getKLScore(rightDepthAvg, leftDepthAvg); // / rightDepthAvg;
}

void ReadDepth::getAvgReadDepth(double& leftAvg, double& rightAvg, TTemplateID templateID, TPosition position, TPosition windowSize, TPosition breakpointSize, BreakpointEvidence::SIDE side)
{
    if (templateID == BreakpointEvidence::NOVEL_TEMPLATE || readDepthInfo[templateID].size() == 0 || windowSize == 0 || breakpointSize == 0)
    {
       leftAvg = 0.0;
       rightAvg = 0.0;
       return;
    }

    if (position > readDepthInfo[templateID].size() - 1)
        position = readDepthInfo[templateID].size() - 1;

    TPosition begin, end;
    if (side == BreakpointEvidence::SIDE::LEFT)    
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
        if (KLTable[p0][p1] < 0.0)
            KLTable[p0][p1] = 0.0;
    }
    return KLTable[p0][p1];
}

double ReadDepth::getPoissonP(unsigned k, double l)
{
    return (pow(l,k) * exp(-l)) / std::tgamma(k+1);
}