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
#include "mergedcandidate.hpp"

void MergedCandidate::doAdditionalJobAfterMerge(Breakpoint* destBp, Breakpoint* srcBp)
{
    ReadSupportInfo *destReadInfo, *srcReadInfo;

    destReadInfo = this->getReadSupport(destBp);
    srcReadInfo = this->getReadSupport(srcBp);

    destReadInfo->splitReadSupport += srcReadInfo->splitReadSupport;
    destReadInfo->pairedEndSupport += srcReadInfo->pairedEndSupport;
    destReadInfo->clippedReadSupport += srcReadInfo->clippedReadSupport;
    destReadInfo->leftReadDepthDiffScore = srcReadInfo->leftReadDepthDiffScore;
    destReadInfo->rightReadDepthDiffScore = srcReadInfo->rightReadDepthDiffScore;
    destReadInfo->leftReadDepth = srcReadInfo->leftReadDepth;
    destReadInfo->rightReadDepth = srcReadInfo->rightReadDepth;
}

ReadSupportInfo* MergedCandidate::initReadSupport(Breakpoint* bp)
{
    this->readSupportInfo[bp].splitReadSupport = 0;
    this->readSupportInfo[bp].pairedEndSupport = 0;
    this->readSupportInfo[bp].clippedReadSupport = 0;

    return getReadSupport(bp);
}

void MergedCandidate::setReadSupport(Breakpoint* bp, ReadSupportInfo& i)
{
    initReadSupport(bp);
    addReadSupport(bp, i);
}

void MergedCandidate::addReadSupport(Breakpoint* bp, ReadSupportInfo& i)
{
    this->readSupportInfo[bp].splitReadSupport = i.splitReadSupport;
    this->readSupportInfo[bp].pairedEndSupport = i.pairedEndSupport;
    this->readSupportInfo[bp].clippedReadSupport = i.clippedReadSupport;
    this->readSupportInfo[bp].leftReadDepthDiffScore = i.leftReadDepthDiffScore;
    this->readSupportInfo[bp].rightReadDepthDiffScore = i.rightReadDepthDiffScore;
    this->readSupportInfo[bp].leftReadDepth = i.leftReadDepth;
    this->readSupportInfo[bp].rightReadDepth = i.rightReadDepth;
}

ReadSupportInfo* MergedCandidate::getReadSupport(Breakpoint* bp)
{
    return &(this->readSupportInfo[bp]);
}