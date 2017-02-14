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
#ifndef APP_MERGEDCANDIDATE_H_
#define APP_MERGEDCANDIDATE_H_

#include "candidate.hpp"

struct ReadSupportInfo
{
    static uint16_t const MAX_VAL = 99;
    uint16_t    splitReadSupport = 0;
    uint16_t    pairedEndSupport = 0;
    uint16_t    clippedReadSupport = 0;
    double      leftReadDepth = 0.0;
    double      rightReadDepth = 0.0;
    double      avgReadDepth = 0.0;
    double      leftReadDepthDiffScore = 0.0;
    double      rightReadDepthDiffScore = 0.0;
    bool        leftReadDepthSelected = true;
    bool        pseudoDeletion = false;
};

class MergedCandidate : public BreakpointCandidate
{    
    private:
        std::map<Breakpoint*, ReadSupportInfo> readSupportInfo;
    
    public:
        void doAdditionalJobAfterMerge(Breakpoint*, Breakpoint*);

        void setReadSupport(Breakpoint*, ReadSupportInfo&);
        void addReadSupport(Breakpoint*, ReadSupportInfo&);
        ReadSupportInfo* getReadSupport(Breakpoint*);
};

#endif // APP_MERGEDCANDIDATE_H_