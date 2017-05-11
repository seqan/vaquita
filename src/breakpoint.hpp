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
#ifndef APP_BREAKPOINT_H_
#define APP_BREAKPOINT_H_

#include "option.hpp"
#include "alignment.hpp"
#include "splitread.hpp"
#include "pairedend.hpp"
#include "clippedread.hpp"
#include "readdepth.hpp"
#include "mergedcandidate.hpp"

struct FinalBreakpointInfo
{
    TTemplateID leftTemplateID;
    TPosition leftPosition;
    bool isLeftReverse;

    TTemplateID rightTemplateID;
    TPosition rightPosition;
    bool isRightReverse;

    double score = 0;
    double gcContent = 0.0;
    double sequenceComplexity = 0.0;
    char vote = 0;

    bool filtered = false;
    bool imprecise = false;
};

typedef std::map<CharString, unsigned>  TKmerSet;
class BreakpointManager
{
    private :
        OptionManager*      optionManager;
        AlignmentManager*   alignmentManager;
        SplitRead           splitReadBreakpoints;
        PairedEndRead       pairedEndBreakpoints;
        ClippedRead         clippedBreakpoints;
        ReadDepth           readDepthBreakpoints;
        MergedCandidate     mergedBreakpoints;

        std::map<Breakpoint*, FinalBreakpointInfo> finalBreakpoints;
        double averageReadDepth = 0.0;

        bool mergeSplitRead();
        bool mergePairedEnd();
        bool mergeClippedRead();
        bool findFinalBreakpoints(void);

        bool filterByEvidenceSumAndVote(void);
        
        bool priByRankAgg(void);
        bool priByEvidenceSum(void);

        bool calculateReadDepth();
        bool addImpreciseBreakpoints();
        void addNewPositionsByClippedSequence(TFoundPosition&, Breakpoint*, bool isLeftClip);
       
    public :
    	BreakpointManager(AlignmentManager& aln) { init(aln); }
        ~BreakpointManager() {};

    	void init(AlignmentManager&);
        bool merge(void);
        bool find(void);
        bool applyFilter(void);
        bool applyPrioritization(void);
        bool applyNormalization(void);
        bool rescueByCombinedEvidence(void);
        void writeBreakpoint(void);
        bool getSequenceFeature(void);
        void getNTCount(CharString&, unsigned&, unsigned&, unsigned&, unsigned&);

        AlignmentManager* getAlignmentManager(void) { return this->alignmentManager; }
        OptionManager* getOptionManager(void) { return this->optionManager; }
        BreakpointCandidate* getSplitRead(void) { return &this->splitReadBreakpoints; }
        BreakpointCandidate* getPairedEndRead(void) { return &this->pairedEndBreakpoints; }
        BreakpointCandidate* getClippedRead(void) { return &this->clippedBreakpoints; }
        BreakpointCandidate* getReadDepth(void) { return &this->readDepthBreakpoints; }
        MergedCandidate* getMergedBreakpoint(void) { return &this->mergedBreakpoints; }
        FinalBreakpointInfo* getFinalBreakpointInfo(Breakpoint* bp) { return &this->finalBreakpoints[bp]; }
        double getReTH(void) { return readDepthBreakpoints.getReTH(); }
        double getDepthTH(void) { return readDepthBreakpoints.getDepthTH(); }
        double getDepthMedian(void) { return readDepthBreakpoints.getMedianDepth(); }
};

#endif // APP_BREAKPOINT_H_
