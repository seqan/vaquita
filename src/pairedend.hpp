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
#ifndef APP_PAIRED_END_READ_H_
#define APP_PAIRED_END_READ_H_

#include <map>
#include "candidate.hpp"

class PairedEndRead : public BreakpointCandidate
{		
	private:
		std::map<TReadName, BreakpointEvidence> pairInfo;
		std::map<Breakpoint*, bool> breakpointIsUsedInfo;
        BreakpointCandidate tempBreakpoints;

		void getCandidateRegion(TPosition&, TPosition&, BreakpointCandidate::SIDE);

    public:
        void prepAfterHeaderParsing(BamHeader& header, BamFileIn& fileIn)
        {
            this->tempBreakpoints.setOptionManager(this->getOptionManager()); 
            //this->setPositionalAdj(this->getOptionManager()->getPairedEndSearchSize()); 
        }
        bool analyze(void);
        void parseReadRecord(TReadName&, BamAlignmentRecord&);
		
        Breakpoint* updateWithCandidateRegion(BreakpointEvidence&, bool);
        bool isBreakpointUsed(Breakpoint* bp) { return this->breakpointIsUsedInfo[bp]; }
        void setBreakpointUsed(Breakpoint* bp, bool b) { this->breakpointIsUsedInfo[bp] = b; }
};
#endif // APP_PAIRED_END_READ_H_
