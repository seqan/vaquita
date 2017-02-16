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
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include "vaquita.hpp"
#include "misc.hpp"
#include "sv.hpp"

int main(int argc, char const ** argv)
{
    // get options
    OptionManager oMgr;
    time_t startTime, endTime;
    bool result;

    // init. (TODO: messages)
    oMgr.init();
    if ( !oMgr.parseCommandLine(argc, argv) ) 
        return 1;
    AlignmentManager alnMgr(oMgr);
    BreakpointManager bpMgr(alnMgr);
    SVManager svMgr(bpMgr);

    // Start
    time(&startTime);
    oMgr.printUserInput();

    // Loading & extraction
    RUN(result, "EVIDENCE EXTRACTION", alnMgr.load()); // segmentation fault if fails
    printTimeMessage("Found evidences");
    printTimeMessage(std::to_string(alnMgr.getSplitReadCount()) + " from split-reads.");
    printTimeMessage(std::to_string(alnMgr.getPairedReadCount()) + " from abnormal read-pairs.");
    printTimeMessage(std::to_string(alnMgr.getClippedReadCount()) + " from soft-clipped reads.");

    // Identification
    RUN(result,"BREAKPOINT IDENTIFICATION", bpMgr.find());
    printTimeMessage("Found breakpoints");
    printTimeMessage(std::to_string(bpMgr.getSplitRead()->getBreakpointCount())  + " from split-read evidences.");
    printTimeMessage(std::to_string(bpMgr.getPairedEndRead()->getBreakpointCount())  + " from read-pair evidences.");
    printTimeMessage(std::to_string(bpMgr.getClippedRead()->getBreakpointCount())  + " from soft-clipped evidences.");

    // Merging
    RUN(result,"BREAKPOINT MERGING", bpMgr.merge());
    printTimeMessage("Breakpoints after merging: " + std::to_string(bpMgr.getMergedBreakpoint()->getBreakpointCount()));

    // Filtering
    RUN(result,"BREAKPOINT FILTERING", bpMgr.applyFilter());
    if (oMgr.getWriteBreakpoint() == true)
        RUN(result,"Write Breakpoints", bpMgr.writeBreakpoint());
    printTimeMessage("Breakpoints after filtering: " + std::to_string(bpMgr.getMergedBreakpoint()->getBreakpointCount()));

    // SV Classificatin
    RUN(result,"SV IDENTIFICATION", svMgr.findSV());
    printTimeMessage("Found SVs");
    printTimeMessage(std::to_string(svMgr.getDeletionCount()) + " deletions (including pseudo-deletions).");
    printTimeMessage(std::to_string(svMgr.getInversionCount()) + " inversions.");
    printTimeMessage(std::to_string(svMgr.getDuplicationCount()) + " duplications.");
    printTimeMessage(std::to_string(svMgr.getTranslocationCount()) + " translocations.");
    printTimeMessage(std::to_string(svMgr.getBreakendCount()) + " breakends.");

    // Result
    RUN(result, "WRITE RESULT", svMgr.writeVCF());

    time(&endTime);
    printMessage("Total elapsed time : " + std::to_string((int)difftime(endTime,startTime)) + " seconds.");

    return 0;
}