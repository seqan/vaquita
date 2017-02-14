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

    // start
    printMessage("======================================================="); 
    time(&startTime);

    // file loading
    RUN(result, "File Loading", alnMgr.load()); // segmentation fault if fails
    
    printMessage("Split-read evidences : " + std::to_string(alnMgr.getSplitReadCount()));
    printMessage("Paired-ends evidences : " + std::to_string(alnMgr.getPairedReadCount()));
    printMessage("Clipped-read evidences : " + std::to_string(alnMgr.getClippedReadCount()));
  
    // find breakpoints (TODO: messages)
    RUN(result,"Find breakpoints", bpMgr.find());

    // merge breakpoints (TODO: messages)
    RUN(result,"Merge Breakpoints", bpMgr.merge());

    // filter breakpoints (TODO: messages)
    RUN(result,"Filtering Breakpoints", bpMgr.applyFilter());
    //RUN(result,"Filtering", bpMgr.filterByScore());
    if (oMgr.getWriteBreakpoint() == true)
        RUN(result,"Write Breakpoints", bpMgr.writeBreakpoint());

    // classification
    RUN(result,"Find structurial variations", svMgr.findSV());

    // write output
    RUN(result, "Write result", svMgr.writeSV());

    time(&endTime);
    printMessage("======================================================="); 
    printMessage("Elapsed time : " + std::to_string((int)difftime(endTime,startTime)) + " seconds.");

    return 0;
}