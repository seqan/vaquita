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
#include <seqan/sequence.h>
#include "vaquita.hpp"
#include "option.hpp"
#include "calloption.hpp"
#include "mergeoption.hpp"
#include "misc.hpp"
#include "sv.hpp"
#include "svmerge.hpp"

#include <seqan3/argument_parser/all.hpp>

#include <sviper/sviper.h>
#include <sviper/io.h>
#include <sviper/auxiliary.h>

int callMain(int argc, char const ** argv)
{
    time_t startTime, endTime;
    bool result;

    // Get options
    CallOptionManager oMgr("call", argc, argv);
    oMgr.init_options();
    if(oMgr.parseCommandLine() != 0)
    {
        return 2;
    }

    bool doShortReads{!oMgr.getInputFile().empty()};
    bool doLongReads{!oMgr.getInputFile(true).empty()};

    // Init.
    AlignmentManager alnMgr(oMgr, false);
    AlignmentManager alnMgrLR(oMgr, true);
    BreakpointManager bpMgr(alnMgr, false);
    BreakpointManager bpMgrLR(alnMgrLR, true);

    // Start
    time(&startTime);
    oMgr.printUserInput();

    // Loading & extraction
    RUN_IF(doShortReads, result, "SHORT READ EVIDENCE EXTRACTION", alnMgr.load()); // TODO: segmentation fault if fails (need fix)
    if (!result && !doLongReads) return 3;
    RUN_IF(doLongReads, result, "LONG READ EVIDENCE EXTRACTION", alnMgrLR.load());
    if (!result) return 3;
    printTimeMessage("Found evidences");
    if (doShortReads)
    {
        printTimeMessage(std::to_string(alnMgr.getSplitReadCount()) + " from short split-reads.");
        printTimeMessage(std::to_string(alnMgr.getPairedReadCount()) + " from short discordant read-pairs.");
        printTimeMessage(std::to_string(alnMgr.getClippedReadCount()) + " from short soft-clipped reads.");
    }
    if (doLongReads)
    {
        printTimeMessage(std::to_string(alnMgrLR.getSplitReadCount()) + " from long split-reads.");
        printTimeMessage(std::to_string(alnMgrLR.getClippedReadCount()) + " from long soft-clipped reads.");
    }

    // Identification
    RUN_IF(doShortReads, result, "SHORT READ BREAKPOINT IDENTIFICATION", bpMgr.find());
    RUN_IF(doLongReads, result, "LONG READ BREAKPOINT IDENTIFICATION", bpMgrLR.find());
    if (!result) return 3;
    printTimeMessage("Found breakpoints");
    if (doShortReads)
    {
        printTimeMessage(std::to_string(bpMgr.getSplitRead()->getBreakpointCount())  + " from short split-read evidences.");
        printTimeMessage(std::to_string(bpMgr.getPairedEndRead()->getBreakpointCount())  + " from short read-pair evidences.");
        printTimeMessage(std::to_string(bpMgr.getClippedRead()->getBreakpointCount())  + " from short soft-clipped evidences.");
    }
    if (doLongReads)
    {
        printTimeMessage(std::to_string(bpMgrLR.getSplitRead()->getBreakpointCount())  + " from long split-read evidences.");
        printTimeMessage(std::to_string(bpMgrLR.getClippedRead()->getBreakpointCount())  + " from long soft-clipped evidences.");
    }

    // Combine long read breakpoints into short reads.
    RUN_IF((doShortReads && doLongReads), result, "COMBINING LONG AND SHORT READ BREAKPOINTS", bpMgr.addLongBP(bpMgrLR));
    // Merging
    RUN_IF(doShortReads, result, "BREAKPOINT MERGING", bpMgr.merge());
    if (!result) return 3;
    RUN_IF(!doShortReads, result, "BREAKPOINT MERGING", bpMgrLR.merge());

    BreakpointManager & finalBpMgr = doShortReads ? bpMgr : bpMgrLR;

    // #1
    // RUN(result,"BREAKPOINT MERGING", bpMgrShortRead.mergeFromLongRead(bpMgrLongRead));

    // #2
    // BreakpointManager combinedBpMgr(bpMgrShortRead, bpMgrLongRead);
    SVManager svMgr(finalBpMgr);

    if (!result) return 3;
    printTimeMessage("Breakpoints after merging: " + std::to_string(finalBpMgr.getMergedBreakpoint()->getBreakpointCount()));

    // Filtering
    RUN(result,"BREAKPOINT FILTERING", finalBpMgr.applyFilter());
    if (!result) return 3;
    int bpCnt = finalBpMgr.getMergedBreakpoint()->getBreakpointCount() - finalBpMgr.getMergedBreakpoint()->getFilteredBreakpointCount();
    printTimeMessage("Breakpoints after filtering: " + std::to_string(bpCnt));

    // SV Classification
    RUN(result,"SV CLASSIFICATION", svMgr.findSV());
    if (!result) return 3;
    printTimeMessage("Found SVs");
    printTimeMessage(std::to_string(svMgr.getDeletionCount()) + " deletions.");
    printTimeMessage(std::to_string(svMgr.getInversionCount()) + " inversions.");
    printTimeMessage(std::to_string(svMgr.getDuplicationCount()) + " duplications.");
    printTimeMessage(std::to_string(svMgr.getTranslocationCount()) + " translocations.");
    //printTimeMessage(std::to_string(svMgr.getBreakendCount()) + " breakends.");

    // SV ordering
    RUN(result,"SV PRIORITIZATION", svMgr.orderSV());

    if (doShortReads && doLongReads)
    {
        sviper::CmdOptions options(oMgr.getThreadCount(),
                                   oMgr.getInputFile(true),
                                   oMgr.getInputFile(),
                                   oMgr.getOutputFile(),
                                   oMgr.getReferenceGenome(),
                                   oMgr.getOutputFile() + "_polished");
        sviper::input_output_information info{options};

        // Check files
        if (!sviper::open_file_success(info.long_read_bai, (info.cmd_options.long_read_file_name + ".bai").c_str()) ||
            !sviper::open_file_success(info.short_read_bai, (info.cmd_options.short_read_file_name + ".bai").c_str()) ||
            !sviper::open_file_success(info.log_file, (info.cmd_options.output_prefix + ".log").c_str()))
            return 1;

        // Prepare file hangles for parallel computing
        // -------------------------------------------------------------------------
        RUN(result, "SViper: PREP PARALLEL PROCESSING", sviper::prep_file_handles(info));

        // Polish variants
        // -------------------------------------------------------------------------
        info.log_file  << "======================================================================" << std::endl
                        << "START polishing variants in of file " << info.cmd_options.candidate_file_name << std::endl
                        << "======================================================================" << std::endl;

        std::vector<VcfRecordEnhanced> & deletions_sviper = svMgr.getVcfBySVType()[svMgr.SVTYPE_DELETION()];
        startTimeMessage("SViper: POLISHING VARIANTS");
        #pragma omp parallel for schedule(guided)
        for (unsigned vidx = 0; vidx < deletions_sviper.size(); ++vidx)
        {
            if (oMgr.getReportFilteredResult() == true || deletions_sviper[vidx].status == VcfRecordEnhanced::STATUS::PASS)
            {
                sviper::Variant tmp{deletions_sviper[vidx]};
                tmp.ref_chrom = std::string{seqan::toCString(alnMgr.getRefName(std::stoi(tmp.ref_chrom)))};
                sviper::polish_variant(tmp, info);
                deletions_sviper[vidx].update(tmp);
            }
        } // parallel for loop
        endTimeMessage("SViper: POLISHING VARIANTS");

        info.log_file  << "======================================================================" << std::endl
                  << "                                 DONE"  << std::endl
                  << "======================================================================" << std::endl;
    }

    // Result
    RUN(result, "WRITE RESULT", svMgr.writeVCF());
    if (oMgr.getWriteBreakpoint()) finalBpMgr.writeBreakpoint();

    time(&endTime);
    printMessage("Total elapsed time : " + std::to_string((int)difftime(endTime,startTime)) + " seconds.");

    return 0;
}

int callMerge(int argc, char const ** argv)
{
    time_t startTime, endTime;
    bool result;

    // Get options
    MergeOptionManager oMgr("merge", argc, argv);
    oMgr.init_options();
    if (oMgr.parseCommandLine() != 0)
        return 1;

    // Init.
    SVMerge svMerge(oMgr);
    RUN(result, "LOADING VCFs", svMerge.loadVcfs());
    RUN(result, "MERGING VCFs", svMerge.merge());
    RUN(result, "WRITING VCF", svMerge.writeVCF());

    return 0;
}

int main(int argc, char const ** argv)
{
    bool callOrMerge = false;
    if (argc > 1)
    {
        std::string cmd(argv[1]);
        if (cmd == "call")
        {
            callMain(argc-1, argv+1);
            callOrMerge = true;
        }
        else if (cmd == "merge")
        {
            callMerge(argc-1, argv+1);
            callOrMerge = true;
        }
    }
    if (!callOrMerge)
    {
        // Either not enough arguments passed, or first argument was not call/merge.
        OptionManager oMgr("vaquita", argc, argv);
        oMgr.init_options();
        if(oMgr.parseCommandLine() != 0)
        {
            return 1;
        }
    }

    return 0;
}
