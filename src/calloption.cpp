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
#include <climits>
#include <cfloat>

#include "vaquita.hpp"
#include "calloption.hpp"
#include "misc.hpp"

// ============================================================================
// Functions
// ============================================================================

void CallOptionManager::init_options()
{
    seqan3::arithmetic_range_validator int_val{0, INT32_MAX};
    seqan3::arithmetic_range_validator dbl_val{0, DBL_MAX};

    (*this).info.short_description = "Identification mode";

    // version & date
    (*this).info.version = SEQAN_APP_VERSION;
    (*this).info.date = SEQAN_DATE;

    // description
    (*this).info.description.push_back(std::string(APP_NAME) + std::string(": ") + std::string(APP_TITLE));
    (*this).info.description.push_back(std::string(APP_AUTHOR_INFO));
    (*this).info.url = std::string(APP_WEBSITE_INFO);

    // usage line
    (*this).info.synopsis.push_back("call [\\fIOPTIONS\\fP] -r [\\fIreference.fa\\fP] [\\fIalignment.bam\\fP] > [\\fIout.vcf\\fP]");

    // mandatory arguments
    (*this).add_positional_option(inputFile, "ALIGNMENT(.bam)", seqan3::input_file_validator({std::string("bam")}));

    // Options
    (*this).add_section("General");
    (*this).add_option(referenceGenome, 'r', "referenceGenome", "Genome sequence file(.fa).", seqan3::option_spec::REQUIRED, seqan3::input_file_validator({std::string("fa")}));
    (*this).add_option(cutoff, 'c', "cutoff", "Minimum number of supporting read-pairs and split-reads.", seqan3::option_spec::DEFAULT, int_val);
    (*this).add_option(minVote, 'v', "minVote", "Minimum number of evidence types(=vote) that support SVs for rescue. -1: Supported by all evidence types.", seqan3::option_spec::DEFAULT, seqan3::arithmetic_range_validator{-1, INT32_MAX});
    (*this).add_option(minMapQual, 'q', "minMapQual", "Mapping quaility cutoff.", seqan3::option_spec::DEFAULT, int_val);
    (*this).add_option(minSVSize, 'm', "minSVSize", "Structural varation size cutoff.", seqan3::option_spec::DEFAULT, int_val);
    (*this).add_option(adjTol, 'a', "adjTol", "Positional adjacency in nucleotide resolution.", seqan3::option_spec::DEFAULT, int_val);
    (*this).add_flag(reportFilteredResult, '\0', "report-filtered", "Report filtered result");
    (*this).add_flag(skipPairedEndRead, '\0', "no-pe", "Do not use read-pair evidence.");
    (*this).add_flag(skipClippedRead, '\0', "no-ce", "Do not use soft-clipped evidence.");
    (*this).add_flag(skipReadDepth, '\0', "no-re", "Do not use read-depth evidence.");
    (*this).add_flag(skipRankAggregation, '\0', "no-rank-aggregation", "Do not use rank-aggregation for prioritization.");
    (*this).add_flag(writeBreakpoint, 'w', "write-breakpoint", "Write breakpoint informtation in a tab-separated format (breakpoints.tsv).");

    (*this).add_section("Split-read evidence");
    (*this).add_option(maxSplit, 's', "maxSplit", "Maximum number of segments in a single read.", seqan3::option_spec::DEFAULT, int_val);
    (*this).add_option(maxOverlap, 'o', "maxOverlap", "Maximum allowed overlaps between segements.", seqan3::option_spec::DEFAULT, int_val);
    (*this).add_option(minSplitReadSupport, 'b', "minSplitReadSupport", "SVs supported by >= b get a vote.", seqan3::option_spec::DEFAULT, dbl_val);

    (*this).add_section("Read-pair evidence");
    (*this).add_option(pairedEndSearchSize, 'p', "pairedEndSearchSize", "Size of breakpoint candidate regions.", seqan3::option_spec::DEFAULT, int_val);
    (*this).add_option(abInsParam, 'i', "abInsParam", "Discordant insertion size: median +/- (MAD * i)", seqan3::option_spec::DEFAULT, dbl_val);
    (*this).add_option(depthOutlier, 'd', "depthOutlier", "Depth outlier: {Q3 + (IQR * d)}", seqan3::option_spec::DEFAULT, dbl_val);
    (*this).add_option(minPairSupport, 'e', "minPairSupport", "SVs supported by >= e get a vote.", seqan3::option_spec::DEFAULT, dbl_val);

    (*this).add_section("Soft-clipped evidence");
    (*this).add_option(minClippedSeqSize, 'l', "minClippedSeqSize", "Minimum size of clipped sequence to be considered.", seqan3::option_spec::DEFAULT, int_val);
    (*this).add_option(clippedSeqErrorRate, 't', "clippedSeqErrorRate", "Maximum edit distance: floor{length of clipped sequence * (1 - t)}.", seqan3::option_spec::DEFAULT, dbl_val);
    // //addOption(*this, ArgParseOption("", "use-assembly", "Use consitency based sequence assembly(deprecated)."));
    // //setDefaultValue(*this, "use-assembly", "false");

    (*this).add_section("Read-Depth evidence");
    (*this).add_option(samplingNum, 'n', "samplingNum", "Number of random sample to estimate the background distribution(Q3, IQR, ..) of read-depth evidence.", seqan3::option_spec::DEFAULT, int_val);
    (*this).add_option(readDepthWindowSize, 'f', "readDepthWindowSize", "Window size to caclulate average read-depth around breakpoints.", seqan3::option_spec::DEFAULT, int_val);
    (*this).add_flag(useREforBalancedSV, '\0', "use-re-for-bs", "Use RE for balanced SVs(eg. inverison).");
    (*this).add_option(reThreshold, 'j', "reThreshold", "SVs satisfy read-depth evidence >= {Q3 + (IQR * h)} get a vote.", seqan3::option_spec::DEFAULT, dbl_val);

    // output options
    /* not supported yet
    addSection(*this, "Output");
    addOption(*this, ArgParseOption("o", "output-file", "Specify an output file. Default: write the file to standard output.", ArgParseOption::OUTPUT_FILE));
    setValidValues(*this, "output-file", ".vcf");
    */
}

int CallOptionManager::parseCommandLine()
{
    try
    {
        (*this).parse();
    }
    catch(seqan3::parser_invalid_argument const & ext)
    {
        std::cerr << "[PARSER ERROR] " << ext.what() << "\n\n";
        const char * dummy_argv[] = {"./dummy", "-h"};
        CallOptionManager parser_dummy("call", 0, dummy_argv);
        parser_dummy.init_options();
        try
        {
            parser_dummy.parse();
        }
        catch(seqan3::parser_invalid_argument const &)
        {

        }
        return 1;
    }
    return 0;
}

std::string CallOptionManager::getBooleanString(bool b)
{
    if (b) return "True";
    else return "False";
}

void CallOptionManager::printUserInput(void)
{
    printMessage("==============================");
    printMessage(std::string(APP_NAME) + std::string(" ") + std::string(SEQAN_APP_VERSION));
    printMessage("==============================");
    printMessage("[Input File]");
    printMessage(this->inputFile);
    printMessage("[General options]");
    printMessage("- referenceGenome: " + this->referenceGenome);
    printMessage("- minMapQual: " + std::to_string(this->minMapQual));
    printMessage("- minSVSize: " + std::to_string(this->minSVSize));
    printMessage("- adjTol: " + std::to_string(this->adjTol));
    printMessage("- cutoff: " + std::to_string(this->cutoff));
    printMessage("- minVote: " + std::to_string(this->minVote));
    printMessage("- reportFilteredResult: " + getBooleanString(this->reportFilteredResult));
    printMessage("- no-pe: " + getBooleanString(this->skipPairedEndRead));
    printMessage("- no-ce: " + getBooleanString(this->skipClippedRead));
    printMessage("- no-re: " + getBooleanString(this->skipReadDepth));
    printMessage("- no-rank-aggregation: " + getBooleanString(this->skipRankAggregation));
    printMessage("[Split-read options]");
    printMessage("- maxSplit: " + std::to_string(this->maxSplit));
    printMessage("- maxOverlap: " + std::to_string(this->maxOverlap));
    printMessage("- minSplitReadSupport: " + std::to_string(this->minSplitReadSupport));
    printMessage("[Read-pair options]");
    printMessage("- abInsParam: " + std::to_string(this->abInsParam));
    printMessage("- depthOutlier: " + std::to_string(this->depthOutlier));
    printMessage("- pairedEndSearchSize: " + std::to_string(this->pairedEndSearchSize));
    printMessage("- minPairSupport: " + std::to_string(this->minPairSupport));
    printMessage("[Soft-clipped-read options]");
    printMessage("- minClippedSeqSize: " + std::to_string(this->minClippedSeqSize));
    printMessage("- clippedSeqErrorRate: " + std::to_string(this->clippedSeqErrorRate));
    //printMessage("- use-assembly: " + std::to_string(this->useAssembly));
    printMessage("[Read-depth options]");
    printMessage("- samplingNum: " + std::to_string(this->samplingNum));
    printMessage("- readDepthWindowSize: " + std::to_string(this->readDepthWindowSize));
    printMessage("- useREforBalancedSV: " + getBooleanString(this->useREforBalancedSV));
    printMessage("- reThreshold: " + std::to_string(this->reThreshold));
    printMessage("==============================");
}
