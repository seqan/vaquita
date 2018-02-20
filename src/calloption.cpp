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
#include "vaquita.hpp"
#include "calloption.hpp"
#include "misc.hpp"

// ============================================================================
// Functions
// ============================================================================

void CallOptionManager::init()
{
    setAppName(*this, APP_NAME);
    setShortDescription(*this, "Identification mode");

    // version & date
    setVersion(*this, SEQAN_APP_VERSION);
    setDate(*this, SEQAN_DATE);

    // description
    addDescription(*this, std::string(APP_NAME) + std::string(": ") + std::string(APP_TITLE));
    addDescription(*this, std::string(APP_AUTHOR_INFO));
    addDescription(*this, std::string(APP_WEBSITE_INFO));

    // usage line
    addUsageLine(*this, "[\\fIOPTIONS\\fP] -r [\\fIreference.fa\\fP] [\\fIalignment.bam\\fP] > [\\fIout.vcf\\fP]");

    // mandatory arguments
    addArgument(*this, ArgParseArgument(ArgParseArgument::INPUT_FILE, "ALIGNMENT(.bam)")); 
    setValidValues(*this, 0, ".bam");
    setHelpText(*this, 0, "Support : .bam");

    // Options
    addSection(*this, "General");
    addOption(*this, ArgParseOption("r", "referenceGenome", "Genome sequence file(.fa).", ArgParseOption::INPUT_FILE));
    setValidValues(*this, "referenceGenome", ".fa");
    setRequired(*this, "r");
    addOption(*this, ArgParseOption("c", "cutoff", "Minimum number of supporting read-pairs and split-reads.", ArgParseOption::INTEGER, "INT"));
    setDefaultValue(*this, "cutoff", "4");
    addOption(*this, ArgParseOption("v", "minVote", "Minimum number of evidence types(=vote) that support SVs for rescue. -1: Supported by all evidence types.", ArgParseOption::INTEGER, "INT"));
    setDefaultValue(*this, "minVote", "-1");
    addOption(*this, ArgParseOption("q", "minMapQual", "Mapping quaility cutoff.", ArgParseOption::INTEGER, "INT"));
    setDefaultValue(*this, "minMapQual", "20");
    addOption(*this, ArgParseOption("m", "minSVSize", "Structural varation size cutoff.", ArgParseOption::INTEGER, "INT"));
    setDefaultValue(*this, "minSVSize", "50");
    addOption(*this, ArgParseOption("a", "adjTol", "Positional adjacency in nucleotide resolution.", ArgParseOption::INTEGER, "INT"));
    setDefaultValue(*this, "adjTol", "50");
    addOption(*this, ArgParseOption("", "report-filtered", "Report filtered result"));
    setDefaultValue(*this, "report-filtered", "false");
    addOption(*this, ArgParseOption("", "no-pe", "Do not use read-pair evidence."));
    setDefaultValue(*this, "no-pe", "false");
    addOption(*this, ArgParseOption("", "no-ce", "Do not use soft-clipped evidence."));
    setDefaultValue(*this, "no-ce", "false");
    addOption(*this, ArgParseOption("", "no-re", "Do not use read-depth evidence."));
    setDefaultValue(*this, "no-re", "false");
    addOption(*this, ArgParseOption("", "no-rank-aggregation", "Do not use rank-aggregation for prioritization."));
    setDefaultValue(*this, "no-rank-aggregation", "false");
    /* this option is for debugging
    addOption(*this, ArgParseOption("", "write-breakpoint", "Write breakpoint informtation in a tab-separated format (breakpoints.tsv)."));
    setDefaultValue(*this, "write-breakpoint", "false");
    */

    addSection(*this, "Split-read evidence");
    addOption(*this, ArgParseOption("ss", "maxSplit", "Maximum number of segments in a single read.", ArgParseOption::INTEGER, "INT"));
    setDefaultValue(*this, "maxSplit", "2");
    addOption(*this, ArgParseOption("so", "maxOverlap", "Maximum allowed overlaps between segements.", ArgParseOption::INTEGER, "INT"));
    setDefaultValue(*this, "maxOverlap", "20");
    addOption(*this, ArgParseOption("se", "minSplitReadSupport", "SVs supported by >= se get a vote.", ArgParseOption::DOUBLE, "FLOAT"));
    setDefaultValue(*this, "minSplitReadSupport", "1");

    addSection(*this, "Read-pair evidence");
    addOption(*this, ArgParseOption("ps", "pairedEndSearchSize", "Size of breakpoint candidate regions.", ArgParseOption::INTEGER, "INT"));
    setDefaultValue(*this, "pairedEndSearchSize", "500");
    addOption(*this, ArgParseOption("pi", "abInsParam", "Discordant insertion size: median +/- (MAD * pi)", ArgParseOption::DOUBLE, "FLOAT"));
    setDefaultValue(*this, "abInsParam", "9.0");
    addOption(*this, ArgParseOption("pd", "depthOutlier", "Depth outlier: {Q3 + (IQR * pd)}", ArgParseOption::DOUBLE, "FLOAT"));
    setDefaultValue(*this, "depthOutlier", "1.0");
    addOption(*this, ArgParseOption("pe", "minPairSupport", "SVs supported by >= pe get a vote.", ArgParseOption::DOUBLE, "FLOAT"));
    setDefaultValue(*this, "minPairSupport", "1");
 
    addSection(*this, "Soft-clipped evidence");
    addOption(*this, ArgParseOption("cs", "minClippedSeqSize", "Minimum size of clipped sequence to be considered.", ArgParseOption::INTEGER, "INT"));
    setDefaultValue(*this, "minClippedSeqSize", "20");
    addOption(*this, ArgParseOption("ce", "clippedSeqErrorRate", "Maximum edit distance: floor{length of clipped sequence * (1 - ce)}.", ArgParseOption::DOUBLE, "FLOAT"));
    setDefaultValue(*this, "clippedSeqErrorRate", "0.1");
    //addOption(*this, ArgParseOption("", "use-assembly", "Use consitency based sequence assembly(deprecated)."));
    //setDefaultValue(*this, "use-assembly", "false");
    
    addSection(*this, "Read-depth evidence");
    addOption(*this, ArgParseOption("rs", "samplingNum", "Number of random sample to estimate the background distribution(Q3, IQR, ..) of read-depth evidence.", ArgParseOption::INTEGER, "INT"));
    setDefaultValue(*this, "samplingNum", "100000");
    addOption(*this, ArgParseOption("rw", "readDepthWindowSize", "Window size to caclulate average read-depth around breakpoints.", ArgParseOption::INTEGER, "INT"));
    setDefaultValue(*this, "readDepthWindowSize", "20");
    addOption(*this, ArgParseOption("", "use-re-for-bs", "Use RE for balanced SVs(eg. inverison)."));
    setDefaultValue(*this, "use-re-for-bs", "false");
    addOption(*this, ArgParseOption("re", "reThreshold", "SVs satisfy read-depth evidence >= {Q3 + (IQR * re)} get a vote.", ArgParseOption::DOUBLE, "FLOAT"));
    setDefaultValue(*this, "reThreshold", "1.0");

    // output options
    /* not supported yet
    addSection(*this, "Output");
    addOption(*this, ArgParseOption("o", "output-file", "Specify an output file. Default: write the file to standard output.", ArgParseOption::OUTPUT_FILE));
    setValidValues(*this, "output-file", ".vcf");
    */
}

bool CallOptionManager::parseCommandLine(int argc, char const ** argv)
{
    ArgumentParser::ParseResult res = parse(*this, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return false;

    // alignment files
    getArgumentValue(this->inputFile, *this, 0);
    getOptionValue(this->cutoff, *this, "cutoff");
    getOptionValue(this->minVote, *this, "minVote");
    getOptionValue(this->minMapQual, *this, "minMapQual"); 
    getOptionValue(this->minSVSize, *this, "minSVSize");
    getOptionValue(this->adjTol, *this, "adjTol");
    this->reportFilteredResult = isSet(*this, "report-filtered");
    this->useRankAggregation = !isSet(*this, "no-rank-aggregation");

    getOptionValue(this->minSplitReadSupport, *this, "minSplitReadSupport");
    getOptionValue(this->maxSplit, *this, "maxSplit");
    getOptionValue(this->maxOverlap, *this, "maxOverlap");

    this->doPairedEndRead = !isSet(*this, "no-pe");
    getOptionValue(this->abInsParam, *this, "abInsParam");
    getOptionValue(this->depthOutlier, *this, "depthOutlier");
    getOptionValue(this->minPairSupport, *this, "minPairSupport");
    getOptionValue(this->pairedEndSearchSize, *this, "pairedEndSearchSize");

    this->doClippedRead = !isSet(*this, "no-ce");
    getOptionValue(this->minClippedSeqSize, *this, "minClippedSeqSize");
    getOptionValue(this->referenceGenome, *this, "referenceGenome");
    getOptionValue(this->clippedSeqErrorRate, *this, "clippedSeqErrorRate");
    //this->useAssembly = isSet(*this, "use-assembly");
    this->useAssembly = false;

    this->doReadDepth = !isSet(*this, "no-re");
    getOptionValue(this->readDepthWindowSize, *this, "readDepthWindowSize");
    getOptionValue(this->reThreshold, *this, "reThreshold");
    getOptionValue(this->samplingNum, *this, "samplingNum");
    this->useREforBalancedSV = isSet(*this, "use-re-for-bs");

    return true;
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
    printMessage(std::string(toCString(this->inputFile)));
    printMessage("[General options]");
    printMessage("- referenceGenome: " + std::string(toCString(this->referenceGenome)));
    printMessage("- minMapQual: " + std::to_string(this->minMapQual));
    printMessage("- minSVSize: " + std::to_string(this->minSVSize));
    printMessage("- adjTol: " + std::to_string(this->adjTol));
    printMessage("- cutoff: " + std::to_string(this->cutoff));
    printMessage("- minVote: " + std::to_string(this->minVote));
    printMessage("- reportFilteredResult: " + getBooleanString(this->reportFilteredResult));
    printMessage("- no-pe: " + getBooleanString(!this->doPairedEndRead));
    printMessage("- no-ce: " + getBooleanString(!this->doClippedRead));
    printMessage("- no-re: " + getBooleanString(!this->doReadDepth));
    printMessage("- no-rank-aggregation: " + getBooleanString(!this->useRankAggregation));
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