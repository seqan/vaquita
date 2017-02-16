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
#include "vaquita.hpp"
#include "option.hpp"
#include "misc.hpp"

// ============================================================================
// Functions
// ============================================================================

void OptionManager::init()
{
    setAppName(*this, APP_NAME_TEMP);
    setShortDescription(*this, "For Structural Variation");
    setCategory(*this, "For Structural Variation");

    // version & date
    setVersion(*this, SEQAN_APP_VERSION);
    setDate(*this, SEQAN_DATE);

    // synopsis
    addUsageLine(*this, "[\\fIOPTIONS\\fP] <\\fIALIGNMENT(.bam)\\fP>");

    // description
    addDescription(*this, "Vaquita: Fast and Accurate Identification of Structural Variation using Combined Evidence (Kim, 2017)");
    addDescription(*this, "See \\fIhttps://github.com/xenigmax/vaquita\\fP");
    addDescription(*this, "Copyright 2016 by Jongkyu Kim.");

    // mandatory arguments
    addArgument(*this, ArgParseArgument(ArgParseArgument::INPUT_FILE, "ALIGNMENT(.bam)")); 
    setValidValues(*this, 0, ".bam");
    setHelpText(*this, 0, "Support : .bam");

    // Options
    addSection(*this, "General");
    addOption(*this, ArgParseOption("q", "minMapQual", "Mapping quaility cutoff.", ArgParseOption::INTEGER, "INT"));
    setDefaultValue(*this, "minMapQual", "20");
    addOption(*this, ArgParseOption("m", "minSVSize", "Structural varation size cutoff.", ArgParseOption::INTEGER, "INT"));
    setDefaultValue(*this, "minSVSize", "50");
    addOption(*this, ArgParseOption("a", "adjTol", "Positional adjacency in nucleotide resolution.", ArgParseOption::INTEGER, "INT"));
    setDefaultValue(*this, "adjTol", "10");
    addOption(*this, ArgParseOption("v", "minVote", "Number of evidence types that supports SV.", ArgParseOption::INTEGER, "INT"));
    setDefaultValue(*this, "minVote", "2");
    addOption(*this, ArgParseOption("vb", "voteBound", "Parameter for categorization of SVs by depth.", ArgParseOption::DOUBLE, "FLOAT"));
    setDefaultValue(*this, "voteBound", "0.9");
    addOption(*this, ArgParseOption("", "write-breakpoint", "Write breakpoint informtation in a tab-separated format (breakpoints.tsv)."));
    setDefaultValue(*this, "write-breakpoint", "false");
    addOption(*this, ArgParseOption("", "use-global-th", "Use global threshold. (temporary)"));
    setDefaultValue(*this, "use-global-th", "false");

    addSection(*this, "Split-read evidence");
    addOption(*this, ArgParseOption("sr", "minSplitReadSupport", "Minimum number of supporting reads.", ArgParseOption::DOUBLE, "FLOAT"));
    setDefaultValue(*this, "minSplitReadSupport", "2");

    addSection(*this, "Read-pair evidence");
    addOption(*this, ArgParseOption("pr", "minPairSupport", "Minimum number of supporting pairs.", ArgParseOption::DOUBLE, "FLOAT"));
    setDefaultValue(*this, "minPairSupport", "2");
    addOption(*this, ArgParseOption("ps", "pairedEndSearchSize", "Size of the candidate regions.", ArgParseOption::INTEGER, "INT"));
    setDefaultValue(*this, "pairedEndSearchSize", "500");
    addOption(*this, ArgParseOption("pi", "abInsParam", "Abnormal insertion size : median +/- (standard_deviation * abInsParam).", ArgParseOption::DOUBLE, "FLOAT"));
    setDefaultValue(*this, "abInsParam", "5.0");
    addOption(*this, ArgParseOption("", "no-pe", "Do not use information from read-pair evidence."));
    setDefaultValue(*this, "no-pe", "false");
    
    addSection(*this, "Soft-clipped evidence");
    addOption(*this, ArgParseOption("cs", "minClippedSeqSize", "Minimum size of clipped sequence to be considered.", ArgParseOption::INTEGER, "INT"));
    setDefaultValue(*this, "minClippedSeqSize", "20");
    addOption(*this, ArgParseOption("ce", "clippedSeqErrorRate", "Edit distance will be decided by floor{length of clipped sequence * (1 - editDistanceParam)}.", ArgParseOption::DOUBLE, "FLOAT"));
    setDefaultValue(*this, "clippedSeqErrorRate", "0.1");
    addOption(*this, ArgParseOption("cg", "referenceGenome", "Genome sequence file(.fa).", ArgParseOption::INPUT_FILE));
    setValidValues(*this, "referenceGenome", ".fa");
    addOption(*this, ArgParseOption("", "use-assembly", "Use consitency based sequence assembly(deprecated)."));
    setDefaultValue(*this, "use-assembly", "false");
    addOption(*this, ArgParseOption("", "no-ce", "Do not use soft-clipped evidence."));
    setDefaultValue(*this, "no-ce", "false");
    
    addSection(*this, "Read-depth evidence");
    addOption(*this, ArgParseOption("rw", "readDepthWindowSize", "Window size to caclulate average read-depth around breakpoints.", ArgParseOption::INTEGER, "INT"));
    setDefaultValue(*this, "readDepthWindowSize", "100");
    addOption(*this, ArgParseOption("rh", "ddsHigh", "Minimum depth discrepancy score for potential deletions : (rh * depth around breakpoints).", ArgParseOption::DOUBLE, "FLOAT"));
    setDefaultValue(*this, "ddsHigh", "0.05");
    addOption(*this, ArgParseOption("rm", "ddsMid", "Minimum depth discrepancy score for potential duplications and translocations : (rm * depth around breakpoints).", ArgParseOption::DOUBLE, "FLOAT"));
    setDefaultValue(*this, "ddsMid", "0.001");
    addOption(*this, ArgParseOption("rl", "ddsLow", "Minimum depth discrepancy score for insertion : (rl * depth around breakpoints).", ArgParseOption::DOUBLE, "FLOAT"));
    setDefaultValue(*this, "ddsLow", "0.001");
    addOption(*this, ArgParseOption("", "no-re", "Do not use read-depth evidence."));
    setDefaultValue(*this, "no-re", "false");

    // output options
    addSection(*this, "Output");
    addOption(*this, ArgParseOption("o", "output-file", "Specify an output file. Default: write the file to standard output.", ArgParseOption::OUTPUT_FILE));
    setValidValues(*this, "output-file", ".vcf");
}

bool OptionManager::parseCommandLine(int argc, char const ** argv)
{
    ArgumentParser::ParseResult res = parse(*this, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return false;

    // alignment files
    getArgumentValue(this->inputFile, *this, 0);
    getOptionValue(this->minMapQual, *this, "minMapQual"); 
    getOptionValue(this->minSVSize, *this, "minSVSize");
    getOptionValue(this->adjTol, *this, "adjTol");
    getOptionValue(this->minVote, *this, "minVote");
    getOptionValue(this->voteBound, *this, "voteBound");

    getOptionValue(this->minSplitReadSupport, *this, "minSplitReadSupport");

    getOptionValue(this->abInsParam, *this, "abInsParam");
    getOptionValue(this->minPairSupport, *this, "minPairSupport");
    getOptionValue(this->pairedEndSearchSize, *this, "pairedEndSearchSize");

    getOptionValue(this->minClippedSeqSize, *this, "minClippedSeqSize");
    getOptionValue(this->referenceGenome, *this, "referenceGenome");
    getOptionValue(this->clippedSeqErrorRate, *this, "clippedSeqErrorRate");
    this->useAssembly = isSet(*this, "use-assembly");

    getOptionValue(this->ddsHigh, *this, "ddsHigh");
    getOptionValue(this->ddsMid, *this, "ddsMid");
    getOptionValue(this->ddsLow, *this, "ddsLow");
    getOptionValue(this->readDepthWindowSize, *this, "readDepthWindowSize");

    this->writeBreakpoint = isSet(*this, "write-breakpoint");
    this->useGlobalTh = isSet(*this, "use-global-th");
    this->doPairedEndRead = !isSet(*this, "no-pe");
    this->doClippedRead = !isSet(*this, "no-ce");
    this->doReadDepth = !isSet(*this, "no-re");

    return true;
}

void OptionManager::printUserInput(void)
{
    printMessage("==============================");
    printMessage("[General options]");
    printMessage("- inputFile: " + std::string(toCString(this->inputFile)));
    printMessage("- minMapQual: " + std::to_string(this->minMapQual));
    printMessage("- minSVSize: " + std::to_string(this->minSVSize));
    printMessage("- adjTol: " + std::to_string(this->adjTol));
    printMessage("- minVote: " + std::to_string(this->minVote));
    printMessage("- voteBound: " + std::to_string(this->voteBound));
    printMessage("- no-pe: " + std::to_string(this->doPairedEndRead));
    printMessage("- no-ce: " + std::to_string(this->doClippedRead));
    printMessage("- no-re: " + std::to_string(this->doReadDepth));
    printMessage("[Split-read option]");
    printMessage("- minSplitReadSupport: " + std::to_string(this->minSplitReadSupport));
    printMessage("- [Read-pair options]");
    printMessage("- abInsParam: " + std::to_string(this->abInsParam));
    printMessage("- minPairSupport: " + std::to_string(this->minPairSupport));
    printMessage("- pairedEndSearchSize: " + std::to_string(this->pairedEndSearchSize));
    printMessage("[Soft-clipped-read options]");
    printMessage("- minClippedSeqSize: " + std::to_string(this->minClippedSeqSize));
    printMessage("- referenceGenome: " + std::string(toCString(this->referenceGenome)));
    printMessage("- clippedSeqErrorRate: " + std::to_string(this->clippedSeqErrorRate));
    printMessage("- use-assembly: " + std::to_string(this->useAssembly));
    printMessage("[Read-depth options]");
    printMessage("- ddsHigh: " + std::to_string(this->ddsHigh));
    printMessage("- ddsMid: " + std::to_string(this->ddsMid));
    printMessage("- ddsLow: " + std::to_string(this->ddsLow));
    printMessage("- readDepthWindowSize: " + std::to_string(this->readDepthWindowSize));
    printMessage("==============================");
}       
