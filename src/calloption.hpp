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
#ifndef APP_CALLOPTION_H_
#define APP_CALLOPTION_H_


#include <seqan3/argument_parser/all.hpp>

// ==========================================================================
// Tags, Classes, Enums
// ==========================================================================

class CallOptionManager : seqan3::argument_parser
{
    private :
        std::string inputFile;
        std::string outputFile;
        std::string outputFormat;

        // general options
        int32_t minMapQual{20};
        int32_t minSVSize{50};
        int32_t adjTol{50};
        int32_t cutoff{4};
        int32_t minVote{-1};
        int32_t priMethod;
        int32_t filterMethod;
        unsigned int threadCount;
        bool writeBreakpoint{false};
        bool reportFilteredResult{false};
        bool skipRankAggregation{true};

        // for split-read
        double minSplitReadSupport{1};
        int32_t maxSplit{2};
        int32_t maxOverlap{20};

        // for paired-end
        double minPairSupport{1.0};
        double abInsParam{9.0};
        double depthOutlier{1.0};
        int32_t pairedEndSearchSize{500};
        bool skipPairedEndRead{true};

        // for clipped read
        double minClippedReadSupport;
        int32_t minClippedSeqSize{20};
        std::string referenceGenome;
        double clippedSeqErrorRate{0.1};
        bool skipClippedRead{true};
        bool useAssembly;

        // for read-depth analysis
        double reThreshold{1.0};
        int32_t samplingNum{100000};
        int32_t readDepthWindowSize{20};
        bool skipReadDepth{true};
        bool useREforBalancedSV{false};

    public :
        using seqan3::argument_parser::argument_parser;
//     	CallOptionManager()
//     	{
//     		threadCount = 0;
// #ifdef _OPENMP
// 	    	threadCount = std::thread::hardware_concurrency();
// #endif
//     	}

    	void init_options();
      	int parseCommandLine();
        std::string getInputFile(void) { return inputFile; }
        void printUserInput(void);
        std::string getBooleanString(bool);

        // general options
        int32_t getMinMapQual(void) { return minMapQual; }
        int32_t getMinSVSize(void) { return minSVSize; }
        int32_t getAdjTol(void) { return adjTol; }
        int32_t getCutoff(void) { return cutoff; }
        int32_t getMinVote(void) { return minVote; }
        int32_t getPriMethod(void) { return priMethod; }
        int32_t getFilterMethod(void) { return filterMethod; }
        bool getWriteBreakpoint(void) { return writeBreakpoint; }
        bool getReportFilteredResult(void) { return reportFilteredResult; }
        bool getUseRankAggregation(void) { return !skipRankAggregation; }
        void setMinVote(int v) { this->minVote = v;}

        // for split-read
        double getMinSplitReadSupport(void) { return minSplitReadSupport; }
        int32_t getMaxSplit(void) { return maxSplit; }
        int32_t getMaxOverlap(void) { return maxOverlap; }

        // for paired-end
        bool doPairedEndAnalysis(void) { return !skipPairedEndRead; }
        double getAbInsParam(void) { return abInsParam; }
        double getMinPairSupport(void) { return minPairSupport; }
        double getPairedEndSearchSize(void) { return pairedEndSearchSize;}
        double getDepthOutlier(void) { return depthOutlier;}

        // for clipped read
        bool doClippedReadAnalysis(void) { return !skipClippedRead; }
        int32_t getMinClippedSeqSize(void) { return minClippedSeqSize; }
        std::string getReferenceGenome(void) { return referenceGenome; }
        double getClippedSeqErrorRate(void) { return clippedSeqErrorRate; }
        int32_t isUsingAssembly(void) { return useAssembly; }

        // for read-depth analysis
        bool doReadDepthAnalysis(void) { return !skipReadDepth; }
        int32_t getReadDepthWindowSize(void) { return readDepthWindowSize; }
        double getReOutlierCutoff(void) { return reThreshold; }
        int32_t getSamplingNum(void) { return samplingNum; }
        bool getUseREforBalancedSV(void) { return useREforBalancedSV; }
};
#endif // APP_CALLOPTION_H_
