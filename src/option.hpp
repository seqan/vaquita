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
#ifndef APP_OPTION_H_
#define APP_OPTION_H_

#include <seqan/sequence.h>
#include <seqan/arg_parse.h>

using namespace seqan;

// ==========================================================================
// Tags, Classes, Enums
// ==========================================================================
class OptionManager : public ArgumentParser
{
    private :
        CharString inputFile;
        CharString outputFile;
        CharString outputFormat;

        // general options
        int32_t minMapQual;
        int32_t minSVSize;
        int32_t adjTol;
        int32_t minVote;
        double voteBound;
        unsigned int threadCount;
        bool writeBreakpoint;
        bool useGlobalTh;

        // for split-read
        double minSplitReadSupport;

        // for paired-end
        double minPairSupport;
        double abInsParam;
        int32_t pairedEndSearchSize;
        bool doPairedEndRead;

        // for clipped read
        double minClippedReadSupport;
        int32_t minClippedSeqSize;
        CharString referenceGenome;
        double clippedSeqErrorRate;
        bool doClippedRead;
        bool useAssembly;

        // for read-depth analysis
        double ddsHigh;
        double ddsMid;
        double ddsLow;
        double avgReadDepth;
        int32_t readDepthWindowSize;
        bool doReadDepth;

    public :
    	OptionManager()
    	{
    		threadCount = 0;
#ifdef _OPENMP
	    	threadCount = std::thread::hardware_concurrency();
#endif
    	}

    	void init(void);
      	bool parseCommandLine(int argc, char const ** argv);
        CharString getInputFile(void) { return inputFile; }
        
        // general options
        int32_t getMinMapQual(void) { return minMapQual; } 
        int32_t getMinSVSize(void) { return minSVSize; }
        int32_t getAdjTol(void) { return adjTol; }
        int32_t getMinVote(void) { return minVote; }
        double getVoteBound(void) { return voteBound; }
        bool getWriteBreakpoint(void) { return writeBreakpoint; }
        bool useGlobalThreshold(void) { return useGlobalTh; }

        // for split-read
        double getMinSplitReadSupport(void) { return minSplitReadSupport; }        

        // for paired-end
        double getAbInsParam(void) { return abInsParam; }
        double getMinPairSupport(void) { return minPairSupport; }
        double getPairedEndSearchSize(void) { return pairedEndSearchSize;}
        bool doPairedEndAnalysis(void) { return doPairedEndRead; }

        // for clipped read
        int32_t getMinClippedSeqSize(void) { return minClippedSeqSize; }
        CharString getReferenceGenome(void) { return referenceGenome; }
        double getClippedSeqErrorRate(void) { return clippedSeqErrorRate; }
        int32_t isUsingAssembly(void) { return useAssembly; }
        bool doClippedReadAnalysis(void) { return doClippedRead; }

        // for read-depth analysis
        void setAverageReadDepth(double x) { avgReadDepth = x; }
        double getAverageReadDepth(void) { return avgReadDepth; }
        bool doReadDepthAnalysis(void) { return doReadDepth; }
        int32_t getReadDepthWindowSize(void) { return readDepthWindowSize; }
        double getDDSHigh(void) { return ddsHigh; }
        double getDDSMid(void) { return ddsMid; }
        double getDDSLow(void) { return ddsLow; }
        //double getReadDepthDiffTh(double x) { return (avgReadDepth*sigMax) / (1.0 + exp(-sigSteepness * (x-(avgReadDepth * (1+sigMid))))); }
        //double getReadDepthGlobalTh(void) { return sigMax; }
};
#endif // APP_OPTION_H_
