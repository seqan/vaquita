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
#ifndef APP_ALIGNMENT_H_
#define APP_ALIGNMENT_H_

#include "option.hpp"
#include "candidate.hpp"

// ==========================================================================
// Tags, Classes, Enums
// ==========================================================================
typedef std::pair<int32_t, BamAlignmentRecord> TInsBamRecordPair;
class AlignmentManager
{
    private :
        static unsigned const INS_SIZE_ESTIMATION_SAMPLE_SIZE = 10000;
        double K = 1.4826;

        // option
        OptionManager* optionManager;

        // bam record
        BamFileIn       bamFileIn;
        BamHeader       bamHeader;
        BamIndex<Bai>   baiIndex;
        BamFileOut*     pBamFileOut;

        // breakpoint candidates
        BreakpointCandidate* splitRead;
        BreakpointCandidate* pairedEndRead;
        BreakpointCandidate* clippedRead;
        BreakpointCandidate* readDepth;

        // information of a dataset
        std::map<TTemplateID, TPosition> templateLength;
        std::vector<TInsBamRecordPair> readsForMedInsEst;
        std::vector<BamAlignmentRecord> readsBuffer;
        double insertMedian;
        double insertDev; // median absolute deviation
        double minAbInsSize;
        double maxAbInsSize;
        int32_t totalRecordNum;
        int32_t splitReadCount;
        int32_t pairedReadCount;
        int32_t clippedReadCount;

        void calcInsertSize();
        //bool isAbnormalInsertion(unsigned ins) { return ( (ins < minAbInsSize) || (ins > maxAbInsSize)); }
        bool isAbnormalInsertion(unsigned ins) { return (ins > maxAbInsSize); }
        static bool pairCompare(const TInsBamRecordPair& e1, const TInsBamRecordPair& e2) { return e1.first < e2.first; }

    public :
        AlignmentManager() : optionManager(NULL) {}
    	AlignmentManager(OptionManager & op) { init(op); }
        ~AlignmentManager() 
        {
            // TODO : release memory! (TBamRecords)
            close(bamFileIn);
            close(*pBamFileOut);
        }
   	
        void init(OptionManager & op)
        { 
            optionManager = &op; 
            insertMedian = 0; 
            insertDev = 0.0; 
            minAbInsSize = std::numeric_limits<double>::min();
            maxAbInsSize = std::numeric_limits<double>::max();
            splitReadCount = 0;
            pairedReadCount = 0;
            clippedReadCount = 0;
        }

        void setBreakpointCandidate(BreakpointCandidate* s, BreakpointCandidate* p, BreakpointCandidate* e, BreakpointCandidate* r)
        {
            splitRead = s;
            pairedEndRead = p;
            clippedRead = e;
            readDepth = r;
        }

        // [start, end)
        void getSequenceAndDepth(CharString&, std::vector<int32_t>&, TTemplateID, TPosition, TPosition);
        void getSequence(CharString&, TTemplateID, TPosition, TPosition);
        void getDepth(std::vector<int32_t>&, TTemplateID, TPosition, TPosition);

        int32_t getSplitReadCount(void) { return splitReadCount; }
        int32_t getPairedReadCount(void) { return pairedReadCount; }
        int32_t getClippedReadCount(void) { return clippedReadCount; }
        int32_t getTotalRecordNum() { return totalRecordNum; }
        double getInsMedian() { return insertMedian; }
        double getInsSD() { return insertDev; }
        double getAbInsParam() { return optionManager->getAbInsParam(); }
        CharString getRefName(int32_t id);
        int32_t getRefCount();
        std::map<TTemplateID, TPosition>* getTemplateLength(void) { return &templateLength; }

        OptionManager* getOptionManager(void) { return this->optionManager; }
        void printRecord(BamAlignmentRecord &);
        bool load(void);
};
#endif // APP_ALIGNMENT_H_