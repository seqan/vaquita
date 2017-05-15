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
#ifndef APP_CLIPPEDREAD_H_
#define APP_CLIPPEDREAD_H_

#include <seqan/score.h>  // The module score.
#include <seqan/index.h>
#include "candidate.hpp"
#include "alignment.hpp"

struct ClippedSequenceSegment
{
    TTemplateID matchedRightTemplateID;
    TTemplateID matchedLeftTemplateID;
    Breakpoint* matchedBp;
    SequenceSegment sequenceSegment;
    SequenceSegment querySegment;
    bool isReverseComplemented;
    BreakpointEvidence::ORIENTATION  orientation;
};

struct ClippedSequenceSegmentCompareByPosition {
    bool operator()(const ClippedSequenceSegment &lhs, const ClippedSequenceSegment &rhs) {
        if (lhs.sequenceSegment.templateID == rhs.sequenceSegment.templateID)
        {
            if (lhs.sequenceSegment.beginPos != rhs.sequenceSegment.beginPos)
                return lhs.sequenceSegment.beginPos < rhs.sequenceSegment.beginPos;
            else
                return lhs.sequenceSegment.endPos < rhs.sequenceSegment.endPos;
        }
        else
            return false;
    }
};

struct SearchRegion
{
    SequenceSegment segment;
    Dna5String sequence;
};

struct CandidateRegionPair
{
    SearchRegion leftRegion;
    SearchRegion rightRegion;
};

typedef std::set<ClippedSequenceSegment, ClippedSequenceSegmentCompareByPosition> TFoundPositionSet;
typedef std::vector<ClippedSequenceSegment>                                       TFoundPosition;
class ClippedRead : public BreakpointCandidate
{
    private :
        static int32_t const MAX_SCORE = MaxValue<int32_t>::VALUE;

        FaiIndex faiIndex;
        std::map<Breakpoint*, CandidateRegionPair> candidateRegion;
        std::map<TTemplateID, TPosition> templateSize;
       
    public:
        void prepAfterHeaderParsing(BamHeader&, BamFileIn&);
        void setOptionManager(OptionManager* op);
        void parseReadRecord(TReadName&, BamAlignmentRecord&);
        void setSearchRegionByOrientation(const BreakpointEvidence::ORIENTATION, const BreakpointCandidate::SIDE, Breakpoint&, TTemplateID&, TPosition&, TPosition&);

        void getReferenceSequence(CharString& seq, CharString chr, TPosition start, TPosition end);
        void getReferenceSequence(CharString& seq, TTemplateID templateID, TPosition start, TPosition end);

        bool searchPairRegion(TFoundPosition&, Breakpoint*, int32_t&, CharString&, SIDE, bool, bool, BreakpointEvidence::ORIENTATION );
        bool searchTwilightZone(TFoundPosition&, Breakpoint*, int32_t&, CharString&, SIDE, bool, BreakpointEvidence::ORIENTATION );
        bool alignByMyersBitVector(TFoundPosition&, CharString&, CharString&, int32_t&);
        bool alignByLocal(TFoundPosition&, CharString&, CharString&, int32_t&);
        bool onlineSearchBySegment(TFoundPosition&, CharString&, CharString&, int32_t&);

        void getConsensusSequence(CharString&, Breakpoint*);
        static int32_t getMaxScore(void) { return ClippedRead::MAX_SCORE; }
};

#endif // APP_CLIPPEDREAD_H_