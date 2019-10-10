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
#ifndef APP_BREAKPOINTCANDIDATE_H_
#define APP_BREAKPOINTCANDIDATE_H_

#define GENOMIC_BIN(x,y) x / y

#include <set>
#include <map>
#include <vector>
#include <utility>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include "calloption.hpp"
#include "intervalindex.hpp"

// ==========================================================================
// Tags, Classes, Enums
// ==========================================================================
typedef uint8_t                         TTemplateID;
typedef uint32_t                        TPosition;
typedef uint16_t                        TRelativePosition;
typedef seqan::CharString               TReadName;
typedef uint32_t                        TReadID;
typedef bool                            TStrand;
typedef std::vector<TPosition>          TPositionList;
typedef std::vector<TRelativePosition>  TRelativePositionList;

struct SequenceSegment
{
    TTemplateID templateID; // eg. chromosome, readID
    TPosition beginPos, endPos;
    TStrand isReverse = false;
};
typedef std::set<SequenceSegment>                                   TSequenceSegmentVector;

struct BreakpointEvidence
{
    static TPosition const INVALID_POS = seqan::MaxValue<TPosition>::VALUE;
    static TTemplateID const NOVEL_TEMPLATE = seqan::MaxValue<TTemplateID>::VALUE;

    // P.: -><-, S.: <-->, I.: ->->, LEFT/RIGHT CLIP : temporary use.
    enum ORIENTATION {PROPERLY_ORIENTED_LARGE=0, PROPERLY_ORIENTED_SMALL=1, SWAPPED=2, INVERTED=3, NOT_DECIDED=4, CLIPPED=5, NUM_OF_ORIENTATION=6};
    enum SIDE {LEFT, RIGHT, BOTH};

    SequenceSegment leftSegment, rightSegment;
    TReadID suppRead;
    seqan::DnaString sequence; // eg. clipped sequence
    ORIENTATION orientation;
};

struct AlignmentInfo
{
    SequenceSegment refSegment, querySegment;
    TPosition querySize;
    std::vector<BreakpointEvidence> indelList, clippedList;
};

typedef std::vector<std::pair<TPosition, seqan::DnaString> >              TClippedSequences;
class Breakpoint
{
    public :
        TTemplateID leftTemplateID, rightTemplateID;
        TStrand leftReverseFlag, rightReverseFlag;
        BreakpointEvidence::ORIENTATION orientation;

        // candidate positions
        TPosition       minLeftPos, maxLeftPos, minRightPos, maxRightPos;
        TPositionList   leftPos, rightPos;

        // supporting reads
        uint16_t suppReads;

        // clipped sequences
        TClippedSequences clippedSequences;
        TPosition clippedConsensusSequenceSize = 0;

        bool needLeftIndexUpdate = false;
        bool needRightIndexUpdate = false;
        bool bFoundExactPosition = false;
};

typedef seqan::String<Breakpoint*>                                  TBreakpointList;
typedef std::set<Breakpoint*>                                       TBreakpointSet;
typedef std::vector<Breakpoint*>                                    TBreakpointVector;
typedef std::pair<std::pair<TPosition, TPosition>, Breakpoint*>     TPosBreakpoint;
typedef std::vector<TPosBreakpoint>                                 TPosBreakpointVector;
typedef std::map<TPosition, TPosBreakpointVector*>                  TBreakpointBinIndex;
typedef std::map<TTemplateID, TBreakpointBinIndex>                  TBreakpointIndex;

typedef IntervalIndex<Breakpoint*>  TBreakpointIntervalIndex;
typedef std::map<TTemplateID, TBreakpointIntervalIndex*>     TBreakpointIntervalIndexMap;

class BreakpointCandidate
{
    private :
        static TPosition const GENOMIC_BIN_SIZE = 1000;
        TReadID currentReadID = 0;

        int32_t posAdj;
        double insertMedian;
        double insertDev;
        double maxAbInsSize;
        double minAbInsSize;

        CallOptionManager* op;
        TBreakpointIntervalIndexMap leftIndexMap, rightIndexMap;
        TBreakpointSet     breakpoints;

        void setPositionWithAdj(TPosition &, TPosition &);

        TBreakpointIntervalIndex* getIndex(BreakpointEvidence::SIDE, TTemplateID);
        bool addIndex(TBreakpointIntervalIndex*, TPosition, TPosition, Breakpoint*);
        bool removeIndex(TBreakpointIntervalIndex*, TPosition, TPosition, Breakpoint*);

    public :
        BreakpointCandidate(CallOptionManager*);
        ~BreakpointCandidate();

        void setPositionalAdj(int32_t a) { this->posAdj = a; }
        int32_t getPositionalAdj(void) { return this->posAdj; }

        // insertion information
        void setInsertionInfo(double, double, double, double);
        double getInsMedian() { return insertMedian; }
        double getInsSD() { return insertDev; }
        double getMaxAbInsSize() { return maxAbInsSize; }
        double getMinAbInsSize() { return minAbInsSize; }

        // options
        void setOptionManager(CallOptionManager* op) { this->op = op; setPositionalAdj(op->getAdjTol()); }
        CallOptionManager* getOptionManager(void) { return this->op; }
        int32_t getBreakpointCount(void) { return this->breakpoints.size(); }

        // operations for breakpoints
        TReadID getCurrentReadID(void) { return this->currentReadID; }
        TReadID getNextReadID(void) { return this->currentReadID++; }

        std::set<Breakpoint*>* getCandidateSet() { return &this->breakpoints; }     // GET candidate set

        void addNewBreakpoint(Breakpoint*); // addition of a new breakpoint
        void findBreakpoint(TBreakpointSet&, TBreakpointSet&, TBreakpointSet&, Breakpoint*);

        Breakpoint* copyAndUpdateBreakpoint(Breakpoint*, bool&);                    // ADD breakpoints by copying
        Breakpoint* moveAndUpdateBreakpoint(Breakpoint*, bool&);                    // ADD breakpoints by copying
        TBreakpointSet::iterator mergeBreakpoint(Breakpoint*, Breakpoint*);         // MERGE two breakpoints
        TBreakpointSet::iterator removeBreakpoint(TBreakpointSet::iterator);        // REMOVE breakpoints
        Breakpoint* updateBreakpoint(Breakpoint*, bool&);  // UPDATE breakpoints including adjustment of indicies
        Breakpoint* updateBreakpoint(BreakpointEvidence&, bool&);
        void updateBreakpointIndex(Breakpoint*);                                    // UPDATE index only
        TBreakpointSet::iterator removeBreakpoint(Breakpoint*);

        bool isOverlap(SequenceSegment&, SequenceSegment&);                         // checking for overlap of 2 intervals
        bool isAdjacent(SequenceSegment&, SequenceSegment&);                        // checking for adjancy of 2 intervals

        // defined by derived functions
        virtual void prepAfterHeaderParsing(seqan::BamHeader& header, seqan::BamFileIn& fileIn) { return; }
        virtual void parseReadRecord(TReadName&, seqan::BamAlignmentRecord&) { return; }
        virtual void checkReadRecord(TReadName&, seqan::BamAlignmentRecord&) { return; }
        virtual bool analyze(void) { return true; }
        virtual void doAdditionalJobAfterMerge(Breakpoint*, Breakpoint*) { return; }
        virtual bool isNew(TReadName&) { return true; }

        // static functions
        static bool isOverlap(TPosition, TPosition, TPosition, TPosition);
        static bool isAdjacent(TPosition, TPosition, TPosition, TPosition, TPosition);
        static bool isAdjacent(TPosition, TPosition, TPosition);
        static void setPositionWithAdj(TPosition &, TPosition &, TPosition);        // make extended intervals using positional adjancy
        static void printBreakpoint(Breakpoint*);                                   // print breakpoints to stderr
        static void printBreakpoint(Breakpoint&);
        static void copyBreakpoint(Breakpoint&, Breakpoint&);                       // copy breakpoint
        static void moveBreakpoint(Breakpoint&, Breakpoint&);

        static void parseCIGAR(AlignmentInfo&, TReadID&, seqan::BamAlignmentRecord&, bool, bool, TPosition, TPosition); // get Alignment info from CIGAR
        static bool compareByQueryPos(AlignmentInfo&, AlignmentInfo&);              // for ordering with sorting function
        static bool compareByChrmAndPos(Breakpoint&, Breakpoint&);                  // for ordering with sorting function
        static bool isMatchedBreakpoint(Breakpoint*, Breakpoint*);
        static void clearPosInfo(Breakpoint*);
        static double PREVENT_DIV_BY_ZERO(void) { return 0.0000000001; }

        static void updateMinMaxPos(Breakpoint*);
        static void updateLeftMinMaxPos(Breakpoint*);
        static void updateRightMinMaxPos(Breakpoint*);
};

#endif // APP_BREAKPOINTCANDIDATE_H_
