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
#ifndef APP_SV_H_
#define APP_SV_H_

#include <string>
#include <seqan/vcf_io.h>
#include "breakpoint.hpp"

class VcfRecordEnhanced : public VcfRecord
{
    public :
        TPosition endPos;
        TPosition targetPos; // for duplications and translocations
        char additionalInfo;
        Breakpoint* breakpoint;

        // evidences        
        unsigned se;
        unsigned pe;
        unsigned ce;
        double re;
        double sc;
        double rd;
        double gc;
        double cp;
        char vt;

        char status = STATUS::PASS;
        bool imprecise = false; // do not have exact positions

        enum STATUS {PASS=0, FILTERED=1, MERGED=2};
        static std::string STATUS_PASS(void) { return "PASS"; }
        static std::string STATUS_FILTERED(void) { return "FILTERED"; }
        static std::string STATUS_MERGED(void) { return "MERGED"; }
};  

struct less_than_vcf
{
    inline bool operator() (const VcfRecordEnhanced& v1, const VcfRecordEnhanced& v2)
    {
        if (v1.rID != v2.rID)
            return v1.rID < v2.rID;
        if (v1.beginPos != v2.beginPos)
            return v1.beginPos < v2.beginPos;
        if (v1.endPos != v2.endPos)
            return v1.endPos < v2.endPos;
        return v1.id < v2.id;
    }
};

class SVManager
{
    private :
        BreakpointManager* bpManager;
        OptionManager* opManager;
        AlignmentManager*   alnManager;

        std::map<std::string, std::vector<VcfRecordEnhanced> > sv;

        bool findDeletion(void);
        bool findInversion(void);
        bool findBreakend(void);
        bool findDuplication(void);
        bool findTranslocation(void);
        bool addTranslocation(VcfRecordEnhanced&);
        bool findTranslocation_old(void);
        uint32_t getSVCount(std::string, bool);

    public :
        SVManager(BreakpointManager& bpm) 
        { 
            this->bpManager = &bpm;
            this->alnManager = bpm.getAlignmentManager();
            this->opManager = bpm.getOptionManager();
        }

        bool findSV(void);
        bool filterSV(void);
        bool rescueSV(void);
        bool orderSV(void);
        bool orderSVByEvidenceSum(void);
        bool orderSVByRankAgg(void);

        bool writeVCF(void);
        uint32_t getDeletionCount(void) { return getSVCount(SVTYPE_DELETION(), false); }
        uint32_t getInversionCount(void) { return getSVCount(SVTYPE_INVERSION(), false);  }
        uint32_t getDuplicationCount(void) { return getSVCount(SVTYPE_DUPLICATION(), false);  }
        uint32_t getTranslocationCount(void) { return getSVCount(SVTYPE_TRANSLOCATION(), false);  }
        uint32_t getBreakendCount(void) { return getSVCount(SVTYPE_BREAKEND(), false) / 2; } // always pairs.
        bool filterImpreciseDel(void);

        static std::string SVTYPE_DELETION(void) { return "DEL"; }
        static std::string SVTYPE_INVERSION(void) { return "INV"; }
        static std::string SVTYPE_DUPLICATION(void) { return "DUP"; }
        static std::string SVTYPE_TRANSLOCATION(void) { return "TRA"; }
        static std::string SVTYPE_BREAKEND(void) { return "BND"; }
};

#endif // APP_SV_H_