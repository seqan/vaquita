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
#ifndef APP_SVMERGE_H_
#define APP_SVMERGE_H_

#include "sv.hpp"
#include "mergeoption.hpp"


class VcfRecordMultiSample : public VcfRecordEnhanced
{
    public :
        TPosition   inBeginPos, outBeginPos;
        TPosition   inEndPos, outEndPos;
        TPosition   inTargetPos, outTargetPos;
        double      minSC, maxSC;

        std::vector<unsigned> recordNum;
        std::vector<std::string> gt;
        std::vector<unsigned> dp;
};  


class SVMerge
{
    private :
        std::map<std::string, std::vector<VcfRecordMultiSample> >  mergedSV;
        std::vector<SVManager>  svSet;
        MergeOptionManager*     optionManager;
        std::vector<std::string>   fileNames;
        std::vector<std::string>   refNames;

    public :
        SVMerge(MergeOptionManager& om)
        {
            this->optionManager = &om;
        }

        bool loadVcf(std::string&);
        bool loadVcfs(void);
        bool loadVcfs(std::vector<std::string>&);
        bool merge(void);
        bool isReciOverlap(TPosition, TPosition, TPosition, TPosition, double);
        bool writeVCF(void);
};

#endif // APP_SVMERGE_H_