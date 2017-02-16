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
// ===========================================/===============================
#include "vaquita.hpp"
#include "sv.hpp"
#include "misc.hpp"

bool SVManager::findSV(void)
{
    findDeletion();
    findInversion();
    findDuplication(); // required : deletion
    findTranslocation(); // required : duplication
    findBreakend(); // the others
    return true;
}

bool SVManager::writeVCF(void)
{
    // merge to a single list
    std::vector<VcfRecordEnhanced> vcfRecords;
    for (auto itSVType = this->sv.begin(); itSVType != this->sv.end(); ++itSVType)
    {    
        std::sort(itSVType->second.begin(), itSVType->second.end(), less_than_vcf());
        int32_t nID = 1;
        for (auto itSV = itSVType->second.begin(); itSV != itSVType->second.end(); ++itSV)
        {
            // skip pseudo deletions
            if (itSVType->first == SVTYPE_DELETION() && itSV->isPseudoDeletion == true)
                continue;

            /* for 6 lines of representation of translocations
            if (itSVType->first != SVTYPE_TRANSLOCATION())
                itSV->id = itSVType->first + "_" + std::to_string(nID++);
            */

            vcfRecords.push_back(*itSV);
        }
    }
    // sort by chromosome & position
    std::sort(vcfRecords.begin(), vcfRecords.end(), less_than_vcf());


    // init.
    VcfHeader vcfHeader;
    VcfFileOut vcfOut;
    open(vcfOut, std::cout, Vcf());

    // headers
    appendValue(vcfHeader, VcfHeaderRecord("fileformat", "VCFv4.1"));
    appendValue(vcfHeader, VcfHeaderRecord("source", APP_NAME));
    appendValue(vcfHeader, seqan::VcfHeaderRecord("INFO", "<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variation\">"));
    appendValue(vcfHeader, seqan::VcfHeaderRecord("INFO", "<ID=SVLEN,Number=1,Type=Integer,Description=\"Size of structural variation compared to reference\">"));
    appendValue(vcfHeader, seqan::VcfHeaderRecord("INFO", "<ID=TARGETPOS,Number=1,Type=Integer,Description=\"Position of the newly inserted sequence in duplication or translocations\">"));
    appendValue(vcfHeader, seqan::VcfHeaderRecord("INFO", "<ID=RT,Number=1,Type=Float,Description=\"Number of split-read and read-pairs supporting the structural variation\">"));
    appendValue(vcfHeader, seqan::VcfHeaderRecord("INFO", "<ID=SE,Number=1,Type=Integer,Description=\"Number of split-reads supporting the structural variation\">"));
    appendValue(vcfHeader, seqan::VcfHeaderRecord("INFO", "<ID=PE,Number=1,Type=Integer,Description=\"Number of read-pairs supporting the structural variation\">"));
    appendValue(vcfHeader, seqan::VcfHeaderRecord("INFO", "<ID=RE,Number=1,Type=Float,Description=\"Number of read-pairs supporting the structural variation\">"));
    appendValue(vcfHeader, seqan::VcfHeaderRecord("INFO", "<ID=RD,Number=1,Type=Float,Description=\"Read-depth around structural variation\">"));
    appendValue(vcfHeader, seqan::VcfHeaderRecord("ALT", "<ID=DEL,Description=\"Deletion\">"));
    appendValue(vcfHeader, seqan::VcfHeaderRecord("ALT", "<ID=INV,Description=\"Inversion\">"));
    appendValue(vcfHeader, seqan::VcfHeaderRecord("ALT", "<ID=DUP,Description=\"Duplication\">"));
    appendValue(vcfHeader, seqan::VcfHeaderRecord("ALT", "<ID=DUP:TANDEM,Description=\"Tandem duplication\">"));
    appendValue(vcfHeader, seqan::VcfHeaderRecord("ALT", "<ID=TRA,Description=\"Translocation\">"));
    appendValue(vcfHeader, seqan::VcfHeaderRecord("ALT", "<ID=BND,Description=\"Breakend\">"));
    appendValue(vcfHeader, seqan::VcfHeaderRecord("FORMAT", "<ID=GT,Number=1,Type=String,Description=\"Genotype\">"));

    // reference & sample
    for (int i=0; i < this->alnManager->getRefCount(); ++i)
        appendValue(contigNames(context(vcfOut)), this->alnManager->getRefName(i));
    appendValue(sampleNames(context(vcfOut)), "SAMPLE");

    // write
    writeHeader(vcfOut, vcfHeader);
    for (auto itSV = vcfRecords.begin(); itSV != vcfRecords.end(); ++itSV)
        writeRecord(vcfOut, *itSV);

    return true;
}

bool SVManager::addTranslocation(int32_t nID, VcfRecordEnhanced& orgRecord)
{
    FinalBreakpointInfo* finalBreakpoint = this->bpManager->getFinalBreakpointInfo(orgRecord.breakpoint);
    ReadSupportInfo* info = this->bpManager->getMergedBreakpoint()->getReadSupport(orgRecord.breakpoint);
    TPosition svLen = orgRecord.endPos - orgRecord.beginPos + 1;
    VcfRecordEnhanced record = orgRecord;

    std::string recordID = SVTYPE_TRANSLOCATION() + "_" + std::to_string(nID);
    record.alt = "<" + SVTYPE_TRANSLOCATION() + ">";
    record.info  = "RT=" + std::to_string(finalBreakpoint->score) + ";";
    record.info += "SE=" + std::to_string(orgRecord.se) + ";";
    record.info += "PE=" + std::to_string(orgRecord.pe) + ";";
    record.info += "CE=" + std::to_string(orgRecord.ce) + ";";
    record.info += "RE=" + std::to_string(orgRecord.re) + ";";
    record.info += "RD=" + std::to_string(info->avgReadDepth)+ ";";
    record.info += "SVTYPE=" + SVTYPE_TRANSLOCATION() + ";";
    record.info += "SVLEN=" + std::to_string(svLen);
            
    // TODOs (is this the best way to represent translocations?)
    record.ref = "N"; // before SV
    record.qual = VcfRecord::MISSING_QUAL();
    record.filter = "PASS";       
    record.format = "GT";
    appendValue(record.genotypeInfos, "./.");

    record.id = recordID;
    record.beginPos = orgRecord.beginPos;
    record.endPos = orgRecord.endPos;
    record.targetPos = orgRecord.targetPos;
    record.info += ";TARGETPOS=" + std::to_string(record.targetPos);

    sv[SVTYPE_TRANSLOCATION()].push_back(record);

    /* 6 lines of representation. Legacy of Trappe 2015
    record.id = recordID + "_1";
    record.beginPos = orgRecord.beginPos - 1;
    sv[SVTYPE_TRANSLOCATION()].push_back(record);
    record.id = recordID + "_2";
    record.beginPos = orgRecord.beginPos;
    sv[SVTYPE_TRANSLOCATION()].push_back(record);

    record.id = recordID + "_3";
    record.beginPos = orgRecord.endPos - 1;
    sv[SVTYPE_TRANSLOCATION()].push_back(record);
    record.id = recordID + "_4";
    record.beginPos = orgRecord.endPos;
    sv[SVTYPE_TRANSLOCATION()].push_back(record);

    record.id = recordID + "_5";
    record.beginPos = orgRecord.targetPos - 1;
    sv[SVTYPE_TRANSLOCATION()].push_back(record);
    record.id = recordID + "_6";
    record.beginPos = orgRecord.targetPos;
    sv[SVTYPE_TRANSLOCATION()].push_back(record);
    */
}

bool SVManager::findTranslocation(void)
{
    int32_t nID = 1;
    int32_t adjTol = this->opManager->getAdjTol();  

    // true : remove it from the duplication list
    std::vector<bool> dupRemoveList;
    for (auto i = 0; i < this->sv[SVTYPE_DUPLICATION()].size(); ++i)
        dupRemoveList.push_back(false);

    // sort duplications according to chromosome & positions
    std::sort(this->sv[SVTYPE_DUPLICATION()].begin(), this->sv[SVTYPE_DUPLICATION()].end(), less_than_vcf());    
    for (auto itDup = this->sv[SVTYPE_DUPLICATION()].begin(); itDup != this->sv[SVTYPE_DUPLICATION()].end(); ++itDup)
    {
        VcfRecordEnhanced record = *itDup;
        bool foundMatch = false;

        // find matched duplication
        for (auto itDupComp = itDup + 1; itDupComp != this->sv[SVTYPE_DUPLICATION()].end(); ++itDupComp)
        {
            // L,L or R,R
            if (itDup->additionalInfo == itDupComp->additionalInfo)
                continue;

            // check matched duplications to register them to translocations
            if ( (itDup->rID == itDupComp->rID) && \
                  BreakpointCandidate::isAdjacent(itDup->beginPos, itDup->beginPos, itDupComp->targetPos, itDupComp->targetPos, adjTol) && \
                  BreakpointCandidate::isAdjacent(itDupComp->endPos, itDupComp->endPos, itDup->targetPos, itDup->targetPos, adjTol))
            {
                foundMatch = true;
                record = *itDup;
                record.se += itDup->se;
                record.pe += itDup->pe;
                record.ce += itDup->ce;
                record.re += itDup->re;
                addTranslocation(nID, record);
                ++nID;

                dupRemoveList[itDup - this->sv[SVTYPE_DUPLICATION()].begin()] = true;
                dupRemoveList[itDupComp - this->sv[SVTYPE_DUPLICATION()].begin()] = true;
            }
        }
    }

    // true duplications
    std::vector<VcfRecordEnhanced> dupAfter;
    for (auto i = 0; i < this->sv[SVTYPE_DUPLICATION()].size(); ++i)
        if (dupRemoveList[i] == false)
            dupAfter.push_back(this->sv[SVTYPE_DUPLICATION()][i]); 
    this->sv[SVTYPE_DUPLICATION()].swap(dupAfter);

    return true;
}

bool SVManager::findDuplication(void)
{
    MergedCandidate* mergedBreakpoint = this->bpManager->getMergedBreakpoint();
    TBreakpointSet* candidateSet = mergedBreakpoint->getCandidateSet();
    auto itBreakpoint = candidateSet->begin();
    int32_t nID = 1;
    int32_t adjTol = this->opManager->getAdjTol();

    while (itBreakpoint != candidateSet->end())
    {
        Breakpoint* bp = *itBreakpoint;
        ReadSupportInfo* info =  mergedBreakpoint->getReadSupport(bp);
        FinalBreakpointInfo* finalBreakpoint = this->bpManager->getFinalBreakpointInfo(bp);

        bool found = true;

        // not in the same template
        if(bp->leftTemplateID == BreakpointEvidence::NOVEL_TEMPLATE || bp->leftTemplateID != bp->rightTemplateID)
            found = false;

        // not an inversion
        if (bp->orientation != BreakpointEvidence::ORIENTATION::SWAPPED)
            found = false;
        
        int32_t svLen = finalBreakpoint->rightPosition - finalBreakpoint->leftPosition + 1;
        if (svLen < this->opManager->getMinSVSize())
            found = false;

        if (found == true)
        {
            // 'L' : left-side, 'R' : right-side, 'I' : imprecise
            std::vector<std::pair<char, VcfRecordEnhanced> > svDel;
            svDel.clear();
            auto itDel = this->sv[SVTYPE_DELETION()].begin();
            while( itDel != this->sv[SVTYPE_DELETION()].end() )
            {
                bool matchFound = false;

                // left-side overlap
                if ( itDel->rID == finalBreakpoint->leftTemplateID && \
                     BreakpointCandidate::isAdjacent(itDel->beginPos, itDel->beginPos, finalBreakpoint->leftPosition, finalBreakpoint->leftPosition, adjTol))
                {
                    if (itDel->endPos < finalBreakpoint->rightPosition)
                    {
                        svDel.push_back(std::make_pair('L', *itDel));
                        matchFound = true;
                    }
                }

                // right-side overlap
                if ( itDel->rID == finalBreakpoint->rightTemplateID && \
                     BreakpointCandidate::isAdjacent(itDel->endPos, itDel->endPos, finalBreakpoint->rightPosition, finalBreakpoint->rightPosition, adjTol))
                {
                    if (finalBreakpoint->leftPosition < itDel->beginPos)
                    {
                        svDel.push_back(std::make_pair('R', *itDel));
                        matchFound = true;
                    }
                }

                if (matchFound)
                {
                    itDel = this->sv[SVTYPE_DELETION()].erase(itDel);
                    // break
                }
                else
                    ++itDel;
            }

            // imprecise duplication
            if (svDel.size() == 0)
            {
                VcfRecordEnhanced tempSvDel;
                tempSvDel.se = 0;
                tempSvDel.pe = 0;
                tempSvDel.ce = 0;
                tempSvDel.re = 999999999.0;
                svDel.push_back(std::make_pair('I', tempSvDel));
            }

            for (auto itSvDel = svDel.begin(); itSvDel != svDel.end(); ++itSvDel)
            {
                char dupType = itSvDel->first;
                auto currentSvDel = itSvDel->second;

                VcfRecordEnhanced record;
                record.breakpoint = bp;
                record.additionalInfo = dupType;
                if (record.additionalInfo == 'L')
                {
                    record.beginPos = currentSvDel.endPos;
                    record.endPos = finalBreakpoint->rightPosition;
                    record.targetPos = finalBreakpoint->leftPosition;
                }
                else if (record.additionalInfo == 'R')
                {
                    record.beginPos = finalBreakpoint->leftPosition;
                    record.endPos = currentSvDel.beginPos;
                    record.targetPos = finalBreakpoint->rightPosition;
                }
                else // 'I'
                {
                    record.beginPos = finalBreakpoint->leftPosition;
                    record.endPos = finalBreakpoint->rightPosition;
                }
                TPosition svLen = record.endPos - record.beginPos + 1;

                record.rID = finalBreakpoint->leftTemplateID;
                record.id = SVTYPE_DUPLICATION() + "_" + std::to_string(nID++);
                record.se = currentSvDel.se + info->splitReadSupport;
                record.pe = currentSvDel.pe + info->pairedEndSupport;
                record.ce = currentSvDel.ce + info->clippedReadSupport;
                if (info->leftReadDepthSelected == true)
                    record.re = info->leftReadDepthDiffScore;
                else
                    record.re = info->rightReadDepthDiffScore;
                    
                record.info  = "RT=" + std::to_string(finalBreakpoint->score) + ";";
                record.info += "SE=" + std::to_string(record.se) + ";";
                record.info += "PE=" + std::to_string(record.pe) + ";";
                record.info += "CE=" + std::to_string(record.ce) + ";";
                record.info += "RE=" + std::to_string(record.re) + ";";
                record.info += "RD=" + std::to_string(info->avgReadDepth)+ ";";
                record.info += "SVTYPE=" + SVTYPE_DUPLICATION() + ";";
                record.info += "SVLEN=" + std::to_string(svLen);

                if (record.additionalInfo == 'I')
                {
                    record.alt = "<" + SVTYPE_DUPLICATION() + ":TANDEM>";
                }
                else
                {
                    record.alt = "<" + SVTYPE_DUPLICATION() + ">";
                    record.info += ";TARGETPOS=" + std::to_string(record.targetPos);
                }

                // TODOs
                record.ref = "N"; // before SV
                record.qual = VcfRecord::MISSING_QUAL();
                record.filter = "PASS";       
                record.format = "GT";
                appendValue(record.genotypeInfos, "./.");

                sv[SVTYPE_DUPLICATION()].push_back(record);
            }
            itBreakpoint = mergedBreakpoint->removeBreakpoint(bp);
        }
        else
            ++itBreakpoint;
    }

    return true;
}

bool SVManager::findInversion(void)
{
    MergedCandidate* mergedBreakpoint = this->bpManager->getMergedBreakpoint();
    TBreakpointSet* candidateSet = mergedBreakpoint->getCandidateSet();
    auto itBreakpoint = candidateSet->begin();
    int32_t nID = 1;

    while (itBreakpoint != candidateSet->end())
    {
        Breakpoint* bp = *itBreakpoint;
        ReadSupportInfo* info =  mergedBreakpoint->getReadSupport(bp);
        FinalBreakpointInfo* finalBreakpoint = this->bpManager->getFinalBreakpointInfo(bp);

        bool found = true;
        /* DEBUG
        bool bFound = false;
        for (auto itit = bp->suppReads.begin(); itit != bp->suppReads.end(); ++itit)
        {
            if (*itit == "simulated.19042/2")
            {
                std::cerr << "findInversion" << "\n";
                BreakpointCandidate::printBreakpoint(bp);
                std::cerr << bp->bFoundExactPosition << "\n";
                bFound = true;
            }
        }
        //*/

        // not in the same template
        if(bp->leftTemplateID == BreakpointEvidence::NOVEL_TEMPLATE || bp->leftTemplateID != bp->rightTemplateID)
            found = false;

        // not an inversion
        if (bp->orientation != BreakpointEvidence::ORIENTATION::INVERSED)
            found = false;
        
        // not an inversion
        if (finalBreakpoint->leftPosition >= finalBreakpoint->rightPosition)
            found = false;

        int32_t svLen = finalBreakpoint->rightPosition - finalBreakpoint->leftPosition + 1;
        if (svLen < this->opManager->getMinSVSize())
            found = false;

        if (found == true)
        {
            VcfRecordEnhanced record;
            double re;
            if (info->leftReadDepthSelected == true)
                re = info->leftReadDepthDiffScore;
            else
                re = info->rightReadDepthDiffScore;

            record.rID = finalBreakpoint->leftTemplateID;
            record.beginPos = finalBreakpoint->leftPosition;
            record.endPos = finalBreakpoint->rightPosition;
            record.id = SVTYPE_INVERSION() + "_" + std::to_string(nID++);
            
            record.alt = "<" + SVTYPE_INVERSION() + ">";
            record.info  = "RT=" + std::to_string(finalBreakpoint->score) + ";";
            record.info += "SE=" + std::to_string(info->splitReadSupport) + ";";
            record.info += "PE=" + std::to_string(info->pairedEndSupport) + ";";
            record.info += "CE=" + std::to_string(info->clippedReadSupport) + ";";
            record.info += "RE=" + std::to_string(re) + ";";
            record.info += "RD=" + std::to_string(info->avgReadDepth)+ ";";
            record.info += "SVTYPE=" + SVTYPE_INVERSION() + ";";
            record.info += "SVLEN=" + std::to_string(svLen);

            // TODOs
            record.ref = "N"; // before SV
            record.qual = VcfRecord::MISSING_QUAL();
            record.filter = "PASS";       
            record.format = "GT";
            appendValue(record.genotypeInfos, "./.");

            sv[SVTYPE_INVERSION()].push_back(record);
            itBreakpoint = mergedBreakpoint->removeBreakpoint(bp);
        }
        else
            ++itBreakpoint;
    }
}

bool SVManager::findDeletion(void)
{
    MergedCandidate* mergedBreakpoint = this->bpManager->getMergedBreakpoint();
    TBreakpointSet* candidateSet = mergedBreakpoint->getCandidateSet();
    auto itBreakpoint = candidateSet->begin();
    int32_t nID = 1;

    while (itBreakpoint != candidateSet->end())
    {
        Breakpoint* bp = *itBreakpoint;
        ReadSupportInfo* info =  mergedBreakpoint->getReadSupport(bp);
        FinalBreakpointInfo* finalBreakpoint = this->bpManager->getFinalBreakpointInfo(bp);

        bool found = true;

        // not in the same template
        if(bp->leftTemplateID == BreakpointEvidence::NOVEL_TEMPLATE || bp->leftTemplateID != bp->rightTemplateID)
            found = false;

        // not a deletion
        if (bp->orientation != BreakpointEvidence::ORIENTATION::PROPERLY_ORIENTED)
            found = false;
        
        // not a deletion
        if (finalBreakpoint->leftPosition >= finalBreakpoint->rightPosition)
            found = false;

        int32_t svLen = finalBreakpoint->rightPosition - finalBreakpoint->leftPosition + 1;
        if (svLen < this->opManager->getMinSVSize())
            found = false;

        if (found == true)
        {
            VcfRecordEnhanced record;

            record.isPseudoDeletion = info->pseudoDeletion;
            record.beginPos = finalBreakpoint->leftPosition;
            record.endPos = finalBreakpoint->rightPosition;
            record.breakpoint = bp;
            record.se = info->splitReadSupport;
            record.pe = info->pairedEndSupport;
            record.ce = info->clippedReadSupport;
            if (info->leftReadDepthSelected == true)
                record.re = info->leftReadDepthDiffScore;
            else
                record.re = info->rightReadDepthDiffScore;
           
            record.rID = finalBreakpoint->leftTemplateID;
            record.id = SVTYPE_DELETION() + "_" + std::to_string(nID++);
            record.alt = "<" + SVTYPE_DELETION() + ">";
            record.info  = "RT=" + std::to_string(finalBreakpoint->score) + ";";
            record.info += "SE=" + std::to_string(record.se) + ";";
            record.info += "PE=" + std::to_string(record.pe) + ";";
            record.info += "CE=" + std::to_string(record.ce) + ";";
            record.info += "RE=" + std::to_string(record.re) + ";";
            record.info += "RD=" + std::to_string(info->avgReadDepth)+ ";";
            record.info += "SVTYPE=" + SVTYPE_DELETION() + ";";
            record.info += "SVLEN=-" + std::to_string(svLen);

            // TODOs
            record.ref = "N"; // before SV
            record.qual = VcfRecord::MISSING_QUAL();
            record.filter = "PASS";       
            record.format = "GT";
            appendValue(record.genotypeInfos, "./.");

            sv[SVTYPE_DELETION()].push_back(record);
            itBreakpoint = mergedBreakpoint->removeBreakpoint(bp);
        }
        else
            ++itBreakpoint;
    }


    /*
    std::cout << this->alnManager->getRefName(finalBreakpoint->leftTemplateID) << "\t";
    if (bp->leftReverseFlag == bp->rightReverseFlag)
       std::cout << "=";
    else
        std::cout << "!=";
    std::cout << "\t";

    std::cout << finalBreakpoint->leftPosition + 1<< "\t";
    std::cout << finalBreakpoint->rightPosition + 1<< "\t";

    std::cout << finalBreakpoint->score << "\t";        
    std::cout << info->splitReadSupport.size() << "\t";
    std::cout << info->pairedEndSupport.size() << "\t";
    std::cout << info->clippedReadSupport.size() << "\t";
    std::cout << info->readDepth << "\t";
    std::cout << info->readDepthDiffScore << "\t";

    std::set<CharString>::iterator it;
    for(it = bp->suppReads.begin(); it != (bp->suppReads.end()); ++it)
        std::cout << *it << "|";
    std::cout << "#\t";

    for(it = info->splitReadSupport.begin(); it != (info->splitReadSupport.end()); ++it)
        std::cout << *it << "|";
    std::cout << "#\t";

    for(it = info->pairedEndSupport.begin(); it != (info->pairedEndSupport.end()); ++it)
        std::cout << *it << "|,";
    std::cout << "#\t";

    for(it = info->clippedReadSupport.begin(); it != (info->clippedReadSupport.end()); ++it)
        std::cout << *it << "|,";
    std::cout << "#\n";
    */
    

    return true;
}


bool SVManager::findBreakend(void)
{
    MergedCandidate* mergedBreakpoint = this->bpManager->getMergedBreakpoint();
    TBreakpointSet* candidateSet = mergedBreakpoint->getCandidateSet();
    auto itBreakpoint = candidateSet->begin();
    int32_t nID = 1;

    while (itBreakpoint != candidateSet->end())
    {
        Breakpoint* bp = *itBreakpoint;
        ReadSupportInfo* info =  mergedBreakpoint->getReadSupport(bp);
        FinalBreakpointInfo* finalBreakpoint = this->bpManager->getFinalBreakpointInfo(bp);

        if (bp->orientation != BreakpointEvidence::ORIENTATION::NOT_DECIDED)
        {
            ++itBreakpoint;
            continue;
        }

        VcfRecordEnhanced record;
        record.breakpoint = bp;
        record.se = info->splitReadSupport;
        record.pe = info->pairedEndSupport;
        record.ce = info->clippedReadSupport;
        if (info->leftReadDepthSelected == true)
            record.re = info->leftReadDepthDiffScore;
        else
            record.re = info->rightReadDepthDiffScore;

        record.info  = "RT=" + std::to_string(finalBreakpoint->score) + ";";
        record.info += "SE=" + std::to_string(record.se) + ";";
        record.info += "PE=" + std::to_string(record.pe) + ";";
        record.info += "CE=" + std::to_string(record.ce) + ";";
        record.info += "RE=" + std::to_string(record.re) + ";";
        record.info += "RD=" + std::to_string(info->avgReadDepth)+ ";";
        record.info += "SVTYPE=" + SVTYPE_BREAKEND() + ";";

        record.qual = VcfRecord::MISSING_QUAL();
        record.filter = "PASS";     
        record.format = "GT";
        appendValue(record.genotypeInfos, "./.");

        record.ref = "N"; // TODO : get the real nt.
        std::string ntAtSV;
        if (finalBreakpoint->isLeftReverse == finalBreakpoint->isRightReverse)
        {
            ntAtSV = "N"; // TODO : get the real nt.
            record.rID = finalBreakpoint->leftTemplateID;
            record.beginPos = finalBreakpoint->leftPosition;
            record.id = SVTYPE_BREAKEND() + "_" + std::to_string(nID) + "_" + "1";
            record.alt = ntAtSV + "[" + CharStringToStdString(this->alnManager->getRefName(finalBreakpoint->rightTemplateID));
            record.alt += ":" + std::to_string(finalBreakpoint->rightPosition + 1) + "[";
            sv[SVTYPE_BREAKEND()].push_back(record);

            ntAtSV = "N";
            record.rID = finalBreakpoint->rightTemplateID;
            record.beginPos = finalBreakpoint->rightPosition;
            record.id = SVTYPE_BREAKEND() + "_" + std::to_string(nID) + "_" + "2";
            record.alt = "]" + CharStringToStdString(this->alnManager->getRefName(finalBreakpoint->leftTemplateID));
            record.alt += ":" + std::to_string(finalBreakpoint->leftPosition + 1) + "]" + ntAtSV;
            sv[SVTYPE_BREAKEND()].push_back(record);
        }
        else
        {
            ntAtSV = "N"; // TODO : get the real nt.
            record.rID = finalBreakpoint->leftTemplateID;
            record.beginPos = finalBreakpoint->leftPosition;
            record.id = SVTYPE_BREAKEND() + "_" + std::to_string(nID) + "_" + "1";
            record.alt = ntAtSV + "]" + CharStringToStdString(this->alnManager->getRefName(finalBreakpoint->rightTemplateID));
            record.alt += ":" + std::to_string(finalBreakpoint->rightPosition + 1) + "]";
            sv[SVTYPE_BREAKEND()].push_back(record);

            ntAtSV = "N"; // TODO : get the real nt.
            record.rID = finalBreakpoint->rightTemplateID;
            record.beginPos = finalBreakpoint->rightPosition;
            record.id = SVTYPE_BREAKEND() + "_" + std::to_string(nID) + "_" + "2";
            record.alt = "[" + CharStringToStdString(this->alnManager->getRefName(finalBreakpoint->leftTemplateID));
            record.alt += ":" + std::to_string(finalBreakpoint->leftPosition + 1) + "[" + ntAtSV;
            sv[SVTYPE_BREAKEND()].push_back(record);
        }

        ++nID;
        itBreakpoint = mergedBreakpoint->removeBreakpoint(bp);
    }
}