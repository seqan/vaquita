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
// ===========================================/===============================
#include "vaquita.hpp"
#include "sv.hpp"
#include "misc.hpp"

bool SVManager::findSV(void)
{
    bool result;
    RUN(result, "FIND INVERSIONS", findInversion());
    RUN(result, "FIND DELETIONS", findDeletion());
    RUN(result, "FIND DUPLICATIONS", findDuplication());
    RUN(result, "FIND TRANSLOCATION", findTranslocation());
    RUN(result, "FIND BREAKEND", findBreakend());

    filterImpreciseDel();
    return true;
}

bool SVManager::filterImpreciseDel(void)
{
    double cutoff = this->bpManager->getDepthTH();
    auto itDel = this->sv[SVTYPE_DELETION()].begin();
    unsigned filtered = 0;
    while( itDel != this->sv[SVTYPE_DELETION()].end() )
    {
        if (itDel->imprecise == true)
        {
            if ( itDel->status == VcfRecordEnhanced::STATUS::PASS && itDel->rd > cutoff)
            {
                ++filtered;
                itDel->status = VcfRecordEnhanced::STATUS::FILTERED;
            }
        }
        ++itDel;
    }
    printTimeMessage("Depth outliers: >" + std::to_string(cutoff));
    printTimeMessage("Filtered imprecise deletions: " + std::to_string(filtered));

    return true;
}

bool SVManager::loadVcf(std::string& fileName, bool useAll)
{
    VcfFileIn vcfIn(toCString(fileName));

    // header
    VcfHeader header;
    readHeader(header, vcfIn);

    // eg. RD=0.0;SVTYPE=DEL
    std::string deli1(";");
    std::string deli2("=");

    while (!atEnd(vcfIn))
    {
        VcfRecordEnhanced record;
        readRecord(record, vcfIn);

        // split by ";"
        std::vector<std::string> info;
        std::string s = CharStringToStdString(record.info);
        splitString(info, s, deli1);

        // split by "="
        std::map<std::string, std::string> mapInfo;
        for (auto it = info.begin(); it != info.end(); ++it)
        {
            std::vector<std::string> info2;
            splitString(info2, *it, deli2);

            mapInfo[info2[0]] = "";
            if (info2.size() > 1)
                mapInfo[info2[0]] = info2[1];
        }

        // fill information
        if (mapInfo.find("SVLEN") != mapInfo.end())
            record.endPos = record.beginPos + abs(std::stoll(mapInfo["SVLEN"])) - 1;
        else
            record.endPos = 0;

        if (mapInfo.find("TARGETPOS") != mapInfo.end())
            record.targetPos = std::stoul(mapInfo["TARGETPOS"]);
        if (mapInfo.find("END") != mapInfo.end())
            record.endPos = std::stod(mapInfo["END"]);
        if (mapInfo.find("SE") != mapInfo.end())
            record.se = std::stoul(mapInfo["SE"]);
        if (mapInfo.find("PE") != mapInfo.end())
            record.pe = std::stoul(mapInfo["PE"]);
        if (mapInfo.find("CE") != mapInfo.end())
            record.ce = std::stoul(mapInfo["CE"]);
        if (mapInfo.find("RE") != mapInfo.end())
            record.re = std::stod(mapInfo["RE"]);
        if (mapInfo.find("SC") != mapInfo.end())
            record.sc = std::stod(mapInfo["SC"]);
        if (mapInfo.find("RD") != mapInfo.end())
            record.rd = std::stod(mapInfo["RD"]);
        if (mapInfo.find("GC") != mapInfo.end())
            record.gc = std::stod(mapInfo["GC"]);
        if (mapInfo.find("CP") != mapInfo.end())
            record.cp = std::stod(mapInfo["CP"]);
        if (mapInfo.find("VT") != mapInfo.end())
            record.vt = std::stoi(mapInfo["VT"]);

        if (CharStringToStdString(record.filter) == VcfRecordEnhanced::STATUS_PASS())
            record.filter = VcfRecordEnhanced::STATUS::PASS;
        else if (CharStringToStdString(record.filter) == VcfRecordEnhanced::STATUS_FILTERED())
            record.filter = VcfRecordEnhanced::STATUS::FILTERED;
        else if (CharStringToStdString(record.filter) == VcfRecordEnhanced::STATUS_MERGED())
            record.filter = VcfRecordEnhanced::STATUS::MERGED;

        record.imprecise = (mapInfo.find("IMPRECISE") != mapInfo.end());
        record.chrName = CharStringToStdString(contigNames(context(vcfIn))[record.rID]);

        // filter
        if (useAll == false && record.filter != VcfRecordEnhanced::STATUS::PASS)
            continue;

        // store it
        sv[mapInfo["SVTYPE"]].push_back(record);
    }

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
            // filtered result
            if (this->opManager->getReportFilteredResult() == false && itSV->status != VcfRecordEnhanced::STATUS::PASS)
                continue;

            // addtional informations
            itSV->info += "END=" + std::to_string(itSV->endPos) + ";";
            itSV->info  = "SC=" + std::to_string(itSV->sc) + ";";
            itSV->info += "VT=" + std::to_string(itSV->vt) + ";";
            itSV->info += "SE=" + std::to_string(itSV->se) + ";";
            itSV->info += "PE=" + std::to_string(itSV->pe) + ";";
            itSV->info += "CE=" + std::to_string(itSV->ce) + ";";
            itSV->info += "RE=" + std::to_string(itSV->re) + ";";
            itSV->info += "RD=" + std::to_string(itSV->rd) + ";";
            itSV->info += "GC=" + std::to_string(itSV->gc) + ";";
            itSV->info += "CP=" + std::to_string(itSV->cp) + ";";
            if (itSV->targetPos != BreakpointEvidence::INVALID_POS)
                itSV->info += "TARGETPOS=" + std::to_string(itSV->targetPos)  + ";";
            if (itSV->endPos != BreakpointEvidence::INVALID_POS)
            {
                itSV->info += "SVLEN=";
                if (itSVType->first == SVTYPE_DELETION())
                    itSV->info += "-";
                itSV->info += std::to_string(itSV->endPos - itSV->beginPos + 1)  + ";";
            }
            itSV->info += "SVTYPE=" + itSVType->first;

            // TODOs
            itSV->ref = "N"; // before SV
            itSV->format = "GT";
            itSV->qual = VcfRecord::MISSING_QUAL();

            if (itSV->status == VcfRecordEnhanced::STATUS::PASS)
                itSV->filter = VcfRecordEnhanced::STATUS_PASS();
            else if (itSV->status == VcfRecordEnhanced::STATUS::FILTERED)
                itSV->filter = VcfRecordEnhanced::STATUS_FILTERED();
            else if (itSV->status == VcfRecordEnhanced::STATUS::MERGED)
                itSV->filter = VcfRecordEnhanced::STATUS_MERGED();
            else
                itSV->filter = VcfRecordEnhanced::STATUS_FILTERED();

            // hetero & not phased
            appendValue(itSV->genotypeInfos, "1/0");

            // add
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
    appendValue(vcfHeader, seqan::VcfHeaderRecord("INFO", "<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variation (SV)\">"));
    appendValue(vcfHeader, seqan::VcfHeaderRecord("INFO", "<ID=SVLEN,Number=1,Type=Integer,Description=\"Size of structural variation compared to reference\">"));
    appendValue(vcfHeader, seqan::VcfHeaderRecord("INFO", "<ID=TARGETPOS,Number=1,Type=Integer,Description=\"Position of the newly inserted sequence in duplication or translocations\">"));
    appendValue(vcfHeader, seqan::VcfHeaderRecord("INFO", "<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">"));
    appendValue(vcfHeader, seqan::VcfHeaderRecord("INFO", "<ID=SC,Number=1,Type=Float,Description=\"Overall score\">"));
    appendValue(vcfHeader, seqan::VcfHeaderRecord("INFO", "<ID=VT,Number=1,Type=Float,Description=\"Number of evidences types supporting the SV\">"));
    appendValue(vcfHeader, seqan::VcfHeaderRecord("INFO", "<ID=SE,Number=1,Type=Integer,Description=\"Number of split-reads supporting the SV\">"));
    appendValue(vcfHeader, seqan::VcfHeaderRecord("INFO", "<ID=PE,Number=1,Type=Integer,Description=\"Number of read-pairs supporting the SV\">"));
    appendValue(vcfHeader, seqan::VcfHeaderRecord("INFO", "<ID=RE,Number=1,Type=Float,Description=\"Read depth descrepancy around the SV\">"));
    appendValue(vcfHeader, seqan::VcfHeaderRecord("INFO", "<ID=RD,Number=1,Type=Float,Description=\"Read-depth around structural variation\">"));
    appendValue(vcfHeader, seqan::VcfHeaderRecord("INFO", "<ID=GC,Number=1,Type=Float,Description=\"GC content around the SV\">"));
    appendValue(vcfHeader, seqan::VcfHeaderRecord("INFO", "<ID=CP,Number=1,Type=Float,Description=\"Shannon entropy around the SV\">"));
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

bool SVManager::addTranslocation(VcfRecordEnhanced& orgRecord)
{
    // copy
    VcfRecordEnhanced record = orgRecord;

    int nID = sv[SVTYPE_TRANSLOCATION()].size();
    record.id = SVTYPE_TRANSLOCATION() + "_" + std::to_string(nID);
    record.alt = "<" + SVTYPE_TRANSLOCATION() + ">";
    sv[SVTYPE_TRANSLOCATION()].push_back(record);

    return true;
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
            if ( itDup->rID == itDupComp->rID)
            {
                // check matched duplications to register them to translocations
                if (  BreakpointCandidate::isAdjacent(itDup->beginPos, itDupComp->targetPos, adjTol) && \
                      BreakpointCandidate::isAdjacent(itDupComp->endPos, itDup->targetPos, adjTol))
                {
                    // merge evidences
                    //record.se = std::max(itDup->se, itDupComp->se);
                    //record.pe = std::max(itDup->pe, itDupComp->pe);
                    //record.ce = std::max(itDup->ce, itDupComp->ce);
                    record.se = itDup->se + itDupComp->se;
                    record.pe = itDup->pe + itDupComp->pe;
                    record.ce = itDup->ce + itDupComp->ce;
                    record.re = std::max(itDup->re, itDupComp->re);
                    record.vt = std::max(itDup->vt, itDupComp->vt);
                    if (itDup->status == VcfRecordEnhanced::STATUS::PASS || itDupComp->status == VcfRecordEnhanced::STATUS::PASS)
                        record.status = VcfRecordEnhanced::STATUS::PASS;

                    // marking to remove
                    //dupRemoveList[itDup - this->sv[SVTYPE_DUPLICATION()].begin()] = true;
                    //dupRemoveList[itDupComp - this->sv[SVTYPE_DUPLICATION()].begin()] = true;
                    itDup->status = VcfRecordEnhanced::STATUS::MERGED;
                    itDupComp->status = VcfRecordEnhanced::STATUS::MERGED;
                    addTranslocation(record);

                }
            }
        }
    }

    // merge overlapped duplications
    for (auto itTra = this->sv[SVTYPE_TRANSLOCATION()].begin(); itTra != this->sv[SVTYPE_TRANSLOCATION()].end(); ++itTra)
    {
        if (itTra->status == VcfRecordEnhanced::STATUS::PASS)
        {
            for (auto itDup = this->sv[SVTYPE_DUPLICATION()].begin(); itDup != this->sv[SVTYPE_DUPLICATION()].end(); ++itDup)
            {
                if (itTra->rID == itDup->rID && itDup->status != VcfRecordEnhanced::STATUS::MERGED)
                {
                    if (BreakpointCandidate::isOverlap(itTra->beginPos, itTra->endPos, itDup->beginPos, itDup->endPos))
                    {
                        itTra->se += itDup->se;
                        itTra->pe += itDup->pe;
                        itTra->ce += itDup->ce;
                        //dupRemoveList[itDup - this->sv[SVTYPE_DUPLICATION()].begin()] = true;
                        itDup->status = VcfRecordEnhanced::STATUS::MERGED;
                    }
                }
            }
        }
    }

    // remove marked duplications (they are translocations now)
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

        // not an duplication
        if (bp->orientation != BreakpointEvidence::ORIENTATION::SWAPPED)
            found = false;

        // skip this
        if (finalBreakpoint->filtered == true)
            found = false;

        if (found == true)
        {
            std::vector<TPosition> leftExactPositions, centerExactPositions, rightExactPositions;
            std::vector<TPosition> leftImprecisePositions, centerImprecisePositions, rightImprecisePositions;
            bool leftMatchFound = false, rightMatchFound = false;

            // set initial positions
            if (finalBreakpoint->imprecise == false)
            {
                leftExactPositions.push_back(finalBreakpoint->leftPosition);
                rightExactPositions.push_back(finalBreakpoint->rightPosition);

            }
            else
            {
                leftImprecisePositions.push_back(finalBreakpoint->leftPosition);
                rightImprecisePositions.push_back(finalBreakpoint->rightPosition);
            }

            // find matched deletions
            std::vector<VcfRecordEnhanced*> vDel;
            auto itDel = this->sv[SVTYPE_DELETION()].begin();
            while( itDel != this->sv[SVTYPE_DELETION()].end() )
            {
                bool matchFound = false;
                uint svLen = 0;

                if (itDel->rID == finalBreakpoint->leftTemplateID)
                {
                    bool leftMatched = BreakpointCandidate::isAdjacent(itDel->beginPos, finalBreakpoint->leftPosition, adjTol);
                    bool rightMatched = BreakpointCandidate::isAdjacent(itDel->endPos, finalBreakpoint->rightPosition, adjTol);

                    if (leftMatched && rightMatched)
                    {
                        if (itDel->imprecise == false)
                        {
                            leftExactPositions.push_back(itDel->beginPos);
                            rightExactPositions.push_back(itDel->endPos);
                        }
                        else
                        {
                            leftImprecisePositions.push_back(itDel->beginPos);
                            rightImprecisePositions.push_back(itDel->endPos);
                        }
                        vDel.push_back(&(*itDel));
                    }
                    else if (leftMatched)
                    {
                        if (itDel->imprecise == false)
                        {
                            leftExactPositions.push_back(itDel->beginPos);
                            centerImprecisePositions.push_back(itDel->endPos);
                        }
                        else
                        {
                            leftImprecisePositions.push_back(itDel->beginPos);
                            centerImprecisePositions.push_back(itDel->endPos);
                        }
                        vDel.push_back(&(*itDel));
                    }
                    else if (rightMatched)
                    {
                        if (itDel->imprecise == false)
                        {
                            centerExactPositions.push_back(itDel->beginPos);
                            rightExactPositions.push_back(itDel->endPos);
                        }
                        else
                        {
                            centerImprecisePositions.push_back(itDel->beginPos);
                            rightImprecisePositions.push_back(itDel->endPos);
                        }
                        vDel.push_back(&(*itDel));
                    }
                    /*
                    else if (finalBreakpoint->imprecise == false && IS_OVERLAP(itDel->beginPos, itDel->endPos, finalBreakpoint->leftPosition, finalBreakpoint->rightPosition))
                    {
                        vDel.push_back(&(*itDel));
                    }
                    */

                    leftMatchFound = leftMatchFound || leftMatched;
                    rightMatchFound = rightMatchFound || rightMatched;
                }
                ++itDel;
            }

            // update supporting reads
            VcfRecordEnhanced record;
            record.se = info->splitReadSupport;
            record.pe = info->pairedEndSupport;
            record.ce = info->clippedReadSupport;
            record.re = std::max(info->leftReadDepthDiffScore, info->rightReadDepthDiffScore);
            record.rd = info->avgReadDepth;
            record.vt = finalBreakpoint->vote;
            record.gc = finalBreakpoint->gcContent;
            record.cp = finalBreakpoint->sequenceComplexity;

            if (finalBreakpoint->filtered == false)
            {
                record.status = VcfRecordEnhanced::STATUS::PASS;
                record.imprecise = finalBreakpoint->imprecise;

                // merged to this duplication
                for (auto it=vDel.begin(); it != vDel.end(); ++it)
                    (*it)->status = VcfRecordEnhanced::STATUS::MERGED;

                // sorting
                std::sort(leftExactPositions.begin(), leftExactPositions.end());
                std::sort(leftImprecisePositions.begin(), leftImprecisePositions.end());
                std::sort(centerExactPositions.begin(), centerExactPositions.end());
                std::sort(centerImprecisePositions.begin(), centerImprecisePositions.end());
                std::sort(rightExactPositions.begin(), rightExactPositions.end());
                std::sort(rightImprecisePositions.begin(), rightImprecisePositions.end());

                // fill records
                record.id = SVTYPE_DUPLICATION() + "_" + std::to_string(nID++);
                record.rID = finalBreakpoint->leftTemplateID;
                record.breakpoint = bp;

                // get positions
                TPosition leftPos, rightPos;
                if (leftExactPositions.size() > 0)
                {
                    leftPos = MID_ELEMENT(leftExactPositions);
                }
                else
                {
                    leftPos = MID_ELEMENT(leftImprecisePositions);
                    record.imprecise = true;
                }

                if (rightExactPositions.size() > 0)
                {
                    rightPos = MID_ELEMENT(rightExactPositions);
                }
                else
                {
                    rightPos = MID_ELEMENT(rightImprecisePositions);
                    record.imprecise = true;
                }

                // tandem duplication
                if (centerExactPositions.size() == 0 && centerImprecisePositions.size() == 0)
                {
                    record.alt = "<" + SVTYPE_DUPLICATION() + ":TANDEM>";
                    record.beginPos = leftPos;
                    record.endPos = rightPos;
                    record.targetPos = BreakpointEvidence::INVALID_POS;

                    // size check
                    if ((record.endPos - record.beginPos + 1) >= this->opManager->getMinSVSize())
                        sv[SVTYPE_DUPLICATION()].push_back(record);
                }
                else
                {
                    record.alt = "<" + SVTYPE_DUPLICATION() + ">";

                    // center
                    TPosition centerPos;
                    if (centerExactPositions.size() > 0)
                    {
                        centerPos = MID_ELEMENT(centerExactPositions);
                    }
                    else
                    {
                        centerPos = MID_ELEMENT(centerImprecisePositions);
                        record.imprecise = true;
                    }

                    if (leftMatchFound && rightMatchFound) // translocation
                    {
                        record.beginPos = leftPos;
                        record.endPos = centerPos;
                        record.targetPos = rightPos;
                        addTranslocation(record);

                        record.status = VcfRecordEnhanced::STATUS::MERGED;
                        sv[SVTYPE_DUPLICATION()].push_back(record);
                    }
                    else if (leftMatchFound)
                    {
                        record.beginPos = centerPos;
                        record.endPos = rightPos;
                        record.targetPos = leftPos;

                        // size check
                        if ((record.endPos - record.beginPos + 1) >= this->opManager->getMinSVSize())
                            sv[SVTYPE_DUPLICATION()].push_back(record);
                    }
                    else if (rightMatchFound)
                    {
                        record.beginPos = leftPos;
                        record.endPos = centerPos;
                        record.targetPos = rightPos;

                        // size check
                        if ((record.endPos - record.beginPos + 1) >= this->opManager->getMinSVSize())
                            sv[SVTYPE_DUPLICATION()].push_back(record);
                    }
                }

            }
            itBreakpoint = mergedBreakpoint->removeBreakpoint(bp);
        }
        else
            ++itBreakpoint;
    } // end while

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
        // not in the same template
        if(bp->leftTemplateID == BreakpointEvidence::NOVEL_TEMPLATE || bp->leftTemplateID != bp->rightTemplateID)
            found = false;

        // not an inversion
        if (bp->orientation != BreakpointEvidence::ORIENTATION::INVERTED)
            found = false;

        // not an inversion
        if (finalBreakpoint->leftPosition >= finalBreakpoint->rightPosition)
            found = false;

        // skip this
        if (finalBreakpoint->filtered == true)
            found = false;;

        // size check
        int32_t svLen = finalBreakpoint->rightPosition - finalBreakpoint->leftPosition + 1;
        if (svLen < this->opManager->getMinSVSize())
            found = false;

        if (found == true)
        {
            VcfRecordEnhanced record;

            record.id = SVTYPE_INVERSION() + "_" + std::to_string(nID++);
            record.alt = "<" + SVTYPE_INVERSION() + ">";

            record.rID = finalBreakpoint->leftTemplateID;
            record.beginPos = finalBreakpoint->leftPosition;
            record.endPos = finalBreakpoint->rightPosition;
            record.targetPos = BreakpointEvidence::INVALID_POS;

            record.se = info->splitReadSupport;
            record.pe = info->pairedEndSupport;
            record.ce = info->clippedReadSupport;
            record.re = std::max(info->leftReadDepthDiffScore, info->rightReadDepthDiffScore);
            record.rd = info->avgReadDepth;
            record.vt = finalBreakpoint->vote;
            record.gc = finalBreakpoint->gcContent;
            record.cp = finalBreakpoint->sequenceComplexity;

            record.breakpoint = bp;
            record.status = (finalBreakpoint->filtered) ? VcfRecordEnhanced::STATUS::FILTERED : VcfRecordEnhanced::STATUS::PASS;
            record.imprecise = finalBreakpoint->imprecise;

            sv[SVTYPE_INVERSION()].push_back(record);
            itBreakpoint = mergedBreakpoint->removeBreakpoint(bp);
        }
        else
            ++itBreakpoint;
    }

    return true;
}

bool SVManager::findDeletion(void)
{
    MergedCandidate* mergedBreakpoint = this->bpManager->getMergedBreakpoint();
    TBreakpointSet* candidateSet = mergedBreakpoint->getCandidateSet();
    auto itBreakpoint = candidateSet->begin();
    int32_t nID = 1;

    TTemplateID rID_before = 0;
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
        if (bp->orientation != BreakpointEvidence::ORIENTATION::PROPERLY_ORIENTED_LARGE)
            found = false;

        // not a deletion
        if (finalBreakpoint->leftPosition >= finalBreakpoint->rightPosition)
            found = false;

        // size check
        int32_t svLen = finalBreakpoint->rightPosition - finalBreakpoint->leftPosition + 1;
        if (svLen < this->opManager->getMinSVSize())
            found = false;

        if (found == true)
        {
            VcfRecordEnhanced record;

            record.id = SVTYPE_DELETION() + "_" + std::to_string(nID++);
            record.alt = "<" + SVTYPE_DELETION() + ">";

            record.rID = finalBreakpoint->leftTemplateID;
            record.beginPos = finalBreakpoint->leftPosition;
            record.endPos = finalBreakpoint->rightPosition;
            record.targetPos = BreakpointEvidence::INVALID_POS;

            record.se = info->splitReadSupport;
            record.pe = info->pairedEndSupport;
            record.ce = info->clippedReadSupport;
            record.re = std::max(info->leftReadDepthDiffScore, info->rightReadDepthDiffScore);
            record.rd = info->avgReadDepth;
            record.vt = finalBreakpoint->vote;
            record.gc = finalBreakpoint->gcContent;
            record.cp = finalBreakpoint->sequenceComplexity;

            record.breakpoint = bp;
            record.status = (finalBreakpoint->filtered) ? VcfRecordEnhanced::STATUS::FILTERED : VcfRecordEnhanced::STATUS::PASS;
            record.imprecise = finalBreakpoint->imprecise;

            sv[SVTYPE_DELETION()].push_back(record);
            itBreakpoint = mergedBreakpoint->removeBreakpoint(bp);
        }
        else
            ++itBreakpoint;
    }

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

        // skip this
        if (finalBreakpoint->filtered == true)
        {
            ++itBreakpoint;
            continue;
        }

        VcfRecordEnhanced record;
        record.breakpoint = bp;
        record.targetPos = BreakpointEvidence::INVALID_POS;
        record.endPos = BreakpointEvidence::INVALID_POS;
        record.se = info->splitReadSupport;
        record.pe = info->pairedEndSupport;
        record.ce = info->clippedReadSupport;
        record.re = std::max(info->leftReadDepthDiffScore, info->rightReadDepthDiffScore);
        record.rd = info->avgReadDepth;
        record.gc = finalBreakpoint->gcContent;
        record.cp = finalBreakpoint->sequenceComplexity;

        record.status = (finalBreakpoint->filtered) ? VcfRecordEnhanced::STATUS::FILTERED : VcfRecordEnhanced::STATUS::PASS;
        record.imprecise = finalBreakpoint->imprecise;

        // represented by 2 records
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

    return true;
}

bool SVManager::orderSV(void)
{
    if (this->opManager->getUseRankAggregation())
        return this->orderSVByRankAgg();
    else
        return this->orderSVByEvidenceSum();

    return true;
}

bool SVManager::orderSVByEvidenceSum(void)
{
    // for each type and SVs
    for (auto itSVType = this->sv.begin(); itSVType != this->sv.end(); ++itSVType)
        for (auto itSV = itSVType->second.begin(); itSV != itSVType->second.end(); ++itSV)
            itSV->sc = itSV->se + itSV->pe + itSV->ce;

    return true;
}

bool SVManager::orderSVByRankAgg(void)
{
    typedef std::pair<std::pair<double, double>, VcfRecordEnhanced*> TOrderingInfo;
    typedef std::map<VcfRecordEnhanced*, std::vector<unsigned int> > TRankInfo;
    std::vector< TOrderingInfo > vSe, vPe, vRe;
    std::vector< std::pair<unsigned int, VcfRecordEnhanced*> > aggRank;

    // for each type
    for (auto itSVType = this->sv.begin(); itSVType != this->sv.end(); ++itSVType)
    {
        // init.
        /*
        vSe.clear();
        vPe.clear();
        vRe.clear();
        aggRank.clear();
        */

        // fill list
        for (auto itSV = itSVType->second.begin(); itSV != itSVType->second.end(); ++itSV)
        {
            if (itSV->status == VcfRecordEnhanced::STATUS::PASS)
            {
                vSe.push_back(std::make_pair(std::make_pair(itSV->se + itSV->ce, itSV->rd), &(*itSV)));
                vPe.push_back(std::make_pair(std::make_pair(itSV->pe, itSV->rd), &(*itSV)));
                vRe.push_back(std::make_pair(std::make_pair(itSV->re, itSV->rd), &(*itSV)));

                /*
                double se = (double) (itSV->se + itSV->ce + BreakpointCandidate::PREVENT_DIV_BY_ZERO()) / itSV->rd;
                double pe = (double) (itSV->pe + BreakpointCandidate::PREVENT_DIV_BY_ZERO()) / itSV->rd;
                double re = (double) (itSV->re + BreakpointCandidate::PREVENT_DIV_BY_ZERO()) / itSV->rd;
                vSe.push_back(std::make_pair(std::make_pair(se, se), &(*itSV)));
                vPe.push_back(std::make_pair(std::make_pair(pe, pe), &(*itSV)));
                vRe.push_back(std::make_pair(std::make_pair(re, re), &(*itSV)));
                */
            }
            else
            {
                itSV->sc = -1;
            }
        }
    }

    {
        // get individual rank.
        auto sorter = [](TOrderingInfo l, TOrderingInfo r)->bool{
        if (l.first.first != r.first.first)
            return l.first.first > r.first.first;
        else
            return l.first.second < r.first.second;
        };
        std::sort(vSe.begin(), vSe.end(), sorter);
        std::sort(vPe.begin(), vPe.end(), sorter);
        std::sort(vRe.begin(), vRe.end(), sorter);

        // aggregate ranks
        TRankInfo mapRanks;
        for (auto rank=0; rank < vSe.size(); ++rank)
        {
            if(mapRanks.find(vSe[rank].second) == mapRanks.end())
                mapRanks.insert(std::make_pair(vSe[rank].second, std::vector<unsigned int>()));
            mapRanks[vSe[rank].second].push_back(rank);

            if (this->opManager->doPairedEndAnalysis())
            {
                if(mapRanks.find(vPe[rank].second) == mapRanks.end())
                    mapRanks.insert(std::make_pair(vPe[rank].second, std::vector<unsigned int>()));
                 mapRanks[vPe[rank].second].push_back(rank);
            }

            if (this->opManager->doReadDepthAnalysis())
            {
                if(mapRanks.find(vRe[rank].second) == mapRanks.end())
                    mapRanks.insert(std::make_pair(vRe[rank].second, std::vector<unsigned int>()));
                mapRanks[vRe[rank].second].push_back(rank);
            }
        }

        // get median rank
        auto it = mapRanks.begin();
        while (it != mapRanks.end())
        {
            std::sort(it->second.begin(), it->second.end());
            aggRank.push_back(std::make_pair(MID_ELEMENT(it->second), it->first));
            ++it;
        }

        // final
        std::sort(aggRank.begin(), aggRank.end(), [](auto &left, auto &right) {return left.first < right.first;} );
        int correctedRank = -1, prevOriginalRank = -1, tieCount = 0;
        for (unsigned i=0; i < aggRank.size(); ++i)
        {
            // tie
            if (prevOriginalRank == aggRank[i].first)
            {
                ++tieCount;
            }
            else
            {
                correctedRank += tieCount + 1;
                prevOriginalRank = aggRank[i].first;
                tieCount = 0;
            }
            aggRank[i].second->sc = (1.0 - ( (double)correctedRank  /  (double) (aggRank.size() - 1))) * 100;
        }
    }

    return true;
}

uint32_t SVManager::getSVCount(std::string svType, bool countFilteredResult = false)
{
    int nCnt = 0;
    for (auto itSV = sv[svType].begin(); itSV != sv[svType].end(); ++itSV)
    {
        if (countFilteredResult == true || itSV->status == VcfRecordEnhanced::STATUS::PASS)
            ++nCnt;
    }
    return nCnt;
}

bool SVManager::addSV(std::string& svType, VcfRecordEnhanced& record)
{
    sv[svType].push_back(record);
    return true;
}
