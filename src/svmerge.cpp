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
#include "svmerge.hpp"
#include "misc.hpp"
#include "intervalindex.hpp"

bool SVMerge::loadVcf(std::string& fileName)
{
    printTimeMessage("Loading: " + fileName);
    
    SVManager newSV(fileName, optionManager->getUseAll());
    printTimeMessage(std::to_string(newSV.getDeletionCount()) + " deletions.");
    printTimeMessage(std::to_string(newSV.getInversionCount()) + " inversions.");
    printTimeMessage(std::to_string(newSV.getDuplicationCount()) + " duplications.");
    printTimeMessage(std::to_string(newSV.getTranslocationCount()) + " translocations.");
    this->svSet.push_back(newSV);
    
    return true;
}

bool SVMerge::loadVcfs(std::vector<std::string>& inputFileNames)
{
    for(auto it = inputFileNames.begin(); it != inputFileNames.end(); ++it)
    {
        if (this->loadVcf(*it) == false)
            return false;
        this->fileNames.push_back(*it);
    }
    return true;
}


bool SVMerge::loadVcfs(void)
{
    std::vector<std::string>& fileNames = optionManager->getVcfs();
    return loadVcfs(fileNames);
}

bool SVMerge::merge(void)
{
    if (this->svSet.size() == 0)
        return false;

    std::map<std::string, std::vector<std::vector<VcfRecordEnhanced>>> clusteredSV;
    std::map<std::string, unsigned> refNameMap;

    // for each dataset
    for (auto itSet = svSet.begin(); itSet != svSet.end(); ++itSet)
    {
        unsigned datasetIdx = (itSet - svSet.begin());

       // for sach SV type
        TVcfBySVType& vcfBySVType = itSet->getVcfBySVType();
        for (auto itSVType = vcfBySVType.begin(); itSVType != vcfBySVType.end(); ++itSVType)
        {
            if (itSVType->first == SVManager::SVTYPE_BREAKEND())
                continue;

            // for each SV
            for (auto itRecord = itSVType->second.begin(); itRecord != itSVType->second.end(); ++itRecord)
            {
                bool bFound = false;
                itRecord->additionalInfo = datasetIdx;
                refNameMap[itRecord->chrName] = 0;

                // for each cluster
                for (auto itClustSV = clusteredSV[itSVType->first].begin(); itClustSV != clusteredSV[itSVType->first].end(); ++itClustSV)
                {
                    // check overlap
                    for (auto itClustSVRecord = itClustSV->begin(); itClustSVRecord != itClustSV->end(); ++itClustSVRecord)
                    {
                        // diff. chromosomes.
                        if (itClustSVRecord->chrName != itRecord->chrName)
                            break;

                        // check overlap & add
                        if (isReciOverlap(itClustSVRecord->beginPos, itClustSVRecord->endPos, itRecord->beginPos, itRecord->endPos, 0.8))
                        {
                            itClustSV->push_back(*itRecord);
                            bFound = true;
                            break;
                        }
                    }

                    if (bFound == true)
                        break;
                }

                // new one
                if (bFound == false)
                {
                    std::vector<VcfRecordEnhanced> firstSV;
                    firstSV.push_back(*itRecord);
                    clusteredSV[itSVType->first].push_back(firstSV);
                }
            }
        }
    }

    // ref names
    for (auto it = refNameMap.begin(); it != refNameMap.end(); ++it)
        this->refNames.push_back(it->first);
    std::sort(this->refNames.begin(), this->refNames.end());
    for (auto it = this->refNames.begin(); it != this->refNames.end(); ++it)
        refNameMap[*it] = it - this->refNames.begin();

    for (auto itType = clusteredSV.begin(); itType != clusteredSV.end(); ++itType)
    {
        for (auto itCluster = itType->second.begin(); itCluster != itType->second.end(); ++itCluster)
        {
            VcfRecordMultiSample record;
            record.rID = refNameMap[itCluster->begin()->chrName];

            // init.
            for (auto itSet = svSet.begin(); itSet != svSet.end(); ++itSet)
            {
                record.gt.push_back("0/0");
                record.dp.push_back(0);
                record.recordNum.push_back(0);
            }

            // extract info.
            std::vector<TPosition> beginVector, endVector, targetVector;
            std::vector<double> sc;
            for (auto itSV = itCluster->begin(); itSV != itCluster->end(); ++itSV)
            {
                // positions
                beginVector.push_back(itSV->beginPos);
                endVector.push_back(itSV->endPos);
                targetVector.push_back(itSV->targetPos);

                // genotype
                if (record.gt[itSV->additionalInfo] == "0/0")
                    record.gt[itSV->additionalInfo] = CharStringToStdString(itSV->genotypeInfos[0]);

                // depth
                record.dp[itSV->additionalInfo] += (unsigned) itSV->rd;
                ++record.recordNum[itSV->additionalInfo];
                sc.push_back(itSV->sc);
            }

            // merged record
            std::sort(beginVector.begin(), beginVector.end());
            std::sort(endVector.begin(), endVector.end());
            std::sort(targetVector.begin(), targetVector.end());
            std::sort(sc.begin(), sc.end());

            record.beginPos = MID_ELEMENT(beginVector); // representative pos.
            record.endPos = MID_ELEMENT(endVector);
            record.inBeginPos = LAST_ELEMENT(beginVector); // inside
            record.inEndPos = FIRST_ELEMENT(endVector);
            record.outBeginPos = FIRST_ELEMENT(beginVector); // outside
            record.outEndPos = LAST_ELEMENT(endVector);
            record.minSC = FIRST_ELEMENT(sc);
            record.maxSC = LAST_ELEMENT(sc);

            if (targetVector.size() > 0)
            {
                record.targetPos = MID_ELEMENT(targetVector);
                record.inTargetPos = FIRST_ELEMENT(targetVector);
                record.outTargetPos = LAST_ELEMENT(targetVector);
            }
            else
            {
                record.targetPos = BreakpointEvidence::INVALID_POS;
                record.inTargetPos = BreakpointEvidence::INVALID_POS;
                record.outTargetPos = BreakpointEvidence::INVALID_POS;
            }

            // average depth
            for (unsigned i = 0; i < svSet.size(); ++i)
            {
                if (record.recordNum[i] > 0)
                    record.dp[i] = (unsigned) (record.dp[i] / record.recordNum[i]);
            }

            this->mergedSV[itType->first].push_back(record);
        }
    }

    return true;
}

bool SVMerge::isReciOverlap(TPosition beginPos1, TPosition endPos1, TPosition beginPos2, TPosition endPos2, double tol)
{
    if (IntervalIndex<TPosition>::isOverlap(beginPos1, endPos1, beginPos2, endPos2))
    {
        TPosition size1 = endPos1 - beginPos1 + 1;
        TPosition size2 = endPos2 - beginPos2 + 1;

        std::vector<TPosition> positions;
        positions.push_back(beginPos1);
        positions.push_back(endPos1);
        positions.push_back(beginPos2);
        positions.push_back(endPos2);
        std::sort(positions.begin(), positions.end());
        TPosition matchedSize = positions[2] - positions[1] + 1;

        return (((double) matchedSize >= (size1 * tol)) && ((double) matchedSize >= (size2 * tol)));
    }

    return false;
}

bool SVMerge::writeVCF(void)
{
    // merge to a single list
    std::vector<VcfRecordEnhanced> vcfRecords;
    for (auto itSVType = this->mergedSV.begin(); itSVType != this->mergedSV.end(); ++itSVType)
    {    
        std::sort(itSVType->second.begin(), itSVType->second.end(), less_than_vcf());
        int32_t nID = 1;
        for (auto itSV = itSVType->second.begin(); itSV != itSVType->second.end(); ++itSV)
        {
            itSV->id = itSVType->first + "_" + std::to_string(nID++);
            itSV->alt = "<" + itSVType->first + ">";

            // addtional informations
            itSV->info  = "MIN_SC=" + std::to_string(itSV->minSC) + ";";
            itSV->info += "MAX_SC=" + std::to_string(itSV->maxSC) + ";";
            itSV->info += "IN_BEGINPOS=" + std::to_string(itSV->inBeginPos) + ";";
            itSV->info += "IN_ENDPOS=" + std::to_string(itSV->inEndPos) + ";";
            itSV->info += "OUT_BEGINPOS=" + std::to_string(itSV->outBeginPos) + ";";
            itSV->info += "OUT_ENDPOS=" + std::to_string(itSV->outEndPos) + ";";
            if (itSV->targetPos != BreakpointEvidence::INVALID_POS)
            {
                itSV->info += "TARGETPOS=" + std::to_string(itSV->targetPos)  + ";";
                itSV->info += "IN_TARGETPOS=" + std::to_string(itSV->inTargetPos)  + ";";
                itSV->info += "OUTTARGETPOS=" + std::to_string(itSV->outTargetPos)  + ";";
            }

            if (itSV->endPos != BreakpointEvidence::INVALID_POS)
            {
                itSV->info += "SVLEN=";
                if (itSVType->first == SVManager::SVTYPE_DELETION())
                    itSV->info += "-";
                itSV->info += std::to_string(itSV->endPos - itSV->beginPos + 1)  + ";";
            }
            itSV->info += "SVTYPE=" + itSVType->first;

            itSV->ref = "N"; // before SV
            itSV->format = "GT:DP";
            itSV->qual = VcfRecord::MISSING_QUAL();
            itSV->filter = VcfRecordEnhanced::STATUS_PASS();

            // genotype info
            for (unsigned datasetID = 0; datasetID < this->svSet.size(); ++datasetID)
                appendValue(itSV->genotypeInfos, itSV->gt[datasetID] + ":" +  std::to_string(itSV->dp[datasetID]));

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
    appendValue(vcfHeader, seqan::VcfHeaderRecord("INFO", "<ID=MIN_SC,Number=1,Type=Float,Description=\"Minimum score\">"));
    appendValue(vcfHeader, seqan::VcfHeaderRecord("INFO", "<ID=MAX_SC,Number=1,Type=Float,Description=\"Maximum score\">"));
    appendValue(vcfHeader, seqan::VcfHeaderRecord("ALT", "<ID=DEL,Description=\"Deletion\">"));
    appendValue(vcfHeader, seqan::VcfHeaderRecord("ALT", "<ID=INV,Description=\"Inversion\">"));
    appendValue(vcfHeader, seqan::VcfHeaderRecord("ALT", "<ID=DUP,Description=\"Duplication\">"));
    appendValue(vcfHeader, seqan::VcfHeaderRecord("ALT", "<ID=DUP:TANDEM,Description=\"Tandem duplication\">"));
    appendValue(vcfHeader, seqan::VcfHeaderRecord("ALT", "<ID=TRA,Description=\"Translocation\">"));
    appendValue(vcfHeader, seqan::VcfHeaderRecord("ALT", "<ID=BND,Description=\"Breakend\">"));
    appendValue(vcfHeader, seqan::VcfHeaderRecord("FORMAT", "<ID=GT,Number=1,Type=String,Description=\"Genotype\">"));

    // reference & sample
    for (auto it = this->refNames.begin(); it != this->refNames.end(); ++it)
        appendValue(contigNames(context(vcfOut)), *it);

    for(auto it = this->fileNames.begin(); it != this->fileNames.end(); ++it)
        appendValue(sampleNames(context(vcfOut)), *it);

    // write
    writeHeader(vcfOut, vcfHeader);
    for (auto itSV = vcfRecords.begin(); itSV != vcfRecords.end(); ++itSV)
        writeRecord(vcfOut, *itSV);

    return true;
}