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
#include "alignment.hpp"
#include "misc.hpp"

void AlignmentManager::getSequenceAndDepth(CharString& seq, std::vector<int32_t>& depth, TTemplateID rID, TPosition beginPos, TPosition endPos)
{
    // init
    seq = "";
    std::vector<std::map<char,int32_t> > profile;
    int32_t seqSize = (endPos - beginPos);
    profile.resize(seqSize);

    depth.clear();
    int32_t depthSize = (endPos - beginPos);
    for (int i=0; i < depthSize; ++i)
        depth.push_back(0);

        
    // Jump the BGZF stream to this position.
    bool hasAlignments = false;
    if (!jumpToRegion(this->bamFileIn, hasAlignments, rID, beginPos, endPos, baiIndex))
    {
        std::cerr << "ERROR: Could not jump to " << rID << ":" << beginPos << "-" << endPos << "\n";
        return;
    }
    if (!hasAlignments)
        return;
    
    // search
    //int32_t clippedReadEditDistance = optionManager->getClippedReadEditDistance();
    BamAlignmentRecord record;
    while ( !atEnd(this->bamFileIn) )
    {
        readRecord(record, this->bamFileIn);

         // stop retrival
        if (record.rID == -1 || record.rID > rID || record.beginPos >= endPos)
            break;

        // skip
        int32_t len = getAlignmentLengthInRef(record);
        if ( (record.beginPos + len) <= beginPos)
          continue;

        // calc. depth
        int32_t startIdx = record.beginPos - beginPos;
        int32_t endIdx = startIdx + len;

        if (startIdx < 0)
            startIdx = 0;

        if (endIdx > depthSize)
            endIdx = depthSize;

        for(int i=startIdx; i < endIdx; ++i)
            ++depth[i];

        // edit distance based filter
        /*
        BamTagsDict tagsDict(record.tags);
        int32_t editDistance, tagIdx;
        if ( findTagKey(tagIdx, tagsDict, "NM") )   
        {
            extractTagValue(editDistance, tagsDict, tagIdx);
            if (editDistance > clippedReadEditDistance)
                continue;
        }
        */

       // if (hasFlagRC(record)) // wrong
       //     reverseComplement(record.seq);
       
        // calc. seq
        int32_t seqStart = beginPos - record.beginPos;

        if (seqStart < 0)
            seqStart = 0;
        if (record.cigar[0].operation == 'S')
            seqStart += record.cigar[0].count;       

        for (int i=startIdx; i < endIdx && seqStart < length(record.seq); ++i, ++seqStart)
            ++profile[i][record.seq[seqStart]];
    }

    // conesnsus sequence
    std::map<char,int32_t>::iterator topIt;
    std::map<char,int32_t>::iterator it;
    int topCnt;
    for (int i=0; i < profile.size(); ++i)
    {
        if (profile[i].size() == 0)
            appendValue(seq, 'N');
        else if (profile[i].size() == 1)
            appendValue(seq, profile[i].begin()->first);
        else
        {
            it = profile[i].begin();
            topIt = it;
            topCnt = 0;
            while( it != profile[i].end() )
            {
                if (it->second > topCnt)
                {
                    topIt = it;
                    topCnt = it->second;
                }
                ++it;
            }
            appendValue(seq, topIt->first);
        }
    }
}

void AlignmentManager::getSequence(CharString& seq, TTemplateID rID, TPosition beginPos, TPosition endPos)
{
    // init
    clear(seq);
    std::vector<std::map<char,int32_t> > profile;
    int32_t seqSize = (endPos - beginPos);
    profile.resize(seqSize);
        
    // Jump the BGZF stream to this position.
    bool hasAlignments = false;
    if (!jumpToRegion(this->bamFileIn, hasAlignments, rID, beginPos, endPos, baiIndex))
    {
        std::cerr << "ERROR: Could not jump to " << rID << ":" << beginPos << "-" << endPos << "\n";
        return;
    }
    if (!hasAlignments)
        return;
    
    // search
    //int32_t clippedReadEditDistance = optionManager->getClippedReadEditDistance();
    BamAlignmentRecord record;
    while ( !atEnd(this->bamFileIn) )
    {
        readRecord(record, this->bamFileIn);

         // stop retrival
        if (record.rID == -1 || record.rID > rID || record.beginPos >= endPos)
            break;

        // skip
        int32_t len = getAlignmentLengthInRef(record);
        if ( (record.beginPos + len) <= beginPos)
          continue;
       
        // calc. seq
        int32_t startIdx = record.beginPos - beginPos;
        int32_t endIdx = startIdx + len;
        int32_t seqStart;

        if (startIdx < 0)
            startIdx = 0;

        if (endIdx > seqSize)
            endIdx = seqSize;

        seqStart = beginPos - record.beginPos;
        if (seqStart < 0)
            seqStart = 0;
        if (record.cigar[0].operation == 'S')
            seqStart += record.cigar[0].count;       

        for (int i=startIdx; i < endIdx && seqStart < length(record.seq); ++i, ++seqStart)
            ++profile[i][record.seq[seqStart]];
    }

    // conesnsus sequence
    std::map<char,int32_t>::iterator topIt;
    std::map<char,int32_t>::iterator it;
    int topCnt;
    for (int i=0; i < profile.size(); ++i)
    {
        if (profile[i].size() == 0)
            appendValue(seq, 'N');
        else if (profile[i].size() == 1)
            appendValue(seq, profile[i].begin()->first);
        else
        {
            it = profile[i].begin();
            topIt = it;
            topCnt = 0;
            while( it != profile[i].end() )
            {
                if (it->second > topCnt)
                {
                    topIt = it;
                    topCnt = it->second;
                }
                ++it;
            }
            appendValue(seq, topIt->first);
        }
    }
}

void AlignmentManager::getDepth(std::vector<int32_t>& depth, TTemplateID rID, TPosition beginPos, TPosition endPos)
{
    // init
    depth.clear();
    int32_t depthSize = (endPos - beginPos);
    for (int i=0; i < depthSize; ++i)
        depth.push_back(0);

    // Jump the BGZF stream to this position.
    bool hasAlignments = false;
    if (!jumpToRegion(this->bamFileIn, hasAlignments, rID, beginPos, endPos, baiIndex))
    {
        std::cerr << "ERROR: Could not jump to " << rID << ":" << beginPos << "-" << endPos << "\n";
        return;
    }
    if (!hasAlignments)
        return;

    // search
    BamAlignmentRecord record;
    while ( !atEnd(this->bamFileIn) )
    {
        readRecord(record, this->bamFileIn);

         // stop retrival
        if (record.rID == -1 || record.rID > rID || record.beginPos >= endPos)
            break;

        // skip
        int32_t len = getAlignmentLengthInRef(record);
        if ( (record.beginPos + len) <= beginPos)
            continue;

        // calc. depth
        int32_t startIdx = record.beginPos - beginPos;
        int32_t endIdx = startIdx + len;

        if (startIdx < 0)
            startIdx = 0;

        if (endIdx > depthSize)
            endIdx = depthSize;

        for(int i=startIdx; i < endIdx; ++i)
            ++depth[i];
    }
}

bool AlignmentManager::load(void)
{
    if( optionManager == NULL )
        return false;

    if ( !open(this->bamFileIn, toCString(optionManager->getInputFile())) )
    {
        std::cerr << "ERROR: Could not open " << optionManager->getInputFile() << std::endl;
        return false;
    }

    CharString baiFileName = optionManager->getInputFile();
    baiFileName += ".bai";
    if (!open(this->baiIndex, toCString(baiFileName)))
    {
        std::cerr << "ERROR: Could not read BAI index file " << baiFileName << "\n";
        return false;
    }
    this->pBamFileOut = new BamFileOut(context(this->bamFileIn), std::cerr, Sam());

    this->totalRecordNum = 0;
    try
    {
        readHeader(this->bamHeader, this->bamFileIn);

        // read headers
        for (unsigned i=0; i < length(this->bamHeader); ++i)
        {
            if (this->bamHeader[i].type == BamHeaderRecordType::BAM_HEADER_FIRST)
            {
                CharString tempStr;
                getTagValue(tempStr, "SO", this->bamHeader[i]);
                if (tempStr != "coordinate")
                {
                    printMessage("[ERROR] : BAM files must be sorted by coordinate.");
                    exit(1);
                }
            }
            break;
        }

        this->splitRead->prepAfterHeaderParsing(this->bamHeader, this->bamFileIn);
        if ( optionManager->doPairedEndAnalysis() )
            this->pairedEndRead->prepAfterHeaderParsing(this->bamHeader, this->bamFileIn);
        if ( optionManager->doClippedReadAnalysis() )
            this->clippedRead->prepAfterHeaderParsing(this->bamHeader, this->bamFileIn);
        if ( optionManager->doReadDepthAnalysis() )
            this->readDepth->prepAfterHeaderParsing(this->bamHeader, this->bamFileIn);

        BamAlignmentRecord record;
        CharString qNameWithPairInfo;
        int32_t minMapQual = optionManager->getMinMapQual();
        int32_t minSVSize = optionManager->getMinSVSize();
        int32_t minClipSeqSize = optionManager->getMinClippedSeqSize();
        bool checkClippedSequence, checkSplitRead, checkPairedEndRead;        
        
        while (!atEnd(this->bamFileIn))
        {
            readRecord(record, this->bamFileIn);
            ++this->totalRecordNum;

            if (this->totalRecordNum % AlignmentManager::PRINT_READ_NUMBER_PER == 0)
            {
                printTimeMessage(std::to_string(this->totalRecordNum) + " records were parsed.");
                break;
            }
            
            // discards low quality reads, secondary mappings
            //if (hasFlagUnmapped(record) || record.mapQ < minMapQual || hasFlagQCNoPass(record) || hasFlagSecondary(record))
            if (hasFlagUnmapped(record) || hasFlagQCNoPass(record) || hasFlagSecondary(record))
                continue;

            // calculate depth
            if ( optionManager->doReadDepthAnalysis() )
                this->readDepth->parseReadRecord(record.qName, record);

            checkClippedSequence = false;
            checkSplitRead = false;
            checkPairedEndRead = false;

            // primary alignment
            if ( hasFlagSupplementary(record) == false)
            {
                BamTagsDict tagsDict(record.tags);
                int32_t tagIdx;

                // split-read
                if (findTagKey(tagIdx, tagsDict, "SA"))
                {
                    //if (record.mapQ >= minMapQual) // filtering by mapping quality
                        checkSplitRead = true;
                }
                else
                {
                    // clipped-reads, filtering by sequence size
                    int lastCigar = length(record.cigar) - 1;               
                    if ( ((record.cigar[0].operation == 'H'  || record.cigar[0].operation == 'S') && \
                           record.cigar[0].count >= minClipSeqSize) || \
                         ((record.cigar[lastCigar].operation == 'H' || record.cigar[lastCigar].operation == 'S') && \
                           record.cigar[lastCigar].count >= minClipSeqSize) )
                    {
                        checkClippedSequence = true;
                    }
                    else // indels, filtering by sequence size
                    {
                        for (unsigned i=1; i < lastCigar; ++i)
                        {
                            if ((record.cigar[i].count >= minSVSize) && \
                                (record.cigar[i].operation == 'P'  || record.cigar[i].operation == 'I' || \
                                 record.cigar[i].operation == 'D'  || record.cigar[i].operation == 'N'))
                            {
                                checkSplitRead = true;
                            }
                        }
                    }

                    // paired-end reads, filterirng my mapping quaility
                    if (record.mapQ >= minMapQual)
                    {
                        // calculate insertion size
                        if (hasFlagAllProper(record) == true) 
                        {
                            // didn't calculate the median yet
                            if (this->maxAbInsSize == std::numeric_limits<double>::max()) 
                            {                            
                                // save the insertion length for median calculation
                                this->readsForMedInsEst.push_back(std::make_pair(abs(record.tLen), record));

                                // calculate the median insertion size
                                if(this->readsForMedInsEst.size() == AlignmentManager::INS_SIZE_ESTIMATION_SAMPLE_SIZE)
                                {
                                    this->calcInsertSize();
                                    this->splitRead->setInsertionInfo(this->getInsMedian(), this->getInsSD());
                                    this->pairedEndRead->setInsertionInfo(this->getInsMedian(), this->getInsSD());
                                    this->clippedRead->setInsertionInfo(this->getInsMedian(), this->getInsSD());
                                }
                            }
                        }

                        // abnormal pairs
                        if (hasFlagNextUnmapped(record) == false) // both reads have to be mapped
                        {
                            // size abnormality
                            if (this->isAbnormalInsertion(abs(record.tLen)))
                                checkPairedEndRead = true;
                            else if (checkPairedEndRead == false) // orientation abnormality
                            {
                                if (hasFlagRC(record) == hasFlagNextRC(record))
                                    checkPairedEndRead = true; // inverted
                                else if ( (record.tLen > 0 && hasFlagRC(record)) || \
                                          (record.tLen < 0 && hasFlagRC(record) == false) )
                                    checkPairedEndRead = true; // swapped (<-- -->)
                            }
                        }
                    }
                }
            }
            else // suppplementray alignment("chimeric" reads)
                checkSplitRead = true;

            // do not consider this record
            if( checkClippedSequence == false && checkSplitRead == false && checkPairedEndRead == false)
                continue;

            // abnormal pair
            qNameWithPairInfo = record.qName;
            if ( optionManager->doPairedEndAnalysis() && checkPairedEndRead == true )
            {
                this->pairedEndRead->parseReadRecord(qNameWithPairInfo, record);
                ++pairedReadCount;
            }

            if ( hasFlagFirst(record) )
                qNameWithPairInfo += "/1";
            if ( hasFlagLast(record) )
                qNameWithPairInfo += "/2";

            // contains clipped sequences
            if ( optionManager->doClippedReadAnalysis() && checkClippedSequence == true)
            {
                this->clippedRead->parseReadRecord(qNameWithPairInfo, record);
                ++clippedReadCount;
            }

            // contains split-reads
            if ( checkSplitRead == true )
            {
                this->splitRead->parseReadRecord(qNameWithPairInfo, record);
                ++splitReadCount;
            }
        }

        // dataset contains small number of reads
        if (this->maxAbInsSize < 0 && this->readsForMedInsEst.size() != 0)
            this->calcInsertSize();

        // check records used in median estimation
        if (optionManager->doPairedEndAnalysis())
        {
            for (unsigned i=0; i < this->readsForMedInsEst.size(); ++i)
            {
                BamAlignmentRecord& record = this->readsForMedInsEst[i].second;
                
                // just check insertion size. there is no chance to get orientation abnormality in these reads.
                if (this->isAbnormalInsertion(abs(record.tLen)))
                {
                    qNameWithPairInfo = record.qName;
                    this->pairedEndRead->parseReadRecord(qNameWithPairInfo, record);
                    ++pairedReadCount;
                }
            }
            this->readsForMedInsEst.clear();
        }
    }
    catch (Exception const & e)
    {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return false;
    }

    printTimeMessage(std::to_string(this->totalRecordNum) + " records were parsed.");
    return true;
}

void AlignmentManager::calcInsertSize(void)
{
    // median
    sort(this->readsForMedInsEst.begin(), this->readsForMedInsEst.end(), AlignmentManager::pairCompare);
    this->insertMedian = MID_ELEMENT(readsForMedInsEst).first;

    // median absolute deviation
    for (auto it=this->readsForMedInsEst.begin(); it != this->readsForMedInsEst.end(); ++it)
        it->first = abs(it->first - this->insertMedian);
    sort(this->readsForMedInsEst.begin(), this->readsForMedInsEst.end(), AlignmentManager::pairCompare);

    // approximated standard deviation
    this->insertDev = (double) (MID_ELEMENT(readsForMedInsEst).first) * AlignmentManager::K;
    this->minAbInsSize = this->insertMedian - (this->insertDev) * optionManager->getAbInsParam();
    this->maxAbInsSize = this->insertMedian + (this->insertDev) * optionManager->getAbInsParam(); 

    printTimeMessage("Estimated insertion size (median,s.d.) : " + \
                      std::to_string((int)this->insertMedian) + \
                      "," + std::to_string((int)this->insertDev));
}

void AlignmentManager::printRecord(BamAlignmentRecord & record) 
{ 
    writeRecord(*this->pBamFileOut, record); 
}

CharString AlignmentManager::getRefName(int32_t id)
{
    return contigNames(context(this->bamFileIn))[id];
}

int32_t AlignmentManager::getRefCount()
{
    return length(contigNames(context(this->bamFileIn)));
}