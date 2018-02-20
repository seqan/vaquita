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
#include "vaquita.hpp"
#include "mergeoption.hpp"

// ============================================================================
// Functions
// ============================================================================

void MergeOptionManager::init()
{
    setAppName(*this, APP_NAME);
    setShortDescription(*this, "Merging mode");

    // version & date
    setVersion(*this, SEQAN_APP_VERSION);
    setDate(*this, SEQAN_DATE);

    addDescription(*this, std::string(APP_NAME) + std::string(": ") + std::string(APP_TITLE));
    addDescription(*this, std::string(APP_AUTHOR_INFO));
    addDescription(*this, std::string(APP_WEBSITE_INFO));

    // synopsis
    addUsageLine(*this, "merge [\\fIOPTIONS\\fP] [\\fIin_1.vcf\\fP, .. ,\\fIin_N.vcf\\fP] > [\\fIout.vcf\\fP]");

    // options
    addOption(*this, ArgParseOption("", "use-all", "Consider all records. (Default: Consider SVs with 'PASS' in the filter column.)"));
    setDefaultValue(*this, "use-all", "false");

    // files
    ArgParseArgument arg(ArgParseArgument::STRING, "INPUT");
    addArgument(*this, arg);
}

bool MergeOptionManager::parseCommandLine(int argc, char const ** argv)
{
    ArgumentParser::ParseResult res = parse(*this, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return false;

    CharString vcfFileString;    
    getArgumentValue(vcfFileString, *this, 0);

    // split strings
    std::string s = CharStringToStdString(vcfFileString), d = ",";
    splitString(this->vcfFiles, s, d);

    // extention check
    for (unsigned i=0; i < this->vcfFiles.size(); ++i)
    {
        std::string ext = this->vcfFiles[i].substr(this->vcfFiles[i].length()-4);
        if (ext != ".vcf")
        {
            std::cerr << "Unknown file extention: " << ext << " (" << this->vcfFiles[i] << ")\n";
            return false;
        }
    }

    this->useAll = isSet(*this, "use-all");

    return true;
}