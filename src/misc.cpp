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
#include <iostream>
#include "misc.hpp"

// ============================================================================
// Functions
// ============================================================================
time_t _t1, _t2;
void printMessage(std::string msg)
{
    std::cerr << msg << std::endl;
}

void printTimeMessage(std::string msg)
{
    time_t t = time(0);
    struct tm* now = localtime(&t);
    printMessage(
        "[" + std::to_string(now->tm_year + 1900) + '-' \
         + std::to_string(now->tm_mon + 1) + '-' \
         + std::to_string(now->tm_mday) + " " \
         + std::to_string(now->tm_hour) + ":" \
         + std::to_string(now->tm_min) + ":" \
         + std::to_string(now->tm_sec) + "] " \
         + msg );
}

void startTimeMessage(std::string msg)
{
    time(&_t1);
    printTimeMessage("[START] " + msg);
}

void endTimeMessage(std::string msg)
{
    time(&_t2);
    printTimeMessage("[END] " + msg + " (" + std::to_string((int)difftime(_t2,_t1)) + " seconds.)");
}


void splitString(std::vector<std::string>& v, std::string& s, std::string& d)
{
    // split strings
    size_t last = 0, next = 0;
    v.clear();
    while ((next = s.find(d, last)) != std::string::npos) 
    { 
        v.push_back(s.substr(last, next-last));
        last = next + d.length();
    }
    v.push_back(s.substr(last));
}
