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
#ifndef APP_MISC_H_
#define APP_MISC_H_

#include <vector>
#include <string>
#include <ctime>

#define FIRST_ELEMENT(x)                    x[0]
#define LAST_ELEMENT(x)                     x[x.size()-1]
#define MID_ELEMENT(x)                      x[x.size()/2]
#define GET_MIN_MAX_ELEMENT(min,max,list)   min=FIRST_ELEMENT(list);max=LAST_ELEMENT(list)
#define RUN(x, y, z)                        {startTimeMessage(y);(x=z);endTimeMessage(y);}
#define RUN_IF(c, x, y, z)                  {if (c) { startTimeMessage(y);(x=z);endTimeMessage(y); } else { x = true; }}
#define CharStringToStdString(x)            std::string(seqan::toCString(x))

// ==========================================================================
// Functions
// ==========================================================================
void printMessage(std::string);
void printTimeMessage(std::string);
void startTimeMessage(std::string);
void endTimeMessage(std::string);
void splitString(std::vector<std::string>&, std::string&, std::string&);
#endif // APP_MISC_H_
