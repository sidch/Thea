//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2009, Siddhartha Chaudhuri/Stanford University
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holders nor the names of contributors
// to this software may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
//============================================================================

#include "Serializable.hpp"

#ifdef THEA_EXTERN_TEMPLATES

THEA_INSTANTIATE_SMART_POINTERS(Thea::Serializable)
THEA_INSTANTIATE_SMART_POINTERS(Thea::SerializableFactory)

#endif // THEA_EXTERN_TEMPLATES

namespace Thea {

TextInputStream::Settings
Serializable::initConfigReadSettings()
{
  TextInputStream::Settings settings;
  settings.caseSensitive = false;
  settings.cppBlockComments = true;
  settings.cppLineComments = true;
  settings.escapeSequencesInStrings = true;
  settings.falseSymbols.insert("false");
  settings.falseSymbols.insert("off");
  settings.falseSymbols.insert("0");
  settings.generateCommentTokens = false;
  settings.generateNewlineTokens = false;
  settings.otherCommentCharacter = '#';
  settings.otherCommentCharacter2 = '\'';
  settings.otherLineComments = true;
  settings.proofSymbols = false;
  settings.signedNumbers = true;
  settings.simpleFloatSpecials = true;
  settings.singleQuotedStrings = true;
  settings.falseSymbols.insert("true");
  settings.falseSymbols.insert("on");
  settings.falseSymbols.insert("1");

  return settings;
}

TextOutputStream::Settings
Serializable::initConfigWriteSettings()
{
  TextOutputStream::Settings settings;
  settings.allowWordWrapInsideDoubleQuotes = true;
  settings.convertNewlines = true;
  settings.falseSymbol = "false";

#ifdef THEA_WINDOWS
  settings.newlineStyle = TextOutputStream::Settings::NEWLINE_WINDOWS;
#else
  settings.newlineStyle = TextOutputStream::Settings::NEWLINE_UNIX;
#endif

  settings.numColumns = 8;
  settings.spacesPerIndent = 4;
  settings.trueSymbol = "true";
  settings.wordWrap = TextOutputStream::Settings::WRAP_NONE;

  return settings;
}

} // namespace Thea
