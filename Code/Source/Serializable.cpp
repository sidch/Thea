//============================================================================
//
// This file is part of the Thea toolkit.
//
// This software is distributed under the BSD license, as detailed in the
// accompanying LICENSE.txt file. Portions are derived from other works:
// their respective licenses and copyright information are reproduced in
// LICENSE.txt and/or in the relevant source files.
//
// Author: Siddhartha Chaudhuri
// First version: 2009
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
  settings.newlineStyle = TextOutputStream::NewlineStyle::WINDOWS;
#else
  settings.newlineStyle = TextOutputStream::NewlineStyle::UNIX;
#endif

  settings.numColumns = 8;
  settings.spacesPerIndent = 4;
  settings.trueSymbol = "true";
  settings.wordWrap = TextOutputStream::WordWrap::NONE;

  return settings;
}

} // namespace Thea
