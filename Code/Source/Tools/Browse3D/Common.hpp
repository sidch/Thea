//============================================================================
//
// This file is part of the Browse3D project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2011, Siddhartha Chaudhuri/Stanford University
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

#ifndef __Browse3D_Common_hpp__
#define __Browse3D_Common_hpp__

#define QT_USE_FAST_CONCATENATION
#define QT_USE_FAST_OPERATOR_PLUS

#include "../../Common.hpp"
#include "../../Colors.hpp"
#include "../../FilePath.hpp"
#include "../../IOStream.hpp"
#include <QDebug>
#include <QtDebug>
#include <QString>
#include <string>
#include <iostream>

/** Allow a std::string to be piped to a Qt debug stream. */
inline QDebug
operator<<(QDebug dbg, std::string const & str)
{
  dbg << str.c_str();
  return dbg;
}

/** Allow a QString to be piped to a standard output stream. */
inline std::ostream &
operator<<(std::ostream & out, QString const & s)
{
  return out << s.toAscii().data();
}

/** Convert a std::string to a Qt string. */
inline QString
toQString(std::string const s)
{
  return QString::fromAscii(s.data(), (int)s.size());
}

/** Convert a Qt string to a std::string (use this function and <b>NOT</b> QString::toStdString()!). */
inline std::string
toStdString(QString const s)
{
  return std::string(s.toAscii().data(), (int)s.length());
}

/** Namespace for data-driven texturing project. */
namespace Browse3D {

using namespace Thea;

// Put quotes around the result of a macro expansion.
#define BROWSE3D_STRINGIFY_(x) #x
#define BROWSE3D_STRINGIFY(x) BROWSE3D_STRINGIFY_(x)

/** Construct a fully qualified path for a file, given the name of the file and the path to its parent directory. */
inline std::string getFullPath(std::string const & dir, std::string const & filename)
{
  return FilePath::concat(dir, filename);
}

/** Construct a fully qualified path for a file, given the name of the file and the path to its parent directory. */
QString getFullPath(QString const & dir, QString const & filename);

} // namespace Browse3D

#endif
