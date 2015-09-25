//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2013, Siddhartha Chaudhuri/Princeton University
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

#include "FilePath.hpp"
#include <boost/filesystem.hpp>

namespace Thea {

std::string
FilePath::baseName(std::string const & path)
{
  std::string node = objectName(path);
  size_t first_dot = node.find_first_of('.');
  return (first_dot == std::string::npos ? node : node.substr(0, first_dot));
}

std::string
FilePath::completeBaseName(std::string const & path)
{
  std::string node = objectName(path);
  size_t last_dot = node.find_last_of('.');
  return (last_dot == std::string::npos ? node : node.substr(0, last_dot));
}

std::string
FilePath::extension(std::string const & path)
{
  std::string node = objectName(path);
  size_t last_dot = node.find_last_of('.');
  return (last_dot == std::string::npos ? std::string() : node.substr(last_dot + 1));
}

std::string
FilePath::completeExtension(std::string const & path)
{
  std::string node = objectName(path);
  size_t first_dot = node.find_first_of('.');
  return (first_dot == std::string::npos ? std::string() : node.substr(first_dot + 1));
}

std::string
FilePath::objectName(std::string const & path)
{
  boost::filesystem::path p(path);
  boost::filesystem::path n = p.filename();

  while (n.string() == "." && p.string() != ".")
  {
    p.remove_filename();
    n = p.filename();
  }

  return n.string();
}

std::string
FilePath::parent(std::string const & path)
{
  boost::filesystem::path p(path);
  boost::filesystem::path n = p.filename();

  while (n.string() == "." && p.string() != ".")
  {
    p.remove_filename();
    n = p.filename();
  }

  return p.parent_path().string();
}

std::string
FilePath::concat(std::string const & parent_name, std::string const & child_name)
{
  boost::filesystem::path p(parent_name);
  p /= child_name;
  return p.string();
}

std::string
FilePath::changeExtension(std::string const & path, std::string const & new_ext)
{
  return concat(parent(path), completeBaseName(path) + '.' + new_ext);
}

std::string
FilePath::changeCompleteExtension(std::string const & path, std::string const & new_ext)
{
  return concat(parent(path), baseName(path) + '.' + new_ext);
}

} // namespace Thea
