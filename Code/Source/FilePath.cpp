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
// First version: 2013
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

bool
FilePath::isAbsolute(std::string const & path)
{
  return boost::filesystem::path(path).is_absolute();
}

bool
FilePath::isRelative(std::string const & path)
{
  return boost::filesystem::path(path).is_relative();
}

} // namespace Thea
