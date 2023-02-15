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

#include "FileSystem.hpp"
#include "StringAlg.hpp"
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <filesystem>

namespace Thea {

bool
FileSystem::exists(std::string const & path)
{
  std::error_code err;
  return std::filesystem::exists(path, err);
}

bool
FileSystem::fileExists(std::string const & path)
{
  std::error_code err;
  return std::filesystem::is_regular_file(path, err);
}

bool
FileSystem::directoryExists(std::string const & path)
{
  std::error_code err;
  return std::filesystem::is_directory(path, err);
}

int64
FileSystem::fileSize(std::string const & path)
{
  std::error_code err;
  std::uintmax_t size = std::filesystem::file_size(path, err);
  if (size == static_cast<std::uintmax_t>(-1))
    return -1;
  else
    return static_cast<int64>(size);
}

std::string
FileSystem::resolve(std::string const & path)
{
  std::error_code err;
  auto a = std::filesystem::absolute(path, err);
  return err ? std::string() : a.string();
}

bool
FileSystem::createDirectory(std::string const & path)
{
  if (directoryExists(path))
    return true;

  std::error_code err;
  return std::filesystem::create_directories(path, err);
}

bool
FileSystem::readWholeFile(std::string const & path, std::string & ret)
{
  if (!fileExists(path))
  {
    THEA_ERROR << "FileSystem: File '" << path << "' not found";
    return false;
  }

  int64 length = fileSize(path);
  if (length <= 0)
  {
    ret.clear();
    return true;
  }

  char * buffer = (char *)std::malloc((size_t)length);
  if (!buffer)
  {
    THEA_ERROR << "FileSystem: Could not allocate buffer to hold " << length << " bytes from file '" << path << '\'';
    return false;
  }

  FILE * f = std::fopen(path.c_str(), "rb");
  if (!f)
  {
    THEA_ERROR << "FileSystem: Couldn't open file '" << path << "' for reading";
    return false;
  }

  size_t num_read = std::fread(buffer, 1, length, f);
  if ((int64)num_read != length)
  {
    THEA_ERROR << "FileSystem: Error reading from file '" << path << '\'';
    return false;
  }

  std::fclose(f);

  ret.assign(buffer, (size_t)length);
  std::free(buffer);

  return true;
}

namespace FileSystemInternal {

bool
objectSatisfiesConstraints(std::filesystem::directory_entry const & object, int types,
                           Array<std::string> const & patterns, bool ignore_case)
{
  if (types > 0 && types != FileSystem::ObjectType::ALL)
  {
    std::filesystem::file_status status = object.symlink_status();
    if (!std::filesystem::is_symlink(status))
      status = object.status();

    bool ok = false;

    if (!ok && (types & FileSystem::ObjectType::FILE) && std::filesystem::is_regular_file(status))
      ok = true;

    if (!ok && (types & FileSystem::ObjectType::DIRECTORY) && std::filesystem::is_directory(status))
      ok = true;

    if (!ok && (types & FileSystem::ObjectType::SYMLINK) && std::filesystem::is_symlink(status))
      ok = true;

    if (!ok)
      return false;
  }

  if (!patterns.empty())
  {
    std::string name = object.path().filename().string();
    if (ignore_case)
      name = toLower(name);

    bool ok = false;
    for (size_t i = 0; !ok && i < patterns.size(); ++i)
      if (patternMatch(patterns[i], name))
        ok = true;

    if (!ok)
      return false;
  }

  return true;
}

} // namespace FileSystemInternal

intx
FileSystem::getDirectoryContents(std::string const & dir, Array<std::string> & objects, int types, std::string const & patterns,
                                 int flags)
{
  if (!directoryExists(dir))
    return -1;

  bool ignore_case  =  (flags & Flags::CASE_INSENSITIVE);
  bool recursive    =  (flags & Flags::RECURSIVE);
  bool sorted       =  (flags & Flags::SORTED);

  Array<std::string> patlist;
  if (!patterns.empty())
  {
    stringSplit(patterns, ' ', patlist, true);

    if (ignore_case)
    {
      for (size_t i = 0; i < patlist.size(); ++i)
        patlist[i] = toLower(patlist[i]);
    }
  }

  objects.clear();

  if (recursive)
  {
    std::filesystem::recursive_directory_iterator objects_end;
    for (std::filesystem::recursive_directory_iterator iter(dir); iter != objects_end; ++iter)
      if (FileSystemInternal::objectSatisfiesConstraints(*iter, types, patlist, ignore_case))
        objects.push_back(iter->path().string());
  }
  else
  {
    std::filesystem::directory_iterator objects_end;
    for (std::filesystem::directory_iterator iter(dir); iter != objects_end; ++iter)
      if (FileSystemInternal::objectSatisfiesConstraints(*iter, types, patlist, ignore_case))
        objects.push_back(iter->path().string());
  }

  if (sorted)
    std::sort(objects.begin(), objects.end());

  return (intx)objects.size();
}

bool
FileSystem::remove(std::string const & path, bool recursive)
{
  std::error_code err;
  if (recursive)
    return std::filesystem::remove_all(path, err) != static_cast<std::uintmax_t>(-1);
  else
    return std::filesystem::remove(path, err);
}

bool
FileSystem::copyFile(std::string const & from, std::string const & to)
{
  std::error_code err;
  return std::filesystem::copy_file(from, to, err);
}

} // namespace Thea
