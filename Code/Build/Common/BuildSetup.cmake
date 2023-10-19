#===============================================================================================================================
#
# Common build setup instructions for the Thea library.
#
# Copyright (C) 2019, Siddhartha Chaudhuri
#
# Before including this file, set the following variables:
#   - Thea_PROJECT_ROOT: path to root directory of project
#   - Thea_FIND_<packagename>: set to true for each package to be found as a dependency (or set Thea_FIND_ALL to true)
#
# Remember to link the target with Thea_DEPS_LIBRARIES and Thea_DEPS_LDFLAGS in the calling script.
#
#===============================================================================================================================

# Set the minimum required CMake version
CMAKE_MINIMUM_REQUIRED(VERSION 3.0...3.27.7)

# See cmake --help-policy CMP0003 for details
IF(POLICY CMP0003)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(POLICY CMP0003)

# See cmake --help-policy CMP0042 for details
IF(POLICY CMP0042)
  CMAKE_POLICY(SET CMP0042 NEW)
ENDIF(POLICY CMP0042)

# See cmake --help-policy CMP0074 for details
IF(POLICY CMP0074)
  CMAKE_POLICY(SET CMP0074 NEW)
ENDIF(POLICY CMP0074)

# If you don't want the full compiler output, remove the following line
SET(CMAKE_VERBOSE_MAKEFILE ON)

# Avoid having to repeat condition after ELSE and ENDIF statements
SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE)

# Postfix for debug builds
SET(CMAKE_DEBUG_POSTFIX "d")

# Convert ProjectRoot to absolute path, if not so already
GET_FILENAME_COMPONENT(ProjectRoot ${Thea_PROJECT_ROOT} ABSOLUTE)

# Path for build products
SET(OutputRoot ${ProjectRoot}/Build/Output)

# Path to put executables in
SET(EXECUTABLE_OUTPUT_PATH ${OutputRoot}/bin)

# Path to put libraries in
SET(LIBRARY_OUTPUT_PATH ${OutputRoot}/lib)

# Path for customized CMake modules
IF(NOT CMAKE_MODULE_PATH)
  SET(CMAKE_MODULE_PATH ${ProjectRoot}/Build/Common/CMake/Modules)
ENDIF()
GET_FILENAME_COMPONENT(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ABSOLUTE)

# Path to root folder for source code
SET(SourceRoot ${ProjectRoot}/Source)

# Path to folder for third-party packages that will be directly built/included as part of the source tree
SET(ThirdPartyRoot ${SourceRoot}/ThirdParty)

# Path to folder with installations of external dependencies
IF(NOT THEA_DEPS_ROOT)
  SET(THEA_DEPS_ROOT ${CMAKE_INSTALL_PREFIX})
ENDIF()
SET(THEA_DEPS_ROOT ${THEA_DEPS_ROOT} CACHE PATH "Path to folder with installations of dependencies")

# Add some additional system search paths which may be omitted by default
SET(CMAKE_SYSTEM_PREFIX_PATH ${CMAKE_SYSTEM_PREFIX_PATH} /usr/local /opt/local)

# Make sure the LITE flag has a non-empty value since it will be cached in BuildFlags
IF(NOT DEFINED LITE)
  SET(LITE FALSE)
ENDIF()

# Locate dependencies
INCLUDE(${ProjectRoot}/Build/Common/FindTheaDependencies.cmake)

# Definitions, compiler switches etc.
STRING(REPLACE ";" " " EXTRA_DEBUG_CFLAGS "${Thea_DEPS_DEBUG_CFLAGS}")
STRING(REPLACE ";" " " EXTRA_RELEASE_CFLAGS "${Thea_DEPS_RELEASE_CFLAGS}")
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${EXTRA_DEBUG_CFLAGS}")
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${EXTRA_RELEASE_CFLAGS}")
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ${EXTRA_RELEASE_CFLAGS}")
SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${EXTRA_DEBUG_CFLAGS}")
SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${EXTRA_RELEASE_CFLAGS}")
SET(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} ${EXTRA_RELEASE_CFLAGS}")

INCLUDE(CheckCXXCompilerFlag)

IF(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  CHECK_CXX_COMPILER_FLAG("-std=c++17" THEA_COMPILER_SUPPORTS_CXX17)
  IF(THEA_COMPILER_SUPPORTS_CXX17)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
  ELSE()
    MESSAGE(FATAL_ERROR "Compiler ${CMAKE_CXX_COMPILER} has no C++17 support")
  ENDIF()

  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Thea_DEPS_CFLAGS} -Wall -fno-strict-aliasing -fPIC")
  SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g2")
  SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DNDEBUG -O2")
  SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DNDEBUG -g2 -O2")

  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${Thea_DEPS_CFLAGS} -Wall -fno-strict-aliasing -fPIC")
  SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g2")
  SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -DNDEBUG -O2")
  SET(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} -DNDEBUG -g2 -O2")

ELSEIF(MSVC)

  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Thea_DEPS_CFLAGS}")
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${Thea_DEPS_CFLAGS}")
  SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2")
  ADD_DEFINITIONS(-D_SCL_SECURE_NO_WARNINGS)

ENDIF()

# Specify that Triangle should be built as a library, use ANSI function declarations and not use a timer
ADD_DEFINITIONS(-DTRILIBRARY -DANSI_DECLARATORS -DNO_TIMER)

# Use extern templates?
IF(NOT DEFINED EXTERN_TEMPLATES)
  SET(EXTERN_TEMPLATES FALSE)
ENDIF()
SET(EXTERN_TEMPLATES ${EXTERN_TEMPLATES} CACHE BOOL "Use extern templates?")

IF(EXTERN_TEMPLATES)
  ADD_DEFINITIONS(-DTHEA_EXTERN_TEMPLATES=1)
ENDIF()

# Include directories
INCLUDE_DIRECTORIES(BEFORE ${Thea_DEPS_INCLUDE_DIRS})

# Link directories
LINK_DIRECTORIES(${Thea_DEPS_LIBRARY_DIRS})
