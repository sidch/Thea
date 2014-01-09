@echo off

echo !!! Please ensure all support libraries have already been installed !!!
echo.
echo To specify a installation directory for Thea other than the directory
echo containing this script, pass this directory as the first argument to
echo the script. This should also be the directory containing the support
echo libraries.
echo.
echo Installing Thea...

pushd %~dp0

if "%1"=="" (
  set THEA_PREFIX="%CD%"
) else (
  set THEA_PREFIX="%1"
)

cd ..\Build
cmake "-DCMAKE_INSTALL_PREFIX=%THEA_PREFIX%" -G "Visual Studio 10 Win64" .

msbuild Thea.sln /p:Configuration="Release"
msbuild INSTALL.vcxproj /p:Configuration="Release"

msbuild Thea.sln /p:Configuration="Debug"
msbuild INSTALL.vcxproj /p:Configuration="Debug"

set THEA_PREFIX=

popd
