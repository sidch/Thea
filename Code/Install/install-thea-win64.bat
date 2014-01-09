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
rem cmake "-DCMAKE_INSTALL_PREFIX=%THEA_PREFIX%" -DBOOST_ROOT=C:\Boost32 -DQT_QMAKE_EXECUTABLE=C:\Qt32\4.8.2\bin\qmake.exe -DCGAL_ROOT=C:\CGAL-3.7-32 -G "Visual Studio 10 Win64" .
cmake "-DCMAKE_INSTALL_PREFIX=%THEA_PREFIX%" -G "Visual Studio 10 Win64" .

devenv /Build Release Thea.sln
devenv /Build Release /Project INSTALL Thea.sln

devenv /Build Debug Thea.sln
devenv /Build Debug /Project INSTALL Thea.sln

set THEA_PREFIX=

popd
