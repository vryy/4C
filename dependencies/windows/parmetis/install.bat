@echo off
rem This file is part of 4C multiphysics licensed under the
rem GNU Lesser General Public License v3.0 or later.
rem
rem See the LICENSE.md file in the top-level for license information.
rem
rem SPDX-License-Identifier: LGPL-3.0-or-later

setlocal enabledelayedexpansion

set "INSTALL_DIR=%~1"

rem Number of procs for building (default 4 if not already set)
if "%NPROCS%"=="" set "NPROCS=4"
set "VERSION=4.0.3"
set "CHECKSUM=f2d9a231b7cf97f1fee6e8c9663113ebf6c240d407d3c118c55b3633d6be6e5f"
set "ARCHIVE=parmetis-%VERSION%.tar.gz"

rem Download parmetis
curl -s -L -o "%ARCHIVE%" "https://ftp.mcs.anl.gov/pub/pdetools/spack-pkgs/%ARCHIVE%"
if errorlevel 1 (
    echo Failed to download %ARCHIVE%
    exit /b 1
)

rem Verify checksum using certutil
for /f "skip=1 tokens=* delims=" %%A in ('certutil -hashfile "%ARCHIVE%" SHA256') do (
    if not defined FILE_HASH (
        set "FILE_HASH=%%A"
        set "FILE_HASH=!FILE_HASH: =!"
    )
)

if /i "%FILE_HASH%"=="%CHECKSUM%" (
    echo Checksum matches
) else (
    echo Checksum does not match
    exit /b 1
)

tar -xzf "%ARCHIVE%"
if errorlevel 1 (
    echo Failed to extract archive
    exit /b 1
)

rem apply patch

ren parmetis-%VERSION% parmetis-%VERSION%-orig

cd parmetis-%VERSION%-orig

set "GIT_DIR=none"
git apply --no-index -p2 -v ..\fixes.patch

cd ..

ren parmetis-%VERSION%-orig parmetis-%VERSION%

rem compiling

mkdir "parmetis-%VERSION%-build"

cd "parmetis-%VERSION%-build"

set "PARMETIS_HOME=%~dp0parmetis-%VERSION%"

set "METIS_HOME=%PARMETIS_HOME%\metis"

cmake ^
  -G "Ninja" ^
  -D CMAKE_BUILD_TYPE:STRING="Release" ^
  -D CMAKE_INSTALL_PREFIX="%INSTALL_DIR%" ^
  -D CMAKE_C_COMPILER="cl.exe" ^
  -D CMAKE_CXX_COMPILER="cl.exe" ^
  -D CMAKE_CXX_FLAGS:STRING="%CMAKE_CXX_FLAGS% -D_WIN32 /EHsc /MP" ^
  -D CMAKE_POLICY_VERSION_MINIMUM="3.5" ^
  -D SHARED:BOOL=OFF ^
  -D METIS_INSTALL:BOOL=ON ^
  %METIS_HOME%

ninja install
ninja clean

del CMakeCache.txt

cmake ^
  -G "Ninja" ^
  -D CMAKE_BUILD_TYPE:STRING="Release" ^
  -D CMAKE_INSTALL_PREFIX="%INSTALL_DIR%" ^
  -D CMAKE_C_COMPILER="cl.exe" ^
  -D CMAKE_CXX_COMPILER="cl.exe" ^
  -D CMAKE_CXX_FLAGS:STRING="%CMAKE_CXX_FLAGS% -D_WIN32 /EHsc /MP" ^
  -D CMAKE_POLICY_VERSION_MINIMUM="3.5" ^
  -D SHARED:BOOL=OFF ^
  -D METIS_PATH="%INSTALL_DIR%" ^
  -D MPI_INCLUDE_PATH="C:/Program Files (x86)/Microsoft SDKs/MPI/Include" ^
  -D MPI_LIBRARIES="C:/Program Files (x86)/Microsoft SDKs/MPI/Lib/x64/msmpi.lib" ^
  %PARMETIS_HOME%

ninja install
ninja clean

rem

cd ..

rem Clean up downloaded and extracted artifacts
del /f /q parmetis*.tar.gz 2>nul
for /d %%D in (parmetis*) do rmdir /s /q "%%D"

if %ERRORLEVEL% neq 0 (
    echo ERROR: install.bat failed with exit code %ERRORLEVEL%
    exit /b %ERRORLEVEL%
)

endlocal
