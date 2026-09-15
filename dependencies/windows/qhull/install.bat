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
set "VERSION=8.1-alpha6"
set "CHECKSUM=033a155cfed37f811f881d4db16668157027e3cb6c9bf2ac5a036cd44c731dc1"
set "ARCHIVE=v%VERSION%.tar.gz"

rem Download qhull
curl -s -L -o "%ARCHIVE%" "https://github.com/qhull/qhull/archive/refs/tags/%ARCHIVE%"
if !ERRORLEVEL! neq 0 (
    echo Failed to download %ARCHIVE%
    exit /b !ERRORLEVEL!
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
if !ERRORLEVEL! neq 0 (
    echo Failed to extract archive
    exit /b !ERRORLEVEL!
)

rem compiling

set "QHULL_HOME=%~dp0qhull-%VERSION%"

mkdir "qhull-%VERSION%-build"

cd "qhull-%VERSION%-build"

cmake ^
  -G "Ninja" ^
  -D CMAKE_BUILD_TYPE:STRING="Release" ^
  -D CMAKE_INSTALL_PREFIX="%INSTALL_DIR%" ^
  -D CMAKE_C_COMPILER="cl.exe" ^
  -D CMAKE_CXX_COMPILER="cl.exe" ^
  -D BUILD_SHARED_LIBS:BOOL=OFF ^
  -D BUILD_STATIC_LIBS:BOOL=ON ^
  -D WITH_LFS:BOOL=ON ^
  %QHULL_HOME%

if !ERRORLEVEL! neq 0 (
    echo ERROR: installing qhull failed at configuration
    exit /b !ERRORLEVEL!
)

ninja install -j%NPROCS%

if !ERRORLEVEL! neq 0 (
    echo ERROR: installing qhull failed at compilation
    exit /b !ERRORLEVEL!
)

rem

cd ..

rem Clean up downloaded and extracted artifacts
del /f /q *.tar.gz 2>nul
for /d %%D in (qhull*) do rmdir /s /q "%%D"

endlocal
