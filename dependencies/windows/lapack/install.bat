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
set "VERSION=3.12.1"
set "CHECKSUM=2ca6407a001a474d4d4d35f3a61550156050c48016d949f0da0529c0aa052422"
set "ARCHIVE=v%VERSION%.tar.gz"

rem Download lapack
curl -s -L -o "%ARCHIVE%" "https://github.com/Reference-LAPACK/lapack/archive/refs/tags/%ARCHIVE%"
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

rem compiling

mkdir "lapack-%VERSION%-build"

cd "lapack-%VERSION%-build"

set "LAPACK_HOME=%~dp0lapack-%VERSION%"

cmake ^
  -G "Ninja" ^
  -D CMAKE_BUILD_TYPE:STRING="Release" ^
  -D CMAKE_INSTALL_PREFIX="%INSTALL_DIR%" ^
  -D CMAKE_C_COMPILER="cl.exe" ^
  -D CMAKE_CXX_COMPILER="cl.exe" ^
  -D CMAKE_Fortran_COMPILER="ifx.exe" ^
  -D CMAKE_Fortran_FLAGS="/names:lowercase /assume:underscore" ^
  -D CMAKE_VERBOSE_MAKEFILE:BOOL=OFF ^
  -D CMAKE_COLOR_MAKEFILE:BOOL=ON ^
  -D BUILD_SHARED_LIBS:BOOL=OFF ^
  -D CBLAS:BOOL=OFF ^
  -D LAPACKE:BOOL=OFF ^
  -D BUILD_SINGLE:BOOL=ON ^
  -D BUILD_DOUBLE:BOOL=ON ^
  -D BUILD_COMPLEX:BOOL=ON ^
  -D BUILD_COMPLEX16:BOOL=ON ^
  -D BUILD_DEPRECATED:BOOL=ON ^
  %LAPACK_HOME%

ninja install

rem

cd ..

rem Clean up downloaded and extracted artifacts
del /f /q lapack*.tar.gz 2>nul
for /d %%D in (lapack*) do rmdir /s /q "%%D"

endlocal
