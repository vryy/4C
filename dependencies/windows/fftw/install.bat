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
set "VERSION=3.3.11"
set "CHECKSUM=5630c24cdeb33b131612f7eb4b1a9934234754f9f388ff8617458d0be6f239a1"
set "ARCHIVE=fftw-%VERSION%.tar.gz"

rem Download fftw
curl -s -L -o "%ARCHIVE%" "https://www.fftw.org/%ARCHIVE%"
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

set "FFTW_HOME=%~dp0fftw-%VERSION%"

mkdir "fftw-%VERSION%-build"

cd "fftw-%VERSION%-build"

cmake ^
  -G "Ninja" ^
  -D CMAKE_POLICY_VERSION_MINIMUM="3.5" ^
  -D CMAKE_BUILD_TYPE:STRING="Release" ^
  -D CMAKE_INSTALL_PREFIX="%INSTALL_DIR%" ^
  -D CMAKE_C_COMPILER="cl.exe" ^
  -D CMAKE_CXX_COMPILER="cl.exe" ^
  -D BUILD_SHARED_LIBS:BOOL=OFF ^
  -D BUILD_TESTS:BOOL=OFF ^
  -D ENABLE_AVX:BOOL=ON ^
  -D ENABLE_AVX2:BOOL=ON ^
  -D ENABLE_SSE:BOOL=ON ^
  -D ENABLE_SSE2:BOOL=ON ^
  -D ENABLE_OPENMP:BOOL=OFF ^
  %FFTW_HOME%

if !ERRORLEVEL! neq 0 (
    echo ERROR: installing fftw failed at configuration
    exit /b !ERRORLEVEL!
)

ninja install -j%NPROCS%

if !ERRORLEVEL! neq 0 (
    echo ERROR: installing fftw failed at compilation
    exit /b !ERRORLEVEL!
)

rem

cd ..

rem Clean up downloaded and extracted artifacts
del /f /q *.tar.gz 2>nul
for /d %%D in (fftw*) do rmdir /s /q "%%D"

endlocal
