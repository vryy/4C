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
set "VERSION=1.2.13"
set "CHECKSUM=1525952a0a567581792613a9723333d7f8cc20b87a81f920fb8bc7e3f2251428"
set "ARCHIVE=v%VERSION%.tar.gz"

rem Download zlib
curl -s -L -o "%ARCHIVE%" "https://github.com/madler/zlib/archive/refs/tags/%ARCHIVE%"
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

set "ZLIB_HOME=%~dp0zlib-%VERSION%

mkdir "zlib-%VERSION%-build"

cd "zlib-%VERSION%-build"

cmake ^
  -G "Ninja" ^
  -D CMAKE_BUILD_TYPE:STRING="Release" ^
  -D CMAKE_INSTALL_PREFIX="%INSTALL_DIR%" ^
  -D CMAKE_POLICY_VERSION_MINIMUM="3.5" ^
  -D CMAKE_C_COMPILER="cl.exe" ^
  -D BUILD_SHARED_LIBS:BOOL=OFF ^
  %ZLIB_HOME%

ninja install -j%NPROCS%

rem

cd ..

rem Clean up downloaded and extracted artifacts
del /f /q zlib*.tar.gz 2>nul
for /d %%D in (zlib*) do rmdir /s /q "%%D"

if !ERRORLEVEL! neq 0 (
    echo ERROR: install.bat failed with exit code !ERRORLEVEL!
    exit /b !ERRORLEVEL!
)

endlocal
