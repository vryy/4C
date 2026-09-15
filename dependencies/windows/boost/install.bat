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
set "VERSION=1_92_0"
set "CHECKSUM=c4a3b310ddd2472416e091067166b0713be97c63f38c212c484ada022fd296ce"
set "ARCHIVE=boost_%VERSION%.tar.gz"

rem Download boost
curl -s -L -o "%ARCHIVE%" "https://archives.boost.io/release/1.92.0/source/%ARCHIVE%"
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

cd "boost_%VERSION%"

set "CC=cl.exe"
set "CXX=cl.exe"
call .\bootstrap.bat

.\b2.exe install -j%NPROCS% --prefix="%INSTALL_DIR%" toolset=msvc variant=release address-model=64 threading=multi

rem

cd ..

rem Clean up downloaded and extracted artifacts
del /f /q boost*.tar.gz 2>nul
for /d %%D in (boost*) do rmdir /s /q "%%D"

if !ERRORLEVEL! neq 0 (
    echo ERROR: install.bat failed with exit code !ERRORLEVEL!
    exit /b !ERRORLEVEL!
)

endlocal
