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
set "VERSION=2.2.0"
set "CHECKSUM=5b8d75125ae7b4fef55d3b39e8d3dbfd238cdf63b6518afd67fa69ac80f06542"
set "ARCHIVE=%VERSION%.tar.gz"

rem Download hdf5
curl -s -L -o "%ARCHIVE%" "https://github.com/HDFGroup/hdf5/archive/refs/tags/%ARCHIVE%"
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

set "HDF5_HOME=%~dp0hdf5-%VERSION%"

mkdir "hdf5-%VERSION%-build"

cd "hdf5-%VERSION%-build"

cmake ^
  -G "Ninja" ^
  -D CMAKE_BUILD_TYPE:STRING="Release" ^
  -D CMAKE_INSTALL_PREFIX="%INSTALL_DIR%" ^
  -D CMAKE_C_COMPILER="cl.exe" ^
  -D CMAKE_CXX_COMPILER="cl.exe" ^
  -D BUILD_SHARED_LIBS:BOOL=OFF ^
  -D BUILD_STATIC_LIBS:BOOL=ON ^
  -D HDF5_BUILD_CPP_LIB:BOOL=ON ^
  -D HDF5_BUILD_HL_LIB:BOOL=ON ^
  -D HDF5_BUILD_TOOLS:BOOL=ON ^
  -D HDF5_BUILD_UTILS:BOOL=ON ^
  -D HDF5_DEFAULT_API_VERSION="v200" ^
  -D HDF5_ENABLE_ZLIB_SUPPORT:BOOL=ON ^
  -D HDF5_USE_ZLIB_STATIC:BOOL=ON ^
  -D ZLIB_ROOT="%USERPROFILE%\opt\zlib" ^
  %HDF5_HOME%

ninja install -j%NPROCS%

rem

cd ..

rem Clean up downloaded and extracted artifacts
del /f /q *.tar.gz 2>nul
for /d %%D in (hdf5*) do rmdir /s /q "%%D"

if !ERRORLEVEL! neq 0 (
    echo ERROR: install.bat failed with exit code !ERRORLEVEL!
    exit /b !ERRORLEVEL!
)

endlocal
