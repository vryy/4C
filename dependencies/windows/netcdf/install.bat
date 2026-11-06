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
set "VERSION=4.10.1"
set "CHECKSUM=33c27231c478c3b35da7c7758fbdd02da1fe407abcb16ddfe195f69d164f930d"
set "ARCHIVE=v%VERSION%.tar.gz"

rem Download netcdf
curl -s -L -o "%ARCHIVE%" "https://github.com/Unidata/netcdf-c/archive/refs/tags/%ARCHIVE%"
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

set "NETCDF_HOME=%~dp0netcdf-c-%VERSION%"

rem Do not use the Zlib find module of netcdf
if exist %NETCDF_HOME%\cmake\modules\FindZLIB.cmake del %NETCDF_HOME%\cmake\modules\FindZLIB.cmake

mkdir "netcdf-c-%VERSION%-build"

cd "netcdf-c-%VERSION%-build"

cmake ^
  -G "Ninja" ^
  -D CMAKE_BUILD_TYPE:STRING="Release" ^
  -D CMAKE_INSTALL_PREFIX="%INSTALL_DIR%" ^
  -D CMAKE_C_COMPILER="cl.exe" ^
  -D CMAKE_CXX_COMPILER="cl.exe" ^
  -D CMAKE_CXX_FLAGS:STRING="%CMAKE_CXX_FLAGS% -D_WIN32 /EHsc /MP" ^
  -D CMAKE_VERBOSE_MAKEFILE:BOOL=OFF ^
  -D CMAKE_COLOR_MAKEFILE:BOOL=ON ^
  -D BUILD_SHARED_LIBS:BOOL=OFF ^
  ^
  -D NETCDF_ENABLE_HDF5=ON ^
  -D HDF5_DIR="%USERPROFILE%/opt/hdf5/cmake" ^
  -D HDF5_USE_STATIC_LIBRARIES=ON ^
  -D ZLIB_ROOT="%USERPROFILE%/opt/zlib" ^
  -D ZLIB_USE_STATIC_LIBS=ON ^
  -D NETCDF_ENABLE_DAP=OFF ^
  -D NETCDF_ENABLE_DAP2=OFF ^
  -D NETCDF_ENABLE_DAP4=OFF ^
  -D NETCDF_BUILD_UTILITIES=OFF ^
  %NETCDF_HOME%

ninja install -j%NPROCS%

rem

cd ..

rem Clean up downloaded and extracted artifacts
del /f /q netcdf*.tar.gz 2>nul
for /d %%D in (netcdf*) do rmdir /s /q "%%D"

if !ERRORLEVEL! neq 0 (
    echo ERROR: install.bat failed with exit code !ERRORLEVEL!
    exit /b !ERRORLEVEL!
)

endlocal
