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
set "VERSION=6.4.0"
set "CHECKSUM=cb9c0b2ba4c28e5ed5817718ba19ae1dd63ccd30bc44c8b8252b54f5f04a44cc"
set "ARCHIVE=v%VERSION%.tar.gz"

rem Download superlu_dist
curl -s -L -o "%ARCHIVE%" "https://github.com/xiaoyeli/superlu_dist/archive/refs/tags/%ARCHIVE%"
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

rem apply patch

ren superlu_dist-%VERSION% superlu_dist-%VERSION%-orig

cd superlu_dist-%VERSION%-orig

set "GIT_DIR=none"
git apply --no-index -p2 -v ..\fixes.patch

cd ..

ren superlu_dist-%VERSION%-orig superlu_dist-%VERSION%

rem compiling

set "SUPERLU_HOME=%~dp0superlu_dist-%VERSION%"

mkdir "superlu_dist-%VERSION%-build"

cd "superlu_dist-%VERSION%-build"

cmake ^
  -G "Ninja" ^
  -D CMAKE_POLICY_VERSION_MINIMUM="3.5" ^
  -D CMAKE_BUILD_TYPE:STRING="Release" ^
  -D CMAKE_INSTALL_PREFIX="%INSTALL_DIR%" ^
  -D CMAKE_C_COMPILER="cl.exe" ^
  -D CMAKE_CXX_COMPILER="cl.exe" ^
  -D CMAKE_C_FLAGS:STRING="%CMAKE_C_FLAGS% /std:c11" ^
  -D CMAKE_CXX_FLAGS:STRING="%CMAKE_CXX_FLAGS%" ^
  -D CMAKE_SHARED_LINKER_FLAGS="/NODEFAULTLIB:libcmt" ^
  -D BUILD_SHARED_LIBS:BOOL=OFF ^
  -D enable_openmp:BOOL=OFF ^
  -D XSDK_ENABLE_Fortran:BOOL=OFF ^
  -D BUILD_TESTING:BOOL=OFF ^
  -D BLAS_LIBRARIES="%USERPROFILE%\opt\lapack\lib\libblas.lib" ^
  -D LAPACK_LIBRARIES="%USERPROFILE%\opt\lapack\lib\liblapack.lib" ^
  -D TPL_BLAS_LIBRARIES="%USERPROFILE%\opt\lapack\lib\libblas.lib" ^
  -D TPL_PARMETIS_INCLUDE_DIRS="%USERPROFILE%\opt\parmetis\include" ^
  -D TPL_PARMETIS_LIBRARIES="%USERPROFILE%\opt\parmetis\lib\parmetis.lib;%USERPROFILE%\opt\parmetis\lib\metis.lib" ^
  %SUPERLU_HOME%

if !ERRORLEVEL! neq 0 (
    echo ERROR: installing superlu_dist failed at configuration
    exit /b !ERRORLEVEL!
)

ninja install -j%NPROCS%

if !ERRORLEVEL! neq 0 (
    echo ERROR: installing superlu_dist failed at compilation
    exit /b !ERRORLEVEL!
)

rem

cd ..

rem Clean up downloaded and extracted artifacts
del /f /q *.tar.gz 2>nul
for /d %%D in (superlu_dist*) do rmdir /s /q "%%D"

endlocal
