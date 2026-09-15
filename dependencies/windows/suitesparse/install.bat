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
set "VERSION=7.14.0"
set "CHECKSUM=c552c4b4bb7d0978796e57263a73295bca0c6b41ad137b45b4f264cfe9300fcb"
set "ARCHIVE=v%VERSION%.tar.gz"
set "USER_DIR=%USERPROFILE:\=\\%"
set "LIB_DIR=%USER_DIR%\\temp"

rem Download suitesparse
curl -s -L -o "%ARCHIVE%" "https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/refs/tags/%ARCHIVE%"
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

set "SUITESPARSE_HOME=%~dp0SuiteSparse-%VERSION%"

mkdir "suitesparse-%VERSION%-build"

cd "suitesparse-%VERSION%-build"

cmake ^
  -G "Ninja" ^
  -D CMAKE_BUILD_TYPE:STRING="Release" ^
  -D CMAKE_INSTALL_PREFIX="%INSTALL_DIR%" ^
  -D CMAKE_C_COMPILER="cl.exe" ^
  -D CMAKE_CXX_COMPILER="cl.exe" ^
  -D CMAKE_SHARED_LINKER_FLAGS="/FORCE:MULTIPLE" ^
  -D CMAKE_EXE_LINKER_FLAGS="/FORCE:MULTIPLE" ^
  -D BUILD_SHARED_LIBS:BOOL=OFF ^
  -D BUILD_STATIC_LIBS:BOOL=ON ^
  -D SUITESPARSE_USE_OPENMP:BOOL=ON ^
  -D SUITESPARSE_USE_FORTRAN:BOOL=ON ^
  -D BLAS_LIBRARIES="%LIB_DIR%\\lapack\\lib\\libblas.lib" ^
  -D BLA_STATIC:BOOL=ON ^
  -D BLA_VENDOR="Generic" ^
  -D LAPACK_LIBRARIES="%LIB_DIR%\\lapack\\lib\\liblapack.lib" ^
  -D SUITESPARSE_ENABLE_PROJECTS="suitesparse_config;amd;colamd;cholmod;umfpack" ^
  -D SUITESPARSE_DEMOS:BOOL=ON ^
  -D BUILD_TESTING:BOOL=ON ^
  -D SUITESPARSE_USE_FORTRAN:BOOL=OFF ^
  -D SUITESPARSE_C_TO_FORTRAN="(name,NAME) name##_" ^
  %SUITESPARSE_HOME%

if !ERRORLEVEL! neq 0 (
    echo ERROR: installing suitesparse failed at configuration
    exit /b !ERRORLEVEL!
)

ninja install -j%NPROCS%

if !ERRORLEVEL! neq 0 (
    echo ERROR: installing suitesparse failed at compilation
    exit /b !ERRORLEVEL!
)

rem

cd ..

rem Clean up downloaded and extracted artifacts
del /f /q suitesparse*.tar.gz 2>nul
for /d %%D in (suitesparse*) do rmdir /s /q "%%D"

rem rename the installed file (from*_static.lib to *.lib)

cd %INSTALL_DIR%\lib
for %%F in (*_static.lib) do (
    set "filename=%%~nF"
    ren "%%F" "!filename:~0,-7!.lib"
)

if !ERRORLEVEL! neq 0 (
    echo ERROR: install.bat failed with exit code !ERRORLEVEL!
    exit /b !ERRORLEVEL!
)

endlocal
