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

rem clone the repository
git clone https://codeberg.org/vryy/cln.git

rem compiling

set "CLN_HOME=%~dp0cln"

mkdir "cln-build"

cd "cln-build"

cmake ^
  -G "Ninja" ^
  -D CMAKE_BUILD_TYPE:STRING="Release" ^
  -D CMAKE_INSTALL_PREFIX="%INSTALL_DIR%" ^
  -D CMAKE_C_COMPILER="cl.exe" ^
  -D CMAKE_CXX_COMPILER="cl.exe" ^
  -D BUILD_SHARED_LIBS:BOOL=OFF ^
  -D CLN_USE_GMP:BOOL=OFF ^
  %CLN_HOME%

ninja install -j%NPROCS%

rem

cd ..

rem Clean up downloaded and extracted artifacts
for /d %%D in (cln*) do rmdir /s /q "%%D"

if !ERRORLEVEL! neq 0 (
    echo ERROR: install.bat failed with exit code !ERRORLEVEL!
    exit /b !ERRORLEVEL!
)

endlocal
