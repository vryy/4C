@echo off
setlocal

REM This file is part of 4C multiphysics licensed under the
REM GNU Lesser General Public License v3.0 or later.
REM
REM See the LICENSE.md file in the top-level for license information.
REM
REM SPDX-License-Identifier: LGPL-3.0-or-later

REM Check that Python is available.
python --version >nul 2>&1
if errorlevel 1 (
    echo Python was not found on PATH.
    exit /b 1
)

REM Run the platform-independent implementation.
python "%~dp0..\..\four_c_python\src\four_c_precommit\check_non_ascii_chars.py" %*

exit /b %ERRORLEVEL%
