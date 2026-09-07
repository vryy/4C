@echo off
setlocal

REM This file is part of 4C multiphysics licensed under the
REM GNU Lesser General Public License v3.0 or later.
REM
REM See the LICENSE.md file in the top-level for license information.
REM
REM SPDX-License-Identifier: LGPL-3.0-or-later

REM Run the platform-independent Python implementation.
python "%~dp0..\..\four_c_python\src\four_c_precommit\commit-msg.py" %*

REM Return the exit code from Python.
exit /b %ERRORLEVEL%
