@echo off
setlocal enabledelayedexpansion

REM This file is part of 4C multiphysics licensed under the
REM GNU Lesser General Public License v3.0 or later.
REM
REM See the LICENSE.md file in the top-level for license information.
REM
REM SPDX-License-Identifier: LGPL-3.0-or-later

REM Install the virtual Python environment necessary for code development in 4C.
REM Call the script from the root directory of the repository:
REM     utilities\set_up_dev_env.bat
REM Optionally, you can specify the path to the python executable:
REM     utilities\set_up_dev_env.bat C:\path\to\python.exe

if not exist "utilities\set_up_dev_env.bat" (
    echo Please run this script from the root directory of the repository.
    exit /b 1
)

REM Path to the python virtual environment (sibling of this script).
set "PYTHON_VENV=%~dp0python-venv"

REM If the virtual environment already exists, delete it.
if exist "%PYTHON_VENV%" (
    rmdir /s /q "%PYTHON_VENV%"
)

REM Path to python (default: python; override by passing an argument).
if "%~1"=="" (
    set "PYTHON_PATH=python"
) else (
    set "PYTHON_PATH=%~1"
)

REM Check that the Python version meets the minimum requirement (>=3.12).
"%PYTHON_PATH%" -c "import sys; exit(0 if sys.version_info >= (3, 12) else 1)" 2>nul
if %ERRORLEVEL% neq 0 (
    echo Provided Python version "%PYTHON_PATH%" does not meet the minimum requirement ^(^>=3.12^).
    echo Please provide a compatible Python executable as an argument to this script.
    exit /b 1
)

REM Setup the virtual environment.
"%PYTHON_PATH%" -m venv "%PYTHON_VENV%"
if %ERRORLEVEL% neq 0 (
    echo Failed to create virtual environment.
    exit /b 1
)

REM Activate the virtual environment.
call "%PYTHON_VENV%\Scripts\activate.bat"
if %ERRORLEVEL% neq 0 (
    echo Failed to activate virtual environment.
    exit /b 1
)

REM Install all the modules defined in requirements.txt.
python -m pip install --upgrade pip
if %ERRORLEVEL% neq 0 exit /b 1

pip install wheel
if %ERRORLEVEL% neq 0 exit /b 1

pip install -e utilities/four_c_python[development]
if %ERRORLEVEL% neq 0 exit /b 1

REM Additionally store the hash of the ingredients for the virtual environment.
REM check_venv is installed as a console-script entry point by the package above;
REM if it is not on PATH yet, call it explicitly via python.
call "%~dp0code_checks\windows\check_venv.bat" --update
if %ERRORLEVEL% neq 0 exit /b 1

REM We used to copy the `commit-msg` hook to `.git\hooks\` manually, but now we use
REM pre-commit to manage it.  Thus remove the old hook if it exists.
if exist ".git\hooks\commit-msg" (
    del /f /q ".git\hooks\commit-msg"
)

REM Install the pre-commit hooks.
COPY /Y .pre-commit-config.windows.yaml .pre-commit-config.yaml
pre-commit install
if %ERRORLEVEL% neq 0 exit /b 1

echo.
echo Virtual environment set up successfully.
endlocal
