@echo off
setlocal EnableExtensions EnableDelayedExpansion

REM This file is part of 4C multiphysics licensed under the
REM GNU Lesser General Public License v3.0 or later.
REM
REM See the LICENSE.md file in the top-level for license information.
REM
REM SPDX-License-Identifier: LGPL-3.0-or-later

REM ---------------------------------------------------------------------------
REM check_venv.bat
REM
REM Checks whether the Python virtual environment is up to date.
REM
REM Usage:
REM     utilities\code_checks\check_venv.bat
REM     utilities\code_checks\check_venv.bat --update
REM
REM Must be called from the repository root.
REM ---------------------------------------------------------------------------

if not exist "utilities\code_checks\windows\check_venv.bat" (
    echo Please run this script from the root directory of the repository.
    exit /b 1
)

set "STORED_HASH_FILE=utilities\python-venv\_venv_hash.txt"
set "HASH_LIST_FILE=%TEMP%\4c_venv_hashes_%RANDOM%.txt"
set "HASH_FILE=%TEMP%\4c_venv_hash_%RANDOM%.txt"

REM Make sure temporary files do not already exist.
del /q "%HASH_LIST_FILE%" "%HASH_FILE%" 2>nul

REM ---------------------------------------------------------------------------
REM Hash utilities\four_c_python recursively.
REM
REM This corresponds approximately to:
REM
REM find ./utilities/four_c_python/ \
REM      -not -wholename '*.egg-info*' \
REM      -not -wholename '*__pycache__*' \
REM      -type f \
REM      -exec sha256sum {} \;
REM
REM ---------------------------------------------------------------------------

for /r "utilities\four_c_python" %%F in (*) do (
    if exist "%%F" (
        set "FILE=%%F"

        REM Exclude .egg-info paths.
        echo !FILE! | findstr /i /c:".egg-info" >nul
        if errorlevel 1 (
            REM Exclude __pycache__ paths.
            echo !FILE! | findstr /i /c:"__pycache__" >nul
            if errorlevel 1 (

                set "FILEHASH="

                REM certutil prints:
                REM
                REM SHA256 hash of file:
                REM abcdef...
                REM CertUtil: -hashfile command completed successfully.
                REM
                REM Extract only the actual hash.
                for /f "skip=1 tokens=1" %%H in (
                    'certutil -hashfile "!FILE!" SHA256 2^>nul'
                ) do (
                    if not defined FILEHASH set "FILEHASH=%%H"
                )

                if not defined FILEHASH (
                    echo Failed to calculate SHA-256 for:
                    echo !FILE!
                    del /q "%HASH_LIST_FILE%" "%HASH_FILE%" 2>nul
                    exit /b 1
                )

                REM Convert the path to a repository-relative path.
                set "RELATIVE=!FILE:%CD%\=!"

                REM Bash sha256sum produces:
                REM
                REM <hash>  <filename>
                echo !FILEHASH!  ./!RELATIVE!>>"%HASH_LIST_FILE%"
            )
        )
    )
)

REM ---------------------------------------------------------------------------
REM Include utilities\set_up_dev_env.bat.
REM
REM The original Bash script included set_up_dev_env.sh. Since this is the
REM Windows version, we include set_up_dev_env.bat instead.
REM ---------------------------------------------------------------------------

set "FILE=utilities\set_up_dev_env.bat"
set "FILEHASH="

for /f "skip=1 tokens=1" %%H in (
    'certutil -hashfile "!FILE!" SHA256 2^>nul'
) do (
    if not defined FILEHASH set "FILEHASH=%%H"
)

if not defined FILEHASH (
    echo Failed to calculate SHA-256 for:
    echo !FILE!
    del /q "%HASH_LIST_FILE%" "%HASH_FILE%" 2>nul
    exit /b 1
)

echo !FILEHASH!  ./utilities/set_up_dev_env.bat>>"%HASH_LIST_FILE%"

REM ---------------------------------------------------------------------------
REM Sort the hash list.
REM
REM The Bash version relies on find's output order. We explicitly sort here
REM so that the Windows hash is deterministic regardless of filesystem order.
REM ---------------------------------------------------------------------------

%SystemRoot%\System32\sort.exe "%HASH_LIST_FILE%" /o "%HASH_LIST_FILE%.sorted" >nul
if errorlevel 1 (
    echo Failed to sort the hash list.
    del /q "%HASH_LIST_FILE%" "%HASH_LIST_FILE%.sorted" "%HASH_FILE%" 2>nul
    exit /b 1
)

move /y "%HASH_LIST_FILE%.sorted" "%HASH_LIST_FILE%" >nul

REM ---------------------------------------------------------------------------
REM Hash the complete list of file hashes.
REM
REM This corresponds to the final:
REM
REM     ... | sha256sum | cut -d ' ' -f 1
REM ---------------------------------------------------------------------------

set "HASH="

for /f "skip=1 tokens=1" %%H in (
    'certutil -hashfile "%HASH_LIST_FILE%" SHA256 2^>nul'
) do (
    if not defined HASH set "HASH=%%H"
)

if not defined HASH (
    echo Failed to calculate the final SHA-256 hash.
    del /q "%HASH_LIST_FILE%" "%HASH_FILE%" 2>nul
    exit /b 1
)

REM ---------------------------------------------------------------------------
REM Parse command-line arguments.
REM ---------------------------------------------------------------------------

if "%~1"=="--update" goto update_hash

if not "%~1"=="" (
    echo Usage: check_venv.bat [--update]
    del /q "%HASH_LIST_FILE%" "%HASH_FILE%" 2>nul
    exit /b 1
)

REM ---------------------------------------------------------------------------
REM Check that the stored hash exists.
REM ---------------------------------------------------------------------------

if not exist "%STORED_HASH_FILE%" (
    echo The hash file does not exist. Please run utilities\set_up_dev_env.bat.
    del /q "%HASH_LIST_FILE%" "%HASH_FILE%" 2>nul
    exit /b 1
)

set /p EXPECTED_HASH=<"%STORED_HASH_FILE%"

REM ---------------------------------------------------------------------------
REM Compare hashes.
REM ---------------------------------------------------------------------------

if /i not "!HASH!"=="!EXPECTED_HASH!" (
    echo Your virtual environment is out of date.
    echo Please run the following command in your source directory:
    echo utilities\set_up_dev_env.bat
    del /q "%HASH_LIST_FILE%" "%HASH_FILE%" 2>nul
    exit /b 1
)

del /q "%HASH_LIST_FILE%" "%HASH_FILE%" 2>nul

exit /b 0

:update_hash

if not "%~2"=="" (
    echo Usage: check_venv.bat [--update]
    del /q "%HASH_LIST_FILE%" "%HASH_FILE%" 2>nul
    exit /b 1
)

> "%STORED_HASH_FILE%" echo !HASH!

echo Updated the hash of the virtual environment.

del /q "%HASH_LIST_FILE%" "%HASH_FILE%" 2>nul

exit /b 0
