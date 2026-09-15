$ErrorActionPreference = "Stop"

# ---------------------------------------------------------------------------
# Intel oneAPI Toolkit version
# ---------------------------------------------------------------------------

$Version = "2026.1.0.191"

$FileName = "intel-oneapi-toolkit-${Version}_offline.exe"

$Url = "https://registrationcenter-download.intel.com/akdlm/IRC_NAS/4144bec3-82ce-4672-bd71-5c93a79cd5e7/intel-oneapi-toolkit-2026.1.0.191_offline.exe"

$DownloadDir = Join-Path $env:RUNNER_TEMP "intel-oneapi"

$Installer = Join-Path $DownloadDir $FileName

$OneApiRoot = "C:\Program Files (x86)\Intel\oneAPI"

# ---------------------------------------------------------------------------
# Download
# ---------------------------------------------------------------------------

New-Item -ItemType Directory -Force -Path $DownloadDir | Out-Null

Write-Host "Downloading Intel oneAPI Toolkit $Version..."
Write-Host "URL: $Url"
Write-Host "Destination: $Installer"

Invoke-WebRequest `
    -Uri $Url `
    -OutFile $Installer

if (!(Test-Path $Installer)) {
    throw "Intel oneAPI installer was not downloaded."
}

Write-Host "Download complete."

# ---------------------------------------------------------------------------
# Install
# ---------------------------------------------------------------------------

Write-Host "Installing Intel oneAPI Toolkit $Version..."

$p = Start-Process `
    -FilePath $Installer `
    -ArgumentList @(
        "-s",
        "-a",
        "--silent",
        "--eula", "accept"
    ) `
    -Wait `
    -PassThru

if ($p.ExitCode -ne 0) {
    throw "Intel oneAPI installation failed with exit code $($p.ExitCode)"
}

Write-Host "Intel oneAPI installation completed."

# ---------------------------------------------------------------------------
# Check installation
# ---------------------------------------------------------------------------

$SetVars = Join-Path $OneApiRoot "setvars.bat"

if (!(Test-Path $SetVars)) {
    throw "setvars.bat was not found: $SetVars"
}

Write-Host "Found: $SetVars"

# ---------------------------------------------------------------------------
# Initialize oneAPI environment
#
# setvars.bat modifies the environment of its cmd.exe process.
# Capture the resulting environment and import it into PowerShell.
# ---------------------------------------------------------------------------

Write-Host "Initializing Intel oneAPI environment..."

$envDump = & cmd.exe /c "`"$SetVars`" >nul 2>&1 && set"

if ($LASTEXITCODE -ne 0) {
    throw "Intel oneAPI setvars.bat failed."
}

foreach ($line in $envDump) {
    if ($line -match '^([^=]+)=(.*)$') {
        [Environment]::SetEnvironmentVariable(
            $matches[1],
            $matches[2],
            "Process"
        )
    }
}

# ---------------------------------------------------------------------------
# Verify ifx
# ---------------------------------------------------------------------------

$Ifx = Get-Command ifx.exe -ErrorAction SilentlyContinue

if ($null -eq $Ifx) {
    throw "ifx.exe was not found after initializing oneAPI."
}

Write-Host ""
Write-Host "=========================================="
Write-Host "Intel Fortran Compiler"
Write-Host "=========================================="

& ifx --version

if ($LASTEXITCODE -ne 0) {
    throw "ifx --version failed."
}

Write-Host ""
Write-Host "ifx found at:"
Write-Host $Ifx.Source

# ---------------------------------------------------------------------------
# Export environment to subsequent GitHub Actions steps
# ---------------------------------------------------------------------------

if ($env:GITHUB_ENV) {

    Write-Host ""
    Write-Host "Exporting oneAPI environment to GITHUB_ENV..."

    foreach ($line in $envDump) {
        Add-Content `
            -Path $env:GITHUB_ENV `
            -Value $line
    }
}

Write-Host ""
Write-Host "Intel oneAPI / ifx setup complete."
