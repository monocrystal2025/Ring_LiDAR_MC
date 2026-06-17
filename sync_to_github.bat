@echo off
setlocal EnableExtensions EnableDelayedExpansion
chcp 65001 >nul

cd /d "%~dp0"
title Sync to GitHub

echo.
echo ========================================
echo   Sync all local changes to GitHub
echo ========================================
echo.

where git >nul 2>nul
if errorlevel 1 (
    echo [ERROR] Git was not found in PATH.
    echo Please install Git for Windows or add Git to PATH.
    echo.
    pause
    exit /b 1
)

git rev-parse --is-inside-work-tree >nul 2>nul
if errorlevel 1 (
    echo [ERROR] This folder is not a Git repository.
    echo Folder: %CD%
    echo.
    pause
    exit /b 1
)

for /f "usebackq delims=" %%B in (`git branch --show-current`) do set "BRANCH=%%B"
if not defined BRANCH (
    echo [ERROR] Could not determine the current Git branch.
    echo.
    pause
    exit /b 1
)

for /f "usebackq delims=" %%R in (`git config --get remote.origin.url`) do set "REMOTE=%%R"
if not defined REMOTE (
    echo [ERROR] No remote named origin is configured.
    echo.
    pause
    exit /b 1
)

echo Repository: %CD%
echo Branch:     %BRANCH%
echo Remote:     %REMOTE%
echo.

echo [1/5] Checking local changes...
git status --short
echo.

echo [2/5] Staging all changed files...
git add -A
if errorlevel 1 (
    echo [ERROR] Failed to stage changes.
    echo.
    pause
    exit /b 1
)

git diff --cached --quiet
if not errorlevel 1 (
    echo No local changes to commit.
    echo.
    pause
    exit /b 0
)

for /f "usebackq delims=" %%T in (`powershell -NoProfile -Command "Get-Date -Format 'yyyy-MM-dd HH:mm:ss'"`) do set "NOW=%%T"

echo [3/5] Creating commit...
git commit -m "Auto sync: %NOW%"
if errorlevel 1 (
    echo [ERROR] Failed to create commit.
    echo.
    pause
    exit /b 1
)

echo.
echo [4/5] Pulling latest remote changes...
git pull --rebase origin %BRANCH%
if errorlevel 1 (
    echo.
    echo [ERROR] Pull/rebase failed. Resolve the conflict in Git, then run this file again.
    echo.
    pause
    exit /b 1
)

echo.
echo [5/5] Pushing to GitHub...
git push origin %BRANCH%
if errorlevel 1 (
    echo.
    echo [ERROR] Push failed. Check your network connection or GitHub SSH authentication.
    echo.
    pause
    exit /b 1
)

echo.
echo Done. All committed changes were pushed to GitHub.
echo.
pause
