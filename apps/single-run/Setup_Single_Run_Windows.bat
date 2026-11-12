@echo off
setlocal
cd /d "%~dp0"
echo Mironova Single Run - one-time setup
echo This application uses 64-bit Python 3.12. Internet is needed for setup.
echo Setup details will be saved to setup-log.txt in this folder.
py -3.12 -c "import sys,struct; sys.exit(0 if sys.version_info[:2] == (3,12) and struct.calcsize('P') == 8 else 1)" > setup-log.txt 2>&1
if errorlevel 1 goto pythonneeded
if exist ".venv312\Scripts\python.exe" goto checkenv
echo Creating the application's Python 3.12 environment...
py -3.12 -m venv .venv312 >> setup-log.txt 2>&1
if errorlevel 1 goto failed
:checkenv
".venv312\Scripts\python.exe" -c "import sys,struct; print(sys.version); sys.exit(0 if sys.version_info[:2] == (3,12) and struct.calcsize('P') == 8 else 1)" >> setup-log.txt 2>&1
if errorlevel 1 goto wrongenv
set "PIP_ONLY_BINARY=:all:"
echo Preparing installation tools...
".venv312\Scripts\python.exe" -m pip install --disable-pip-version-check --upgrade pip setuptools wheel >> setup-log.txt 2>&1
if errorlevel 1 goto failed
echo Installing prebuilt dependencies. No C++ compiler is needed...
".venv312\Scripts\python.exe" -m pip install --disable-pip-version-check -r requirements_windows.txt >> setup-log.txt 2>&1
if errorlevel 1 goto failed
echo Installing the application...
".venv312\Scripts\python.exe" -m pip install --disable-pip-version-check --no-deps --no-build-isolation . >> setup-log.txt 2>&1
if errorlevel 1 goto failed
echo Checking the installation...
".venv312\Scripts\python.exe" -m pip check >> setup-log.txt 2>&1
if errorlevel 1 goto failed
".venv312\Scripts\python.exe" -m unittest discover -s tests -v >> setup-log.txt 2>&1
if errorlevel 1 goto failed
".venv312\Scripts\python.exe" -m pip freeze > installed-versions.txt
echo Setup complete. Double-click Start_Single_Run.bat.
pause
exit /b 0
:pythonneeded
echo.
echo Please install Python 3.12 for Windows, 64-bit, with the Python launcher enabled.
echo Python 3.15 can stay installed. This app selects Python 3.12 explicitly.
echo Official installer: https://www.python.org/ftp/python/3.12.10/python-3.12.10-amd64.exe
echo Then run Setup_Windows.bat again.
pause
exit /b 1
:wrongenv
echo The existing .venv312 folder is not a working 64-bit Python 3.12 environment.
echo Extract this application into a NEW folder and run setup there.
goto failed
:failed
echo.
type setup-log.txt
echo.
echo Setup did not complete. Share setup-log.txt from this folder.
pause
exit /b 1
