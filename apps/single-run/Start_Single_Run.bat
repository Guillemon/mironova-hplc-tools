@echo off
setlocal
cd /d "%~dp0"
if not exist ".venv312\Scripts\python.exe" goto setupneeded
".venv312\Scripts\python.exe" -c "import sys,struct,streamlit; from rainbow.agilent.chemstation import parse_ch,parse_ms,parse_uv; sys.exit(0 if sys.version_info[:2] == (3,12) and struct.calcsize('P') == 8 else 1)" >nul 2>&1
if errorlevel 1 goto setupneeded
set "AUTOQUANT_DATA=%~dp0single_run_results"
set "MPLCONFIGDIR=%~dp0single_run_results\plot_cache"
echo Single Run UV / Mass Spec is starting.
echo If the browser does not open, use http://127.0.0.1:8502.
echo Keep this window open while using the application.
".venv312\Scripts\python.exe" -m streamlit run single_run_app.py --server.address 127.0.0.1 --server.port 8502 --browser.serverAddress 127.0.0.1 --browser.gatherUsageStats false
if errorlevel 1 pause
exit /b
:setupneeded
echo Run Setup_Windows.bat once, then start this application again.
pause
exit /b 1
