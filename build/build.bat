@::!/dos/rocks
@echo off
goto :parse

:parse
if "%~1"=="" goto :end

    if /i "%~1"=="/debug"       set "OptDebug=yes"  & shift & goto :parse
    if /i "%~1"=="-debug"       set "OptDebug=yes"  & shift & goto :parse
    if /i "%~1"=="--debug"      set "OptDebug=yes"  & shift & goto :parse

    if /i "%~1"=="/build"       goto :build
    if /i "%~1"=="-build"       goto :build
    if /i "%~1"=="--build"      goto :build
    
    if /i "%~1"=="/clean"       goto :clean
    if /i "%~1"=="-clean"       goto :clean
    if /i "%~1"=="--clean"      goto :clean

    shift

:build
    call %HOME%\Anaconda3\Scripts\activate.bat mitsuba_build
    cd ..
    if defined OptDebug (
        scons --parallelize --cfg=build\config-win64-msvc2017-debug.py
    ) else (
        scons --parallelize --cfg=build\config-win64-msvc2017.py
    )
    goto :end

:clean
    call %HOME%\Anaconda3\Scripts\activate.bat mitsuba_build
    cd ..
    if defined OptDebug (
        scons --parallelize --cfg=build\config-win64-msvc2017-debug.py -c
    ) else (
        scons --parallelize --cfg=build\config-win64-msvc2017.py -c
    )
    goto :end

:end
    exit /B