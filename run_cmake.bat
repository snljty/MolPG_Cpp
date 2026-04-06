@echo off

set CMAKE_PREFIX_PATH=C:/eigen-3.4.1;C:/fmt-12.1.0
cmake . -B build -G "MinGW Makefiles" -D CMAKE_INSTALL_PREFIX=%CD% -LH
rem cmake --build build -j
rem cmake --install build
