:: made in BOP s.c., bop@bop.com.pl; public domain software
:: history: 20.08.1999: ver. 0.1, pre-release
::          05.01.2000: ver. 0.2, first release
@echo off
if "%2"=="" goto USAGE
gawk -f epsincl.awk %1 > %2
goto DONE
:USAGE
echo USAGE:
echo   epsincl inputfile outputfile
:DONE
