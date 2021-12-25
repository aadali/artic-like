@echo off

:: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: ::
:: PREPARE:
:: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: ::

:: setting the environment variable texfonts is needed by dvitomp:

:: mp by Piet van Oostrum and emTeX
:: set texfonts=c:/emtex/tfms!
:: set tex_call=call tex &plain \input
:: set mp_call=mp &plain
:: set dvitomp_call=c:\emtex\mp\mplib\dvitomp

:: fpTeX, MiKTeX, DJ Delorie
set mp_call=mpost &mpost 
set dvitomp_call=dvitomp

:: COMMON
set tex_call=call tex \input
set dvips_call=dvips -F -K
:: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: ::
:: DO JOB:
:: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: ::

%tex_call% storyx
%dvitomp_call% storyx.dvi storyx.mpx
:: it is advisable to avoid quotes -- hence the `str' operation:
%mp_call% \vardef ifi = str storyx.mpx enddef; vardef ofi = str storyx.mpy enddef; input fix_brks.mp
%mp_call% test-bop.mp
%tex_call% test-bop.tex
%dvips_call% test-bop.dvi

:: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: ::
:: CLEAN:
:: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: :: ::

del storyx.dvi
del storyx.log
:: del storyx.mpx
:: del storyx.mpy
del test-bop.dvi
del test-bop.log
del test-bop.10*
