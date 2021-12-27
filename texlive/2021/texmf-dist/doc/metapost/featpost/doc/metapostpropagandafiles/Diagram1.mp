% Metapost TeX macro
% Title: /home/nobre-bu/mydirs/www/nobre/down/featpost/doc/metapostpropagandafiles/Diagram1.dia
% Creator: Dia v0.97.1
% CreationDate: Mon Jan  9 10:37:50 2012
% For: nobre


beginfig(1);
verbatimtex
%&latex
\documentclass{minimal}
\begin{document}
etex
% Define macro for horizontal centering.
vardef hcentered primary P =
  P shifted -(xpart center P, 0)
enddef;
% Define macro for right justification.
vardef rjust primary P =
  P shifted -(xpart (lrcorner P - llcorner P), 0)
enddef;
  % picture(1.950000,-12.600000)(13.850000,-1.350000)
  x = 1.000000cm; y = -1.000000cm;

  t = 1.000000;

% set_line_color 0.000000, 0.000000, 0.000000
% set_line_color 1.000000, 1.000000, 1.000000
  path p;
  p = (6.850000x,3.425000y)..(4.425000x,5.450000y)..(2.000000x,3.425000y)..(4.425000x,1.400000y)..cycle;
  fill p withcolor (1.000000,1.000000,1.000000);
% set_linewidth 0.100000
% set_line_color 0.000000, 0.000000, 0.000000
  draw (6.850000x,3.425000y)..(4.425000x,5.450000y)..(2.000000x,3.425000y)..(4.425000x,1.400000y)..cycle
    withpen pencircle scaled 0.1000x;
% set_linewidth 0.100000
% set_linewidth 0.100000
  path p;
  p = (8.300000x,8.758333y)..controls (9.400000x,8.189583y) and (9.950000x,8.000000y)
    ..(11.050000x,8.000000y)..controls (12.150000x,8.000000y) and (12.700000x,8.189583y)
    ..(13.800000x,8.758333y)--(13.800000x,11.791667y)..controls (12.700000x,12.360417y) and (12.150000x,12.550000y)
    ..(11.050000x,12.550000y)..controls (9.950000x,12.550000y) and (9.400000x,12.360417y)
    ..(8.300000x,11.791667y)--(8.300000x,8.758333y)
    ..cycle;
  fill p withcolor (1.000000,1.000000,1.000000);
% set_line_color 0.000000, 0.000000, 0.000000
  draw (8.300000x,8.758333y)..controls (9.400000x,8.189583y) and (9.950000x,8.000000y)
    ..(11.050000x,8.000000y)..controls (12.150000x,8.000000y) and (12.700000x,8.189583y)
    ..(13.800000x,8.758333y)--(13.800000x,11.791667y)..controls (12.700000x,12.360417y) and (12.150000x,12.550000y)
    ..(11.050000x,12.550000y)..controls (9.950000x,12.550000y) and (9.400000x,12.360417y)
    ..(8.300000x,11.791667y)--(8.300000x,8.758333y)
    withpen pencircle scaled 0.1000x;
% set_line_color 0.000000, 0.000000, 0.000000
  draw (8.300000x,8.758333y)..controls (9.400000x,9.327083y) and (9.950000x,9.516667y)
    ..(11.050000x,9.516667y)..controls (12.150000x,9.516667y) and (12.700000x,9.327083y)
    ..(13.800000x,8.758333y)
    withpen pencircle scaled 0.1000x;
% set_line_color 0.000000, 0.000000, 0.000000
  draw hcentered btex {\usefont{OT1}{phv}{m}{n} } etex scaled 1.64t shifted (11.050000x,10.854167y);
% set_linewidth 0.100000
% set_line_color 0.000000, 0.000000, 0.000000
  draw (6.850000x,3.425000y)--(11.050000x,3.425000y)--(11.050000x,7.513197y)
    withpen pencircle scaled 0.1000x;
% set_line_color 0.000000, 0.000000, 0.000000
% fill_polygon
  path p;
  p = (11.050000x,7.888197y)--(10.800000x,7.388197y)--(11.050000x,7.513197y)--(11.300000x,7.388197y)--cycle;
  fill p 
    withpen pencircle scaled 0.1000x;
% set_line_color 0.000000, 0.000000, 0.000000
  draw (11.050000x,7.888197y)--(10.800000x,7.388197y)--(11.050000x,7.513197y)--(11.300000x,7.388197y)--cycle
    withpen pencircle scaled 0.1000x;
% set_line_color 0.000000, 0.000000, 0.000000
  draw hcentered btex {\usefont{OT1}{phv}{m}{n} Hello} etex scaled 1.64t shifted (4.425000x,3.646250y);
% set_line_color 0.000000, 0.000000, 0.000000
  draw hcentered btex {\usefont{OT1}{phv}{m}{n} world!} etex scaled 1.64t shifted (11.050000x,10.496250y);
endfig;
end;
