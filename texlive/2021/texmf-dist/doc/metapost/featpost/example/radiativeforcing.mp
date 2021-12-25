%
% fig2dev (version 3.2.5-alpha7) -L (m)mp version 0.05 --- Preamble
%

%
% mp output driver options:
%  -L mp
%

% +MP-ADDITIONAL-HEADER
verbatimtex
%&latex
\documentclass{article}
\begin{document}
etex
% -MP-ADDITIONAL-HEADER


%SetFigFont macros for latex
verbatimtex
%
\begingroup\makeatletter\ifx\SetFigFontNFSS\undefined%
\gdef\SetFigFontNFSS#1#2#3#4#5{%
  \reset@font\fontsize{#1}{#2pt}%
  \fontfamily{#3}\fontseries{#4}\fontshape{#5}%
  \selectfont}%
\fi\endgroup%
\ifx\SetFigFontSize\undefined%
\gdef\SetFigFontSize#1#2{%
  \fontsize{#1}{#2pt}%
  \selectfont}%
\fi%
etex

% Now draw the figure
beginfig(0)
% Some reasonable defaults
  labeloffset:=0;
  truecorners:=1;
  bboxmargin:=0;
% Begin ellipse object
  pickup pencircle scaled 0.94;
  path c;
  c = fullcircle scaled 1638.80 yscaled 1.00
         rotated 0.00 shifted (428.03,-242.65);
 fill c withcolor ((0.53,0.81,1.00) + 0.00white);
  draw c withcolor (0.53,0.81,1.00);
% End ellipse object
% Begin ellipse object
  pickup pencircle scaled 0.94;
  path c;
  c = fullcircle scaled 1126.68 yscaled 1.00
         rotated 0.00 shifted (428.03,-242.65);
 fill c withcolor ((0.00,0.69,0.00) + 0.00white);
  draw c withcolor (0.00,0.69,0.00);
% End ellipse object
% Begin ellipse object
  pickup pencircle scaled 0.94;
  path c;
  c = fullcircle scaled 141.73 yscaled 1.00
         rotated 0.00 shifted (147.40,650.27);
 fill c withcolor ((1.00,1.00,0.00) + 0.00white);
  draw c withcolor (1.00,1.00,0.00);
% End ellipse object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 2.83;
  path p;
  p = (155.91, 568.06)
    --(238.11, 301.61);
  draw p withcolor (1.00,1.00,0.00);
% Draw arrowhead type 1
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 2.83;
  path arr;
  arr = (231.94, 311.87)
      --(237.17, 304.57)
      --(237.35, 313.57)
      --(231.94, 311.87) -- cycle;
  fill arr withcolor (1.00,1.00,0.00);
  draw arr withcolor (1.00,1.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 2.83;
  path p;
  p = (184.25, 573.73)
    --(260.79, 315.78);
  draw p withcolor (1.00,1.00,0.00);
% Draw arrowhead type 1
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 2.83;
  path arr;
  arr = (254.74, 326.17)
      --(259.91, 318.80)
      --(260.16, 327.81)
      --(254.74, 326.17) -- cycle;
  fill arr withcolor (1.00,1.00,0.00);
  draw arr withcolor (1.00,1.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 2.83;
  path p;
  p = (212.60, 596.41)
    --(286.30, 321.45);
  draw p withcolor (1.00,1.00,0.00);
% Draw arrowhead type 1
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 2.83;
  path arr;
  arr = (280.50, 331.97)
      --(285.48, 324.47)
      --(285.98, 333.48)
      --(280.50, 331.97) -- cycle;
  fill arr withcolor (1.00,1.00,0.00);
  draw arr withcolor (1.00,1.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 2.83;
  path p;
  p = (121.89, 562.39)
    --(204.09, 293.10);
  draw p withcolor (1.00,1.00,0.00);
% Draw arrowhead type 1
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 2.83;
  path arr;
  arr = (197.98, 303.37)
      --(203.21, 296.06)
      --(203.40, 305.07)
      --(197.98, 303.37) -- cycle;
  fill arr withcolor (1.00,1.00,0.00);
  draw arr withcolor (1.00,1.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 2.83;
  path p;
  p = (351.50, 344.13)
    --(606.61, 774.99);
  draw p withcolor (1.00,0.00,0.00);
% Draw arrowhead type 1
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 2.83;
  path arr;
  arr = (603.09, 763.53)
      --(605.04, 772.28)
      --(598.24, 766.43)
      --(603.09, 763.53) -- cycle;
  fill arr withcolor (1.00,0.00,0.00);
  draw arr withcolor (1.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 2.83;
  path p;
  p = (311.81, 332.79)
    --(572.60, 766.49);
  draw p withcolor (1.00,0.00,0.00);
% Draw arrowhead type 1
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 2.83;
  path arr;
  arr = (568.94, 755.15)
      --(570.96, 763.84)
      --(564.09, 758.05)
      --(568.94, 755.15) -- cycle;
  fill arr withcolor (1.00,0.00,0.00);
  draw arr withcolor (1.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 2.83;
  path p;
  p = (391.18, 346.96)
    --(513.07, 545.39)
    --(632.13, 352.63);
  draw p withcolor (1.00,0.00,0.00);
% Draw arrowhead type 1
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 2.83;
  path arr;
  arr = (623.56, 361.07)
      --(630.49, 355.28)
      --(628.41, 364.03)
      --(623.56, 361.07) -- cycle;
  fill arr withcolor (1.00,0.00,0.00);
  draw arr withcolor (1.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 2.83;
  path p;
  p = (428.03, 349.80)
    --(549.92, 536.88)
    --(615.12, 434.83);
  draw p withcolor (1.00,0.00,0.00);
% Draw arrowhead type 1
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 2.83;
  path arr;
  arr = (606.43, 443.15)
      --(613.42, 437.48)
      --(611.21, 446.17)
      --(606.43, 443.15) -- cycle;
  fill arr withcolor (1.00,0.00,0.00);
  draw arr withcolor (1.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 2.83;
  path p;
  p = (2.83, 789.17)
    --(839.06, 789.17)
    --(839.06, 202.39)
    --(2.83, 202.39)--cycle;
  picture tmppic;
  tmppic = currentpicture;
  currentpicture := nullpicture;
  clip tmppic to p;
  draw tmppic;
% End polyline object
endfig;
end
