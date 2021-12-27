% labelconstruct.mp
% L. Nobre G. 
% 2004

input featpost3Dplus2D;

% bash: export TEX=latex
% tcsh: setenv TEX latex

f := 1.1*(2,1,0.5);
ParallelProj := true;

verbatimtex
\documentclass{article} 
\usepackage{beton,concmath,ccfonts} 
\begin{document} 
etex

beginfig(1);
  kindofcube(false,true,(0,-0.5,0),90,0,0,1.2,0.1,0.4);
  kindofcube(false,true,(0,0,0),0,0,0,0.5,0.1,0.8);
  labelinspace(false,(0.45,0.1,0.65),(-0.4,0,0),(0,0,0.1),
    btex \framebox{\textsc{Label}} etex);
endfig;

verbatimtex 
\end{document} 
etex

end;
