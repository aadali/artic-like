% embroncordada.mp
% L. Nobre G. 
% 2004

input featpost3Dplus2D;

% bash: export TEX=latex
% tcsh: setenv TEX latex

verbatimtex 
\documentclass{article}
%\usepackage{pandora}
\usepackage{palatino}
\begin{document} 
etex

beginfig(1);
  picture embroncordada;
  numeric u, v;
  path pa;
  u = 1mm;
  v = 1.1mm;
  z1 = (105u,250v);
  z2 = (105u,220v);
  z3 = (125u,200v);
  z4 = (85u,180v);
  z5 = (25u,150v);
  z6 = (105u,120v);
  z7 = (185u,150v);
  z8 = (110u,170v);
  z9 = (45u,150v);
  z10= (100u,105v);
  z11= (170u,150v);
  z12= (90u,195v);
  z13= (35u,150v);
  z14= (120u,135v);
  z15= (80u,100v);
  z16= (115u,75v);
  z17= (105u,50v);
  pa = z1{down}..z2..z3..z4..z5..z6..z7..z8..z9..z10..
       z11..z12..z13..z14..z15..z16..{down}z17;
  ropepattern( pa, 8u, 35 );
  embroncordada = thelabel( btex
        \begin{tabular}{c}
          \textsf{EMBRON} \\
          \textsf{CORDADA}
	\end{tabular}
    etex scaled 8, (105u,150v) );
  draw embroncordada;
%  dotlabels(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17);
endfig;

verbatimtex 
\end{document} 
etex

end;
