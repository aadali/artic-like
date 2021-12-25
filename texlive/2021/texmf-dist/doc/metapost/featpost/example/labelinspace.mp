% labelinspace.mp
% L. Nobre G. 
% 2003

input featpost3Dplus2D;

% bash: export TEX=latex
% tcsh: setenv TEX latex

f := (2,3,0.5);
ParallelProj := true;

verbatimtex
\documentclass{article} 
%\usepackage{beton,concmath,ccfonts} 
\begin{document} 
etex

beginfig(1);
  cartaxes( 1, 1.3, 1 );
  labelinspace( false, black, 0.8*blue, 0.7*red,
    btex \framebox{$\frac{\sqrt{5}-1}{2}$} etex);
  labelinspace( false, 0.8*blue, -0.8*blue, 0.6*green,
    btex \framebox{$\phi$} etex);
  labelinspace( true, 0.7*red, 0.6*green, -0.7*red,
    btex \framebox{$\frac{1+\sqrt{5}}{2}$} etex);
endfig;

verbatimtex 
\end{document} 
etex

end;
