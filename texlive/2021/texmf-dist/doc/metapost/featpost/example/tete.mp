
verbatimtex 
\documentclass{article}
\usepackage{palatino}
\begin{document} 
etex

beginfig(1);
  picture embroncordada;
  embroncordada = thelabel( btex
        \begin{tabular}{c}
          \TeX \\
	  \& \\
	  \LaTeX
	\end{tabular}
    etex scaled 8, (105mm,150mm) );
  draw embroncordada;
endfig;

verbatimtex 
\end{document} 
etex

end;
