% logofontest.mp
% L. Nobre G.
% 2004

input boxes;

verbatimtex
\documentclass{article}
%\usepackage{charter}% no bold
\usepackage{newcent}
%\usepackage{palatino}
\pagestyle{empty}
\begin{document}
etex

beginfig(0);
	picture thestatement;
	thestatement = thelabel( btex
		     \begin{tabular}{c}
		       \textbf{\texttt{{\Large T}$\!$ecnumbria}} \\
		       \texttt{\scriptsize a sombra da tecnologia}
		     \end{tabular}
	etex, origin);
	circleit.state( thestatement );
	state.dx = state.dy;
	pickup pencircle scaled 1pt;
	fill bpath.state withcolor (0.25,0.33,0.75);
	draw bpath.state;
	draw thestatement;
endfig;

verbatimtex \end{document} etex

end;

