verbatimtex
\documentclass{article}
\usepackage{bookman}
\pagestyle{empty}
\begin{document}
etex

input boxes;

beginfig(0);
	picture thestatement;
	thestatement = thelabel( btex
		     \begin{tabular}{c}
			Pictures \\
			from \\
			Out There
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
