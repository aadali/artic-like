% featpostlogo.mp
% L. Nobre G.
% 2002

input featpost3Dplus2D;

verbatimtex
\documentclass{article}
\newfont{\normald}{logod10 scaled 5000}
\pagestyle{empty}
\begin{document}
etex

input boxes;

beginfig(1);
	picture logo;
	defaultdx := 3mm;
	defaultdy := 2mm;
	logo = thelabel.urt( btex {\normald FEATPOST} etex, origin );
	boxit.thelogo( logo );
	drawboxed( thelogo );
endfig;

verbatimtex
\end{document}
etex

end;
