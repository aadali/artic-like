% tiposdetranst.mp
% L. Nobre G.
% 2003

input expressg;

verbatimtex
\documentclass{article}
\usepackage{palatino}
\begin{document}
etex

beginfig(1);
	picture pi[];
	numeric altura, lagura, unid, i, alt, lag, magura, RaD;
	unid = 5mm;
	altura = 7unid;
	lagura = 14unid;
	lag = 5.6unid;
	smoothrad := 3unid;
	alt = lag;
	RaD = lag;
	magura = 0.5lagura;
	pi1 = btex \shortstack{{\Huge U}\\ uniforme} etex;
	pi2 = btex \shortstack{{\Huge H}\\ homog\'enea} etex;
	pi3 = btex \shortstack{{\Huge B}\\ bandas} etex;
	z1 = ( -magura, altura );
	z2 = ( magura, altura );
	z3 = ( 0, 0 );
	pickup thickpen;
	drawcirclebox( 1, RaD )( pi1 );
	drawcirclebox( 2, RaD )( pi2 );
	drawcirclebox( 3, RaD )( pi3 );
	drawnormalOA( 1br, 3tl );
	drawnormalOA( 3tr, 2bl );
	drawnormalOA( 1mr, 2ml );
endfig;

verbatimtex
\end{document}
etex

end.
