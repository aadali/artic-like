% tiposdetransb.mp
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
	magura = 0.5*lagura;
	pi1 = btex \shortstack{{\Huge U}\\ uniforme} etex;
	pi2 = btex \shortstack{{\Huge H}\\ homog\'enea (0)} etex;
	pi3 = btex \shortstack{{\Huge B}\\ bandas (1)} etex;
	pi4 = btex \shortstack{{\Huge D}\\ defeitos} etex;
	pi5 = btex \shortstack{{\Huge P}\\ paredes} etex;
	pi6 = btex \shortstack{{\Huge O}\\ outra (2)} etex;
	z1 = ( -lagura, altura );
	z2 = ( lagura, altura );
	z3 = ( -magura, 0 );
	z4 = ( magura, 0 );
	z5 = ( 0, -altura );
	z6 = ( -lagura, -altura );
	pickup thickpen;
	drawcirclebox( 1, RaD )( pi1 );
	drawcirclebox( 2, RaD )( pi2 );
	drawcirclebox( 3, RaD )( pi3 );
	drawcirclebox( 4, RaD )( pi4 );
	drawcirclebox( 5, RaD )( pi5 );
	drawcirclebox( 6, RaD )( pi6 );
	drawnormalOA( 1br, 3tl );
	drawnormalOA( 3tr, 2ml );
	drawnormalOA( 4tr, 2bl );
	drawnormalOA( 3br, 5tl );
	drawnormalOA( 5tr, 4bl );
	drawnormalOA( 1mr, 2ml );
	drawnormalOA( 3mr, 4ml );
	drawnormalOA( 1bm, 6tm );
endfig;

verbatimtex
\end{document}
etex

end.
