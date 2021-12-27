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
	numeric altura, lagura, unid, i, alt, lag, magura;
	unid = 5mm;
	altura = 7unid;
	lagura = 14unid;
	lag = 5.6unid;
	smoothrad := 3unid;
	alt = lag;
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
	drawPREF( 1, lag, alt )( pi1 );
	drawPREF( 2, lag, alt )( pi2 );
	drawPREF( 3, lag, alt )( pi3 );
	drawPREF( 4, lag, alt )( pi4 );
	drawPREF( 5, lag, alt )( pi5 );
	drawPREF( 6, lag, alt )( pi6 );
	z8 = ( xpart z1c, ypart z3c );
	z9 = ( xpart z2c, ypart z4c );
	z7 = ( xpart z3c, ypart z2c );
	z10= ( xpart z3c, ypart z5c );
	z11= ( xpart z4c, ypart z5c );
	drawnormalthree( 1bm, 8, 3ml );
	drawnormalthree( 3tm, 7, 2ml );
	drawnormalthree( 4mr, 9, 2bm );
	drawnormalthree( 3bm, 10,5ml );
	drawnormalthree( 5mr, 11,4bm );
	smooth( 1bm, 8, 3ml );
	smooth( 3tm, 7, 2ml );
	smooth( 4mr, 9, 2bm );
	smooth( 3bm, 10,5ml );
	smooth( 5mr, 11,4bm );
	drawnormalOA( 1mr, 2ml );
	drawnormalOA( 3mr, 4ml );
	drawnormalOA( 1bm, 6tm );
	drawOA( 8, 3ml );
	drawOA( 7, 2ml );
	drawOA( 9, 2bm );
	drawOA( 10,5ml );
	drawOA( 11,4bm );
endfig;

verbatimtex
\end{document}
etex

end.
