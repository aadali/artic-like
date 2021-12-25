% troncoconedef.mp
% L. Nobre G.
% 2002

input featpost3Dplus2D;

verbatimtex 
\documentclass{article} 
\usepackage{beton} 
\usepackage{concmath} 
\usepackage{ccfonts}
\begin{document} 
etex

f := (5,0,3);

def twoclicpathesinapen( expr CycleA, CycleB )=
    begingroup
	numeric TheLengthOfA, TheLengthOfB, TheMargin;
	path SubPathA, SubPathB;
	TheMargin = 0.02;
	TheLengthOfA = ( length CycleA ) - TheMargin;
 	TheLengthOfB = ( length CycleB ) - TheMargin;
	SubPathA = subpath ( 0, TheLengthOfA ) of CycleA;
	SubPathB = subpath ( 0, TheLengthOfB ) of CycleB;
	( makepath makepen ( SubPathA--SubPathB--cycle ) )
    endgroup
enddef;

beginfig(0);
	numeric rayo, rayt;
	path onec, twoc, botc;
	color axec, baso, bast;
	rayo = 1;
	rayt = 0.7;
	axec = red;
	baso = black;
	bast = (1,1,1);
	onec = goodcirclepath( baso, axec, rayo );
 	twoc = goodcirclepath( bast, axec, rayt ); 
 	botc = twoclicpathesinapen( onec, twoc );
	fill botc withcolor green;
	draw onec;
	draw twoc;
	draw botc;
	dotlabel.bot( btex onec0 etex, point 0 of onec ); 
	dotlabel.bot( btex twoc0 etex, point 0 of twoc ); 
endfig;

verbatimtex 
\end{document} 
etex

end.
