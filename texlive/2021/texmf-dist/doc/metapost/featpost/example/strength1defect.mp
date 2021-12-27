% strength1defect.mp
% L. Nobre G.
% 2002

verbatimtex 
\documentclass{article} 
\usepackage{palatino} 
\usepackage{concmath} 
\usepackage{ccfonts} 
\begin{document} 
etex

color col[],cyan, magenta, yellow;
cyan=blue+green; magenta=red+blue;yellow=green+red;
col[-2]:=col[-1]:=col0:=black;
col1:=red;
col2:=(.2,.2,1); %blue
col3:=(1,.66,0); %orange
col4:=.85*green;
col5:=.9*magenta;
col6:=0.85*cyan;
col7:=.85*yellow;

beginfig(1);
	numeric ymax, xmax, i, j, noise, twist, grid, size, u, actsize;
	numeric extral, yposa, yposb, yposc;
	pair actpos, direct, one, two;
	path linea, lineb, linec, linex, liney;
	color cola, colb, colc, colx;
	picture pica, picb, picc;
	u = 0.3mm;
	xmax = 9;
	ymax = 9;
	noise = 1u;
	extral = 1;
	size = 5u;
	yposa = 8;
	yposb = 5;
	yposc = 2;
%	cola = (0.78,0.32,0.12);
%	colb = (0.18,0.82,0.12);
%	colc = (0.18,0.22,0.92);
	colx = (0.78,0.75,0.72);
	cola = col5;
	colb = col4;
	colc = col3;
	grid = size + 3*noise;
	linea = (grid*(-xmax-extral,yposa))--(grid*(xmax+extral,yposa));
	lineb = (grid*(-xmax-extral,yposb))--(grid*(xmax+extral,yposb));
	linec = (grid*(-xmax-extral,yposc))--(grid*(xmax+extral,yposc));
	linex = (grid*(-xmax-2*extral,0))--(grid*(xmax+2*extral,0));
	liney = (grid*(0,-ymax-2*extral))--(grid*(0,ymax+2*extral));
	pica = thelabel.rt( btex A etex, grid*(xmax+extral,yposa) );
	picb = thelabel.rt( btex B etex, grid*(xmax+extral,yposb) );
	picc = thelabel.rt( btex C etex, grid*(xmax+extral,yposc) );
	draw pica withcolor cola;
	draw picb withcolor colb;
	draw picc withcolor colc;
	label.llft( btex $x$ etex, grid*(xmax+2*extral,0) );
	label.llft( btex $y$ etex, grid*(0,ymax+2*extral) );
	pickup pencircle scaled 2u;
	draw linea withcolor cola;
	draw lineb withcolor colb;
	draw linec withcolor colc;
	pickup pencircle scaled 0.5u;
	drawarrow linex withcolor colx;
	drawarrow liney withcolor colx;
	pickup pencircle scaled u;
	for i=-xmax upto xmax:					
	    for j=-ymax upto ymax:
		actpos:=(grid*i+noise*normaldeviate,grid*j+noise*normaldeviate);
		direct:=unitvector(actpos);
		actsize:=size/2;
		one:=actpos-actsize*direct;
		two:=actpos+actsize*direct;
		draw one--two;
	    endfor; 
	endfor;
endfig;

verbatimtex 
\end{document} 
etex

end;
