input boxes;

verbatimtex 
\documentclass{article} 
\usepackage{palatino} 
\usepackage{concmath} 
\usepackage{ccfonts} 
\begin{document} 
etex

defaultfont:="cmssu30";

beginfig(1);
	numeric ymax, xmax, ymed,ylow;
	ymax = 2in;
	xmax = 3in;
	ymed = 1.25in;
	ylow = 0.2in;
	z1 = (xmax,0);
	z2 = (xmax,ymed);
	z5 = (xmax,ymax);
	z3 = z2+0.4*(z5-z2);
	z4 = z2+0.7*(z5-z2);
	z7 = (0,ymax);
	z6 = z5+0.7*(z7-z5);
	z9 = (0,ymed);
	z8 = z7+0.8*(z9-z7);
	z10= z9+0.5*(origin-z9); 
	z11= (xmax/12,ymed);
	z13= (xmax/3,ymed);
	z14= (5*xmax/12,ymed);
	z15= (7*xmax/12,0.8*ymed);
	z16= (7*xmax/12,0.6*ymed);
	z17= (3*xmax/4,ymed);
	z18= (5*xmax/6,ymed);
	path a,b,c,d;
	a = z6{down}..z14..z15..z17..z4;
	b = z9--z2;
	c = z10..{right}z13..z2;
	d = z8..z11..z16..z18..z3;
	z12 = c intersectionpoint d;
	pickup pencircle scaled 2pt;
	draw origin--z1--z5--z7--cycle;
	draw a;
	draw b;
	draw c;
	draw d;
%	dotlabels.urt(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18);
	label.top(btex $h^-_{I\! I}$ etex,z15);
	label.bot(btex $h^+_{I\! I}$ etex,z16);
	label.lrt(btex $h_{I}$ etex,z10);
	label.urt(btex $h_0$ etex,(x12,ymed));
	label.bot(btex $K_3/K_1$ etex,0.5*(z1-origin));
	label.lft(btex $h$ etex,0.5*(z7-origin));
	pickup pencircle scaled 0pt;
	draw (x11,0)--(x11,ymax);
	draw (x12,0)--(x12,ymax);
	draw (x13,0)--(x13,ymax);
	draw (x14,0)--(x14,ymax);
	draw (x17,0)--(x17,ymax);
	draw (x18,0)--(x18,ymax);
	label(btex $C_1$ etex,((x7+x11)/2,ylow));
	label(btex $C_3$ etex,((x12+x11)/2,ylow));
	label(btex $C_2$ etex,((x13+x12)/2,ylow));
	label(btex $B$ etex,((x14+x13)/2,ylow));
	label(btex $A_1$ etex,((x17+x14)/2,ylow));
	label(btex $B$ etex,((x18+x17)/2,ylow));
	label(btex $A_2$ etex,((x1+x18)/2,ylow));
endfig;

beginfig(2);
	numeric ymax, xmax, i, j, noise, twist, grid, size, u, actsize;
	pair actpos, direct, one, two;
	u = 0.5mm;
	xmax = 9;
	ymax = 9;
	noise = 1u;
	size = 5u;
	pickup pencircle scaled u;
	grid = size + 3*noise;
	twist = 70/(grid*(xmax++ymax));
	label.lft( btex tor\c{c}\~{a}o etex, grid*(-xmax-1,ymax) );
	for i=-xmax upto -1:					% twist
	    for j=1 upto ymax:
		actpos:=(grid*i+noise*normaldeviate,grid*j+noise*normaldeviate);
		direct:=dir(90);
		actsize:=size*sind(3*(xpart actpos)*twist)/2;
%		direct:=dir(angle(actpos)+90);
%		actsize:=size*sind(abs(actpos)*twist)/2;
		one:=actpos-actsize*direct;
		two:=actpos+actsize*direct;
		draw one--two;
	    endfor; 
	endfor;
	label.rt( btex flex\~{a}o etex, grid*(xmax+1,ymax) );
	for i=1 upto xmax:					% bend
	    for j=1 upto ymax:
		actpos:=(grid*i+noise*normaldeviate,grid*j+noise*normaldeviate);
		direct:=dir(angle(actpos)+90);
		actsize:=size/2;
		one:=actpos-actsize*direct;
		two:=actpos+actsize*direct;
		draw one--two;
	    endfor; 
	endfor;
	label.rt( btex sela etex, grid*(xmax+1,-ymax) );
	for i=1 upto xmax:					% saddle-splay
	    for j=-ymax upto -1:
		actpos:=(grid*i+noise*normaldeviate,grid*j+noise*normaldeviate);
		direct:=dir(-angle(actpos));
		actsize:=size*sind(abs(actpos)*twist)/2;
		one:=actpos-actsize*direct;
		two:=actpos+actsize*direct;
		draw one--two;
	    endfor; 
	endfor;
	label.lft( btex afunilamento etex, grid*(-xmax-1,-ymax) );
	for i=-xmax upto -1:					% splay
	    for j=-ymax upto -1:
		actpos:=(grid*i+noise*normaldeviate,grid*j+noise*normaldeviate);
		direct:=unitvector(actpos);
%		actsize:=size/2;
		actsize:=size*sind(abs(actpos)*twist)/2;
		one:=actpos-actsize*direct;
		two:=actpos+actsize*direct;
		draw one--two;
	    endfor; 
	endfor;
endfig;

beginfig(3);
	numeric vertj;
	path ab, at, tb;
	vertj = 1cm;
	boxit.afun( btex afunilamento etex );
	boxit.bend( btex flex\~{a}o etex );
	boxit.torc( btex tor\c{c}\~{a}o etex );
	afun.se=bend.sw-vertj*up; 
        afun.ne=bend.nw-vertj*up;
	bend.se=torc.sw+2*vertj*up; 
        bend.ne=torc.nw+2*vertj*up;
	drawunboxed( afun, bend, torc ); 
	ab = afun.n{up}..{right}bend.w;
	at = afun.s{down}..{right}torc.w;
	tb = torc.n{up}..{left}bend.e;
	pickup pencircle scaled 1pt;
	drawarrow ab;
	drawarrow at;
	drawarrow tb;
	label.top( btex $\beta$ etex, point 0.5*length ab of ab ); 
	label.bot( btex $\psi$ etex, point 0.5*length at of at ); 
	label.rt( btex $\beta$ etex, point 0.5*length tb of tb ); 
endfig;

verbatimtex 
\end{document} 
etex

end;
