% avalzero.mp
% L. Nobre G. 
% 1999

verbatimtex \documentclass{article} \usepackage{palatino} \usepackage{concmath} \usepackage{ccfonts} \begin{document} etex

beginfig(1); 						% crescente
	numeric u;
	u = 1in;
	draw fullcircle scaled 2u;
	draw (left--right) scaled u;
	draw (down--up) scaled u;
	draw (left--up) scaled u;
	draw (left--down) scaled u;
	path int, cresc, out;
	int = halfcircle rotated -90 scaled (2*1.414u) shifted ((-0.999,0)*u);
	out = halfcircle rotated -90 scaled 2u;
	cresc = buildcycle( int, reverse out ); 
	fill cresc withcolor 0.5white;
	draw cresc withpen pencircle scaled 1pt;
endfig;

beginfig(2); 						% trapezio
	numeric u, ind, esq[], rig[], h[], ch;
	pair hlab[];
	u = 0.55in;
	ch = 0;
	for ind=0 upto 4:
		esq[ind] = -abs(((normaldeviate)-1.5)*u); 
		rig[ind] = abs(((normaldeviate)+1.5)*u); 
		h[ind] = (abs((normaldeviate)+0.25)*u/4)+2mm+ch;
		ch := ch+h[ind];
	endfor; 
	draw (esq0,h0)--(rig0,h0);
	for ind=1 upto 4:
		draw (esq[ind-1],h[ind-1])--(esq[ind],h[ind]);
		draw (rig[ind-1],h[ind-1])--(rig[ind],h[ind]);
		draw ((esq[ind],h[ind])--(rig[ind],h[ind])) dashed evenly;
		hlab[ind] = (0,(h[ind]+h[ind-1])/2);
	endfor;
	draw (esq4,h4)--(rig4,h4);
	draw ((0,h0)--(0,h4)) dashed evenly;
	label.rt( btex $h_1$ etex, hlab1 );
	label.rt( btex $h_2$ etex, hlab2 );
	label.rt( btex $h_3$ etex, hlab3 );
	label.rt( btex $h_4$ etex, hlab4 );
	label.top( btex $l_0$ etex, ((rig0+esq0)/2,h0) );
	label.top( btex $l_1$ etex, ((rig1+esq1)/2,h1) );
	label.top( btex $l_2$ etex, ((rig2+esq2)/2,h2) );
	label.top( btex $l_3$ etex, ((rig3+esq3)/2,h3) );
	label.top( btex $l_4$ etex, ((rig4+esq4)/2,h4) );
endfig;

beginfig(3);						% olimpiad
	numeric u;
	u=1in;
	path squr, arcir, vert, horiz, area;
	squr = unitsquare scaled u;
	arcir = halfcircle scaled 2u;
	vert = (origin--up) scaled u shifted (u/2,0);
	horiz = (origin--right) scaled u shifted (0,u/2);\
	area = buildcycle( vert, arcir, horiz );
	fill area withcolor 0.5white;
	draw area withpen pencircle scaled 1pt;
	draw squr;
	draw arcir cutafter squr;
	draw vert;
	draw horiz;
endfig;

beginfig(4);						% ortoproj
	numeric u;
	u=1in;
	draw unitsquare scaled u;
	path one, two;
	one = (0.5,0)--(0.5,0.5)--(1,1);
	two = (0,1)--(0.5,0.5);
	draw one scaled u;
	draw two scaled u;
	draw currentpicture shifted (0,1.5u);
endfig;

verbatimtex \end{document} etex

end;
