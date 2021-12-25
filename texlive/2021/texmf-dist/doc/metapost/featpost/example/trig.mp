verbatimtex 
\documentclass{article} 
\usepackage{palatino} 
\usepackage{concmath} 
\usepackage{ccfonts} 
\begin{document} 
etex

%prologues := 2;


vardef trigonometry(expr ct)=

	save	z, p, m, pic;
	numeric m;
	pair	z[];
	path	p[];
	picture pic;

	pic := nullpicture;
	m := if (xpart ct) > (ypart ct): ypart ct else: xpart ct fi;

	p0 = unitsquare xscaled (xpart ct) yscaled (ypart ct);
	p3 = origin--ct;
	p2 = (m,0)--(m,m)--(0,m);
	p1 = (m,0){up}..{left}(0,m);

	z0 = origin;
	z2 = point 0 of p2;
	z5 = point 2 of p2;
	z6 = ct;
	z3 = p1 intersectionpoint p3;
	z4 = p2 intersectionpoint p3;
	z1 = (xpart z3, 0);

	p4 = z1--z3;
	m := xpart (p2 intersectiontimes p3);
	if m <= 1:
		p5 = z2--z4;
		p6 = z5--z6;
	else:
		p5 = z5--z4;
		p6 = z2--z6;
	fi;

	m := top(0)-bot(0);

	for i=0 upto 4:	
		addto pic doublepath p[i] withpen pencircle scaled (m/3);
	endfor;
	for i=3 upto 6:	
		addto pic doublepath p[i] withpen pencircle scaled m;	
	endfor;

	for i=0 upto 6:	
		addto pic doublepath z[i] withpen pencircle scaled 3m;	
	endfor;

	addto pic also thelabel.llft(btex \strut $a$ etex , z0);
	addto pic also thelabel.bot(btex \strut $b$ etex, z1);
	addto pic also thelabel.bot(btex \strut $c$ etex, z2);
	addto pic also thelabel.lft(btex $d$ etex, z3);
	addto pic also thelabel.ulft(btex $e$ etex, z4);
	addto pic also thelabel.lft(btex \strut $f$ etex, z5);
	addto pic also thelabel.rt(btex \strut $g$ etex, z6);

	pic
enddef;




beginfig(1);

pickup pencircle scaled .6pt;

draw trigonometry((6cm,3cm));
draw trigonometry((3cm,6cm))shifted (7cm,0);

endfig;
verbatimtex 
\end{document} 
etex
end
