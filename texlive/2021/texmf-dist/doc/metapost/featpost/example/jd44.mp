% jd44.mp
% J. Dinis
% long ago
% BUG: With MetaPost, version 1.208 (kpathsea version 5.0.0) points z5 and z7
% (the "elbows") produce wrong figures !!!!

verbatimtex 
\documentclass{article} 
\usepackage{palatino} 
\begin{document} 
etex

% prologues := 2;


def drawpoints(expr a, b)=
  for i=a upto b:
	label(decimal i, z[i]);
  endfor
enddef;


def drawtopic expr p = addto pic doublepath p withpen currentpen enddef;


vardef climber(expr sc) =
	clearxy;
	save a,b,c,d,e,f, pic;
	picture pic;
	pic := nullpicture;
	
	a := .8 sc;
	b := 1 sc;
	c := 1.2 sc;
	d := .8 sc;
	e := .3 sc;
	f := .6 sc;
	
	z4 = z3 + (e,0) rotated 45;
	z1 = z4 + (c,0) rotated 150;
	z1 - z2 = z3 - z4;
	z13 = .5[z1,z2];
	z5 = z13 + (0,a);
	z6 = z5 + (0,a);
	z7 = z2 + (f,0) rotated 45;
	z8 = z7 + (a,0) rotated 45;
	z9 = z3 + (b,0);
	z10 = z9 + (b,0) rotated -45;
	z11 = z4 + (a,0) rotated 70;
	z12 = z11 + (b,0) rotated -30;
	z0 = z1 + (d,0) rotated (100+angle(z2-z1));
	
	x4-(.25pt*sc) = 0;
	y4-(.25pt*sc) = 0;
	
	pickup pencircle scaled (.5pt*sc);
	drawtopic z13--z6;
	drawtopic z3--z9--z10;
	
	pickup pencircle scaled (.6pt*sc);
	drawtopic z2--z8 withcolor white;
	drawtopic z4--z11--z12  withcolor white;
	drawtopic z1--z2--z4--z3--cycle  withcolor white;
	
	pickup pencircle scaled (.5pt*sc);
	drawtopic z1--z2--z4--z3--cycle;
	drawtopic z2--z7--z8;
	drawtopic z4--z11--z12;
	
	pickup pencircle scaled (1pt*sc);
	drawtopic z0;
	
%	drawpoints(0,12);
	pic
enddef;


vardef windsurfer(expr sc) =
	clearxy;
	save a,b,c,d,e,f, pic;
	picture pic;
	pic := nullpicture;
	
	a := .9 sc;
	b := .6 sc;
	c := 1.2 sc;
	d := .8 sc;
	e := .3 sc;
	f := .6 sc;
	
	z3 = z4 + (e,0) rotated 170;
	z2 = z4 + (c,0) rotated 80;
	z1-z2 = z3-z4;
	
	z7 = z2 + (b,0) rotated 200;
	z8 = z7 + (a,0) rotated 200;
	z5 = z1 + (b,0) rotated 200;
	z6 = z5 + (a,0) rotated 200;
	
	z9 = z3 + (d,0) rotated 195;
	z10 = z9 + (a,0) rotated 240;
	z11 = z4 + (d,0) rotated 195;
	z12 = z11 + (a,0) rotated 240;
	
	z0 = z2 + (d,0) rotated 90;
	
	x10-(.25pt*sc) = 0;
	y12-(.25pt*sc) = 0;
	
	pickup pencircle scaled (.5pt*sc);
	drawtopic z1--z5--z6;
	drawtopic z3--z9--z10;

	pickup pencircle scaled (.6pt*sc);
	drawtopic z2--z7--z8 withcolor white;
	drawtopic z4--z11--z12  withcolor white;
	drawtopic z1--z2--z4--z3--cycle  withcolor white;
	
	pickup pencircle scaled (.5pt*sc);
	drawtopic z2--z7--z8;
	drawtopic z4--z11--z12;
	drawtopic z1--z2--z4--z3--cycle;
	
	pickup pencircle scaled (1pt*sc);
	drawtopic z0;
	
 %  drawpoints(0,12);
   pic
enddef;


vardef fourtyfour(expr sc)=
	save f, p;
	picture f, p;
	f := btex 4 etex scaled sc;
	f := f yscaled 1.5 rotated 35 yscaled .85 rotated -39.4803;
	p := nullpicture;
	addto p also (f shifted -.44(lrcorner f - llcorner f));
	addto p also (f shifted .44(lrcorner f - llcorner f));
	p
enddef;


beginfig(1);
  draw climber(100) shifted (288,144);
endfig;


beginfig(2);
  draw windsurfer(100) shifted (72,72);
endfig;


beginfig(3);
	draw fourtyfour(40) shifted (144,144);
endfig;


% e à quarta é de vez.
%
beginfig(4);
	draw fourtyfour(40);
	draw climber(50) shifted llcorner currentpicture;
	draw windsurfer(50) shifted (0.8*(lrcorner currentpicture));
%	currentpicture := currentpicture shifted (80mm,120mm);
endfig;

verbatimtex 
\end{document} 
etex

end
