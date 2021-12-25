% cinemwork.mp
% L. Nobre G. 
% 2013

verbatimtex \documentclass{article}\begin{document} etex

def wmg = withcolor mygray enddef;

def resetpic = 
    begingroup
	tmpicture := currentpicture;
	currentpicture := nullpicture
    endgroup
enddef;

def roundedcornersbox( expr swc, nec, cornerrad ) =
    begingroup
	pair swe, swn, nes, new, nwe, nws, sen, sew;
	numeric lowx, higx, lowy, higy;
	path returnpath;
	lowx = xpart swc;
	higx = xpart nec;
	lowy = ypart swc;
	higy = ypart nec;
	swe = ( lowx+cornerrad, lowy );
	nwe = ( lowx+cornerrad, higy );
	new = ( higx-cornerrad, higy );
	sew = ( higx-cornerrad, lowy );
	swn = ( lowx, lowy+cornerrad );
	nws = ( lowx, higy-cornerrad );
	nes = ( higx, higy-cornerrad );
	sen = ( higx, lowy+cornerrad );
	returnpath = nws{up}..{right}nwe--
		     new{right}..{down}nes--
		     sen{down}..{left}sew--
		     swe{left}..{up}swn--cycle;
	( returnpath )
    endgroup
enddef;

def rectracarea( expr unit, largnumtracs, espesfrac, tangle )=
    begingroup
	picture tmpicture, returnpicture;
	numeric largarea, espesarea, i;
	pair esquarea, direarea, espessvec;
	path boxarea, traco;
	tmpicture = currentpicture;
	currentpicture := nullpicture;
	largarea = largnumtracs*unit;
	espesarea = espesfrac*largarea;
	esquarea = (-largarea/2,0);
	direarea = (largarea/2,0);
	espessvec = (0,espesarea);
	boxarea = esquarea--direarea--
		  (direarea+espessvec)--(esquarea+espessvec)--cycle;
	traco = ( origin--( dir( tangle ) ) ) scaled 6in;
	for i=-largnumtracs upto largnumtracs:
	    draw traco shifted (2*i*unit,0);
	endfor;
	clip currentpicture to boxarea;
	returnpicture = currentpicture;
	currentpicture := nullpicture;
	draw tmpicture;
	( returnpicture )
    endgroup	
enddef;

def picketfence =
    begingroup
	picture storepicture, returnpicture, minor, major, mtics;
	path borderbox;
	numeric rule;
	storepicture := currentpicture;
	currentpicture := nullpicture;
	pickup pencircle scaled 0;
	borderbox = origin--(136mm,0)--(136mm,45mm)--(0,45mm)--cycle;
	linecap := butt;
	pickup pencircle scaled 5mm;
	draw origin--(0,1cm);
	resetpic;
	minor = tmpicture;
	draw origin--(0,1.5cm);
	resetpic;
	mtics = tmpicture;
	draw origin--(0,3cm);
	resetpic;
	major = tmpicture;
	for rule = -5cm step 1cm until 5cm:
	    draw minor shifted (rule+68mm,0);
	endfor;
	for rule = -5cm step 5cm until 5cm:
	    draw major shifted (rule+68mm,0);
	endfor;
	for rule = -1 step 2 until 1:
	    draw mtics shifted (rule*5mm+68mm,3cm);
	endfor;
	linecap := rounded;
	pickup pencircle scaled 0;
	clip currentpicture to borderbox;
	draw borderbox;
	returnpicture := currentpicture rotated 180 shifted (136mm,45mm);
	currentpicture := storepicture;
	( returnpicture scaled (1/(136mm)) )
    endgroup
enddef;

def woodblock =
    begingroup
	picture storepicture, returnpicture;
	path borderbox;
	storepicture := currentpicture;
	currentpicture := nullpicture;
	borderbox = origin--(136mm,0)--(136mm,45mm)--(0,45mm)--cycle;
	fill borderbox wmg;
	returnpicture := currentpicture;
	currentpicture := storepicture;
	( returnpicture scaled (1/(136mm)) )
    endgroup
enddef;

beginfig(1);
	numeric u, ma, mb, mc, md, me, mf, mg, mh, mi, mj, mk, ml, mq;
	numeric mn, mo, mp, rangle, xzero, ystart, yzero, yend, mr;
	numeric magn, xthree, ymed;	
	picture bancada, tmpicture, storepicture;
	path carro, lwheel, rwheel, zoom;
	color mygray;
	ahangle := 30;
	ahlength := 10pt;
	u := 1mm;
	ma := 2u;
	mb := 150u;
	mc := 0.25mb;
	md := mb-mc;
	me := 8u;
	mf := 16u;
	rangle := -angle(md-mc,mf-me);
	mg := 3u;
	mh := 7u;
	mi := 1.2u;
	mj := 4mi;
	mk := 0.7mh;
	ml := 0.16mb;
	mn := mi+mj;
	mo := (0.5*mg)+(1.5*ma);
	mq := 7u;
	mr := (mf*md-me*mc)/(mf-me);
	mygray := 0.8white;
	magn := 2.5;

%	draw picketfence scaled 0.8ml shifted (mn+0.1ml,mo+mk);
	draw woodblock scaled 0.8ml shifted (mn+0.1ml,mo+mk);
  
        linecap := butt;
	pickup pencircle scaled ma;
	draw (0,0)--(mb,0) wmg;
	draw (0,0.5mg)--(mb,0.5mg);
	draw (0,-0.5mg)--(mb,-0.5mg);
	draw (mi,0)--(mi,mh);	
	draw (mb-mi,0)--(mb-mi,mh);	
	draw (mc,0)--(mc,-mf);
	draw (md,0)--(md,-me);
	linecap := rounded;
	draw (-0.5mg,0)--(-mq-2u,0) dashed evenly withpen pencircle scaled 0pt;
	path angarc;
	angarc = fullcircle cutafter (origin--1.1*dir(-rangle));
	angarc := (angarc rotated 180) scaled (2*(mq+mr)) shifted (mr,0);
	draw angarc dashed evenly withpen pencircle scaled 0pt;
	label.rt( btex $\theta$ etex scaled magn, point 0.5 of angarc );
	pickup pencircle scaled 0.67ma;
	carro := roundedcornersbox( (mn,mo),(mn+ml,mo+mk), 2u );
	draw carro;
	lwheel := fullcircle rotated 90 scaled mk shifted (mn+0.2ml,0.5*(mk+mg));
	rwheel := fullcircle rotated 90 scaled mk shifted (mn+0.8ml,0.5*(mk+mg));
	pickup pencircle scaled 0.3ma;
	draw lwheel cutafter carro cutbefore carro; 
	draw rwheel cutbefore carro cutafter carro; 

%	xzero := (0.8ml)*52.5/136+mn+(0.5ml);
	xzero := (1.0ml)*52.5/136+mn+(0.5ml);
	ystart := 3.5*(mo+mk);
	yzero := 1.05*(mo+mk)+(0.8ml)*45/136;
	pickup pencircle scaled 0.3ma;
	drawarrow (xzero,ystart)--(xzero,yzero);
	label.top( btex $x_0$ etex scaled magn, (xzero,ystart) );
	%zoom := fullcircle scaled (0.5ml) shifted (xzero,mo+mk+(0.8ml)*45/136);
	%storepicture := currentpicture;
	%clip currentpicture to zoom;
	%draw zoom;
	%resetpic;
	%draw storepicture;
	%draw tmpicture scaled 2 shifted (-5cm,5mm);	

	xthree := 0.5*(mb+mn+ml);
	yend := 2mo;
	drawarrow (xthree,ystart)--(xthree,yend);
	label.top( btex $x_3$ etex scaled magn, (xthree,ystart) );
	mp := 0.8*(mb-0.5*(mb+mn+ml));
	drawarrow (xthree-mp,ystart)--(xthree-mp,yend);
	label.top( btex $x_1$ etex scaled magn, (xthree-mp,ystart) );
	drawarrow (xthree+mp,ystart)--(xthree+mp,yend);
	label.top( btex $x_2$ etex scaled magn, (xthree+mp,ystart) );
	ymed := 0.5*(ystart+yend);
	drawdblarrow (xthree-mp,ymed)--(xthree,ymed);
	drawdblarrow (xthree,ymed)--(xthree+mp,ymed);
	label.top( btex $D/2$ etex scaled magn, (xthree-(0.5*mp),ymed) );
	label.top( btex $D/2$ etex scaled magn, (xthree+(0.5*mp),ymed) );
	resetpic;
	draw tmpicture shifted (-md,me) rotated rangle shifted (0.25mb,0);
	bancada := rectracarea( u, 1.2*mb/u, 0.06, 45 ) rotated 180;
	draw (-0.6mb,0)--(0.6mb,0);
	draw bancada;
	resetpic; 
	draw tmpicture shifted (10.5cm, 3cm );
endfig;

beginfig(2);
	numeric u, ma, mb, mc, md, me, mf, mg, mh, mi, mj, mk, ml;
	numeric mn, mo, rangle, xzero, ystart, yzero, yend;
	numeric magn, xthree, ymed, xphot;	
	picture bancada, tmpicture, storepicture;
	path carro, lwheel, rwheel, zoom;
	color mygray;
	ahangle := 30;
	ahlength := 10pt;
	u := 1mm;
	ma := 2u;
	mb := 150u;
	mc := 0.25mb;
	md := 3mc;
	me := 8u;
	mf := 11u;
	rangle := -angle(md-mc,mf-me);
	mg := 3u;
	mh := 7u;
	mi := 1.2u;
	mj := 4mi;
	mk := 0.7mh;
	ml := 0.16mb;
	mn := mi+mj;
	mo := (0.5*mg)+(1.5*ma);
	drawoptions( withcolor 0.6white );
	linecap := rounded;
	pickup pencircle scaled 0.67ma;
	carro := roundedcornersbox( (mn,mo),(mn+ml,mo+mk), 2u );
	draw carro;
	lwheel := fullcircle rotated 90 scaled mk shifted (mn+0.2ml,0.5*(mk+mg));
	rwheel := fullcircle rotated 90 scaled mk shifted (mn+0.8ml,0.5*(mk+mg));
	pickup pencircle scaled 0.3ma;
	draw lwheel cutafter carro cutbefore carro; 
	draw rwheel cutbefore carro cutafter carro; 
	draw picketfence scaled 0.8ml shifted (mn+0.1ml,mo+mk);
	drawoptions();
	mygray := 0.8white;
	linecap := butt;
	pickup pencircle scaled ma;
	draw (0,0)--(mb,0) wmg;
	draw (0,0.5mg)--(mb,0.5mg);
	draw (0,-0.5mg)--(mb,-0.5mg);
	draw (mi,0)--(mi,mh);	
	draw (mb-mi,0)--(mb-mi,mh);	
	draw (mc,0)--(mc,-mf);
	draw (md,0)--(md,-me);

	xzero := (0.8ml)*52.5/136+mn+(0.5ml);
	ystart := 3.5*(mo+mk);
	yzero := 1.05*(mo+mk)+(0.8ml)*45/136;
	yend := 2mo;
	ymed := 0.5*(ystart+yend);
	magn := 2.5;
	pickup pencircle scaled 0.3ma;
	drawarrow (xzero,ystart)--(xzero,yzero);
	label.top( btex $x_0$ etex scaled magn, (xzero,ystart) );

	pickup pencircle scaled 0.67ma;
	carro := roundedcornersbox( (mn+0.5*mb,mo),(mn+ml+0.5*mb,mo+mk), 2u );
	draw carro;
	lwheel := fullcircle rotated 90 scaled mk;
        lwheel := lwheel shifted (mn+0.2ml+0.5*mb,0.5*(mk+mg));
	rwheel := fullcircle rotated 90 scaled mk;
	rwheel := rwheel shifted (mn+0.8ml+0.5*mb,0.5*(mk+mg));
	pickup pencircle scaled 0.3ma;
	draw lwheel cutafter carro cutbefore carro; 
	draw rwheel cutbefore carro cutafter carro; 
	draw picketfence scaled 0.8ml shifted (mn+0.1ml+0.5*mb,mo+mk);

	pickup pencircle scaled 0;
	ahlength := 2pt;
	xzero := (0.8ml)*47.5/136+mn+(0.5ml)+0.5*mb;
	xphot := (0.8ml)*42.5/136+mn+(0.5ml)+0.5*mb;
	magn := 0.625;
	drawarrow (xzero,ymed)--(xzero,yzero);
	label.rt( btex $x$ etex scaled magn, (xzero,ymed) );
	drawarrow (xphot,ymed)--(xphot,yzero);
	label.lft( btex $x_1$ etex scaled magn, (xphot,ymed) );
	drawarrow (xphot-0.1ml,0.5[ymed,yzero])--(xphot,0.5[ymed,yzero]);
	draw (xphot,0.5[ymed,yzero])--(xzero,0.5[ymed,yzero]);
	drawarrow (xzero+0.1ml,0.5[ymed,yzero])--(xzero,0.5[ymed,yzero]);
	label.rt( btex $D/2$ etex scaled magn, (xzero+0.1ml,0.5[ymed,yzero]) );

	zoom := fullcircle scaled (0.7ml);
	zoom := zoom shifted (xzero,mo+mk+(0.8ml)*45/136);
	storepicture := currentpicture;
	clip currentpicture to zoom;
	draw zoom;
	resetpic;
	draw storepicture;
	draw tmpicture scaled 4 shifted (-31cm,5mm);

	resetpic;
	draw tmpicture shifted (-md,me) rotated rangle shifted (0.25mb,0);
	pickup pencircle scaled 0.3ma;
	bancada := rectracarea( u, 1.2*mb/u, 0.06, 45 ) rotated 180;
	draw (-0.6mb,0)--(0.6mb,0);
	draw bancada;
	resetpic;
	draw tmpicture shifted (10.5cm, 3cm );
endfig;

verbatimtex \end{document} etex

end;
