% experimental2Dsetup.mp
% L. Nobre G. 
% 2002

verbatimtex 
\documentclass{article} 
\usepackage{palatino} 
\usepackage{concmath} 
\usepackage{ccfonts} 
\begin{document} 
etex

beginfig(1);
%	linecap := squared;

	numeric salto, shif, rotang, cfac, posang;
	numeric innerray, outerray, midleray, thebiray, otherray;
	numeric midmidcir, midotharc, lblfac, wavecray, midwavarc;
	path loricir, roricir, baserec, sampleh, wavecarc;
	path innercirc, outercirc, midlecirc, otherarch;
	path oridirh, oridirn, oridirp, oridirv, oridire;
	path dirhoriz, dirnegat, dirposit, dirverti, direff;
	path lowcut, higcut, glassre, wavecarch, field;
	pair midmidcirpos, midotharcpos, midwavarcpos, mideffcirpos;
	pen thickpen;
	color samplecolor, magnetcolor, samplebordercolor;
	picture lblpiplu, lblalpha, lblwavec, lblqu, lblno, lblfield;
	numeric hei, wid, thi, pro, pol, mar;
	path themagnet, effcirc, diry, oridiry;
	numeric circleray, squareray, glasssize, effang, effraya, effrayb;
	picture lblx, lbly, lbleff;
	numeric mideffcir, yray;
	pair lastoridire;
	
	pickup pencircle scaled 0.15mm;
	thickpen = pencircle scaled 1.7bp;

	samplecolor = 0.88white;
	magnetcolor = 0.37white;
	samplebordercolor = 0.1white;
%	samplecolor = (0.85,0.9,0.8);
%	magnetcolor = (0.2,0.27,0.58);
%	samplebordercolor = blue;
	salto = 1mm;
	rotang = 7;
	posang = -28;
	effang = -20;
	shif = 50salto;
	thi = 5salto;
	pol = 25salto;
	pro = 8salto;
	wid = 58salto;
	cfac = 0.5;
	lblfac = 1.5;
	innerray = 27salto;
	circleray = 0.93innerray;
	squareray = 0.93circleray;
	glasssize = 0.42squareray;
	outerray = 60salto;
	midleray = 55salto;
	yray = 53salto;
	thebiray = 145salto;
	otherray = 115salto;
	wavecray = 0.5*(thebiray+otherray);
	effraya = 0.45*(midleray+otherray);
	effrayb = 0.56*(midleray+otherray);
	innercirc = fullcircle scaled innerray;
	outercirc = fullcircle scaled outerray;
	midlecirc = (halfcircle cutbefore (origin--dir(90+posang))) 
		    scaled midleray;
	midmidcir = 0.5*length midlecirc;
	midmidcirpos = point midmidcir of midlecirc;
	effcirc = ((halfcircle cutbefore (origin--dir(90+effang)))
	            cutafter (origin--up))
		    scaled effraya;
	mideffcir = 0.5*length effcirc;
	mideffcirpos = point mideffcir of effcirc;
	otherarch = ((halfcircle cutbefore (origin--dir(90+posang)))
		    cutafter (origin--up))
		    scaled otherray;
	midotharc = 0.5*length otherarch;
	midotharcpos = point midotharc of otherarch;
	wavecarch = ((halfcircle cutafter (origin--dir(90+rotang)))
		    cutbefore (origin--up))
		    scaled wavecray;
	midwavarc = 0.5*length wavecarch;
	midwavarcpos = point midwavarc of wavecarch;
	oridirh = origin--0.5*outerray*left;
	oridirv = origin--0.5*thebiray*up;
	oridirn = origin--0.5*thebiray*dir(90+rotang);
	oridirp = origin--0.5*wavecray*dir(90+posang);
	oridiry = origin--0.5*yray*dir(180+posang);
	lastoridire = 0.5*effrayb*dir(90+effang);
	oridire = origin--lastoridire;
	dirhoriz = oridirh cutbefore innercirc;
	dirverti = oridirv cutbefore innercirc;
	dirnegat = oridirn cutbefore innercirc;
	dirposit = oridirp cutbefore innercirc;
	direff = oridire cutbefore innercirc;
	diry = oridiry cutbefore innercirc;
	field = (0.5innerray,0)--(0.5wid-salto,0);
	loricir = halfcircle scaled circleray rotated 90;
	roricir = halfcircle scaled circleray rotated 270;
        baserec = (unitsquare shifted (-0.5*(1,1))) 
                  xscaled 2
                  scaled squareray;
	glassre = (unitsquare shifted (-0.5*(1,1))) 
		  yscaled (13.0/6.0)
		  scaled glasssize;
	lowcut = subpath (0,1) of baserec;
	higcut = subpath (2,3) of baserec;
	sampleh = buildcycle( higcut, loricir, reverse lowcut, roricir );
	draw glassre rotated 90 dashed evenly withcolor (cfac*white);
%	draw glassre rotated rotang dashed evenly withcolor (cfac*(white-blue));
	fill glassre rotated posang withcolor samplecolor;
	draw glassre rotated posang withcolor samplebordercolor;
	draw sampleh rotated posang withpen thickpen withcolor blue;
	drawarrow reverse midlecirc; 
	drawarrow reverse otherarch; 
	drawarrow wavecarch;
	drawarrow reverse effcirc;
	draw dirhoriz dashed evenly scaled 2;
	draw dirverti dashed evenly scaled 2;
	drawarrow dirnegat withpen thickpen;
	drawarrow dirposit withpen thickpen;
	drawarrow direff withpen thickpen;
	drawarrow reverse field withpen thickpen;
	drawarrow direff withpen thickpen;
	drawarrow diry withpen thickpen;
	lblpiplu = btex $\alpha = \frac{\pi}{2}+\beta$ etex scaled lblfac;
	lblalpha = btex $\beta$ etex scaled lblfac;
	lblwavec = btex $\Psi$ etex scaled lblfac;
	lblqu = btex $\vec{q}\,$ etex scaled lblfac;
	lblno = btex $\vec{n}_0$ etex scaled lblfac;
	lblfield = btex $\vec{B}$ etex scaled lblfac;
	lblx = btex $x$ etex scaled lblfac;
	lbly = btex $y$ etex scaled lblfac;
	lbleff = btex $\beta_{\mbox{\small eff}}$ etex scaled lblfac;
	label.ulft( lblpiplu, midmidcirpos );
	label.top( lblalpha, midotharcpos );
	label.top( lblwavec, midwavarcpos );
	label.llft( lblqu, point 0.7 of dirnegat );
	label.lrt( lblx, point 0.7 of dirposit );
	label.top( lblfield, point 0.5 of field );
	label.llft( lbly, point 0.5 of diry );
	label.top( lbleff, mideffcirpos );
	label.top( lblno, lastoridire );
	hei = pol+4thi;
	mar = (hei-2*thi-pol)*0.5;
	z100 = (pro,mar);
	z0 = (wid,pol)*0.5;
	z1 = z0 + z100;
	z2 = z1 + (thi,0);
	z3 = (xpart z2, pol*0.5+mar-hei);
	z4 = (-xpart z2, ypart z3);
	z5 = (-xpart z2, ypart z2);
	z6 = (-xpart z1, ypart z1);
	z7 = (-xpart z0, ypart z0);
	z8 = (xpart z7, -pol*0.5);
	z9 = z8 - z100;
	z10 = z9 - (0,thi);
	z11= (xpart z1, ypart z10);
	z12= z11+ (0,thi);
	z13= (xpart z0, -ypart z0);
	themagnet = z0--z1--z2--z3--z4--z5--z6--z7--
		    z8--z9--z10--z11--z12--z13--cycle;
	fill themagnet withcolor magnetcolor;
	draw themagnet withpen thickpen;
endfig;

verbatimtex 
\end{document} 
etex

end;
