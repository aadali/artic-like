% probtodooterr.mp
% L. Nobre G.
% 2003

verbatimtex 
\documentclass{article} 
\begin{document} 
etex

beginfig(1);
  ahangle := 30;
  ahlength := 5mm;
  currentpen := pencircle scaled 1mm;
	numeric u, ray, delta, theta, tray, vlen;
	numeric cray, clen, chig, cfra, cang;
	path cpat, fall, vpat, tpat, apat, rpat;
	picture cpic;
	u = 1mm;
	cray = 4u;
	clen = 27u;
	chig = 12u;
	cfra = 1.3;
	cang = 15;
	theta = 38;
	ray = 65u;
	delta = 12u;
	vlen = 39u;
	tray = 0.3ray;
	z1 = (-clen/2,0);
	z2 = (-clen/2+cray,chig);
	z3 = (cray,chig);
	z4 = (clen/2,1.5*cray);
	z5 = (clen/2,0);
	z6 = (-clen/2+cfra*cray,0);
	z7 = (clen/2-cfra*cray,0);
	z8 = (ray,0);
	z9 = (ray+delta,0);
	z10 = (0,ray);
	z11 = ray*(sind(theta),cosd(theta));
	z12 = z11+vlen*(cosd(theta),-sind(theta));	
	cpat = z1{up}..{right}z2--z3{right}..z4..{down}z5{left}...{left}z1;
%	dotlabels( 1,2,3,4,5);
	fall = z11{dir(-theta)}..{dir(-85)}z9;
	vpat = z11--z12;
	rpat = origin--z11;
	tpat = (halfcircle cutafter (origin--up)) scaled (2ray);
	apat = ((halfcircle cutafter (origin--up))
	  cutbefore rpat) scaled (2tray);	
	draw cpat;
	filldraw fullcircle scaled cray shifted z6;
	filldraw fullcircle scaled cray shifted z7;
	cpic = currentpicture;
	currentpicture := nullpicture;
	draw z9--tpat;
	draw fall dashed evenly;
	drawarrow vpat;
	drawarrow rpat;
	draw apat;
	draw cpic rotated (-cang)
	  shifted ((ray+0.5cray)*(sind(cang),cosd(cang)));
	draw origin--z10 dashed evenly;
	label.top(btex $\theta$ etex scaled 2, point 0.5 of apat );
	label.top(btex $l$ etex scaled 2, 0.5[z8,z9] );
	label.lrt(btex $R$ etex scaled 2, point 0.5 of rpat );
	label.urt(btex $\vec{v}_0$ etex scaled 2, point 0.5 of vpat );

endfig;

verbatimtex 
\end{document} 
etex

end;


