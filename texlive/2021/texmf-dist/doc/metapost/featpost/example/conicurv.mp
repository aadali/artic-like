% conicurv.mp 
% L. Nobre G.
% Y2K, 2008

input featpost3Dplus2D;

verbatimtex 
\documentclass{article} 
\usepackage{beton} 
\usepackage{concmath} 
\usepackage{ccfonts} 
\begin{document} 
etex

beginfig(1);
	f := (10,-5,5.44);
	Spread := 30;
	numeric shortradius, longradius, theta, width, height;
	numeric updiff, refsize, vecsize, anglar, bord, totup;
	color rn, fa, fc, pg, centre, speed, refx, refy, refz, refY, refX;
	color iplow, iphig, oplow,ophig, anglebase, central, refo;
	color eplow, ephig, aplow,aphig, anglebase, central, refo;
	theta = 30;
	width = 3;
	shortradius = 2;
	bord = 2;
	refsize = 1;
	vecsize = 2;
	height = 0.3;
	anglar = 1.75;
	totup = 3;
	longradius = shortradius + width*cosd(theta);
	updiff = width*sind(theta);
	iplow = (0,shortradius,0);
	iphig = (0,longradius,updiff);
	oplow = (-shortradius,0,0);
	ophig = (-longradius,0,updiff);
	aplow = -iplow;
	aphig = (0,-longradius,updiff);
	eplow = -oplow;
	ephig = (longradius,0,updiff);
	anglebase = (0,longradius,0);
	centre = 0.5[iplow,iphig]+(0,0,height);
	central = (0,0,Z(centre));
%	refo = (0,0,-shortradius*sind(theta)/cosd(theta));
	refo = (0,0.5*Y(centre),Z(centre));
	refx = refsize*(0,cosd(theta),sind(theta));
	refy = refsize*(0,-sind(theta),cosd(theta));
	refz = refsize*(1,0,0);
	refX = refsize*green;
	refY = refsize*blue;
%	anglinen( iplow, oplow, (0,0,0), shortradius, "", 0 );
%	anglinen( iphig, ophig, (0,0,updiff), longradius, "", 0 );
	angline(iphig,anglebase,iplow,width/anglar,btex $\theta$ etex,rt);
	draw rp(central) withpen pencircle scaled 2pt;
	draw rp(iplow)--rp(iphig);
	draw rp(oplow)--rp(ophig);
	draw rp(aplow)--rp(aphig);
	draw rp(eplow)--rp(ephig);
	draw rp(iphig)--rp(anglebase)--rp(aplow) dashed evenly;
	draw rp(oplow)--rp(eplow) dashed evenly;
	draw rp(central)--rp(centre) dashed evenly;
%	draw rp((0,-bord,-bord))--rp((0,longradius+bord,-bord))
%	   --rp((0,longradius+bord,totup))--rp((0,-bord,totup))--cycle
%	     withpen pencircle scaled 2pt;
	draw rp((0,0,-bord))--rp((0,longradius+bord,-bord))
	   --rp((0,longradius+bord,totup))--rp((0,0,totup))--cycle;
	drawarrow rp(refo)..rp(refo+refy);
	drawarrow rp(refo)..rp(refo+refx);
	drawarrow rp(black)..rp(refX);
	drawarrow rp(black)..rp(refY);
	drawarrow rp(refo)..rp(refo+refz);
	label.llft(btex $y$ etex, rp(refo+refy));
	label.lrt(btex $x$ etex, rp(refo+refx));
	label.lft(btex $Y$ etex, rp(refY));
	label.bot(btex $X$ etex, rp(refX));
	label.llft(btex $z$ etex, rp(refo+refz));
	rn = centre+vecsize*refy;
	fa = centre+vecsize*refx;
	fc = centre+vecsize*(0,1,0);
	pg = centre+vecsize*(0,0,-1);
	speed = centre+vecsize*(1,0,0);
	drawarrow rp(centre)..rp(rn);
	drawarrow rp(centre)..rp(fa);
	drawarrow rp(centre)..rp(fc);
	drawarrow rp(centre)..rp(pg);
	drawarrow rp(centre)..rp(speed);
	label.rt(btex $\vec{R}_N$ etex, rp(rn));
	label.top(btex $\vec{F}_a$ etex, rp(fa));
	label.urt(btex $\vec{F}_c$ etex, rp(fc));
	label.rt(btex $\vec{P}$ etex, rp(pg));
	label.rt(btex $\vec{v}$ etex, rp(speed));
	draw rp(centre) withpen pencircle scaled 10pt withcolor blue;
	pickup pencircle scaled 2pt;
	draw rigorouscircle( (0,0,updiff), (0,0,1), longradius );
	draw rigorouscircle( (0,0,0), (0,0,1), shortradius );
endfig;

verbatimtex 
\end{document} 
etex

end;
