% hiddenlinegraph.mp
% L. Nobre G.
% 2003
% The application of the big nice trick to draw 3D hidden line surfaces 

input featpost3Dplus2D;

verbatimtex 
\documentclass{article} 
\usepackage{beton} 
\usepackage{concmath} 
\usepackage{ccfonts} 
\begin{document} 
etex

def zsurfaci( expr xc, yc ) = 
	begingroup
		numeric zout;
		zout = yc**4-xc*yc*(1-yc**2);
		( ( xc, yc, zout ) )
	endgroup 
enddef;

def zsurface( expr xc, yc ) =
	( xc, yc, cosd(xc*57)*cosd(yc*57)+4*mexp(-(xc**2+yc**2)*6.4) )
enddef;

def zsurfaca( expr xc, yc ) =
	( xc, yc, 1.5-0.2*xc )
enddef;

%f := (6,1,0);
f := (5,1.5,3);
Spread := 147;
NF := 0;
LightSource := (1,-1,2);

beginfig(1); 
	numeric nx, ny, lowx, higx, lowy, higy, 
	        facz, lowz, higz, base, topo;
	path chair;
	nx = 22;
	ny = 22;
	lowy = -7;
	higy = 7;
	lowx = -7;
	higx = 7;
	facz = 2;
	lowz = -1;
	higz = 5;
%	lowy = -2.5;
%	higy = 1;
%	lowx = 1.5;
%	higx = 3;
%	facz = 3;
%	lowz = -1.5;
%	higz = 1.5;
	partrimesh( nx,ny,lowx,higx,lowy,higy,
	                  lowx,higx,lowy,higy,lowz,higz,facz,zsurface);
	partrimesh( nx,ny,lowx,higx,lowy,higy,
	                  lowx,higx,lowy,higy,lowz,higz,facz,zsurfaca);
	base = -1/facz;
	topo = 1/facz;
	chair = rp((-1,-1,base))--rp((1,-1,base))--rp((1,1,base))
	      --rp((-1,1,base))--rp((-1,-1,base))--rp((-1,-1,topo))
	      --rp((-1,1,topo))--rp((-1,1,base))--cycle;
	draw chair;
%	fill_faces( withcolor blue );
	draw_invisible( false, false, green, blue );
endfig;

verbatimtex \end{document} etex

end;
