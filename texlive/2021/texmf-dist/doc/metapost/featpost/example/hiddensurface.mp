% hiddensurface.mp
% L. Nobre G.
% 2004

input featpost3Dplus2D;

% The big nice trick to draw 3D hidden line surfaces

verbatimtex 
\documentclass{article} 
\usepackage{beton} 
\usepackage{concmath} 
\usepackage{ccfonts} 
\begin{document} 
etex

numeric toplevel;
toplevel = 2.2;

def zsurface( expr xc, yc ) = 
	begingroup
		numeric zout;
		zout = yc**4 -xc*yc*(1-yc**2);
		if zout > toplevel:
			zout := toplevel;
		fi;
		( zout )
	endgroup 
enddef;

def zline( expr yc ) = 
	begingroup
		numeric zout;
		zout = 4*yc*yc*yc/(1-3*yc*yc); %% I don't get it...
		if zout > 3:
			zout := 3;
		fi;
		( zout )
	endgroup 
enddef;

f := (6,1,2);
Spread := 67;

beginfig(1); 
	numeric i, j, nx, ny, lowx, higx, lowy, higy, refz;
	numeric posx, posy, posz, facz, lowl, higl, pensize;
	color gridx[], pos, begcont, endcont, begconf, endconf;
	path facex, descontjump;
	pensize = 3pt;
	nx = 20;
	ny = 50;
	lowy = -2.5;
	higy = 1;
	lowx = 1.5;
	higx = 3;
	refz = -2.5;
	facz = 5;
	for i=0 upto nx:
	    for j=0 upto ny:
		posx := lowx + i*( higx - lowx )/nx;
		posy := lowy + j*( higy - lowy )/ny;
		posz := zsurface( posx, posy );
		posx := bracket(lowx,posx,higx);
		posy := bracket(lowy,posy,higy);
		posz := posz/facz;
		gridx[j] := ( posx, posy, posz );
	    endfor;
	    pos := ( posx, -1, refz/facz );
	    facex := rp( pos )--rp( gridx0 );
	    for j=1 upto ny:
		facex := facex--rp( gridx[j] );
	    endfor;
	    facex:=facex--rp((posx,1,refz/facz))--cycle;
	    unfill facex;
	    draw facex;
	endfor;
	lowl = -2.0;
%	higl = -0.7;
	higl = -1.0;
	for j=0 upto ny:
		posy := lowl + j*( higl - lowl )/ny;
		posx := zline( posy );
		posz := zsurface( posx, posy );
		posx := bracket(lowx,posx,higx);
		posy := bracket(lowy,posy,higy);
		posz := posz/facz;
		gridx[j] := ( posx, posy, posz );
	endfor;
	facex := rp( gridx0 );
	for j=1 upto ny:
		facex := facex--rp( gridx[j] );
	endfor;
	draw facex withcolor red; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	lowl := -1.4;
	higl := -0.8;
	for j=0 upto ny:
		posy := lowl + j*( higl - lowl )/ny;
		posx := -2*posy;
		posz := zsurface( posx, posy );
		posx := bracket(lowx,posx,higx);
		posy := bracket(lowy,posy,higy);
		posz := posz/facz;
		gridx[j] := ( posx, posy, posz );
	endfor;
	facex := rp( gridx0 );
	for j=1 upto ny:
		facex := facex--rp( gridx[j] );
	endfor;
	draw facex withcolor green;
	lowl := -1.0;
	higl := 0.5;
	for j=0 upto 1:
		posy := lowl + j*( higl - lowl );
		posx := 2;
		posz := zsurface( posx, posy );
		posx := bracket(lowx,posx,higx);
		posy := bracket(lowy,posy,higy);
		posz := posz/facz;
		gridx[j] := ( posx, posy, posz );
	endfor;
%	facex := rp( gridx0 );
%	for j=1 upto ny:
%		facex := facex--rp( gridx[j] );
%	endfor;
%	draw facex withcolor blue;	
	pickup pencircle scaled pensize;
	draw rp(gridx0) withcolor red;
	draw rp(gridx1) withcolor red;	
	posx := 2.83;
	posy := -1.933;
	posz := zsurface( posx, posy );
	posx := bracket(lowx,posx,higx);
	posy := bracket(lowy,posy,higy);
	posz := posz/facz;
	draw rp((posx,posy,posz)) withcolor blue;
	posy := 0.517;
	posy := bracket(lowy,posy,higy);
	draw rp((posx,posy,posz)) withcolor blue;
	begconf = (1.15,-1,refz/facz);
	endconf = (1.15, 1,refz/facz);
	begcont = ( 1,-1,1.15*toplevel/facz);
	endcont = (-1,-1,1.15*toplevel/facz);
	pickup pencircle scaled 0.5pensize;
	drawarrow rp(begconf)--rp(endconf);
	drawarrow rp(begcont)--rp(endcont);
	label.ulft( btex $\Lambda$ etex, rp(0.5*(begcont+endcont)) );
	label.bot( btex $\Upsilon$ etex, rp(0.5*(begconf+endconf)) );
endfig;

verbatimtex \end{document} etex

end;
