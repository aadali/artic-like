% hiddensurfaces.mp
% L. Nobre G.
% 2003

input featpost3Dplus2D;

def zsurface( expr xc, yc ) = 
	begingroup
		numeric zout;
		zout = yc**4-xc*yc*(1-yc**2);
		if zout > 1.5:
			zout := 1.5;
		fi;
		( zout )
	endgroup 
enddef;

f := (2.5,2,1.25);
Spread := 127;

beginfig(1); 
	numeric i, j, nx, ny, lowx, higx, lowy, higy, refz
	        posx, posy, posz, facz, counter;
	color poi[][];
	nx = 30;
	ny = 30;
	lowy = -2.5;
	higy = 1;
	lowx = 1.5;
	higx = 3;
	refz = 0;
	facz = 5;
	counter = 0;
	for i=0 upto nx:
	    for j=0 upto ny:
		posx := lowx + i*( higx - lowx )/nx;
		posy := lowy + j*( higy - lowy )/ny;
		posz := zsurface( posx, posy );
		posx := bracket(lowx,posx,higx);
		posy := bracket(lowy,posy,higy);
		posz := posz/facz + refz;
		poi[i][j] := ( posx, posy, posz );
	    endfor;
	endfor;
	for i=1 upto nx:
	    for j=1 upto ny:
		counter := incr(counter);
		npf[counter] := 3;
		F[counter]p1 := poi[i-1][j-1];
		F[counter]p2 := poi[i][j-1];
		F[counter]p3 := poi[i][j];
		counter := incr(counter);
		npf[counter] := 3;
		F[counter]p1 := poi[i-1][j-1];
		F[counter]p2 := poi[i-1][j];
		F[counter]p3 := poi[i][j];
	    endfor;
	endfor;
	NF := counter;
	drawarrow rp((-1,-1,refz))--rp((1,-1,refz));
	drawarrow rp((-1,-1,refz))--rp((-1,1,refz));
	fill_faces( withcolor blue );
endfig;

end;
