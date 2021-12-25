% hexagonaltrimesh.mp
% L. Nobre G.
% 2003

input featpost3Dplus2D;

def zsurface( expr xc, yc ) =
	cosd(xc*57)*cosd(yc*57)+4*mexp(-(xc**2+yc**2)*6.4)
enddef;

f := 7*(4,1,5);
Spread := 35;
LightSource := 10*(4,-3,4);
SubColor := 0.4background;

beginfig(1); 
	numeric np, ssize;
	path chair;
	np = 40;
	ssize = 15;

	hexagonaltrimesh( true,np,ssize,zsurface);
endfig;

end;
