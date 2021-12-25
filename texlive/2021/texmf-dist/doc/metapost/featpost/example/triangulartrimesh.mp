% triangulartrimesh.mp
% L. Nobre G.
% IYP (2005)

input featpost3Dplus2D;

def zsurface( expr xc, yc ) =
	cosd(xc*57)*cosd(yc*57)+4*mexp(-(xc**2+yc**2)*6.4)
enddef;

f := 16*(4,1,1);
Spread := 35;
LightSource := 10*(4,-3,4);
defaultfont := "cmvtt10";
beginfig(1); 
  numeric np, ssize, lx, ly, topl, ystep, i, zl;
  string L[];
  L0 = "P";
  L1 = "H";
  L2 = "Y";
  L3 = "S";
  L4 = "I";
  L5 = "C";
  L6 = "S";
  np = 72;
  ssize = 15;
  ly = -0.5*ssize;
  lx = -sqrt(3)*ssize/6;
  ystep = ssize/7;
  topl = 2*zsurface( 0, 0 );
  for i=0 upto 6:
    zl := zsurface(lx,i*ystep+ly);
    labelinspace(false,(lx,i*ystep+ly,zl),green*ystep,blue*(topl-zl),L[i]);
  endfor;
  hexagonaltrimesh( false,np,ssize,zsurface);
  label.ulft("matagalatlante.org",lrcorner currentpicture); 
endfig;

end;
