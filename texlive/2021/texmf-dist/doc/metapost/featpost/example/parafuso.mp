% parafuso.mp
% L. Nobre G.
% 2005

input featpost3Dplus2D;

def theline( expr TheVal ) = 
  begingroup
    numeric par[], cred, cgre, cblu, param;
    par0 = 0;
    par1 = 6*360;
    par2 = 0.3;
    par3 = 850;
    param = par0+TheVal*(par1-par0);
    cred = -par2*cosd( param );
    cblu = par2*sind( param );
    cgre = param/par3;
    ( (cred,cgre,cblu) )
  endgroup 
enddef;

def thexyline( expr TheVal ) = 
  begingroup
    numeric par[], cred, cgre, param;
    par0 = 0;
    par1 = 6*360;
    par2 = 0.3;
    par3 = 850;
    param = par0+TheVal*(par1-par0);
    cred = -par2*cosd( param );
    cgre = param/par3;
    ( (cred,cgre,0) )
  endgroup 
enddef;

def loxodrome( expr TheVal ) = 
  begingroup
    numeric cred, cgre, cblu, alpha, theta, ratio, radio;
    ratio = 22;
    radio = 2;
    theta = ratio*180*TheVal;
    alpha = 170*TheVal-85;
    cred = radio*cosd(theta)*cosd(alpha);
    cgre = radio*sind(theta)*cosd(alpha);
    cblu = radio*sind(alpha);
    ( (cred,cgre,cblu) )
  endgroup 
enddef;

f := (2,4,1.8);
Spread := 77;
ShadowOn := false;

beginfig(1); 
  numeric axsize, zaxpos, zaxlen;
  color xbeg, xend, ybeg, yend, zbeg, zend;
  axsize = 0.85;
  zaxpos = 0.55;
  zaxlen = 2.1;
  pickup pencircle scaled 1.5pt;
  xbeg = (axsize,0,0);
  xend = (-axsize,0,0);
  ybeg = (0,0,-axsize);
  yend = (0,0,axsize);
  zbeg = (zaxpos,-zaxpos,0);
  zend = (zaxpos,zaxlen,0);
  drawarrow rp( xbeg )--rp( xend );
  drawarrow rp( ybeg )--rp( yend );
%	defaultfont := "putr8r";
  defaultscale := 1.95;
  label.rt( "A", rp( xend ) );
  label.lft( "B", rp( yend ) );
  emptyline( false, 1, black, 0.5background, 1000, 0.82, 2, theline );
  drawarrow rp( zbeg )--rp( zend );
  label.bot( "C", rp( zend ) );
endfig;

beginfig(2); 
  numeric axsize, zaxpos, zaxlen;
  color xbeg, xend, ybeg, yend, zbeg, zend;
  axsize = 0.85;
  zaxpos = 0.55;
  zaxlen = 2.1;
  pickup pencircle scaled 1.5pt;
  xbeg = (axsize,0,0);
  xend = (-axsize,0,0);
  ybeg = (0,0,-axsize);
  yend = (0,0,axsize);
  zbeg = (zaxpos,-zaxpos,0);
  zend = (zaxpos,zaxlen,0);
  drawarrow rp( xbeg )--rp( xend );
  drawarrow rp( ybeg )--rp( yend );
%	defaultfont := "putr8r";
  defaultscale := 0.65;
  label.rt( "X", rp( xend ) );
  label.lft( "Y", rp( yend ) );
  closedline( false, 70, 0.5, 0.95, theline );
  drawarrow rp( zbeg )--rp( zend );
  label.bot( "Z", rp( zend ) );
  doitnow;
endfig;

f := (2,4,2.8);
beginfig(3); 
  NCL := 0;
  Nobjects := 0;
  closedline( false, 270, 0.5, 0.85, loxodrome );
  doitnow;
endfig;

beginfig(4); 
  numeric axsize, zaxpos, zaxlen;
  color xbeg, xend, ybeg, yend, zbeg, zend;
  axsize = 0.85;
  zaxpos = 0.55;
  zaxlen = 2.1;
  pickup pencircle scaled 1.5pt;
  xbeg = (axsize,0,0);
  xend = (-axsize,0,0);
  ybeg = (0,0,-axsize);
  yend = (0,0,axsize);
  zbeg = (zaxpos,-zaxpos,0);
  zend = (zaxpos,zaxlen,0);
  drawarrow rp( xbeg )--rp( xend );
  drawarrow rp( ybeg )--rp( yend );
%	defaultfont := "putr8r";
  defaultscale := 0.65;
  label.rt( "X", rp( xend ) );
  label.lft( "Y", rp( yend ) );
  closedline( false, 70, 0.5, 0.95, thexyline );
  closedline( false, 70, 0.5, 0.95, theline );
  drawarrow rp( zbeg )--rp( zend );
  label.bot( "Z", rp( zend ) );
  doitnow;
endfig;

end;



