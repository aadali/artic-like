% parafusoprojxy.mp
% L. Nobre G.
% 2009

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

f := (3,5,15);
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
  defaultscale := 0.65;
  label.rt( "X", rp( xend ) );
  label.lft( "Y", rp( yend ) );
  emptyline( false, 1, black, 0.7background, 700, 0.8, 2, thexyline );
  drawarrow rp( zbeg )--rp( zend );
  label.bot( "Z", rp( zend ) );
endfig;

end;



