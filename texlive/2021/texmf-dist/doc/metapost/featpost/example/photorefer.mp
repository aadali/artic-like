% photorefer.mp
% L. Nobre G.
% 2006

input featpost3Dplus2D;

%\usepackage{beton,concmath,ccfonts} 

verbatimtex \documentclass{article} \begin{document} etex

f := ( 6, 5, 2.5 );
Spread := 160;

beginfig(2);
  color a, x, b, c, D, y, w, g, h;
  color auxa, auxb, auxc, auxd, auxe, auxf, azz, azy, bzz, bzy;
  pair appear, rpc, rpch, rpcv;
  numeric wid, hei, alfa, beta, forw, back, labscal, framang, lx, ly, verscal;
  pen thedotpen, borderpen, normalpen;
  picture dasha, dashb, labx, laba, labb, labc, labD, labh, labv;
  path outerb, innerb, apath, bpath;
  framang = 0;
  alfa = 0;
  beta = 0;
%  framang = -15;
%  alfa = -20;
%  beta = -15;

  wid = 1.5;
  hei = 1;
  forw = 1.5;
  back = 1.5;
  lx = 1.15;
  ly = 0.4;
  verscal = 0.47;
  labscal = 2.8;
  labeloffset := 8pt;
  ahangle := 30;
  ahlength := 10pt;
  b = ( -0.8*back, 0.25*wid, 0.9*hei );
  c = ( 0, 0.85*wid, 0.35*hei );
  D = ( cosd( alfa )*cosd( beta ), sind( alfa )*cosd( beta ), sind( beta ) );
  a = c+0.95*D;
%  a = c+0.7*D;
  labx = btex $\vec{p}$ etex scaled labscal;
  laba = btex $\vec{f}$ etex scaled labscal;
  labb = btex $\vec{R}$ etex scaled labscal;
  labc = btex $\vec{o}$ etex scaled labscal;
  labh = btex $\hat{h}$ etex scaled labscal;
  labv = btex $\hat{v}$ etex scaled labscal;  
  thedotpen = pencircle scaled 4pt;
  borderpen = pencircle scaled 2pt;
  normalpen = pencircle scaled 1.5pt;
  dasha = evenly scaled 2;
%  dashb = dashpattern( on 7pt off 5pt on 7pt );
  dashb = withdots;
  azz = (X(a),Y(a),0);
  azy = (X(a),0,Z(a));
  bzz = (X(b),Y(b),0);
  bzy = (X(b),0,Z(b));
  apath = rp(azz)--rp(a)--rp(azy);
  bpath = rp(bzz)--rp(b)--rp(bzy);
  auxa = ( forw, 0, 0 );
  auxb = ( forw, 0, hei );
  auxc = ( -back, 0, hei );
  auxd = ( -back, 0, 0 );
  auxe = ( -back, wid, 0 );
  auxf = ( forw, wid, 0 );
  outerb = rp(auxa)--rp(auxb)--rp(auxc)--rp(auxd)--rp(auxe)--rp(auxf)--cycle;
  draw outerb withpen borderpen;
  draw rp(auxa)--rp(auxd) withpen borderpen;

  draw bpath dashed dasha withpen normalpen;
  draw rp(b) withpen thedotpen;
  label.rt( labb, rp(b) );
  
  draw rp(a)--rp(b) withpen normalpen;
  
  cdotprod( y-c, D ) = 0;
  cdotprod( y-black, green ) = 0;
  cdotprod( y-black, blue ) = 0;
  
  cdotprod( w-c, D ) = 0;
  cdotprod( w-black, green ) = 0;
  cdotprod( w-hei*blue, blue ) = 0;

  cdotprod( g-c, D ) = 0;
  cdotprod( g-wid*green, green ) = 0;
  cdotprod( g-hei*blue, blue ) = 0;

  cdotprod( h-c, D ) = 0;
  cdotprod( h-wid*green, green ) = 0;
  cdotprod( h-black, blue ) = 0;

  innerb = rp(y)--rp(w)--rp(g)--rp(h)--cycle;
  unfill innerb;
  draw innerb withpen borderpen;

  draw apath dashed dasha withpen normalpen;
  draw rp(a) withpen thedotpen;
  label.llft( laba, rp(a) );

  x = whatever[a,b];
  cdotprod( x-c, D ) = 0;
  draw rp(x) withpen thedotpen;
  label.lft( labx, rp(x) );

  draw rp(a)--rp(x) withpen normalpen;
  draw rp(x)--rp(b) withpen normalpen dashed dashb;

  draw rp(c) withpen thedotpen;
  draw rp(a)--rp(c) withpen borderpen;
  label.bot( labc, rp(c) );

  numeric theang;
  color perpo, perpa, perph, perpv;
  perpo = N( x-c );
  perpa = ncrossprod( D, c-x );
  theang = -getangle( w-y, perpa );
  perph = perpo*cosd( framang+theang ) + perpa*sind( framang+theang );
  perpv =-perpo*sind( framang+theang ) + perpa*cosd( framang+theang );
  color corn, corm;
  numeric coox, cooy;
  coox = cdotprod( x-c, perph );
  cooy = cdotprod( x-c, perpv );
  draw rp(c)--rp(c+lx*perph)--rp(c+lx*perph+ly*perpv)
            --rp(c+ly*perpv)--cycle dashed evenly;
  draw rp(c)--rp(c+coox*perph)--rp(c+coox*perph+cooy*perpv)
            --rp(c+cooy*perpv)--cycle dashed evenly;  
  rpc = rp(c);
  rpch = rp(c+verscal*perph);
  rpcv = rp(c+verscal*perpv);
  drawarrow rpc--rpch withpen normalpen;
  drawarrow rpc--rpcv withpen normalpen;
  label.rt( labv, rpcv );
%  label.top( labh, rpch );
  label.lrt( labh, rpch );
endfig;

f := ( 5, 3, 1 );
Spread := 70;

beginfig(1);
  labeloffset := 3pt;
  ahangle := 30;
  ahlength := 5pt;
  numeric alpha, beta, gamma, axeslen, blen, angray;
  alpha = 65;
  beta  = 70;
  gamma = 75;
  blen  = 1.4;
  axeslen = 1.2;
  angray = 0.4;
  orthaxes( axeslen, btex $E$ etex, axeslen, btex $N$ etex, axeslen, btex $Z$ etex );
  path vo, vh, vv, vhp, vop;
  z0 = rp( black );
  color po, ph, pv, php, pop;
  pop = ( cosd(alpha), sind(alpha), 0 );
  po = blen*( cosd(beta)*pop + sind(beta)*blue );
  php = ( sind(alpha), -cosd(alpha), 0 );
  numeric hx, hy, hv, vx, vy, vz;
  hx = sind(alpha)*cosd(gamma)-sind(beta)*cosd(alpha)*sind(gamma);
  hy =-cosd(alpha)*cosd(gamma)-sind(beta)*sind(alpha)*sind(gamma);
  hz = cosd(beta)*sind(gamma);
  vx =-sind(alpha)*sind(gamma)-sind(beta)*cosd(alpha)*cosd(gamma);
  vy = cosd(alpha)*sind(gamma)-sind(beta)*sind(alpha)*cosd(gamma);
  vz = cosd(beta)*cosd(gamma);
  ph = ( hx, hy, hz );
  pv = ( vx, vy, vz );
  z1 = rp( po );
  z2 = rp( ph );
  z3 = rp( pv );
  z4 = rp( php );
  z5 = rp( pop );
  vo = z0--z1;
  vh = z0--z2;
  vv = z0--z3;
  vhp = z0--z4;
  vop = z0--z5;
  drawarrow vo;
  drawarrow vh;
  drawarrow vv;
  draw vhp dashed evenly;
  draw vop dashed evenly;
  squareangline( ph, pv, black, 0.1 );
  angline( pop, po, black, angray, btex $\beta$ etex, rt );
  angline( php, ph, black, angray, btex $\gamma$ etex, lft );
  angline( red, pop, black, angray, btex $\alpha$ etex, bot );
  label.top( btex $\vec{f}-\vec{o}$ etex, z1 );
  label.top( btex $\hat{h}$ etex, z2 );
  label.top( btex $\hat{v}$ etex, z3 );
endfig;

verbatimtex 
\end{document} 
etex

end.
