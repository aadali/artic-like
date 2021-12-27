% integratorfigures.mp
% L. Nobre G.
% 2007

input featpost3Dplus2D;

verbatimtex
 \documentclass{article} 
 \begin{document} 
etex

beginfig(1);
  f := (5,3,2);
  Spread := 90;
  numeric axlen, disclen, discpos, dradius, dxo, dy;
  color dcenter, xlim, dxr, dxco, dxro, vdy, dxa, dxb;
  path xtip, neadx, fardx, dx;
  axlen = 1.65;
  disclen = 0.15;
  discpos = 0.85;
  dy = 0.02;
  vdy = -dy*green;
  dradius = 1 +-+ discpos;
  dxo = dradius + 0.5;
  cartaxes( axlen, axlen, axlen );
  dcenter = discpos*red;
  xlim = axlen*red;
  xtip = rp(dcenter)--rp(xlim);
  drawoptions( withpen pencircle scaled 1.2pt );
  rigorousdisc( 0, true, discpos*red, dradius, -disclen*red );
  drawoptions();
%  draw goodcirclepath( black, blue, 1 ) dashed evenly;
%  draw goodcirclepath( black, red, 1 ) dashed evenly;
%  draw goodcirclepath( black, green, 1 ) dashed evenly;
%  drawoptions( withpen pencircle scaled 1.2pt );
  draw spatialhalfcircle( black, blue, 1, true );
  draw spatialhalfcircle( black, red, 1, true );
  draw spatialhalfcircle( black, green, 1, true );
  drawoptions();
  drawarrow xtip;
%  draw rp(black)--rp(dcenter) dashed evenly;
  draw rigorousfearpath( black, 1 );
  dxco = (discpos,dxo,0);
  dxr = (discpos-disclen,dradius,0);
  dxro = (discpos-disclen,dxo,0);
  neadx = rp(dcenter)--rp(dxco);
  fardx = rp(dxr)--rp(dxro);
  dxa = dxco+vdy;
  dxb = dxro+vdy;
  draw fardx;
  draw neadx;
  dx = rp(dxa)--rp(dxb);
  drawdblarrow dx;
%  label.rt( btex $dx$ etex, rp(0.5[dxa,dxb]) );
  label.rt( "dx", rp(0.5[dxa,dxb]) );
  label.bot( "x" ,rp(xlim));
endfig;

input boxes;

beginfig(2);
  f := (3,2.05,2.35);
  Spread := 140;
  numeric labscal, bray, sray, angb, phi, phib;
  picture labrb, labrs, plb, pls;
  color cenb, momb, lpos, hpos;
  path vaa, vb, ha, hb, hapart;
  labscal = 2;
  ahlength := 8pt;
  ahangle := 30;
  bray = 0.95;
  sray = 0.4;
  angb = 9;
  phi = 55;
  phib = 83;
  labrb= btex $(b+r\cos\phi)d\theta$ etex scaled labscal;
  labrs= btex $r\,d\phi$ etex scaled labscal;
  drawoptions( withpen pencircle scaled 1.7pt withcolor 0.65blue );
  smoothtorus( black, blue, bray, sray );
%  drawoptions( withpen pencircle scaled 1.3pt dashed evenly );
  drawoptions( );
  vaa = goodcirclepath( bray*red, green, sray );
  draw vaa withcolor red;
  cenb = bray*(red*cosd(angb)+green*sind(angb));
  momb = green*cosd(angb)-red*sind(angb);
  lpos = cenb*(1+sray*cosd(phi)/bray)+blue*sray*sind(phi);
  hpos = cenb*(1+sray*cosd(phib)/bray)+blue*sray*sind(phib);
  vb = goodcirclepath( cenb, momb, sray );
  draw vb withcolor (blue+green);
  ha = goodcirclepath( blue*sray*sind(phi), blue, bray+sray*cosd(phi) );
  draw ha;
  hb = goodcirclepath( blue*sray*sind(phib), blue, bray+sray*cosd(phib) );
  draw hb withcolor green;
%  show length ha;
  hapart = subpath (3,6) of ha;
  drawarrow hapart withpen pencircle scaled 1.3pt;
  fill buildcycle( vaa, hb, vb, reverse hapart ) withcolor 0.5white;
  plb = thelabel.llft( labrb, rp(lpos) );
  pls = thelabel.lrt( labrs, rp(hpos) );
  boxit.bb( plb );
  boxit.bs( pls );
  unfill bpath bb;
  unfill bpath bs;
  drawboxed( bb, bs );
endfig;
  
verbatimtex 
 \end{document} 
etex

end;
  