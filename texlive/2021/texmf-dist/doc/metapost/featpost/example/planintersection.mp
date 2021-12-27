% planintersection.mp
% L. Nobre G. 
% 2004

input featpost3Dplus2D;

verbatimtex \documentclass{article}\begin{document}{\sffamily etex

beginfig(1);
  f := 0.5*(10,14,5.5);
  color a, x, b, c, D, y, w, g, h;
  color auxa, auxb, auxc, auxd, auxe, auxf;
  pair appear;
  numeric wid, hei, alfa, beta, forw, back, labscal;
  pen thedotpen, borderpen, normalpen;
  picture dasha, dashb, labx, laba, labb, labc, labD;
  path outerb, innerb, apath, bpath;
  wid = 1.5;
  hei = 1;
  forw = 1.5;
  back = 1.5;
  alfa = -10;
  beta = -30;
  labscal = 2.8;
  labeloffset := 8pt;
  ahangle := 30;
  ahlength := 10pt;
  a = ( 0.7*forw, 0.15*wid, 0.3*hei );
  b = ( -0.7*back, 0.5*wid, 0.8*hei );
  c = ( 0, 0.7*wid, 0.5*hei );
  D = ( cosd( alfa )*cosd( beta ), sind( alfa )*cosd( beta ), sind( beta ) );
  labx = btex $\vec{x}$ etex scaled labscal;
  laba = btex $\vec{a}$ etex scaled labscal;
  labb = btex $\vec{b}$ etex scaled labscal;
  labc = btex $\vec{c}$ etex scaled labscal;
  labD = btex $\vec{D}$ etex scaled labscal;
  thedotpen = pencircle scaled 4pt;
  borderpen = pencircle scaled 2pt;
  normalpen = pencircle scaled 1.5pt;
  dasha = evenly scaled 2;
%  dashb = dashpattern( on 7pt off 5pt on 7pt );
  dashb = withdots;
  apath = rp((X(a),Y(a),0))--rp(a)--rp((X(a),0,Z(a)));
  bpath = rp((X(b),Y(b),0))--rp(b)--rp((X(b),0,Z(b)));
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
  label.top( labb, rp(b) );
  
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
  label.top( laba, rp(a) );

  x = whatever[a,b];
  cdotprod( x-c, D ) = 0;
  draw rp(x) withpen thedotpen;
  label.top( labx, rp(x) );

  draw rp(a)--rp(x) withpen normalpen;
  draw rp(x)--rp(b) withpen normalpen dashed dashb;

  draw rp(c) withpen thedotpen;
  drawarrow rp(c)--rp(c+0.7*D) withpen borderpen;
  label.top( labc, rp(c) );
  
  label.bot( labD, rp(c+0.5*D) );
  
  draw rp(x)--rp(c) dashed dasha withpen normalpen;
  drawoptions( withpen normalpen );
  squareangline( x, c+D, c, 0.08*wid );
  
endfig;

verbatimtex }\end{document} etex

end;
