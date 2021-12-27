% revolvers.mp
% L. Nobre G. 
% 2012

input featpost3Dplus2D;

verbatimtex \documentclass{article}\begin{document}{\sffamily etex

def posang( expr Angul ) =
  rp( green*cosd( Angul ) + blue*sind( Angul ) )
enddef;
  
beginfig(1);
  f := (3,5,4);
  numeric labscal;
  picture labr, labc, laba, labum, labdo, labtr, labqu, labci;
  labscal = 2;
  ahlength := 8pt;
  ahangle := 30;
  labc = btex $\vec{c}$ etex scaled labscal;
  laba = btex $\vec{A}$ etex scaled labscal;
  labr = btex $R$ etex scaled labscal;
  labum= btex 1 etex scaled labscal;
  labdo= btex 2 etex scaled labscal;
  labtr= btex 3 etex scaled labscal;
  labqu= btex 4 etex scaled labscal;
  labci= btex 5 etex scaled labscal;
  drawoptions( withpen pencircle scaled 1.7pt withcolor 0.65blue );
  tropicalglobe( 5, black, 1, blue );
  drawoptions( withpen pencircle scaled 1.3pt );
  drawdblarrow rp(black)--rp(red);
  label.top( labr, rp((0.5,0,0)) );
  drawarrow rp(black)--rp(blue);
  label.rt( laba, rp((0,0,0.5)) );
  dotlabel.lrt( labc, rp(black) );
  dotlabel.rt( labum, posang(-72) );
  dotlabel.top( labdo, posang(-36) );
  dotlabel.top( labtr, posang(0) );
  dotlabel.top( labqu, posang(36) );
  dotlabel.top( labci, posang(72) );
endfig;

beginfig(2);
  f := (3,5,5.7);
  numeric labscal;
  picture labri, labc, laba, labro;
  labscal = 2;
  ahlength := 8pt;
  ahangle := 30;
  labc = btex $\vec{c}$ etex scaled labscal;
  laba = btex $\vec{A}$ etex scaled labscal;
  labri= btex $R_i$ etex scaled labscal;
  labro= btex $R_o$ etex scaled labscal;
  drawoptions( withpen pencircle scaled 1.7pt withcolor 0.65blue );
  rigorousdisc( 0.5, true, black, 1, 0.85*blue );
  drawoptions( withpen pencircle scaled 1.3pt );
  drawdblarrow rp(black)--rp(red);
  label.lrt( labro, rp(0.5*red) );
  drawdblarrow rp(black+0.85*blue)--rp(0.5*green+0.85*blue);
  label.urt( labri, rp(0.25*green+0.85*blue) );
  drawarrow rp(black)--rp(0.85*blue);
  label.ulft( laba, rp((0,0,0.5)) );
  dotlabel.rt( labc, rp(black) );
endfig;

beginfig(3);
  f := (3,2.05,2.65);
  numeric labscal, bray, sray;
  picture labrb, labc, laba, labrs;
  labscal = 2;
  ahlength := 8pt;
  ahangle := 30;
  bray = 0.7;
  sray = 0.4;
  labc = btex $\vec{c}$ etex scaled labscal;
  laba = btex $\vec{A}$ etex scaled labscal;
  labrb= btex $R_b$ etex scaled labscal;
  labrs= btex $R_s$ etex scaled labscal;
  drawoptions( withpen pencircle scaled 1.7pt withcolor 0.65blue );
  smoothtorus( black, blue, bray, sray );
  drawoptions( withpen pencircle scaled 1.3pt );
  draw goodcirclepath( bray*red, green, sray ) dashed evenly;
  drawdblarrow rp(black)--rp(bray*red);
  label.ulft( labrb, rp(0.5*bray*red) );
  drawdblarrow rp(bray*red+sray*blue)--rp(bray*red);
  label.llft( labrs, rp(bray*red+0.5*sray*blue) );
  drawarrow rp(black)--rp(sray*blue);
  label.rt( laba, rp((0,0,0.5*sray)) );
  dotlabel.lrt( labc, rp(black) );
endfig;

beginfig(4);
  f := (3,5,4);
  numeric labscal;
  picture labr, labc, laba;
  labscal = 2;
  ahlength := 8pt;
  ahangle := 30;
  labc = btex $\vec{c}$ etex scaled labscal;
  laba = btex $\vec{A}$ etex scaled labscal;
  labr = btex $R$ etex scaled labscal;
  drawoptions( withpen pencircle scaled 1.7pt withcolor 0.65blue );
  spatialhalfsfear( black, blue, 1 );
  drawoptions( withpen pencircle scaled 1.3pt );
  drawdblarrow rp(black)--rp(red);
  label.ulft( labr, rp(0.5*red) );
  drawarrow rp(black)--rp(0.4*blue);
  label.rt( laba, rp((0,0,0.2)) );
  dotlabel.lrt( labc, rp(black) );
endfig;

beginfig(5);
  f := (3,5,7);
  numeric labscal;
  picture labr, labc, laba;
  labscal = 2;
  ahlength := 10pt;
  ahangle := 30;
  labc = btex $\vec{c}$ etex scaled labscal;
  laba = btex $\vec{S}$ etex scaled labscal;
  labr = btex $R$ etex scaled labscal;
  drawoptions( withpen pencircle scaled 1.7pt withcolor 0.65blue );
  spheroid( black, 2*blue, 1 );
  drawoptions( withpen pencircle scaled 1.3pt );
  draw ellipticpath( black, red, 2*blue ) dashed evenly;
  draw ellipticpath( black, green, 2*blue ) dashed evenly;
  drawdblarrow rp(black)--rp(red);
  label.ulft( labr, rp((0.5,0,0)) );
  drawarrow rp(black)--rp(2*blue);
  label.rt( laba, rp((0,0,1)) );
  dotlabel.lrt( labc, rp(black) );
endfig;

beginfig(6);
  f := (3,5,4);
  numeric labscal;
  picture labl, labm, labr, labc;
  labscal = 2;
  ahlength := 10pt;
  ahangle := 30;
  labc = btex $\vec{c}$ etex scaled labscal;
  labm = btex $\vec{M}$ etex scaled labscal;
  labl = btex $L$ etex scaled labscal;
  labr = btex $R$ etex scaled labscal;
  draw rp((cosd(145),sind(-145),0))--rp(black)--rp((cosd(145),sind(145),0));
  drawarrow rp(black)--rp(red) withpen pencircle scaled 1.3pt;
  drawoptions( withpen pencircle scaled 1.7pt withcolor 0.65blue );
  banana( black, 1, black, 0.3, 145 );
  drawoptions( withpen pencircle scaled 1.3pt );
  drawdblarrow rp(black)--rp(0.3*blue);
  label.lft( labl, rp((0,0,0.15)) );
  drawdblarrow rp(0.3*blue)--rp(0.3*blue-green);
  label.llft( labr, rp(0.3*blue-0.5*green) );
  label.ulft( labm, rp(0.5*red) );
  dotlabel.bot( labc, rp(black) );
endfig;

beginfig(7);
  numeric labscal, ray;
  picture laba, labb, labr, labc;
  labscal = 2;
  ray = 0.25;
  ahlength := 10pt;
  ahangle := 30;
  labc = btex $\vec{c}$ etex scaled labscal;
  laba = btex $\vec{A}$ etex scaled labscal;
  labb = btex $\vec{B}$ etex scaled labscal;
  labr = btex $R$ etex scaled labscal;
  drawoptions( withpen pencircle scaled 1.7pt withcolor 0.65blue );
  quartertorus( black, -red, red-green, ray );
  drawoptions( withpen pencircle scaled 1.3pt );
  drawarrow rp(black)--rp(red-green) withpen pencircle scaled 1.3pt;
  label.bot( labb, rp((red-green)/sqrt(2)) );
  drawarrow rp(black)--rp(-red) withpen pencircle scaled 1.3pt;
  label.bot( laba, rp(-0.5*red) );
  drawdblarrow rp(-red)--rp(ray*blue-red);
  label.rt( labr, rp(0.5*ray*blue-red) );
  dotlabel.bot( labc, rp(black) );
  
endfig;

verbatimtex }\end{document} etex

end;
