% cone.mp
% L. Nobre G.
% 2004

input featpost3Dplus2D;

verbatimtex \documentclass{article}\begin{document}{\sffamily etex

beginfig(1);
  f := (5,3,4);
  numeric labscal;
  picture labr, labc, laba, labv;
  labscal = 2;
  ahlength := 8pt;
  ahangle := 30;
  labc = btex $\vec{c}$ etex scaled labscal;
  laba = btex $\vec{A}$ etex scaled labscal;
  labr = btex $R$ etex scaled labscal;
  labv = btex $\vec{v}$ etex scaled labscal;
  drawoptions( withpen pencircle scaled 1.7pt withcolor 0.65blue );
  verygoodcone(true,black,blue,0.8,blue+green);
  drawoptions( withpen pencircle scaled 1.3pt );
  draw rp(black)--rp(green)--rp(blue+green) dashed evenly;
  drawdblarrow rp(black)--rp(0.8green);
  label.urt( labr, rp(0.4green) );
  drawarrow rp(black)--rp(0.45*blue);
  label.rt( laba, rp((0,0,0.25)) );
  dotlabel.llft( labc, rp(black) );
  dotlabel.rt( labv, rp(blue+green) );

endfig;
  
verbatimtex 
\end{document} 
etex

end.

  
