% ellipticprism.mp
% L. Nobre G.
% 2004

input featpost3Dplus2D;

verbatimtex \documentclass{article}\begin{document}{\sffamily etex

beginfig(1);
  f := (5,3,4);
  numeric labscal;
  picture labsa, labba, labsb, labd, labc;
  labscal = 2;
  ahlength := 8pt;
  ahangle := 30;
  labc = btex $\vec{c}$ etex scaled labscal;
  labsa = btex $\vec{S}_1$ etex scaled labscal;
  labsb = btex $\vec{S}_2$ etex scaled labscal;
  labba = btex $\vec{B}_1$ etex scaled labscal;
  labd= btex $D$ etex scaled labscal;
  drawoptions( withpen pencircle scaled 1.7pt withcolor 0.65blue );
  whatisthis( black, 0.5red, green, 0.85, 0.8 );
  drawoptions( withpen pencircle scaled 1.3pt );
  drawdblarrow rp(black)--rp(0.85blue);
  label.rt( labd, rp(0.4blue) );
  drawarrow rp(black)--rp(0.5red);
  label.ulft( labsa, rp(0.25red) );
  drawarrow rp(black)--rp(green);
  label.top( labba, rp((0,0.5,0)) );
  drawarrow rp(0.85blue)--rp(0.8*0.5*red+0.85blue);
  label.ulft( labsb, rp(0.2red+0.85blue) );
  dotlabel.bot( labc, rp(black) );
endfig;
  
verbatimtex 
\end{document} 
etex

end.

  
