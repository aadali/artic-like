% intersection2D.mp
% L. Nobre G. 
% 2004

input featpost3Dplus2D;

verbatimtex \documentclass{article}\begin{document}{\sffamily etex

beginfig(1);
  numeric u, r, d;
  pair Aa, Ab, Ba, Bb, Ca, Cb, Da, Db, Ea, Eb, Or, Kn;
  pen dotpen, penpen;
  color knowncolour, notknowncol;
  picture kn, un, rl;
  u = 4mm;
  Aa = (5u,4u);
  Ab = (3u,2u);
  Ba = (13u,5u);
  Bb = 0.5[Aa,Ab];
  r = abs( Aa - Ab );
  d = abs( Aa - Ba ) +-+ r;
  Ea = radialcross( Aa, r, Ba, d, true );
  Eb = radialcross( Aa, r, Ba, d, false );
  Ca = 1.4[Ba,Ea];
  Cb = 1.15[Ea,Ba];
  Da = 1.5[Ba,Eb];
  Db = 1.1[Eb,Ba];
  Kn = (0,1u);
  Or = (10u,1u);
  dotpen = pencircle scaled 0.5u;
  penpen = pencircle scaled 0.2u;
  knowncolour = blue;
  notknowncol = red;
  kn = btex given etex;
  un = btex unknown etex;
  rl = btex $R$ etex; 
  draw fullcircle scaled 2r shifted Aa;
  draw Ca--Cb;
  draw Da--Db dashed evenly;
  drawarrow Aa--Ab;
  drawoptions( withpen dotpen withcolor knowncolour );
  draw Aa;
  draw Ba;
  draw Kn+Or;
  drawoptions( withcolor knowncolour );
  label.ulft( rl, Bb );
  drawoptions( withpen dotpen withcolor notknowncol );
  draw Ea;
  draw Eb;
  draw Or;
  drawoptions();
  label.rt( kn, Kn+Or );
  label.rt( un, Or );
endfig;

verbatimtex }\end{document} etex

end;
