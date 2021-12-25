% intersection2D.mp
% L. Nobre G. 
% 2004

verbatimtex \documentclass{article}\begin{document}{\sffamily etex

beginfig(1);
  numeric u;
  pair Aa, Ab, Ba, Bb, Un, Kn, Or;
  pen dotpen, penpen;
  color knowncolour, notknowncol;
  picture kn, un;
  u = 4mm;
  Aa = (0,8u);
  Ab = (15u,0);
  Ba = (12u,0.5u);
  Bb = (13u,6u);
  Kn = (0,1u);
  Or = (1.5u,1.5u);
  dotpen = pencircle scaled 0.4u;
  penpen = pencircle scaled 0.2u;
  knowncolour = blue;
  notknowncol = red;
  kn = btex given etex;
  un = btex unknown etex;
  drawoptions( withpen penpen );
  draw Aa--Ab;
  draw Ba--Bb;
  drawoptions( withpen dotpen withcolor knowncolour );
  draw Aa;
  draw Ab;
  draw Ba;
  draw Bb;
  draw Kn+Or;
  drawoptions( withpen dotpen withcolor notknowncol );
  Un = whatever[Aa,Ab];
  Un = whatever[Ba,Bb];
  draw Un;
  draw Or;
  drawoptions();
  label.rt( kn, Kn+Or );
  label.rt( un, Or );
  label.lft( btex z1 etex, Aa );
  label.rt( btex z2 etex, Ab );
  label.lft( btex z3 etex, Ba );
  label.rt( btex z4 etex, Bb );
  label.urt( btex z5 etex, Un );
endfig;

verbatimtex }\end{document} etex

end;
