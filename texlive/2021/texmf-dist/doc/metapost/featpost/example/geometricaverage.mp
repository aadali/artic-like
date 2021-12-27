% geometricaverage.mp
% L. Nobre G.
% 2014

input featpost3Dplus2D;

verbatimtex \documentclass{article}\begin{document} etex

%f := (5,3.5,1);
%f := 0.4*(5,3.5,3);
%f := 0.4*(5,1.5,1.5);
f := 0.4*(2.5,5,1.5);
beginfig(1);
  cartaxes(1.5,2,1);
  color va, vb, vc, vd, ve, vf;
  va = black;
  vb = (0,1.5,0);
  vc = (1,1.5,0);
  vd = red;
  ve = (1,0.75,0);
  vf = (0.5,1.5,0);
  z3 = rp(va);
  z4 = rp(vb);
  z5 = rp(vc);
  z6 = rp(vd);
  z7 = rp(ve);
  z9 = rp(vf);
  draw z3--z4--z5--z6--cycle dashed evenly;
  z1 = whatever[z3,z4];
  z1 = whatever[z6,z5];
  z2 = whatever[z6,z3];
  z2 = whatever[z5,z4];
  draw z1--z3 dashed evenly;
  draw z1--z6 dashed evenly;
  draw z2--z3 dashed evenly;
  draw z2--z4 dashed evenly;
  draw z2--z7 withcolor red;
  draw z1--z9 withcolor red;
  numeric laa, lba, lca;
  laa = abs(z3-z1);
  lba = abs(z4-z1);
  lca = sqrt( laa )*sqrt( lba );
  z8 = (lca/lba)[z1,z4];
  draw z2--z8;
  
  laa := 1/abs(z3-z2);
  lba := 1/abs(z6-z2);
  lca := 0.5*( laa + lba ); % What works is the arithmetic average of
                            % reciprocal distances!!
  z10 = (lba/lca)[z2,z6];
  draw z1--z10;
  labels.top(1,2,3,4,5,6,7,8,9,10);

  z11 = vp( vf-vd );
  drawarrow z9--z11;
endfig;

verbatimtex \end{document} etex

end;
