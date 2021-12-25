% cartaxes.mp
% L. Nobre G.
% 2006

input featpost3Dplus2D;

verbatimtex \documentclass{article}\begin{document} etex

beginfig(1);
  cartaxes(1,1,1);
endfig;

f := (5,3.5,1);
beginfig(2);
  cartaxes(1,1,1);
  color va, vb, vc, vd;
  va = (0.29,0.7,1.0);
  vb = (X(va),Y(va),0);
  vc = N((-Y(va),X(va),0));
  vd = (0,Y(vc),0);
  drawarrow rp(black)--rp(va);
  draw rp(black)--rp(vb)--rp(va) dashed evenly;
  draw rp(vc)--rp(vd) dashed evenly;
  drawarrow rp(black)--rp(vc);
  squareangline( va, vc, black, 0.15 );
  angline( va, red, black, 0.75, decimal getangle( va, red ), lft );
endfig;

verbatimtex \end{document} etex

end;
