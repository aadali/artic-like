% parabolaonthestage.mp
% L. Nobre G.
% 2014

input featpost3Dplus2D;
Spread := 30;
f := (5,3,2);
beginfig(1);
  numeric i, li;
  li = 50;
  setthestage( 18, 10 );
  draw for i=-li step 0.25 until li-0.25:
      rp((5-(i**2),i,0))...
  endfor rp((5-(li**2),li,0)) withcolor blue withpen pencircle scaled 3pt;
endfig;

end.