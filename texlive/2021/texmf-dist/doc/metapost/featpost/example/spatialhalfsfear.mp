% spatialhalfsfear.mp
% L. Nobre G.
% 2004

input featpost3Dplus2D;

f := (3,-0.2,3);
beginfig(1);
  spatialhalfsfear( black, blue, 1 );
  drawoptions( withpen pencircle scaled 3pt withcolor (blue+green) );
  spatialhalfsfear( black, blue+green, 0.8 );
  drawoptions( withpen pencircle scaled 3pt withcolor red );
  spatialhalfsfear( black, red, 0.6 );
  drawoptions( withpen pencircle scaled 3pt withcolor green );
  spatialhalfsfear( black, green, 0.7 );
endfig;
end.
