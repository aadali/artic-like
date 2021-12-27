% iplaca.mp
% L. Nobre G.
% 2012

prologues := 1;

numeric u;
u = 52mm;

def placaminho( expr yv ) =
  begingroup
    pair p[];
    numeric xh;
    xh = u-yv;
    p1 = (xh,yv);
    p2 = (-xh,yv);
    p3 = (-xh,-yv);
    p4 = (0,-yv);
    p5 = (-yv,0);
    p6 = (0,yv);
    p7 = (yv,0);
    p8 = (xh,-yv);
    ( p1--p2--p3--p4--p5--p6--p7--p4--p8--cycle )
  endgroup
enddef;

beginfig(1);
  path a, b, c;
  color g;
  pen j;
  j = pencircle scaled 0.02u;
  g = 0.5white;
  a = placaminho( 0.15u ) shifted (left*u*1.9);
  fill a withcolor g;
  draw a withpen j;
  b = placaminho( 0.3u );
  fill b withcolor g;
  draw b withpen j;
  c = placaminho( 0.45u ) shifted (right*u*1.7);
  fill c withcolor g;
  draw c withpen j;
endfig;

end.
    