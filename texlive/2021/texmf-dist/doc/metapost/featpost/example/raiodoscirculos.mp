% raiodoscirculos.mp
% L. Nobre G.
% 2014

beginfig(1);
  numeric u, r, p;
  u = 2cm;
  r = 2/(3+sqrt(7));
  p = u/40;
  pen mypen;
  mypen = pencircle scaled p;
  drawoptions( withpen mypen );
  path mysquare, halcir, fulcir;
  mysquare = origin--right--(right+up)--up--cycle;
  halcir = halfcircle scaled 2r shifted (r*right);
  fulcir = fullcircle scaled 2r shifted ((1-r)*(right+up));
  draw mysquare scaled u;
  draw halcir scaled u;
  draw fulcir scaled u;
endfig;

beginfig(2);
  drawoptions( withpen mypen );
  draw mysquare scaled u;
  r := (sqrt(5)-1)/2;
  halcir := halfcircle shifted (right/2);
  draw halcir scaled u;
  path quarcirc;
  quarcirc = (halfcircle cutafter (origin--up)) rotated 180;
  draw quarcirc scaled 2r shifted (right+up) scaled u;
endfig;

beginfig(3);
  drawoptions( withpen mypen );
  draw mysquare scaled u;
  draw halcir scaled u;
  r := 2-sqrt(3);
  draw fullcircle scaled 2r shifted ((1-r)*(right+up)) scaled u;
endfig;

beginfig(4);
  drawoptions( withpen mypen );
  draw mysquare scaled u;
  draw halcir scaled u;
  r := 1/3;
  path secondhalfcircle;
  secondhalfcircle = halfcircle scaled 2r rotated 90 shifted (1,1-r);
  draw secondhalfcircle scaled u;
endfig;

beginfig(5);
  drawoptions( withpen mypen );
  draw mysquare scaled u;
  r := sqrt(2)-1;
  draw quarcirc scaled 2 shifted (right+up) scaled u;
  draw quarcirc rotated 180 scaled (2r*u);
endfig;

beginfig(6);
  drawoptions( withpen mypen );
  draw mysquare scaled u;
  r := 1/4;
  halcir := halfcircle scaled 0.5 shifted (r*right);
  draw halcir scaled u;
  draw quarcirc scaled 2 shifted (right+up) scaled u;
endfig;

beginfig(7);
  drawoptions( withpen mypen );
  path thisquare;
  thisquare = mysquare scaled u;
  draw thisquare;
  r := (sqrt(7)-1)/3;
  halcir := (halfcircle scaled 2r shifted (r*right)) cutbefore mysquare;
  draw halcir scaled u;
  draw quarcirc scaled 2r shifted (right+up) scaled u;
endfig;

beginfig(8);
  drawoptions( withpen mypen );
  draw mysquare scaled u;
  r := 1/(1+sqrt(2));
  draw quarcirc rotated 180 scaled (2r*u);
  draw fullcircle scaled 2r shifted ((1-r)*(right+up)) scaled u;
endfig;

beginfig(9);
  drawoptions( withpen mypen );
  draw mysquare scaled u;
  r := 1/(1+sqrt(2));
  halcir := halfcircle scaled 2r shifted (r*right);
  draw halcir scaled u;
  secondhalfcircle := halfcircle scaled 2r rotated 90 shifted (1,1-r);
  draw secondhalfcircle scaled u;
endfig;

end.