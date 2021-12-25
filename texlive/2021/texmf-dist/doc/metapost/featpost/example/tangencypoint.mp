% tangencypoint.mp
% L. Nobre G.
% IYP (2005)

input featpost3Dplus2D;

beginfig(1);
  numeric u, i;
  u = 5mm;
  pen a, b, c;
  a = pencircle scaled 3pt;
  b = pencircle scaled 5pt;
  c = pencircle scaled 1pt;
  z1 = (1u,1u);
  z2 = (4u,4u);
  z3 = (4u,5u);
  z4 = (3u,5u);
  z5 = (3u,6u);
  z6 = (4u,7u);
  z7 = (6u,1u);
  path cp;
  cp = z1{up}..z2..z3..z4..z5..{up}z6;
  draw cp withpen c;
  for i=1 upto 6:
    draw z[i] withpen a withcolor 0.5*(red+green);
  endfor;
  z8 = firsttangencypoint( cp, z7, 5 );
  draw z7 withpen b withcolor green;
  draw z7--z8 withpen c withcolor blue;
endfig;

end.






