% jdbug.mp
% L. Nobre G.
% 2012

prologues := 1;

beginfig(1);
  numeric N, i, u, c, s;
  pen thepen;
  picture lastpic;
  u = 5mm;
  thepen = pencircle scaled 0.2u;
  N = 5;
  s = 180/N;
  c = u*cosd(s);
  draw for i=0 step 2s until 360-2s:
      (u*dir(i))--(c*dir(i+s))--
  endfor cycle withpen thepen;
  lastpic = currentpicture;
endfig;

beginfig(2);
  draw lastpic rotated s;
endfig;

end.