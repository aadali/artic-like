% interdigitrodes.mp
% L. Nobre G.
% 2012

prologues := 1;

beginfig(1);
  numeric u, t, i, j, len;
  u = 1mm;
  t = 4u;
  i = 1u;
  len = 2*(t+i);
  path basepath;
  basepath = ((0,0.5t)---(3len,0.5t)--(3len,-0.5t)---(0,-0.5t)..cycle)
    shifted (-len,0.5*(t+i));
  fill basepath;
  fill basepath shifted (0,len);
  fill basepath shifted (0,-len);
  fill basepath shifted (0,-2len);
  fill (2len,-2len)--(3len,-2len)--(3len,2len)--(2len,2len)--cycle;
  draw currentpicture rotated 180 withcolor 0.333white;
endfig;

end.