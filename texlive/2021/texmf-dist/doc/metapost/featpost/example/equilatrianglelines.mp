% equilatrianglelines.mp
% L. Nobre G.
% 2009

prologues := 1;

beginfig(1);
  numeric u, factordown, ang, angl, angh, ydown, xl, xh, astep;
  pen pendot;
  path uline, rline, dline;
  pair auxp;
  u = 2cm;
  factordown = 2.34;
  astep = 2.75;
  xl = -2u;
  xh = 10u;
  angl = 7;
  angh = 38;
  pendot = pencircle scaled (2.5mm);
  ydown = -u*factordown;
  uline = (xl,u)--(xh,u);
  rline = uline rotated (-60);
  dline = (xl,ydown)--(xh,ydown);
  draw uline;
  draw (xl,0)--(xh,0);
  draw dline;
  draw rline;
  z1 = rline intersectionpoint dline;
  draw z1 withpen pendot;
  draw origin withpen pendot;
  for ang=angl step astep until angh:
    pair uapex;
    ypart uapex = u;
    uapex = whatever*dir(ang);
    auxp := uapex+abs(uapex)*dir(ang-120);
    draw origin--uapex--auxp--cycle;
  endfor;
endfig;

end.