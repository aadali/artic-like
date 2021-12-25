% reuleaux.mp
% L. Nobre G.
% 2012

prologues := 1;

beginfig(1);
  numeric u, desv, ray;
  path rouletri, pathpart, pa, pb, pc, border;
  u = 150mm;
  desv = u*sqrt(3)/3;
  ray = (sqrt(3)/3-0.5)*u/2;
  pathpart = halfcircle scaled u cutafter (origin--u*dir(60));
  pathpart := pathpart shifted (desv*dir(-150)/2);
  pa = pathpart;
  pb = pathpart rotated 120;
  pc = pathpart rotated -120;
  rouletri = pa--pb--pc--cycle;
  draw rouletri;
  border=(u/4,u/4)--(-u/4,u/4)--(-u/4,-u/4)--(u/4,-u/4)--cycle;
  draw border shifted (ray*up);
endfig;

beginfig(2);
  draw border withcolor red;
  numeric i, N, astep;
  N = 4;
  astep = 3*N;
  for i=astep step astep until 360:
    draw rouletri rotated (-i/3) shifted (ray*dir(i-90));
  endfor;
endfig;

end.
