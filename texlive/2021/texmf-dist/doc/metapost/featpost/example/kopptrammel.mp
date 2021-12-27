% kopptrammel.mp
% L. Nobre G.
% IYP (2005)

beginfig(1);
  numeric i, u, astep, pfrac, ra, rb, totray, adiff, astart;
  pair ca, cb, pa, pb, pc;
  path cir;
  ra = 0.3;
  rb = 0.9;
  astep = 10;
  adiff = -20;
  pfrac = 3;
  astart = 20;
  u = 2cm;
  ca = (2u,4u);
  cb = (4u,4u);
  totray = 0.5*abs(cb-ca);
  cir = fullcircle scaled (2*totray);
  draw cir shifted ca;
  draw cir shifted cb;
  pickup pencircle scaled (u/20);
  for i=0 step astep until 360:
    pa := ca + ra*totray*dir( astart+i );
    pb := cb + rb*totray*dir( astart-i-adiff );
    pc := pfrac[pa,pb];
    draw pa;
    draw pb;
    draw pc withcolor red;
  endfor;
  pickup pencircle scaled (u/40);
  draw pa--pc withcolor blue;
  draw pb--pc;
  draw ca;
  draw cb;
  draw pfrac[ca,cb];
endfig;

end.