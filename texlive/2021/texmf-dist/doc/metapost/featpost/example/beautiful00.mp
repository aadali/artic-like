% beautiful00.mp
% L. Nobre G.
% 2014

prologues := 1;

beginfig(1);
  numeric u, firstep, cylray, plangle, i, a, as, speed, dist, accel, li, aux;
  pair poi, tip, ref;
  pen one;
  accel = 0.27;
  u = 1.2cm;
  firstep = 6u;
  cylray = 2u;
  plangle = 12;
  as = 5;
  ahangle := 30;
  ahlength := 0.035cylray;
  one = pencircle scaled 0.015cylray;
  drawoptions( withpen one );
  li = 3;
  aux = li*accel*firstep*(li+1);
  fill (-cylray,0)--(aux,0)--(aux,-ahlength)--(-cylray,-ahlength)--cycle withcolor 0.5white;
  for i=0 upto li:
    speed := i*accel;
    ref := (speed*firstep*i,0);
    for a=as-90 step as until 270-as:
      poi := ref+cylray*(dir(a)+up);
      dist := abs( poi-ref );
      tip := poi+speed*dist*dir(0.5*a-45);
      drawarrow poi--tip;
    endfor;
  endfor;
  picture tmpicture;
  tmpicture = currentpicture;
  currentpicture := nullpicture;
  draw tmpicture rotated -plangle;
endfig;

end.