% muslimpattern.mp
% L. Nobre G.
% 2012

prologues := 1;
background := 0.9white;

numeric u, pa, pb, hexf;
u = 35mm;
pa = 0.12u;
pb = 0.05u;
hexf = sqrt(3)/2;

vardef markpath( expr Pat ) =
  draw Pat withpen pencircle scaled pa;
  undraw Pat withpen pencircle scaled pb;
enddef;

vardef hc( expr Xindex, Yindex ) =
  2*u*hexf*(Xindex*right+Yindex*dir(60))
enddef;

beginfig(1);
  numeric i, ang;
  pair c, b, a, d;
  a = (hexf,0);
  b = (0.33hexf,-0.4); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHANGE only this.
  c = (0,1.5);
  d = 1.1[b,a];
  path seg[], aux, hex, auxa, auxb;
  for i=0 upto 5:
    seg[i] = (d--b--c) rotated (i*60) scaled u;
  endfor;
  seg6 = seg0;
  ang = angle(c-b);
  aux = ((a-dir(ang))--(a+dir(ang))) scaled u;
  z1 = seg0 intersectionpoint seg5;
  z2 = seg0 intersectiontimes seg5;
  z3 = unitvector( direction xpart z2 of seg0 );
  for i=1 upto 6:
    markpath( aux rotated (i*60) );
    markpath( seg[i-1] );
  endfor;
  auxa = ((a-pb*dir(ang)/u)--(a+pb*dir(ang)/u)) scaled u;
  auxb = ((a-pa*dir(ang)/u)--(a+pa*dir(ang)/u)) scaled u;
  for i=1 upto 6:
    draw auxa rotated ((i-1)*60) withpen pencircle scaled pa;
    undraw auxb rotated ((i-1)*60) withpen pencircle scaled pb;
  endfor;
  auxa := ((z1-pb*z3)--(z1+pb*z3));
  auxb := ((z1-pa*z3)--(z1+pa*z3));
  draw auxa withpen pencircle scaled pa;
  undraw auxb withpen pencircle scaled pb;
  hex = for i=30 step 60 until 360: u*dir(i)-- endfor cycle;
  clip currentpicture to hex;
  picture tile;
  tile = currentpicture;
endfig;

beginfig(2);
  numeric i, j;
  path rect;
  rect = hc(0,-1)--hc(1,-1)--hc(0,1)--hc(-1,1)--cycle;
  fill rect withcolor (1,0.986,0.993);
  for i=-1 upto 1:
    for j=-1 upto 1:
      draw tile shifted hc(i,j);
    endfor;
  endfor;
  clip currentpicture to rect;
endfig;

end.