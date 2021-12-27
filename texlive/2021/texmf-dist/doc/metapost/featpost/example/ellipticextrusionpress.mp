% ellipticextrusionpress.mp
% L. Nobre G.
% 2011

prologues := 1;

beginfig(1);
  numeric u, a, b, i, j, n, angstep, halflen;
  pair vec[];
  u = 1cm;
  a = 8u;
  b = 2u;
  n = 40;
  angstep = 180/n;
  halflen = a*0.5*3.14159;
  draw (-halflen,u)--(halflen,u);
  draw (-u,u)--(-u,0)--(u,0)--(u,u);
  for i=0 upto n:
    vec[i] = (a*sind(i*angstep),-b*cosd(i*angstep));
    draw ((-halflen,a)+vec[i])
      for j=1 upto n:
      --((a*j*3.14159/n-halflen,a)+(vec[i] rotated (-j*angstep)))
    endfor;
  endfor;
  draw ((0,a)+(vec[0] rotated (-90)))
    for j=1 upto n:
    ..((0,a)+(vec[j] rotated (-90)))
  endfor;
endfig;

end.