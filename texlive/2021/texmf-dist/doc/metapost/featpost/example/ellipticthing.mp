% ellipticthing.mp
% L. Nobre G.
% 2011

prologues := 1;

beginfig(1);
  numeric u, a, b, i, n, astep;
  u = 1cm;
  a = 5u;
  b = 3u;
  n = 24;
  astep = 360/n;
  path elli, arcb;
  elli = (a,0) for i=1 upto (n-1):
      ..(a*cosd(i*astep),b*sind(i*astep))
  endfor ..cycle;
  draw elli;
  draw fullcircle scaled (2*b);
  draw fullcircle scaled (a-b) shifted (0.5*(a+b),0);
  draw fullcircle scaled (a-b) shifted (-0.5*(a+b),0);
  pair refpoint, refslope, centerpo, focus;
  refpoint = (a,b)/sqrt(2);
  refslope = unitvector(b,a);
  centerpo = refpoint + whatever*refslope;
  xpart centerpo = 0;
  draw centerpo--refpoint;
  draw centerpo--(-xpart refpoint,ypart refpoint);
  numeric len, ang, lenb, angb;
  len = abs(refpoint-centerpo);
  ang = angle(refpoint-centerpo);
  draw refpoint for i=1 upto n:
      ..(centerpo+len*dir(2*i*(90-ang)/n+ang))
  endfor;
  focus = (4u,0);
  lenb = abs(focus-centerpo);
  angb = angle(focus-centerpo);
  arcb = focus for i=1 upto n:
      ..(centerpo+lenb*dir(2*i*(90-angb)/n+angb))
  endfor;
  draw arcb;
  draw arcb rotated 180;
  draw centerpo--focus;
  draw centerpo--(-xpart focus,ypart focus);
endfig;

end.