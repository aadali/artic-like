% solardish.mp
% L. Nobre G.
% 2013

prologues := 1;

vardef reflectrayr(expr s,a,p,l)=
  save tI,tn,ia,I,J;
  pair I,J;
  tI=xpart(p intersectiontimes ((s+1mm*dir(a))--(s+30cm*dir(a))));
  if tI>0:
    I=point tI of p;
    draw s--I;
    tn=angle(direction tI of p)+90;
    ia:=tn-angle(s-I);
    reflectrayr(I,tn+ia,p,l);
  else:
    J=s+dir(a)*l;
    drawarrow s--J;
  fi;
enddef;

beginfig(1); % parabolic dish
  numeric u, i, num, calx, incidang;
  pair sunpoi;
  u = 8cm;
  num = 12;
  incidang = -60;
  path p;
  p=(-2u,u) for i=1 upto 2num:
      hide( calx:=-2+2*i/num; )
      ..(calx*u,u*(calx**2)/4)
  endfor;
  for i=1 upto 2num-1:
    sunpoi := ((-2+2*i/num,1)*u)-u*dir(incidang);
    reflectrayr(sunpoi,incidang,p,14cm);
  endfor;
  draw p withcolor red;
endfig;

beginfig(2); % circle evolute dish
  numeric u, i, num, calc, incidang, panray;
  pair sunpoi, crosa, crosb;
  u = 8cm; % this value is related with "14cm" inside reflectray
  num = 12;
  incidang = -60;
  panray = 2cm;
  path p, q, w;
  w = fullcircle scaled 2panray shifted (0,u-panray);
  p=(0,u-2panray) for i=1 upto num:
      hide(
	calc:=(i/num)*panray*3.14159;
        sunpoi := (0,u-panray)+dir(i*180/num-90)*panray-calc*dir(i*180/num);
	)
      ..sunpoi
  endfor;
  q = reverse (p xscaled -1) & p;
  sunpoi := (xpart point 2 of q,u)-u*dir(incidang);
  crosa = sunpoi+whatever*dir(incidang);
  crosa = (0,u-panray)-whatever*dir(incidang+90);
  crosb = ((0,u-panray)--crosa) intersectionpoint w;
  for i=0 upto num:
    sunpoi := (i/num)[crosa,crosb]-u*dir(incidang);
    reflectrayr(sunpoi,incidang,q,5cm);
  endfor;
  pair crosa;
  sunpoi := (xpart point ((length q)-2) of q,u)-u*dir(incidang);
  crosa = sunpoi+whatever*dir(incidang);
  crosa = (0,u-panray)-whatever*dir(incidang+90);
  crosb := ((0,u-panray)--crosa) intersectionpoint w;
  for i=0 upto num:
    sunpoi := (i/num)[crosa,crosb]-u*dir(incidang);
    reflectrayr(sunpoi,incidang,q,5cm);
  endfor;
  draw q withcolor red;
  fill w withcolor 0.5white;
endfig;

beginfig(3); % spheric dish
  numeric u, i, num, incidang;
  pair sunpoi;
  u = 8cm;
  num = 12;
  incidang = -60;
  path p;
  p = halfcircle scaled 4u rotated 180 shifted (0,u);
  for i=1 upto 2num-1:
    sunpoi := ((-2+2*i/num,1)*u)-u*dir(incidang);
    reflectrayr(sunpoi,incidang,p,14cm);
  endfor;
  draw p withcolor red;
endfig;

beginfig(4); % elliptic dish
  numeric u, i, num, cal, incidang;
  pair sunpoi;
  u = 8cm;
  num = 12;
  incidang = -55;
  path p;
  p=(-2u,u) for i=1 upto 2num:
      hide( cal:=90*i/num; )
      ..(-2*cosd(cal)*u,u-u*sind(cal))
  endfor;
  for i=1 upto 2num-1:
    sunpoi := ((-2+2*i/num,1)*u)-u*dir(incidang);
    reflectrayr(sunpoi,incidang,p,14cm);
  endfor;
  draw p withcolor red;
endfig;

end.