% revolutionparaboloid.mp
% L. Nobre G.
% 2014

%outputformat := "png";
%outputformatoptions := "format=rgb antialias=best";
%hppp := 0.5;
%vppp := 0.5;

input featpost3Dplus2D;

%prologues := 1;

f := 0.7*(3,3,-5);
  
beginfig(1);
  numeric u, r;
  color b, t;
  u = 2mm;
  b = (-1.0,0.0,0.0);
  t = (-0.4,0.3,0.9);
  r = 0.6;
  drawoptions(withpen pencircle scaled u);
  revolparab( b, t, r );
endfig;

def parab( expr xis ) = (0.25*(xis**2)) enddef;
def oloid( expr xis ) = (1-(xis/2)**2) enddef;
def deriv( expr xis ) = (0.5*xis) enddef;

beginfig(2);
  pair pzero, pone;
  pzero = (-0.5,-5.5);
  numeric u, delta, n, ste, minx, maxx, itex, tanx, marg, slope, sa, aux;
  n = 50;
  marg = 1;
  u = 5mm;  
  delta = sqrt((xpart pzero)**2-4*(ypart pzero));
  minx = (xpart pzero)-delta-marg;
  tanx = (xpart pzero)+delta;
  maxx = -minx;
  ste = (maxx-minx)/(n+1);
  z0 = (minx,parab(minx));
  z1 = ((xpart pzero)-delta,parab((xpart pzero)-delta));
  z2 = (tanx,parab(tanx));
  z3 = 0.5[z1,z2];
  draw z0 for itex=minx+ste step ste until maxx:
      --(itex,parab(itex)) endfor withcolor blue;
  draw pzero--z1--z2--cycle;
  draw pzero--z3;
  sa = -angle(z1-z2);
  slope = sind(sa)/cosd(sa);
  aux = 2*(-slope+(slope++1));
  z4 = (aux,parab(aux));
  z6 = z4+whatever*dir(angle(1,deriv(aux)));
  aux := 2*(-slope-(slope++1));
  z5 = (aux,parab(aux));
  z6 = z5+whatever*dir(angle(1,deriv(aux)));
  draw z6--z4--z5--cycle withcolor red;
  %draw up withcolor blue;
  picture scalpic;
  scalpic = currentpicture;
  currentpicture := nullpicture;
  draw scalpic scaled u withpen pencircle scaled 0;
endfig;

beginfig(3);
%  ParallelProj := true;
  f := (7, 6, 5);
%  f := (2, 1, 8);
%  f := (7,6,8);
%  f := (0.05,0,9.5);
%  f := (0,9,0);
  numeric ray, coneray, angrotparab, angrotfoc, steprotfoc, itex, xzer, yzer;
  numeric tmpxzero, tmpdelta, tmpxtplu, tmpxtmin, stepx;
  color focus, tmpxaxis, auxca, auxcb, auxcc, auxcd;
  pen penone;
  pair auxa, auxb;
  penone = pencircle scaled 3mm;
  xzer = 1;
  yzer = -0.5;
  stepx = 0.05;
  ray = 3;
  revolparab(blue*parab(ray),black,ray);
  coneray = 4;
  steprotfoc = 15;
  focus = (xzer,0,yzer);
  for angrotfoc=0 step steprotfoc until 360-steprotfoc:
    angrotparab := 0.5*angrotfoc;
    tmpxzero := xzer*cosd(angrotparab);
    tmpxaxis := (cosd(angrotparab),sind(angrotparab),0);
      z0 = (minx,parab(minx));
    draw rp(-ray*tmpxaxis+blue*parab(-ray))
      for itex=-ray+stepx step stepx until ray:
      --rp(itex*tmpxaxis+blue*parab(itex)) endfor withcolor blue;
    tmpdelta := sqrt(tmpxzero**2-4*yzer);
    tmpxtplu := tmpxzero+tmpdelta;
    tmpxtmin := tmpxzero-tmpdelta;
    auxca := tmpxtplu*tmpxaxis+blue*parab(tmpxtplu);
    auxcb := tmpxtmin*tmpxaxis+blue*parab(tmpxtmin);
    auxa := rp(auxca);
    auxb := rp(auxcb);
    auxcc := auxca-focus;
    auxcd := auxcb-focus;
    draw auxa withpen penone;
    draw auxb withpen penone;    
  endfor;
  draw rp(focus) withpen penone withcolor green;
  produce_auto_scale;
endfig;

beginfig(4);
  numeric major, minor, focdi, rotaroun, steparoun, ra;
  color majorc, minorc, focdic, garbac, vertec;
  major = 4;
  minor = 2.5;
  f := (-5,0.5,0.4);
  steparoun = 20;
  ra = 0.35;
  focdi = major +-+ minor;
  majorc = major*red;
  minorc = minor*green;
  focdic = -focdi*red;
  for rotaroun=0.5*steparoun step steparoun until 180-0.5*steparoun:
    garbac := planarrotation(majorc,minorc,rotaroun);
    vertec :=0.5[garbac,focdic];
    revolparab( garbac, vertec, ra );
    revolparab( (X(garbac),-Y(garbac),0), (X(vertec),-Y(vertec),0), ra );
  endfor;
endfig;

beginfig(5);
  f := (9,6,3);
  numeric u, slo, delt, higx, lowx, cenx, rai, focz;
  color cen, foc, low, base;
  path ell, cir, both;
  u = 2mm;
  slo = 0.35;
  delt =  slo ++ 1;
  lowx = 2*(slo-delt);
  higx = 2*(slo+delt);
  rai = 2*delt;
  drawoptions(withpen pencircle scaled u withcolor blue);
  revolparab( (0,0,-3), blue, 4 );
  drawoptions(withpen pencircle scaled u);  
  z4 = (lowx,oloid(lowx));
  z6 = z4+whatever*dir(angle(1,-deriv(lowx)));
  z5 = (higx,oloid(higx));
  z6 = z5+whatever*dir(angle(1,-deriv(higx)));
  focz = ypart z6;
  cen = (2*slo,0,ypart 0.5[z4,z5]);
  foc = (2*slo,0,focz);
  low = (lowx,0,ypart z4);
  ell = ellipticpath( cen, low-cen, rai*green );
  cir = fullcircle scaled u;
  both = twocyclestogether( ell, cir shifted rp(foc) );
  draw both;
  draw ell;
  draw rp(2*cen-low)--rp(low)--rp(foc)--cycle;
  %draw rp(cen+rai*green)--rp(cen-rai*green);
  fill cir shifted rp(black);
  base = (2*slo,0,-3);
  draw rigorouscircle( base, blue, rai );
  draw rp(base)--rp(foc) dashed withdots withcolor red;
endfig;

end.