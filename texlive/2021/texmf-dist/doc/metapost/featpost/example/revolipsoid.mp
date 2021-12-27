% revolipsoid.mp
% L. Nobre G.
% 2011

input featpost3Dplus2D;

prologues := 1;

Spread := 20;
f := 7.1*(3,5,4);
ShadowOn := true;
HoriZon := -5;
LightSource := 10*(4,-2,4);

beginfig(1);
  setthestage(30,45);
  drawoptions( withpen pencircle scaled 1mm );
  numeric num, r, ste, i, j, aux;
  r = 1.6;
  num = 5;
  aux = 0.5*(num-1);
  ste = 90/aux;
  for i=-aux upto aux:
    for j=-aux upto aux:
      spheroid((red*j+green*i)*2*(r**aux),(cosd(i*ste),0,sind(i*ste)),r**(-j));
    endfor;
  endfor;
endfig;

ShadowOn := false;

beginfig(2);
  drawoptions( withpen pencircle scaled 1mm );
  numeric num, r, ste, i, j, aux, fac;
  color aim;
  r = 0.55;
  fac = 0.6;
  num = 6;
  aux = 0.5*(num-1);
  ste = 90/aux;
  for i=-aux upto aux:
    for j=-aux upto aux:
      aim := N(f-viewcentr)+fac*(i*red+j*green)/aux;
      spheroid((red*j+green*i)*3*r,aim,r);
    endfor;
  endfor;
endfig;

end.
  
