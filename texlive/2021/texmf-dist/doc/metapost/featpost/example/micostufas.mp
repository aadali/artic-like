% micostufas.mp
% L. Nobre G.
% 2013

prologues := 1;

beginfig(1);
  numeric u, basecant, altucnat, diagcant, baseestu, altuestu, diagestu;
  numeric aux, rayestu;
  u = 1cm;
  basecant = 0.5u;
  diagcant = 2.6u;
  altucant = diagcant +-+ 0.5basecant;
  baseestu = 6u;
  altuestu = 4u;
  diagestu = altuestu ++ 0.5baseestu;
  numeric angalpha, scalefac, cantdesv;
  angalpha = 180-2*angle(0.5baseestu,altuestu); %show angalpha;
  rayestu = altuestu/sind(angalpha);
  message "Raio de curvatura: " & decimal(rayestu/u);
  cantdesv = 0.25*(baseestu+basecant);
  path auxpa, cant, goticarch;
  cant = (-0.5basecant,0)--(0,altucant)--(0.5basecant,0);
  draw cant;
  draw cant shifted (-cantdesv,0);
  draw cant shifted (cantdesv,0);
  z1 = (-0.5*cosd(angalpha),0);
  auxpa = (halfcircle cutafter (origin--dir(angalpha))) shifted z1;
  aux = xpart point 0 of auxpa;
  scalefac = 0.5*baseestu/aux;
  goticarch = auxpa scaled scalefac;
  draw goticarch;
  draw goticarch xscaled -1;
  z0 = z1 scaled scalefac;
  drawarrow z0--(z0+rayestu*dir(angalpha/2));
  numeric comprimarco;
  comprimarco = rayestu*3.14159*angalpha/(180*u);
  message "Comprimento de cada arco: " & decimal(comprimarco);
endfig;

end.