% ellipticproperties.mp
% L. Nobre G.
% 2011

prologues := 1;

beginfig(1);
  numeric b, k, kt, theta, i, ste, N, aux, cp, cm;
  numeric gamma, gammb, alpha, rhoa, alphb, rhob, fx, fy;
  pair X, Y, A, B, U, L;
  path ell;
  b = 35;
  k = 0.85;
  theta = 30;
  ste = 10;
  N = 5;
  X = k*b*right;
  Y = b*dir(90-theta);
  ell = for i=ste step ste until 360: (X*cosd(i)+Y*sind(i)).. endfor cycle;
  draw ell;
  draw Y--origin--X;
  alpha = 0.5*angle((k**2)-1,2*k*sind(theta));
  alphb = 0.5*angle(1-(k**2),2*k*sind(-theta));
  A = X*cosd(alpha)+Y*sind(alpha);
  B = X*cosd(alphb)+Y*sind(alphb);
  draw B--origin--A;
  rhoa = abs(A);
  gamma = angle(A);
  rhob = abs(B);
  gammb = angle(B);
  kt = rhoa/rhob;
  for i=1 upto N:
    aux := rhoa+(i-0.5)*k*b;
    fx := aux*cosd(gamma);
    fy := aux*sind(-gamma);
    cp := (rhob*fx/kt+fy*((fx/kt)++fy+-+rhob))/(((fx/kt)**2)+(fy**2));
    cm := (rhob*fx/kt-fy*((fx/kt)++fy+-+rhob))/(((fx/kt)**2)+(fy**2));
    U := A*cm-B*(1+-+cm);
    L := A*cp+B*(1+-+cp);
    draw U--L--(aux,0)--cycle withcolor blue;
    draw 0.5[U,L];
  endfor;
endfig;

end.