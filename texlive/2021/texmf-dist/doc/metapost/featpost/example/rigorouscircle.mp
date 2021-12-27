% rigorouscircle.mp
% L. Nobre G. 
% 2003

input featpost3Dplus2D;

f := (2,0.2,0.1);


beginfig(1);
  numeric ind, cray, sray, cnum, ang, aux, auxs;
  color sidevec, rayvec, angulmomv, vecaux;
  path paux;
  pen pendot;
  pendot = pencircle scaled 2pt;
  sray = 1;
  cnum = 8;
  ang = 45;
  auxs = 0.15;
  ray = 3*cosd(ang)*sray/cnum;
  sidevec = ( sray*sind(ang), 0, 0 );
  for ind=1 upto cnum:
    aux := ind*360/cnum;
%    draw rigorouscircle( black, (0,sind(aux),cosd(aux)), sray ) dashed evenly;
    rayvec := sray*cosd(ang)*(0,cosd(aux),sind(aux));

    vecaux := rayvec+sidevec;
    drawarrow rp(black)--rp(vecaux+auxs*N(vecaux));
    paux := rigorouscircle( vecaux, vecaux, ray );
    draw paux;
    drawoptions( withpen pendot );
    drawarrow (subpath (0,1) of paux);
    draw rp(vecaux);
    drawoptions();

    vecaux := rayvec-sidevec;
    drawarrow rp(black)--rp(vecaux+auxs*N(vecaux));
    paux := rigorouscircle( vecaux, vecaux, ray );
    draw paux;
    drawoptions( withpen pendot );
    drawarrow (subpath (0,1) of paux);
    draw rp(vecaux);
    drawoptions();
  endfor;
endfig;

end;


