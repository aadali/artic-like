% pathofstraightline.mp
% L. Nobre G.
% 2003

input featpost3Dplus2D;

SphericalDistortion := true;
Spread := 100;
f := 0.4*(1.5,0.5,1);

beginfig(1);
  numeric gridstep, sidenumber, i, coord;
  color pa, pb, pc, pd;
  gridstep = 0.1;
  sidenumber = 10;
  coord = 0.5*sidenumber*gridstep;
  for i=0 upto sidenumber:
    pa := (-coord,-coord+i*gridstep,0);
    pb := (coord,-coord+i*gridstep,0);
    pc := (-coord+i*gridstep,-coord,0);
    pd := (-coord+i*gridstep,coord,0);
    draw pathofstraightline( pa, pb );
    draw pathofstraightline( pc, pd );
  endfor;
endfig;

end.

  