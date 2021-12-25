% ultraimprovertex.mp
% L. Nobre G.
% 2012

input featpost3Dplus2D;

beginfig(1);
  numeric r, i, ste;
  r = 0.65;
  ste = 4;
  color cc, cl, sc, sl, pdir, inip;
  cl = (3*r,0,0);
  cc = (-1.5*r,0,-0.5*r);
  sl = (0,1.7*r,0);
  sc = (0,0,0.5*r);
  rigorousdisc( 0, true, cc, r, cl);
  spheroid( sc, sl, r );
  drawoptions( withpen pencircle scaled 3pt );
  for i=ste step ste until 360:
    pdir := (-sind(i),cosd(i),0);
    inip := r*(cosd(i),sind(i),0);
    draw rp(ultraimprovertex( black, pdir, cc, r, cl, sc, sl, r, inip ));
  endfor;
endfig;

beginfig(2);
  numeric r, i, ste;
  r = 0.65;
  ste = 4;
  color cc, cl, sc, sl, pdir, inip;
  cl = (3.5*r,0,0);
  cc = (-2*r,0,-0.5*r);
  sl = (0,1.9*r,0);
  sc = (0,0,0.95*r);
  spheroid( cc, cl, r );
  spheroid( sc, sl, r );
  drawoptions( withpen pencircle scaled 3pt );
  draw for i=ste step ste until 360:
    hide( pdir := (-sind(i),cosd(i),0); inip := r*(cosd(i),sind(i),0); )
    rp(necplusimprovertex( black, pdir, cc, cl, r, sc, sl, r, inip ))..
  endfor cycle;
endfig;

beginfig(3);
  numeric r, i, ste;
  r = 0.65;
  ste = 4;
  color cc, cl, sc, sl, pdir, inip;
  cl = (3*r,0,0);
  cc = (-1.5*r,0,-0.5*r);
  sl = (0.7*r,0,0);
  sc = (0,0,0.5*r);
  rigorousdisc( 0, true, cc, r, cl);
  spheroid( sc, sl, r );
  drawoptions( withpen pencircle scaled 3pt );
  for i=ste step ste until 360:
    pdir := (-sind(i),cosd(i),0);
    inip := r*(cosd(i),sind(i),0);
    draw rp(ultraimprovertex( black, pdir, cc, r, cl, sc, sl, r, inip ));
  endfor;
endfig;

end.
