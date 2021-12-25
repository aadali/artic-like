% vanallenbelt.mp
% L. Nobre G.
% 2004

input featpost3Dplus2D;

def magnfield( expr pos ) =
  begingroup
    save dist, anga, angb, angp, bv, bh, br, bx, by, bz, fin;
    numeric dist, anga, angb, bv, bh, br, bx, by, bz;
    pair angp;
    color fin;
    dist = conorm( pos );
    angp = getanglepair( pos );
    anga = xpart( angp );
    angb = ypart( angp );
    bh = 0.5*cosd( angb )/dist**3;
    bv = sind( angb )/dist**3;
    br = bv*cosd( angb ) - bh*sind( angb );
    bz = bv*sind( angb ) + bh*cosd( angb );
    bx = br*cosd( anga );
    by = br*sind( anga );
    fin = ( bx, by, bz );
    ( 25*fin )
  endgroup
enddef;

def pmf( expr pos ) =
  (2*red/(((greenpart pos) ++ (bluepart pos))**2))
enddef;

beginfig(1);
  color vel;
  vel = 1.1*(-0.25,0.0,0.1);
  cartaxes(1,1,1);
  draw magnetictrajectorypath(200, red, vel, 0.22, magnfield );
%  draw magnetictrajectorypath(400, red,-vel, 0.22, magnfield );
endfig;

beginfig(2);
  f := 2*red;
  Spread := 200;
  numeric ang;
  color vel;
  ang = 145;
  vel = 0.8*( 0, cosd( ang ), sind( ang ) );
  draw magnetictrajectorypath(200, green, vel, 0.1, pmf );
  drawarrow reverse (rp(blue)--rp(black)--rp(green));
endfig;

end.