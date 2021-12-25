% cubicstructurefcc.mp
% L. Nobre G. 
% 2011

input featpost3Dplus2D;

prologues := 1;

f := 2*(5,3,3.5);
Spread := 20;

beginfig(1);
  color Ve[], va, vb, vc, vu;
  numeric i, j, k, e[]p[];
  string comm, comn, como, a, b;
  Ve1  = (2,1,-1);
  Ve2  = (0,1,-1);
  Ve3  = (-2,1,-1);
  Ve4  = (-2,-1,-1);
  Ve5  = (0,-1,-1);
  Ve6  = (2,-1,-1);
  Ve7  = (1,1,0);
  Ve8  = (-1,1,0);
  Ve9  = (-2,0,0);
  Ve10 = (-1,-1,0);
  Ve11 = (1,-1,0);
  Ve12 = (2,0,0);
  Ve13 = (2,1,1);
  Ve14 = (0,1,1);
  Ve15 = (-2,1,1);
  Ve16 = (-2,-1,1);
  Ve17 = (0,-1,1);
  Ve18 = (2,-1,1);
  Ve19 = (1,0,1);
  Ve20 = (0,0,0);
  Ve21 = (-1,0,1);
  Ve22 = (1,0,-1);
  Ve23 = (-1,0,-1);
  for i=1 upto 23:
    comm := "unfill rigorousfearpath(" & cstr(Ve[i]) & ",0.08);";
    comn := "draw rigorousfearpath(" & cstr(Ve[i]) & ",0.08);";
    getready( comm & comn, Ve[i] );
  endfor; 
  e1p1 = 1; e2p1 = 2; e3p1 = 3; e4p1 = 4; e5p1 = 5; e6p1 = 6;
  e1p2 = 2; e2p2 = 3; e3p2 = 4; e4p2 = 5; e5p2 = 6; e6p2 = 1;
  e7p1 = 13; e8p1 = 14; e9p1 = 15; e10p1 = 16; e11p1 = 17; e12p1 = 18;
  e7p2 = 14; e8p2 = 15; e9p2 = 16; e10p2 = 17; e11p2 = 18; e12p2 = 13;
  e13p1 = 1; e14p1 = 2; e15p1 = 3; e16p1 = 4; e17p1 = 5; e18p1 = 6;
  e13p2 = 13; e14p2 = 14; e15p2 = 15; e16p2 = 16; e17p2 = 17; e18p2 = 18;
  e19p1 = 2; e19p2 = 5; e20p1 = 14; e20p2 = 17;
  for i=1 upto 20:
    vc := Ve[e[i]p2]-Ve[e[i]p1];
    vu := N(vc);
    va := Ve[e[i]p1]+0.2*vu;
    vb := Ve[e[i]p2]-0.2*vu;
    a := cstr( va );
    b := cstr( vb-va );
    comm := "rigorousdisc( 0, true," & a & ",0.08," & b & ");";
    getready( comm, 0.5[va,vb] );
  endfor;
  doitnow;
endfig;

end;
