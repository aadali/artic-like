% acmaglev.mp
% L. Nobre G.
% IYP (2005)

input featpost3Dplus2D;

Spread := 25;
f := (0,-2,0);
ParallelProj := true;
ShiftV := ( 0, 1 );

def vecfunc( expr posit )=
  begingroup
    save a, c, field;
    color field;
    numeric a, c;
    a = (sind(X(posit))-Z(posit)*0.95)/40;
    c = -Y( posit ) ;
    field = ( 1, c, a );
    ( field )
  endgroup
enddef;
  
beginfig(1);
  numeric i, tot;
  path dephased, sinusoid, zeroaxis;
  transform T;
  tot = 900;
  T = identity yscaled 2 xscaled 1/36 shifted (0,Shifts);
  zeroaxis = rp( black )--rp( (tot,0,0) );
  draw zeroaxis transformed T withcolor green;
  sinusoid = rp( black ) for i=1 upto tot: ...rp( (i,0,sind(i)) ) endfor;
  draw sinusoid transformed T;
  dephased = fieldlinepath(tot,(0,0,1),1,vecfunc);
  draw dephased transformed T withcolor red;
endfig;

end.
