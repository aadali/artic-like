% fakehole.mp
% L. Nobre G.
% IYP (2005)

input featpost3Dplus2D;

Spread := 23;
f := (13,6,4.5);
ShadowOn := true;
LightSource := 10*(4,-1.7,2);
beginfig(1);
  numeric bigray, smaray, posray, posang, disang, aresiz, dislen;
  numeric numl, numr, holang, i;
  color pos, disaxe, holaxe;
  aresiz := 13;
  bigray := 2;
  smaray := 0.25;
  posray := 4;
  posang := -10;
  disang := 15;
  dislen := 0.4;
  holang := 115;
  numl := 45;
  disaxe := dislen*(cosd(disang),sind(disang),0);
  pos := posray*(cosd(posang),sind(posang),0);
  holaxe := (cosd(disang+90)*cosd(holang),
             sind(disang+90)*cosd(holang),
                             sind(holang));
  HoriZon := -bigray;
  setthestage( 0.5*numl, aresiz );
  rigorousdisc( smaray, true, pos, bigray, disaxe );
  numr := floor(bigray/smaray);
  for i=2-numr step 3 until -3:
    fakehole( pos+i*smaray*holaxe, disaxe, smaray );
    fakehole( pos-i*smaray*holaxe, disaxe, smaray );
  endfor;
endfig;

ShadowOn := false;
Spread := 25;
ShiftV := 105.00mm*(0.5,1);
beginfig(2);
  numeric len, ray;
  len = 10;
  ray = 0.9;
  f := (5,2*ray,1+1.5*ray); 
  viewcentr := (0,2.5*ray,1+0.5*ray);
  drawoptions( withpen pencircle scaled 2.3pt );
  kindofcube(false,true,black,0,0,0,1,len*ray,2);
  numeric i;                          
  for i=2 step 2 until len-1:
    fakehole( (0,i*ray,1), red, ray );
  endfor;
endfig;
  
end.



