% closingbox.mp
% L. Nobre G.
% 2006

input featpost3Dplus2D;

beginfig(1);
  f := 1.5*(5,-6,2);
  Spread := 80;
  numeric ang;
  color auxa, auxb, auxc, auxd;
  ang = 48;
  V1 := (-1,1,0);
  V3 := (-1,-1,0);
  V2 := (-3,1,0);
  V2 := V1 + rotvecaroundanother( ang, V2-V1, V1-V3 );
  V4 := (-3,-1,0);
  V4 := V3 + rotvecaroundanother( ang, V4-V3, V1-V3 ); 
  V6 := (-1,-3,0);
  V6 := V3 + rotvecaroundanother( ang, V6-V3, V4-V3 );
  V7 := V4 + V3 - V1;
  V7 := V4 + rotvecaroundanother( ang, V7-V4, V4-V3 );
  auxa = N( V4-V3 );
  V5 := V4 + auxa;
  V5 := V4 + rotvecaroundanother( ang, auxa, V4-V7 );
  auxb = N( V7-V6 );
  V8 := V7 + auxb;
  V8 := V7 + rotvecaroundanother( ang, auxb, V4-V7 ); 

  V11:= (1,-1,0);
  V13:= (1,1,0);
  V12 := (3,-1,0);
  V12 := V11 + rotvecaroundanother( ang, V12-V11, V11-V13 );
  V14 := (3,1,0);
  V14 := V13 + rotvecaroundanother( ang, V14-V13, V11-V13 ); 
  V16 := (1,3,0);
  V16 := V13 + rotvecaroundanother( ang, V16-V13, V14-V13 );
  V17 := V14 + V13 - V11;
  V17 := V14 + rotvecaroundanother( ang, V17-V14, V14-V13 ); 
  auxc = N( V14-V13 );
  V15 := V14 + auxc;
  V15 := V14 + rotvecaroundanother( ang, auxc, V14-V17 );
  auxd = N( V17-V16 );
  V18 := V17 + auxd;
  V18 := V17 + rotvecaroundanother( ang, auxd, V14-V17 ); 
  
  makeface1(1,3,11,13);
  makeface2(1,2,4,3);
  makeface3(3,4,7,6);
  makeface4(4,5,8,7);
  makeface5(11,12,14,13);
  makeface6(13,14,17,16);
  makeface7(14,15,18,17);

  draw_invisible( false, false, TableC3, black );

endfig;

end.