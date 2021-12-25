% sharpraytrace.mp
% L. Nobre G.
% IYP (2005)

input featpost3Dplus2D;

f := (5,1.5,4.5);
Spread := 60;
PrintStep := 5;

numeric phi;
phi = 0.5*(1+sqrt(5));

V1 := ( 1, phi,0);
V2 := (-1, phi,0);
V3 := (-1,-phi,0);
V4 := ( 1,-phi,0);
V5 := (0, 1, phi);
V6 := (0,-1, phi);
V7 := (0,-1,-phi);
V8 := (0, 1,-phi);
V9 := ( phi,0, 1);
V10:= ( phi,0,-1);
V11:= (-phi,0,-1);
V12:= (-phi,0, 1);
makeface1(1,2,3,4);
makeface2(5,6,7,8);
makeface3(9,10,11,12);

beginfig(1);
  pickup pencircle scaled 2pt;
%  numeric asize;
%  asize = 0.5*(phi+1);
%  cartaxes(asize,asize,asize);
  sharpraytrace;
endfig;

Spread := 120;
f:=3*(1,0.5,0.75);
PrintStep := 0.15;

numeric a, b, c, d, e;
a = 0.1;
b = 0.3;
c = 1.1;
d = 2.9;
e = 0.5;

V1 := ( a ,  b , c);
V2 := ( a-1, b , 0);
V3 := ( a-1-e, b-d,0);
V4 := (-b ,  a , c);
V5 := (-b , a-1, 0);
V6 := (-b+d,a-1-e, 0);
V7 := (-a , -b , c);
V8 := (-a+1,-b , 0);
V9 := (-a+1+e,-b+d,0);
V10:= ( b , -a , c);
V11:= ( b ,-a+1, 0);
V12:= (b-d,-a+1+e, 0);
makeface1(1,2,3);
makeface2(4,5,6);
makeface3(7,8,9);
makeface4(10,11,12);

beginfig(2);
  pickup pencircle scaled 2pt;
  sharpraytrace;
endfig;

Spread := 120;
f:=(3,5,-0.6);
PrintStep := 5;

numeric a, b, c, d;
a = 0.2;
b = 0.3;
c = 1.3;
d = 1.4;

V1 := ( b, 2.5*a, b);
V2 := ( b,-2.5*a, b);
V3 := (-c,-2.5*a,-3*b);
V4 := (-c, 2.5*a,-3*b);
V5 := (-b, 4*a, b);
V6 := (-b,-4*a, b);
V7 := ( c,-4*a,-3*b);
V8 := ( c, 4*a,-3*b);
V9 := ( d, 2*a,-b);
V10:= ( d,-2*a,-b);
V11:= (-d,-2*a,-b);
V12:= (-d, 2*a,-b);
makeface1(1,2,3,4);
makeface2(5,6,7,8);
makeface3(9,10,11,12);

beginfig(3);
  pickup pencircle scaled 2pt;
  sharpraytrace;
endfig;

Spread := 60;
f:=7*(1,0.45,0.95);
PrintStep := 0.25;

numeric a, b, c;
a = 0.0;
b = 0.6;
c = 2.0;

V1 := ( 3,-1, a);
V2 := (-3,-1, b);
V3 := (-3,-1, c);
V4 := ( 1, 3, a);
V5 := ( 1,-3, b);
V6 := ( 1,-3, c);
V7 := (-3, 1, a);
V8 := ( 3, 1, b);
V9 := ( 3, 1, c);
V10:= (-1,-3, a);
V11:= (-1, 3, b);
V12:= (-1, 3, c);
makeface1(1,2,3);
makeface2(4,5,6);
makeface3(7,8,9);
makeface4(10,11,12);

beginfig(4);
  pickup pencircle scaled 2pt;
  sharpraytrace;
endfig;

end;
