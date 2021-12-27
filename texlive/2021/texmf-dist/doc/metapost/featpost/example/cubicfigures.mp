% cubicfigures.mp
% L. Nobre G. 
% 2003

input featpost3Dplus2D;

f := ( 1.2 , 2.0 , 1.6 );
Spread := 75;

V1 := (1,1,1);
V2 := (-1,1,1);
V3 := (-1,-1,1);
V4 := (1,-1,1);
V5 := (1,1,-1);
V6 := (-1,1,-1);
V7 := (-1,-1,-1);
V8 := (1,-1,-1);
makeface1(1,2,3,4);makeface2(5,6,7,8);
makeface3(1,2,6,5);makeface4(2,3,7,6);
makeface5(3,4,8,7);makeface6(4,1,5,8);
makeline1(1,7);makeline2(2,8);
makeline3(3,5);makeline4(4,6);

beginfig(1);
    ParallelProj := true;
    SphericalDistortion := false;
    pickup pencircle scaled 1pt;
    draw_all_test(true);
endfig;

beginfig(2);
    ParallelProj := false;
    SphericalDistortion := false;
    draw_all_test(true);
endfig;

beginfig(3);
    ParallelProj := false;
    SphericalDistortion := true;
    PrintStep := 5;
    draw_all_test(true);
endfig;

end;
